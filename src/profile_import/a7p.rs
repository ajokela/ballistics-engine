//! ArcherBC2 `.a7p` file parser.
//!
//! File layout: 32 ASCII hex characters (MD5 of the remainder) followed by a
//! proto3 `Payload { Profile profile = 1; }` message. Field numbers and
//! fixed-point scale factors below were confirmed empirically against a real
//! ArcherBC2 export; they are interoperability facts, not vendored schema
//! (the upstream a7p project is LGPL-3.0; nothing from it is copied here).

use super::md5::md5_hex;
use super::wire::{
    collect_repeated_i32, parse_message, varint_to_i32, WireError, WireValue,
};

// Fixed-point scale factors (raw integer -> physical value = raw / SCALE).
const SCALE_TWIST: f64 = 100.0; // r_twist -> inches/turn
const SCALE_VELOCITY: f64 = 10.0; // c_muzzle_velocity, coef mv -> m/s
const SCALE_DIMENSION: f64 = 1000.0; // b_diameter, b_length -> inches
const SCALE_WEIGHT: f64 = 10.0; // b_weight -> grains
const SCALE_PRESSURE: f64 = 10.0; // c_zero_air_pressure -> hPa
const SCALE_COEF: f64 = 10_000.0; // coef bc_cd -> BC or Cd
const SCALE_TCOEFF: f64 = 1000.0; // c_t_coeff -> % per 15 C
const SCALE_DISTANCE: f64 = 100.0; // distances -> meters
// CUSTOM coef rows reuse the `bc_cd` field for Cd (same SCALE_COEF as BC/Cd) but the `mv`
// field means Mach number, not m/s, and uses its own scale (empirically confirmed, same
// footing as the other scale factors above — see the module doc).
const SCALE_MACH: f64 = 10_000.0; // coef mv (CUSTOM only) -> Mach number

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum A7pBcType {
    G1,
    G7,
    Custom,
    /// Forward compatibility: an enum value this build does not know.
    Other(i32),
}

#[derive(Debug, Clone, PartialEq)]
pub enum EnvelopeStatus {
    Verified,
    Mismatch { expected: String, actual: String },
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UnknownField {
    pub context: &'static str,
    pub number: u32,
}

#[derive(Debug, Clone)]
pub struct A7pProfile {
    pub profile_name: String,
    pub cartridge_name: String,
    pub bullet_name: String,
    pub short_name_top: String,
    pub short_name_bot: String,
    pub user_note: String,
    pub caliber: String,
    pub device_uuid: String,
    pub zero_x_raw: i32,
    pub zero_y_raw: i32,
    pub w_pitch_raw: i32,
    pub sight_height_mm: f64,
    pub twist_in_per_turn: f64,
    pub twist_right: bool,
    pub muzzle_velocity_mps: f64,
    pub zero_temperature_c: f64,
    pub temp_coeff_pct_per_15c: f64,
    pub powder_temperature_c: f64,
    pub air_temperature_c: f64,
    pub air_pressure_hpa: f64,
    pub air_humidity_pct: f64,
    pub bullet_diameter_in: f64,
    pub bullet_weight_gr: f64,
    pub bullet_length_in: f64,
    pub bc_type: A7pBcType,
    /// Raw (bc_cd, mv) integer pairs; scaling depends on bc_type (see
    /// `bc_rows()` / `custom_rows()`).
    pub coef_rows_raw: Vec<(i32, i32)>,
    pub distances_m: Vec<f64>,
    pub zero_distance_m: Option<f64>,
    pub switches_count: usize,
}

impl Default for A7pProfile {
    fn default() -> Self {
        A7pProfile {
            profile_name: String::new(),
            cartridge_name: String::new(),
            bullet_name: String::new(),
            short_name_top: String::new(),
            short_name_bot: String::new(),
            user_note: String::new(),
            caliber: String::new(),
            device_uuid: String::new(),
            zero_x_raw: 0,
            zero_y_raw: 0,
            w_pitch_raw: 0,
            sight_height_mm: 0.0,
            twist_in_per_turn: 0.0,
            twist_right: true, // proto3 default TwistDir::RIGHT = 0
            muzzle_velocity_mps: 0.0,
            zero_temperature_c: 0.0,
            temp_coeff_pct_per_15c: 0.0,
            powder_temperature_c: 0.0,
            air_temperature_c: 0.0,
            air_pressure_hpa: 0.0,
            air_humidity_pct: 0.0,
            bullet_diameter_in: 0.0,
            bullet_weight_gr: 0.0,
            bullet_length_in: 0.0,
            bc_type: A7pBcType::G1,
            coef_rows_raw: Vec::new(),
            distances_m: Vec::new(),
            zero_distance_m: None,
            switches_count: 0,
        }
    }
}

impl A7pProfile {
    /// G1/G7 interpretation: (BC, velocity m/s), file order preserved.
    pub fn bc_rows(&self) -> Vec<(f64, f64)> {
        self.coef_rows_raw
            .iter()
            .map(|&(bc, mv)| (f64::from(bc) / SCALE_COEF, f64::from(mv) / SCALE_VELOCITY))
            .collect()
    }

    /// CUSTOM interpretation: (Cd, Mach), file order preserved. Only meaningful when
    /// `bc_type == A7pBcType::Custom` — for G1/G7 files the same raw rows mean
    /// (BC, velocity m/s), see `bc_rows()`.
    pub fn custom_rows(&self) -> Vec<(f64, f64)> {
        self.coef_rows_raw
            .iter()
            .map(|&(cd, mv)| (f64::from(cd) / SCALE_COEF, f64::from(mv) / SCALE_MACH))
            .collect()
    }
}

#[derive(Debug)]
pub enum A7pError {
    /// Shorter than the 32-byte envelope prefix.
    TooShort,
    /// Envelope prefix is not ASCII hex.
    BadPrefix,
    /// Payload did not contain a Profile message.
    MissingProfile,
    /// Malformed protobuf payload. Carries a pre-formatted message rather
    /// than the crate-internal `WireError` type, which is `pub(crate)` to
    /// the `profile_import` subtree and must not leak into this public enum
    /// (would otherwise trip `private_interfaces`).
    Wire(String),
}

impl std::fmt::Display for A7pError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            A7pError::TooShort => write!(f, "file too short to be a .a7p (needs 32-byte checksum prefix)"),
            A7pError::BadPrefix => write!(f, "checksum prefix is not ASCII hex — not a .a7p file"),
            A7pError::MissingProfile => write!(f, "payload contains no profile message"),
            A7pError::Wire(e) => write!(f, "malformed protobuf payload: {e}"),
        }
    }
}

impl std::error::Error for A7pError {}

impl From<WireError> for A7pError {
    fn from(e: WireError) -> Self {
        A7pError::Wire(e.to_string())
    }
}

#[derive(Debug)]
pub struct A7pDocument {
    pub profile: A7pProfile,
    pub envelope: EnvelopeStatus,
    pub unknown_fields: Vec<UnknownField>,
}

/// Wrap a serialized Payload in the .a7p envelope (MD5 hex prefix). The
/// inverse of the envelope check in [`parse_a7p`]; used by tests today and by
/// .a7p export later.
pub fn wrap_payload(payload: &[u8]) -> Vec<u8> {
    let mut out = md5_hex(payload).into_bytes();
    out.extend_from_slice(payload);
    out
}

pub fn parse_a7p(bytes: &[u8]) -> Result<A7pDocument, A7pError> {
    if bytes.len() < 32 {
        return Err(A7pError::TooShort);
    }
    let (prefix, payload) = bytes.split_at(32);
    let expected = std::str::from_utf8(prefix)
        .ok()
        .filter(|s| s.chars().all(|c| c.is_ascii_hexdigit()))
        .ok_or(A7pError::BadPrefix)?
        .to_ascii_lowercase();
    let actual = md5_hex(payload);
    let envelope = if expected == actual {
        EnvelopeStatus::Verified
    } else {
        EnvelopeStatus::Mismatch { expected, actual }
    };

    let mut unknown_fields = Vec::new();
    let payload_fields = parse_message(payload)?;
    let mut profile_bytes: Option<&[u8]> = None;
    for field in &payload_fields {
        match (field.number, &field.value) {
            (1, WireValue::Bytes(b)) => profile_bytes = Some(b),
            _ => unknown_fields.push(UnknownField {
                context: "Payload",
                number: field.number,
            }),
        }
    }
    let profile_bytes = profile_bytes.ok_or(A7pError::MissingProfile)?;

    let fields = parse_message(profile_bytes)?;
    let mut p = A7pProfile::default();
    let mut zero_idx: i32 = 0;

    let get_str = |v: &WireValue| -> String {
        match v {
            WireValue::Bytes(b) => String::from_utf8_lossy(b).into_owned(),
            _ => String::new(),
        }
    };

    for field in &fields {
        let int = |v: &WireValue| -> i32 {
            match v {
                WireValue::Varint(raw) => varint_to_i32(*raw),
                _ => 0,
            }
        };
        match field.number {
            1 => p.profile_name = get_str(&field.value),
            2 => p.cartridge_name = get_str(&field.value),
            3 => p.bullet_name = get_str(&field.value),
            4 => p.short_name_top = get_str(&field.value),
            5 => p.short_name_bot = get_str(&field.value),
            6 => p.user_note = get_str(&field.value),
            7 => p.zero_x_raw = int(&field.value),
            8 => p.zero_y_raw = int(&field.value),
            9 => p.sight_height_mm = f64::from(int(&field.value)),
            10 => p.twist_in_per_turn = f64::from(int(&field.value)) / SCALE_TWIST,
            11 => p.muzzle_velocity_mps = f64::from(int(&field.value)) / SCALE_VELOCITY,
            12 => p.zero_temperature_c = f64::from(int(&field.value)),
            13 => p.temp_coeff_pct_per_15c = f64::from(int(&field.value)) / SCALE_TCOEFF,
            14 => zero_idx = int(&field.value),
            15 => p.air_temperature_c = f64::from(int(&field.value)),
            16 => p.air_pressure_hpa = f64::from(int(&field.value)) / SCALE_PRESSURE,
            17 => p.air_humidity_pct = f64::from(int(&field.value)),
            18 => p.w_pitch_raw = int(&field.value),
            19 => p.powder_temperature_c = f64::from(int(&field.value)),
            20 => p.bullet_diameter_in = f64::from(int(&field.value)) / SCALE_DIMENSION,
            21 => p.bullet_weight_gr = f64::from(int(&field.value)) / SCALE_WEIGHT,
            22 => p.bullet_length_in = f64::from(int(&field.value)) / SCALE_DIMENSION,
            23 => p.twist_right = int(&field.value) == 0,
            24 => {
                p.bc_type = match int(&field.value) {
                    0 => A7pBcType::G1,
                    1 => A7pBcType::G7,
                    2 => A7pBcType::Custom,
                    other => A7pBcType::Other(other),
                }
            }
            25 => p.switches_count += 1,
            26 => {} // handled below via collect_repeated_i32 (packed or not)
            27 => {
                if let WireValue::Bytes(b) = &field.value {
                    let row = parse_message(b)?;
                    let mut bc_cd = 0i32;
                    let mut mv = 0i32;
                    for rf in &row {
                        match (rf.number, &rf.value) {
                            (1, WireValue::Varint(raw)) => bc_cd = varint_to_i32(*raw),
                            (2, WireValue::Varint(raw)) => mv = varint_to_i32(*raw),
                            _ => unknown_fields.push(UnknownField {
                                context: "CoefRow",
                                number: rf.number,
                            }),
                        }
                    }
                    p.coef_rows_raw.push((bc_cd, mv));
                }
            }
            28 => p.caliber = get_str(&field.value),
            29 => p.device_uuid = get_str(&field.value),
            other => unknown_fields.push(UnknownField {
                context: "Profile",
                number: other,
            }),
        }
    }

    p.distances_m = collect_repeated_i32(&fields, 26)?
        .into_iter()
        .map(|d| f64::from(d) / SCALE_DISTANCE)
        .collect();
    p.zero_distance_m = usize::try_from(zero_idx)
        .ok()
        .and_then(|idx| p.distances_m.get(idx).copied());
    // proto3: absent c_zero_distance_idx means 0 == first distance, but an
    // empty distances list means there is nothing to point at.

    Ok(A7pDocument {
        profile: p,
        envelope,
        unknown_fields,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    // Test-only encoders, spec-derived (same rationale as wire::tests).
    fn enc_varint(mut v: u64, out: &mut Vec<u8>) {
        loop {
            let byte = (v & 0x7f) as u8;
            v >>= 7;
            if v == 0 {
                out.push(byte);
                return;
            }
            out.push(byte | 0x80);
        }
    }
    fn enc_i32(number: u32, value: i64, out: &mut Vec<u8>) {
        enc_varint(u64::from(number) << 3, out);
        enc_varint(value as u64, out);
    }
    fn enc_str(number: u32, s: &str, out: &mut Vec<u8>) {
        enc_varint((u64::from(number) << 3) | 2, out);
        enc_varint(s.len() as u64, out);
        out.extend_from_slice(s.as_bytes());
    }
    fn enc_bytes(number: u32, payload: &[u8], out: &mut Vec<u8>) {
        enc_varint((u64::from(number) << 3) | 2, out);
        enc_varint(payload.len() as u64, out);
        out.extend_from_slice(payload);
    }

    /// A synthetic profile using the values of a real-world .338 LM export
    /// (300 gr OTM) — the values are facts about the FORMAT's scaling, not
    /// copied file bytes.
    fn synthetic_profile_bytes() -> Vec<u8> {
        let mut p = Vec::new();
        enc_str(1, "338LM 300GR", &mut p);
        enc_str(2, "LAPUA MAGNUM", &mut p);
        enc_str(3, "300GR OTM", &mut p);
        enc_i32(9, 90, &mut p); // sc_height 90 mm
        enc_i32(10, 1000, &mut p); // r_twist 10.00 in
        enc_i32(11, 7920, &mut p); // 792.0 m/s
        enc_i32(12, 15, &mut p);
        enc_i32(13, 1000, &mut p); // 1.000 %/15C
        enc_i32(14, 1, &mut p); // zero at distances[1]
        enc_i32(15, 15, &mut p);
        enc_i32(16, 10000, &mut p); // 1000.0 hPa
        enc_i32(17, 50, &mut p);
        enc_i32(20, 338, &mut p); // 0.338 in
        enc_i32(21, 3000, &mut p); // 300.0 gr
        enc_i32(22, 1800, &mut p); // 1.800 in
        enc_i32(23, 1, &mut p); // LEFT twist (non-default on purpose)
        enc_i32(24, 1, &mut p); // G7 (non-default on purpose)
        // distances, packed: 100.00 m and 200.00 m
        let mut packed = Vec::new();
        enc_varint(10_000, &mut packed);
        enc_varint(20_000, &mut packed);
        enc_bytes(26, &packed, &mut p);
        // two BC rows
        let mut row1 = Vec::new();
        enc_i32(1, 3810, &mut row1); // bc 0.381
        enc_i32(2, 7920, &mut row1); // at 792.0 m/s
        enc_bytes(27, &row1, &mut p);
        let mut row2 = Vec::new();
        enc_i32(1, 3600, &mut row2);
        enc_i32(2, 5000, &mut row2);
        enc_bytes(27, &row2, &mut p);
        enc_str(28, ".338 Lapua Magnum", &mut p);
        // an unknown future field the parser must survive AND report
        enc_i32(63, 42, &mut p);

        let mut payload = Vec::new();
        enc_bytes(1, &p, &mut payload); // Payload.profile = 1
        payload
    }

    #[test]
    fn parses_synthetic_profile_with_verified_envelope() {
        let file = wrap_payload(&synthetic_profile_bytes());
        let doc = parse_a7p(&file).expect("parse");
        assert!(matches!(doc.envelope, EnvelopeStatus::Verified));
        let p = &doc.profile;
        assert_eq!(p.profile_name, "338LM 300GR");
        assert_eq!(p.bullet_name, "300GR OTM");
        assert!((p.sight_height_mm - 90.0).abs() < 1e-9);
        assert!((p.twist_in_per_turn - 10.0).abs() < 1e-9);
        assert!(!p.twist_right); // LEFT
        assert!((p.muzzle_velocity_mps - 792.0).abs() < 1e-9);
        assert!((p.temp_coeff_pct_per_15c - 1.0).abs() < 1e-9);
        assert!((p.air_pressure_hpa - 1000.0).abs() < 1e-9);
        assert!((p.air_humidity_pct - 50.0).abs() < 1e-9);
        assert!((p.bullet_diameter_in - 0.338).abs() < 1e-9);
        assert!((p.bullet_weight_gr - 300.0).abs() < 1e-9);
        assert!((p.bullet_length_in - 1.8).abs() < 1e-9);
        assert!(matches!(p.bc_type, A7pBcType::G7));
        assert_eq!(p.coef_rows_raw, vec![(3810, 7920), (3600, 5000)]);
        assert_eq!(p.distances_m, vec![100.0, 200.0]);
        // zero index 1 -> 200 m
        assert_eq!(p.zero_distance_m, Some(200.0));
        // the unknown field 63 was reported, not silently dropped
        assert!(doc
            .unknown_fields
            .iter()
            .any(|u| u.context == "Profile" && u.number == 63));
    }

    #[test]
    fn proto3_defaults_apply_when_fields_are_absent() {
        // Empty profile message: everything defaults (G1, RIGHT twist, zeros).
        let mut payload = Vec::new();
        enc_bytes(1, &[], &mut payload);
        let doc = parse_a7p(&wrap_payload(&payload)).expect("parse");
        let p = &doc.profile;
        assert!(matches!(p.bc_type, A7pBcType::G1));
        assert!(p.twist_right);
        assert_eq!(p.muzzle_velocity_mps, 0.0);
        assert_eq!(p.zero_distance_m, None); // no distances to index into
        assert!(p.coef_rows_raw.is_empty());
    }

    #[test]
    fn corrupted_envelope_is_reported_not_fatal() {
        let mut file = wrap_payload(&{
            let mut payload = Vec::new();
            enc_bytes(1, &[], &mut payload);
            payload
        });
        file[0] = if file[0] == b'0' { b'1' } else { b'0' };
        let doc = parse_a7p(&file).expect("parse succeeds with warning");
        assert!(matches!(doc.envelope, EnvelopeStatus::Mismatch { .. }));
    }

    #[test]
    fn short_files_and_garbage_are_clean_errors() {
        assert!(matches!(parse_a7p(&[]), Err(A7pError::TooShort)));
        assert!(matches!(parse_a7p(&[0u8; 20]), Err(A7pError::TooShort)));
        // 32-byte prefix of non-hex garbage + non-proto payload
        let mut junk = vec![b'z'; 32];
        junk.extend_from_slice(&[0xff, 0xff, 0xff]);
        assert!(parse_a7p(&junk).is_err());
    }

    /// A synthetic CUSTOM (bc_type=2) profile: two coef rows whose raw ints are chosen so
    /// the G1/G7 scale (mv / 10) and the CUSTOM scale (mv / 10_000) would decode to visibly
    /// different numbers — catching a copy-paste of `bc_rows()`'s scale factor into
    /// `custom_rows()`.
    fn synthetic_custom_profile_bytes() -> Vec<u8> {
        let mut p = Vec::new();
        enc_str(1, "CUSTOM CURVE", &mut p);
        enc_i32(24, 2, &mut p); // CUSTOM
        let mut row1 = Vec::new();
        enc_i32(1, 5000, &mut row1); // Cd 0.5000
        enc_i32(2, 5000, &mut row1); // Mach 0.5000
        enc_bytes(27, &row1, &mut p);
        let mut row2 = Vec::new();
        enc_i32(1, 2300, &mut row2); // Cd 0.2300
        enc_i32(2, 30000, &mut row2); // Mach 3.0000
        enc_bytes(27, &row2, &mut p);
        let mut payload = Vec::new();
        enc_bytes(1, &p, &mut payload);
        payload
    }

    #[test]
    fn custom_rows_uses_the_mach_scale_not_the_velocity_scale() {
        let file = wrap_payload(&synthetic_custom_profile_bytes());
        let doc = parse_a7p(&file).expect("parse");
        let p = &doc.profile;
        assert!(matches!(p.bc_type, A7pBcType::Custom));
        assert_eq!(p.coef_rows_raw, vec![(5000, 5000), (2300, 30000)]);

        let rows = p.custom_rows();
        assert_eq!(rows.len(), 2);
        assert!((rows[0].0 - 0.5).abs() < 1e-9); // Cd
        assert!((rows[0].1 - 0.5).abs() < 1e-9); // Mach 0.5, NOT 500.0 m/s
        assert!((rows[1].0 - 0.23).abs() < 1e-9);
        assert!((rows[1].1 - 3.0).abs() < 1e-9); // Mach 3.0, NOT 3000.0 m/s

        // bc_rows() on the same raw data uses the OTHER scale (m/s /10, not Mach /10_000) —
        // demonstrates the two interpretations are genuinely distinct, not aliases.
        let bc_interpretation = p.bc_rows();
        assert!((bc_interpretation[0].1 - 500.0).abs() < 1e-9);
        assert!((bc_interpretation[1].1 - 3000.0).abs() < 1e-9);
    }

    #[test]
    fn custom_rows_is_empty_when_no_coef_rows_present() {
        let mut p = Vec::new();
        enc_i32(24, 2, &mut p); // CUSTOM, no rows
        let mut payload = Vec::new();
        enc_bytes(1, &p, &mut payload);
        let doc = parse_a7p(&wrap_payload(&payload)).expect("parse");
        assert!(doc.profile.custom_rows().is_empty());
    }

    #[test]
    fn out_of_range_zero_index_yields_none() {
        let mut p = Vec::new();
        enc_i32(14, 7, &mut p); // index 7, but only 1 distance
        let mut packed = Vec::new();
        enc_varint(10_000, &mut packed);
        enc_bytes(26, &packed, &mut p);
        let mut payload = Vec::new();
        enc_bytes(1, &p, &mut payload);
        let doc = parse_a7p(&wrap_payload(&payload)).expect("parse");
        assert_eq!(doc.profile.zero_distance_m, None);
        assert_eq!(doc.profile.distances_m, vec![100.0]);
    }
}
