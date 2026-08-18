//! Saved ballistic profile data model.
//!
//! `ProfileData` is the schema of `~/.ballistics/profiles/<name>.json` — the CLI's saved
//! profiles — and, since the bridge grew `profile.validate`/`profile.normalize`/
//! `profile.import_a7p`, also the profile document the JSON command bridge exchanges with
//! embedding apps. It lived in `main.rs` from its introduction; it moved here VERBATIM so
//! the bridge and the CLI share one definition. The serde wire shape is a compatibility
//! contract: every stored profile on disk and every fixture in the test suite must keep
//! loading, so field names, defaults, `skip_serializing_if` decisions, and the ABSENCE of
//! `deny_unknown_fields` (unknown keys are tolerated and dropped, the documented
//! forward-compat behavior) are all deliberately unchanged. See the round-trip test at the
//! bottom of this file before touching any attribute.
//!
//! fs-free by design: this module must compile for wasm32, so file persistence
//! (`save_profile`/`load_profile`) and unit conversion of loaded profiles (`converted_to`,
//! which rides on the CLI's `UnitConverter`) stay in `main.rs`.

use serde::{Deserialize, Serialize};

use crate::adjustment::parse_click_value;
use crate::cli_api::UnitSystem;
use crate::optic::{HoldBounds, OpticProfile, TravelLimits, TurretState};
use crate::reticle::ReticleDescription;
use crate::truing_dsf::DsfPoint;

/// Saved ballistic profile for quick recall
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProfileData {
    pub name: String,
    pub velocity: f64,
    pub bc: f64,
    pub mass: f64,
    pub diameter: f64,
    pub drag_model: String,
    #[serde(default)]
    pub twist_rate: Option<f64>,
    #[serde(default)]
    pub sight_height: Option<f64>,
    #[serde(default)]
    pub zero_distance: Option<f64>,
    #[serde(default = "default_unit_system")]
    pub units: String,
    #[serde(default = "default_temperature")]
    pub temperature: f64,
    #[serde(default = "default_pressure")]
    pub pressure: f64,
    #[serde(default = "default_humidity")]
    pub humidity: f64,
    #[serde(default)]
    pub altitude: f64,
    #[serde(default)]
    pub bullet_name: Option<String>,
    #[serde(default)]
    pub created: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub wind_speed: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub wind_direction: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub shooting_angle: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub auto_zero: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub twist_right: Option<bool>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub use_bc_segments: Option<bool>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub bullet_length: Option<f64>,
    /// Turret elevation click graduation for `--adjustment-unit clicks` (MBA-1355), e.g. "0.1mil" or
    /// "0.25moa" — parsed by `parse_click_value` at both save-time (validation) and
    /// resolve-time (`resolve_click_values`). Unit-invariant (an angular graduation, not a
    /// linear measurement), so `converted_to` leaves it untouched — see the `bc_segments`/
    /// `drag_curve` comment below for why.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub elevation_click: Option<String>,
    /// Turret windage click graduation for `--adjustment-unit clicks` (MBA-1355). Falls back to the
    /// resolved elevation graduation when unset — see `resolve_click_values`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub windage_click: Option<String>,
    /// Velocity-banded BC breakpoints (MBA-1323 Phase 2: multi-row `.a7p` G1/G7 import).
    /// `velocity_mps` in each entry is ALWAYS meters/second regardless of this profile's
    /// `units` field — see [`ProfileBcSegment`]. The scalar `bc` field above is retained as
    /// the highest-velocity row for tools that only understand a single BC; this list is the
    /// authoritative full schedule once present. `None` (the pre-Phase-2 shape) means "no
    /// velocity-banded schedule was captured" and callers fall back to the scalar `bc`.
    ///
    /// FORWARD-COMPAT WARNING (one-way): `#[serde(default)]` means this field round-trips
    /// safely through readers that predate Phase 2, but "safely" only means *deserialization*
    /// doesn't error — a pre-Phase-2 (or otherwise un-updated, e.g. stale WASM/bindings) reader
    /// silently drops this key and solves with only the scalar `bc` above. That is a materially
    /// different, unwarned answer whenever the schedule's non-fastest bands matter (empirically
    /// confirmed: ~639 m/s vs. ~411 m/s impact velocity for the same imported profile — see
    /// CLI_USAGE.md's "Multi-BC and CUSTOM drag curves" section). There is no sentinel trick
    /// available here the way there is for `drag_curve`/CUSTOM below (a real, plausible-looking
    /// `bc` value is unavoidable for back-compat with single-BC tools), so this direction of
    /// version skew degrades silently by design and must stay documented rather than "fixed".
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub bc_segments: Option<Vec<ProfileBcSegment>>,
    /// Full Mach/Cd drag curve (MBA-1323 Phase 2: `.a7p` `bc_type == CUSTOM` import). When
    /// present, the scalar `bc`/`drag_model` fields are not physically meaningful for the
    /// solve (`drag_model` reads "CUSTOM"; see `map_a7p_to_profile`'s CUSTOM handling for why
    /// `bc` is an inert `0.0` sentinel rather than a real coefficient).
    ///
    /// FORWARD-COMPAT NOTE: unlike `bc_segments` above, a reader that predates Phase 2 (or
    /// otherwise doesn't consume this field) is safe by construction, not just by omission: it
    /// still sees `bc == 0.0` and `drag_model == "CUSTOM"`, so `BallisticInputs::validate_for_solve`
    /// rejects the solve loudly ("bc_value must be finite and greater than zero") instead of
    /// silently running physics under a fabricated coefficient.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub drag_curve: Option<Vec<ProfileDragPoint>>,
    /// Mach-keyed drop-scale-factor table (MBA-1357), accumulated one point at a time by
    /// the `dsf` verb. `None` for every profile with no DSF calibration yet — including
    /// every profile saved before this field existed, which loads clean and solves
    /// untrued (same `#[serde(default)]` forward-compat pattern as `bc_segments`/
    /// `drag_curve` above: an old reader that predates this field silently drops it on
    /// re-save, degrading to untrued drop with no error).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub dsf_points: Option<Vec<DsfPoint>>,
    /// Which standard atmosphere `bc`/`bc_segments` are referenced to (MBA-1365): `None`
    /// (the omitted-field default, and every profile saved before this field existed) or
    /// `"icao"` mean ICAO; `"army-standard-metro"` declares the older Army Standard Metro
    /// reference some vendor-published BCs use instead. Parsed by
    /// `parse_bc_reference_profile_field`, written by `bc_reference_profile_field` (which
    /// never writes `"icao"` — it stays the omitted default so an untouched profile
    /// round-trips with no new key). Unit-invariant, like `bc_segments`/`drag_curve` above,
    /// so `converted_to` leaves it untouched.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub bc_reference: Option<String>,
    /// Whether `pressure` is absolute station pressure or a sea-level-corrected altimeter
    /// setting (QNH), mirroring `bc_reference` (MBA-1397): `None` (the omitted-field default,
    /// and every profile saved before this field existed) or `"absolute"` mean absolute;
    /// `"qnh"` declares a QNH pressure that must be reduced to station pressure before use.
    /// Parsed by `parse_pressure_reference_profile_field`, written by
    /// `pressure_reference_profile_field` (which never writes `"absolute"` -- it stays the
    /// omitted default so an untouched profile round-trips with no new key).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pressure_reference: Option<String>,
    /// Density altitude (MBA-1366), feet imperial / meters metric per `units` (same convention
    /// as `altitude`). `None` (the omitted-field default, and every profile saved before this
    /// field existed) means no density-altitude override is stored; a saved value supersedes
    /// `altitude`/`pressure` when the profile is loaded (see `trajectory`'s
    /// `--density-altitude` for the full precedence rule). `converted_to` rescales it exactly
    /// like `altitude` since it shares the same feet/meters convention.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub density_altitude: Option<f64>,
    /// Deliberate vertical POI offset AT the zero range (MBA-1359, Kestrel "zero height"):
    /// positive = the rifle is deliberately zeroed to impact HIGH by this much at the zero
    /// distance. ALWAYS meters regardless of this profile's `units` field (same unit-fixed
    /// convention as [`ProfileBcSegment::velocity_mps`]), so `converted_to` leaves it
    /// untouched. `None` (the omitted-field default, and every profile saved before this
    /// field existed) means no offset; an old reader silently drops it on re-save (the
    /// `bc_segments` forward-compat pattern).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub zero_poi_up_m: Option<f64>,
    /// Deliberate horizontal POI offset AT the zero range (MBA-1359, Kestrel "zero
    /// offset"): positive = impacts RIGHT. ALWAYS meters, like `zero_poi_up_m` above.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub zero_poi_right_m: Option<f64>,
    /// Lateral sight-to-bore mount offset (MBA-1396, offset-mounted optics): positive =
    /// sight RIGHT of bore. ALWAYS meters (unit-fixed like the zero POI fields above), so
    /// `converted_to` leaves it untouched; same `#[serde(default)]` forward-compat
    /// pattern (an old reader silently drops it on re-save).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sight_offset_lateral_m: Option<f64>,
    /// Elevation-axis scope tracking correction factor from a tall-target test
    /// (MBA-1358, Litz), stored as the published ratio `actual measured travel /
    /// dialed travel` (0.95 = the scope under-tracks by 5%). Elevation dial-unit
    /// outputs (mil/MOA/SMOA/IPHY/clicks) are DIVIDED by this factor — an
    /// under-tracking scope needs more dial — and dialed truing observations are
    /// MULTIPLIED by it (scope-dial -> true angular); raw drop inches never scale.
    /// NOTE on conventions: Kestrel's "Scope Cal" MULTIPLIES its factor into the
    /// solution because it stores the reciprocal (dialed/actual); we divide because we
    /// store the published actual/dialed — same physics, opposite bookkeeping.
    /// Dimensionless, so `converted_to` leaves it untouched. Validated on load: must
    /// be strictly between 0.5 and 1.5 (a factor outside that band means the
    /// tall-target test went wrong, not that the scope does). `None` = 1.0 = no
    /// correction. Derive with `ballistics tall-target`; overridden by
    /// `--elevation-cf`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub elevation_cf: Option<f64>,
    /// Windage-axis scope tracking correction factor (MBA-1358), same contract and
    /// direction as `elevation_cf`: windage-axis dial-unit outputs (including mover
    /// lead/ring) are divided by it; overridden by `--windage-cf`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub windage_cf: Option<f64>,
    /// Named zero sets (MBA-1360): alternate zero distances and per-load dial
    /// corrections (Lapua Sight-In POI / ATrag zero zones / Strelok multi-zero class).
    /// Managed by `profile zero-set add|remove|list`; selected at solve time with
    /// `--zero-set NAME`. Nothing here applies unless a set is explicitly selected —
    /// the profile's own `zero_distance`/`auto_zero` remain the master zero.
    ///
    /// FORWARD-COMPAT (the `bc_segments` pattern, deliberately): `#[serde(default)]`
    /// means a reader that predates this field loads the profile cleanly and solves
    /// with the master zero — which is exactly what a CURRENT reader does when no
    /// `--zero-set` is selected, so an old reader can never silently produce a
    /// different default answer. Requesting an alternate set on an old binary fails
    /// loudly at the flag (`--zero-set` is an unknown argument there). The one-way
    /// skew is re-SAVING: an old reader that rewrites the profile silently drops this
    /// key (documented, like `bc_segments`; there is no sentinel trick available that
    /// wouldn't corrupt the master-zero fields old readers rely on).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub zero_sets: Option<Vec<ProfileZeroSet>>,
    /// The optic's reticle (MBA-1361), so `reticle hold --profile NAME` can place a
    /// firing solution without being handed a description every time. Set with
    /// `profile save --reticle-json <file>`; carried forward untouched by a re-save.
    ///
    /// Angular data only (milliradians from the optical center), so
    /// `ProfileData::converted_to` (in `main.rs`) leaves it alone for the same reason it leaves
    /// `elevation_click` alone — a subtension is not a linear measurement.
    ///
    /// FORWARD-COMPAT (the `bc_segments` pattern): `#[serde(default)]` means a reader that
    /// predates this field loads the profile cleanly. Nothing about a trajectory depends
    /// on it — it is a display/hold aid consumed only by the `reticle` command — so an old
    /// reader cannot produce a different ballistic answer because of it; it simply has no
    /// `reticle` verb. The one-way skew is re-SAVING, which drops the key, exactly as
    /// documented for `bc_segments` and `zero_sets`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reticle: Option<ReticleDescription>,
    /// Turret mechanics and reticle hold bounds (MBA-1348): twelve fields (this one
    /// through `hold_bound_right_mil` below) assembled by `ProfileData::optic_profile`
    /// into a `ballistics_engine::optic::OpticProfile` for the dial/hold/hybrid
    /// correction planner. Every one is independently `Option` and, like
    /// `elevation_click`/`windage_click` above, ALWAYS stored in mil (or its own natural
    /// type for `clicks_per_revolution`/`zero_stop`) regardless of this profile's `units`
    /// field — angular turret/reticle geometry, not a linear measurement — so
    /// `converted_to` leaves all twelve untouched. Set with `profile save
    /// --clicks-per-rev/--zero-stop/--travel-up/--travel-down/--windage-travel-left/
    /// --windage-travel-right/--turret-elev/--turret-wind/--hold-up/--hold-down/
    /// --hold-left/--hold-right`; validated at save time via `optic_profile()` +
    /// `OpticProfile::validate()` (a profile can never be saved with, say, a dialed
    /// turret state outside its own declared travel, or a non-positive click size).
    ///
    /// FORWARD-COMPAT (the `bc_segments` pattern, deliberately): `#[serde(default)]`
    /// means a reader that predates these fields loads the profile cleanly with every one
    /// of them absent — identical to what a CURRENT reader does for a profile that never
    /// set them, so an old reader can never silently produce a different ballistic
    /// answer because of them (nothing about a trajectory solve reads them; only a later
    /// dial/hold planner does). The one-way skew is re-SAVING on an old binary, which
    /// silently drops all twelve keys, exactly as documented for `bc_segments`/
    /// `zero_sets`/`reticle` above.
    ///
    /// This field specifically: click detents per full turret revolution, for turrets
    /// whose cap marks revolutions at all (many hunting turrets do not). `None` means
    /// unknown/not applicable, never a specific count.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub clicks_per_revolution: Option<u32>,
    /// Whether the elevation turret hard-stops at its lowest travel so it cannot be
    /// dialed below zero (MBA-1348) — purely descriptive metadata, never read by
    /// `plan_corrections` (see `OpticProfile::zero_stop`'s own doc comment for why).
    /// `None` (the omitted-field default) means not recorded; `optic_profile()` treats
    /// that the same as `Some(false)`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub zero_stop: Option<bool>,
    /// Elevation travel remaining UP from the current zero (not the turret's mechanical
    /// bottom), DIAL-space mil (MBA-1348). Required together with
    /// `elevation_travel_down_mil` — `optic_profile()` returns a named-field `Err` if
    /// only one of the pair is set, rather than silently treating the unset half as zero
    /// travel (a specific, likely-false physical claim, not an honest "unknown").
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub elevation_travel_up_mil: Option<f64>,
    /// Elevation travel remaining DOWN from the current zero, DIAL-space mil (MBA-1348).
    /// See `elevation_travel_up_mil`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub elevation_travel_down_mil: Option<f64>,
    /// Windage travel remaining LEFT from the current zero, DIAL-space mil (MBA-1348) —
    /// `TravelLimits::down_mil` on the windage axis (see that type's doc comment for the
    /// left/down convention). Required together with `windage_travel_right_mil`, like
    /// `elevation_travel_up_mil`/`_down_mil` above.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub windage_travel_left_mil: Option<f64>,
    /// Windage travel remaining RIGHT from the current zero, DIAL-space mil (MBA-1348) —
    /// `TravelLimits::up_mil` on the windage axis. See `windage_travel_left_mil`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub windage_travel_right_mil: Option<f64>,
    /// The elevation turret's current dialed offset from zero, DIAL-space mil, signed:
    /// positive is dialed UP (MBA-1348). Required together with
    /// `turret_windage_dialed_mil` — `optic_profile()` returns a named-field `Err` if
    /// only one axis of the pair is set, rather than silently assuming the other reads
    /// zero.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub turret_elevation_dialed_mil: Option<f64>,
    /// The windage turret's current dialed offset from zero, DIAL-space mil, signed:
    /// positive is dialed RIGHT (MBA-1348). See `turret_elevation_dialed_mil`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub turret_windage_dialed_mil: Option<f64>,
    /// The reticle's usable hold extent ABOVE center, TRUE angular mil (MBA-1348) — see
    /// `HoldBounds`, an explicit spec input (manufacturer spec sheet or bench
    /// measurement), never derived from this profile's own `reticle` field. Required
    /// together with `hold_bound_down_mil`/`hold_bound_left_mil`/`hold_bound_right_mil`
    /// — all four or none.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub hold_bound_up_mil: Option<f64>,
    /// The reticle's usable hold extent BELOW center, TRUE angular mil (MBA-1348). See
    /// `hold_bound_up_mil`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub hold_bound_down_mil: Option<f64>,
    /// The reticle's usable hold extent LEFT of center, TRUE angular mil (MBA-1348). See
    /// `hold_bound_up_mil`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub hold_bound_left_mil: Option<f64>,
    /// The reticle's usable hold extent RIGHT of center, TRUE angular mil (MBA-1348). See
    /// `hold_bound_up_mil`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub hold_bound_right_mil: Option<f64>,
}

/// One named zero condition / per-load dial correction (MBA-1360).
///
/// `zero_distance` uses the SAME display-unit convention as [`ProfileData::zero_distance`]
/// (yards imperial / meters metric per the profile's `units` field), and
/// `ProfileData::converted_to` (in `main.rs`) rescales it identically. `poi_up_mil`/`poi_right_mil`
/// are constant ANGULAR dial corrections in MILs (unit-invariant, untouched by
/// `converted_to`), ADDED to the dial outputs (elevation/windage adjustments) when the
/// set is selected — positive = dial UP/RIGHT more; a load that impacts high/right
/// relative to the master zero therefore stores negative values. This is deliberately
/// dial-side (a constant angular correction per the ticket), unlike the MBA-1359
/// linear-at-zero-range POI offsets, which bias the solved zero itself; the two compose.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct ProfileZeroSet {
    pub name: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub zero_distance: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub poi_up_mil: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub poi_right_mil: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub notes: Option<String>,
}

/// One velocity-banded BC breakpoint (profile schema v2, MBA-1323 Phase 2). Stored as a raw
/// breakpoint, NOT a pre-computed `velocity_min`/`velocity_max` band — banding into the
/// engine's [`crate::BCSegmentData`] shape happens at solve time (`bc_segments_from_profile`), so the
/// stored JSON stays a simple, order-independent list that round-trips cleanly.
///
/// `velocity_mps` is ALWAYS meters/second, independent of [`ProfileData::units`]. This is
/// intentional (matches the engine's internal BC-segment plumbing, which is also always in
/// engine units) and is why `ProfileData::converted_to` (in `main.rs`) leaves this field untouched — see
/// the comment there.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct ProfileBcSegment {
    pub bc: f64,
    pub velocity_mps: f64,
}

/// One Mach/Cd point of a saved custom drag curve (profile schema v2, MBA-1323 Phase 2).
/// Both fields are unit-invariant (Mach is dimensionless; Cd is dimensionless), so
/// `ProfileData::converted_to` (in `main.rs`) leaves `drag_curve` untouched too.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct ProfileDragPoint {
    pub mach: f64,
    pub cd: f64,
}

fn default_unit_system() -> String {
    "imperial".to_string()
}
fn default_temperature() -> f64 {
    59.0
}
fn default_pressure() -> f64 {
    29.92
}
fn default_humidity() -> f64 {
    50.0
}

/// Assembles an all-or-nothing pair of angular profile fields — `OpticProfile`'s
/// `TravelLimits`/`TurretState` axes — from two independently-`Option` values (MBA-1348).
/// `None` only when NEITHER is set; `Some` only when BOTH are; a named-field `Err` when
/// exactly one is. A silently-defaulted `0.0` for the missing half would assert a
/// specific, likely-false physical fact (e.g. "zero down travel from zero") rather than
/// "not recorded" — exactly the silent fabrication `OpticProfile`'s own doc comments
/// warn against — so an incomplete pair is rejected by name instead of guessed at.
///
/// `a_flag`/`b_flag` name the CLI FLAGS (e.g. `"--travel-up"`), not the internal
/// `ProfileData` struct fields (MBA-1348 review fix I2): the field↔flag mapping is not
/// mechanical (`--travel-up` drops the axis prefix `elevation_travel_up_mil` carries;
/// `--turret-elev` abbreviates `turret_elevation_dialed_mil`), and a user who typed
/// `--travel-up` has no reason to recognize `elevation_travel_up_mil` in an error. This
/// matches the convention two call sites above already use
/// (`"--elevation-click '{v}' is invalid"`).
pub fn require_angular_pair(
    a_flag: &'static str,
    a: Option<f64>,
    b_flag: &'static str,
    b: Option<f64>,
) -> Result<Option<(f64, f64)>, String> {
    match (a, b) {
        (None, None) => Ok(None),
        (Some(a), Some(b)) => Ok(Some((a, b))),
        (Some(_), None) => Err(format!(
            "{a_flag} is set but {b_flag} is not — both are required together (or neither)"
        )),
        (None, Some(_)) => Err(format!(
            "{b_flag} is set but {a_flag} is not — both are required together (or neither)"
        )),
    }
}

/// The `require_angular_pair` all-or-nothing rule extended to `HoldBounds`' four fields
/// (MBA-1348): all four or none, any other combination is a named-field `Err` (same
/// no-silent-fabrication rationale — see `require_angular_pair`). Names the four CLI
/// FLAGS in any error, not the internal struct fields (MBA-1348 review fix I2) — see
/// `require_angular_pair`'s own doc comment for why that distinction matters.
pub fn require_hold_bounds(
    up: Option<f64>,
    down: Option<f64>,
    left: Option<f64>,
    right: Option<f64>,
) -> Result<Option<HoldBounds>, String> {
    let flags = [
        ("--hold-up", up),
        ("--hold-down", down),
        ("--hold-left", left),
        ("--hold-right", right),
    ];
    let set_count = flags.iter().filter(|(_, v)| v.is_some()).count();
    if set_count == 0 {
        return Ok(None);
    }
    if set_count < 4 {
        let missing: Vec<&str> = flags
            .iter()
            .filter(|(_, v)| v.is_none())
            .map(|(name, _)| *name)
            .collect();
        return Err(format!(
            "reticle hold bounds are incomplete — missing {} (all four of --hold-up/\
             --hold-down/--hold-left/--hold-right are required together, or none)",
            missing.join(", ")
        ));
    }
    Ok(Some(HoldBounds {
        up_mil: up.unwrap(),
        down_mil: down.unwrap(),
        left_mil: left.unwrap(),
        right_mil: right.unwrap(),
    }))
}

/// MBA-1358: a scope tracking correction factor must be strictly between 0.5 and 1.5.
/// Values outside that band mean the tall-target test (or a hand edit) went wrong — a
/// real scope does not mistrack by 50% — so they are a hard error naming the offending
/// field/flag rather than a silently applied scale. The band itself lives in the shared
/// `crate::adjustment::tracking_cf_in_range` so the WASM terminal enforces
/// the identical bound.
pub fn validate_tracking_cf(value: f64, source: &str) -> Result<(), String> {
    if crate::adjustment::tracking_cf_in_range(value) {
        Ok(())
    } else {
        Err(format!(
            "{source} must be a tracking correction factor strictly between 0.5 and 1.5 \
             (got {value}); derive it from a tall-target test with `ballistics tall-target`"
        ))
    }
}

impl ProfileData {
    /// Assembles this profile's `OpticProfile` (MBA-1348) from `elevation_click`/
    /// `windage_click` (parsed via `parse_click_value`, windage falling back to the
    /// resolved elevation graduation — the same precedence `resolve_click_values` uses)
    /// plus the twelve turret/hold fields declared alongside `reticle` above.
    ///
    /// Returns `Ok(None)` only when NONE of the twelve fields are set: the profile has
    /// never been given any turret model at all. Every other combination either succeeds
    /// (`Ok(Some(..))`) or is a named-field `Err` — including `elevation_click` itself
    /// being unset while ANY of the other eleven fields IS set, since those eleven are
    /// meaningless without a click graduation (`OpticProfile` cannot be constructed
    /// without one, structurally) and silently ignoring them would let a save persist
    /// turret/hold data no downstream code could ever use.
    ///
    /// Does NOT call `OpticProfile::validate()` itself — callers that need a
    /// pre-validated profile (`profile save`) call it explicitly, so a validation
    /// failure is attributed to the operation asking for it rather than to assembly.
    pub fn optic_profile(&self) -> Result<Option<OpticProfile>, String> {
        // Shape/pairing checks run FIRST and UNCONDITIONALLY -- before the elevation_click
        // gate below -- so an incomplete travel/turret-state/hold-bound pair is rejected by
        // name even on a profile that never set elevation_click at all (MBA-1348 review
        // fix: a save with, e.g., only --travel-up and no --elevation-click previously
        // skipped every one of these checks via the early `Ok(None)` return, silently
        // persisting an incomplete pair no downstream code could ever have used anyway).
        let elevation_travel = require_angular_pair(
            "--travel-up",
            self.elevation_travel_up_mil,
            "--travel-down",
            self.elevation_travel_down_mil,
        )?
        .map(|(up, down)| TravelLimits { up_mil: up, down_mil: down });

        let windage_travel = require_angular_pair(
            "--windage-travel-left",
            self.windage_travel_left_mil,
            "--windage-travel-right",
            self.windage_travel_right_mil,
        )?
        .map(|(left, right)| TravelLimits { down_mil: left, up_mil: right });

        let turret_state = require_angular_pair(
            "--turret-elev",
            self.turret_elevation_dialed_mil,
            "--turret-wind",
            self.turret_windage_dialed_mil,
        )?
        .map(|(elevation_mil, windage_mil)| TurretState { elevation_mil, windage_mil });

        let reticle_hold_bounds = require_hold_bounds(
            self.hold_bound_up_mil,
            self.hold_bound_down_mil,
            self.hold_bound_left_mil,
            self.hold_bound_right_mil,
        )?;

        let Some(elev_str) = &self.elevation_click else {
            // MBA-1348 review fix: every one of the other eleven fields is meaningless
            // without a click graduation -- `OpticProfile` cannot even be constructed
            // without one, structurally -- so silently accepting them here (the way a
            // bare `Ok(None)` would) lets a save persist turret/hold data that no
            // downstream code, including this profile's own future re-saves, could ever
            // use. A profile that has genuinely never touched any of the twelve fields
            // still returns Ok(None), unchanged.
            if self.clicks_per_revolution.is_some()
                || self.zero_stop.is_some()
                || elevation_travel.is_some()
                || windage_travel.is_some()
                || turret_state.is_some()
                || reticle_hold_bounds.is_some()
            {
                return Err(
                    "turret/hold fields are set but --elevation-click is not -- a turret \
                     model needs a click graduation to be usable at all; supply \
                     --elevation-click, or also pass --clear-turret to clear the other \
                     fields (e.g. after --clear-click removed the click alone)"
                        .to_string(),
                );
            }
            return Ok(None);
        };
        let elevation_click = parse_click_value(elev_str)?;
        let windage_click = match &self.windage_click {
            Some(s) => parse_click_value(s)?,
            None => elevation_click,
        };

        Ok(Some(OpticProfile {
            elevation_click,
            windage_click,
            clicks_per_revolution: self.clicks_per_revolution,
            zero_stop: self.zero_stop.unwrap_or(false),
            elevation_travel,
            windage_travel,
            turret_state,
            reticle_hold_bounds,
        }))
    }

    /// Parse this profile's `units` field into the CLI unit system it names.
    ///
    /// Factored out of the CLI's `converted_to` (which now calls this) so the bridge's
    /// `profile.validate` applies the identical check with the identical message.
    pub fn unit_system(&self) -> Result<UnitSystem, String> {
        match self.units.trim().to_ascii_lowercase().as_str() {
            "imperial" => Ok(UnitSystem::Imperial),
            "metric" => Ok(UnitSystem::Metric),
            other => Err(format!(
                "Profile '{}' has unsupported units '{}'; expected 'imperial' or 'metric'",
                self.name, other
            )),
        }
    }

    /// The cheap invariants the CLI already applies to a saved profile, collected instead
    /// of short-circuited: the `units` string (`converted_to`'s gate), the MBA-1358
    /// tracking-CF band (`load_profile`'s gate), and the MBA-1348 turret/optic assembly +
    /// validation (`profile save`'s gate, including `parse_click_value` on the stored
    /// click graduations). No new physics checks — this is exactly the existing load/save
    /// surface, aggregated for the bridge's `profile.validate`. Empty means the profile
    /// passes every one of those gates.
    pub fn validation_warnings(&self) -> Vec<String> {
        let mut warnings = Vec::new();
        if let Err(err) = self.unit_system() {
            warnings.push(err);
        }
        if let Some(cf) = self.elevation_cf {
            if let Err(err) = validate_tracking_cf(cf, "profile field elevation_cf") {
                warnings.push(err);
            }
        }
        if let Some(cf) = self.windage_cf {
            if let Err(err) = validate_tracking_cf(cf, "profile field windage_cf") {
                warnings.push(err);
            }
        }
        match self.optic_profile() {
            Err(err) => warnings.push(format!("turret/optic model: {err}")),
            Ok(Some(optic)) => {
                if let Err(err) = optic.validate() {
                    warnings.push(format!("turret/optic model: {err}"));
                }
            }
            Ok(None) => {}
        }
        warnings
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A literal on-disk profile document exercising the big wire surface: the five
    /// mandatory fields, bc_segments, drag_curve, dsf_points, tracking CFs, named zero
    /// sets with POI corrections, the click graduations plus all twelve turret/hold
    /// fields, and a reticle — PLUS a top-level key this build has never heard of
    /// ("future_field") and an unknown key nested inside a zero set, because tolerating
    /// unknown keys (no `deny_unknown_fields`) is the documented forward-compat behavior
    /// this struct has always had and the move into the library must not change.
    const FIXTURE: &str = r#"{
        "name": "wire-shape-fixture",
        "velocity": 762.0,
        "bc": 0.243,
        "mass": 11.33980925,
        "diameter": 7.8232,
        "drag_model": "G7",
        "twist_rate": 254.0,
        "sight_height": 50.8,
        "zero_distance": 91.44,
        "units": "metric",
        "temperature": 15.0,
        "pressure": 1013.207888,
        "humidity": 55.0,
        "altitude": 304.8,
        "bullet_name": "175gr SMK",
        "created": "1755400000",
        "wind_speed": 4.4704,
        "wind_direction": 90.0,
        "shooting_angle": -5.0,
        "auto_zero": 91.44,
        "twist_right": false,
        "use_bc_segments": true,
        "bullet_length": 30.48,
        "elevation_click": "0.1mil",
        "windage_click": "0.25moa",
        "bc_segments": [
            {"bc": 0.243, "velocity_mps": 792.0},
            {"bc": 0.230, "velocity_mps": 400.0}
        ],
        "drag_curve": [
            {"mach": 0.5, "cd": 0.23},
            {"mach": 1.2, "cd": 0.45},
            {"mach": 3.0, "cd": 0.28}
        ],
        "dsf_points": [
            {"mach": 0.9, "dsf": 1.04}
        ],
        "bc_reference": "army-standard-metro",
        "pressure_reference": "qnh",
        "density_altitude": 500.0,
        "zero_poi_up_m": 0.05,
        "zero_poi_right_m": -0.02,
        "sight_offset_lateral_m": 0.01,
        "elevation_cf": 0.97,
        "windage_cf": 1.02,
        "zero_sets": [
            {"name": "suppressed", "zero_distance": 200.0, "poi_up_mil": -0.3,
             "poi_right_mil": 0.1, "notes": "suppressed load",
             "zero_set_future_field": true},
            {"name": "match", "poi_up_mil": 0.25}
        ],
        "reticle": {
            "name": "mil-grid 0.5/10",
            "focal_plane": "ffp",
            "reference_magnification": 10.0,
            "marks": [
                {"down_mil": 0.0, "right_mil": 0.0, "kind": "center"},
                {"down_mil": 1.0, "right_mil": 0.0, "kind": "hash", "label": "1.0"}
            ]
        },
        "clicks_per_revolution": 100,
        "zero_stop": true,
        "elevation_travel_up_mil": 26.0,
        "elevation_travel_down_mil": 4.0,
        "windage_travel_left_mil": 12.0,
        "windage_travel_right_mil": 12.0,
        "turret_elevation_dialed_mil": 5.4,
        "turret_windage_dialed_mil": -0.2,
        "hold_bound_up_mil": 5.0,
        "hold_bound_down_mil": 10.0,
        "hold_bound_left_mil": 6.0,
        "hold_bound_right_mil": 6.0,
        "future_field": {"nested": [1, 2, 3]}
    }"#;

    /// deserialize -> serialize -> deserialize is the identity on every stored field,
    /// and the two serializations are byte-identical to each other (the engine's own
    /// output is stable). Unknown input keys are tolerated on the way in — the
    /// pre-existing behavior — and are NOT preserved on the way out (serde drops them;
    /// also pre-existing, documented in the module doc).
    #[test]
    fn wire_shape_round_trips_field_for_field() {
        let first: ProfileData = serde_json::from_str(FIXTURE).expect("fixture must load");
        let serialized = serde_json::to_string(&first).expect("serialize");
        let second: ProfileData = serde_json::from_str(&serialized).expect("reload");

        // The engine re-serializes its own output identically.
        assert_eq!(
            serde_json::to_value(&first).unwrap(),
            serde_json::to_value(&second).unwrap()
        );

        // Field-level equality across the whole surface (ProfileData does not derive
        // PartialEq, deliberately — nothing in production compares whole profiles).
        assert_eq!(second.name, "wire-shape-fixture");
        assert_eq!(second.velocity.to_bits(), 762.0f64.to_bits());
        assert_eq!(second.bc.to_bits(), 0.243f64.to_bits());
        assert_eq!(second.mass.to_bits(), 11.33980925f64.to_bits());
        assert_eq!(second.diameter.to_bits(), 7.8232f64.to_bits());
        assert_eq!(second.drag_model, "G7");
        assert_eq!(second.twist_rate, Some(254.0));
        assert_eq!(second.sight_height, Some(50.8));
        assert_eq!(second.zero_distance, Some(91.44));
        assert_eq!(second.units, "metric");
        assert_eq!(second.temperature.to_bits(), 15.0f64.to_bits());
        assert_eq!(second.pressure.to_bits(), 1013.207888f64.to_bits());
        assert_eq!(second.humidity.to_bits(), 55.0f64.to_bits());
        assert_eq!(second.altitude.to_bits(), 304.8f64.to_bits());
        assert_eq!(second.bullet_name.as_deref(), Some("175gr SMK"));
        assert_eq!(second.created.as_deref(), Some("1755400000"));
        assert_eq!(second.wind_speed, Some(4.4704));
        assert_eq!(second.wind_direction, Some(90.0));
        assert_eq!(second.shooting_angle, Some(-5.0));
        assert_eq!(second.auto_zero, Some(91.44));
        assert_eq!(second.twist_right, Some(false));
        assert_eq!(second.use_bc_segments, Some(true));
        assert_eq!(second.bullet_length, Some(30.48));
        assert_eq!(second.elevation_click.as_deref(), Some("0.1mil"));
        assert_eq!(second.windage_click.as_deref(), Some("0.25moa"));
        assert_eq!(
            second.bc_segments,
            Some(vec![
                ProfileBcSegment { bc: 0.243, velocity_mps: 792.0 },
                ProfileBcSegment { bc: 0.230, velocity_mps: 400.0 },
            ])
        );
        assert_eq!(
            second.drag_curve,
            Some(vec![
                ProfileDragPoint { mach: 0.5, cd: 0.23 },
                ProfileDragPoint { mach: 1.2, cd: 0.45 },
                ProfileDragPoint { mach: 3.0, cd: 0.28 },
            ])
        );
        assert_eq!(second.dsf_points, Some(vec![DsfPoint { mach: 0.9, dsf: 1.04 }]));
        assert_eq!(second.bc_reference.as_deref(), Some("army-standard-metro"));
        assert_eq!(second.pressure_reference.as_deref(), Some("qnh"));
        assert_eq!(second.density_altitude, Some(500.0));
        assert_eq!(second.zero_poi_up_m, Some(0.05));
        assert_eq!(second.zero_poi_right_m, Some(-0.02));
        assert_eq!(second.sight_offset_lateral_m, Some(0.01));
        assert_eq!(second.elevation_cf, Some(0.97));
        assert_eq!(second.windage_cf, Some(1.02));
        let sets = second.zero_sets.as_deref().expect("zero_sets kept");
        assert_eq!(sets.len(), 2);
        assert_eq!(sets[0].name, "suppressed");
        assert_eq!(sets[0].zero_distance, Some(200.0));
        assert_eq!(sets[0].poi_up_mil, Some(-0.3));
        assert_eq!(sets[0].poi_right_mil, Some(0.1));
        assert_eq!(sets[0].notes.as_deref(), Some("suppressed load"));
        assert_eq!(sets[1].name, "match");
        assert_eq!(sets[1].zero_distance, None);
        assert_eq!(sets[1].poi_up_mil, Some(0.25));
        assert_eq!(sets[1].poi_right_mil, None);
        let reticle = second.reticle.as_ref().expect("reticle kept");
        assert_eq!(reticle.name, "mil-grid 0.5/10");
        assert_eq!(reticle.marks.len(), 2);
        assert_eq!(reticle.marks[1].label.as_deref(), Some("1.0"));
        assert_eq!(second.clicks_per_revolution, Some(100));
        assert_eq!(second.zero_stop, Some(true));
        assert_eq!(second.elevation_travel_up_mil, Some(26.0));
        assert_eq!(second.elevation_travel_down_mil, Some(4.0));
        assert_eq!(second.windage_travel_left_mil, Some(12.0));
        assert_eq!(second.windage_travel_right_mil, Some(12.0));
        assert_eq!(second.turret_elevation_dialed_mil, Some(5.4));
        assert_eq!(second.turret_windage_dialed_mil, Some(-0.2));
        assert_eq!(second.hold_bound_up_mil, Some(5.0));
        assert_eq!(second.hold_bound_down_mil, Some(10.0));
        assert_eq!(second.hold_bound_left_mil, Some(6.0));
        assert_eq!(second.hold_bound_right_mil, Some(6.0));

        // Unknown keys were tolerated (this test got this far) and dropped on re-save
        // (pre-existing serde behavior, documented in the module doc): no invented keys.
        let reserialized: serde_json::Value = serde_json::from_str(&serialized).unwrap();
        assert!(reserialized.get("future_field").is_none());

        // This fixture also passes every load/save invariant `validation_warnings`
        // aggregates (units, CF band, optic assembly + validation).
        assert_eq!(second.validation_warnings(), Vec::<String>::new());
    }

    /// The minimal legacy (pre-v2) document: only the five mandatory fields. The four
    /// defaulted scalars fill in (imperial / 59 F / 29.92 inHg / 50%), every `Option`
    /// stays `None`, and — load-bearing for on-disk compatibility — re-serialization
    /// keeps the non-`skip_serializing_if` optional keys (`twist_rate`, `sight_height`,
    /// `zero_distance`, `bullet_name`, `created`) present as JSON `null` while every
    /// `skip_serializing_if` key stays absent.
    #[test]
    fn legacy_minimal_document_defaults_and_reserializes_with_stable_keys() {
        let legacy = r#"{
            "name": "legacy",
            "velocity": 2500.0,
            "bc": 0.475,
            "mass": 175.0,
            "diameter": 0.308,
            "drag_model": "G1"
        }"#;
        let profile: ProfileData = serde_json::from_str(legacy).expect("legacy must load");
        assert_eq!(profile.units, "imperial");
        assert_eq!(profile.temperature.to_bits(), 59.0f64.to_bits());
        assert_eq!(profile.pressure.to_bits(), 29.92f64.to_bits());
        assert_eq!(profile.humidity.to_bits(), 50.0f64.to_bits());
        assert_eq!(profile.altitude.to_bits(), 0.0f64.to_bits());
        assert!(profile.bc_segments.is_none());
        assert!(profile.zero_sets.is_none());

        let out: serde_json::Value =
            serde_json::from_str(&serde_json::to_string(&profile).unwrap()).unwrap();
        // Plain-`default` options serialize as explicit null (the historical shape)...
        assert!(out.get("twist_rate").is_some_and(serde_json::Value::is_null));
        assert!(out.get("bullet_name").is_some_and(serde_json::Value::is_null));
        assert!(out.get("created").is_some_and(serde_json::Value::is_null));
        // ...while `skip_serializing_if` options stay absent entirely.
        for absent in [
            "wind_speed", "bc_segments", "drag_curve", "dsf_points", "zero_sets",
            "reticle", "elevation_click", "elevation_cf", "clicks_per_revolution",
        ] {
            assert!(out.get(absent).is_none(), "{absent} must stay absent");
        }
        assert_eq!(profile.validation_warnings(), Vec::<String>::new());
    }

    /// `validation_warnings` reports exactly the CLI's own gates: an out-of-band
    /// tracking CF (load gate), a bad units string (`converted_to` gate), and an
    /// incomplete turret pair / unparseable click (save gate) — and nothing for a
    /// profile that merely leaves everything optional unset.
    #[test]
    fn validation_warnings_mirror_the_cli_gates() {
        let mut profile: ProfileData = serde_json::from_str(FIXTURE).unwrap();
        profile.elevation_cf = Some(0.4); // outside the (0.5, 1.5) band
        profile.units = "nautical".to_string();
        profile.elevation_travel_down_mil = None; // breaks the all-or-nothing pair
        let warnings = profile.validation_warnings();
        assert_eq!(warnings.len(), 3, "{warnings:?}");
        assert!(warnings.iter().any(|w| w.contains("elevation_cf")), "{warnings:?}");
        assert!(warnings.iter().any(|w| w.contains("unsupported units")), "{warnings:?}");
        assert!(warnings.iter().any(|w| w.contains("--travel-down")), "{warnings:?}");
    }
}
