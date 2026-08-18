//! `.a7p` -> [`ProfileData`] mapping and the import report (moved verbatim from `main.rs`
//! so the CLI's `profile import` and the bridge's `profile.import_a7p` share ONE mapping).
//!
//! The CLI keeps the presentation layer (`render_import_report`, the `--strict` refusal,
//! `--name` sanitization notes, saving); everything that decides WHAT a `.a7p` becomes —
//! field conversions, the honest unmapped/warning accounting, `--zero-click` click-count
//! conversion — lives here, behind the same `profile-import` feature as the parser.

use crate::adjustment::ClickValue;
use crate::profile::{ProfileBcSegment, ProfileData, ProfileDragPoint, ProfileZeroSet};

/// Get a timestamp string without chrono (same implementation as the CLI's own
/// `timestamp_string`; duplicated rather than exported because a seconds-since-epoch
/// formatter is not engine API).
fn timestamp_string() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();
    format!("{}", secs)
}

/// Everything `profile import` produced: the profile to save plus the honest
/// account of what mapped, what did not, and why.
#[derive(Debug)]
pub struct ImportReport {
    /// (source field, raw value, converted value, destination field)
    pub mapped: Vec<[String; 4]>,
    /// (source field, human explanation) — data the profile store cannot hold.
    pub unmapped: Vec<(String, String)>,
    pub warnings: Vec<String>,
}

#[derive(Debug)]
pub struct A7pImportOutcome {
    pub profile: ProfileData,
    pub report: ImportReport,
}

// Re-export of the shared constant under the name this module's a7p-mapping
// call sites already use (MBA-1327: single source of truth for grain<->gram).
use crate::constants::GRAMS_PER_GRAIN as GRAIN_TO_GRAM;
const IN_TO_MM: f64 = 25.4;

/// Restrict imported profile names to characters that are safe as file names
/// in the profile store (`~/.ballistics/profiles/<name>.json`).
pub fn sanitize_profile_name(raw: &str) -> String {
    let cleaned: String = raw
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || matches!(c, ' ' | '.' | '_' | '-') {
                c
            } else {
                '_'
            }
        })
        .collect();
    let trimmed = cleaned.trim().to_string();
    if trimmed.is_empty() {
        "imported-a7p".to_string()
    } else {
        trimmed
    }
}

/// `ClickValue`'s canonical suffixed profile-field string (MBA-1348), e.g. `"0.1mil"` —
/// exactly what `parse_click_value` parses back, produced via the type's own
/// `Serialize` impl rather than a second, hand-rolled suffix match that could drift
/// from it.
fn click_value_to_profile_string(click: ClickValue) -> String {
    serde_json::to_string(&click)
        .expect("ClickValue serialization is infallible")
        .trim_matches('"')
        .to_string()
}

pub fn map_a7p_to_profile(
    doc: &super::A7pDocument,
    name_override: Option<&str>,
    // MBA-1359: the source device's scope click graduation (`--zero-click`). The .a7p
    // format stores zeroing state as raw device click counts WITHOUT the click size, so
    // this is the only way to convert `zero_x`/`zero_y` into a linear POI offset. `None`
    // keeps the historical behavior (the click counts are reported as unmapped).
    zero_click: Option<ClickValue>,
) -> Result<A7pImportOutcome, String> {
    use super::{A7pBcType, EnvelopeStatus};
    let src = &doc.profile;

    let mut report = ImportReport {
        mapped: Vec::new(),
        unmapped: Vec::new(),
        warnings: Vec::new(),
    };
    if let EnvelopeStatus::Mismatch { expected, actual } = &doc.envelope {
        report.warnings.push(format!(
            "checksum mismatch (file says {expected}, payload hashes to {actual}) — file may be corrupted"
        ));
    }

    let name = match name_override {
        Some(n) => n.to_string(),
        None => sanitize_profile_name(&src.profile_name),
    };

    let mut push = |field: &str, raw: String, converted: String, dest: &str| {
        report
            .mapped
            .push([field.to_string(), raw, converted, dest.to_string()]);
    };
    push(
        "profile_name",
        src.profile_name.clone(),
        name.clone(),
        "name",
    );
    push(
        "c_muzzle_velocity",
        format!("{:.1} m/s", src.muzzle_velocity_mps),
        format!("{:.1} m/s", src.muzzle_velocity_mps),
        "velocity (muzzle velocity)",
    );

    // Resolve drag model + BC-related profile fields. Branches by bc_type because G1/G7 and
    // CUSTOM disagree on what `coef_rows_raw` even means (velocity-BC rows vs Mach-Cd rows —
    // see A7pProfile::bc_rows()/custom_rows()) and on what the scalar `bc` field should hold.
    let (drag_model, bc, bc_segments, drag_curve): (
        &str,
        f64,
        Option<Vec<ProfileBcSegment>>,
        Option<Vec<ProfileDragPoint>>,
    ) = match src.bc_type {
        A7pBcType::G1 | A7pBcType::G7 => {
            let drag_model = if matches!(src.bc_type, A7pBcType::G1) {
                "G1"
            } else {
                "G7"
            };
            let rows = src.bc_rows();
            // The row measured at the highest velocity is the muzzle-regime BC, retained as
            // the scalar `bc` for back-compat with tools that only understand one BC.
            let (bc, bc_row_velocity) = rows.iter().copied().max_by(|a, b| a.1.total_cmp(&b.1)).ok_or_else(
                || "no BC rows in file — cannot build a profile without a BC".to_string(),
            )?;
            push(
                "coef_rows[fastest]",
                format!("BC {bc:.3} @ {bc_row_velocity:.0} m/s"),
                format!("{bc:.3} ({drag_model})"),
                "bc + drag_model",
            );

            let bc_segments = if rows.len() > 1 {
                // Descending by velocity: matches bc_segments_from_profile's/the engine's
                // "fastest row governs the muzzle regime" convention, and puts the back-compat
                // scalar `bc` above (= sorted[0].bc) in visible agreement with this list.
                let mut sorted = rows.clone();
                sorted.sort_by(|a, b| b.1.total_cmp(&a.1));
                push(
                    "coef_rows[all]",
                    format!("{} row(s), fastest {bc:.3} @ {bc_row_velocity:.0} m/s", rows.len()),
                    format!("{} bc_segments (velocity-banded, descending)", sorted.len()),
                    "bc_segments",
                );
                Some(
                    sorted
                        .into_iter()
                        .map(|(bc, velocity_mps)| ProfileBcSegment { bc, velocity_mps })
                        .collect(),
                )
            } else {
                None
            };
            (drag_model, bc, bc_segments, None)
        }
        A7pBcType::Custom => {
            let mut pairs = src.custom_rows(); // (Cd, Mach), file order
            pairs.sort_by(|a, b| a.1.total_cmp(&b.1)); // ascending by Mach (DragTable requirement)
            let mach_values: Vec<f64> = pairs.iter().map(|&(_, mach)| mach).collect();
            let cd_values: Vec<f64> = pairs.iter().map(|&(cd, _)| cd).collect();
            // Validate now (at import time) rather than only at first solve, so a malformed
            // curve is a clear import-time error instead of a confusing later failure.
            crate::drag::DragTable::try_new(mach_values.clone(), cd_values.clone())
                .map_err(|e| format!("CUSTOM drag curve is invalid: {e}"))?;
            push(
                "coef_rows[custom]",
                format!("{} Cd/Mach row(s)", pairs.len()),
                format!(
                    "{} drag_curve point(s), Mach {:.3}-{:.3}",
                    pairs.len(),
                    mach_values.first().copied().unwrap_or(0.0),
                    mach_values.last().copied().unwrap_or(0.0)
                ),
                "drag_curve",
            );
            // No scalar BC applies to a full drag curve — see map_a7p_to_profile's module-level
            // rationale in the report below. `bc: 0.0` is an intentionally-invalid sentinel:
            //   * it is physically inert once drag_curve is consumed (custom_drag_table divides
            //     by sectional density, not `bc_value` — BallisticInputs::custom_drag_denominator);
            //   * commands that do not YET consume drag_curve (see CLI_USAGE.md's a7p import
            //     section) will fail loudly (`bc_value must be finite and greater than zero`)
            //     instead of silently running the wrong physics under an assumed G1 model. That
            //     loud failure is the honest outcome for an unwired path, not a bug to paper over.
            report.warnings.push(
                "bc_type CUSTOM: no scalar BC applies to a full drag curve; the profile's 'bc' \
                 field is set to 0.0 as an inert sentinel. It is unused once drag_curve is \
                 consumed (a custom drag table replaces the BC-based retardation model \
                 entirely); commands that do not yet consume drag_curve will fail loudly \
                 (bc_value must be > 0) rather than silently assuming a G1 model — see \
                 CLI_USAGE.md."
                    .to_string(),
            );
            let drag_curve = Some(
                pairs
                    .into_iter()
                    .map(|(cd, mach)| ProfileDragPoint { mach, cd })
                    .collect(),
            );
            ("CUSTOM", 0.0, None, drag_curve)
        }
        A7pBcType::Other(v) => {
            return Err(format!("unknown bc_type {v} — file newer than this importer"))
        }
    };

    push(
        "b_weight",
        format!("{:.1} gr", src.bullet_weight_gr),
        format!("{:.4} g", src.bullet_weight_gr * GRAIN_TO_GRAM),
        "mass",
    );
    push(
        "b_diameter",
        format!("{:.3} in", src.bullet_diameter_in),
        format!("{:.3} mm", src.bullet_diameter_in * IN_TO_MM),
        "diameter",
    );
    push(
        "b_length",
        format!("{:.3} in", src.bullet_length_in),
        format!("{:.2} mm", src.bullet_length_in * IN_TO_MM),
        "bullet_length",
    );
    push(
        "r_twist / twist_dir",
        format!(
            "{:.2} in/turn, {}",
            src.twist_in_per_turn,
            if src.twist_right { "RIGHT" } else { "LEFT" }
        ),
        format!("{:.1} mm/turn", src.twist_in_per_turn * IN_TO_MM),
        "twist_rate + twist_right",
    );
    push(
        "sc_height",
        format!("{:.0} mm", src.sight_height_mm),
        format!("{:.0} mm", src.sight_height_mm),
        "sight_height",
    );
    if let Some(zd) = src.zero_distance_m {
        push(
            "distances[c_zero_distance_idx]",
            format!("{zd:.1} m"),
            format!("{zd:.1} m"),
            "zero_distance + auto_zero",
        );
    }
    push(
        "c_zero_air_temperature",
        format!("{:.1} C", src.air_temperature_c),
        format!("{:.1} C", src.air_temperature_c),
        "temperature",
    );
    push(
        "c_zero_air_pressure",
        format!("{:.1} hPa", src.air_pressure_hpa),
        format!("{:.1} hPa", src.air_pressure_hpa),
        "pressure",
    );
    push(
        "c_zero_air_humidity",
        format!("{:.0} %", src.air_humidity_pct),
        format!("{:.0} %", src.air_humidity_pct),
        "humidity",
    );
    if !src.bullet_name.is_empty() {
        push(
            "bullet_name",
            src.bullet_name.clone(),
            src.bullet_name.clone(),
            "bullet_name",
        );
    }

    // Honest non-mapping: things the profile store cannot hold today.
    let tcoeff_mps_per_c =
        src.muzzle_velocity_mps * (src.temp_coeff_pct_per_15c / 100.0) / 15.0;
    report.unmapped.push((
        "c_t_coeff".to_string(),
        format!(
            "{:.3} %/15C = {:.3} m/s per C powder sensitivity — profile schema does not model \
             powder sensitivity",
            src.temp_coeff_pct_per_15c, tcoeff_mps_per_c
        ),
    ));
    report.unmapped.push((
        "c_zero_p_temperature".to_string(),
        format!("{:.0} C powder temperature at zeroing", src.powder_temperature_c),
    ));
    report.unmapped.push((
        "c_zero_temperature".to_string(),
        format!(
            "{:.0} C ambient temperature when the zero was established",
            src.zero_temperature_c
        ),
    ));
    if src.w_pitch_raw != 0 {
        report.unmapped.push((
            "c_zero_w_pitch".to_string(),
            format!(
                "zeroing pitch value ({}) — base pitch is rifle-mount geometry outside \
                 the turret model (OpticProfile has no base-pitch concept); its effect \
                 is already absorbed by zeroing, not a turret setting left to store",
                src.w_pitch_raw
            ),
        ));
    }
    // MBA-1359: `zero_x`/`zero_y` are the device's zeroing state in CLICK counts x 1000
    // (upstream a7p spec: "zeroing h-clicks / v-clicks for specific device"). The file does
    // NOT carry the device's click size, so conversion is only possible when the user
    // supplies it (`--zero-click`). Axis conventions, confirmed against upstream tooling
    // (a7p's own CLI negates X on entry: `zero_x += round(x_offset * -1000)`,
    // `zero_y += round(y_offset * 1000)`): user-facing right-offset clicks = -zero_x/1000,
    // up-offset clicks = zero_y/1000. Linear offset at the zero range =
    // clicks x (click size / adjustment_factor(base)) [radians] x zero distance [m].
    // MBA-1348: whenever the caller supplies the device's click graduation
    // (--zero-click), that is itself a modeled fact regardless of whether the file's
    // zero_x/zero_y counts can ALSO be converted to a POI offset below (which
    // additionally needs a zero distance) -- so it is recorded as the profile's turret
    // click graduation independently of that offset conversion. "only when not already
    // set" is always satisfied here since map_a7p_to_profile always builds a fresh
    // ProfileData (elevation_click starts unset), stated for the same defensive reason
    // resolve_click_values documents its own precedence.
    let click_from_zero_click: Option<String> = zero_click.map(click_value_to_profile_string);
    if let Some(click_str) = &click_from_zero_click {
        push(
            "--zero-click",
            format!("{click_str} (device click size, supplied on the command line)"),
            click_str.clone(),
            "elevation_click + windage_click",
        );
    }
    let mut zero_poi_up_m: Option<f64> = None;
    let mut zero_poi_right_m: Option<f64> = None;
    let mut zero_sets: Option<Vec<ProfileZeroSet>> = None;
    if src.zero_x_raw != 0 || src.zero_y_raw != 0 {
        match (zero_click, src.zero_distance_m) {
            (Some(click), Some(zero_distance_m)) if zero_distance_m > 0.0 => {
                let click_rad =
                    click.size / crate::adjustment::adjustment_factor(click.base);
                let up_clicks = f64::from(src.zero_y_raw) / 1000.0;
                let right_clicks = -f64::from(src.zero_x_raw) / 1000.0;
                let up_m = up_clicks * click_rad * zero_distance_m;
                let right_m = right_clicks * click_rad * zero_distance_m;
                zero_poi_up_m = Some(up_m);
                zero_poi_right_m = Some(right_m);
                // MBA-1360: ALSO record the click state as a zero set named "a7p-zero",
                // in DIAL-CORRECTION convention (the negated angular POI offset: a zero
                // state that impacts high/right needs less up/right dial). The engine-
                // field path above stays the primary consumer; see CLI_USAGE for why
                // selecting this set only makes sense on a profile whose zero_poi_*
                // fields have been cleared (both applied at once double-counts).
                zero_sets = Some(vec![ProfileZeroSet {
                    name: "a7p-zero".to_string(),
                    zero_distance: None,
                    poi_up_mil: Some(-(up_clicks * click_rad * 1000.0)),
                    poi_right_mil: Some(-(right_clicks * click_rad * 1000.0)),
                    notes: Some("imported .a7p zero_x/zero_y click state".to_string()),
                }]);
                push(
                    "zero_x / zero_y",
                    format!(
                        "({}, {}) raw = {:.2} right / {:.2} up device clicks",
                        src.zero_x_raw, src.zero_y_raw, right_clicks, up_clicks
                    ),
                    format!(
                        "POI {:.2} cm up / {:.2} cm right at {zero_distance_m:.0} m \
                         (--zero-click {}{})",
                        up_m * 100.0,
                        right_m * 100.0,
                        click.size,
                        match click.base {
                            crate::adjustment::ClickBase::Mil => "mil",
                            crate::adjustment::ClickBase::Moa => "moa",
                            crate::adjustment::ClickBase::Smoa => "smoa",
                        },
                    ),
                    "zero_poi_up_m + zero_poi_right_m + zero_sets[a7p-zero]",
                );
            }
            (Some(_), _) => {
                report.unmapped.push((
                    "zero_x / zero_y".to_string(),
                    format!(
                        "scope zeroing click offsets ({}, {}) — the file stores no zero \
                         distance, so --zero-click cannot convert them to a POI offset",
                        src.zero_x_raw, src.zero_y_raw
                    ),
                ));
            }
            // No --zero-click: the historical report line, byte-identical (the .a7p
            // format itself carries no click size to convert with).
            (None, _) => {
                report.unmapped.push((
                    "zero_x / zero_y".to_string(),
                    format!(
                        "scope zeroing click offsets ({}, {}) — device click size not \
                         supplied; pass --zero-click to record the profile's turret \
                         graduation and convert this offset",
                        src.zero_x_raw, src.zero_y_raw
                    ),
                ));
            }
        }
    }
    if !src.distances_m.is_empty() {
        report.unmapped.push((
            "distances".to_string(),
            format!("{} range-card entries (device UI list)", src.distances_m.len()),
        ));
    }
    if src.switches_count > 0 {
        report.unmapped.push((
            "switches".to_string(),
            format!("{} device UI switch entries", src.switches_count),
        ));
    }
    for (field, value) in [
        ("cartridge_name", &src.cartridge_name),
        ("caliber", &src.caliber),
        ("short_name_top", &src.short_name_top),
        ("short_name_bot", &src.short_name_bot),
        ("device_uuid", &src.device_uuid),
    ] {
        if !value.is_empty() {
            report
                .unmapped
                .push((field.to_string(), format!("\"{}\"", value.trim())));
        }
    }
    if !src.user_note.trim().is_empty() {
        report
            .unmapped
            .push(("user_note".to_string(), format!("\"{}\"", src.user_note.trim())));
    }
    for unknown in &doc.unknown_fields {
        report.unmapped.push((
            format!("{} field #{}", unknown.context, unknown.number),
            "unknown field (file newer than this importer)".to_string(),
        ));
    }

    let profile = ProfileData {
        name,
        velocity: src.muzzle_velocity_mps,
        bc,
        mass: src.bullet_weight_gr * GRAIN_TO_GRAM,
        diameter: src.bullet_diameter_in * IN_TO_MM,
        drag_model: drag_model.to_string(),
        twist_rate: Some(src.twist_in_per_turn * IN_TO_MM),
        sight_height: Some(src.sight_height_mm),
        zero_distance: src.zero_distance_m,
        units: "metric".to_string(),
        temperature: src.air_temperature_c,
        pressure: src.air_pressure_hpa,
        humidity: src.air_humidity_pct,
        altitude: 0.0,
        bullet_name: if src.bullet_name.is_empty() {
            None
        } else {
            Some(src.bullet_name.clone())
        },
        created: Some(timestamp_string()),
        wind_speed: None,
        wind_direction: None,
        shooting_angle: None,
        auto_zero: src.zero_distance_m,
        twist_right: Some(src.twist_right),
        // Left unset (not forced to Some(true)): consumers already treat "bc_segments
        // present" as implying velocity-segment behavior is active (the existing
        // `effective_use_bc_segments = use_bc_segments || bc_segments_data.is_some()`
        // pattern), so this boolean keeps its original meaning rather than being
        // overloaded by import.
        use_bc_segments: None,
        bullet_length: Some(src.bullet_length_in * IN_TO_MM),
        // MBA-1348: populated only when --zero-click supplied a device click size (see
        // the mapping above); otherwise left for `profile save
        // --elevation-click/--windage-click` (or a later edit) to fill in.
        elevation_click: click_from_zero_click.clone(),
        windage_click: click_from_zero_click,
        bc_segments,
        drag_curve,
        // .a7p carries no DSF concept; that's the `dsf` verb's job post-import.
        dsf_points: None,
        // .a7p carries no BC-reference-standard concept; imported BCs are treated as
        // ICAO-referenced (the omitted-field default) unless the user later edits the
        // saved profile with `--bc-reference`.
        bc_reference: None,
        // .a7p carries no QNH/pressure-reference concept; `air_pressure_hpa` is treated as
        // absolute station pressure (the omitted-field default) unless the user later edits
        // the saved profile with `--pressure-type`.
        pressure_reference: None,
        // .a7p carries no density-altitude concept; the imported temperature/pressure/altitude
        // are used as-is unless the user later edits the saved profile with
        // `--density-altitude`.
        density_altitude: None,
        // MBA-1359: populated only when --zero-click supplied a device click size to
        // convert the file's zero_x/zero_y click counts with (see the mapping above).
        zero_poi_up_m,
        zero_poi_right_m,
        // MBA-1396: .a7p has no lateral sight-offset concept; left for a hand edit.
        sight_offset_lateral_m: None,
        // MBA-1358: .a7p has no scope tracking-CF concept; derive with `tall-target`.
        elevation_cf: None,
        windage_cf: None,
        // MBA-1360: Some only when --zero-click converted zero_x/zero_y (see above).
        zero_sets,
        // MBA-1361: .a7p carries no reticle description; attach one after import with
        // `profile save --reticle-json`.
        reticle: None,
        // MBA-1348: .a7p carries no turret-mechanics or reticle-hold-bounds concept
        // beyond the click graduation above; these are for `profile save` (or a hand
        // edit) to fill in after import.
        clicks_per_revolution: None,
        zero_stop: None,
        elevation_travel_up_mil: None,
        elevation_travel_down_mil: None,
        windage_travel_left_mil: None,
        windage_travel_right_mil: None,
        turret_elevation_dialed_mil: None,
        turret_windage_dialed_mil: None,
        hold_bound_up_mil: None,
        hold_bound_down_mil: None,
        hold_bound_left_mil: None,
        hold_bound_right_mil: None,
    };

    Ok(A7pImportOutcome { profile, report })
}
