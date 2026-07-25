// WASM bindings for the ballistics engine with full CLI feature parity
use serde_json;
use wasm_bindgen::prelude::*;

use crate::bc_table_5d::Bc5dTable;
use crate::cli_api::{
    calculate_zero_angle_with_conditions, estimate_bc_fit, run_monte_carlo_with_direction_std_dev,
    AtmosphericConditions, BallisticInputs as InternalBallisticInputs, BcFitMode, MonteCarloParams,
    TrajectorySolver, WindConditions,
};
use crate::constants::GRAINS_PER_GRAM;
use crate::drag_model::DragModel;
use crate::moving_target::{calculate_lead, mover_ring};
use crate::truing_dsf::{apply_dsf, DsfPoint, DsfTable};
use std::cell::RefCell;

#[wasm_bindgen]
pub struct WasmBallistics {
    /// Optional BC5D correction table loaded from an in-memory `.bin`. When
    /// present, trajectory runs with `--use-bc-segments` synthesize
    /// velocity-dependent BC segments from it (offline parity with the online
    /// solver's ClusterBCDegradation + BC-segment + weather corrections).
    bc5d_table: RefCell<Option<Bc5dTable>>,
    /// Optional custom Mach:Cd drag table loaded from an in-memory CSV (MBA-1328;
    /// mirrors `bc5d_table` above). When present, every `trajectory`, `zero`, `lead`,
    /// and `monte-carlo` run applies it automatically as a full physical substitute
    /// for the BC + G-model drag (see `calculate_drag_coefficient`) — no `--use-*`
    /// gate flag needed, unlike BC5D's `--use-bc-segments`.
    drag_table: RefCell<Option<crate::drag::DragTable>>,
}

// Unit system for conversions
#[derive(Debug, Clone, Copy, PartialEq)]
enum UnitSystem {
    Imperial,
    Metric,
}

impl UnitSystem {
    fn from_str(s: &str) -> Self {
        match s.to_lowercase().as_str() {
            "metric" => UnitSystem::Metric,
            _ => UnitSystem::Imperial,
        }
    }
}

// Output format
#[derive(Debug, Clone, Copy, PartialEq)]
enum OutputFormat {
    Table,
    Json,
    Csv,
}

impl OutputFormat {
    fn from_str(s: &str) -> Self {
        match s.to_lowercase().as_str() {
            "json" => OutputFormat::Json,
            "csv" => OutputFormat::Csv,
            _ => OutputFormat::Table,
        }
    }
}

/// SMOA (a.k.a. IPHY) per MIL — the crate's locked MBA-724 dial-constant ratio
/// (`adjustment::adjustment_factor(Smoa) / adjustment_factor(Mil)` == 3.6 exact).
/// Mirrors native `main.rs`'s `smoa_per_mil()` helper so the WASM terminal's SMOA/IPHY
/// figures agree with the native CLI bit-for-bit (MBA-1355).
fn smoa_per_mil() -> f64 {
    crate::adjustment::adjustment_factor(crate::adjustment::ClickBase::Smoa)
        / crate::adjustment::adjustment_factor(crate::adjustment::ClickBase::Mil)
}

/// How the mover Ring table column (trajectory's `--target-speed`) turns its raw mil
/// angle into a display value for the active `--adjustment-unit` (MBA-1355), mirroring
/// native `main.rs`'s `RingUnit`: every unit except Clicks is a constant multiply on the
/// mil angle; Clicks instead rounds to a whole click count via a resolved `ClickValue`
/// (the elevation graduation — Ring isn't cleanly an elevation or windage axis, so it
/// reuses elevation, same choice the native CLI's PDF dope card Drop column makes).
#[derive(Debug, Clone, Copy)]
enum RingDisplayUnit {
    Factor(f64, &'static str),
    Clicks(crate::adjustment::ClickValue),
}

/// Parse a `--powder-temp-curve` "TEMP:VEL,..." value into sorted SI
/// `(temperature_c, velocity_m_s)` points — ONE parser shared by the trajectory and
/// lead command handlers so their validation cannot drift (review fix M2: the two
/// inline copies had already drifted from native `parse_powder_temp_curve` by
/// accepting duplicate temperatures, which make the interpolation order-dependent).
/// Validation mirrors the native parser: TEMP:VELOCITY shape, positive velocity,
/// at least 2 points, and no duplicate temperatures after sorting.
fn parse_powder_temp_curve_str(
    curve_str: &str,
    units: UnitSystem,
) -> Result<Vec<(f64, f64)>, JsValue> {
    let mut pts: Vec<(f64, f64)> = Vec::new();
    for raw in curve_str.split(',') {
        let part = raw.trim();
        if part.is_empty() {
            continue;
        }
        let (t_str, v_str) = part
            .split_once(':')
            .ok_or_else(|| JsValue::from_str("--powder-temp-curve point must be TEMP:VELOCITY"))?;
        let t: f64 = t_str
            .trim()
            .parse()
            .map_err(|_| JsValue::from_str("Invalid temperature in --powder-temp-curve"))?;
        let v: f64 = v_str
            .trim()
            .parse()
            .map_err(|_| JsValue::from_str("Invalid velocity in --powder-temp-curve"))?;
        if !(v > 0.0) {
            return Err(JsValue::from_str(
                "--powder-temp-curve velocity must be positive",
            ));
        }
        let (t_c, v_ms) = match units {
            UnitSystem::Imperial => ((t - 32.0) * 5.0 / 9.0, v * 0.3048),
            UnitSystem::Metric => (t, v),
        };
        pts.push((t_c, v_ms));
    }
    if pts.len() < 2 {
        return Err(JsValue::from_str(
            "--powder-temp-curve needs at least 2 TEMP:VELOCITY points",
        ));
    }
    pts.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    for w in pts.windows(2) {
        if (w[0].0 - w[1].0).abs() < f64::EPSILON {
            return Err(JsValue::from_str(
                "--powder-temp-curve has duplicate temperatures",
            ));
        }
    }
    Ok(pts)
}

/// Self-describing units/axes metadata for [`format_trajectory_json`]'s `-o json` output
/// (MBA-1315), additive: appended as a third top-level key alongside the pre-existing
/// `trajectory`/`summary` keys, which are unchanged. The WASM `-o json` fields already carry
/// their unit in the key name (`range_yards`, `drop_inches`, `drift_inches`, ...), so this
/// block's `units` restates those labels generically (matching the native CLI legacy JSON's
/// `legend.units`, MBA-1315) and its `axes` states the sign convention those field names don't
/// spell out on their own.
///
/// Sign conventions verified against `crate::wind::wind_vector`'s documented McCoy frame
/// (x downrange, y up, z right) and the native CLI's crosswind recon (MBA-1315): `drop` is
/// LOS-relative (`los_height - position.y`), positive when the bullet is below the line of
/// sight; `drift` is `position.z` directly, positive to the shooter's right, identical in
/// sign and source to the native CLI legacy JSON's `x`.
fn trajectory_json_legend(units: UnitSystem) -> serde_json::Value {
    let (distance, drop, drift, velocity, energy) = match units {
        UnitSystem::Imperial => ("yd", "in", "in", "fps", "ft-lb"),
        UnitSystem::Metric => ("m", "cm", "cm", "m/s", "J"),
    };
    serde_json::json!({
        "units": {
            "range": distance,
            "drop": drop,
            "drift": drift,
            "velocity": velocity,
            "energy": energy
        },
        "axes": {
            "range": "downrange distance from the muzzle; zero at the muzzle, always increasing.",
            "drop": "vertical miss from the line of sight; positive means the bullet is below \
                     the line of sight (has fallen below the aim point). Not the same reference \
                     as the native CLI legacy JSON's world-frame `y`.",
            "drift": "lateral miss from the line of sight; positive means to the shooter's \
                      right (e.g. a crosswind FROM the left, --wind-direction 270, drifts it \
                      positive; FROM the right, --wind-direction 90, drifts it negative). Same \
                      sign and source as the native CLI legacy JSON's `x`."
        }
    })
}

/// Map the WASM terminal's session unit system onto the engine-side
/// [`crate::cli_api::UnitSystem`] shared by the truing and WEZ cores (MBA-1343).
fn engine_units(units: UnitSystem) -> crate::cli_api::UnitSystem {
    match units {
        UnitSystem::Imperial => crate::cli_api::UnitSystem::Imperial,
        UnitSystem::Metric => crate::cli_api::UnitSystem::Metric,
    }
}

/// Fetch the value token for the value-taking flag at `args[i]`, or error when the
/// command line ends right after the flag. Mirrors native clap's "a value is
/// required" failure; the hand-rolled parsers previously skipped such a dangling
/// flag silently, which corrupted fits (a dangling `--observed` silently fell
/// back to single-observation truing — MBA-1343 review).
fn require_value<'a>(args: &[&'a str], i: usize) -> Result<&'a str, JsValue> {
    if i + 1 < args.len() {
        Ok(args[i + 1])
    } else {
        Err(JsValue::from_str(&format!(
            "Error: a value is required for '{}'",
            args[i]
        )))
    }
}

/// MBA-1356: pairing-requirement text for `--cd-scale` without a loaded drag table, shared by
/// every WASM terminal command that accepts it (trajectory/zero/monte-carlo). Mirrors the
/// native CLI's `resolve_cd_scale` text — except the WASM terminal has no `--bc-adjustment`
/// flag on ANY surface (unlike native, where `trajectory` has one), so unlike native's
/// text this never suggests it (0.28.1 sweep: surface-accurate suggestions only).
const CD_SCALE_REQUIRES_DRAG_TABLE: &str = "--cd-scale requires --drag-table";

/// MBA-1356: nudge for a `--cd-scale` far outside the documented typical truing range
/// (0.90-1.10), mirroring the native CLI's `cd_scale_range_warning` (src/main.rs) text exactly.
/// `None` inside `[0.5, 2.0]` — the engine's own gate (`require_positive` in
/// `validate_for_solve`) already covers finite && > 0. The trailing `\n\n` matches the shape of
/// the BC5D coercion warning (MBA-1386, bcdd213) since both are prepended, table-only, ahead of
/// the run's normal output.
fn cd_scale_range_warning(value: f64) -> Option<String> {
    if (0.5..=2.0).contains(&value) {
        None
    } else {
        Some(format!(
            "warning: --cd-scale {value} is far outside the typical truing range (0.90-1.10)\n\n"
        ))
    }
}

/// `"{n} points, Mach {lo:.2}-{hi:.2}"` summary of a DSF table's points — mirrors native
/// `main.rs`'s `dsf_table_summary` exactly (byte-identical wording), so the WASM
/// terminal's "DSF table active (...)" note reads the same as the native CLI's
/// (MBA-1411). Kept as a local copy since the native version is private to the binary
/// crate. `points` empty reports Mach 0.00-0.00 rather than panicking.
fn dsf_table_summary(points: &[DsfPoint]) -> String {
    let lo = points.first().map(|p| p.mach).unwrap_or(0.0);
    let hi = points.last().map(|p| p.mach).unwrap_or(0.0);
    format!("{} points, Mach {:.2}-{:.2}", points.len(), lo, hi)
}

/// Parse every `--dsf-point MACH:DSF` occurrence into a validated [`DsfTable`] (MBA-1411:
/// WASM has no saved-profile storage, so the DSF table is a per-call argument, mirroring
/// the MBA-1343 per-call pattern used by `true-velocity`'s multi-`--observed` and
/// `monte-carlo --wez`). `None` when no `--dsf-point` flags were given — the caller then
/// skips `apply_dsf` entirely, an exact no-op. Bounds/cap validation (0 < mach < 1.2,
/// 0.5 < dsf < 2.0, <= 6 points) is NOT reimplemented here — `DsfTable::from_points`'s
/// error strings ARE the user-facing errors, identical to the native `dsf`/saved-profile
/// paths.
fn parse_dsf_points(raw: &[String]) -> Result<Option<DsfTable>, JsValue> {
    if raw.is_empty() {
        return Ok(None);
    }
    let mut points = Vec::with_capacity(raw.len());
    for s in raw {
        let parts: Vec<&str> = s.split(':').collect();
        if parts.len() != 2 {
            return Err(JsValue::from_str(&format!(
                "--dsf-point expects MACH:DSF (e.g. 0.95:1.08), got '{}'",
                s
            )));
        }
        let mach: f64 = parts[0]
            .trim()
            .parse()
            .map_err(|_| JsValue::from_str(&format!("--dsf-point: invalid MACH in '{}'", s)))?;
        let dsf: f64 = parts[1]
            .trim()
            .parse()
            .map_err(|_| JsValue::from_str(&format!("--dsf-point: invalid DSF in '{}'", s)))?;
        points.push(DsfPoint { mach, dsf });
    }
    DsfTable::from_points(points)
        .map(Some)
        .map_err(|e| JsValue::from_str(&e))
}

/// Local mirror of the native CLI's `MonteCarloOutput` modes for the WEZ path
/// (MBA-1343 Phase C). The WASM `-o` spellings map onto these in
/// `run_monte_carlo_wez`.
enum WezOutputMode {
    Summary,
    Statistics,
    Full,
}

/// Dominant-bucket label shared by the WEZ summary-table and statistics-CSV
/// renderers — replicates the native CLI's `wez_dominant_label` (main.rs) exactly.
fn wez_dominant_label(row: &crate::wez::WezRow) -> String {
    if row.attribution_unavailable {
        "n/a".to_string()
    } else {
        row.dominant_error_source
            .map(|b| b.to_string())
            .unwrap_or_else(|| "none".to_string())
    }
}

/// Meters -> user distance units (native `UnitConverter::distance_from_metric`).
fn distance_from_metric(val: f64, units: UnitSystem) -> f64 {
    match units {
        UnitSystem::Metric => val,
        UnitSystem::Imperial => val / 0.9144, // meters to yards
    }
}

/// m/s -> user wind units (native `UnitConverter::wind_from_metric`).
fn wind_from_metric(val: f64, units: UnitSystem) -> f64 {
    match units {
        UnitSystem::Metric => val,
        UnitSystem::Imperial => val / 0.44704, // m/s to mph
    }
}

/// `-o full`: the entire [`crate::wez::WezResult`] as pretty JSON — replicates the
/// native `print_wez_full` (the 0.25.0 JSON contract) byte-for-byte, including the
/// trailing newline native's `println!` emits.
fn format_wez_full(result: &crate::wez::WezResult) -> Result<String, JsValue> {
    let json = serde_json::to_string_pretty(result)
        .map_err(|e| JsValue::from_str(&format!("Error serializing JSON: {e}")))?;
    Ok(format!("{json}\n"))
}

/// `-o statistics`: one CSV row per range step — replicates the native
/// `print_wez_statistics` byte-for-byte.
fn format_wez_statistics(result: &crate::wez::WezResult, units: UnitSystem) -> String {
    let mut out = String::new();
    out.push_str("range,p_hit,dominant_error_source,wind_call_share,mv_sd_share,other_share\n");
    for row in &result.rows {
        out.push_str(&format!(
            "{:.2},{:.4},{},{:.4},{:.4},{:.4}\n",
            distance_from_metric(row.range_m, units),
            row.p_hit,
            wez_dominant_label(row),
            row.wind_call_share,
            row.mv_sd_share,
            row.other_share
        ));
    }
    out
}

/// Default `-o summary`: the human-readable sweep table — replicates the native
/// `print_wez_summary` byte-for-byte.
fn format_wez_summary(result: &crate::wez::WezResult, units: UnitSystem) -> String {
    let dist_unit = match units {
        UnitSystem::Imperial => "yd",
        UnitSystem::Metric => "m",
    };
    let wind_unit = match units {
        UnitSystem::Imperial => "mph",
        UnitSystem::Metric => "m/s",
    };
    let mut out = String::new();
    out.push_str(&format!(
        "WEZ sweep: {} sims/step, wind call {:.2} {wind_unit} + wind std {:.2} {wind_unit} (quadrature) = {:.2} {wind_unit} effective\n",
        result.num_sims_per_step,
        wind_from_metric(result.wind_call_error_mps, units),
        wind_from_metric(result.wind_speed_std_mps, units),
        wind_from_metric(result.combined_wind_speed_std_mps, units),
    ));
    out.push_str(
        "┌────────────┬──────────┬───────────────┬───────────┬───────────┬───────────┐\n",
    );
    out.push_str(&format!(
        "│ Range ({dist_unit:>3}) │  P(hit)  │ Dominant      │ Wind call │  MV SD    │ Other/grp │\n"
    ));
    out.push_str(
        "├────────────┼──────────┼───────────────┼───────────┼───────────┼───────────┤\n",
    );
    for row in &result.rows {
        out.push_str(&format!(
            "│ {:>10.1} │ {:>7.1}% │ {:<13} │ {:>8.1}% │ {:>8.1}% │ {:>8.1}% │\n",
            distance_from_metric(row.range_m, units),
            row.p_hit * 100.0,
            wez_dominant_label(row),
            row.wind_call_share * 100.0,
            row.mv_sd_share * 100.0,
            row.other_share * 100.0,
        ));
    }
    out.push_str(
        "└────────────┴──────────┴───────────────┴───────────┴───────────┴───────────┘\n",
    );
    out
}

/// Render a multi-observation truing report as table / JSON / CSV — replicates the
/// native CLI's `display_multi_truing_result` (main.rs) byte-for-byte so terminal
/// parity can be tested against the native binary (MBA-1343 Phase C).
fn format_multi_truing_result(
    report: &crate::truing::MultiTruingReport,
    drop_unit: crate::truing::DropUnit,
    units: UnitSystem,
    chrono_fps: Option<f64>,
    output: OutputFormat,
) -> String {
    let vel_unit = match units {
        UnitSystem::Imperial => "fps",
        UnitSystem::Metric => "m/s",
    };
    let range_unit = match units {
        UnitSystem::Imperial => "yd",
        UnitSystem::Metric => "m",
    };
    let drop_label = drop_unit.label();
    let mv_display = match units {
        UnitSystem::Imperial => report.fitted_mv_fps,
        UnitSystem::Metric => report.fitted_mv_fps * 0.3048,
    };
    let range_display = |range_yd: f64| match units {
        UnitSystem::Imperial => range_yd,
        UnitSystem::Metric => range_yd * 0.9144,
    };
    // Chrono comparison (chrono_fps is already in fps).
    let (adj_display, adj_pct) = match chrono_fps {
        Some(c) => {
            let adj_fps = report.fitted_mv_fps - c;
            let adj_disp = match units {
                UnitSystem::Imperial => adj_fps,
                UnitSystem::Metric => adj_fps * 0.3048,
            };
            let pct = if c != 0.0 { adj_fps / c * 100.0 } else { 0.0 };
            (Some(adj_disp), Some(pct))
        }
        None => (None, None),
    };

    // MBA-1405 Task 2: MV-calibration window (90-100% of the downward Mach 1.2
    // crossing) — same additive JSON fields and table lines as the native CLI's
    // `display_multi_truing_result`. The native out-of-window WARNING is a stderr
    // diagnostic with no WASM equivalent (this terminal has no stderr channel and
    // never reproduces native's stderr text — see the "Fitting N observations..."
    // progress line, also absent here), so it is deliberately not mirrored.
    let mv_window = crate::truing::mv_calibration_window(report.mach_1_2_distance_m);

    let mut out = String::new();
    match output {
        OutputFormat::Json => {
            let obs_json: Vec<serde_json::Value> = report
                .observations
                .iter()
                .enumerate()
                .map(|(i, o)| {
                    serde_json::json!({
                        format!("range_{range_unit}"): range_display(o.range_yd),
                        format!("observed_drop_{drop_label}"): o.drop,
                        format!("predicted_drop_{drop_label}"): report.predicted[i],
                        format!("residual_{drop_label}"): report.residuals[i],
                    })
                })
                .collect();

            let mut json_output = serde_json::json!({
                "mode": if report.bc_fitted { "joint_mv_bc" } else { "mv_only" },
                "fitted_muzzle_velocity": mv_display,
                "velocity_unit": vel_unit,
                "bc_fitted": report.bc_fitted,
                "fitted_bc": report.fitted_bc,
                "input_bc": report.bc_input,
                "observations": obs_json,
                format!("rms_residual_{drop_label}"): report.rms,
                "iterations": report.iterations,
                "converged": report.converged,
                "bc_sensitivity_ratio": report.sensitivity_ratio,
                "condition_number": if report.condition_number.is_finite() {
                    serde_json::json!(report.condition_number)
                } else {
                    serde_json::Value::Null
                },
                "quality": report.quality,
                "legend": {
                    "units": {
                        "range": range_unit,
                        "drop": drop_label,
                        "velocity": vel_unit,
                    },
                },
            });
            if !report.reason.is_empty() {
                json_output["bc_hold_reason"] = serde_json::json!(report.reason);
            }
            if let Some(adj) = adj_display {
                json_output["velocity_adjustment"] = serde_json::json!(adj);
            }
            if let Some(pct) = adj_pct {
                json_output["adjustment_percent"] = serde_json::json!(pct);
            }
            // MBA-1405 Task 2: additive-only fields (meters, unit-invariant); null when
            // the trajectory never crosses Mach 1.2 — no note text ever accompanies
            // these in JSON (purity rule), matching the native CLI's JSON exactly.
            json_output["mv_window_start_m"] = match mv_window {
                Some((lo_m, _)) => serde_json::json!(lo_m),
                None => serde_json::Value::Null,
            };
            json_output["mv_window_end_m"] = match mv_window {
                Some((_, hi_m)) => serde_json::json!(hi_m),
                None => serde_json::Value::Null,
            };
            // Native prints a serialization error to stderr and leaves stdout empty;
            // the returned string mirrors stdout exactly either way.
            if let Ok(s) = serde_json::to_string_pretty(&json_output) {
                out.push_str(&s);
                out.push('\n');
            }
        }
        OutputFormat::Csv => {
            out.push_str(&format!(
                "range_{range_unit},observed_drop_{drop_label},predicted_drop_{drop_label},residual_{drop_label}\n"
            ));
            for (i, o) in report.observations.iter().enumerate() {
                out.push_str(&format!(
                    "{:.1},{:.4},{:.4},{:+.4}\n",
                    range_display(o.range_yd),
                    o.drop,
                    report.predicted[i],
                    report.residuals[i]
                ));
            }
            out.push('\n');
            out.push_str(&format!(
                "fitted_muzzle_velocity_{vel_unit},bc_fitted,fitted_bc,input_bc,rms_residual_{drop_label},iterations,converged,bc_sensitivity_ratio,condition_number\n"
            ));
            out.push_str(&format!(
                "{:.1},{},{:.4},{:.4},{:.4},{},{},{:.4},{}\n",
                mv_display,
                report.bc_fitted,
                report.fitted_bc,
                report.bc_input,
                report.rms,
                report.iterations,
                report.converged,
                report.sensitivity_ratio,
                if report.condition_number.is_finite() {
                    format!("{:.1}", report.condition_number)
                } else {
                    "inf".to_string()
                },
            ));
        }
        OutputFormat::Table => {
            out.push('\n');
            out.push_str("=== VELOCITY + BC TRUING (multi-observation) ===\n");
            out.push('\n');
            out.push_str(&format!(
                "  Fitted muzzle velocity: {:>9.1} {}\n",
                mv_display, vel_unit
            ));
            if report.bc_fitted {
                out.push_str(&format!(
                    "  Fitted BC:              {:>9.4}  (input {:.4})\n",
                    report.fitted_bc, report.bc_input
                ));
            } else {
                out.push_str(&format!(
                    "  BC:                     {:>9.4}  (held; not fitted)\n",
                    report.fitted_bc
                ));
            }
            if let Some(adj) = adj_display {
                out.push_str(&format!(
                    "  Adjustment from chrono: {:>+9.1} {}\n",
                    adj, vel_unit
                ));
                if let Some(pct) = adj_pct {
                    out.push_str(&format!("  Adjustment percentage:  {:>+9.2}%\n", pct));
                }
            }
            out.push('\n');
            out.push_str(&format!(
                "  {:>10}  {:>14}  {:>14}  {:>12}\n",
                format!("Range ({range_unit})"),
                format!("Observed ({drop_label})"),
                format!("Predicted ({drop_label})"),
                format!("Resid ({drop_label})"),
            ));
            out.push_str(&format!("  {}\n", "-".repeat(56)));
            for (i, o) in report.observations.iter().enumerate() {
                out.push_str(&format!(
                    "  {:>10.1}  {:>14.3}  {:>14.3}  {:>+12.3}\n",
                    range_display(o.range_yd),
                    o.drop,
                    report.predicted[i],
                    report.residuals[i]
                ));
            }
            out.push_str(&format!("  {}\n", "-".repeat(56)));
            out.push_str(&format!(
                "  RMS residual: {:.3} {}   |   iterations: {}{}\n",
                report.rms,
                drop_label,
                report.iterations,
                if report.converged {
                    ""
                } else {
                    " (not fully converged)"
                }
            ));
            out.push('\n');
            out.push_str(&format!("  {}\n", report.quality));
            if !report.reason.is_empty() {
                out.push_str(&format!("  Note: {}\n", report.reason));
            }
            out.push_str(&format!(
                "  Diagnostics: BC sensitivity ratio {:.4}, conditioning {}\n",
                report.sensitivity_ratio,
                if report.condition_number.is_finite() {
                    format!("{:.0}", report.condition_number)
                } else {
                    "inf".to_string()
                }
            ));
            // MBA-1405 Task 2: MV-calibration window / no-window note, table only —
            // mirrors the native CLI's `display_multi_truing_result` exactly.
            match mv_window {
                Some((lo_m, hi_m)) => {
                    let lo = distance_from_metric(lo_m, units);
                    let hi = distance_from_metric(hi_m, units);
                    out.push_str(&format!(
                        "  MV-calibration window: {lo:.1}-{hi:.1} {range_unit} (90-100% of the Mach 1.2 distance)\n"
                    ));
                }
                None if report.muzzle_mach >= 1.2 => {
                    let max_display = distance_from_metric(report.window_solved_range_m, units);
                    out.push_str(&format!(
                        "  note: no MV window: trajectory is supersonic through {max_display:.1} \
                         {range_unit}; MV is identifiable at any range\n"
                    ));
                }
                None => {
                    out.push_str(
                        "  note: no MV window: trajectory never reaches Mach 1.2; calibrate \
                         muzzle velocity with a chronograph, then collect DSF points\n",
                    );
                }
            }
            out.push_str("  for optimal observation ranges run: ballistics plan-truing\n");
            out.push('\n');
        }
    }
    out
}

/// Render a single-observation truing result as table / JSON / CSV — replicates the
/// native CLI's `display_true_velocity_result` (main.rs) byte-for-byte.
#[allow(
    clippy::too_many_arguments,
    reason = "flat arguments mirror the native display_true_velocity_result signature"
)]
fn format_true_velocity_result(
    effective_vel: f64,
    vel_unit: &str,
    velocity_adjustment: Option<f64>,
    adjustment_percent: Option<f64>,
    confidence: &str,
    iterations: i32,
    final_error_mil: f64,
    calculated_drop_mil: f64,
    measured_drop: f64,
    units: UnitSystem,
    output: OutputFormat,
    used_bc_table: bool,
) -> Result<String, JsValue> {
    let mut out = String::new();
    match output {
        OutputFormat::Json => {
            let mut json_output = serde_json::json!({
                "effective_velocity": effective_vel,
                "velocity_unit": vel_unit,
                "confidence": confidence,
                "iterations": iterations,
                "final_error_mil": final_error_mil,
                "calculated_drop_mil": calculated_drop_mil,
                "measured_drop_mil": measured_drop,
                "used_bc_table": used_bc_table,
            });
            if let Some(adj) = velocity_adjustment {
                let adj_display = match units {
                    UnitSystem::Imperial => adj,
                    UnitSystem::Metric => adj * 0.3048,
                };
                json_output["velocity_adjustment"] = serde_json::json!(adj_display);
            }
            if let Some(pct) = adjustment_percent {
                json_output["adjustment_percent"] = serde_json::json!(pct);
            }
            let s = serde_json::to_string_pretty(&json_output)
                .map_err(|e| JsValue::from_str(&format!("Error serializing JSON: {e}")))?;
            out.push_str(&s);
            out.push('\n');
        }
        OutputFormat::Csv => {
            out.push_str("effective_velocity,unit,adjustment,adjustment_pct,confidence,iterations,final_error_mil,calculated_drop_mil,used_bc_table\n");
            out.push_str(&format!(
                "{:.1},{},{},{},{},{},{:.4},{:.2},{}\n",
                effective_vel,
                vel_unit,
                velocity_adjustment
                    .map(|v| {
                        let adj = match units {
                            UnitSystem::Imperial => v,
                            UnitSystem::Metric => v * 0.3048,
                        };
                        format!("{:.1}", adj)
                    })
                    .unwrap_or_default(),
                adjustment_percent
                    .map(|v| format!("{:.2}", v))
                    .unwrap_or_default(),
                confidence,
                iterations,
                final_error_mil,
                calculated_drop_mil,
                used_bc_table,
            ));
        }
        OutputFormat::Table => {
            out.push('\n');
            out.push_str("╔════════════════════════════════════════════════════════════╗\n");
            out.push_str("║           VELOCITY TRUING RESULTS                          ║\n");
            out.push_str("╠════════════════════════════════════════════════════════════╣\n");
            out.push_str(&format!(
                "║  Effective Muzzle Velocity: {:>8.1} {:>4}                 ║\n",
                effective_vel, vel_unit
            ));
            if let Some(adj) = velocity_adjustment {
                let adj_display = match units {
                    UnitSystem::Imperial => adj,
                    UnitSystem::Metric => adj * 0.3048,
                };
                out.push_str(&format!(
                    "║  Adjustment from Chrono:    {:>+8.1} {:>4}                 ║\n",
                    adj_display, vel_unit
                ));
                if let Some(pct) = adjustment_percent {
                    out.push_str(&format!(
                        "║  Adjustment Percentage:     {:>+8.2}%                      ║\n",
                        pct
                    ));
                }
            }
            out.push_str("╠════════════════════════════════════════════════════════════╣\n");
            out.push_str(&format!(
                "║  Confidence:                {:>8}                        ║\n",
                confidence
            ));
            out.push_str(&format!(
                "║  Iterations:                {:>8}                        ║\n",
                iterations
            ));
            out.push_str(&format!(
                "║  Final Error:               {:>8.4} MIL                  ║\n",
                final_error_mil
            ));
            out.push_str(&format!(
                "║  Calculated Drop:           {:>8.2} MIL                  ║\n",
                calculated_drop_mil
            ));
            out.push_str(&format!(
                "║  Measured Drop:             {:>8.2} MIL                  ║\n",
                measured_drop
            ));
            if used_bc_table {
                out.push_str("╠════════════════════════════════════════════════════════════╣\n");
                out.push_str("║  BC5D Table:                     Yes                       ║\n");
            }
            out.push_str("╚════════════════════════════════════════════════════════════╝\n");
            out.push('\n');
        }
    }
    Ok(out)
}

#[wasm_bindgen]
impl WasmBallistics {
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        WasmBallistics {
            bc5d_table: RefCell::new(None),
            drag_table: RefCell::new(None),
        }
    }

    /// Load a BC5D correction table from the raw bytes of a `bc5d_<caliber>.bin`
    /// file. The host (browser `fetch()` or Node `fs`/`fetch`) is responsible
    /// for retrieving the file — WASM has no filesystem or network — and passes
    /// the bytes here.
    ///
    /// Once loaded, any `trajectory` run that includes `--use-bc-segments` will
    /// apply velocity-dependent BC segments synthesized from this table. Load a
    /// table matching the bullet's caliber (e.g. `bc5d_308.bin` for a .308).
    ///
    /// Returns a short human-readable summary of the loaded table. Replaces any
    /// previously loaded table.
    #[wasm_bindgen(js_name = loadBc5dTable)]
    pub fn load_bc5d_table(&self, bytes: &[u8]) -> Result<String, JsValue> {
        let table = Bc5dTable::from_bytes(bytes)
            .map_err(|e| JsValue::from_str(&format!("Failed to parse BC5D table: {}", e)))?;
        let summary = format!(
            "Loaded BC5D table: caliber {:.3}\", {} cells, api_version {}",
            table.caliber(),
            table.total_cells(),
            table.api_version()
        );
        *self.bc5d_table.borrow_mut() = Some(table);
        Ok(summary)
    }

    /// Report whether a BC5D table is currently loaded.
    #[wasm_bindgen(js_name = hasBc5dTable)]
    pub fn has_bc5d_table(&self) -> bool {
        self.bc5d_table.borrow().is_some()
    }

    /// Load a custom Mach:Cd drag table (MBA-1328) — a measured or manufacturer-published
    /// drag curve (Hornady CDM data, a Lapua/Doppler-radar-derived deck, or your own) — from
    /// the raw bytes of a CSV file. The host (browser `fetch()` or Node `fs`/`fetch`) is
    /// responsible for retrieving the file — WASM has no filesystem or network — and passes
    /// the bytes here, mirroring [`Self::load_bc5d_table`].
    ///
    /// The bytes must be valid UTF-8 CSV text; parsing is delegated to
    /// [`crate::drag::DragTable::from_csv_str`] — the SAME parser native `--drag-table`
    /// uses — so the accepted format is identical: two columns `mach,cd` per line, blank
    /// lines and `#` comments ignored, a single leading textual header row tolerated, Mach
    /// strictly ascending with at least 2 points, Cd finite and > 0.
    ///
    /// Once loaded, EVERY subsequent `trajectory`, `zero`, `lead`, and `monte-carlo` run
    /// applies the table automatically — no `--use-*` gate flag is needed (unlike a loaded
    /// BC5D table, which only takes effect with `--use-bc-segments`). A loaded table is a
    /// full physical substitute for the G-model + BC (see `calculate_drag_coefficient`):
    /// `-b`/`--bc` may still be supplied but is ignored for drag once a table is active
    /// (matching the native `--drag-table` CLI semantics documented in CLI_USAGE.md).
    ///
    /// Returns a short human-readable summary of the loaded table (point count + Mach
    /// range). Replaces any previously loaded table.
    ///
    /// MBA-1409: also accepts `.drg` vendor drag-curve text (the same format the native
    /// `--drag-table` CLI accepts by `.drg` file extension) as a fallback. WASM has no
    /// filesystem and thus no extension to dispatch on, so the bytes are tried as CSV first
    /// (exactly as before this change); only on CSV failure, if the text
    /// [`crate::drag_file::looks_like_drg`], it is retried through
    /// [`crate::drag_file::parse_drg`]. If both fail, the returned error names both formats.
    #[wasm_bindgen(js_name = loadDragTable)]
    pub fn load_drag_table(&self, bytes: &[u8]) -> Result<String, JsValue> {
        let text = std::str::from_utf8(bytes).map_err(|e| {
            JsValue::from_str(&format!("Drag table bytes are not valid UTF-8: {}", e))
        })?;
        let table = match crate::drag::DragTable::from_csv_str(text) {
            Ok(t) => t,
            Err(csv_err) => {
                if !crate::drag_file::looks_like_drg(text) {
                    return Err(JsValue::from_str(&format!(
                        "Failed to parse drag table as CSV ({csv_err}); text does not look like \
                         a .drg file either"
                    )));
                }
                match crate::drag_file::parse_drg(text) {
                    Ok(curve) => {
                        let (mach, cd): (Vec<f64>, Vec<f64>) = curve.points.into_iter().unzip();
                        crate::drag::DragTable::try_new(mach, cd).map_err(|e| {
                            JsValue::from_str(&format!("Failed to parse drag table: {}", e))
                        })?
                    }
                    Err(drg_err) => {
                        return Err(JsValue::from_str(&format!(
                            "Failed to parse drag table as CSV ({csv_err}) or as .drg ({drg_err})"
                        )));
                    }
                }
            }
        };
        let summary = format!(
            "Loaded drag table: {} points, Mach {:.3}-{:.3}",
            table.mach_values.len(),
            table.mach_values.first().copied().unwrap_or(0.0),
            table.mach_values.last().copied().unwrap_or(0.0),
        );
        *self.drag_table.borrow_mut() = Some(table);
        Ok(summary)
    }

    /// Report whether a custom drag table is currently loaded.
    #[wasm_bindgen(js_name = hasDragTable)]
    pub fn has_drag_table(&self) -> bool {
        self.drag_table.borrow().is_some()
    }

    /// Unload any custom drag table previously installed via [`Self::load_drag_table`],
    /// reverting every subsequent `trajectory`, `zero`, `lead`, and `monte-carlo` run to the
    /// standard G-model + BC drag (the `-b`/`--bc` value with the selected G1/G7 curve).
    ///
    /// The inverse of `loadDragTable`. It lets a single engine instance alternate between a
    /// measured-curve (CDM) solve and a plain G7-BC solve without constructing a second
    /// instance: `load → solve CDM → clear → solve G7`, with the G7 half uncontaminated —
    /// the same load/clear pattern `clearWindSegments` provides for segmented wind.
    ///
    /// Returns `true` if a table was loaded (and is now cleared), `false` if none was loaded
    /// (a harmless no-op). Idempotent.
    #[wasm_bindgen(js_name = clearDragTable)]
    pub fn clear_drag_table(&self) -> bool {
        self.drag_table.borrow_mut().take().is_some()
    }

    /// Run a command and return the output
    #[wasm_bindgen(js_name = runCommand)]
    pub fn run_command(&self, command: &str) -> Result<String, JsValue> {
        // Handle quoted arguments properly
        let mut args: Vec<String> = Vec::new();
        let mut current_arg = String::new();
        let mut in_quotes = false;
        let mut quote_char = ' ';

        for c in command.chars() {
            if !in_quotes && (c == '\'' || c == '"') {
                in_quotes = true;
                quote_char = c;
            } else if in_quotes && c == quote_char {
                in_quotes = false;
                quote_char = ' ';
            } else if !in_quotes && c.is_whitespace() {
                if !current_arg.is_empty() {
                    args.push(current_arg.clone());
                    current_arg.clear();
                }
            } else {
                current_arg.push(c);
            }
        }

        if !current_arg.is_empty() {
            args.push(current_arg);
        }

        let args: Vec<&str> = args.iter().map(|s| s.as_str()).collect();

        if args.is_empty() {
            return Ok(self.show_help());
        }

        // Skip "ballistics" if present
        let args = if !args.is_empty() && (args[0] == "ballistics" || args[0] == "./ballistics") {
            &args[1..]
        } else {
            &args[..]
        };

        if args.is_empty() || args[0] == "help" || args[0] == "--help" || args[0] == "-h" {
            return Ok(self.show_help());
        }

        // Parse global unit system first
        let mut units = UnitSystem::Imperial;
        for i in 0..args.len() {
            if args[i] == "--units" || args[i] == "-u" {
                if i + 1 < args.len() {
                    units = UnitSystem::from_str(args[i + 1]);
                }
                break;
            }
        }

        // MBA-1289: `--units <system>` may appear anywhere — including BEFORE the command
        // word, the order the help text advertises. Units were already parsed from the
        // full list above, so strip the flag and its value here; otherwise the dispatch
        // below would read `--units` as the command and fail with "Unknown command".
        let mut stripped: Vec<&str> = Vec::with_capacity(args.len());
        let mut i = 0;
        while i < args.len() {
            if args[i] == "--units" || args[i] == "-u" {
                i += 2; // skip the flag and its value (a dangling flag skips harmlessly)
            } else {
                stripped.push(args[i]);
                i += 1;
            }
        }
        let args = stripped;
        if args.is_empty() {
            return Ok(self.show_help());
        }

        // Route to appropriate command handler
        match args[0] {
            "version" => Ok(format!(
                "Ballistics Engine v{}\nWASM Build\n",
                env!("CARGO_PKG_VERSION")
            )),
            "trajectory" => self.handle_trajectory_command(&args[1..], units),
            "zero" => self.handle_zero_command(&args[1..], units),
            "monte-carlo" | "montecarlo" => self.handle_monte_carlo_command(&args[1..], units),
            "true-velocity" => self.handle_true_velocity_command(&args[1..], units),
            "estimate-bc" => self.handle_estimate_bc_command(&args[1..], units),
            "lead" => self.handle_lead_command(&args[1..], units),
            "powder" => self.handle_powder_command(&args[1..], units),
            _ => Ok(format!(
                "Error: Unknown command '{}'\n\n{}",
                args[0],
                self.show_help()
            )),
        }
    }

    fn handle_trajectory_command(
        &self,
        args: &[&str],
        units: UnitSystem,
    ) -> Result<String, JsValue> {
        // Default values based on unit system
        let (default_velocity, default_mass, default_diameter, default_temp, default_pressure) =
            match units {
                UnitSystem::Imperial => (2700.0, 168.0, 0.308, 59.0, 29.92),
                UnitSystem::Metric => (823.0, 10.9, 7.82, 15.0, 1013.25),
            };

        // Initialize all parameters with defaults
        let mut velocity = default_velocity;
        let mut angle = 0.0;
        let mut bc = 0.475;
        let mut mass = default_mass;
        let mut diameter = default_diameter;
        let mut drag_model = "G1";
        // MBA-1356: whole-curve drag scale for a custom drag table (see cd_scale_range_warning
        // / CD_SCALE_REQUIRES_DRAG_TABLE above). None until parsed; resolved against
        // self.drag_table once the arg-parse loop is done.
        let mut cd_scale: Option<f64> = None;
        let mut max_range = if units == UnitSystem::Imperial {
            1000.0
        } else {
            914.4
        };
        let mut time_step = 0.001;
        let mut wind_speed = 0.0;
        // f64 is required: `wind_direction.to_radians()` below needs a known receiver
        // type (method resolution can't infer it from the field assignment alone).
        let mut wind_direction: f64 = 0.0;
        // Vertical wind (mph imperial / m/s metric); positive = updraft (raises POI).
        let mut wind_vertical: f64 = 0.0;
        // Raw "SPEED:ANGLE:UNTIL" strings; every --wind-segment occurrence is collected
        // (the parse loop visits all args, so repeats accumulate here).
        let mut wind_segment_strs: Vec<String> = Vec::new();
        // Raw "VMIN:VMAX:BC" strings from every --bc-segment occurrence (manual
        // velocity-keyed BC segments; take precedence over a loaded BC5D table).
        let mut bc_segment_strs: Vec<String> = Vec::new();
        // Raw "MACH:DSF" strings from every --dsf-point occurrence (MBA-1411: per-call
        // DSF table — WASM has no saved-profile storage to carry one).
        let mut dsf_point_strs: Vec<String> = Vec::new();
        let mut temperature = default_temp;
        let mut pressure = default_pressure;
        let mut humidity = 50.0;
        let mut altitude = 0.0;
        let mut output_format = OutputFormat::Table;
        let mut full = false;
        // MBA-1337 p3: native --plot parity. None = no chart; Some(style) appends the
        // two terminal charts after the table, exactly like the native CLI.
        let mut plot: Option<crate::terminal_plot::CanvasStyle> = None;
        let mut auto_zero: Option<f64> = None;
        let mut sight_height = if units == UnitSystem::Imperial {
            2.0
        } else {
            50.0
        }; // inches or mm
        let mut muzzle_height = if units == UnitSystem::Imperial {
            60.0
        } else {
            1500.0
        }; // inches or mm (standing)
        let mut target_height = 0.0; // inches or mm (ground level)

        // Advanced physics flags
        let mut enable_magnus = false;
        let mut enable_coriolis = false;
        let mut use_euler = false;
        let mut use_rk4_fixed = false; // Use fixed-step RK4 instead of adaptive RK45
        let mut enable_spin_drift = false;
        let mut enable_wind_shear = false;
        let mut sample_trajectory = false;
        let mut sample_interval = 10.0;
        let mut enable_pitch_damping = false;
        let mut enable_precession = false;
        let mut use_bc_segments = false;
        let mut use_powder_sensitivity = false;
        // Bero-feedback (0.23.0): let a tester keep a full-range trajectory without abusing
        // --muzzle-height to defeat ground truncation (see the muzzle-height sanity warning
        // below — that abuse silently thins the density the whole flight sees).
        let mut ignore_ground_impact = false;
        // Bero-feedback (0.23.0): surface the velocity-keyed BC ladder actually applied
        // (manual --bc-segment or BC5D-synthesized), since the WASM has no stderr the way
        // the native CLI's --print-bc-segments does.
        let mut print_bc_segments = false;
        // Mover ring (MBA-1325): mph imperial / m/s metric, same convention as `lead
        // --target-speed`. 0.0 (default) disables the per-point Ring column/fields.
        let mut target_speed: f64 = 0.0;
        // Angular unit for the ring TABLE column (review fix M3) — native trajectory
        // already exposes --adjustment-unit (PDF dope card), so the WASM surface honors
        // it too. Only the table reads it; CSV keeps ring_mil and JSON keeps
        // mover_ring_m/mover_ring_mil (unit-in-the-name contract).
        let mut adjustment_unit = "mil";
        // MBA-1355: turret click graduations for `--adjustment-unit clicks`. Only
        // elevation is ever consumed (the Ring column, like the native PDF dope card's
        // Drop column, reuses the elevation graduation) — windage is parsed/validated
        // for CLI-flag parity with native trajectory but is otherwise inert here (WASM's
        // trajectory table has no PDF dope card / windage-driven column).
        let mut elevation_click_value: Option<&str> = None;
        let mut windage_click_value: Option<&str> = None;

        // Additional parameters
        let mut twist_rate: Option<f64> = None;
        let mut twist_right = true;
        let mut latitude: Option<f64> = None;
        let mut shot_direction: Option<f64> = None; // compass bearing, degrees, 0=N (Coriolis)
        let mut shooting_angle = 0.0;
        let mut cant_angle_deg = 0.0;
        let mut powder_temp_sensitivity = if units == UnitSystem::Imperial {
            1.0
        } else {
            0.3048 / (5.0 / 9.0)
        };
        let mut powder_temp = if units == UnitSystem::Imperial {
            crate::constants::DEFAULT_POWDER_REFERENCE_TEMP_F
        } else {
            crate::constants::DEFAULT_POWDER_REFERENCE_TEMP_C
        };
        // Optional measured powder-temperature -> velocity curve ("TEMP:VEL,..."),
        // parsed after unit resolution. Supersedes the linear sensitivity model.
        let mut powder_temp_curve_str: Option<String> = None;
        // Track whether --powder-temp was explicitly given. When a curve is present it
        // becomes the powder temperature the curve is interpolated at (decoupled from the
        // air temperature); when not given, the curve falls back to the air temperature.
        let mut powder_temp_provided = false;

        // Zero-day condition overrides. When --auto-zero is used, these let the zero
        // ANGLE be solved under the conditions the rifle was actually zeroed in (a
        // different day: air temperature, pressure, humidity, altitude, and — via the
        // caller's own powder-temp/velocity table — muzzle velocity), while the
        // trajectory itself runs under the current shot-day conditions. Omitting all zero-day
        // flags reuses the shot-day values exactly; coupled powder fallbacks are resolved below.
        let mut zero_velocity: Option<f64> = None;
        let mut zero_temperature: Option<f64> = None;
        let mut zero_pressure: Option<f64> = None;
        let mut zero_humidity: Option<f64> = None;
        let mut zero_altitude: Option<f64> = None;
        let mut zero_powder_temp: Option<f64> = None;

        // Parse arguments
        let mut i = 0;
        while i < args.len() {
            match args[i] {
                "-v" | "--velocity" => {
                    if i + 1 < args.len() {
                        velocity = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid velocity"))?;
                        i += 1;
                    }
                }
                "-a" | "--angle" => {
                    if i + 1 < args.len() {
                        angle = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid angle"))?;
                        i += 1;
                    }
                }
                "-b" | "--bc" => {
                    if i + 1 < args.len() {
                        bc = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid BC"))?;
                        i += 1;
                    }
                }
                "-m" | "--mass" => {
                    if i + 1 < args.len() {
                        mass = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid mass"))?;
                        i += 1;
                    }
                }
                "-d" | "--diameter" => {
                    if i + 1 < args.len() {
                        diameter = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid diameter"))?;
                        i += 1;
                    }
                }
                "--drag-model" => {
                    if i + 1 < args.len() {
                        drag_model = args[i + 1];
                        i += 1;
                    }
                }
                "--cd-scale" => {
                    if i + 1 < args.len() {
                        cd_scale = Some(
                            args[i + 1]
                                .parse()
                                .map_err(|_| JsValue::from_str("Invalid cd-scale"))?,
                        );
                        i += 1;
                    }
                }
                "--max-range" => {
                    if i + 1 < args.len() {
                        max_range = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid max range"))?;
                        i += 1;
                    }
                }
                "--time-step" => {
                    if i + 1 < args.len() {
                        time_step = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid time step"))?;
                        i += 1;
                    }
                }
                "--wind-speed" => {
                    if i + 1 < args.len() {
                        wind_speed = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid wind speed"))?;
                        i += 1;
                    }
                }
                "--wind-direction" => {
                    if i + 1 < args.len() {
                        wind_direction = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid wind direction"))?;
                        i += 1;
                    }
                }
                "--wind-vertical" => {
                    if i + 1 < args.len() {
                        wind_vertical = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid wind vertical"))?;
                        i += 1;
                    }
                }
                "--wind-segment" => {
                    if i + 1 < args.len() {
                        wind_segment_strs.push(args[i + 1].to_string());
                        i += 1;
                    }
                }
                "--bc-segment" => {
                    if i + 1 < args.len() {
                        bc_segment_strs.push(args[i + 1].to_string());
                        i += 1;
                    }
                }
                "--dsf-point" => {
                    if i + 1 < args.len() {
                        dsf_point_strs.push(args[i + 1].to_string());
                        i += 1;
                    }
                }
                "--temperature" => {
                    if i + 1 < args.len() {
                        temperature = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid temperature"))?;
                        i += 1;
                    }
                }
                "--pressure" => {
                    if i + 1 < args.len() {
                        pressure = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid pressure"))?;
                        i += 1;
                    }
                }
                "--humidity" => {
                    if i + 1 < args.len() {
                        humidity = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid humidity"))?;
                        i += 1;
                    }
                }
                "--altitude" => {
                    if i + 1 < args.len() {
                        altitude = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid altitude"))?;
                        i += 1;
                    }
                }
                "-o" | "--output" => {
                    if i + 1 < args.len() {
                        output_format = OutputFormat::from_str(args[i + 1]);
                        i += 1;
                    }
                }
                "--plot" => {
                    // clap parity (num_args 0..=1, default_missing_value "braille"):
                    // bare --plot before another flag = braille; an explicit next
                    // token is the style value and anything but braille/ascii is a
                    // hard error, exactly like the native CLI.
                    plot = Some(crate::terminal_plot::CanvasStyle::Braille);
                    if i + 1 < args.len() && !args[i + 1].starts_with('-') {
                        match args[i + 1] {
                            "braille" => i += 1,
                            "ascii" => {
                                plot = Some(crate::terminal_plot::CanvasStyle::Ascii);
                                i += 1;
                            }
                            other => {
                                return Err(JsValue::from_str(&format!(
                                    "Invalid value '{other}' for --plot (expected braille or ascii)"
                                )));
                            }
                        }
                    }
                }
                "--full" => full = true,
                "--auto-zero" | "-z" => {
                    if i + 1 < args.len() {
                        auto_zero = Some(
                            args[i + 1]
                                .parse()
                                .map_err(|_| JsValue::from_str("Invalid zero distance"))?,
                        );
                        i += 1;
                    }
                }
                "--sight-height" => {
                    if i + 1 < args.len() {
                        sight_height = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid sight height"))?;
                        i += 1;
                    }
                }
                // MBA-1339: --bore-height is the native CLI's name for this same parameter,
                // now on identical inches/mm units — accept both names on both surfaces.
                "--muzzle-height" | "--bore-height" => {
                    if i + 1 < args.len() {
                        muzzle_height = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid muzzle height"))?;
                        i += 1;
                    }
                }
                "--target-height" => {
                    if i + 1 < args.len() {
                        target_height = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid target height"))?;
                        i += 1;
                    }
                }
                "--enable-magnus" => enable_magnus = true,
                "--enable-coriolis" => enable_coriolis = true,
                "--use-euler" => use_euler = true,
                "--use-rk4-fixed" => use_rk4_fixed = true,
                "--enable-spin-drift" => enable_spin_drift = true,
                "--enable-wind-shear" => enable_wind_shear = true,
                "--sample-trajectory" => sample_trajectory = true,
                "--sample-interval" => {
                    if i + 1 < args.len() {
                        sample_interval = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid sample interval"))?;
                        i += 1;
                    }
                }
                "--enable-pitch-damping" => enable_pitch_damping = true,
                "--enable-precession" => enable_precession = true,
                "--use-bc-segments" => use_bc_segments = true,
                "--use-powder-sensitivity" => use_powder_sensitivity = true,
                "--ignore-ground-impact" => ignore_ground_impact = true,
                "--print-bc-segments" => print_bc_segments = true,
                "--target-speed" => {
                    if i + 1 < args.len() {
                        target_speed = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid target speed"))?;
                        // Same bound the native flags enforce via f64_range(0.0, 300.0)
                        // (review fix M1) — reject, never silently clamp.
                        if !(0.0..=300.0).contains(&target_speed) {
                            return Err(JsValue::from_str(&format!(
                                "--target-speed must be between 0 and 300 (mph/m/s), got {}",
                                target_speed
                            )));
                        }
                        i += 1;
                    }
                }
                "--adjustment-unit" => {
                    if i + 1 < args.len() {
                        adjustment_unit = args[i + 1];
                        i += 1;
                    }
                }
                "--elevation-click-value" => {
                    if i + 1 < args.len() {
                        elevation_click_value = Some(args[i + 1]);
                        i += 1;
                    }
                }
                "--windage-click-value" => {
                    if i + 1 < args.len() {
                        windage_click_value = Some(args[i + 1]);
                        i += 1;
                    }
                }
                "--twist-rate" => {
                    if i + 1 < args.len() {
                        twist_rate = Some(
                            args[i + 1]
                                .parse()
                                .map_err(|_| JsValue::from_str("Invalid twist rate"))?,
                        );
                        i += 1;
                    }
                }
                "--twist-right" => {
                    if i + 1 < args.len() {
                        twist_right = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid twist direction"))?;
                        i += 1;
                    }
                }
                "--latitude" => {
                    if i + 1 < args.len() {
                        latitude = Some(
                            args[i + 1]
                                .parse()
                                .map_err(|_| JsValue::from_str("Invalid latitude"))?,
                        );
                        i += 1;
                    }
                }
                "--shot-direction" => {
                    if i + 1 < args.len() {
                        shot_direction = Some(
                            args[i + 1]
                                .parse()
                                .map_err(|_| JsValue::from_str("Invalid shot-direction"))?,
                        );
                        i += 1;
                    }
                }
                "--shooting-angle" => {
                    if i + 1 < args.len() {
                        shooting_angle = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid shooting angle"))?;
                        i += 1;
                    }
                }
                "--cant" | "--cant-angle" => {
                    if i + 1 < args.len() {
                        cant_angle_deg = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid cant angle"))?;
                        i += 1;
                    }
                }
                "--powder-temp-sensitivity" => {
                    if i + 1 < args.len() {
                        powder_temp_sensitivity = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid powder temp sensitivity"))?;
                        i += 1;
                    }
                }
                "--powder-temp" => {
                    if i + 1 < args.len() {
                        powder_temp = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid powder temp"))?;
                        powder_temp_provided = true;
                        i += 1;
                    }
                }
                "--powder-temp-curve" => {
                    if i + 1 < args.len() {
                        powder_temp_curve_str = Some(args[i + 1].to_string());
                        i += 1;
                    }
                }
                "--zero-velocity" => {
                    if i + 1 < args.len() {
                        zero_velocity = Some(
                            args[i + 1]
                                .parse()
                                .map_err(|_| JsValue::from_str("Invalid zero velocity"))?,
                        );
                        i += 1;
                    }
                }
                "--zero-temperature" => {
                    if i + 1 < args.len() {
                        zero_temperature = Some(
                            args[i + 1]
                                .parse()
                                .map_err(|_| JsValue::from_str("Invalid zero temperature"))?,
                        );
                        i += 1;
                    }
                }
                "--zero-pressure" => {
                    if i + 1 < args.len() {
                        zero_pressure = Some(
                            args[i + 1]
                                .parse()
                                .map_err(|_| JsValue::from_str("Invalid zero pressure"))?,
                        );
                        i += 1;
                    }
                }
                "--zero-humidity" => {
                    if i + 1 < args.len() {
                        zero_humidity = Some(
                            args[i + 1]
                                .parse()
                                .map_err(|_| JsValue::from_str("Invalid zero humidity"))?,
                        );
                        i += 1;
                    }
                }
                "--zero-altitude" => {
                    if i + 1 < args.len() {
                        zero_altitude = Some(
                            args[i + 1]
                                .parse()
                                .map_err(|_| JsValue::from_str("Invalid zero altitude"))?,
                        );
                        i += 1;
                    }
                }
                "--zero-powder-temp" => {
                    if i + 1 < args.len() {
                        zero_powder_temp = Some(
                            args[i + 1]
                                .parse()
                                .map_err(|_| JsValue::from_str("Invalid zero powder temp"))?,
                        );
                        i += 1;
                    }
                }
                // --units/-u (+ its value) is consumed globally in run_command, which
                // pre-scans it to set the unit system before dispatch. Skip it here so
                // it isn't rejected as an unknown flag (this is what blocked metric input).
                "--units" | "-u" => {
                    i += 1;
                }
                // Reject unrecognized flags instead of silently ignoring them, so a
                // typo or a flag that isn't wired into this WASM surface is caught
                // immediately rather than looking like a no-op. (The native CLI's clap
                // parser already does this; the hand-rolled WASM parser did not.)
                other if other.starts_with('-') => {
                    return Err(JsValue::from_str(&format!("Unknown flag: {}", other)));
                }
                _ => {}
            }
            i += 1;
        }

        // MBA-1411: validate the per-call DSF table (if any) before doing any solve work —
        // dimensionless (Mach/DSF ratio), so no unit conversion is needed.
        let dsf_table = parse_dsf_points(&dsf_point_strs)?;

        // Build inputs with unit conversions
        let mut inputs = InternalBallisticInputs::default();

        // Convert units to metric (internal representation)
        match units {
            UnitSystem::Imperial => {
                inputs.muzzle_velocity = velocity * 0.3048; // fps to m/s
                inputs.bullet_mass = mass * crate::constants::GRAINS_TO_KG; // grains to kg
                inputs.bullet_diameter = diameter * 0.0254; // inches to meters
                inputs.sight_height = sight_height * 0.0254; // inches to meters
                inputs.muzzle_height = muzzle_height * 0.0254; // inches to meters
                inputs.target_height = target_height * 0.0254; // inches to meters
            }
            UnitSystem::Metric => {
                inputs.muzzle_velocity = velocity; // already m/s
                inputs.bullet_mass = mass * 0.001; // grams to kg
                inputs.bullet_diameter = diameter * 0.001; // mm to meters
                inputs.sight_height = sight_height * 0.001; // mm to meters
                inputs.muzzle_height = muzzle_height * 0.001; // mm to meters
                inputs.target_height = target_height * 0.001; // mm to meters
            }
        }
        // MBA-1135: mass-based length estimate (mirrors CLI/FFI), replacing the mass-blind
        // 4.5-caliber heuristic. WASM otherwise left it at the struct default regardless of the
        // supplied caliber/weight, skewing the Miller Sg / Litz spin drift / Magnus.
        inputs.bullet_length =
            crate::stability::estimate_bullet_length_m(inputs.bullet_diameter, inputs.bullet_mass);
        if inputs.bullet_length <= 0.0 {
            inputs.bullet_length = inputs.bullet_diameter * 4.5;
        }

        inputs.bc_value = bc;
        let drag_model_parsed = DragModel::from_str(drag_model)
            .ok_or_else(|| JsValue::from_str("Invalid drag model"))?;
        inputs.bc_type = drag_model_parsed;
        // Custom drag table (MBA-1328): a table loaded via loadDragTable() is a full physical
        // substitute for the BC + G-model (see calculate_drag_coefficient) — apply it
        // unconditionally when present, mirroring native --drag-table (main.rs
        // load_drag_table_or_exit). Unlike a loaded BC5D table it needs no --use-* gate flag.
        if let Some(table) = self.drag_table.borrow().as_ref() {
            inputs.custom_drag_table = Some(table.clone());
        }
        // MBA-1356: --cd-scale requires a loaded drag table, mirroring the native CLI's
        // --cd-scale/--drag-table pairing requirement. Validate before any solve.
        if let Some(scale) = cd_scale {
            if self.drag_table.borrow().is_none() {
                return Err(JsValue::from_str(CD_SCALE_REQUIRES_DRAG_TABLE));
            }
            inputs.cd_scale = scale;
        }
        let cd_scale_warning = cd_scale.and_then(cd_scale_range_warning);
        inputs.muzzle_angle = angle * std::f64::consts::PI / 180.0; // degrees to radians
        inputs.shooting_angle = shooting_angle * std::f64::consts::PI / 180.0;
        inputs.cant_angle = cant_angle_deg * std::f64::consts::PI / 180.0;
        inputs.ground_threshold = if ignore_ground_impact {
            f64::NEG_INFINITY
        } else {
            0.0
        };

        // Bero-feedback (0.23.0): trajectory altitude feeds air density (thinning it over a
        // long flight), so a tester using --muzzle-height to defeat ground truncation instead
        // of --ignore-ground-impact silently flies through unrealistically thin air. Warn
        // (don't clamp) once the bore sits implausibly high — matches the native CLI threshold.
        let muzzle_height_warning = if inputs.muzzle_height > 1000.0 {
            let muzzle_unit_label = match units {
                UnitSystem::Imperial => "in",
                UnitSystem::Metric => "mm",
            };
            Some(format!(
                "WARNING: --muzzle-height {muzzle_height}{muzzle_unit_label} puts the bore \
                 {:.0} m above ground. Trajectory altitude feeds air density, so this thins \
                 the air over the whole flight (a 25 km 'muzzle height' flies in ~2% density). \
                 If you meant site elevation use --altitude; to defeat ground truncation use \
                 --ignore-ground-impact.\n\n",
                inputs.muzzle_height
            ))
        } else {
            None
        };

        // Set advanced physics flags. enable_advanced_effects remains the umbrella
        // flag, but Magnus and Coriolis are now gated independently so enabling one
        // does not silently enable the other.
        if enable_magnus || enable_coriolis {
            inputs.enable_advanced_effects = true;
        }
        inputs.enable_magnus = enable_magnus;
        inputs.enable_coriolis = enable_coriolis;
        // Set integration method: Euler < RK4 fixed < RK45 adaptive (default)
        if use_euler {
            inputs.use_rk4 = false;
            inputs.use_adaptive_rk45 = false;
        } else if use_rk4_fixed {
            inputs.use_rk4 = true;
            inputs.use_adaptive_rk45 = false; // Fixed-step RK4
        } else {
            inputs.use_rk4 = true;
            inputs.use_adaptive_rk45 = true; // Default: adaptive RK45
        }
        inputs.use_enhanced_spin_drift = enable_spin_drift;
        inputs.enable_wind_shear = enable_wind_shear;
        inputs.enable_trajectory_sampling = sample_trajectory;
        inputs.sample_interval = sample_interval;
        inputs.enable_pitch_damping = enable_pitch_damping;
        inputs.enable_precession_nutation = enable_precession;
        inputs.use_bc_segments = use_bc_segments;
        inputs.use_powder_sensitivity = use_powder_sensitivity;

        // Velocity-keyed BC segments, in priority order:
        //   1. manual --bc-segment "VMIN:VMAX:BC" pairs (explicit user input)
        //   2. a BC5D table loaded via loadBc5dTable() + --use-bc-segments
        // Velocities are in the command's display units (fps imperial, m/s metric);
        // the solver compares against fps, so convert.
        let vel_to_fps = match units {
            UnitSystem::Imperial => 1.0,
            UnitSystem::Metric => 3.280_839_895,
        };
        let mut manual_bc_segments: Vec<crate::BCSegmentData> = Vec::new();
        for s in &bc_segment_strs {
            let parts: Vec<&str> = s.split(':').collect();
            if parts.len() != 3 {
                return Err(JsValue::from_str(&format!(
                    "--bc-segment expects VMIN:VMAX:BC (e.g. 1500:1800:0.243), got '{}'",
                    s
                )));
            }
            let vmin: f64 = parts[0].trim().parse().map_err(|_| {
                JsValue::from_str(&format!("--bc-segment: invalid VMIN in '{}'", s))
            })?;
            let vmax: f64 = parts[1].trim().parse().map_err(|_| {
                JsValue::from_str(&format!("--bc-segment: invalid VMAX in '{}'", s))
            })?;
            let bcv: f64 = parts[2].trim().parse().map_err(|_| {
                JsValue::from_str(&format!("--bc-segment: invalid BC in '{}'", s))
            })?;
            if !(vmin < vmax) {
                return Err(JsValue::from_str(&format!(
                    "--bc-segment: VMIN must be < VMAX in '{}'",
                    s
                )));
            }
            if bcv <= 0.0 {
                return Err(JsValue::from_str(&format!(
                    "--bc-segment: BC must be > 0 in '{}'",
                    s
                )));
            }
            manual_bc_segments.push(crate::BCSegmentData {
                velocity_min: vmin * vel_to_fps,
                velocity_max: vmax * vel_to_fps,
                bc_value: bcv,
            });
        }

        let mut bc5d_coercion_warning: Option<String> = None;
        if !manual_bc_segments.is_empty() {
            // Manual segments win; imply --use-bc-segments so they're applied.
            inputs.bc_segments_data = Some(manual_bc_segments);
            inputs.use_bc_segments = true;
        } else if use_bc_segments {
            // BC5D offline correction: synthesize velocity-dependent BC segments
            // from a loaded table. The table's native units are grains + fps.
            if let Some(table) = self.bc5d_table.borrow().as_ref() {
                let (weight_grains, muzzle_fps) = match units {
                    UnitSystem::Imperial => (mass, velocity),
                    UnitSystem::Metric => (mass * GRAINS_PER_GRAM, velocity * 3.280839895),
                };
                // MBA-1386: generate_segment_schedule types anything non-G7 as G1
                // (bc_table_5d.rs); surface the coercion like the native aux paths do.
                if !drag_model.eq_ignore_ascii_case("g1") && !drag_model.eq_ignore_ascii_case("g7")
                {
                    bc5d_coercion_warning = Some(format!(
                        "warning: BC5D correction tables support G1/G7 only; treating drag \
                         model '{}' as G1\n\n",
                        drag_model.to_uppercase()
                    ));
                }
                if let Some(schedule) =
                    table.generate_segment_schedule(bc, drag_model, weight_grains, muzzle_fps)
                {
                    inputs.bc_segments_data = Some(schedule.segments);
                    inputs.bc_value = schedule.fallback_bc;
                }
            }
        }

        // Set additional parameters
        let explicit_twist_inches = twist_rate.map(|rate| match units {
            UnitSystem::Imperial => rate,
            UnitSystem::Metric => rate / 25.4,
        });
        inputs.twist_rate = crate::stability::resolve_twist_inches(
            explicit_twist_inches,
            inputs.bullet_diameter,
            inputs.bullet_mass,
            inputs.muzzle_velocity,
        );
        inputs.is_twist_right = twist_right;
        if let Some(lat) = latitude {
            inputs.latitude = Some(lat);
        }
        inputs.shot_azimuth = shot_direction.map(|d| d.to_radians()).unwrap_or(0.0);
        inputs.powder_temp_sensitivity = match units {
            UnitSystem::Imperial => powder_temp_sensitivity * 0.3048 / (5.0 / 9.0),
            UnitSystem::Metric => powder_temp_sensitivity,
        };
        inputs.powder_temp = match units {
            UnitSystem::Imperial => (powder_temp - 32.0) * 5.0 / 9.0,
            UnitSystem::Metric => powder_temp,
        };
        // When --powder-temp was explicitly given, it becomes the POWDER temperature the
        // curve is interpolated at (decoupled from --temperature / air density). When not
        // given, powder_curve_temp_c stays None so the curve falls back to the air temp.
        inputs.powder_curve_temp_c = if powder_temp_provided {
            Some(inputs.powder_temp)
        } else {
            None
        };
        // Parse the optional powder-temperature -> velocity curve into SI points
        // (shared parser — see parse_powder_temp_curve_str).
        if let Some(curve_str) = &powder_temp_curve_str {
            inputs.powder_temp_curve = Some(parse_powder_temp_curve_str(curve_str, units)?);
        }

        // Set wind conditions
        let mut wind = WindConditions::default();
        match units {
            UnitSystem::Imperial => {
                wind.speed = wind_speed * 0.44704; // mph to m/s
                wind.vertical_speed = wind_vertical * 0.44704; // mph to m/s
            }
            UnitSystem::Metric => {
                wind.speed = wind_speed; // already m/s
                wind.vertical_speed = wind_vertical; // already m/s
            }
        }
        // WindConditions.direction is RADIANS (0=North, PI/2=East); --wind-direction is degrees.
        // Convert (matches native CLI); previously a 90-degree crosswind was fed as 90 radians.
        wind.direction = wind_direction.to_radians();

        // Set atmospheric conditions
        let mut atmosphere = AtmosphericConditions::default();
        match units {
            UnitSystem::Imperial => {
                atmosphere.temperature = (temperature - 32.0) * 5.0 / 9.0; // F to C
                atmosphere.pressure = pressure * 33.863886666667; // inHg to hPa
                atmosphere.altitude = altitude * 0.3048; // feet to meters
            }
            UnitSystem::Metric => {
                atmosphere.temperature = temperature;
                atmosphere.pressure = pressure;
                atmosphere.altitude = altitude;
            }
        }
        atmosphere.humidity = humidity;
        inputs.temperature = atmosphere.temperature;
        inputs.pressure = atmosphere.pressure;
        inputs.humidity = (humidity / 100.0).clamp(0.0, 1.0);
        inputs.altitude = atmosphere.altitude;

        // Handle auto-zero if specified
        let mut zero_info = String::new();
        if let Some(zero_distance) = auto_zero {
            let zero_distance_m = match units {
                UnitSystem::Imperial => zero_distance * 0.9144, // yards to meters
                UnitSystem::Metric => zero_distance,
            };

            // Build the condition set the zero ANGLE is solved under. It starts from the
            // shot-day inputs/atmosphere and is overridden only by whichever --zero-*
            // flags the caller supplied, so a rifle zeroed on a different day (different
            // air density and/or muzzle velocity) is modeled correctly while the
            // trajectory below still runs under the current shot-day conditions.
            let mut zero_inputs = inputs.clone();
            // The zero is torn on a LEVEL range: shot-day slope and cant belong to
            // the SHOT, not the zero geometry. The native CLI's zero_inputs literal
            // never carries them (both default to 0); this clone must strip them or
            // an inclined shot (--shooting-angle) makes the zero root-find
            // unbracketable — and where it does converge, it bakes the incline into
            // the zero, which is the wrong physics (Bero's PRS report).
            zero_inputs.shooting_angle = 0.0;
            zero_inputs.cant_angle = 0.0;
            // Fix-half of MBA-1384: the native CLI zero solve never carries
            // Coriolis — its zero_inputs literal ends in Default::default()
            // (enable_coriolis=false, latitude=None, shot_azimuth=0), while this
            // clone runs after the shot's Coriolis fields are set. Strip them so
            // the terminal and the CLI agree on the rifle zero. Zero-day opt-in
            // flags (--zero-latitude/--zero-azimuth) are the deferred feature
            // half of MBA-1384.
            zero_inputs.enable_coriolis = false;
            zero_inputs.latitude = None;
            zero_inputs.shot_azimuth = 0.0;
            let mut zero_atmosphere = atmosphere.clone();
            let mut zero_day_overridden = false;
            if let Some(zv) = zero_velocity {
                zero_inputs.muzzle_velocity = match units {
                    UnitSystem::Imperial => zv * 0.3048, // fps to m/s
                    UnitSystem::Metric => zv,
                };
                // An explicit zero-day velocity takes precedence: disable both velocity
                // adjustment models for the zero solve so neither changes the supplied value.
                // (zero_inputs is a clone of inputs, which may carry the shot-day models.)
                zero_inputs.powder_temp_curve = None;
                zero_inputs.use_powder_sensitivity = false;
                zero_day_overridden = true;
            }
            if let Some(zt) = zero_temperature {
                let t_c = match units {
                    UnitSystem::Imperial => (zt - 32.0) * 5.0 / 9.0, // F to C
                    UnitSystem::Metric => zt,
                };
                zero_atmosphere.temperature = t_c;
                zero_inputs.temperature = t_c;
                zero_day_overridden = true;
            }
            if let Some(zp) = zero_pressure {
                let p_hpa = match units {
                    UnitSystem::Imperial => zp * 33.863886666667, // inHg to hPa
                    UnitSystem::Metric => zp,
                };
                zero_atmosphere.pressure = p_hpa;
                zero_inputs.pressure = p_hpa;
                zero_day_overridden = true;
            }
            if let Some(zh) = zero_humidity {
                zero_atmosphere.humidity = zh;
                zero_inputs.humidity = (zh / 100.0).clamp(0.0, 1.0);
                zero_day_overridden = true;
            }
            if let Some(za) = zero_altitude {
                let a_m = match units {
                    UnitSystem::Imperial => za * 0.3048, // feet to meters
                    UnitSystem::Metric => za,
                };
                zero_atmosphere.altitude = a_m;
                zero_inputs.altitude = a_m;
                zero_day_overridden = true;
            }
            // An explicit zero-day powder temperature wins. Otherwise an explicit zero-day air
            // temperature retains the established "powder follows zero-day air" behavior. With
            // neither override, keep the cloned shot-day powder temperature so a no-override
            // zero solve reproduces the flight conditions exactly.
            if let Some(zpt) = zero_powder_temp {
                zero_inputs.powder_curve_temp_c = Some(match units {
                    UnitSystem::Imperial => (zpt - 32.0) * 5.0 / 9.0,
                    UnitSystem::Metric => zpt,
                });
                zero_day_overridden = true;
            } else if zero_temperature.is_some() {
                zero_inputs.powder_curve_temp_c = None;
            }

            match calculate_zero_angle_with_conditions(
                zero_inputs.clone(),
                zero_distance_m,
                zero_inputs.muzzle_height + zero_inputs.sight_height, // Zero crosses the line of sight (matches CLI)
                wind.clone(),
                zero_atmosphere.clone(),
            ) {
                Ok(zero_angle) => {
                    inputs.muzzle_angle = zero_angle;
                    let moa_adjustment = zero_angle * 180.0 / std::f64::consts::PI * 60.0;
                    let mrad_adjustment = zero_angle * 1000.0;
                    zero_info = format!(
                        "Rifle zeroed at {} {}{} (Adjustment: {:.2} MOA / {:.2} mrad up)\n\n",
                        zero_distance,
                        if units == UnitSystem::Imperial {
                            "yards"
                        } else {
                            "meters"
                        },
                        if zero_day_overridden {
                            " under supplied zero-day conditions"
                        } else {
                            ""
                        },
                        moa_adjustment,
                        mrad_adjustment
                    );
                }
                Err(e) => {
                    return Ok(format!("Error calculating zero: {}\n\nTry a shorter zero distance or check your ballistic parameters.", e));
                }
            }
        }

        // Create solver and calculate
        let mut solver = TrajectorySolver::new(inputs.clone(), wind, atmosphere);
        let max_range_m = match units {
            UnitSystem::Imperial => max_range * 0.9144, // yards to meters
            UnitSystem::Metric => max_range,
        };
        solver.set_max_range(max_range_m);
        solver.set_time_step(time_step);

        // Downrange-segmented wind overrides the scalar wind when present.
        if !wind_segment_strs.is_empty() {
            if enable_wind_shear {
                return Err(JsValue::from_str(
                    "--wind-segment cannot be combined with --enable-wind-shear \
                     (downrange segments + altitude shear is not yet a defined model)",
                ));
            }
            let imperial = matches!(units, UnitSystem::Imperial);
            let mut segments = Vec::with_capacity(wind_segment_strs.len());
            for s in &wind_segment_strs {
                let seg = crate::wind::parse_wind_segment_str(s, imperial)
                    .map_err(|err| JsValue::from_str(&err))?;
                segments.push(seg);
            }
            solver.set_wind_segments(segments);
        }

        // Mover ring (MBA-1325): --target-speed uses the same mph (imperial) / m/s
        // (metric) convention as --wind-speed / `lead --target-speed`. 0 disables the
        // per-point Ring column/field in every output format below (additive).
        let target_speed_mps = match units {
            UnitSystem::Imperial => target_speed * 0.44704, // mph to m/s
            UnitSystem::Metric => target_speed,
        };
        // Validate --adjustment-unit like handle_lead_command does; only the ring
        // table column consumes it on this command. MBA-1355: smoa/iphy join
        // mil/moa as constant-factor units; clicks needs a resolved elevation click
        // graduation (--elevation-click-value or nothing else — WASM has no --profile),
        // mirroring native's resolve_click_values/reject_clicks_out_of_scope split
        // (trajectory is one of the two commands where clicks actually resolves).
        let adjustment_unit_lower = adjustment_unit.to_lowercase();
        let ring_unit = match adjustment_unit_lower.as_str() {
            "mil" => RingDisplayUnit::Factor(1.0, "mil"),
            "moa" => RingDisplayUnit::Factor(
                crate::moving_target::MOA_PER_UNIT_RATIO / crate::moving_target::MIL_PER_UNIT_RATIO,
                "moa",
            ),
            // SMOA/IPHY are numerically identical (exact inches-per-hundred-yards); only
            // the header/cell text differs, so they get their own labels but share the
            // ratio, matching native's Ring(smoa)/Ring(iphy) treatment.
            "smoa" => RingDisplayUnit::Factor(smoa_per_mil(), "smoa"),
            "iphy" => RingDisplayUnit::Factor(smoa_per_mil(), "iphy"),
            "clicks" => {
                let elev_str = elevation_click_value.ok_or_else(|| {
                    JsValue::from_str(
                        "--adjustment-unit clicks requires a turret elevation graduation: pass \
                         --elevation-click-value <SIZE><UNIT> (e.g. 0.25moa or 0.1mil)",
                    )
                })?;
                let elevation = crate::adjustment::parse_click_value(elev_str)
                    .map_err(|e| JsValue::from_str(&e))?;
                // Windage is validated for CLI-flag parity with native (a typo'd
                // --windage-click-value still errors here, same as native trajectory)
                // but is otherwise inert — the Ring column always uses elevation.
                if let Some(w) = windage_click_value {
                    crate::adjustment::parse_click_value(w).map_err(|e| JsValue::from_str(&e))?;
                }
                RingDisplayUnit::Clicks(elevation)
            }
            _ => {
                return Err(JsValue::from_str(&format!(
                    "Invalid --adjustment-unit '{}' (expected mil, moa, smoa, iphy, or clicks)",
                    adjustment_unit
                )));
            }
        };

        match solver.solve() {
            Ok(mut result) => {
                // MBA-1411: apply the per-call DSF table (if any), IN PLACE, before any
                // output format below reads `result` — mirrors native main.rs's single
                // post-solve `apply_dsf` site (drop-only: velocity/energy/TOF/windage and
                // sampled_points' non-drop fields are untouched; see apply_dsf's doc
                // comment for the invariant). Must run identically ahead of every format,
                // including the Ring column (mover_ring, below), which keys off
                // `result.points`' downrange position/time only — untouched by apply_dsf.
                if let Some(table) = dsf_table.as_ref() {
                    apply_dsf(&mut result, table);
                }
                let output = match output_format {
                    OutputFormat::Table => self.format_trajectory_table(
                        &result,
                        auto_zero,
                        units,
                        full,
                        inputs.muzzle_height + inputs.sight_height,
                        target_speed_mps,
                        ring_unit,
                    ),
                    OutputFormat::Json => self.format_trajectory_json(
                        &result,
                        units,
                        inputs.muzzle_height + inputs.sight_height,
                        target_speed_mps,
                    ),
                    OutputFormat::Csv => self.format_trajectory_csv(
                        &result,
                        units,
                        full,
                        inputs.muzzle_height + inputs.sight_height,
                        target_speed_mps,
                    ),
                };
                // JSON/CSV must stay pure machine output: the human-readable "Rifle zeroed
                // at ..." banner is table-only, exactly like the bc-segments and warning
                // blocks below. Prepending it to a JSON/CSV payload makes it unparseable.
                let mut combined = if matches!(output_format, OutputFormat::Table) {
                    // MBA-1386/MBA-1356: table-only, like every human-readable block here —
                    // neither warning may contaminate JSON/CSV payloads.
                    format!(
                        "{}{}{}{}",
                        cd_scale_warning.as_deref().unwrap_or(""),
                        bc5d_coercion_warning.as_deref().unwrap_or(""),
                        zero_info,
                        output
                    )
                } else {
                    output
                };
                // Human-readable append only for the table view: tacking text after a
                // JSON/CSV payload would break any downstream parser of those formats.
                if print_bc_segments && matches!(output_format, OutputFormat::Table) {
                    combined.push_str(&self.format_bc_segments_report(&inputs, units));
                }
                // MBA-1411: table-output-only note that a per-call DSF table was applied —
                // the Drop column above already reflects it (`apply_dsf`, called before
                // this match arm's formatting, corrects every output format identically).
                // JSON/CSV get no equivalent text or top-level field (purity rule, matching
                // native main.rs's saved-profile note).
                if let Some(table) = dsf_table.as_ref() {
                    if matches!(output_format, OutputFormat::Table) {
                        combined.push_str(&format!(
                            "\nDSF table active ({})\n",
                            dsf_table_summary(table.points())
                        ));
                    }
                }
                // MBA-1337 p3: table-only chart append, mirroring the native CLI's
                // --plot block (72x12 cells; drop then lateral drift vs range). JSON
                // and CSV stay pure machine output.
                if let Some(style) = plot {
                    if matches!(output_format, OutputFormat::Table) && !result.points.is_empty() {
                        let (dist_div, range_unit) = match units {
                            UnitSystem::Imperial => (0.9144, "yd"),
                            UnitSystem::Metric => (1.0, "m"),
                        };
                        let drop_label = format!("drop ({})", range_unit);
                        let drop_points: Vec<(f64, f64)> = result
                            .points
                            .iter()
                            .map(|p| (p.position.x / dist_div, p.position.y / dist_div))
                            .collect();
                        combined.push_str("\nDrop vs Range:\n");
                        combined.push_str(&crate::terminal_plot::render_chart(
                            &[(drop_label.as_str(), drop_points.as_slice())],
                            72,
                            12,
                            style,
                        ));
                        combined.push('\n');

                        let drift_label = format!("drift ({})", range_unit);
                        let drift_points: Vec<(f64, f64)> = result
                            .points
                            .iter()
                            .map(|p| (p.position.x / dist_div, p.position.z / dist_div))
                            .collect();
                        combined.push_str("\nLateral Drift vs Range:\n");
                        combined.push_str(&crate::terminal_plot::render_chart(
                            &[(drift_label.as_str(), drift_points.as_slice())],
                            72,
                            12,
                            style,
                        ));
                        combined.push('\n');
                    }
                }
                if let Some(warning) = &muzzle_height_warning {
                    if matches!(output_format, OutputFormat::Table) {
                        combined = format!("{}{}", warning, combined);
                    }
                }
                Ok(combined)
            }
            Err(e) => {
                let mut combined = format!("Error: {}", e);
                // MBA-1386 (0.28.1 sweep): the Ok arm above prepends bc5d_coercion_warning
                // (table-only, like every other human-readable block here); the Err arm was
                // dropping it entirely instead of carrying it the same way.
                if let Some(warning) = &bc5d_coercion_warning {
                    if matches!(output_format, OutputFormat::Table) {
                        combined = format!("{}{}", warning, combined);
                    }
                }
                if let Some(warning) = &muzzle_height_warning {
                    combined = format!("{}{}", warning, combined);
                }
                Ok(combined)
            }
        }
    }

    /// Bero-feedback (0.23.0) `--print-bc-segments`: report the velocity-keyed BC ladder
    /// actually applied to this run (manual `--bc-segment` entries or a BC5D-synthesized
    /// schedule), or a one-line note when none are active. `BCSegmentData` velocities are
    /// always stored in fps (see the parse-time conversion in `handle_trajectory_command`),
    /// so this converts to the command's display units and derives a standard-atmosphere
    /// Mach span alongside them.
    fn format_bc_segments_report(
        &self,
        inputs: &InternalBallisticInputs,
        units: UnitSystem,
    ) -> String {
        match inputs.bc_segments_data.as_ref().filter(|s| !s.is_empty()) {
            Some(segments) => {
                let (fps_to_display, unit_label) = match units {
                    UnitSystem::Imperial => (1.0, "fps"),
                    UnitSystem::Metric => (0.3048, "m/s"),
                };
                let mut block = String::from("\nBC Segments (active)\n=====================\n");
                for seg in segments {
                    let mach_min =
                        seg.velocity_min * 0.3048 / crate::constants::SPEED_OF_SOUND_MPS;
                    let mach_max =
                        seg.velocity_max * 0.3048 / crate::constants::SPEED_OF_SOUND_MPS;
                    block.push_str(&format!(
                        "  {:.1}-{:.1} {} (Mach {:.2}-{:.2}): BC {:.5}\n",
                        seg.velocity_min * fps_to_display,
                        seg.velocity_max * fps_to_display,
                        unit_label,
                        mach_min,
                        mach_max,
                        seg.bc_value,
                    ));
                }
                block
            }
            None => "\nNo BC segments active for this run.\n".to_string(),
        }
    }

    fn handle_zero_command(&self, args: &[&str], units: UnitSystem) -> Result<String, JsValue> {
        // Default values
        let (default_velocity, default_mass, default_diameter) = match units {
            UnitSystem::Imperial => (2700.0, 168.0, 0.308),
            UnitSystem::Metric => (823.0, 10.9, 7.82),
        };

        let mut velocity = default_velocity;
        let mut bc = 0.475;
        let mut mass = default_mass;
        let mut diameter = default_diameter;
        let mut target_distance = if units == UnitSystem::Imperial {
            100.0
        } else {
            100.0
        };
        let mut sight_height = if units == UnitSystem::Imperial {
            2.0
        } else {
            50.0
        };
        // Heights above GROUND (--muzzle-height / --target-height) do NOT change the zero ANGLE —
        // for a same-elevation target they cancel — so they are intentionally not parsed here
        // (silently ignored by the catch-all arm below). The SIGHT height IS honored: the zero
        // targets the line-of-sight height at the zero distance (see the calculate_zero call).
        let mut drag_model = "G1";
        // MBA-1356: whole-curve drag scale for a custom drag table. None until parsed;
        // resolved against self.drag_table once the arg-parse loop is done.
        let mut cd_scale: Option<f64> = None;

        // Parse arguments
        let mut i = 0;
        while i < args.len() {
            match args[i] {
                "-v" | "--velocity" => {
                    if i + 1 < args.len() {
                        velocity = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid velocity"))?;
                        i += 1;
                    }
                }
                "-b" | "--bc" => {
                    if i + 1 < args.len() {
                        bc = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid BC"))?;
                        i += 1;
                    }
                }
                "-m" | "--mass" => {
                    if i + 1 < args.len() {
                        mass = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid mass"))?;
                        i += 1;
                    }
                }
                "-d" | "--diameter" => {
                    if i + 1 < args.len() {
                        diameter = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid diameter"))?;
                        i += 1;
                    }
                }
                "--target-distance" => {
                    if i + 1 < args.len() {
                        target_distance = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid target distance"))?;
                        i += 1;
                    }
                }
                "--sight-height" => {
                    if i + 1 < args.len() {
                        sight_height = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid sight height"))?;
                        i += 1;
                    }
                }
                "--drag-model" => {
                    if i + 1 < args.len() {
                        drag_model = args[i + 1];
                        i += 1;
                    }
                }
                "--cd-scale" => {
                    if i + 1 < args.len() {
                        cd_scale = Some(
                            args[i + 1]
                                .parse()
                                .map_err(|_| JsValue::from_str("Invalid cd-scale"))?,
                        );
                        i += 1;
                    }
                }
                // --units/-u (+ its value) is consumed globally in run_command, which
                // pre-scans it to set the unit system before dispatch. Skip it here so
                // it isn't rejected as an unknown flag (this is what blocked metric input).
                "--units" | "-u" => {
                    i += 1;
                }
                // Reject unrecognized flags instead of silently ignoring them, so a
                // typo or a flag that isn't wired into this WASM surface is caught
                // immediately rather than looking like a no-op. (The native CLI's clap
                // parser already does this; the hand-rolled WASM parser did not.)
                other if other.starts_with('-') => {
                    return Err(JsValue::from_str(&format!("Unknown flag: {}", other)));
                }
                _ => {}
            }
            i += 1;
        }

        // Build inputs
        let mut inputs = InternalBallisticInputs::default();

        // Convert units
        match units {
            UnitSystem::Imperial => {
                inputs.muzzle_velocity = velocity * 0.3048;
                inputs.bullet_mass = mass * crate::constants::GRAINS_TO_KG;
                inputs.bullet_diameter = diameter * 0.0254;
                inputs.sight_height = sight_height * 0.0254;
            }
            UnitSystem::Metric => {
                inputs.muzzle_velocity = velocity;
                inputs.bullet_mass = mass * 0.001;
                inputs.bullet_diameter = diameter * 0.001;
                inputs.sight_height = sight_height * 0.001;
            }
        }
        // MBA-1135: mass-based length estimate (mirrors CLI/FFI); replaces the mass-blind
        // 4.5-caliber heuristic. WASM otherwise left it at the struct default.
        inputs.bullet_length =
            crate::stability::estimate_bullet_length_m(inputs.bullet_diameter, inputs.bullet_mass);
        if inputs.bullet_length <= 0.0 {
            inputs.bullet_length = inputs.bullet_diameter * 4.5;
        }

        inputs.bc_value = bc;
        let drag_model_parsed = DragModel::from_str(drag_model)
            .ok_or_else(|| JsValue::from_str("Invalid drag model"))?;
        inputs.bc_type = drag_model_parsed;
        // Custom drag table (MBA-1328): see handle_trajectory_command for the rationale —
        // applied unconditionally, no gate flag, mirrors native --drag-table on `zero`.
        if let Some(table) = self.drag_table.borrow().as_ref() {
            inputs.custom_drag_table = Some(table.clone());
        }
        // MBA-1356: --cd-scale requires a loaded drag table, mirroring the native CLI's
        // --cd-scale/--drag-table pairing requirement. Validate before any solve.
        if let Some(scale) = cd_scale {
            if self.drag_table.borrow().is_none() {
                return Err(JsValue::from_str(CD_SCALE_REQUIRES_DRAG_TABLE));
            }
            inputs.cd_scale = scale;
        }
        let cd_scale_warning = cd_scale.and_then(cd_scale_range_warning);

        let target_distance_m = match units {
            UnitSystem::Imperial => target_distance * 0.9144,
            UnitSystem::Metric => target_distance,
        };

        // MBA-951: target the line-of-sight height at the zero distance (= sight_height), matching
        // the CLI convention in every zero call. Previously 0.0, which solved a BORE-line zero and
        // ignored sight height entirely (off by the sight-height angle — ~2 MOA at 100 yd).
        let los_height = inputs.sight_height;
        let result = match calculate_zero_angle_with_conditions(
            inputs,
            target_distance_m,
            los_height,
            WindConditions::default(),
            AtmosphericConditions::default(),
        ) {
            Ok(zero_angle) => {
                let zero_degrees = zero_angle * 180.0 / std::f64::consts::PI;
                let moa_adjustment = zero_degrees * 60.0;
                let mrad_adjustment = zero_angle * 1000.0;

                format!(
                    "Zero Calculation Results\n\
                     ========================\n\
                     Target Distance: {} {}\n\
                     Zero Angle: {:.4}°\n\
                     MOA Adjustment: {:.2} MOA up\n\
                     Mrad Adjustment: {:.2} mrad up\n\
                     Sight Height: {} {}\n",
                    target_distance,
                    if units == UnitSystem::Imperial {
                        "yards"
                    } else {
                        "meters"
                    },
                    zero_degrees,
                    moa_adjustment,
                    mrad_adjustment,
                    sight_height,
                    if units == UnitSystem::Imperial {
                        "inches"
                    } else {
                        "mm"
                    }
                )
            }
            Err(e) => format!("Error calculating zero: {}", e),
        };
        // MBA-1356: table-only prepend, same pattern as the trajectory/monte-carlo handlers —
        // this command has no JSON/CSV output mode to protect from contamination.
        Ok(format!("{}{}", cd_scale_warning.as_deref().unwrap_or(""), result))
    }

    fn handle_lead_command(&self, args: &[&str], units: UnitSystem) -> Result<String, JsValue> {
        // Default values
        let (default_velocity, default_mass, default_diameter) = match units {
            UnitSystem::Imperial => (2700.0, 168.0, 0.308),
            UnitSystem::Metric => (823.0, 10.9, 7.82),
        };

        let mut velocity = default_velocity;
        let mut bc = 0.475;
        let mut mass = default_mass;
        let mut diameter = default_diameter;
        let mut drag_model = "G1";
        let mut sight_height = if units == UnitSystem::Imperial {
            2.0
        } else {
            50.0
        };
        let mut target_speed: Option<f64> = None;
        let mut target_angle = 90.0;
        let mut range = 500.0;
        let mut adjustment_unit = "mil";
        let mut lead_output = "table";
        // MBA-1411 (carried from the MBA-1356 review's "WASM lead untrued-curve gap"):
        // this command already applies a loaded custom drag table unconditionally (see
        // the self.drag_table block below) but, unlike trajectory/zero/monte-carlo, never
        // wired up --cd-scale to go with it — a table trued via `trajectory --cd-scale`
        // silently reverted to its untrued curve here. None until parsed; resolved against
        // self.drag_table once the arg-parse loop is done, exactly like handle_trajectory_command.
        let mut cd_scale: Option<f64> = None;
        // Powder temperature plumbing (MBA-1325), identical defaults/parsing to
        // handle_trajectory_command so a curve/sensitivity correction resolves the
        // same muzzle velocity from either command.
        let mut use_powder_sensitivity = false;
        let mut powder_temp_sensitivity = if units == UnitSystem::Imperial {
            1.0
        } else {
            0.3048 / (5.0 / 9.0)
        };
        let mut powder_temp = if units == UnitSystem::Imperial {
            crate::constants::DEFAULT_POWDER_REFERENCE_TEMP_F
        } else {
            crate::constants::DEFAULT_POWDER_REFERENCE_TEMP_C
        };
        let mut powder_temp_curve_str: Option<String> = None;
        // Track whether --powder-temp was explicitly given (mirrors handle_trajectory_command):
        // when a curve is present it becomes the powder temperature the curve is interpolated
        // at (decoupled from the air temperature); when not given, the curve falls back to it.
        let mut powder_temp_provided = false;
        // Environmental conditions (MBA-1325 env-flags addendum): native-lead parity for
        // --temperature/--pressure/--humidity/--altitude/--wind-speed/--wind-direction.
        // Parsed as Option — an ABSENT flag leaves the pre-existing default objects
        // (WindConditions::default(), AtmosphericConditions::default(), BallisticInputs::
        // default()'s env fields) completely untouched, so a no-env-flag run stays
        // byte-identical to the pre-flag build by construction. (A display-unit default
        // like trajectory's 29.92 inHg would NOT: 29.92 x 33.863886666667 = 1013.2075 hPa,
        // whereas this command has always solved at AtmosphericConditions::default()'s
        // 1013.25.) Explicit values convert exactly like handle_trajectory_command's.
        let mut temperature: Option<f64> = None;
        let mut pressure: Option<f64> = None;
        let mut humidity: Option<f64> = None;
        let mut altitude: Option<f64> = None;
        let mut wind_speed: Option<f64> = None;
        let mut wind_direction: Option<f64> = None;

        // Parse arguments
        let mut i = 0;
        while i < args.len() {
            match args[i] {
                "-v" | "--velocity" => {
                    if i + 1 < args.len() {
                        velocity = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid velocity"))?;
                        i += 1;
                    }
                }
                "-b" | "--bc" => {
                    if i + 1 < args.len() {
                        bc = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid BC"))?;
                        i += 1;
                    }
                }
                "-m" | "--mass" => {
                    if i + 1 < args.len() {
                        mass = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid mass"))?;
                        i += 1;
                    }
                }
                "-d" | "--diameter" => {
                    if i + 1 < args.len() {
                        diameter = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid diameter"))?;
                        i += 1;
                    }
                }
                "--drag-model" => {
                    if i + 1 < args.len() {
                        drag_model = args[i + 1];
                        i += 1;
                    }
                }
                "--cd-scale" => {
                    if i + 1 < args.len() {
                        cd_scale = Some(
                            args[i + 1]
                                .parse()
                                .map_err(|_| JsValue::from_str("Invalid cd-scale"))?,
                        );
                        i += 1;
                    }
                }
                "--sight-height" => {
                    if i + 1 < args.len() {
                        sight_height = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid sight height"))?;
                        i += 1;
                    }
                }
                "--temperature" => {
                    if i + 1 < args.len() {
                        temperature = Some(
                            args[i + 1]
                                .parse()
                                .map_err(|_| JsValue::from_str("Invalid temperature"))?,
                        );
                        i += 1;
                    }
                }
                "--pressure" => {
                    if i + 1 < args.len() {
                        pressure = Some(
                            args[i + 1]
                                .parse()
                                .map_err(|_| JsValue::from_str("Invalid pressure"))?,
                        );
                        i += 1;
                    }
                }
                "--humidity" => {
                    if i + 1 < args.len() {
                        humidity = Some(
                            args[i + 1]
                                .parse()
                                .map_err(|_| JsValue::from_str("Invalid humidity"))?,
                        );
                        i += 1;
                    }
                }
                "--altitude" => {
                    if i + 1 < args.len() {
                        altitude = Some(
                            args[i + 1]
                                .parse()
                                .map_err(|_| JsValue::from_str("Invalid altitude"))?,
                        );
                        i += 1;
                    }
                }
                "--wind-speed" => {
                    if i + 1 < args.len() {
                        wind_speed = Some(
                            args[i + 1]
                                .parse()
                                .map_err(|_| JsValue::from_str("Invalid wind speed"))?,
                        );
                        i += 1;
                    }
                }
                "--wind-direction" => {
                    if i + 1 < args.len() {
                        wind_direction = Some(
                            args[i + 1]
                                .parse()
                                .map_err(|_| JsValue::from_str("Invalid wind direction"))?,
                        );
                        i += 1;
                    }
                }
                "--target-speed" => {
                    if i + 1 < args.len() {
                        let speed: f64 = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid target speed"))?;
                        // Same bound native lead enforces via f64_range(0.0, 300.0)
                        // (review fix M1) — reject, never silently clamp.
                        if !(0.0..=300.0).contains(&speed) {
                            return Err(JsValue::from_str(&format!(
                                "--target-speed must be between 0 and 300 (mph/m/s), got {}",
                                speed
                            )));
                        }
                        target_speed = Some(speed);
                        i += 1;
                    }
                }
                "--target-angle" => {
                    if i + 1 < args.len() {
                        target_angle = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid target angle"))?;
                        i += 1;
                    }
                }
                "--range" => {
                    if i + 1 < args.len() {
                        range = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid range"))?;
                        i += 1;
                    }
                }
                "--use-powder-sensitivity" => use_powder_sensitivity = true,
                "--powder-temp-sensitivity" => {
                    if i + 1 < args.len() {
                        powder_temp_sensitivity = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid powder temp sensitivity"))?;
                        i += 1;
                    }
                }
                "--powder-temp" => {
                    if i + 1 < args.len() {
                        powder_temp = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid powder temp"))?;
                        powder_temp_provided = true;
                        i += 1;
                    }
                }
                "--powder-temp-curve" => {
                    if i + 1 < args.len() {
                        powder_temp_curve_str = Some(args[i + 1].to_string());
                        i += 1;
                    }
                }
                "--adjustment-unit" => {
                    if i + 1 < args.len() {
                        adjustment_unit = args[i + 1];
                        i += 1;
                    }
                }
                "-o" | "--output" => {
                    if i + 1 < args.len() {
                        lead_output = args[i + 1];
                        i += 1;
                    }
                }
                // --units/-u (+ its value) is consumed globally in run_command, which
                // pre-scans it to set the unit system before dispatch. Skip it here so
                // it isn't rejected as an unknown flag (this is what blocked metric input).
                "--units" | "-u" => {
                    i += 1;
                }
                // Reject unrecognized flags instead of silently ignoring them, so a
                // typo or a flag that isn't wired into this WASM surface is caught
                // immediately rather than looking like a no-op. (The native CLI's clap
                // parser already does this; the hand-rolled WASM parser did not.)
                other if other.starts_with('-') => {
                    return Err(JsValue::from_str(&format!("Unknown flag: {}", other)));
                }
                _ => {}
            }
            i += 1;
        }

        let target_speed = target_speed
            .ok_or_else(|| JsValue::from_str("--target-speed is required"))?;

        // MBA-1355: smoa/iphy join mil/moa as real display units (mirrors native
        // handle_lead's Smoa|Iphy match arm — sol.lead_mil * smoa_per_mil()). clicks is
        // out-of-scope for `lead` — real click resolution only exists for
        // trajectory/come-ups (native's `reject_clicks_out_of_scope`; WASM has no
        // come-ups command, so only trajectory's Ring column resolves clicks).
        let adjustment_unit_lower = adjustment_unit.to_lowercase();
        match adjustment_unit_lower.as_str() {
            "mil" | "moa" | "smoa" | "iphy" => {}
            "clicks" => {
                return Err(JsValue::from_str(
                    "error: --adjustment-unit clicks is currently supported for trajectory and come-ups only (MBA-1355)",
                ));
            }
            _ => {
                return Err(JsValue::from_str(&format!(
                    "Invalid --adjustment-unit '{}' (expected mil, moa, smoa, iphy, or clicks)",
                    adjustment_unit
                )));
            }
        }

        let lead_output_lower = lead_output.to_lowercase();
        if lead_output_lower != "table" && lead_output_lower != "json" {
            return Err(JsValue::from_str(&format!(
                "Invalid --output '{}' (expected table or json)",
                lead_output
            )));
        }

        // Build inputs (mirrors handle_zero_command's unit conversions)
        let mut inputs = InternalBallisticInputs::default();
        match units {
            UnitSystem::Imperial => {
                inputs.muzzle_velocity = velocity * 0.3048;
                inputs.bullet_mass = mass * crate::constants::GRAINS_TO_KG;
                inputs.bullet_diameter = diameter * 0.0254;
                inputs.sight_height = sight_height * 0.0254;
            }
            UnitSystem::Metric => {
                inputs.muzzle_velocity = velocity;
                inputs.bullet_mass = mass * 0.001;
                inputs.bullet_diameter = diameter * 0.001;
                inputs.sight_height = sight_height * 0.001;
            }
        }
        // MBA-1135: mass-based length estimate (mirrors CLI/FFI); replaces the mass-blind
        // 4.5-caliber heuristic. WASM otherwise left it at the struct default.
        inputs.bullet_length =
            crate::stability::estimate_bullet_length_m(inputs.bullet_diameter, inputs.bullet_mass);
        if inputs.bullet_length <= 0.0 {
            inputs.bullet_length = inputs.bullet_diameter * 4.5;
        }

        inputs.bc_value = bc;
        let drag_model_parsed = DragModel::from_str(drag_model)
            .ok_or_else(|| JsValue::from_str("Invalid drag model"))?;
        inputs.bc_type = drag_model_parsed;
        // Custom drag table (MBA-1328): see handle_trajectory_command for the rationale —
        // applied unconditionally, no gate flag. Native `lead` has no --drag-table flag of
        // its own, but calculate_lead solves through the same TrajectorySolver/BallisticInputs
        // as trajectory/zero, so a loaded table is honored identically here.
        if let Some(table) = self.drag_table.borrow().as_ref() {
            inputs.custom_drag_table = Some(table.clone());
        }
        // MBA-1411: --cd-scale requires a loaded drag table, mirroring the native CLI's
        // --cd-scale/--drag-table pairing requirement and the identical check in
        // handle_trajectory_command/handle_zero_command. Validate before any solve.
        if let Some(scale) = cd_scale {
            if self.drag_table.borrow().is_none() {
                return Err(JsValue::from_str(CD_SCALE_REQUIRES_DRAG_TABLE));
            }
            inputs.cd_scale = scale;
        }
        let cd_scale_warning = cd_scale.and_then(cd_scale_range_warning);

        // Powder temperature (MBA-1325): identical resolution to handle_trajectory_command
        // so a curve/sensitivity muzzle-velocity correction reproduces exactly between the
        // two commands. calculate_lead builds a TrajectorySolver from these inputs, which
        // resolves the correction (see cli_api::TrajectorySolver::new) — no further plumbing
        // needed here.
        inputs.use_powder_sensitivity = use_powder_sensitivity;
        inputs.powder_temp_sensitivity = match units {
            UnitSystem::Imperial => powder_temp_sensitivity * 0.3048 / (5.0 / 9.0),
            UnitSystem::Metric => powder_temp_sensitivity,
        };
        inputs.powder_temp = match units {
            UnitSystem::Imperial => (powder_temp - 32.0) * 5.0 / 9.0,
            UnitSystem::Metric => powder_temp,
        };
        // When --powder-temp was explicitly given, it becomes the POWDER temperature the
        // curve is interpolated at (decoupled from --temperature / air density). When not
        // given, powder_curve_temp_c stays None so the curve falls back to the air temp.
        inputs.powder_curve_temp_c = if powder_temp_provided {
            Some(inputs.powder_temp)
        } else {
            None
        };
        // Parse the optional powder-temperature -> velocity curve into SI points
        // (shared parser — see parse_powder_temp_curve_str).
        if let Some(curve_str) = &powder_temp_curve_str {
            inputs.powder_temp_curve = Some(parse_powder_temp_curve_str(curve_str, units)?);
        }

        // Environmental conditions (MBA-1325 env-flags addendum). Start from the exact
        // objects this command has always solved with and override only what was
        // explicitly flagged — absent flags are structurally byte-identical to the
        // pre-flag build. Each override applies handle_trajectory_command's conversion
        // and mirrors native handle_lead's plumbing: the value lands BOTH on the solver
        // condition object (wind/atmosphere — what the solve consumes) and on the
        // matching BallisticInputs field (native's build_trajectory_components sets
        // both; inputs.humidity takes the documented 0-1 fraction per its field
        // contract, while atmosphere.humidity takes percent — same split as
        // handle_trajectory_command).
        let mut wind = WindConditions::default();
        let mut atmosphere = AtmosphericConditions::default();
        if let Some(t) = temperature {
            let t_c = match units {
                UnitSystem::Imperial => (t - 32.0) * 5.0 / 9.0, // F to C
                UnitSystem::Metric => t,
            };
            atmosphere.temperature = t_c;
            inputs.temperature = t_c;
        }
        if let Some(p) = pressure {
            let p_hpa = match units {
                UnitSystem::Imperial => p * 33.863886666667, // inHg to hPa
                UnitSystem::Metric => p,
            };
            atmosphere.pressure = p_hpa;
            inputs.pressure = p_hpa;
        }
        if let Some(h) = humidity {
            atmosphere.humidity = h; // percent
            inputs.humidity = (h / 100.0).clamp(0.0, 1.0); // fraction
        }
        if let Some(a) = altitude {
            let a_m = match units {
                UnitSystem::Imperial => a * 0.3048, // feet to meters
                UnitSystem::Metric => a,
            };
            atmosphere.altitude = a_m;
            inputs.altitude = a_m;
        }
        if let Some(w) = wind_speed {
            let w_ms = match units {
                UnitSystem::Imperial => w * 0.44704, // mph to m/s
                UnitSystem::Metric => w,
            };
            wind.speed = w_ms;
            inputs.wind_speed = w_ms;
        }
        if let Some(d) = wind_direction {
            // Degrees, wind-FROM (0=headwind, 90=from right) -> radians, matching
            // handle_trajectory_command and WindConditions' documented convention.
            let d_rad = d.to_radians();
            wind.direction = d_rad;
            inputs.wind_angle = d_rad;
        }

        // --target-speed uses the same mph (imperial) / m/s (metric) convention as
        // --wind-speed in handle_trajectory_command.
        let target_speed_mps = match units {
            UnitSystem::Imperial => target_speed * 0.44704, // mph to m/s
            UnitSystem::Metric => target_speed,
        };
        let range_m = match units {
            UnitSystem::Imperial => range * 0.9144, // yards to meters
            UnitSystem::Metric => range,
        };

        match calculate_lead(
            inputs,
            wind,
            atmosphere,
            target_speed_mps,
            target_angle,
            range_m,
        ) {
            Ok(sol) => {
                let (dist_unit, speed_unit) = match units {
                    UnitSystem::Imperial => ("yd", "mph"),
                    UnitSystem::Metric => ("m", "m/s"),
                };
                let lead_disp = match units {
                    UnitSystem::Imperial => sol.lead_m * 1.09361, // m to yards
                    UnitSystem::Metric => sol.lead_m,
                };
                let intercept_disp = match units {
                    UnitSystem::Imperial => sol.corrected_range_m * 1.09361,
                    UnitSystem::Metric => sol.corrected_range_m,
                };
                // Requested --adjustment-unit is listed first; MIL is always shown second
                // (MBA-1355: SMOA/IPHY join MOA as a requestable primary unit, sharing the
                // native smoa_per_mil() conversion off sol.lead_mil).
                let lead_smoa = sol.lead_mil * smoa_per_mil();
                let lead_adj_line = match adjustment_unit_lower.as_str() {
                    "moa" => format!("{:.2} MOA / {:.2} MIL", sol.lead_moa, sol.lead_mil),
                    "smoa" => format!("{:.2} SMOA / {:.2} MIL", lead_smoa, sol.lead_mil),
                    "iphy" => format!("{:.2} IPHY / {:.2} MIL", lead_smoa, sol.lead_mil),
                    _ => format!("{:.2} MIL / {:.2} MOA", sol.lead_mil, sol.lead_moa),
                };

                if lead_output_lower == "json" {
                    let payload = serde_json::json!({
                        "target_speed": target_speed,
                        "target_speed_unit": speed_unit,
                        "target_angle_deg": target_angle,
                        "range": range,
                        "distance_unit": dist_unit,
                        "tof_s": sol.time_of_flight_s,
                        "lead": lead_disp,
                        "lead_mil": sol.lead_mil,
                        "lead_moa": sol.lead_moa,
                        "lead_smoa": lead_smoa,
                        "intercept_range": intercept_disp,
                        "iterations": sol.iterations,
                        "adjustment_unit": adjustment_unit_lower,
                    });
                    return Ok(serde_json::to_string_pretty(&payload)
                        .unwrap_or_else(|_| "Error formatting JSON".to_string()));
                }

                // MBA-1411: table-only prepend, same pattern as
                // handle_trajectory_command/handle_zero_command — the JSON branch above
                // already returned, so everything past this point is the table view.
                Ok(format!(
                    "{}Moving-Target Lead\n\
                     ===================\n\
                     Target: {:.1} {} at {:.0}\u{b0} \
                     (0=away, 90=left-to-right, 180=toward, 270=right-to-left;\n\
                     positive lead = hold in direction of travel)\n\n\
                     Range: {:.0} {}\n\
                     Time of Flight: {:.3} s\n\
                     Lead: {:.2} {} ({})\n\
                     Intercept Range: {:.1} {}\n\
                     Iterations: {}\n",
                    cd_scale_warning.as_deref().unwrap_or(""),
                    target_speed,
                    speed_unit,
                    target_angle,
                    range,
                    dist_unit,
                    sol.time_of_flight_s,
                    lead_disp,
                    dist_unit,
                    lead_adj_line,
                    intercept_disp,
                    dist_unit,
                    sol.iterations
                ))
            }
            Err(e) => Err(JsValue::from_str(&format!("Error calculating lead: {}", e))),
        }
    }

    fn handle_monte_carlo_command(
        &self,
        args: &[&str],
        units: UnitSystem,
    ) -> Result<String, JsValue> {
        // Default values
        let (default_velocity, default_mass, default_diameter) = match units {
            UnitSystem::Imperial => (2700.0, 168.0, 0.308),
            UnitSystem::Metric => (823.0, 10.9, 7.82),
        };

        let mut velocity = default_velocity;
        let mut angle = 0.0;
        let mut bc = 0.475;
        let mut mass = default_mass;
        let mut diameter = default_diameter;
        let mut num_sims = 1000;
        // The two std-dev flags stay unset until parsed: the base (non-WEZ) path keeps
        // this handler's historical defaults (10 fps / 3 m/s velocity std, 2 mph /
        // 0.5 m/s wind std), while the --wez path uses the native CLI's clap defaults
        // (1.0 each, user units) so a WEZ sweep with identical args is byte-identical
        // to the native binary (MBA-1343 Phase C).
        let mut velocity_std: Option<f64> = None;
        let mut angle_std = 0.1;
        let mut bc_std = 0.01;
        let mut wind_speed_std: Option<f64> = None;
        let mut wind_direction_std = 0.0;
        let mut drag_model = "G1";
        // MBA-1356: whole-curve drag scale for a custom drag table. None until parsed;
        // resolved against self.drag_table once the arg-parse loop is done, shared by both
        // the base and --wez dispatch below.
        let mut cd_scale: Option<f64> = None;
        // MBA-1343 Phase C: WEZ (Weapon Employment Zone) sweep mode, mirroring the
        // native `monte-carlo --wez` flag set (MBA-1317). The wez-only flags are kept
        // as Options so using one without --wez can be rejected (native clap's
        // `requires = "wez"`), not silently ignored.
        let mut wez = false;
        let mut target_size: Option<String> = None;
        let mut wind_call_error: Option<f64> = None;
        let mut wez_start: Option<f64> = None;
        let mut wez_end: Option<f64> = None;
        let mut wez_step: Option<f64> = None;
        let mut output = "summary";

        // Parse arguments. Every value-taking flag goes through `require_value`, so a
        // flag dangling at the end of the line is an error, never a silent skip; bare
        // non-flag tokens are rejected like native clap's "unexpected argument".
        let mut i = 0;
        while i < args.len() {
            match args[i] {
                "-v" | "--velocity" => {
                    velocity = require_value(args, i)?
                        .parse()
                        .map_err(|_| JsValue::from_str("Invalid velocity"))?;
                    i += 1;
                }
                "-a" | "--angle" => {
                    angle = require_value(args, i)?
                        .parse()
                        .map_err(|_| JsValue::from_str("Invalid angle"))?;
                    i += 1;
                }
                "-b" | "--bc" => {
                    bc = require_value(args, i)?
                        .parse()
                        .map_err(|_| JsValue::from_str("Invalid BC"))?;
                    i += 1;
                }
                "-m" | "--mass" => {
                    mass = require_value(args, i)?
                        .parse()
                        .map_err(|_| JsValue::from_str("Invalid mass"))?;
                    i += 1;
                }
                "-d" | "--diameter" => {
                    diameter = require_value(args, i)?
                        .parse()
                        .map_err(|_| JsValue::from_str("Invalid diameter"))?;
                    i += 1;
                }
                "-n" | "--num-sims" => {
                    num_sims = require_value(args, i)?
                        .parse()
                        .map_err(|_| JsValue::from_str("Invalid number of simulations"))?;
                    i += 1;
                }
                "--velocity-std" => {
                    velocity_std = Some(
                        require_value(args, i)?
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid velocity std"))?,
                    );
                    i += 1;
                }
                "--angle-std" => {
                    angle_std = require_value(args, i)?
                        .parse()
                        .map_err(|_| JsValue::from_str("Invalid angle std"))?;
                    i += 1;
                }
                "--bc-std" => {
                    bc_std = require_value(args, i)?
                        .parse()
                        .map_err(|_| JsValue::from_str("Invalid BC std"))?;
                    i += 1;
                }
                // --wind-std is the native CLI's name for this flag; accept both so a
                // native command line pastes into the terminal unchanged.
                "--wind-speed-std" | "--wind-std" => {
                    wind_speed_std = Some(
                        require_value(args, i)?
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid wind speed std"))?,
                    );
                    i += 1;
                }
                "--wind-direction-std" | "--wind-dir-std" => {
                    wind_direction_std = require_value(args, i)?
                        .parse()
                        .map_err(|_| JsValue::from_str("Invalid wind direction std"))?;
                    i += 1;
                }
                "--drag-model" => {
                    drag_model = require_value(args, i)?;
                    i += 1;
                }
                "--cd-scale" => {
                    cd_scale = Some(
                        require_value(args, i)?
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid cd-scale"))?,
                    );
                    i += 1;
                }
                "--wez" => {
                    wez = true;
                }
                "--target-size" => {
                    target_size = Some(require_value(args, i)?.to_string());
                    i += 1;
                }
                "--wind-call-error" => {
                    wind_call_error = Some(
                        require_value(args, i)?
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid wind call error"))?,
                    );
                    i += 1;
                }
                "--wez-start" => {
                    wez_start = Some(
                        require_value(args, i)?
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid wez start range"))?,
                    );
                    i += 1;
                }
                "--wez-end" => {
                    wez_end = Some(
                        require_value(args, i)?
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid wez end range"))?,
                    );
                    i += 1;
                }
                "--wez-step" => {
                    wez_step = Some(
                        require_value(args, i)?
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid wez step"))?,
                    );
                    i += 1;
                }
                "-o" | "--output" => {
                    output = require_value(args, i)?;
                    i += 1;
                }
                // --units/-u (+ its value) is consumed globally in run_command, which
                // pre-scans it to set the unit system before dispatch. Skip it here so
                // it isn't rejected as an unknown flag (this is what blocked metric input).
                "--units" | "-u" => {
                    i += 1;
                }
                // Reject unrecognized flags instead of silently ignoring them, so a
                // typo or a flag that isn't wired into this WASM surface is caught
                // immediately rather than looking like a no-op. (The native CLI's clap
                // parser already does this; the hand-rolled WASM parser did not.)
                other if other.starts_with('-') => {
                    return Err(JsValue::from_str(&format!("Unknown flag: {}", other)));
                }
                // Bare non-flag tokens are equally errors (native clap: "unexpected
                // argument"): silently ignoring them made e.g. a mistyped flag value
                // look like a successful run (MBA-1343 review).
                other => {
                    return Err(JsValue::from_str(&format!(
                        "Error: unexpected argument '{other}'"
                    )));
                }
            }
            i += 1;
        }

        // MBA-1356: --cd-scale requires a loaded drag table, mirroring the native CLI's
        // --cd-scale/--drag-table pairing requirement. Validate before any solve — shared by
        // both the --wez and base dispatch below, since both read the same self.drag_table.
        if cd_scale.is_some() && self.drag_table.borrow().is_none() {
            return Err(JsValue::from_str(CD_SCALE_REQUIRES_DRAG_TABLE));
        }

        // MBA-1343 Phase C: WEZ sweep mode. Mirrors the native dispatch
        // (Commands::MonteCarlo with `wez: true` in main.rs): convert to metric, run
        // the shared seeded compute (ballistics_engine::wez::compute_wez), and render
        // summary / statistics / full exactly as the native printers do.
        if wez {
            return self.run_monte_carlo_wez(
                units,
                velocity,
                angle,
                bc,
                mass,
                diameter,
                num_sims,
                velocity_std,
                angle_std,
                bc_std,
                wind_speed_std,
                wind_direction_std,
                wind_call_error,
                target_size.as_deref(),
                wez_start,
                wez_end,
                wez_step,
                drag_model,
                cd_scale,
                output,
            );
        }
        // The wez-only flags mirror native clap's `requires = "wez"`: using one
        // without --wez is an error, never a silent no-op.
        if target_size.is_some() {
            return Err(JsValue::from_str("--target-size requires --wez"));
        }
        if wind_call_error.is_some() {
            return Err(JsValue::from_str("--wind-call-error requires --wez"));
        }
        if wez_start.is_some() {
            return Err(JsValue::from_str("--wez-start requires --wez"));
        }
        if wez_end.is_some() {
            return Err(JsValue::from_str("--wez-end requires --wez"));
        }
        if wez_step.is_some() {
            return Err(JsValue::from_str("--wez-step requires --wez"));
        }
        // Base monte-carlo has a single (summary) output format on the WASM surface;
        // the statistics/full modes native supports there remain a documented gap.
        match output.to_lowercase().as_str() {
            "summary" | "table" => {}
            other => {
                return Err(JsValue::from_str(&format!(
                    "-o {other} requires --wez in the WASM terminal (base monte-carlo has a single summary output)"
                )));
            }
        }
        // Historical base-path defaults (pre-WEZ behavior, unchanged).
        let velocity_std = velocity_std.unwrap_or(if units == UnitSystem::Imperial {
            10.0
        } else {
            3.0
        });
        let wind_speed_std = wind_speed_std.unwrap_or(if units == UnitSystem::Imperial {
            2.0
        } else {
            0.5
        });

        // Build inputs
        let mut inputs = InternalBallisticInputs::default();

        // Convert units
        match units {
            UnitSystem::Imperial => {
                inputs.muzzle_velocity = velocity * 0.3048;
                inputs.bullet_mass = mass * crate::constants::GRAINS_TO_KG;
                inputs.bullet_diameter = diameter * 0.0254;
            }
            UnitSystem::Metric => {
                inputs.muzzle_velocity = velocity;
                inputs.bullet_mass = mass * 0.001;
                inputs.bullet_diameter = diameter * 0.001;
            }
        }

        // MBA-1135: mass-based length estimate (mirrors CLI/FFI); replaces the mass-blind
        // 4.5-caliber heuristic. WASM otherwise left it at the struct default.
        inputs.bullet_length =
            crate::stability::estimate_bullet_length_m(inputs.bullet_diameter, inputs.bullet_mass);
        if inputs.bullet_length <= 0.0 {
            inputs.bullet_length = inputs.bullet_diameter * 4.5;
        }

        inputs.bc_value = bc;
        // Honor --drag-model (mirrors the trajectory/zero handlers); previously the Monte
        // Carlo path silently always used the G1 default even when G7 was intended.
        let drag_model_parsed = DragModel::from_str(drag_model)
            .ok_or_else(|| JsValue::from_str("Invalid drag model"))?;
        inputs.bc_type = drag_model_parsed;
        // Custom drag table (MBA-1328): see handle_trajectory_command for the rationale —
        // applied unconditionally, no gate flag, mirrors native --drag-table on `monte-carlo`.
        // Note: with a table active, --bc-std dispersion becomes a no-op for drag (the table
        // fixes Cd directly) — same documented caveat as native (CLI_USAGE.md).
        if let Some(table) = self.drag_table.borrow().as_ref() {
            inputs.custom_drag_table = Some(table.clone());
        }
        // MBA-1356: pairing already validated (before any solve, above); just apply. The
        // warning is prepended to the result string below, table-only.
        if let Some(scale) = cd_scale {
            inputs.cd_scale = scale;
        }
        let cd_scale_warning = cd_scale.and_then(cd_scale_range_warning);
        inputs.muzzle_angle = angle * std::f64::consts::PI / 180.0;
        inputs.muzzle_height = 1.5;
        inputs.ground_threshold = 0.0;

        // Create Monte Carlo parameters
        let params = MonteCarloParams {
            num_simulations: num_sims,
            velocity_std_dev: velocity_std
                * (if units == UnitSystem::Imperial {
                    0.3048
                } else {
                    1.0
                }),
            angle_std_dev: angle_std * std::f64::consts::PI / 180.0,
            bc_std_dev: bc_std,
            wind_speed_std_dev: wind_speed_std
                * (if units == UnitSystem::Imperial {
                    0.44704
                } else {
                    1.0
                }),
            base_wind_speed: 0.0,
            base_wind_direction: 0.0,
            azimuth_std_dev: 0.001,
            target_distance: None,
        };

        let body = match run_monte_carlo_with_direction_std_dev(
            inputs,
            params,
            wind_direction_std * std::f64::consts::PI / 180.0,
        ) {
            Ok(results) => {
                // Calculate statistics
                let mean_range: f64 =
                    results.ranges.iter().sum::<f64>() / results.ranges.len() as f64;
                let mean_velocity: f64 = results.impact_velocities.iter().sum::<f64>()
                    / results.impact_velocities.len() as f64;

                let range_std = {
                    let variance: f64 = results
                        .ranges
                        .iter()
                        .map(|r| (r - mean_range).powi(2))
                        .sum::<f64>()
                        / results.ranges.len() as f64;
                    variance.sqrt()
                };

                let velocity_std = {
                    let variance: f64 = results
                        .impact_velocities
                        .iter()
                        .map(|v| (v - mean_velocity).powi(2))
                        .sum::<f64>()
                        / results.impact_velocities.len() as f64;
                    variance.sqrt()
                };

                // Convert back to display units
                let (range_unit, velocity_unit) = match units {
                    UnitSystem::Imperial => ("yards", "fps"),
                    UnitSystem::Metric => ("meters", "m/s"),
                };

                let display_mean_range = match units {
                    UnitSystem::Imperial => mean_range * 1.09361,
                    UnitSystem::Metric => mean_range,
                };

                let display_range_std = match units {
                    UnitSystem::Imperial => range_std * 1.09361,
                    UnitSystem::Metric => range_std,
                };

                let display_mean_velocity = match units {
                    UnitSystem::Imperial => mean_velocity * 3.28084,
                    UnitSystem::Metric => mean_velocity,
                };

                let display_velocity_std = match units {
                    UnitSystem::Imperial => velocity_std * 3.28084,
                    UnitSystem::Metric => velocity_std,
                };

                format!(
                    "Monte Carlo Simulation Results\n\
                     ==============================\n\
                     Simulations Run: {}\n\n\
                     Range Statistics:\n\
                     - Mean: {:.1} {}\n\
                     - Std Dev: {:.1} {}\n\
                     - Min: {:.1} {}\n\
                     - Max: {:.1} {}\n\n\
                     Impact Velocity Statistics:\n\
                     - Mean: {:.0} {}\n\
                     - Std Dev: {:.0} {}\n\
                     - Min: {:.0} {}\n\
                     - Max: {:.0} {}",
                    num_sims,
                    display_mean_range,
                    range_unit,
                    display_range_std,
                    range_unit,
                    results
                        .ranges
                        .iter()
                        .min_by(|a, b| a.partial_cmp(b).unwrap())
                        .unwrap()
                        * (if units == UnitSystem::Imperial {
                            1.09361
                        } else {
                            1.0
                        }),
                    range_unit,
                    results
                        .ranges
                        .iter()
                        .max_by(|a, b| a.partial_cmp(b).unwrap())
                        .unwrap()
                        * (if units == UnitSystem::Imperial {
                            1.09361
                        } else {
                            1.0
                        }),
                    range_unit,
                    display_mean_velocity,
                    velocity_unit,
                    display_velocity_std,
                    velocity_unit,
                    results
                        .impact_velocities
                        .iter()
                        .min_by(|a, b| a.partial_cmp(b).unwrap())
                        .unwrap()
                        * (if units == UnitSystem::Imperial {
                            3.28084
                        } else {
                            1.0
                        }),
                    velocity_unit,
                    results
                        .impact_velocities
                        .iter()
                        .max_by(|a, b| a.partial_cmp(b).unwrap())
                        .unwrap()
                        * (if units == UnitSystem::Imperial {
                            3.28084
                        } else {
                            1.0
                        }),
                    velocity_unit,
                )
            }
            Err(e) => format!("Error running Monte Carlo simulation: {}", e),
        };
        // MBA-1356: table-only prepend, same pattern as the trajectory handler — this
        // command has no JSON/CSV output mode to protect from contamination.
        Ok(format!("{}{}", cd_scale_warning.as_deref().unwrap_or(""), body))
    }

    /// Native `monte-carlo --wez` parity path (MBA-1343 Phase C): unit conversions
    /// mirror the native dispatch factor-for-factor, the compute is the shared
    /// per-step-seeded [`crate::wez::compute_wez`] core (deterministic, so outputs are
    /// directly comparable to the native binary), and the three renderers replicate
    /// the native printers byte-for-byte.
    #[allow(
        clippy::too_many_arguments,
        reason = "flat arguments mirror the stable Monte Carlo CLI command shape (MBA-1317)"
    )]
    fn run_monte_carlo_wez(
        &self,
        units: UnitSystem,
        velocity: f64,
        angle: f64,
        bc: f64,
        mass: f64,
        diameter: f64,
        num_sims: usize,
        velocity_std: Option<f64>,
        angle_std: f64,
        bc_std: f64,
        wind_speed_std: Option<f64>,
        wind_direction_std: f64,
        wind_call_error: Option<f64>,
        target_size: Option<&str>,
        wez_start: Option<f64>,
        wez_end: Option<f64>,
        wez_step: Option<f64>,
        drag_model: &str,
        cd_scale: Option<f64>,
        output: &str,
    ) -> Result<String, JsValue> {
        // Native clap defaults: --velocity-std 1.0, --wind-std 1.0, --wind-call-error
        // 0.0, --wez-start 200, --wez-end 1000, --wez-step 100 — all in user units.
        // (The base WASM monte-carlo keeps its own historical std defaults; see
        // handle_monte_carlo_command.)
        let velocity_std = velocity_std.unwrap_or(1.0);
        let wind_speed_std = wind_speed_std.unwrap_or(1.0);
        let wind_call_error = wind_call_error.unwrap_or(0.0);
        let wez_start = wez_start.unwrap_or(200.0);
        let wez_end = wez_end.unwrap_or(1000.0);
        let wez_step = wez_step.unwrap_or(100.0);

        // Native MonteCarloOutput mode names first; the WASM terminal's conventional
        // table/csv/json spellings map onto the same three modes (table -> summary,
        // csv -> statistics, json -> full).
        let output = match output.to_lowercase().as_str() {
            "summary" | "table" => WezOutputMode::Summary,
            "statistics" | "csv" => WezOutputMode::Statistics,
            "full" | "json" => WezOutputMode::Full,
            other => {
                return Err(JsValue::from_str(&format!(
                    "Invalid --output '{other}' (expected summary, statistics, or full; \
                     table/csv/json are accepted aliases)"
                )))
            }
        };

        let target_size_spec = target_size.ok_or_else(|| {
            JsValue::from_str("--target-size is required with --wez (e.g. --target-size 18x30)")
        })?;
        let target_size_parsed = crate::wez::parse_target_size(target_size_spec)
            .map_err(|e| JsValue::from_str(&format!("--target-size: {e}")))?;
        let target_size_metric = target_size_parsed.to_metric(engine_units(units));

        // Unit conversions, mirroring the native dispatch (UnitConverter::*_to_metric).
        let (velocity_metric, mass_metric, diameter_metric) = match units {
            UnitSystem::Imperial => (
                velocity * 0.3048,
                mass * crate::constants::GRAINS_TO_KG,
                diameter * 0.0254,
            ),
            UnitSystem::Metric => (velocity, mass * 0.001, diameter * 0.001),
        };
        let velocity_std_metric = match units {
            UnitSystem::Imperial => velocity_std * 0.3048, // fps to m/s
            UnitSystem::Metric => velocity_std,
        };
        let to_wind = |v: f64| match units {
            UnitSystem::Imperial => v * 0.44704, // mph to m/s
            UnitSystem::Metric => v,
        };
        let to_distance = |v: f64| match units {
            UnitSystem::Imperial => v * 0.9144, // yards to meters
            UnitSystem::Metric => v,
        };

        // Custom drag table (MBA-1328): applied unconditionally when loaded, mirroring
        // native --drag-table on `monte-carlo` (see handle_monte_carlo_command).
        let custom_drag_table = self.drag_table.borrow().as_ref().cloned();

        // --drag-model (a WASM-surface extension; the native `monte-carlo` has no such
        // flag and always sweeps G1). Same validation/message as the base path below.
        let drag_model = DragModel::from_str(drag_model)
            .ok_or_else(|| JsValue::from_str("Invalid drag model"))?;

        // MBA-1356: pairing already validated by the caller (handle_monte_carlo_command,
        // before any solve); just resolve the value and the out-of-range nudge here.
        let cd_scale_warning = cd_scale.and_then(cd_scale_range_warning);
        let cd_scale = cd_scale.unwrap_or(1.0);

        // Flags the WASM monte-carlo surface does not (yet) expose keep the native
        // defaults: base wind (--wind-speed / --wind-direction / --wind-vertical) and
        // --cant are all 0.
        let result = crate::wez::compute_wez(
            velocity_metric,
            angle,
            bc,
            mass_metric,
            diameter_metric,
            num_sims,
            velocity_std_metric,
            angle_std,
            bc_std,
            to_wind(wind_speed_std),
            wind_direction_std,
            0.0, // wind_speed
            0.0, // wind_direction
            0.0, // wind_vertical
            to_wind(wind_call_error),
            target_size_metric,
            to_distance(wez_start),
            to_distance(wez_end),
            to_distance(wez_step),
            drag_model,
            custom_drag_table,
            cd_scale,
            0.0, // cant
        )
        .map_err(|e| JsValue::from_str(&e.to_string()))?;

        match output {
            // JSON/CSV must stay pure machine output — the cd_scale nudge is table-only,
            // same convention as every other prepended warning in this file.
            WezOutputMode::Full => format_wez_full(&result),
            WezOutputMode::Statistics => Ok(format_wez_statistics(&result, units)),
            WezOutputMode::Summary => Ok(format!(
                "{}{}",
                cd_scale_warning.as_deref().unwrap_or(""),
                format_wez_summary(&result, units)
            )),
        }
    }

    /// `true-velocity` (MBA-1343 Phase C): the native command's local compute paths on
    /// the WASM terminal. Zero `--observed` flags -> the classic single-observation
    /// binary search ([`crate::truing::calculate_true_velocity_local`]); one or more ->
    /// the MBA-1316 joint MV+BC calibration
    /// ([`crate::truing::run_multi_observation_truing_core`]). Unit conversions and all
    /// three output formats mirror the native dispatch/printers byte-for-byte.
    fn handle_true_velocity_command(
        &self,
        args: &[&str],
        units: UnitSystem,
    ) -> Result<String, JsValue> {
        use crate::truing::{DragModelArg, DropUnit};

        // Default bullet parameters follow this file's house convention (the native
        // command clap-requires -b/-m/-d instead).
        let (default_mass, default_diameter) = match units {
            UnitSystem::Imperial => (168.0, 0.308),
            UnitSystem::Metric => (10.9, 7.82),
        };

        let mut measured_drop: Option<f64> = None;
        let mut range: Option<f64> = None;
        let mut observed: Vec<String> = Vec::new();
        let mut drop_unit = DropUnit::Mil;
        let mut bc = 0.475;
        let mut drag_model = "g1";
        let mut mass = default_mass;
        let mut diameter = default_diameter;
        let mut chrono_velocity: Option<f64> = None;
        let mut zero_distance = 100.0;
        let mut sight_height: Option<f64> = None;
        let mut temperature: Option<f64> = None;
        let mut pressure: Option<f64> = None;
        let mut humidity = 50.0;
        let mut altitude = 0.0;
        let mut output = "table";

        // Every value-taking flag goes through `require_value`, so a flag dangling at
        // the end of the line is an error, never a silent skip (a dangling --observed
        // used to silently fall back to single-observation truing); bare non-flag
        // tokens are rejected like native clap's "unexpected argument" (a space-
        // separated second observation, e.g. `--observed 600:4.8 700:5.9`, used to
        // silently drop the second point and corrupt the fit — MBA-1343 review).
        let mut i = 0;
        while i < args.len() {
            match args[i] {
                "--measured-drop" => {
                    measured_drop = Some(
                        require_value(args, i)?
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid measured drop"))?,
                    );
                    i += 1;
                }
                "--range" => {
                    range = Some(
                        require_value(args, i)?
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid range"))?,
                    );
                    i += 1;
                }
                "--observed" => {
                    observed.push(require_value(args, i)?.to_string());
                    i += 1;
                }
                "--drop-unit" => {
                    drop_unit = DropUnit::parse(require_value(args, i)?)
                        .map_err(|e| JsValue::from_str(&e))?;
                    i += 1;
                }
                "-b" | "--bc" => {
                    bc = require_value(args, i)?
                        .parse()
                        .map_err(|_| JsValue::from_str("Invalid BC"))?;
                    i += 1;
                }
                "--drag-model" => {
                    drag_model = require_value(args, i)?;
                    i += 1;
                }
                "-m" | "--mass" => {
                    mass = require_value(args, i)?
                        .parse()
                        .map_err(|_| JsValue::from_str("Invalid mass"))?;
                    i += 1;
                }
                "-d" | "--diameter" => {
                    diameter = require_value(args, i)?
                        .parse()
                        .map_err(|_| JsValue::from_str("Invalid diameter"))?;
                    i += 1;
                }
                "--chrono-velocity" => {
                    chrono_velocity = Some(
                        require_value(args, i)?
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid chrono velocity"))?,
                    );
                    i += 1;
                }
                "--zero-distance" => {
                    zero_distance = require_value(args, i)?
                        .parse()
                        .map_err(|_| JsValue::from_str("Invalid zero distance"))?;
                    i += 1;
                }
                "--sight-height" => {
                    sight_height = Some(
                        require_value(args, i)?
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid sight height"))?,
                    );
                    i += 1;
                }
                "--temperature" => {
                    temperature = Some(
                        require_value(args, i)?
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid temperature"))?,
                    );
                    i += 1;
                }
                "--pressure" => {
                    pressure = Some(
                        require_value(args, i)?
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid pressure"))?,
                    );
                    i += 1;
                }
                "--humidity" => {
                    humidity = require_value(args, i)?
                        .parse()
                        .map_err(|_| JsValue::from_str("Invalid humidity"))?;
                    i += 1;
                }
                "--altitude" => {
                    altitude = require_value(args, i)?
                        .parse()
                        .map_err(|_| JsValue::from_str("Invalid altitude"))?;
                    i += 1;
                }
                "-o" | "--output" => {
                    output = require_value(args, i)?;
                    i += 1;
                }
                // Native's --offline forces the local calculation; the WASM terminal is
                // always local, so accept it as a harmless no-op for command parity.
                "--offline" => {}
                // --units/-u (+ its value) is consumed globally in run_command, which
                // pre-scans it to set the unit system before dispatch. Skip it here so
                // it isn't rejected as an unknown flag.
                "--units" | "-u" => {
                    i += 1;
                }
                // Reject unrecognized flags instead of silently ignoring them (this
                // includes the native online/BC5D flags this surface does not support).
                other if other.starts_with('-') => {
                    return Err(JsValue::from_str(&format!("Unknown flag: {}", other)));
                }
                // Bare non-flag tokens are equally errors (native clap: "unexpected
                // argument").
                other => {
                    return Err(JsValue::from_str(&format!(
                        "Error: unexpected argument '{other}'"
                    )));
                }
            }
            i += 1;
        }

        let measured_drop =
            measured_drop.ok_or_else(|| JsValue::from_str("--measured-drop is required"))?;
        let range = range.ok_or_else(|| JsValue::from_str("--range is required"))?;

        let drag_model_arg = match drag_model.to_lowercase().as_str() {
            "g1" => DragModelArg::G1,
            "g7" => DragModelArg::G7,
            _ => return Err(JsValue::from_str("Invalid drag model (expected G1 or G7)")),
        };

        let output = match output.to_lowercase().as_str() {
            "table" => OutputFormat::Table,
            "json" => OutputFormat::Json,
            "csv" => OutputFormat::Csv,
            other => {
                return Err(JsValue::from_str(&format!(
                    "Invalid --output '{other}' (expected table, json, or csv)"
                )))
            }
        };

        if !(0.0..=100.0).contains(&humidity) {
            return Err(JsValue::from_str("--humidity must be between 0 and 100"));
        }

        // Native clap's f64_range validators for the bullet parameters (same bounds
        // under either unit system). The hand-rolled parser previously accepted any
        // finite (or NaN) value here.
        if !(0.001..=2.0).contains(&bc) {
            return Err(JsValue::from_str(&format!(
                "Error: invalid value '{bc}' for '--bc': must be in range 0.001..=2"
            )));
        }
        if !(0.1..=2000.0).contains(&mass) {
            return Err(JsValue::from_str(&format!(
                "Error: invalid value '{mass}' for '--mass': must be in range 0.1..=2000"
            )));
        }
        if !(0.01..=60.0).contains(&diameter) {
            return Err(JsValue::from_str(&format!(
                "Error: invalid value '{diameter}' for '--diameter': must be in range 0.01..=60"
            )));
        }

        // Resolve temperature/pressure AFTER units are known — replicates the native
        // UnitConverter::resolve_temperature / resolve_pressure defaults and range
        // checks (identical error strings).
        let temperature = match (temperature, units) {
            (None, UnitSystem::Imperial) => 59.0,
            (None, UnitSystem::Metric) => 15.0,
            (Some(v), UnitSystem::Imperial) => {
                if !(-148.0..=392.0).contains(&v) {
                    return Err(JsValue::from_str(&format!(
                        "--temperature {v} F is out of range (expected ~-148..392 F for imperial units)"
                    )));
                }
                v
            }
            (Some(v), UnitSystem::Metric) => {
                if !(-100.0..=200.0).contains(&v) {
                    return Err(JsValue::from_str(&format!(
                        "--temperature {v} C is out of range (expected ~-100..200 C for metric units)"
                    )));
                }
                v
            }
        };
        let pressure = match (pressure, units) {
            (None, UnitSystem::Imperial) => 29.92,
            (None, UnitSystem::Metric) => 1013.25,
            (Some(v), UnitSystem::Imperial) => {
                if !(8.0..=33.0).contains(&v) {
                    return Err(JsValue::from_str(&format!(
                        "--pressure {v} inHg is out of range (expected ~8..33 inHg for imperial units)"
                    )));
                }
                v
            }
            (Some(v), UnitSystem::Metric) => {
                if !(250.0..=1100.0).contains(&v) {
                    return Err(JsValue::from_str(&format!(
                        "--pressure {v} hPa is out of range (expected ~250..1100 hPa for metric units)"
                    )));
                }
                v
            }
        };

        // Convert to the truing core's internal imperial units — factor-for-factor the
        // native Commands::TrueVelocity dispatch.
        let range_yd = match units {
            UnitSystem::Imperial => range,
            UnitSystem::Metric => range / 0.9144, // meters to yards
        };
        let weight_gr = match units {
            UnitSystem::Imperial => mass,
            UnitSystem::Metric => mass / crate::constants::GRAMS_PER_GRAIN, // grams to grains
        };
        let caliber_in = match units {
            UnitSystem::Imperial => diameter,
            UnitSystem::Metric => diameter / 25.4, // mm to inches
        };
        let chrono_fps = chrono_velocity.map(|v| match units {
            UnitSystem::Imperial => v,
            UnitSystem::Metric => v / 0.3048, // m/s to fps
        });
        let zero_yd = match units {
            UnitSystem::Imperial => zero_distance,
            UnitSystem::Metric => zero_distance / 0.9144,
        };
        let sight_height_default = match units {
            UnitSystem::Imperial => 2.0,
            UnitSystem::Metric => 50.0,
        };
        let sight_height = sight_height.unwrap_or(sight_height_default);
        let sight_in = match units {
            UnitSystem::Imperial => sight_height,
            UnitSystem::Metric => sight_height / 25.4, // mm to inches
        };
        let temp_f = match units {
            UnitSystem::Imperial => temperature,
            UnitSystem::Metric => temperature * 9.0 / 5.0 + 32.0, // C to F
        };
        let press_inhg = match units {
            UnitSystem::Imperial => pressure,
            UnitSystem::Metric => pressure / 33.8639, // hPa to inHg
        };
        let alt_ft = match units {
            UnitSystem::Imperial => altitude,
            UnitSystem::Metric => altitude / 0.3048, // meters to feet
        };

        // No BC5D segment synthesis on this path: native gates it behind
        // --bc-table-dir / --bc-table-auto, which the WASM terminal does not expose
        // (the loaded-BC5D host API only wires into `trajectory --use-bc-segments`).
        let bc_segments: Option<Vec<crate::BCSegmentData>> = None;

        if !observed.is_empty() {
            // MBA-1316: one or more --observed impacts -> joint MV+BC calibration via
            // the shared core. Token parsing happens here (the typed core takes
            // parsed observations); the parse and validation errors (e.g. duplicate
            // ranges) carry the native messages verbatim.
            let mut observations: Vec<crate::truing::TruingObservation> =
                Vec::with_capacity(observed.len() + 1);
            observations.push(crate::truing::TruingObservation {
                range_yd,
                drop: measured_drop,
            });
            for token in &observed {
                observations.push(
                    crate::truing::parse_truing_observation(token, engine_units(units))
                        .map_err(|e| JsValue::from_str(&e))?,
                );
            }
            let report = crate::truing::run_multi_observation_truing_core(
                &observations,
                drop_unit,
                bc,
                drag_model_arg,
                weight_gr,
                caliber_in,
                zero_yd,
                sight_in,
                temp_f,
                press_inhg,
                humidity,
                alt_ft,
                &bc_segments,
            )
            .map_err(|e| JsValue::from_str(&e.to_string()))?;
            Ok(format_multi_truing_result(
                &report, drop_unit, units, chrono_fps, output,
            ))
        } else {
            // Classic single-observation velocity truing (drop is always MIL here).
            let result = crate::truing::calculate_true_velocity_local(
                measured_drop,
                range_yd,
                bc,
                drag_model_arg,
                weight_gr,
                caliber_in,
                zero_yd,
                sight_in,
                temp_f,
                press_inhg,
                humidity,
                alt_ft,
                &bc_segments,
            )
            .map_err(|e| JsValue::from_str(&format!("Local calculation failed: {}", e)))?;

            let velocity_adjustment = chrono_fps.map(|c| result.effective_velocity_fps - c);
            let adjustment_percent =
                chrono_fps.map(|c| (result.effective_velocity_fps - c) / c * 100.0);

            let effective_vel = match units {
                UnitSystem::Imperial => result.effective_velocity_fps,
                UnitSystem::Metric => result.effective_velocity_fps * 0.3048,
            };
            let vel_unit = match units {
                UnitSystem::Imperial => "fps",
                UnitSystem::Metric => "m/s",
            };

            format_true_velocity_result(
                effective_vel,
                vel_unit,
                velocity_adjustment,
                adjustment_percent,
                &result.confidence,
                result.iterations,
                result.final_error_mil,
                result.calculated_drop_mil,
                measured_drop,
                units,
                output,
                false,
            )
        }
    }

    fn handle_estimate_bc_command(
        &self,
        args: &[&str],
        units: UnitSystem,
    ) -> Result<String, JsValue> {
        let (default_velocity, default_mass, default_diameter) = match units {
            UnitSystem::Imperial => (2700.0, 168.0, 0.308),
            UnitSystem::Metric => (823.0, 10.9, 7.82),
        };

        let mut velocity = default_velocity;
        let mut mass = default_mass;
        let mut diameter = default_diameter;
        let mut data_points: Vec<(f64, f64)> = Vec::new();
        let mut vel_points: Vec<(f64, f64)> = Vec::new();
        let mut drag_model_str = String::from("both");
        let mut zero_range: Option<f64> = None;
        let mut sight_height: Option<f64> = None;
        let mut temperature: Option<f64> = None;
        let mut pressure: Option<f64> = None;
        let mut humidity: f64 = 50.0;
        let mut altitude: f64 = 0.0;

        // Parse "d,v;d,v;..." pairs, tolerating surrounding quotes/whitespace.
        fn parse_pairs(raw: &str) -> Result<Vec<(f64, f64)>, JsValue> {
            let s = raw.trim().trim_matches('\'').trim_matches('"');
            let mut out = Vec::new();
            for pair in s.split(';') {
                let pair = pair.trim();
                if pair.is_empty() {
                    continue;
                }
                let parts: Vec<&str> = pair.split(',').collect();
                if parts.len() != 2 {
                    return Err(JsValue::from_str(&format!(
                        "Malformed data pair '{}': expected \"distance,value\"",
                        pair
                    )));
                }
                let d: f64 = parts[0]
                    .trim()
                    .parse()
                    .map_err(|_| JsValue::from_str(&format!("Invalid distance '{}'", parts[0].trim())))?;
                let v: f64 = parts[1]
                    .trim()
                    .parse()
                    .map_err(|_| JsValue::from_str(&format!("Invalid value '{}'", parts[1].trim())))?;
                out.push((d, v));
            }
            Ok(out)
        }

        // Parse arguments
        let mut i = 0;
        while i < args.len() {
            match args[i] {
                "-v" | "--velocity" => {
                    if i + 1 < args.len() {
                        velocity = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid velocity"))?;
                        i += 1;
                    }
                }
                "-m" | "--mass" => {
                    if i + 1 < args.len() {
                        mass = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid mass"))?;
                        i += 1;
                    }
                }
                "-d" | "--diameter" => {
                    if i + 1 < args.len() {
                        diameter = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid diameter"))?;
                        i += 1;
                    }
                }
                "--data" => {
                    // Drop data: distance,drop pairs.
                    if i + 1 < args.len() {
                        data_points = parse_pairs(args[i + 1])?;
                        i += 1;
                    }
                }
                "--velocity-data" => {
                    // Velocity-retention data: distance,velocity pairs.
                    if i + 1 < args.len() {
                        vel_points = parse_pairs(args[i + 1])?;
                        i += 1;
                    }
                }
                "--drag-model" => {
                    if i + 1 < args.len() {
                        drag_model_str = args[i + 1].to_lowercase();
                        i += 1;
                    }
                }
                "--zero-range" => {
                    if i + 1 < args.len() {
                        zero_range = Some(args[i + 1].parse().map_err(|_| JsValue::from_str("Invalid zero-range"))?);
                        i += 1;
                    }
                }
                "--sight-height" => {
                    if i + 1 < args.len() {
                        sight_height = Some(args[i + 1].parse().map_err(|_| JsValue::from_str("Invalid sight-height"))?);
                        i += 1;
                    }
                }
                "--temperature" => {
                    if i + 1 < args.len() {
                        temperature = Some(args[i + 1].parse().map_err(|_| JsValue::from_str("Invalid temperature"))?);
                        i += 1;
                    }
                }
                "--pressure" => {
                    if i + 1 < args.len() {
                        pressure = Some(args[i + 1].parse().map_err(|_| JsValue::from_str("Invalid pressure"))?);
                        i += 1;
                    }
                }
                "--humidity" => {
                    if i + 1 < args.len() {
                        humidity = args[i + 1].parse().map_err(|_| JsValue::from_str("Invalid humidity"))?;
                        i += 1;
                    }
                }
                "--altitude" => {
                    if i + 1 < args.len() {
                        altitude = args[i + 1].parse().map_err(|_| JsValue::from_str("Invalid altitude"))?;
                        i += 1;
                    }
                }
                // --units/-u (+ its value) is consumed globally in run_command, which
                // pre-scans it to set the unit system before dispatch. Skip it here so
                // it isn't rejected as an unknown flag (this is what blocked metric input).
                "--units" | "-u" => {
                    i += 1;
                }
                // Reject unrecognized flags instead of silently ignoring them, so a
                // typo or a flag that isn't wired into this WASM surface is caught
                // immediately rather than looking like a no-op. (The native CLI's clap
                // parser already does this; the hand-rolled WASM parser did not.)
                other if other.starts_with('-') => {
                    return Err(JsValue::from_str(&format!("Unknown flag: {}", other)));
                }
                _ => {}
            }
            i += 1;
        }

        if data_points.is_empty() && vel_points.is_empty() {
            return Ok("Error: No data provided. Use --data \"dist,drop;...\" and/or \
                 --velocity-data \"dist,vel;...\".\nExample: --data \"300,29.0;500,89.9;700,204.6\""
                .to_string());
        }

        // Select drag models to estimate. `both`/`all` fits G1 + G7; any other recognized
        // family name (MBA-1386 widened DragModel past the scalar G1/G7 pair) fits just
        // that one family, matching the native CLI's `parse_drag_models`.
        let models: Vec<DragModel> = match drag_model_str.as_str() {
            "both" | "all" | "g1,g7" | "g1g7" => vec![DragModel::G1, DragModel::G7],
            other => match DragModel::from_str(other) {
                Some(model) => vec![model],
                None => {
                    return Err(JsValue::from_str(&format!(
                        "Unknown --drag-model '{}'; use g1, g7, both, or a specific drag family \
                         (g2, g5, g6, g8, gi, gs, ra4).",
                        other
                    )))
                }
            },
        };

        // Convert scalar inputs to metric.
        let velocity_mps = match units {
            UnitSystem::Imperial => velocity * 0.3048,
            UnitSystem::Metric => velocity,
        };
        let mass_kg = match units {
            UnitSystem::Imperial => mass * crate::constants::GRAINS_TO_KG,
            UnitSystem::Metric => mass * 0.001,
        };
        let diameter_m = match units {
            UnitSystem::Imperial => diameter * 0.0254,
            UnitSystem::Metric => diameter * 0.001,
        };

        // Atmosphere the data was measured at (defaults = ICAO standard). BC only means
        // something relative to air density, so this must match the dope card's conditions.
        let atmosphere = AtmosphericConditions {
            temperature: temperature
                .map(|t| match units {
                    UnitSystem::Imperial => (t - 32.0) * 5.0 / 9.0,
                    UnitSystem::Metric => t,
                })
                .unwrap_or(15.0),
            pressure: pressure
                .map(|p| match units {
                    UnitSystem::Imperial => p * 33.8639,
                    UnitSystem::Metric => p,
                })
                .unwrap_or(1013.25),
            humidity,
            altitude: match units {
                UnitSystem::Imperial => altitude * 0.3048,
                UnitSystem::Metric => altitude,
            },
        };
        let zero_m = zero_range.map(|z| match units {
            UnitSystem::Imperial => z * 0.9144,
            UnitSystem::Metric => z,
        });
        let sight_m = sight_height
            .map(|s| match units {
                UnitSystem::Imperial => s * 0.0254,
                UnitSystem::Metric => s / 1000.0,
            })
            .unwrap_or(0.05);

        // Convert both series to metric (drop -> m, velocity -> m/s; metric drop input is mm).
        let drop_metric: Vec<(f64, f64)> = data_points
            .iter()
            .map(|(dist, drop)| match units {
                UnitSystem::Imperial => (*dist * 0.9144, *drop * 0.0254),
                UnitSystem::Metric => (*dist, *drop * 0.001),
            })
            .collect();
        let vel_metric: Vec<(f64, f64)> = vel_points
            .iter()
            .map(|(dist, v)| match units {
                UnitSystem::Imperial => (*dist * 0.9144, *v * 0.3048),
                UnitSystem::Metric => (*dist, *v),
            })
            .collect();

        // Data bases in a stable order (drop first, then velocity).
        let mut bases: Vec<(BcFitMode, &[(f64, f64)])> = Vec::new();
        if !drop_metric.is_empty() {
            bases.push((BcFitMode::Drop, &drop_metric));
        }
        if !vel_metric.is_empty() {
            bases.push((BcFitMode::Velocity, &vel_metric));
        }

        let (vu, mu, du) = match units {
            UnitSystem::Imperial => ("fps", "grains", "inches"),
            UnitSystem::Metric => ("m/s", "grams", "mm"),
        };
        let mut lines = vec![
            "BC Estimation Results".to_string(),
            "=====================".to_string(),
            format!(
                "Inputs: v={} {}, m={} {}, d={} {}",
                velocity, vu, mass, mu, diameter, du
            ),
            String::new(),
            format!(
                "  {:<6} {:<20} {:>12}   {}",
                "Model", "Fit basis", "Estimated BC", "Fit RMS"
            ),
        ];
        // Guard the common mistake: a zeroed dope card (a point with ~0 drop) fed without
        // --zero-range, which makes the drop fit bore-referenced and returns a wrong BC.
        if zero_m.is_none()
            && drop_metric.iter().any(|(_, dr)| dr.abs() < 0.05)
            && drop_metric.iter().any(|(_, dr)| dr.abs() > 0.25)
        {
            let zd = drop_metric
                .iter()
                .min_by(|a, b| a.1.abs().partial_cmp(&b.1.abs()).unwrap())
                .map(|(d, _)| match units {
                    UnitSystem::Imperial => d / 0.9144,
                    UnitSystem::Metric => *d,
                })
                .unwrap_or(0.0);
            let du2 = if units == UnitSystem::Imperial { "yd" } else { "m" };
            lines.insert(0, format!(
                "⚠ Data looks zeroed near {zd:.0} {du2} but --zero-range not given; drop is being"
            ));
            lines.insert(1, format!(
                "  treated as bore-referenced. For a dope card, pass --zero-range {zd:.0}."
            ));
            lines.insert(2, String::new());
        }

        let mut any_unreliable = false;
        for &model in &models {
            for &(mode, pts) in &bases {
                // Zero range only shapes a drop fit; a velocity fit is frame-independent.
                let zr = match mode {
                    BcFitMode::Drop => zero_m,
                    BcFitMode::Velocity => None,
                };
                let est = estimate_bc_fit(
                    velocity_mps, mass_kg, diameter_m, pts, model, mode,
                    atmosphere.clone(), zr, sight_m,
                )
                .map_err(|e| JsValue::from_str(&format!("Error estimating BC: {}", e)))?;
                if est.at_bound {
                    any_unreliable = true;
                }
                let (rms_user, unit) = match mode {
                    BcFitMode::Drop => match units {
                        UnitSystem::Imperial => (est.rms_error / 0.0254, "in"),
                        UnitSystem::Metric => (est.rms_error * 1000.0, "mm"),
                    },
                    BcFitMode::Velocity => match units {
                        UnitSystem::Imperial => (est.rms_error / 0.3048, "fps"),
                        UnitSystem::Metric => (est.rms_error, "m/s"),
                    },
                };
                // 0.28.1 sweep: label every family by its actual name (DragModel's Display
                // is exactly its variant name), not a generic "G?" for anything past G1/G7.
                let model_name = model.to_string();
                let basis = match mode {
                    BcFitMode::Drop => "drop",
                    BcFitMode::Velocity => "velocity",
                };
                lines.push(format!(
                    "  {:<6} {:<20} {:>12.3}   {:>6.2} {:<4}{}",
                    model_name,
                    format!("{} ({} pts)", basis, pts.len()),
                    est.bc,
                    rms_user,
                    unit,
                    if est.at_bound { " ⚠ UNRELIABLE (hit BC limit)" } else { "" }
                ));
            }
        }
        if any_unreliable {
            lines.push(String::new());
            lines.push("⚠ A fit ran to the BC search limit — the data did not determine a real".to_string());
            lines.push("  value. Add more/longer-range points and check --zero-range / --temperature.".to_string());
        }
        Ok(lines.join("\n"))
    }

    fn format_trajectory_table(
        &self,
        result: &crate::cli_api::TrajectoryResult,
        zero_distance: Option<f64>,
        units: UnitSystem,
        full: bool,
        los_height_m: f64,
        target_speed_mps: f64,
        ring_unit: RingDisplayUnit,
    ) -> String {
        // Mover ring (MBA-1325): additive "Ring" column, only when --target-speed > 0.
        let ring_enabled = target_speed_mps > 0.0;
        // MBA-1355: header carries the unit — Ring(mil)/Ring(moa)/Ring(smoa)/Ring(iphy)/
        // Ring(clicks) — matching native's `ring_hdr` convention exactly (CLI_USAGE.md:
        // "the header/column suffix follows the same (mil)/(moa) convention as every
        // other unit — e.g. the mover Ring column reads Ring(clicks)").
        let ring_hdr_unit = match ring_unit {
            RingDisplayUnit::Factor(_, label) => label,
            RingDisplayUnit::Clicks(_) => "clicks",
        };
        let mut output = String::new();
        output.push_str("Trajectory Calculation Results\n");
        output.push_str("==============================\n\n");
        if ring_enabled {
            output.push_str(&format!(
                "Range | Drop | Drift | Velocity | Energy | Time | Ring({})\n",
                ring_hdr_unit
            ));
            output.push_str("------|------|-------|----------|--------|------|------\n");
        } else {
            output.push_str("Range | Drop | Drift | Velocity | Energy | Time\n");
            output.push_str("------|------|-------|----------|--------|------\n");
        }

        // Determine sampling interval based on max range and full flag
        let max_range_display = match units {
            UnitSystem::Imperial => result.max_range * 1.09361, // m to yards
            UnitSystem::Metric => result.max_range,
        };

        let sample_interval = if full {
            if max_range_display < 100.0 {
                10.0
            } else {
                25.0
            }
        } else {
            if max_range_display < 500.0 {
                50.0
            } else {
                100.0
            }
        };

        let mut current_range = 0.0;

        // The line of sight sits muzzle_height + sight_height above ground at EVERY
        // cant angle (the scope stays on the fixed LOS; the bore rotates around it).
        // Deriving it from points[0].y + sight_height double-counted the canted
        // bore's vertical rise sh*(1-cos c) — Bero's 90-degree report (MBA-1297):
        // reported drop was true miss + one full sight height.
        let los_height = los_height_m;

        for (idx, point) in result.points.iter().enumerate() {
            // McCoy coordinate system: X=downrange, Y=vertical, Z=lateral
            let range_display = match units {
                UnitSystem::Imperial => point.position.x * 1.09361, // X is downrange (McCoy), m to yards
                UnitSystem::Metric => point.position.x,             // X is downrange (McCoy)
            };

            let is_last_point = idx == result.points.len() - 1;

            // Show point if it's at the sampling interval OR if it's the last point OR if it's the zero distance
            let should_show = range_display >= current_range
                || is_last_point
                || (zero_distance.is_some()
                    && (range_display - zero_distance.unwrap()).abs() < 1.0);

            if should_show {
                let drop = los_height - point.position.y;
                let drift = point.position.z; // Z is lateral (windage, McCoy)
                let velocity = point.velocity_magnitude;

                // Format values based on unit system
                let (range_str, drop_str, drift_str, velocity_str, energy_str) = match units {
                    UnitSystem::Imperial => {
                        let energy_ftlb = 0.5
                            * (result.points[0].kinetic_energy / 0.5)
                            * (velocity / result.points[0].velocity_magnitude).powi(2)
                            * 0.737562149;
                        (
                            format!("{:.0} yd", range_display),
                            format!("{:.1} in", drop * 39.3701),
                            format!("{:.1} in", drift * 39.3701),
                            format!("{:.0} fps", velocity * 3.28084),
                            format!("{:.0} ft-lb", energy_ftlb),
                        )
                    }
                    UnitSystem::Metric => (
                        format!("{:.0} m", range_display),
                        format!("{:.1} cm", drop * 100.0),
                        format!("{:.1} cm", drift * 100.0),
                        format!("{:.0} m/s", velocity),
                        format!("{:.0} J", point.kinetic_energy),
                    ),
                };

                if ring_enabled {
                    let (_, ring_mil) =
                        mover_ring(target_speed_mps, point.time, point.position.x);
                    let ring_str = match ring_mil {
                        Some(mil) => match ring_unit {
                            RingDisplayUnit::Factor(f, label) => format!("{:.2} {}", mil * f, label),
                            // clicks_for(drop_yd, range_yd, click) only needs the
                            // drop_yd/range_yd RATIO — passing (mil, 1000.0) reuses it
                            // directly on an already-computed mil angle (ring_mil is
                            // `ring_m / downrange_m * 1000`), same trick native's
                            // RingUnit::Clicks arm uses (main.rs run_trajectory).
                            RingDisplayUnit::Clicks(click) => {
                                format!("{}", crate::adjustment::clicks_for(mil, 1000.0, &click))
                            }
                        },
                        None => "-".to_string(),
                    };
                    output.push_str(&format!(
                        "{:6} | {:6} | {:7} | {:10} | {:8} | {:.3} s | {}\n",
                        range_str, drop_str, drift_str, velocity_str, energy_str, point.time, ring_str
                    ));
                } else {
                    output.push_str(&format!(
                        "{:6} | {:6} | {:7} | {:10} | {:8} | {:.3} s\n",
                        range_str, drop_str, drift_str, velocity_str, energy_str, point.time
                    ));
                }

                if range_display >= current_range {
                    current_range += sample_interval;
                }
            }
        }

        // Add summary
        output.push_str(&format!(
            "\nMax Range: {:.0} {}\n",
            if units == UnitSystem::Imperial {
                result.max_range * 1.09361
            } else {
                result.max_range
            },
            if units == UnitSystem::Imperial {
                "yards"
            } else {
                "meters"
            }
        ));

        // Check termination reason
        if let Some(last_point) = result.points.last() {
            // Debug info about last point
            // McCoy coordinate system: X=downrange, Y=vertical, Z=lateral
            let last_height = last_point.position.y;
            let last_range = last_point.position.x; // X is downrange (McCoy)
            let last_time = last_point.time;

            // Ground threshold is typically around 0.01m (1cm), so if y is close to or below that, it hit ground
            if last_height <= 0.01 {
                output.push_str(&format!(
                    "Bullet struck ground at {:.0} {}\n",
                    if units == UnitSystem::Imperial {
                        last_range * 1.09361
                    } else {
                        last_range
                    },
                    if units == UnitSystem::Imperial {
                        "yards"
                    } else {
                        "meters"
                    }
                ));
            } else if last_time >= 99.0 {
                output.push_str(&format!(
                    "Trajectory timeout at {:.0} {} (time limit reached)\n",
                    if units == UnitSystem::Imperial {
                        last_range * 1.09361
                    } else {
                        last_range
                    },
                    if units == UnitSystem::Imperial {
                        "yards"
                    } else {
                        "meters"
                    }
                ));
            }
        }

        output.push_str(&format!(
            "Max Height: {:.1} {}\n",
            if units == UnitSystem::Imperial {
                result.max_height * 39.3701
            } else {
                result.max_height * 100.0
            },
            if units == UnitSystem::Imperial {
                "inches"
            } else {
                "cm"
            }
        ));
        output.push_str(&format!(
            "Time of Flight: {:.2} seconds\n",
            result.time_of_flight
        ));
        output.push_str(&format!(
            "Impact Velocity: {:.0} {}\n",
            if units == UnitSystem::Imperial {
                result.impact_velocity * 3.28084
            } else {
                result.impact_velocity
            },
            if units == UnitSystem::Imperial {
                "fps"
            } else {
                "m/s"
            }
        ));

        output
    }

    fn format_trajectory_json(
        &self,
        result: &crate::cli_api::TrajectoryResult,
        units: UnitSystem,
        los_height_m: f64,
        target_speed_mps: f64,
    ) -> String {
        // LOS height is cant-invariant (see format_trajectory_table).
        let los_height = los_height_m;
        // Mover ring (MBA-1325): additive fields, only when --target-speed > 0. Units are
        // in the field names (mover_ring_m is always meters, MBA-1315 hygiene) regardless
        // of --units, matching the native CLI's `-o json` per-point contract.
        let ring_enabled = target_speed_mps > 0.0;
        // McCoy coordinate system: X=downrange, Y=vertical, Z=lateral
        let points: Vec<serde_json::Value> = result
            .points
            .iter()
            .map(|p| {
                let mut point = match units {
                    UnitSystem::Imperial => {
                        serde_json::json!({
                            "range_yards": p.position.x * 1.09361,  // X is downrange (McCoy)
                            "drop_inches": (los_height - p.position.y) * 39.3701,
                            "drift_inches": p.position.z * 39.3701,  // Z is lateral (windage, McCoy)
                            "velocity_fps": p.velocity_magnitude * 3.28084,
                            "energy_ftlb": p.kinetic_energy * 0.737562149,
                            "time_seconds": p.time
                        })
                    }
                    UnitSystem::Metric => {
                        serde_json::json!({
                            "range_meters": p.position.x,  // X is downrange (McCoy)
                            "drop_cm": (los_height - p.position.y) * 100.0,
                            "drift_cm": p.position.z * 100.0,  // Z is lateral (windage, McCoy)
                            "velocity_mps": p.velocity_magnitude,
                            "energy_joules": p.kinetic_energy,
                            "time_seconds": p.time
                        })
                    }
                };
                if ring_enabled {
                    let (ring_m, ring_mil) =
                        mover_ring(target_speed_mps, p.time, p.position.x);
                    if let Some(obj) = point.as_object_mut() {
                        obj.insert("mover_ring_m".to_string(), serde_json::json!(ring_m));
                        if let Some(mil) = ring_mil {
                            obj.insert("mover_ring_mil".to_string(), serde_json::json!(mil));
                        }
                    }
                }
                point
            })
            .collect();

        let summary = match units {
            UnitSystem::Imperial => {
                serde_json::json!({
                    "max_range_yards": result.max_range * 1.09361,
                    "max_height_inches": result.max_height * 39.3701,
                    "time_of_flight_seconds": result.time_of_flight,
                    "impact_velocity_fps": result.impact_velocity * 3.28084
                })
            }
            UnitSystem::Metric => {
                serde_json::json!({
                    "max_range_meters": result.max_range,
                    "max_height_cm": result.max_height * 100.0,
                    "time_of_flight_seconds": result.time_of_flight,
                    "impact_velocity_mps": result.impact_velocity
                })
            }
        };

        let output = serde_json::json!({
            "trajectory": points,
            "summary": summary,
            "legend": trajectory_json_legend(units),
        });

        serde_json::to_string_pretty(&output)
            .unwrap_or_else(|_| "Error formatting JSON".to_string())
    }

    fn format_trajectory_csv(
        &self,
        result: &crate::cli_api::TrajectoryResult,
        units: UnitSystem,
        full: bool,
        los_height_m: f64,
        target_speed_mps: f64,
    ) -> String {
        let mut output = String::new();
        // Mover ring (MBA-1325): extra column, header carries the unit; only emitted
        // when --target-speed enabled it (additive, matches table/JSON).
        let ring_enabled = target_speed_mps > 0.0;

        // Header
        match units {
            UnitSystem::Imperial => {
                output.push_str("Range(yards),Drop(inches),Drift(inches),Velocity(fps),Energy(ft-lb),Time(seconds)");
            }
            UnitSystem::Metric => {
                output.push_str(
                    "Range(meters),Drop(cm),Drift(cm),Velocity(m/s),Energy(joules),Time(seconds)",
                );
            }
        }
        if ring_enabled {
            output.push_str(",Ring(mil)\n");
        } else {
            output.push('\n');
        }

        // Determine sampling interval
        let max_range_display = match units {
            UnitSystem::Imperial => result.max_range * 1.09361,
            UnitSystem::Metric => result.max_range,
        };

        let sample_interval = if full {
            if max_range_display < 100.0 {
                10.0
            } else {
                25.0
            }
        } else {
            if max_range_display < 500.0 {
                50.0
            } else {
                100.0
            }
        };

        let mut current_range = 0.0;
        // LOS height is cant-invariant (see format_trajectory_table).
        let los_height = los_height_m;

        for (idx, point) in result.points.iter().enumerate() {
            // McCoy coordinate system: X=downrange, Y=vertical, Z=lateral
            let range_display = match units {
                UnitSystem::Imperial => point.position.x * 1.09361, // X is downrange (McCoy)
                UnitSystem::Metric => point.position.x,             // X is downrange (McCoy)
            };

            let is_last_point = idx == result.points.len() - 1;

            if range_display >= current_range || is_last_point {
                let drop = los_height - point.position.y;

                let row = match units {
                    UnitSystem::Imperial => {
                        let energy_ftlb = point.kinetic_energy * 0.737562149;
                        format!(
                            "{:.1},{:.2},{:.2},{:.0},{:.0},{:.3}",
                            range_display,
                            drop * 39.3701,
                            point.position.z * 39.3701, // Z is lateral (windage, McCoy)
                            point.velocity_magnitude * 3.28084,
                            energy_ftlb,
                            point.time
                        )
                    }
                    UnitSystem::Metric => format!(
                        "{:.1},{:.2},{:.2},{:.0},{:.0},{:.3}",
                        range_display,
                        drop * 100.0,
                        point.position.z * 100.0, // Z is lateral (windage, McCoy)
                        point.velocity_magnitude,
                        point.kinetic_energy,
                        point.time
                    ),
                };
                if ring_enabled {
                    // Downrange is position.x (McCoy frame) regardless of the CSV's
                    // lateral/downrange column swap above.
                    let (_, ring_mil) = mover_ring(
                        target_speed_mps,
                        point.time,
                        point.position.x,
                    );
                    match ring_mil {
                        Some(mil) => output.push_str(&format!("{},{:.3}\n", row, mil)),
                        None => output.push_str(&format!("{},\n", row)),
                    }
                } else {
                    output.push_str(&format!("{}\n", row));
                }

                if range_display >= current_range {
                    current_range += sample_interval;
                }
            }
        }

        output
    }

    /// MBA-737: standalone powder-temperature velocity resolution — native `powder`
    /// command parity. The physics rides the same shared
    /// `resolve_powder_adjusted_velocity` the trajectory/lead solvers call. NOTE:
    /// this command always applies the linear model, while trajectory/lead gate it
    /// behind --use-powder-sensitivity (a curve applies there unconditionally).
    /// Output (table/json/csv) is returned as the terminal string.
    fn handle_powder_command(&self, args: &[&str], units: UnitSystem) -> Result<String, JsValue> {
        use crate::cli_api::{parse_powder_sweep, resolve_powder_adjusted_velocity};

        let mut velocity: Option<f64> = None;
        let mut temperature: Option<f64> = None;
        let mut powder_temp_sensitivity: Option<f64> = None;
        let mut powder_temp: Option<f64> = None;
        let mut powder_temp_curve_str: Option<String> = None;
        let mut sweep_str: Option<String> = None;
        let mut mass: Option<f64> = None;
        let mut out_fmt = "table";

        // Parse arguments
        let mut i = 0;
        while i < args.len() {
            match args[i] {
                "-v" | "--velocity" => {
                    if i + 1 < args.len() {
                        velocity = Some(
                            args[i + 1]
                                .parse()
                                .map_err(|_| JsValue::from_str("Invalid velocity"))?,
                        );
                        i += 1;
                    }
                }
                "--temperature" => {
                    if i + 1 < args.len() {
                        temperature = Some(
                            args[i + 1]
                                .parse()
                                .map_err(|_| JsValue::from_str("Invalid temperature"))?,
                        );
                        i += 1;
                    }
                }
                "--powder-temp-sensitivity" => {
                    if i + 1 < args.len() {
                        powder_temp_sensitivity = Some(args[i + 1].parse().map_err(|_| {
                            JsValue::from_str("Invalid powder temp sensitivity")
                        })?);
                        i += 1;
                    }
                }
                "--powder-temp" => {
                    if i + 1 < args.len() {
                        powder_temp = Some(
                            args[i + 1]
                                .parse()
                                .map_err(|_| JsValue::from_str("Invalid powder temperature"))?,
                        );
                        i += 1;
                    }
                }
                "--powder-temp-curve" => {
                    if i + 1 < args.len() {
                        powder_temp_curve_str = Some(args[i + 1].to_string());
                        i += 1;
                    }
                }
                "--sweep" => {
                    if i + 1 < args.len() {
                        sweep_str = Some(args[i + 1].to_string());
                        i += 1;
                    }
                }
                "-m" | "--mass" => {
                    if i + 1 < args.len() {
                        mass = Some(
                            args[i + 1]
                                .parse()
                                .map_err(|_| JsValue::from_str("Invalid mass"))?,
                        );
                        i += 1;
                    }
                }
                "-o" | "--output" => {
                    if i + 1 < args.len() {
                        out_fmt = args[i + 1];
                        i += 1;
                    }
                }
                _ => {}
            }
            i += 1;
        }
        if !matches!(out_fmt, "table" | "json" | "csv") {
            return Err(JsValue::from_str(
                "Invalid output format for powder (expected table, json, or csv)",
            ));
        }
        // Same numeric ranges the native clap definition enforces (f64_range).
        if let Some(v) = velocity {
            if !(0.0..=6000.0).contains(&v) {
                return Err(JsValue::from_str("Velocity must be between 0 and 6000"));
            }
        }
        if let Some(m) = mass {
            if !(0.1..=2000.0).contains(&m) {
                return Err(JsValue::from_str("Mass must be between 0.1 and 2000"));
            }
        }

        // SI conversions — identical factors to the native handler.
        let powder_temp_curve_si: Option<Vec<(f64, f64)>> = match powder_temp_curve_str.as_deref()
        {
            Some(s) => Some(parse_powder_temp_curve_str(s, units)?),
            None => None,
        };
        let has_curve = powder_temp_curve_si.as_ref().is_some_and(|c| !c.is_empty());
        let curve_ref = powder_temp_curve_si.as_deref();

        let to_temp_c = |t: f64| match units {
            UnitSystem::Imperial => (t - 32.0) * 5.0 / 9.0,
            UnitSystem::Metric => t,
        };
        let from_temp_c = |t: f64| match units {
            UnitSystem::Imperial => t * 9.0 / 5.0 + 32.0,
            UnitSystem::Metric => t,
        };
        let to_vel_si = |v: f64| match units {
            UnitSystem::Imperial => v * 0.3048,
            UnitSystem::Metric => v,
        };
        let from_vel_si = |v: f64| match units {
            UnitSystem::Imperial => v / 0.3048,
            UnitSystem::Metric => v,
        };

        let sens_display = powder_temp_sensitivity.unwrap_or(match units {
            UnitSystem::Imperial => 1.0,
            UnitSystem::Metric => 0.3048 / (5.0 / 9.0),
        });
        // fps/degF -> m/s/degC (imperial): velocity factor over the delta-degree factor.
        let sens_si = match units {
            UnitSystem::Imperial => sens_display * 0.3048 / (5.0 / 9.0),
            UnitSystem::Metric => sens_display,
        };
        let ref_temp_c = powder_temp
            .map(to_temp_c)
            .unwrap_or(crate::constants::DEFAULT_POWDER_REFERENCE_TEMP_C);
        let velocity_si = velocity.map(to_vel_si);
        let mass_kg: Option<f64> = mass.map(|m| match units {
            UnitSystem::Imperial => m * crate::constants::GRAINS_TO_KG,
            UnitSystem::Metric => m * 0.001,
        });

        // Validation — same rules and messages as the native handler.
        if !has_curve && velocity.is_none() {
            return Err(JsValue::from_str(
                "--velocity is required for the linear model (or supply --powder-temp-curve)",
            ));
        }
        let sweep_temps: Option<Vec<f64>> = match sweep_str.as_deref() {
            Some(s) => Some(parse_powder_sweep(s).map_err(|e| JsValue::from_str(&e))?),
            None => None,
        };
        if sweep_temps.is_none() {
            if !has_curve && temperature.is_none() {
                return Err(JsValue::from_str(
                    "--temperature is required (the shot-day air temperature), or use --sweep",
                ));
            }
            if has_curve && powder_temp.is_none() && temperature.is_none() {
                return Err(JsValue::from_str(
                    "--powder-temp or --temperature is required with --powder-temp-curve, or use --sweep",
                ));
            }
        }

        let resolve_at_c = |temp_c: f64| -> f64 {
            resolve_powder_adjusted_velocity(
                velocity_si.unwrap_or(0.0),
                temp_c,
                !has_curve,
                sens_si,
                ref_temp_c,
                curve_ref,
                None,
            )
        };
        let energy_at = |v_mps: f64| -> Option<f64> {
            mass_kg.map(|m| {
                let e_j = 0.5 * m * v_mps * v_mps;
                match units {
                    UnitSystem::Imperial => e_j * 0.737562, // Joules to ft-lbs
                    UnitSystem::Metric => e_j,
                }
            })
        };
        let (vel_unit, temp_unit, energy_unit) = match units {
            UnitSystem::Imperial => ("fps", "°F", "ft·lb"),
            UnitSystem::Metric => ("m/s", "°C", "J"),
        };

        struct PowderRow {
            temp_display: f64,
            velocity_display: f64,
            shift_display: Option<f64>,
            energy_display: Option<f64>,
        }
        let row_at = |temp_display: f64, temp_c: f64| -> PowderRow {
            let v_mps = resolve_at_c(temp_c);
            let v_display = from_vel_si(v_mps);
            PowderRow {
                temp_display,
                velocity_display: v_display,
                shift_display: velocity.map(|nominal| v_display - nominal),
                energy_display: energy_at(v_mps),
            }
        };

        let rows: Vec<PowderRow> = match &sweep_temps {
            Some(temps) => temps.iter().map(|&t| row_at(t, to_temp_c(t))).collect(),
            None => {
                let (t_display, t_c) = if has_curve {
                    match (powder_temp, temperature) {
                        (Some(td), _) => (td, to_temp_c(td)),
                        (None, Some(td)) => (td, to_temp_c(td)),
                        _ => unreachable!("validated above"),
                    }
                } else {
                    let td = temperature.expect("validated above");
                    (td, to_temp_c(td))
                };
                vec![row_at(t_display, t_c)]
            }
        };

        let model_desc = if has_curve {
            let curve = curve_ref.expect("has_curve");
            format!(
                "measured curve, {} points ({:.1}\u{2013}{:.1} {})",
                curve.len(),
                from_temp_c(curve[0].0),
                from_temp_c(curve[curve.len() - 1].0),
                temp_unit
            )
        } else {
            format!("linear, {:.2} {}/{}", sens_display, vel_unit, temp_unit)
        };

        let mut out = String::new();
        match out_fmt {
            "json" => {
                let round1 = |x: f64| (x * 10.0).round() / 10.0;
                let rows_json: Vec<serde_json::Value> = rows
                    .iter()
                    .map(|r| {
                        let mut o = serde_json::json!({
                            "temperature": round1(r.temp_display),
                            "velocity": round1(r.velocity_display),
                        });
                        if let Some(s) = r.shift_display {
                            o["shift"] = serde_json::json!(round1(s));
                        }
                        if let Some(e) = r.energy_display {
                            o["energy"] = serde_json::json!(e.round());
                        }
                        o
                    })
                    .collect();
                let mut result = serde_json::json!({
                    "command": "powder",
                    "units": match units {
                        UnitSystem::Imperial => "imperial",
                        UnitSystem::Metric => "metric",
                    },
                    "model": if has_curve { "curve" } else { "linear" },
                });
                if !has_curve {
                    result["sensitivity"] = serde_json::json!(sens_display);
                    result["reference_temp"] = serde_json::json!(round1(from_temp_c(ref_temp_c)));
                } else {
                    result["curve_points"] =
                        serde_json::json!(curve_ref.expect("has_curve").len());
                }
                if let Some(v) = velocity {
                    result["nominal_velocity"] = serde_json::json!(v);
                }
                if let Some(m) = mass {
                    result["mass"] = serde_json::json!(m);
                }
                if sweep_temps.is_some() {
                    result["sweep"] = serde_json::json!(rows_json);
                } else {
                    result["resolved"] = rows_json
                        .into_iter()
                        .next()
                        .unwrap_or(serde_json::json!(null));
                }
                out.push_str(
                    &serde_json::to_string_pretty(&result)
                        .map_err(|e| JsValue::from_str(&e.to_string()))?,
                );
                out.push('\n');
            }
            "csv" => {
                let (temp_a, vel_a, energy_a) = match units {
                    UnitSystem::Imperial => ("f", "fps", "ftlb"),
                    UnitSystem::Metric => ("c", "mps", "j"),
                };
                let has_shift = rows.first().is_some_and(|r| r.shift_display.is_some());
                let has_energy = rows.first().is_some_and(|r| r.energy_display.is_some());
                out.push_str(&format!("temperature_{},velocity_{}", temp_a, vel_a));
                if has_shift {
                    out.push_str(&format!(",shift_{}", vel_a));
                }
                if has_energy {
                    out.push_str(&format!(",energy_{}", energy_a));
                }
                out.push('\n');
                for r in &rows {
                    out.push_str(&format!("{:.1},{:.1}", r.temp_display, r.velocity_display));
                    if let Some(s) = r.shift_display {
                        out.push_str(&format!(",{:.1}", s));
                    }
                    if let Some(e) = r.energy_display {
                        out.push_str(&format!(",{:.0}", e));
                    }
                    out.push('\n');
                }
            }
            _ => {
                out.push_str("Powder Temperature Velocity\n");
                out.push_str("===========================\n");
                out.push_str(&format!("  Model:              {}\n", model_desc));
                if !has_curve {
                    out.push_str(&format!(
                        "  Reference temp:     {:.1} {}\n",
                        from_temp_c(ref_temp_c),
                        temp_unit
                    ));
                }
                if let Some(v) = velocity {
                    out.push_str(&format!("  Nominal velocity:   {:.1} {}\n", v, vel_unit));
                }
                out.push('\n');
                if rows.len() == 1 {
                    let r = &rows[0];
                    let temp_label = if has_curve { "Powder temp:" } else { "Shot temp:" };
                    out.push_str(&format!(
                        "  {:<20}{:.1} {}\n",
                        temp_label, r.temp_display, temp_unit
                    ));
                    match r.shift_display {
                        Some(s) => out.push_str(&format!(
                            "  Resolved velocity:  {:.1} {}  ({:+.1})\n",
                            r.velocity_display, vel_unit, s
                        )),
                        None => out.push_str(&format!(
                            "  Resolved velocity:  {:.1} {}\n",
                            r.velocity_display, vel_unit
                        )),
                    }
                    if let Some(e) = r.energy_display {
                        out.push_str(&format!("  Muzzle energy:      {:.0} {}\n", e, energy_unit));
                    }
                } else {
                    let has_shift = rows.first().is_some_and(|r| r.shift_display.is_some());
                    let has_energy = rows.first().is_some_and(|r| r.energy_display.is_some());
                    out.push_str(&format!(
                        "  {:>10}  {:>14}",
                        format!("Temp ({})", temp_unit),
                        format!("Velocity ({})", vel_unit)
                    ));
                    if has_shift {
                        out.push_str(&format!("  {:>12}", format!("Shift ({})", vel_unit)));
                    }
                    if has_energy {
                        out.push_str(&format!("  {:>15}", format!("Energy ({})", energy_unit)));
                    }
                    out.push('\n');
                    for r in &rows {
                        out.push_str(&format!(
                            "  {:>10.1}  {:>14.1}",
                            r.temp_display, r.velocity_display
                        ));
                        if let Some(s) = r.shift_display {
                            out.push_str(&format!("  {:>12.1}", s));
                        }
                        if let Some(e) = r.energy_display {
                            out.push_str(&format!("  {:>15.0}", e));
                        }
                        out.push('\n');
                    }
                }
            }
        }
        Ok(out)
    }

    fn show_help(&self) -> String {
        r#"Ballistics Engine - WebAssembly Version

Usage: ballistics [OPTIONS] <COMMAND>

Commands:
  trajectory      Calculate ballistic trajectory
  zero           Calculate sight adjustment for zero
  monte-carlo    Run Monte Carlo simulation
  true-velocity  Calculate effective muzzle velocity from observed drop
  estimate-bc    Estimate BC from trajectory data
  lead           Calculate moving-target lead (hold)
  powder         Resolve powder-temperature velocity shift (no trajectory)
  help           Show this help message

Global Options:
  -u, --units <SYSTEM>  Unit system (imperial/metric) [default: imperial]

Host API (JavaScript, not a CLI flag):
  loadDragTable(bytes)   Load a custom Mach:Cd drag table (CSV, same format as native
                         --drag-table). Also accepts .drg vendor drag-curve text as a fallback
                         when the bytes aren't valid CSV. Once loaded, EVERY
                         trajectory/zero/lead/monte-carlo run applies it automatically (no
                         --use-* flag needed); -b/--bc is still accepted but ignored for drag
                         while a table is active.
  hasDragTable()         Report whether a drag table is currently loaded.
  clearDragTable()       Unload the drag table, reverting to the selected --drag-model +
                         BC drag. Lets one instance compare CDM vs G7: load -> solve ->
                         clear -> solve.
  loadBc5dTable(bytes)   Load a BC5D correction table (see --use-bc-segments below).
  hasBc5dTable()         Report whether a BC5D table is currently loaded.

Trajectory Command:
  ballistics trajectory [OPTIONS]
  
  Basic Options:
    -v, --velocity <VEL>         Muzzle velocity (fps/m/s)
    -b, --bc <BC>                Ballistic coefficient
    -m, --mass <MASS>            Mass (grains/grams)
    -d, --diameter <DIA>         Diameter (inches/mm)
    -a, --angle <ANGLE>          Launch angle (degrees)
    --drag-model <MODEL>         Drag model (G1/G2/G5/G6/G7/G8/GI/GS/RA4)
    --cd-scale <FACTOR>          Whole-curve drag scale for a loaded drag table
                                 (1.0 = neutral, typical 0.90-1.10). Requires a
                                 drag table (loadDragTable)
    --max-range <RANGE>          Maximum range (yards/meters)
    -z, --auto-zero <DIST>       Auto-zero at distance
    -o, --output <FORMAT>        Output format (table/json/csv)
    --full                       Show all trajectory points
    --plot [STYLE]               Append terminal charts after the table (drop and
                                 drift vs range); STYLE = braille (default) or ascii
                                 (for fonts without braille glyph coverage)
    --target-speed <SPEED>       Mover ring: per-point ring radius = speed x ToF
                                 (mph/m/s, 0-300; 0 = off). Adds a Ring column (table),
                                 mover_ring_m/mover_ring_mil (json), ring_mil (csv)
    --adjustment-unit <UNIT>     Ring table column unit: mil/moa/smoa/iphy/clicks
                                 [default: mil]
    --elevation-click-value <S>  Turret elevation click graduation, e.g. 0.25moa or
                                 0.1mil; required (once) when --adjustment-unit clicks
    --windage-click-value <S>    Turret windage click graduation (accepted for CLI
                                 parity; the Ring column always uses the elevation one)

  Environmental:
    --wind-speed <SPEED>         Wind speed (mph/m/s)
    --wind-direction <DIR>       Wind direction (deg; 0=headwind, 90=from right)
    --wind-vertical <SPEED>      Vertical wind (mph/m/s); positive = updraft (raises POI)
    --wind-segment <S:A:D[:V]>   Downrange wind seg speed:angle:until-dist[:vertical] (repeatable).
                                 Optional 4th field is ALWAYS m/s updraft-positive, unlike
                                 speed which follows --units
    --temperature <TEMP>         Temperature (F/C)
    --pressure <PRESSURE>        Pressure (inHg/hPa)
    --humidity <HUMIDITY>        Humidity (0-100%)
    --altitude <ALT>             Altitude (feet/meters)
    
  Advanced Physics:
    --enable-magnus              Enable Magnus effect
    --enable-coriolis            Enable Coriolis effect
    --enable-spin-drift          Enable empirical Litz spin drift
    --enable-wind-shear          Enable altitude-dependent wind
    --enable-pitch-damping       Enable transonic stability
    --enable-precession          Enable angular motion physics
    --use-euler                  Use Euler integration (less accurate)
    --use-rk4-fixed              Use fixed-step RK4 (default: adaptive RK45)
    --time-step <SECONDS>        Integration time step (seconds) [default: 0.001]
    --sample-trajectory          Enable trajectory sampling
    --sample-interval <DIST>     Trajectory sampling interval (yards/meters) [default: 10]
    --use-bc-segments            Use velocity-based BC (from a loaded BC5D table)
    --bc-segment <VMIN:VMAX:BC>  Manual velocity-keyed BC segment (repeatable; fps/m/s per --units)
    --print-bc-segments          Print the active BC segment ladder (velocity/Mach span + BC)
                                 applied to this run, or a note when none are active
    --dsf-point <MACH:DSF>       Drop-scale-factor truing point (repeatable; up to 6). Scales
                                 predicted drop by DSF at MACH (0 < MACH < 1.2, 0.5 < DSF < 2.0)
                                 — WASM has no saved profile, so pass the table per call
                                 (native CLI equivalent: a saved profile's `dsf` table)
    --use-powder-sensitivity     Enable powder temp sensitivity
    --ignore-ground-impact       Disable ground-impact truncation; trajectory runs to
                                 --max-range regardless of drop below the muzzle

  Additional Parameters:
    --twist-rate <RATE>          Barrel twist (inches/turn imperial, mm/turn metric)
    --twist-right <BOOL>         Right-hand twist (true/false)
    --latitude <LAT>             Latitude for Coriolis (degrees)
    --shot-direction <DEG>       Compass bearing of the shot for Coriolis (0=N, 90=E)
    --shooting-angle <ANGLE>     Uphill/downhill angle (degrees)
    --cant <DEGREES>             Rifle cant angle (degrees); positive = clockwise from the
                                 shooter, moving point of impact right and low
    --sight-height <HEIGHT>      Sight height above bore (inches/mm)
    --muzzle-height <HEIGHT>     Shooter height above ground (inches/mm; also accepts
                                 --bore-height, matching the native CLI). A value above
                                 ~39,370in/1,000,000mm (1000m) triggers a warning: it feeds
                                 air density, thinning it over the whole flight — use
                                 --altitude for site elevation, --ignore-ground-impact to
                                 defeat ground truncation
    --target-height <HEIGHT>     Target height above ground (inches/mm)
    --powder-temp <TEMP>         Powder temperature
    --powder-temp-sensitivity <SENS>  Velocity change per degree

Zero Command:
  ballistics zero [OPTIONS]
  
  Options:
    -v, --velocity <VEL>         Muzzle velocity
    -b, --bc <BC>                Ballistic coefficient
    -m, --mass <MASS>            Mass
    -d, --diameter <DIA>         Diameter
    --drag-model <MODEL>         Drag model (G1/G2/G5/G6/G7/G8/GI/GS/RA4)
    --cd-scale <FACTOR>          Whole-curve drag scale for a loaded drag table
                                 (1.0 = neutral, typical 0.90-1.10). Requires a
                                 drag table (loadDragTable)
    --target-distance <DIST>     Target distance for zero
    --sight-height <HEIGHT>      Sight height above bore

Monte Carlo Command:
  ballistics monte-carlo [OPTIONS]

  Options:
    -v, --velocity <VEL>         Base velocity
    -b, --bc <BC>                Base BC
    -m, --mass <MASS>            Mass
    -d, --diameter <DIA>         Diameter
    --drag-model <MODEL>         Drag model (G1/G2/G5/G6/G7/G8/GI/GS/RA4)
    --cd-scale <FACTOR>          Whole-curve drag scale for a loaded drag table
                                 (1.0 = neutral, typical 0.90-1.10). Requires a
                                 drag table (loadDragTable)
    -n, --num-sims <N>           Number of simulations
    --velocity-std <STD>         Velocity std deviation
    --angle-std <STD>            Angle std deviation
    --bc-std <STD>               BC std deviation
    --wind-speed-std <STD>       Wind speed std deviation (--wind-std also accepted)
    --wind-direction-std <STD>   Wind direction std deviation (degrees)

  WEZ (Weapon Employment Zone) sweep:
    --wez                        WEZ sweep mode: hit probability vs range for a
                                 target size instead of a single summary
    -a, --angle <ANGLE>          Launch (elevation) angle in degrees, held for
                                 every sweep step (from `ballistics zero`,
                                 typically) [default: 0]
    --target-size <WxH|R>        WEZ target size: WIDTHxHEIGHT (inches imperial,
                                 cm metric; e.g. 18x30) for a rectangle centered
                                 on the point of aim, or a single number for a
                                 circular radius. Required with --wez
    --wind-call-error <ERR>      Shooter's wind CALL error (mph/m/s): uncertainty
                                 in the shooter's own wind estimate, composed with
                                 the wind std in quadrature [default: 0]
    --wez-start <DIST>           Sweep start range (yd/m) [default: 200]
    --wez-end <DIST>             Sweep end range, inclusive (yd/m) [default: 1000]
    --wez-step <DIST>            Sweep step (yd/m) [default: 100]
    -o, --output <FORMAT>        WEZ output format: summary/statistics/full
                                 (table/csv/json accepted as aliases)
                                 [default: summary]

  Browser note: a WEZ sweep runs num-sims full trajectory solves per range step
  (deterministic but heavy) — prefer -n 300 or fewer for interactive use.

True Velocity Command:
  ballistics true-velocity --range <DIST> --measured-drop <DROP> [OPTIONS]

  Calculates the effective muzzle velocity that reproduces an observed drop.
  With one or more --observed impacts it fits muzzle velocity and (when the
  observations constrain it) ballistic coefficient jointly via the real
  forward model (joint MV+BC calibration).

  Options:
    --measured-drop <DROP>       Measured drop at --range (MIL by default;
                                 follows --drop-unit in multi-observation mode).
                                 Required
    --range <DIST>               Range where drop was measured (yd/m). Required
    --observed <RANGE:DROP>      Additional observed impact (repeatable; e.g.
                                 600:5.1). RANGE follows --units, DROP follows
                                 --drop-unit. Enables joint MV+BC calibration
    --drop-unit <UNIT>           Drop unit for --measured-drop/--observed in
                                 multi-observation mode: mil/moa/in [default: mil]
    -b, --bc <BC>                Ballistic coefficient (starting value; fitted
                                 when the observations allow)
    --drag-model <MODEL>         Drag model (G1/G7) [default: g1]
    -m, --mass <MASS>            Bullet weight (grains/grams)
    -d, --diameter <DIA>         Bullet diameter (inches/mm)
    --chrono-velocity <VEL>      Chronograph velocity for comparison (fps/m/s)
    --zero-distance <DIST>       Zero distance (yd/m) [default: 100]
    --sight-height <HEIGHT>      Sight height above bore (in/mm) [default: 2/50]
    --temperature <T>            Temperature (°F/°C) [default: 59/15]
    --pressure <P>               Pressure (inHg/hPa) [default: 29.92/1013.25]
    --humidity <H>               Humidity (0-100%) [default: 50]
    --altitude <A>               Altitude (ft/m) [default: 0]
    --offline                    Accepted for native-command parity (the WASM
                                 terminal always calculates locally)
    -o, --output <FORMAT>        Output format (table/json/csv) [default: table]

Estimate BC Command:
  ballistics estimate-bc [OPTIONS]

  Options:
    -v, --velocity <VEL>         Muzzle velocity (fps/m/s)
    -m, --mass <MASS>            Mass (grains/grams)
    -d, --diameter <DIA>         Diameter (inches/mm)
    --data <PAIRS>               Drop data: "dist,drop;..." (yd,in / m,mm)
    --velocity-data <PAIRS>      Velocity data: "dist,vel;..." (yd,fps / m,m/s)
    --drag-model <MODEL>         g1, g7, both, or any single family (g2/g5/g6/g8/gi/gs/ra4) [default: both]
    --zero-range <DIST>          Zero range of the drop data (yd/m). Dope cards are
                                 zeroed — pass this so drop is fit below line of sight.
    --sight-height <H>           Sight height above bore (in/mm) for the zeroed fit
    --temperature <T>            Air temp the data was measured at (°F/°C) [59/15]
    --pressure <P>               Pressure the data was measured at (inHg/hPa) [29.92/1013.25]
    --humidity <H>               Relative humidity (percent) [50]
    --altitude <A>               Altitude the data was measured at (ft/m) [0]

  Prints a BC for each drag model x data basis you supply. For a DOPE CARD, pass
  --zero-range (drop is below line of sight) and set the atmosphere it was made at.
  Without --zero-range, drop is treated as bore-referenced (flat-fire). --velocity-data
  gives a velocity-retention fit (immune to zero/angle). A fit that can't be pinned
  down is flagged UNRELIABLE.

Lead Command:
  ballistics lead --target-speed <SPEED> [OPTIONS]

  Options:
    -v, --velocity <VEL>          Muzzle velocity (fps/m/s)
    -b, --bc <BC>                 Ballistic coefficient
    -m, --mass <MASS>             Mass (grains/grams)
    -d, --diameter <DIA>          Diameter (inches/mm)
    --drag-model <MODEL>          Drag model (G1/G2/G5/G6/G7/G8/GI/GS/RA4)
    --cd-scale <FACTOR>            Whole-curve drag scale for a loaded drag table
                                  (1.0 = neutral, typical 0.90-1.10). Requires a
                                  drag table (loadDragTable)
    --sight-height <HEIGHT>       Sight height above bore (inches/mm)
    --temperature <T>             Air temperature (°F/°C) [default: 15°C standard]
    --pressure <P>                Barometric pressure (inHg/hPa) [default: 1013.25 hPa]
    --humidity <H>                Relative humidity (percent) [default: 50]
    --altitude <A>                Altitude (ft/m) [default: 0]
    --wind-speed <SPEED>          Wind speed (mph/m/s) [default: 0]
    --wind-direction <DEG>        Wind direction, degrees, wind-FROM [default: 0]
                                  0=headwind, 90=from right, 180=tailwind, 270=from left
    --use-powder-sensitivity      Enable linear powder temperature sensitivity
    --powder-temp-sensitivity <S> Sensitivity (fps/°F or m/s/°C) [default: 1.0 fps/°F]
    --powder-temp <T>             Powder temperature (°F/°C); curve lookup temp, or
                                  linear-model reference temp [default: 70°F/21°C]
    --powder-temp-curve <T:V,...> Measured powder-temp->velocity curve (overrides
                                  --powder-temp-sensitivity)
    --target-speed <SPEED>        Target speed (mph/m/s) [required]
    --target-angle <DEG>          Direction of target travel, degrees [default: 90]
                                  0=away, 90=left-to-right, 180=toward, 270=right-to-left;
                                  positive lead = hold in direction of travel
    --range <DIST>                Range to target (yards/meters) [default: 500]
    --adjustment-unit <UNIT>      mil/moa/smoa/iphy [default: mil] (clicks: trajectory
                                  and come-ups only, MBA-1355)
    -o, --output <FORMAT>         Output format (table/json) [default: table]

  Time of flight is solved under the supplied wind/atmosphere (wind-aware lead);
  the lead figure itself is pure target motion — wind hold stays separate.

Powder Command:
  ballistics powder [OPTIONS]

  Resolves the powder-temperature-adjusted muzzle velocity without running a
  trajectory. A curve carries to trajectory/lead as-is; for the linear model those
  commands additionally need --use-powder-sensitivity (powder always applies it).

  Options:
    -v, --velocity <VEL>          Nominal muzzle velocity (fps/m/s) at the reference
                                  temperature. Required for the linear model; with a
                                  curve it only anchors the reported shift.
    --temperature <T>             Shot-day air temperature (°F/°C); curve lookup
                                  fallback when --powder-temp is unset
    --powder-temp-sensitivity <S> Sensitivity (fps/°F or m/s/°C) [default: 1.0 fps/°F]
    --powder-temp <T>             Powder temperature (°F/°C); curve lookup temp
                                  [default: --temperature], or linear-model reference
                                  temp [default: 70°F/21°C]
    --powder-temp-curve <T:V,...> Measured powder-temp->velocity curve (overrides
                                  --powder-temp-sensitivity; clamped at endpoints)
    --sweep <START:END:STEP>      Velocity table across temperatures (°F/°C)
    -m, --mass <MASS>             Bullet mass (grains/grams): adds muzzle energy
    -o, --output <FORMAT>         Output format (table/json/csv) [default: table]

Examples:
  ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308
  ballistics trajectory --auto-zero 200 --enable-spin-drift
  ballistics --units metric trajectory -v 823 -b 0.475 -m 10.9
  ballistics zero --target-distance 300
  ballistics estimate-bc -v 2700 -m 168 -d 0.308 --data "100,2.1;200,9.4;300,22.8"
  ballistics estimate-bc -v 2650 -m 77 -d 0.224 --data "300,29;500,89.9" \
    --velocity-data "300,1980;500,1560" --drag-model both
  ballistics monte-carlo -n 1000 --velocity-std 10
  ballistics monte-carlo -v 2700 -b 0.475 -m 168 -d 0.308 --wez \
    --target-size 18x30 -n 300 --wez-start 200 --wez-end 500 --wez-step 100
  ballistics true-velocity --range 300 --measured-drop 1.8 -b 0.475 -m 168 -d 0.308
  ballistics true-velocity --range 300 --measured-drop 1.8 --observed 600:5.1 \
    --observed 900:10.9 -b 0.475 -m 168 -d 0.308 --chrono-velocity 2700"#
            .to_string()
    }
}

// ============================================================================
// Object-Oriented API
// ============================================================================

/// Object-oriented calculator for programmatic use
/// Provides a type-safe, fluent API alternative to the CLI interface
#[wasm_bindgen]
pub struct Calculator {
    // Projectile properties
    bc: f64,
    velocity_fps: f64,
    mass_grains: f64,
    diameter_inches: f64,
    drag_model: String,

    // Environmental conditions
    wind_speed_mph: f64,
    wind_direction_deg: f64,
    // Downrange wind segments as (speed_mph, angle_deg, until_yards); when non-empty
    // these are emitted as --wind-segment flags and override the scalar wind above.
    wind_segments: Vec<(f64, f64, f64)>,
    temperature_f: f64,
    pressure_inhg: f64,
    humidity_percent: f64,
    altitude_ft: f64,

    // Shooting parameters
    sight_height_inches: f64,
    zero_range_yards: Option<f64>,
    max_range_yards: f64,

    // Advanced options
    enable_spin_drift: bool,
    enable_coriolis: bool,
    twist_rate_inches: Option<f64>,
    latitude_deg: Option<f64>,
}

#[wasm_bindgen]
impl Calculator {
    /// Create a new calculator with default values
    /// Defaults: .308 Winchester 168gr at 2700 fps, standard atmosphere
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        Calculator {
            bc: 0.475,
            velocity_fps: 2700.0,
            mass_grains: 168.0,
            diameter_inches: 0.308,
            drag_model: "G1".to_string(), // G1 matches the G1-scale default BC (0.475) and the CLI default

            wind_speed_mph: 0.0,
            wind_direction_deg: 0.0, // wind-FROM convention: 0 = headwind, matching the CLI/trajectory defaults
            wind_segments: Vec::new(),
            temperature_f: 59.0,
            pressure_inhg: 29.92,
            humidity_percent: 50.0,
            altitude_ft: 0.0,

            sight_height_inches: 2.0,
            zero_range_yards: None,
            max_range_yards: 1000.0,

            enable_spin_drift: false,
            enable_coriolis: false,
            twist_rate_inches: None,
            latitude_deg: None,
        }
    }

    // Projectile property setters (fluent API)

    #[wasm_bindgen(js_name = setBC)]
    pub fn set_bc(mut self, bc: f64) -> Self {
        self.bc = bc;
        self
    }

    #[wasm_bindgen(js_name = setVelocity)]
    pub fn set_velocity(mut self, velocity_fps: f64) -> Self {
        self.velocity_fps = velocity_fps;
        self
    }

    #[wasm_bindgen(js_name = setMass)]
    pub fn set_mass(mut self, mass_grains: f64) -> Self {
        self.mass_grains = mass_grains;
        self
    }

    #[wasm_bindgen(js_name = setDiameter)]
    pub fn set_diameter(mut self, diameter_inches: f64) -> Self {
        self.diameter_inches = diameter_inches;
        self
    }

    #[wasm_bindgen(js_name = setDragModel)]
    pub fn set_drag_model(mut self, model: &str) -> Self {
        self.drag_model = model.to_string();
        self
    }

    // Environmental setters

    #[wasm_bindgen(js_name = setWind)]
    pub fn set_wind(mut self, speed_mph: f64, direction_deg: f64) -> Self {
        self.wind_speed_mph = speed_mph;
        self.wind_direction_deg = direction_deg;
        self
    }

    /// Add a downrange wind segment: `speed_mph` from `direction_deg`
    /// (wind-FROM: 0 = headwind, 90 = from the right) out to `until_yards`.
    /// Each segment applies from the previous boundary to `until_yards`; wind is
    /// zero beyond the last segment. When any segment is added it overrides the
    /// scalar `setWind` value. Repeatable.
    #[wasm_bindgen(js_name = addWindSegment)]
    pub fn add_wind_segment(
        mut self,
        speed_mph: f64,
        direction_deg: f64,
        until_yards: f64,
    ) -> Self {
        self.wind_segments
            .push((speed_mph, direction_deg, until_yards));
        self
    }

    /// Remove all downrange wind segments (reverts to the scalar `setWind`).
    #[wasm_bindgen(js_name = clearWindSegments)]
    pub fn clear_wind_segments(mut self) -> Self {
        self.wind_segments.clear();
        self
    }

    #[wasm_bindgen(js_name = setTemperature)]
    pub fn set_temperature(mut self, temp_f: f64) -> Self {
        self.temperature_f = temp_f;
        self
    }

    #[wasm_bindgen(js_name = setPressure)]
    pub fn set_pressure(mut self, pressure_inhg: f64) -> Self {
        self.pressure_inhg = pressure_inhg;
        self
    }

    #[wasm_bindgen(js_name = setHumidity)]
    pub fn set_humidity(mut self, humidity: f64) -> Self {
        self.humidity_percent = humidity;
        self
    }

    #[wasm_bindgen(js_name = setAltitude)]
    pub fn set_altitude(mut self, altitude_ft: f64) -> Self {
        self.altitude_ft = altitude_ft;
        self
    }

    // Shooting parameter setters

    #[wasm_bindgen(js_name = setSightHeight)]
    pub fn set_sight_height(mut self, height_inches: f64) -> Self {
        self.sight_height_inches = height_inches;
        self
    }

    #[wasm_bindgen(js_name = setZeroRange)]
    pub fn set_zero_range(mut self, range_yards: f64) -> Self {
        self.zero_range_yards = Some(range_yards);
        self
    }

    #[wasm_bindgen(js_name = setMaxRange)]
    pub fn set_max_range(mut self, range_yards: f64) -> Self {
        self.max_range_yards = range_yards;
        self
    }

    // Advanced option setters

    #[wasm_bindgen(js_name = enableSpinDrift)]
    pub fn enable_spin_drift_opt(mut self, enabled: bool, twist_rate: Option<f64>) -> Self {
        self.enable_spin_drift = enabled;
        self.twist_rate_inches = twist_rate;
        self
    }

    #[wasm_bindgen(js_name = enableCoriolis)]
    pub fn enable_coriolis_opt(mut self, enabled: bool, latitude: Option<f64>) -> Self {
        self.enable_coriolis = enabled;
        self.latitude_deg = latitude;
        self
    }

    // Calculation method

    /// Calculate trajectory and return result as JavaScript object
    /// Returns: { range_yards, drop_inches, windage_inches, velocity_fps, energy_ftlb, time_sec }
    #[wasm_bindgen(js_name = calculateTrajectory)]
    pub fn calculate_trajectory(&self, range_yards: f64) -> Result<JsValue, JsValue> {
        // Build CLI command from parameters
        let mut cmd = format!(
            "ballistics trajectory -v {} -b {} -m {} -d {} --drag-model {} --max-range {}",
            self.velocity_fps,
            self.bc,
            self.mass_grains,
            self.diameter_inches,
            self.drag_model,
            range_yards
        );

        // Add environmental parameters
        if self.wind_speed_mph > 0.0 {
            cmd.push_str(&format!(
                " --wind-speed {} --wind-direction {}",
                self.wind_speed_mph, self.wind_direction_deg
            ));
        }
        for (speed_mph, direction_deg, until_yards) in &self.wind_segments {
            cmd.push_str(&format!(
                " --wind-segment {}:{}:{}",
                speed_mph, direction_deg, until_yards
            ));
        }

        cmd.push_str(&format!(
            " --temperature {} --pressure {} --humidity {} --altitude {}",
            self.temperature_f, self.pressure_inhg, self.humidity_percent, self.altitude_ft
        ));

        // Add shooting parameters
        cmd.push_str(&format!(" --sight-height {}", self.sight_height_inches));

        if let Some(zero) = self.zero_range_yards {
            cmd.push_str(&format!(" --auto-zero {}", zero));
        }

        // Add advanced options
        if self.enable_spin_drift {
            cmd.push_str(" --enable-spin-drift");
            if let Some(twist) = self.twist_rate_inches {
                cmd.push_str(&format!(" --twist-rate {}", twist));
            }
        }

        if self.enable_coriolis {
            cmd.push_str(" --enable-coriolis");
            if let Some(lat) = self.latitude_deg {
                cmd.push_str(&format!(" --latitude {}", lat));
            }
        }

        // Use JSON output format for easy parsing
        cmd.push_str(" -o json");

        // Execute command
        let wasm_ballistics = WasmBallistics::new();
        let result_str = wasm_ballistics.run_command(&cmd)?;

        // Strip any text before the JSON (like zero info messages)
        let json_str = if let Some(json_start) = result_str.find('{') {
            &result_str[json_start..]
        } else {
            &result_str
        };

        // Parse JSON result
        let json_result: serde_json::Value = serde_json::from_str(json_str).map_err(|e| {
            JsValue::from_str(&format!(
                "JSON parse error: {}. Result was: {}",
                e,
                &json_str[..json_str.len().min(500)]
            ))
        })?;

        // Find the point closest to the requested range
        if let Some(trajectory) = json_result.get("trajectory").and_then(|t| t.as_array()) {
            let target_point = trajectory
                .iter()
                .min_by_key(|p| {
                    let range = p.get("range_yards").and_then(|r| r.as_f64()).unwrap_or(0.0);
                    ((range - range_yards).abs() * 100.0) as i64
                })
                .ok_or_else(|| JsValue::from_str("No trajectory points found"))?;

            // Manually construct JavaScript object to avoid Map conversion
            let result = js_sys::Object::new();
            js_sys::Reflect::set(
                &result,
                &"range_yards".into(),
                &target_point
                    .get("range_yards")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(0.0)
                    .into(),
            )?;
            js_sys::Reflect::set(
                &result,
                &"drop_inches".into(),
                &target_point
                    .get("drop_inches")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(0.0)
                    .into(),
            )?;
            js_sys::Reflect::set(
                &result,
                &"drift_inches".into(),
                &target_point
                    .get("drift_inches")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(0.0)
                    .into(),
            )?;
            js_sys::Reflect::set(
                &result,
                &"velocity_fps".into(),
                &target_point
                    .get("velocity_fps")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(0.0)
                    .into(),
            )?;
            js_sys::Reflect::set(
                &result,
                &"energy_ftlb".into(),
                &target_point
                    .get("energy_ftlb")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(0.0)
                    .into(),
            )?;
            js_sys::Reflect::set(
                &result,
                &"time_seconds".into(),
                &target_point
                    .get("time_seconds")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(0.0)
                    .into(),
            )?;

            Ok(result.into())
        } else {
            Err(JsValue::from_str("Invalid trajectory data"))
        }
    }

    /// Get full trajectory table as array of points
    /// Returns array of: [{ range_yards, drop_inches, windage_inches, velocity_fps, energy_ftlb, time_sec }, ...]
    #[wasm_bindgen(js_name = getFullTrajectory)]
    pub fn get_full_trajectory(&self) -> Result<JsValue, JsValue> {
        // Build CLI command (similar to calculate_trajectory but get full table)
        let mut cmd = format!(
            "ballistics trajectory -v {} -b {} -m {} -d {} --drag-model {} --max-range {} --full",
            self.velocity_fps,
            self.bc,
            self.mass_grains,
            self.diameter_inches,
            self.drag_model,
            self.max_range_yards
        );

        // Add environmental parameters
        if self.wind_speed_mph > 0.0 {
            cmd.push_str(&format!(
                " --wind-speed {} --wind-direction {}",
                self.wind_speed_mph, self.wind_direction_deg
            ));
        }
        for (speed_mph, direction_deg, until_yards) in &self.wind_segments {
            cmd.push_str(&format!(
                " --wind-segment {}:{}:{}",
                speed_mph, direction_deg, until_yards
            ));
        }

        cmd.push_str(&format!(
            " --temperature {} --pressure {} --humidity {} --altitude {}",
            self.temperature_f, self.pressure_inhg, self.humidity_percent, self.altitude_ft
        ));

        cmd.push_str(&format!(" --sight-height {}", self.sight_height_inches));

        if let Some(zero) = self.zero_range_yards {
            cmd.push_str(&format!(" --auto-zero {}", zero));
        }

        if self.enable_spin_drift {
            cmd.push_str(" --enable-spin-drift");
            if let Some(twist) = self.twist_rate_inches {
                cmd.push_str(&format!(" --twist-rate {}", twist));
            }
        }

        if self.enable_coriolis {
            cmd.push_str(" --enable-coriolis");
            if let Some(lat) = self.latitude_deg {
                cmd.push_str(&format!(" --latitude {}", lat));
            }
        }

        cmd.push_str(" -o json");

        // Execute command
        let wasm_ballistics = WasmBallistics::new();
        let result_str = wasm_ballistics.run_command(&cmd)?;

        // Strip any text before the JSON (like zero info messages)
        let json_str = if let Some(json_start) = result_str.find('{') {
            &result_str[json_start..]
        } else {
            &result_str
        };

        // Parse JSON result and return trajectory array
        let json_result: serde_json::Value = serde_json::from_str(json_str).map_err(|e| {
            JsValue::from_str(&format!(
                "JSON parse error: {}. Result: {}",
                e,
                &json_str[..json_str.len().min(500)]
            ))
        })?;

        if let Some(trajectory) = json_result.get("trajectory") {
            // Convert trajectory array to JavaScript array of plain objects
            let js_array = js_sys::Array::new();

            if let Some(points) = trajectory.as_array() {
                for point in points {
                    let js_point = js_sys::Object::new();

                    // Extract each field and add to JavaScript object
                    if let Some(range) = point.get("range_yards").and_then(|v| v.as_f64()) {
                        js_sys::Reflect::set(&js_point, &"range_yards".into(), &range.into())?;
                    }
                    if let Some(drop) = point.get("drop_inches").and_then(|v| v.as_f64()) {
                        js_sys::Reflect::set(&js_point, &"drop_inches".into(), &drop.into())?;
                    }
                    // The JSON producer emits "drift_inches"/"time_seconds"; read those (the old
                    // "windage_inches"/"time_sec" lookups always missed, dropping both fields).
                    // Keep the public output keys (windage_inches/time_sec) unchanged.
                    if let Some(windage) = point.get("drift_inches").and_then(|v| v.as_f64()) {
                        js_sys::Reflect::set(&js_point, &"windage_inches".into(), &windage.into())?;
                    }
                    if let Some(velocity) = point.get("velocity_fps").and_then(|v| v.as_f64()) {
                        js_sys::Reflect::set(&js_point, &"velocity_fps".into(), &velocity.into())?;
                    }
                    if let Some(energy) = point.get("energy_ftlb").and_then(|v| v.as_f64()) {
                        js_sys::Reflect::set(&js_point, &"energy_ftlb".into(), &energy.into())?;
                    }
                    if let Some(time) = point.get("time_seconds").and_then(|v| v.as_f64()) {
                        js_sys::Reflect::set(&js_point, &"time_sec".into(), &time.into())?;
                    }

                    js_array.push(&js_point);
                }
            }

            Ok(js_array.into())
        } else {
            Err(JsValue::from_str("No trajectory data found"))
        }
    }
}


