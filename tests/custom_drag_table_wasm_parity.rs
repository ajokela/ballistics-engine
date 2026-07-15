// MBA-1328: proves the WASM `loadDragTable(bytes)` assembly path is equivalent to the native
// `--drag-table <FILE>` assembly path.
//
// `src/wasm.rs` is `#[cfg(target_arch = "wasm32")]`-gated (see src/lib.rs), so it cannot be
// compiled into or called from a native test binary. Both platforms nonetheless bottom out in
// the SAME shared library code:
//   - native: `load_drag_table_or_exit` (src/main.rs) calls `DragTable::from_file(path)`, which
//     is itself just `std::fs::read_to_string(path)` followed by `DragTable::from_csv_str` (see
//     src/drag.rs) — i.e. "read the file into a String, then parse the CSV text".
//   - WASM: `loadDragTable` (src/wasm.rs) calls `std::str::from_utf8(bytes)` followed by the
//     SAME `DragTable::from_csv_str` — i.e. "UTF-8-decode the bytes into a &str, then parse the
//     CSV text".
// The only structural difference between the two loaders is the I/O source (disk file vs.
// in-memory bytes); the CSV parser and everything downstream (`BallisticInputs::custom_drag_table`
// -> `TrajectorySolver`) is identical code already shared by both platforms. This test drives
// both loaders against byte-identical CSV content and asserts the parsed tables — and full
// trajectory solves built from them — are bit-for-bit identical, which is the concrete claim
// "the WASM assembly path equals the native assembly path" reduces to.
//
// `tests/custom_drag_table_cli.rs` separately drives the compiled CLI binary end-to-end through
// the real `--drag-table` flag; this test complements it by pinning the underlying loader
// equivalence at the library level.

use ballistics_engine::drag::DragTable;
use ballistics_engine::{
    AtmosphericConditions, BallisticInputs, DragModel, TrajectorySolver, WindConditions,
};
use std::io::Write;

/// A representative measured deck (mirrors the worked example in CLI_USAGE.md's Custom Drag
/// Tables section): a transonic bump a plain G1/G7 lookup would not reproduce, so a solve
/// actually exercises the table rather than coincidentally matching a built-in curve.
const DECK_CSV: &str = "mach,cd\n\
0.5,0.220\n\
0.8,0.230\n\
1.0,0.520\n\
1.2,0.480\n\
1.5,0.400\n\
2.0,0.330\n\
2.5,0.300\n";

fn write_temp_csv(name: &str, contents: &str) -> std::path::PathBuf {
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let path =
        std::env::temp_dir().join(format!("mba1328_parity_{name}_{}_{nonce}.csv", std::process::id()));
    std::fs::File::create(&path)
        .unwrap()
        .write_all(contents.as_bytes())
        .unwrap();
    path
}

/// Build `BallisticInputs` exactly as `handle_trajectory_command` (src/wasm.rs) assembles them
/// for imperial units with all defaults (2700 fps / 168 gr / .308 / G1 / --bc 0.475), plus the
/// SAME hook the WASM auto-apply uses to install a loaded table:
/// `inputs.custom_drag_table = Some(table)`. No advanced physics (Magnus/Coriolis/spin-drift) is
/// enabled, so `twist_rate` is left at its inert struct default — it does not feed the drag
/// acceleration this test compares.
fn wasm_style_inputs(table: DragTable) -> BallisticInputs {
    let mut inputs = BallisticInputs::default();
    inputs.muzzle_velocity = 2700.0 * 0.3048; // fps -> m/s
    inputs.bullet_mass = 168.0 * 0.00006479891; // grains -> kg
    inputs.bullet_diameter = 0.308 * 0.0254; // inches -> meters
    inputs.sight_height = 2.0 * 0.0254; // WASM imperial default: 2 in
    inputs.muzzle_height = 60.0 * 0.0254; // WASM imperial default: 60 in (5 ft)
    inputs.target_height = 0.0;
    inputs.bullet_length =
        ballistics_engine::stability::estimate_bullet_length_m(inputs.bullet_diameter, inputs.bullet_mass);
    if inputs.bullet_length <= 0.0 {
        inputs.bullet_length = inputs.bullet_diameter * 4.5;
    }
    inputs.bc_value = 0.475; // WASM default --bc; physically ignored once custom_drag_table is set
    inputs.bc_type = DragModel::G1;
    inputs.muzzle_angle = 0.0;
    inputs.shooting_angle = 0.0;
    inputs.cant_angle = 0.0;
    inputs.ground_threshold = 0.0; // WASM: ignore_ground_impact=false -> 0.0 (not the struct default)
    inputs.temperature = (59.0 - 32.0) * 5.0 / 9.0; // WASM imperial default: 59 F
    inputs.pressure = 29.92 * 33.863886666667; // WASM imperial default: 29.92 inHg
    inputs.humidity = (50.0_f64 / 100.0).clamp(0.0, 1.0); // WASM default: 50%
    inputs.altitude = 0.0;
    inputs.custom_drag_table = Some(table);
    inputs
}

fn atmosphere_for_wasm_style() -> AtmosphericConditions {
    AtmosphericConditions {
        temperature: (59.0 - 32.0) * 5.0 / 9.0,
        pressure: 29.92 * 33.863886666667,
        humidity: 50.0, // AtmosphericConditions::humidity is PERCENT, not the 0-1 fraction
        altitude: 0.0,
    }
}

fn solve(table: DragTable) -> ballistics_engine::TrajectoryResult {
    let inputs = wasm_style_inputs(table);
    let mut solver =
        TrajectorySolver::new(inputs, WindConditions::default(), atmosphere_for_wasm_style());
    solver.set_max_range(300.0 * 0.9144); // 300 yd, matching WASM's yards->meters conversion
    solver.set_time_step(0.001); // WASM default --time-step
    solver.solve().expect("solve")
}

fn assert_results_bit_identical(
    a: &ballistics_engine::TrajectoryResult,
    b: &ballistics_engine::TrajectoryResult,
) {
    assert_eq!(a.points.len(), b.points.len(), "point count diverged");
    assert!(!a.points.is_empty(), "solve produced no points");
    for (i, (pa, pb)) in a.points.iter().zip(b.points.iter()).enumerate() {
        assert_eq!(pa.time.to_bits(), pb.time.to_bits(), "point {i} time diverged");
        assert_eq!(
            pa.velocity_magnitude.to_bits(),
            pb.velocity_magnitude.to_bits(),
            "point {i} velocity_magnitude diverged"
        );
        assert_eq!(
            pa.kinetic_energy.to_bits(),
            pb.kinetic_energy.to_bits(),
            "point {i} kinetic_energy diverged"
        );
        assert_eq!(
            pa.position.x.to_bits(),
            pb.position.x.to_bits(),
            "point {i} position.x diverged"
        );
        assert_eq!(
            pa.position.y.to_bits(),
            pb.position.y.to_bits(),
            "point {i} position.y diverged"
        );
        assert_eq!(
            pa.position.z.to_bits(),
            pb.position.z.to_bits(),
            "point {i} position.z diverged"
        );
    }
    assert_eq!(a.max_range.to_bits(), b.max_range.to_bits());
    assert_eq!(a.max_height.to_bits(), b.max_height.to_bits());
    assert_eq!(a.time_of_flight.to_bits(), b.time_of_flight.to_bits());
    assert_eq!(a.impact_velocity.to_bits(), b.impact_velocity.to_bits());
    assert_eq!(a.impact_energy.to_bits(), b.impact_energy.to_bits());
}

/// The two loaders must parse byte-identical CSV content into byte-identical `DragTable`s.
#[test]
fn wasm_and_native_loaders_parse_identical_tables() {
    let path = write_temp_csv("loaders", DECK_CSV);

    // Native: exactly what `load_drag_table_or_exit` (src/main.rs) does for --drag-table.
    let native_table = DragTable::from_file(&path).expect("native from_file");
    std::fs::remove_file(&path).ok();

    // WASM: exactly what `loadDragTable` (src/wasm.rs) does with the bytes it's handed.
    let bytes = DECK_CSV.as_bytes();
    let csv = std::str::from_utf8(bytes).expect("utf8 decode");
    let wasm_table = DragTable::from_csv_str(csv).expect("wasm from_csv_str");

    assert_eq!(native_table.mach_values, wasm_table.mach_values);
    assert_eq!(native_table.cd_values, wasm_table.cd_values);
}

/// The two loaders, plugged into otherwise-identical `BallisticInputs`, must produce
/// bit-for-bit identical `TrajectorySolver` results — the concrete "WASM assembly path equals
/// native assembly path" claim.
#[test]
fn wasm_and_native_drag_table_assembly_produce_identical_trajectory() {
    let path = write_temp_csv("trajectory", DECK_CSV);

    let native_table = DragTable::from_file(&path).expect("native from_file");
    std::fs::remove_file(&path).ok();

    let bytes = DECK_CSV.as_bytes();
    let csv = std::str::from_utf8(bytes).expect("utf8 decode");
    let wasm_table = DragTable::from_csv_str(csv).expect("wasm from_csv_str");

    let native_result = solve(native_table);
    let wasm_result = solve(wasm_table);

    assert_results_bit_identical(&native_result, &wasm_result);

    // Sanity precondition: DECK_CSV's table must actually be driving the solve (not silently
    // ignored), otherwise the bit-identical assertion above would be vacuous. Compare against a
    // flat, deliberately different low-drag curve.
    let flat_table =
        DragTable::from_csv_str("mach,cd\n0.0,0.10\n1.0,0.10\n2.5,0.10\n").unwrap();
    let flat_result = solve(flat_table);
    assert!(
        (native_result.impact_velocity - flat_result.impact_velocity).abs() > 1.0,
        "runs with different tables should differ materially; DECK_CSV's table may be unused"
    );
}

/// Malformed CSV must error identically through both loaders, citing the same parser message
/// (both call the same `DragTable::from_csv_str`; only the I/O front-end differs).
#[test]
fn wasm_and_native_loaders_reject_malformed_csv_identically() {
    let bad = "0.5,0.23\n1.0,notanumber\n";
    let path = write_temp_csv("bad", bad);

    let native_err = DragTable::from_file(&path).unwrap_err();
    std::fs::remove_file(&path).ok();

    let wasm_err = DragTable::from_csv_str(std::str::from_utf8(bad.as_bytes()).unwrap())
        .unwrap_err();

    assert_eq!(
        native_err, wasm_err,
        "native (from_file) and WASM (from_csv_str) must report the identical parser error"
    );
    assert!(native_err.contains("line 2"), "got: {native_err}");
}
