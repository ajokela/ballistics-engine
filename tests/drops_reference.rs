//! MBA-1403: `--drops-reference {los|target}` — LOS vs target-plane sampled drops.
//!
//! `los` (the default) keeps the historical LOS-perpendicular drop byte-identically;
//! `target` reports drop as vertical in the target plane: the LOS-perpendicular drop
//! divided by cos(shooting_angle) (JBM's "target plane" checkbox), with the sampled
//! table/CSV drop column relabeled. Pins:
//!   (a) absent flag == explicit `--drops-reference los`, byte-identical (table + CSV);
//!   (b) target mode at a 30 degree incline scales every sampled drop by exactly
//!       1/cos(30 deg) (bit-exact division) and touches nothing else;
//!   (c) target mode at a 0 degree angle leaves drop VALUES identical (cos 0 == 1)
//!       while still relabeling the column;
//!   (d) table/CSV labels carry the active reference;
//!   (e) solve-json wire: the field decodes (not unknown_field), scales `drop_m` by
//!       exactly 1/cos, leaves `windage_m`/summary untouched, and omitted == explicit
//!       "los" byte-identically; an invalid value is a structured error;
//!   (f) target mode at |shooting angle| >= 90 degrees is rejected.

use ballistics_engine::{
    AtmosphericConditions, BallisticInputs, DropsReference, TrajectorySolver, WindConditions,
};
use std::process::Command;

const BIN: &str = env!("CARGO_BIN_EXE_ballistics");

const BASE_ARGS: &[&str] = &[
    "trajectory",
    "--velocity",
    "2700",
    "--bc",
    "0.475",
    "--mass",
    "168",
    "--diameter",
    "0.308",
    "--max-range",
    "600",
    "--auto-zero",
    "100",
    "--sample-trajectory",
    "--sample-interval",
    "100",
];

fn run_ok(args: &[&str]) -> std::process::Output {
    let output = Command::new(BIN).args(args).output().expect("run binary");
    assert!(
        output.status.success(),
        "command failed: {:?}\nstderr: {}",
        args,
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn library_inputs(shooting_angle_rad: f64, reference: DropsReference) -> BallisticInputs {
    BallisticInputs {
        bc_value: 0.475,
        muzzle_velocity: 823.0,
        bullet_mass: 0.01089,
        bullet_diameter: 0.007823,
        sight_height: 0.05,
        shooting_angle: shooting_angle_rad,
        enable_trajectory_sampling: true,
        sample_interval: 100.0,
        drops_reference: reference,
        ..Default::default()
    }
}

fn sampled(
    shooting_angle_rad: f64,
    reference: DropsReference,
) -> Vec<ballistics_engine::trajectory_sampling::TrajectorySample> {
    let mut solver = TrajectorySolver::new(
        library_inputs(shooting_angle_rad, reference),
        WindConditions::default(),
        AtmosphericConditions::default(),
    );
    solver.set_max_range(600.0);
    solver
        .solve()
        .expect("solve")
        .sampled_points
        .expect("sampling enabled")
}

// (a) absent flag == explicit --drops-reference los, byte-identical.
#[test]
fn explicit_los_flag_is_byte_identical() {
    for extra in [
        vec![],
        vec!["-o", "csv", "--full"],
        vec!["-o", "json", "--full"],
    ] {
        let mut with_flag: Vec<&str> = BASE_ARGS.to_vec();
        with_flag.extend(["--drops-reference", "los"]);
        with_flag.extend(extra.iter());
        let mut without: Vec<&str> = BASE_ARGS.to_vec();
        without.extend(extra.iter());
        assert_eq!(
            run_ok(&with_flag).stdout,
            run_ok(&without).stdout,
            "explicit --drops-reference los must not change a single output byte ({extra:?})"
        );
    }
}

// (b) target mode at a 30 degree incline: every sampled drop is exactly drop/cos(30);
// distance/drift/velocity/energy/time are bit-identical.
#[test]
fn target_mode_scales_sampled_drop_by_inverse_cos_exactly() {
    let angle = 30.0_f64.to_radians();
    let los = sampled(angle, DropsReference::Los);
    let target = sampled(angle, DropsReference::Target);
    assert_eq!(los.len(), target.len());
    assert!(los.len() > 3, "expected several sampled rows");

    let cos_theta = angle.cos();
    for (l, t) in los.iter().zip(target.iter()) {
        assert_eq!(
            t.drop_m.to_bits(),
            (l.drop_m / cos_theta).to_bits(),
            "target drop must be exactly los drop / cos(30) at {} m",
            l.distance_m
        );
        assert_eq!(t.distance_m.to_bits(), l.distance_m.to_bits());
        assert_eq!(t.wind_drift_m.to_bits(), l.wind_drift_m.to_bits());
        assert_eq!(t.velocity_mps.to_bits(), l.velocity_mps.to_bits());
        assert_eq!(t.energy_j.to_bits(), l.energy_j.to_bits());
        assert_eq!(t.time_s.to_bits(), l.time_s.to_bits());
    }
}

// (c) target mode at 0 degrees: cos(0) == 1.0, division is bit-exact — values identical.
#[test]
fn target_mode_flat_values_match_los() {
    let los = sampled(0.0, DropsReference::Los);
    let target = sampled(0.0, DropsReference::Target);
    assert_eq!(los.len(), target.len());
    for (l, t) in los.iter().zip(target.iter()) {
        assert_eq!(t.drop_m.to_bits(), l.drop_m.to_bits());
    }
}

// (d) table and CSV labels carry the active reference.
#[test]
fn target_mode_relabels_table_and_csv() {
    let mut table_args: Vec<&str> = BASE_ARGS.to_vec();
    table_args.extend(["--shooting-angle", "30", "--drops-reference", "target"]);
    let table = String::from_utf8(run_ok(&table_args).stdout).unwrap();
    assert!(
        table.contains("Drop (target)(in)"),
        "target-mode sampled table must be labeled 'Drop (target)(in)':\n{table}"
    );

    let mut csv_args: Vec<&str> = table_args.clone();
    csv_args.extend(["-o", "csv", "--full"]);
    let csv = String::from_utf8(run_ok(&csv_args).stdout).unwrap();
    assert!(
        csv.lines()
            .next()
            .is_some_and(|h| h.contains(",drop_target_in,")),
        "target-mode sampled CSV header must carry drop_target_in:\n{csv}"
    );

    // LOS mode keeps the historical labels.
    let mut los_csv_args: Vec<&str> = BASE_ARGS.to_vec();
    los_csv_args.extend(["--shooting-angle", "30", "-o", "csv", "--full"]);
    let los_csv = String::from_utf8(run_ok(&los_csv_args).stdout).unwrap();
    assert!(
        los_csv.lines().next().is_some_and(|h| h.contains(",drop_in,")),
        "LOS-mode sampled CSV header must stay drop_in:\n{los_csv}"
    );
}

// (e) solve-json wire path.
#[test]
fn solve_json_wire_field_decodes_scales_and_defaults() {
    let request_json = |shot_extra: &str| -> String {
        format!(
            r#"{{
                "schema_version": 1,
                "projectile": {{"mass_kg": 0.01089, "diameter_m": 0.007823,
                                "drag_model": "G1", "ballistic_coefficient": 0.475}},
                "rifle": {{"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05}},
                "shot": {{"max_range_m": 600.0, "shooting_angle_rad": 0.5235987755982988{shot_extra}}},
                "atmosphere": {{}},
                "wind": {{}},
                "solver": {{}},
                "effects": {{}},
                "sampling": {{"interval_m": 100.0}}
            }}"#
        )
    };

    let solve = |json: &str| -> ballistics_engine::solve_json::SolveSuccessV1 {
        let request = ballistics_engine::solve_json::decode_solve_request_v1(json)
            .expect("request with drops_reference must decode (not unknown_field)");
        ballistics_engine::solve_v1(request).expect("solve")
    };

    let baseline = solve(&request_json(""));
    let target = solve(&request_json(r#", "drops_reference": "target""#));
    let cos_theta = 0.5235987755982988_f64.cos();

    assert_eq!(baseline.samples.len(), target.samples.len());
    for (b, t) in baseline.samples.iter().zip(target.samples.iter()) {
        // serde_json round-trips f64 exactly (ryu), so the division is still bit-exact
        // on the wire.
        assert_eq!(
            t.drop_m.to_bits(),
            (b.drop_m / cos_theta).to_bits(),
            "wire drop_m must be exactly los drop_m / cos(30) at {} m",
            b.distance_m
        );
        assert_eq!(t.windage_m.to_bits(), b.windage_m.to_bits());
        assert_eq!(t.time_s.to_bits(), b.time_s.to_bits());
    }
    assert_eq!(
        serde_json::to_string(&baseline.summary).unwrap(),
        serde_json::to_string(&target.summary).unwrap(),
        "the summary block is not a drop output and must not change"
    );

    // Omitted == explicit "los", byte-identical response.
    let explicit_los = solve(&request_json(r#", "drops_reference": "los""#));
    assert_eq!(
        serde_json::to_string(&baseline).unwrap(),
        serde_json::to_string(&explicit_los).unwrap(),
        "explicit \"los\" must be byte-identical to omitting the field"
    );

    // Invalid values are structured errors, not silent defaults.
    let error = ballistics_engine::solve_json::decode_solve_request_v1(&request_json(
        r#", "drops_reference": "bogus""#,
    ))
    .expect_err("an invalid drops_reference must be rejected");
    let rendered = serde_json::to_string(&error).unwrap();
    assert!(
        rendered.contains("drops_reference"),
        "error must name the offending field: {rendered}"
    );
}

// (f) target mode at or beyond 90 degrees is rejected (1/cos undefined).
#[test]
fn target_mode_rejects_vertical_shooting_angles() {
    let mut solver = TrajectorySolver::new(
        library_inputs(90.0_f64.to_radians(), DropsReference::Target),
        WindConditions::default(),
        AtmosphericConditions::default(),
    );
    solver.set_max_range(600.0);
    let error = solver
        .solve()
        .expect_err("target drops reference at 90 degrees must be rejected");
    assert!(
        error.to_string().contains("drops reference"),
        "unexpected error: {error}"
    );

    // The same angle in LOS mode still solves (the guard is target-mode only).
    let mut los_solver = TrajectorySolver::new(
        library_inputs(90.0_f64.to_radians(), DropsReference::Los),
        WindConditions::default(),
        AtmosphericConditions::default(),
    );
    los_solver.set_max_range(600.0);
    los_solver.solve().expect("LOS mode is unaffected");
}
