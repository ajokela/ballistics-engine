//! MBA-1395: equivalent horizontal range (BDC shoot-to) for inclined shots.
//!
//! Definition pinned here: the flat-fire range whose ANGULAR elevation correction —
//! against the SAME dialed zero — matches the inclined solution's at the target range
//! (SIG AMR / Leica EHR / Gunwerks style), computed by inverse lookup over one flat
//! re-solve. Explicitly NOT the rifleman's-rule `range * cos(angle)` linear-drop
//! approximation. Pins:
//!   (a) a 30 degree inclined shot's EHR exists and is shorter than the true range;
//!   (b) EHR is monotone in target range;
//!   (c) ill-defined regions return None: target at/inside the zero range, and a
//!       negative correction (bullet above the LOS past a near zero);
//!   (d) the angular-match definition holds (flat correction at the EHR equals the
//!       inclined correction at the target) AND diverges measurably from rifleman's
//!       rule at long range — we must not match the cosine rule;
//!   (e) CLI: flat shots and un-zeroed shots print no EHR line (byte-identical class),
//!       inclined zeroed shots print exactly one;
//!   (f) solve-json v1: summary.equivalent_horizontal_range_m present only for the
//!       inclined zeroed request, absent (not null) otherwise.

use ballistics_engine::{
    calculate_zero_angle, AtmosphericConditions, BallisticInputs, TrajectorySolver,
    WindConditions,
};
use std::process::Command;

const BIN: &str = env!("CARGO_BIN_EXE_ballistics");
const ZERO_100YD_M: f64 = 91.44;

fn base_inputs() -> BallisticInputs {
    BallisticInputs {
        bc_value: 0.475,
        muzzle_velocity: 823.0,
        bullet_mass: 0.01089,
        bullet_diameter: 0.007823,
        sight_height: 0.05,
        ..Default::default()
    }
}

/// Inputs dialed to a solved zero at `zero_distance_m`, then inclined by
/// `shooting_angle_rad` — the same "zero level, shoot inclined" state the CLI builds.
fn zeroed_inclined_inputs(zero_distance_m: f64, shooting_angle_rad: f64) -> BallisticInputs {
    let mut inputs = base_inputs();
    let los_height = inputs.muzzle_height + inputs.sight_height;
    let zero_angle = calculate_zero_angle(inputs.clone(), zero_distance_m, los_height)
        .expect("zero solve");
    inputs.muzzle_angle = zero_angle;
    inputs.shooting_angle = shooting_angle_rad;
    inputs
}

fn solver_for(inputs: BallisticInputs, max_range_m: f64) -> TrajectorySolver {
    let mut solver = TrajectorySolver::new(
        inputs,
        WindConditions::default(),
        AtmosphericConditions::default(),
    );
    solver.set_max_range(max_range_m);
    solver
}

// (a) 30 degree incline: EHR exists, lies past the zero, and is SHORTER than true range.
#[test]
fn inclined_ehr_is_shorter_than_true_range() {
    let angle = 30.0_f64.to_radians();
    let solver = solver_for(zeroed_inclined_inputs(ZERO_100YD_M, angle), 800.0);
    let target = 731.52; // 800 yd
    let ehr = solver
        .equivalent_horizontal_range(target, ZERO_100YD_M)
        .expect("EHR must exist for a normal inclined shot past the zero");
    assert!(
        ehr < target,
        "shoot-to must be shorter than the true range: {ehr} vs {target}"
    );
    assert!(
        ehr > ZERO_100YD_M,
        "shoot-to must lie past the zero range: {ehr}"
    );
}

// (b) EHR is monotone in target range.
#[test]
fn ehr_is_monotone_in_target_range() {
    let angle = 30.0_f64.to_radians();
    let solver = solver_for(zeroed_inclined_inputs(ZERO_100YD_M, angle), 800.0);
    let mut previous = 0.0;
    for target in [274.32, 457.2, 640.08] {
        let ehr = solver
            .equivalent_horizontal_range(target, ZERO_100YD_M)
            .unwrap_or_else(|| panic!("EHR must exist at {target} m"));
        assert!(
            ehr > previous,
            "EHR must increase with target range: {ehr} after {previous}"
        );
        previous = ehr;
    }
}

// (c) ill-defined regions return None instead of a guessed number.
#[test]
fn boundary_and_negative_correction_return_none() {
    let angle = 30.0_f64.to_radians();
    let solver = solver_for(zeroed_inclined_inputs(ZERO_100YD_M, angle), 800.0);

    // Target exactly at, and inside, the zero range.
    assert_eq!(solver.equivalent_horizontal_range(ZERO_100YD_M, ZERO_100YD_M), None);
    assert_eq!(solver.equivalent_horizontal_range(45.0, ZERO_100YD_M), None);

    // Negative correction: zeroed at 25 yd (a NEAR zero), the bullet rides ABOVE the
    // LOS at 100 yd — past the zero range, but no positive BDC correction exists.
    let near_zero_m = 22.86;
    let near_solver = solver_for(zeroed_inclined_inputs(near_zero_m, angle), 800.0);
    assert_eq!(
        near_solver.equivalent_horizontal_range(91.44, near_zero_m),
        None,
        "a bullet above the LOS has no defined shoot-to range"
    );
}

// (d) the angular-match definition is pinned — and it is NOT rifleman's rule.
#[test]
fn angular_match_definition_diverges_from_riflemans_rule() {
    let angle = 45.0_f64.to_radians();
    let target = 731.52; // 800 yd — deep enough for drag to separate the two rules
    let inputs = zeroed_inclined_inputs(ZERO_100YD_M, angle);
    let solver = solver_for(inputs.clone(), 800.0);
    let ehr = solver
        .equivalent_horizontal_range(target, ZERO_100YD_M)
        .expect("EHR");

    // Defining property: the flat solution's angular correction at the EHR equals the
    // inclined solution's at the true range (both against the same zero).
    let inclined_result = solver_for(inputs.clone(), 800.0).solve().expect("inclined");
    let inclined_obs = inclined_result
        .observation_at_range_checked(target)
        .expect("inclined observation");
    let inclined_correction = inclined_obs.drop_m / target;
    assert!(
        inclined_correction > 0.0,
        "test setup: the inclined correction must be positive"
    );

    let mut flat_inputs = inputs;
    flat_inputs.shooting_angle = 0.0;
    let flat_result = solver_for(flat_inputs, 800.0).solve().expect("flat");
    let flat_obs = flat_result
        .observation_at_range_checked(ehr)
        .expect("flat observation at the EHR");
    let flat_correction = flat_obs.drop_m / ehr;
    assert!(
        (flat_correction - inclined_correction).abs() < 1e-6,
        "angular-match definition violated: flat correction {flat_correction} at EHR \
         {ehr} vs inclined correction {inclined_correction} at {target}"
    );

    // And it must NOT be the rifleman's-rule cosine range: at 800 yd / 45 degrees the
    // two disagree by well over the solver's centimeter-level convergence.
    let riflemans = target * angle.cos();
    assert!(
        (ehr - riflemans).abs() > 5.0,
        "EHR {ehr} must not reproduce rifleman's rule {riflemans} — the ticket pins \
         the angular-match definition, not linear drop"
    );
}

fn run_ok(args: &[&str]) -> String {
    let output = Command::new(BIN).args(args).output().expect("run binary");
    assert!(
        output.status.success(),
        "command failed: {:?}\nstderr: {}",
        args,
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).expect("utf8 stdout")
}

const CLI_BASE: &[&str] = &[
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
];

// (e) CLI emission rules.
#[test]
fn cli_emits_line_only_for_inclined_zeroed_shots() {
    // Flat + zeroed: nothing.
    let mut flat: Vec<&str> = CLI_BASE.to_vec();
    flat.extend(["--auto-zero", "100"]);
    assert!(
        !run_ok(&flat).contains("Equivalent horizontal range"),
        "a flat shot must not print an EHR line"
    );

    // Inclined + un-zeroed (bare bore angle): nothing — no zero, no BDC semantics.
    let mut unzeroed: Vec<&str> = CLI_BASE.to_vec();
    unzeroed.extend(["--shooting-angle", "30"]);
    assert!(
        !run_ok(&unzeroed).contains("Equivalent horizontal range"),
        "an un-zeroed shot must not print an EHR line"
    );

    // Inclined + zeroed: exactly one line, value below the solved range, past the zero.
    let mut inclined: Vec<&str> = CLI_BASE.to_vec();
    inclined.extend(["--shooting-angle", "30", "--auto-zero", "100"]);
    let stdout = run_ok(&inclined);
    let lines: Vec<&str> = stdout
        .lines()
        .filter(|l| l.contains("Equivalent horizontal range"))
        .collect();
    assert_eq!(lines.len(), 1, "exactly one EHR summary line:\n{stdout}");
    let line = lines[0];
    assert!(
        line.contains("(shoot-to for BDC)") && line.contains("yd"),
        "unexpected line format: {line}"
    );
    let value: f64 = line
        .trim_start_matches("Equivalent horizontal range:")
        .split_whitespace()
        .next()
        .expect("value")
        .parse()
        .expect("numeric shoot-to range");
    assert!(
        value > 100.0 && value < 600.0,
        "shoot-to must lie between the zero and the true range: {value}"
    );
}

// (f) solve-json v1: field present only for the inclined zeroed request.
#[test]
fn solve_json_summary_field_only_when_defined() {
    let request_json = |shot: &str| -> String {
        format!(
            r#"{{
                "schema_version": 1,
                "projectile": {{"mass_kg": 0.01089, "diameter_m": 0.007823,
                                "drag_model": "G1", "ballistic_coefficient": 0.475}},
                "rifle": {{"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05}},
                "shot": {shot},
                "atmosphere": {{}},
                "wind": {{}},
                "solver": {{}},
                "effects": {{}},
                "sampling": {{"interval_m": 100.0}}
            }}"#
        )
    };
    let solve = |json: &str| -> ballistics_engine::solve_json::SolveSuccessV1 {
        let request =
            ballistics_engine::solve_json::decode_solve_request_v1(json).expect("decode");
        ballistics_engine::solve_v1(request).expect("solve")
    };

    let flat = solve(&request_json(
        r#"{"max_range_m": 600.0, "zero_distance_m": 91.44}"#,
    ));
    assert_eq!(flat.summary.equivalent_horizontal_range_m, None);
    assert!(
        !serde_json::to_string(&flat).unwrap().contains("equivalent_horizontal_range_m"),
        "absent — not null — for a flat solve (byte-identity class)"
    );

    let unzeroed = solve(&request_json(
        r#"{"max_range_m": 600.0, "shooting_angle_rad": 0.5235987755982988}"#,
    ));
    assert_eq!(unzeroed.summary.equivalent_horizontal_range_m, None);

    // NOTE: solve-json v1 zeroes in the WorldVertical frame (MBA-1412), so an inclined
    // zero needs a target height consistent with the slope: the 100 yd LOS point sits at
    // world height 91.44*sin(30) + 0.05*cos(30) for a 0.05 m sight height.
    let inclined = solve(&request_json(
        r#"{"max_range_m": 600.0, "zero_distance_m": 91.44,
            "target_height_m": 45.76330127018922,
            "shooting_angle_rad": 0.5235987755982988}"#,
    ));
    let ehr = inclined
        .summary
        .equivalent_horizontal_range_m
        .expect("inclined zeroed request must carry the shoot-to range");
    assert!(
        ehr > 91.44 && ehr < inclined.summary.actual_range_m,
        "shoot-to must lie between the zero and the true range: {ehr} vs {}",
        inclined.summary.actual_range_m
    );
}
