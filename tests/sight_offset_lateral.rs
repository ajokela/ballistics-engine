//! MBA-1396: lateral sight offset for offset-mounted optics.
//!
//! `--sight-offset` (engine field `sight_offset_lateral_m`, solve-json
//! `rifle.sight_offset_lateral_m`) is PHYSICAL mount geometry: positive = the sight sits
//! right of the bore, so the bullet starts that far LEFT of the sight line
//! (`initial_position` z-term), and a zero solve adds the windage convergence
//! (`offset / zero_distance`) so the trajectory crosses the sight line laterally at the
//! zero range — not a constant parallel offset. Pins:
//!   (a) absent/0.0 flag is byte-identical;
//!   (b) with a zero: drift ~0 at the zero range, left before it, right past it;
//!   (c) without a zero (explicit --angle): the constant physical displacement remains;
//!   (d) composes with cant;
//!   (e) NOT conflated with MBA-1359's zero POI offset — both set => additive, distinct.

use serde_json::Value;
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
    "300",
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

fn full_json(extra: &[&str]) -> Value {
    let mut args: Vec<&str> = BASE_ARGS.to_vec();
    args.extend_from_slice(extra);
    args.extend(["--full", "-o", "json"]);
    serde_json::from_slice(&run_ok(&args).stdout).expect("json output")
}

/// Lateral position (`x`, yards, positive = right) interpolated at downrange `at` yards.
fn lateral_at(doc: &Value, at: f64) -> f64 {
    let points = doc["trajectory"].as_array().expect("trajectory points");
    let mut prev: Option<(f64, f64)> = None;
    for p in points {
        let z = p["z"].as_f64().expect("z");
        let x = p["x"].as_f64().expect("x");
        if let Some((pz, px)) = prev {
            if pz <= at && z >= at && z > pz {
                return px + (x - px) * (at - pz) / (z - pz);
            }
        }
        prev = Some((z, x));
    }
    panic!("no bracketing points around {at} yd");
}

const OFFSET_IN: f64 = 0.5;
const OFFSET_YD: f64 = OFFSET_IN / 36.0;

// (a) explicit 0.0 is byte-identical to a run without the flag.
#[test]
fn zero_valued_flag_is_byte_identical() {
    let mut with_flag: Vec<&str> = BASE_ARGS.to_vec();
    with_flag.extend(["--auto-zero", "100", "--sight-offset", "0", "-o", "json"]);
    let mut without: Vec<&str> = BASE_ARGS.to_vec();
    without.extend(["--auto-zero", "100", "-o", "json"]);
    assert_eq!(
        run_ok(&with_flag).stdout,
        run_ok(&without).stdout,
        "a 0.0 --sight-offset must not change a single output byte"
    );
}

// (b) with a zero: starts left by the offset, crosses the sight line at the zero range,
// and continues right past it (windage-zero convergence, not a parallel offset).
#[test]
fn offset_with_zero_converges_at_zero_range() {
    let baseline = full_json(&["--auto-zero", "100"]);
    let offset = full_json(&["--auto-zero", "100", "--sight-offset", "0.5"]);

    let d_muzzle = lateral_at(&offset, 2.0) - lateral_at(&baseline, 2.0);
    let d_50 = lateral_at(&offset, 50.0) - lateral_at(&baseline, 50.0);
    let d_100 = lateral_at(&offset, 100.0) - lateral_at(&baseline, 100.0);
    let d_200 = lateral_at(&offset, 200.0) - lateral_at(&baseline, 200.0);

    assert!(
        (d_muzzle + OFFSET_YD).abs() < OFFSET_YD * 0.1,
        "near the muzzle the bullet must sit ~{OFFSET_YD} yd LEFT of the sight line, got {d_muzzle}"
    );
    assert!(
        d_50 < -OFFSET_YD * 0.3 && d_50 > -OFFSET_YD,
        "halfway to the zero the displacement must be partially converged, got {d_50}"
    );
    assert!(
        d_100.abs() < OFFSET_YD * 0.05,
        "at the zero range the trajectory must cross the sight line (drift ~0), got {d_100}"
    );
    assert!(
        (d_200 - OFFSET_YD).abs() < OFFSET_YD * 0.1,
        "past the zero the drift must continue RIGHT (~+{OFFSET_YD} yd at 200), got {d_200}"
    );
}

// (c) without a zero solve (explicit --angle), only the physical displacement applies:
// a constant left offset at every range, no convergence.
#[test]
fn offset_without_zero_is_a_constant_displacement() {
    let baseline = full_json(&["--angle", "0.05"]);
    let offset = full_json(&["--angle", "0.05", "--sight-offset", "0.5"]);

    for at in [2.0, 100.0, 250.0] {
        let delta = lateral_at(&offset, at) - lateral_at(&baseline, at);
        assert!(
            (delta + OFFSET_YD).abs() < OFFSET_YD * 0.05,
            "with no zero solve the displacement must stay ~{OFFSET_YD} yd LEFT at {at} yd, got {delta}"
        );
    }
}

// (d) composes with cant: the offset still converges at the zero range on top of the
// cant-induced drift (which the zero deliberately does not correct).
#[test]
fn offset_composes_with_cant() {
    let cant_only = full_json(&["--auto-zero", "100", "--cant", "5"]);
    let cant_and_offset = full_json(&["--auto-zero", "100", "--cant", "5", "--sight-offset", "0.5"]);

    let d_100 = lateral_at(&cant_and_offset, 100.0) - lateral_at(&cant_only, 100.0);
    let d_muzzle = lateral_at(&cant_and_offset, 2.0) - lateral_at(&cant_only, 2.0);
    assert!(
        (d_muzzle + OFFSET_YD).abs() < OFFSET_YD * 0.1,
        "offset must still displace the muzzle position under cant, got {d_muzzle}"
    );
    assert!(
        d_100.abs() < OFFSET_YD * 0.05,
        "offset must still converge at the zero range under cant, got {d_100}"
    );
    // And the cant drift itself must still be present (the two effects are additive):
    // canted fire with an upward zero pushes POI laterally downrange vs an un-canted shot.
    let no_cant = full_json(&["--auto-zero", "100"]);
    let cant_drift_200 = lateral_at(&cant_only, 200.0) - lateral_at(&no_cant, 200.0);
    assert!(
        cant_drift_200.abs() > 1e-4,
        "fixture sanity: cant alone must produce lateral drift at 200 yd"
    );
}

// (e) NOT conflated with MBA-1359: a zero POI right-offset shows AT the zero range, the
// mount offset converges to zero there — both set => additive, distinct effects.
#[test]
fn distinct_from_zero_poi_offset_and_additive() {
    let baseline = full_json(&["--auto-zero", "100"]);
    let poi_only = full_json(&["--auto-zero", "100", "--zero-poi-right", "0.2"]);
    let both = full_json(&[
        "--auto-zero",
        "100",
        "--zero-poi-right",
        "0.2",
        "--sight-offset",
        "0.5",
    ]);

    let poi_yd = 0.2 / 36.0;
    // At the zero range: the mount offset contributes ~0, so `both` shows only the POI bias.
    let both_100 = lateral_at(&both, 100.0) - lateral_at(&baseline, 100.0);
    assert!(
        (both_100 - poi_yd).abs() < poi_yd * 0.15,
        "at the zero range only the POI offset must remain (~{poi_yd} yd right), got {both_100}"
    );
    // At the muzzle: the POI bias contributes ~0, so `both` shows only the mount offset.
    let both_muzzle = lateral_at(&both, 2.0) - lateral_at(&baseline, 2.0);
    assert!(
        (both_muzzle + OFFSET_YD).abs() < OFFSET_YD * 0.1,
        "at the muzzle only the mount displacement must remain (~{OFFSET_YD} yd left), got {both_muzzle}"
    );
    // Additivity: (both - poi_only) reproduces the pure mount-offset signature.
    let mount_part_muzzle = lateral_at(&both, 2.0) - lateral_at(&poi_only, 2.0);
    assert!(
        (mount_part_muzzle + OFFSET_YD).abs() < OFFSET_YD * 0.1,
        "subtracting the POI-only run must leave the mount-offset displacement, got {mount_part_muzzle}"
    );
}

// Wire path (solve-json v1): rifle.sight_offset_lateral_m decodes (NOT unknown_field),
// applies with the convergence, and explicit 0.0 stays byte-identical to omitting it.
#[test]
fn solve_json_wire_field_decodes_and_applies() {
    let request_json = |offset: Option<f64>| -> String {
        let rifle = match offset {
            Some(v) => format!(
                r#"{{"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05,
                    "sight_offset_lateral_m": {v}}}"#
            ),
            None => r#"{"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05}"#.to_string(),
        };
        format!(
            r#"{{
                "schema_version": 1,
                "projectile": {{"mass_kg": 0.01089, "diameter_m": 0.007823,
                                "drag_model": "G1", "ballistic_coefficient": 0.475}},
                "rifle": {rifle},
                "shot": {{"max_range_m": 300.0, "zero_distance_m": 91.44}},
                "atmosphere": {{}},
                "wind": {{}},
                "solver": {{}},
                "effects": {{}},
                "sampling": {{"interval_m": 25.0}}
            }}"#
        )
    };

    let solve = |json: &str| -> ballistics_engine::solve_json::SolveSuccessV1 {
        let request = ballistics_engine::solve_json::decode_solve_request_v1(json)
            .expect("request with sight_offset_lateral_m must decode (not unknown_field)");
        ballistics_engine::solve_v1(request).expect("solve")
    };

    let baseline = solve(&request_json(None));
    let offset = solve(&request_json(Some(0.0127))); // 0.5 in, sight right of bore

    let windage_near = |s: &ballistics_engine::solve_json::SolveSuccessV1, at: f64| {
        s.samples
            .iter()
            .min_by(|a, b| {
                (a.distance_m - at)
                    .abs()
                    .total_cmp(&(b.distance_m - at).abs())
            })
            .map(|s| (s.distance_m, s.windage_m))
            .expect("sample")
    };

    // Early sample: displaced left (negative windage), partially converged.
    let (d_early, w_early) = windage_near(&offset, 25.0);
    let (_, w_early_base) = windage_near(&baseline, 25.0);
    let expected_early = -0.0127 * (1.0 - d_early / 91.44);
    assert!(
        ((w_early - w_early_base) - expected_early).abs() < 0.0015,
        "expected ~{expected_early} m left at {d_early} m, got {}",
        w_early - w_early_base
    );
    // At the zero range: converged onto the sight line.
    let (_, w_zero) = windage_near(&offset, 91.44);
    let (_, w_zero_base) = windage_near(&baseline, 91.44);
    assert!(
        (w_zero - w_zero_base).abs() < 0.0015,
        "expected ~0 windage delta at the zero range, got {}",
        w_zero - w_zero_base
    );

    // Explicit 0.0 == omitted, byte-identical response.
    let explicit_zero = solve(&request_json(Some(0.0)));
    assert_eq!(
        serde_json::to_string(&baseline).unwrap(),
        serde_json::to_string(&explicit_zero).unwrap(),
        "explicit 0.0 sight_offset_lateral_m must be byte-identical to omitting it"
    );
}
