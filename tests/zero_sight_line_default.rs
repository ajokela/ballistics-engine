//! The solve-json v1 zero frame (MBA-1537).
//!
//! `shot.target_height_m` is a world-vertical height above the local ground datum. It used to
//! default to a flat `0` even when the elevation search ran, which aimed the search at the
//! GROUND rather than at the line of sight: a caller who said "zero at 100 yards" and nothing
//! else got a solve sitting a full line-of-sight height LOW at the stated zero, never crossing
//! the line of sight at all. Consequences pinned here: `rifle.sight_height_m` had no effect
//! whatever on the solved elevation, a raised `rifle.muzzle_height_m` moved it (it should not:
//! heights above the ground cancel) and then stopped it converging at all once the muzzle
//! cleared the bullet's drop over the zero distance, and solve-json disagreed with the zero
//! surfaces that solve `ZeroTargetFrame::SightLine` against the sight height --
//! `calculate_zero_angle_with_conditions` and its resolved-conditions twin, and the C ABI's
//! `ballistics_calculate_zero_angle`.
//!
//! The default is now the line of sight, `muzzle_height_m + sight_height_m`, applied exactly
//! when the elevation search runs. An explicitly supplied height still wins, `0.0` included, so
//! a caller who means the ground datum can still say so.

use ballistics_engine::solve_json::{decode_solve_request_v1, SolveSuccessV1};
use ballistics_engine::{
    calculate_zero_angle_with_resolved_conditions, solve_v1, AtmosphericConditions,
    BallisticInputs, DragModel, WindConditions,
};

const ZERO_DISTANCE_M: f64 = 91.44; // 100 yd
const SIGHT_HEIGHT_M: f64 = 0.0381; // 1.5 in
const MUZZLE_VELOCITY_MPS: f64 = 807.72;
const MASS_KG: f64 = 0.01134;
const DIAMETER_M: f64 = 0.00782;
const BC: f64 = 0.243;

/// `find_zero_angle` converges once the trial height is within 1e-4 m of the target, and the
/// reported sample is linearly interpolated onto the zero distance on top of that. 1 mm is a
/// comfortable ceiling over both and still two orders of magnitude under the sight height this
/// test is distinguishing from zero.
const CROSSES_TOLERANCE_M: f64 = 0.001;

const HISTORICAL_GROUND_DATUM_NOTICE: &str =
    "Target height defaulted to 0 m above the local ground datum.";

/// A .308 175 gr G7 in calm sea-level air. `sampling.interval_m` is the zero distance itself so
/// a sample lands exactly on it; `length_m` is omitted so the projectile-length estimate is the
/// same one the C ABI derives (it carries no length field), which the three-surface test below
/// depends on.
fn request_json(rifle_extra: &str, shot_extra: &str) -> String {
    format!(
        r#"{{
            "schema_version": 1,
            "projectile": {{"mass_kg": {MASS_KG}, "diameter_m": {DIAMETER_M},
                            "drag_model": "G7", "ballistic_coefficient": {BC}}},
            "rifle": {{"muzzle_velocity_mps": {MUZZLE_VELOCITY_MPS}{rifle_extra}}},
            "shot": {{"max_range_m": 300.0{shot_extra}}},
            "atmosphere": {{"altitude_m": 0.0, "temperature_k": 288.15, "pressure_pa": 101325.0,
                            "relative_humidity": 0.5}},
            "wind": {{"speed_mps": 0.0, "direction_from_rad": 0.0}},
            "solver": {{"method": "rk4", "time_step_s": 0.001}},
            "effects": {{}},
            "sampling": {{"interval_m": {ZERO_DISTANCE_M}}}
        }}"#
    )
}

fn solve(json: &str) -> SolveSuccessV1 {
    let request = decode_solve_request_v1(json).expect("request must decode");
    solve_v1(request).expect("request must solve")
}

fn sight(height_m: f64) -> String {
    format!(r#", "sight_height_m": {height_m}"#)
}

fn zero_only() -> String {
    format!(r#", "zero_distance_m": {ZERO_DISTANCE_M}"#)
}

/// `drop_m` at the zero distance, positive BELOW the line of sight.
fn drop_at_zero(response: &SolveSuccessV1) -> f64 {
    response
        .samples
        .iter()
        .find(|sample| (sample.distance_m - ZERO_DISTANCE_M).abs() < 0.05)
        .unwrap_or_else(|| panic!("no sample at the zero distance {ZERO_DISTANCE_M} m"))
        .drop_m
}

fn target_height_notices(response: &SolveSuccessV1) -> Vec<&str> {
    response
        .assumptions
        .iter()
        .filter(|notice| notice.path.as_deref() == Some("$.shot.target_height_m"))
        .map(|notice| notice.message.as_str())
        .collect()
}

/// The headline: a request that states only a zero distance produces a trajectory that actually
/// crosses the line of sight there. Checked over a spread of sight heights because the previous
/// behaviour left a residual of exactly the sight height, so a single height could not tell
/// "crosses" from "off by a constant".
#[test]
fn a_zero_distance_solve_with_no_target_height_crosses_the_line_of_sight() {
    for sight_height_m in [0.0, 0.0381, 0.0508, 0.06604, 0.1] {
        let response = solve(&request_json(&sight(sight_height_m), &zero_only()));

        assert_eq!(
            response.resolved_request.shot.target_height_m, sight_height_m,
            "with muzzle_height_m at its 0 default the resolved target height is the sight height"
        );
        let drop_m = drop_at_zero(&response);
        assert!(
            drop_m.abs() < CROSSES_TOLERANCE_M,
            "sight height {sight_height_m} m: the trajectory is {drop_m} m off the line of \
             sight at its own stated zero distance"
        );
    }
}

/// The sharpest statement of the old defect: the solved elevation was bit-identical for every
/// sight height, so `rifle.sight_height_m` did not reach the zero at all. It does now, and by
/// the geometric amount -- the extra height the bullet must be carrying at the zero distance,
/// over the distance it has to be carried.
#[test]
fn sight_height_moves_the_solved_elevation_by_the_geometry() {
    let low = 0.0381;
    let high = 0.06604;

    let low_angle = solve(&request_json(&sight(low), &zero_only()))
        .resolved_request
        .shot
        .muzzle_angle_rad;
    let high_angle = solve(&request_json(&sight(high), &zero_only()))
        .resolved_request
        .shot
        .muzzle_angle_rad;

    let measured = high_angle - low_angle;
    let geometric = (high - low) / ZERO_DISTANCE_M;
    assert!(
        measured > 0.0,
        "a taller sight must need more elevation, not less (low {low_angle}, high {high_angle})"
    );
    // The two angles come out of a bisection with a 1e-7 rad angle floor, and the trajectory is
    // not perfectly linear over the zero distance, so this is a proportionality check rather
    // than an identity.
    assert!(
        (measured - geometric).abs() < 0.05 * geometric,
        "the elevation difference {measured} rad should track the sight-height difference over \
         the zero distance ({geometric} rad)"
    );
}

/// An explicitly supplied height still decides the frame, including the `0.0` that used to be
/// the silent default: a caller who means the ground datum can still ask for it, and gets the
/// old numbers -- a trajectory a full sight height below the line of sight at the zero.
#[test]
fn an_explicit_zero_target_height_still_zeroes_to_the_ground_datum() {
    let ground = solve(&request_json(
        &sight(SIGHT_HEIGHT_M),
        &format!(r#"{}, "target_height_m": 0.0"#, zero_only()),
    ));

    assert_eq!(ground.resolved_request.shot.target_height_m, 0.0);
    assert!(
        target_height_notices(&ground).is_empty(),
        "a supplied value is not a default and must not raise a default notice"
    );

    let drop_m = drop_at_zero(&ground);
    assert!(
        (drop_m - SIGHT_HEIGHT_M).abs() < CROSSES_TOLERANCE_M,
        "a ground-datum zero puts the bullet a sight height ({SIGHT_HEIGHT_M} m) below the line \
         of sight at the zero distance; got {drop_m} m"
    );

    // And it is the ground datum specifically, not "whatever the sight height happens to be":
    // the solved elevation is the same for a different sight height, which is exactly the
    // behaviour that used to apply to every request.
    let taller = solve(&request_json(
        &sight(0.06604),
        &format!(r#"{}, "target_height_m": 0.0"#, zero_only()),
    ));
    assert_eq!(
        ground.resolved_request.shot.muzzle_angle_rad, taller.resolved_request.shot.muzzle_angle_rad,
        "an explicit ground-datum zero is indifferent to sight height, as it always was"
    );
}

/// `muzzle_height_m` is part of the line of sight, not a separate correction. A raised muzzle
/// used to make this request fail to converge outright, because the ground datum it was aiming
/// at was then below the muzzle.
#[test]
fn a_raised_muzzle_composes_with_the_sight_height() {
    let mut angles = Vec::new();

    for muzzle_height_m in [0.0, 0.5, 1.2] {
        let response = solve(&request_json(
            &format!(
                r#"{}, "muzzle_height_m": {muzzle_height_m}"#,
                sight(SIGHT_HEIGHT_M)
            ),
            &zero_only(),
        ));

        assert_eq!(
            response.resolved_request.shot.target_height_m,
            muzzle_height_m + SIGHT_HEIGHT_M,
            "the resolved target height is the line of sight above the ground datum"
        );
        let drop_m = drop_at_zero(&response);
        assert!(
            drop_m.abs() < CROSSES_TOLERANCE_M,
            "muzzle height {muzzle_height_m} m: {drop_m} m off the line of sight at the zero"
        );
        angles.push(response.resolved_request.shot.muzzle_angle_rad);
    }

    // Heights above the ground cancel once the target tracks the line of sight, so the solved
    // elevation must not depend on how high the muzzle is held -- the same invariance
    // `tests/zero_sight_height.rs` pins for the native zero, which the old flat-`0` default
    // broke on this surface (it made the angle a function of muzzle height, then stopped
    // converging altogether).
    assert!(
        angles.windows(2).all(|pair| pair[0] == pair[1]),
        "the zero angle must be invariant to muzzle height; got {angles:?}"
    );
}

/// The default is a LEVEL line of sight in a world-vertical field, so it does not describe an
/// inclined shot's sight line — and `docs/SOLVE_JSON_V1.md` hands inclined callers the
/// projection to supply instead. This pins that formula, and pins that the omitted-field case
/// really is the one the doc warns about rather than quietly working.
#[test]
fn an_inclined_zero_needs_its_target_height_stated() {
    for degrees in [-10.0_f64, 10.0] {
        let radians = degrees.to_radians();
        let projected =
            ZERO_DISTANCE_M * radians.sin() + (0.0 + SIGHT_HEIGHT_M) * radians.cos();

        let stated = solve(&request_json(
            &sight(SIGHT_HEIGHT_M),
            &format!(
                r#"{}, "shooting_angle_rad": {radians}, "target_height_m": {projected}"#,
                zero_only()
            ),
        ));
        let drop_m = drop_at_zero(&stated);
        assert!(
            drop_m.abs() < CROSSES_TOLERANCE_M,
            "{degrees} deg: the documented projection should cross the line of sight at the \
             zero distance; got {drop_m} m"
        );

        // And the omitted case does NOT, which is why the doc tells a caller to state it.
        let omitted = decode_solve_request_v1(&request_json(
            &sight(SIGHT_HEIGHT_M),
            &format!(r#"{}, "shooting_angle_rad": {radians}"#, zero_only()),
        ))
        .expect("request must decode");
        match solve_v1(omitted) {
            // Uphill the level default is unreachable and the search says so.
            Err(_) => assert!(degrees > 0.0, "{degrees} deg: unexpected solve failure"),
            // Downhill it converges, on the world height the default names rather than on the
            // sight line -- a long way off at the stated zero.
            Ok(response) => assert!(
                drop_at_zero(&response).abs() > 1.0,
                "{degrees} deg: an omitted inclined target height is documented as NOT zeroing \
                 to the sight line; if that has changed, the doc and this test must change with it"
            ),
        }
    }
}

/// A supplied `muzzle_angle_rad` runs no elevation search, so there is no zero to frame and
/// nothing about such a request changes: the angle is used verbatim, the resolved target height
/// stays `0`, and the notice is still the historical ground-datum wording. Pinned for the angle
/// alone and for the angle alongside a zero distance, which is what rebuilding a request from a
/// previous `resolved_request` produces.
#[test]
fn a_supplied_muzzle_angle_is_untouched_by_the_new_default() {
    let angle_rad = 0.00123456;

    for shot_extra in [
        format!(r#", "muzzle_angle_rad": {angle_rad}"#),
        format!(r#"{}, "muzzle_angle_rad": {angle_rad}"#, zero_only()),
    ] {
        let response = solve(&request_json(&sight(SIGHT_HEIGHT_M), &shot_extra));

        assert_eq!(
            response.resolved_request.shot.muzzle_angle_rad, angle_rad,
            "{shot_extra}: the supplied angle must be used verbatim"
        );
        assert_eq!(
            response.resolved_request.shot.target_height_m, 0.0,
            "{shot_extra}: no elevation search ran, so the zero-frame default must not apply"
        );
        assert_eq!(
            target_height_notices(&response),
            vec![HISTORICAL_GROUND_DATUM_NOTICE],
            "{shot_extra}: the historical notice text is part of the unchanged response"
        );
    }
}

/// The notice announces the default exactly when it is applied -- a defaulted value that moves
/// every elevation number in the response has to say so -- and stays silent when the caller
/// supplied the height themselves.
#[test]
fn the_line_of_sight_notice_appears_exactly_when_the_default_is_applied() {
    let explicit_height = format!(r#", "target_height_m": {SIGHT_HEIGHT_M}"#);
    let explicit_angle = r#", "muzzle_angle_rad": 0.001"#;

    // (shot fields, does the line-of-sight default apply?)
    let cases: [(String, bool); 6] = [
        (zero_only(), true),
        (format!("{}{explicit_height}", zero_only()), false),
        (format!("{}{explicit_angle}", zero_only()), false),
        (String::new(), false),
        (explicit_height.clone(), false),
        (explicit_angle.to_string(), false),
    ];

    for (shot_extra, expect_default) in cases {
        let response = solve(&request_json(&sight(SIGHT_HEIGHT_M), &shot_extra));
        let notices = target_height_notices(&response);

        if expect_default {
            assert_eq!(notices.len(), 1, "shot {shot_extra:?}: expected one notice");
            let message = notices[0];
            assert!(
                message.contains("line of sight"),
                "shot {shot_extra:?}: the notice must name the frame, got {message:?}"
            );
            assert!(
                message.contains(&SIGHT_HEIGHT_M.to_string()),
                "shot {shot_extra:?}: the notice must name the value it applied, got {message:?}"
            );
        } else if shot_extra.contains("target_height_m") {
            assert!(
                notices.is_empty(),
                "shot {shot_extra:?}: a supplied height is not a default, got {notices:?}"
            );
        } else {
            assert_eq!(
                notices,
                vec![HISTORICAL_GROUND_DATUM_NOTICE],
                "shot {shot_extra:?}: with no elevation search the historical default stands"
            );
        }
    }
}

/// The point of the ticket: a caller moving between the JSON bridge and the native zero helper
/// no longer silently changes rifles. Both are handed the same rifle, the same sea-level station
/// reading (where neither surface's atmosphere resolution alters it), and the same zero
/// distance; solve-json is told nothing about the target height and the native helper is told
/// the sight height, which is what its `ZeroTargetFrame::SightLine` contract means.
#[test]
fn solve_json_and_the_native_zero_helper_agree() {
    let bridge = solve(&request_json(&sight(SIGHT_HEIGHT_M), &zero_only()))
        .resolved_request
        .shot
        .muzzle_angle_rad;

    let native = calculate_zero_angle_with_resolved_conditions(
        native_inputs(),
        ZERO_DISTANCE_M,
        SIGHT_HEIGHT_M,
        WindConditions::default(),
        AtmosphericConditions {
            temperature: 15.0,
            pressure: 1013.25,
            humidity: 50.0,
            altitude: 0.0,
        },
    )
    .expect("the native zero solves");

    assert_eq!(
        bridge, native,
        "solve-json's omitted-target-height zero and the native sight-line zero are the same \
         solve and must return the same angle"
    );
}

/// The same rifle the JSON fixture describes, as engine inputs. `bullet_length` is left at the
/// mass-and-diameter estimate, which is what both the JSON path (no `length_m` supplied) and the
/// C ABI (no length field at all) derive.
fn native_inputs() -> BallisticInputs {
    BallisticInputs {
        muzzle_velocity: MUZZLE_VELOCITY_MPS,
        bc_value: BC,
        bc_type: DragModel::G7,
        bullet_mass: MASS_KG,
        bullet_diameter: DIAMETER_M,
        sight_height: SIGHT_HEIGHT_M,
        use_rk4: true,
        use_adaptive_rk45: false,
        ..BallisticInputs::default()
    }
}

/// The third surface. The C ABI carries no target height at all: it passes the sight height for
/// the caller (`src/ffi.rs`, "we want the bullet to hit at sight height at the zero distance").
/// At sea level the legacy station-atmosphere sentinels do not fire, so the two surfaces are
/// being asked the same physical question and their answers are directly comparable.
#[cfg(feature = "ffi")]
#[test]
fn the_c_abi_zero_agrees_with_solve_json() {
    use ballistics_engine::ffi::{
        ballistics_calculate_zero_angle, FFIAtmosphericConditions, FFIBallisticInputs,
        FFIWindConditions,
    };

    let bridge = solve(&request_json(&sight(SIGHT_HEIGHT_M), &zero_only()))
        .resolved_request
        .shot
        .muzzle_angle_rad;

    let inputs = FFIBallisticInputs {
        muzzle_velocity: MUZZLE_VELOCITY_MPS,
        muzzle_angle: 0.0,
        bc_value: BC,
        bullet_mass: MASS_KG,
        bullet_diameter: DIAMETER_M,
        bc_type: 1, // G7
        sight_height: SIGHT_HEIGHT_M,
        target_distance: ZERO_DISTANCE_M,
        temperature: 15.0,
        twist_rate: BallisticInputs::default().twist_rate,
        is_twist_right: 1,
        shooting_angle: 0.0,
        altitude: 0.0,
        latitude: f64::NAN,
        azimuth_angle: 0.0,
        use_rk4: 1,
        use_adaptive_rk45: 0,
        enable_wind_shear: 0,
        enable_trajectory_sampling: 0,
        sample_interval: 0.0,
        enable_pitch_damping: 0,
        enable_precession_nutation: 0,
        enable_spin_drift: 0,
        enable_magnus: 0,
        enable_coriolis: 0,
        shot_azimuth: 0.0,
        cant_angle: 0.0,
        zero_poi_vertical: 0.0,
        zero_poi_horizontal: 0.0,
        sight_offset_lateral: 0.0,
    };
    let wind = FFIWindConditions {
        speed: 0.0,
        direction: 0.0,
        vertical_speed: 0.0,
    };
    let atmosphere = FFIAtmosphericConditions {
        temperature: 15.0,
        pressure: 1013.25,
        humidity: 50.0,
        altitude: 0.0,
    };

    let c_abi =
        unsafe { ballistics_calculate_zero_angle(&inputs, &wind, &atmosphere, ZERO_DISTANCE_M) };

    assert!(c_abi.is_finite(), "the C ABI zero returned {c_abi}");
    assert_eq!(
        c_abi, bridge,
        "the C ABI and solve-json must resolve the same zero for the same rifle"
    );
}
