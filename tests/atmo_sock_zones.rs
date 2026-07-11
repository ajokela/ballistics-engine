//! MBA-1137: per-zone (downrange-varying) atmosphere affecting drag, composing with the
//! MBA-1136 altitude-lapse density. The CLI also applies zone humidity to local sound speed;
//! the fast path retains its documented dry-sound fallback.
//!
//! Validation gates:
//!  - Identity: a flat AtmoSock whose zones equal the station base → BIT-IDENTICAL trajectory.
//!  - Double-apply guard: flat AtmoSock + steep arc → byte-identical to the altitude-only
//!    (no-AtmoSock) result (proves the y-lapse still runs and a flat zone adds nothing); and a
//!    level shot + strong x temperature gradient → a measurable zone density effect with ~no
//!    altitude contribution.
//!  - Physical direction: hot (low-density) near zone → less drag / more retained velocity than a
//!    cold near zone, matching an offline CIPM hand calc.
//!  - Fast/cli parity: the same per-zone density inputs through cli_api and the fast path agree
//!    loosely on the endpoint (catches an incomplete fast_trajectory.rs density-freeze fix).

use ballistics_engine::atmosphere::{calculate_air_density_cimp, AtmoSock};
use ballistics_engine::fast_trajectory::{fast_integrate, FastIntegrationParams};
use ballistics_engine::wind::WindSock;
use ballistics_engine::{
    AtmosphericConditions, BallisticInputs, DragModel, TrajectorySolver, WindConditions,
};

const TARGET_M: f64 = 914.4; // 1000 yards
const GRAINS_TO_KG: f64 = 0.00006479891;
const IN_TO_M: f64 = 0.0254;

/// Reference .308, 175 gr projectile (G7 BC 0.243), SI-canonical.
fn ref_308_175gr(muzzle_angle_rad: f64) -> BallisticInputs {
    let mut i = BallisticInputs::default();
    i.muzzle_velocity = 792.5; // m/s (~2600 fps)
    i.bullet_mass = 175.0 * GRAINS_TO_KG;
    i.bullet_diameter = 0.308 * IN_TO_M;
    i.bullet_length = 1.24 * IN_TO_M;
    i.caliber_inches = 0.308;
    i.weight_grains = 175.0;
    i.bc_value = 0.243;
    i.bc_type = DragModel::G7;
    i.twist_rate = 10.0;
    i.muzzle_angle = muzzle_angle_rad;
    i.sight_height = 0.0;
    i.muzzle_height = 0.0;
    // Keep advanced lateral effects OFF so byte-identity comparisons are clean (drag-only).
    i.use_enhanced_spin_drift = false;
    i.enable_magnus = false;
    i.enable_coriolis = false;
    i
}

/// Non-default explicit station atmosphere so `resolve_station_conditions` returns it unchanged
/// (explicit temp/pressure are authoritative) — the base against which a "flat" AtmoSock must be
/// byte-identical.
fn station_atmo() -> AtmosphericConditions {
    AtmosphericConditions {
        temperature: 20.0, // Celsius (non-default)
        pressure: 1000.0,  // hPa (non-default)
        humidity: 50.0,    // percent
        altitude: 0.0,
    }
}

fn solve_with(
    inputs: &BallisticInputs,
    atmo: &AtmosphericConditions,
    segments: Option<Vec<(f64, f64, f64, f64)>>,
    max_range: f64,
) -> ballistics_engine::TrajectoryResult {
    let mut solver = TrajectorySolver::new(inputs.clone(), WindConditions::default(), atmo.clone());
    solver.set_max_range(max_range);
    if let Some(segs) = segments {
        solver.set_atmo_segments(segs);
    }
    solver.solve().expect("cli_api solve should succeed")
}

/// Vertical position (m, relative to bore line) at a downrange distance.
fn drop_at(result: &ballistics_engine::TrajectoryResult, range_m: f64) -> f64 {
    result
        .position_at_range(range_m)
        .unwrap_or_else(|| result.points.last().unwrap().position)
        .y
}

/// Impact velocity at the end of the (max_range-bounded) trajectory.
fn impact_velocity(result: &ballistics_engine::TrajectoryResult) -> f64 {
    result.points.last().unwrap().velocity_magnitude
}

// ------------------------------------------------------------------------------------------------
// Identity: a flat AtmoSock == station → BIT-IDENTICAL trajectory.
// ------------------------------------------------------------------------------------------------

#[test]
fn atmo_sock_identity_is_bit_identical() {
    let inputs = ref_308_175gr(0.008); // modest arc, reaches the target
    let atmo = station_atmo();

    let baseline = solve_with(&inputs, &atmo, None, TARGET_M);

    // Flat AtmoSock: every zone equals the station base (temp, pressure, humidity), held to well
    // beyond the target so a single zone covers the whole flight.
    let flat = vec![
        (20.0, 1000.0, 50.0, TARGET_M * 0.5),
        (20.0, 1000.0, 50.0, TARGET_M * 5.0),
    ];
    let with_sock = solve_with(&inputs, &atmo, Some(flat), TARGET_M);

    assert_eq!(
        baseline.points.len(),
        with_sock.points.len(),
        "flat AtmoSock must not change the number of integration points"
    );
    // Full-trajectory bit-identity: position (x,y,z), speed and time at every point.
    for (a, b) in baseline.points.iter().zip(with_sock.points.iter()) {
        assert_eq!(a.time.to_bits(), b.time.to_bits(), "time drifted");
        assert_eq!(a.position.x.to_bits(), b.position.x.to_bits(), "x drifted");
        assert_eq!(a.position.y.to_bits(), b.position.y.to_bits(), "y drifted");
        assert_eq!(a.position.z.to_bits(), b.position.z.to_bits(), "z drifted");
        assert_eq!(
            a.velocity_magnitude.to_bits(),
            b.velocity_magnitude.to_bits(),
            "velocity drifted"
        );
    }
    // And the headline scalars.
    assert_eq!(
        baseline.time_of_flight.to_bits(),
        with_sock.time_of_flight.to_bits(),
        "TOF must be bit-identical"
    );
    assert_eq!(
        baseline.impact_velocity.to_bits(),
        with_sock.impact_velocity.to_bits(),
        "impact velocity must be bit-identical"
    );
}

// ------------------------------------------------------------------------------------------------
// Double-apply guard (part A): flat AtmoSock + steep arc == altitude-only (no-AtmoSock).
// Proves the y-lapse (get_local_atmosphere) still runs and a flat zone contributes nothing extra.
// ------------------------------------------------------------------------------------------------

#[test]
fn atmo_sock_flat_zone_plus_steep_arc_matches_altitude_only() {
    // A steep launch so the bullet climbs tens of meters and the altitude lapse is materially
    // exercised over the flight.
    let inputs = ref_308_175gr(0.30); // ~17 deg
    let atmo = station_atmo();

    let altitude_only = solve_with(&inputs, &atmo, None, TARGET_M);

    // Sanity: the arc really does gain significant altitude, so the y-lapse is non-trivial here.
    let max_y = altitude_only
        .points
        .iter()
        .map(|p| p.position.y)
        .fold(f64::MIN, f64::max);
    assert!(
        max_y > 40.0,
        "steep arc should climb well above the muzzle (got max y = {max_y:.1} m)"
    );

    let flat = vec![(20.0, 1000.0, 50.0, TARGET_M * 5.0)];
    let flat_zone = solve_with(&inputs, &atmo, Some(flat), TARGET_M);

    assert_eq!(altitude_only.points.len(), flat_zone.points.len());
    for (a, b) in altitude_only.points.iter().zip(flat_zone.points.iter()) {
        assert_eq!(a.position.y.to_bits(), b.position.y.to_bits(), "y drifted");
        assert_eq!(
            a.velocity_magnitude.to_bits(),
            b.velocity_magnitude.to_bits(),
            "velocity drifted"
        );
    }
}

// ------------------------------------------------------------------------------------------------
// Double-apply guard (part B): level shot (y≈0) + strong x temperature gradient → a MEASURABLE
// zone density effect with ~no altitude contribution.
// ------------------------------------------------------------------------------------------------

#[test]
fn atmo_sock_level_shot_x_gradient_is_measurable() {
    let inputs = ref_308_175gr(0.0); // pure level fire: y stays within a few meters of the bore
    let atmo = station_atmo();

    // Confirm the altitude excursion is tiny (so the effect below is the X gradient, not the lapse).
    let baseline = solve_with(&inputs, &atmo, None, TARGET_M);
    let max_abs_y = baseline
        .points
        .iter()
        .map(|p| p.position.y.abs())
        .fold(0.0_f64, f64::max);
    assert!(
        max_abs_y < 12.0,
        "level shot should stay near the bore line (|y|max = {max_abs_y:.1} m)"
    );

    // Strong downrange temperature gradient: a hot (thin) near half and a cold (dense) far half,
    // both at the station pressure/humidity so ONLY temperature (hence density) varies with X.
    let gradient = vec![
        (45.0, 1000.0, 50.0, TARGET_M * 0.5), // hot near
        (-20.0, 1000.0, 50.0, TARGET_M * 5.0), // cold far
    ];
    let with_gradient = solve_with(&inputs, &atmo, Some(gradient), TARGET_M);

    let v_base = impact_velocity(&baseline);
    let v_grad = impact_velocity(&with_gradient);
    // A measurable change in retained velocity (well above numerical noise).
    assert!(
        (v_grad - v_base).abs() > 1.0,
        "x temperature gradient should measurably change retained velocity: base {v_base:.2} m/s \
         vs gradient {v_grad:.2} m/s"
    );
}

// ------------------------------------------------------------------------------------------------
// Physical direction: hot near zone (low density) → less drag / more retained velocity than a
// cold near zone. Matches an offline CIPM density hand-calc.
// ------------------------------------------------------------------------------------------------

#[test]
fn atmo_sock_hot_near_retains_more_velocity_than_cold_near() {
    let inputs = ref_308_175gr(0.008);
    let atmo = station_atmo();

    // Offline CIPM sanity: hot air at the same pressure/humidity is less dense than cold air.
    let rho_hot = calculate_air_density_cimp(45.0, 1000.0, 50.0);
    let rho_cold = calculate_air_density_cimp(-20.0, 1000.0, 50.0);
    assert!(
        rho_hot < rho_cold,
        "CIPM: hot air ({rho_hot:.4}) must be less dense than cold air ({rho_cold:.4})"
    );

    // Hot near / standard far vs cold near / standard far. The near zone covers the high-velocity
    // stretch where drag (∝ v²·ρ) dominates, so density there drives retained velocity.
    let hot_near = vec![
        (45.0, 1000.0, 50.0, TARGET_M * 0.5),
        (20.0, 1000.0, 50.0, TARGET_M * 5.0),
    ];
    let cold_near = vec![
        (-20.0, 1000.0, 50.0, TARGET_M * 0.5),
        (20.0, 1000.0, 50.0, TARGET_M * 5.0),
    ];

    let v_hot = impact_velocity(&solve_with(&inputs, &atmo, Some(hot_near), TARGET_M));
    let v_cold = impact_velocity(&solve_with(&inputs, &atmo, Some(cold_near), TARGET_M));

    assert!(
        v_hot > v_cold,
        "hot (thin) near zone should retain MORE velocity than cold (dense): hot {v_hot:.2} m/s \
         vs cold {v_cold:.2} m/s"
    );
}

// ------------------------------------------------------------------------------------------------
// Fast/cli parity: the same per-zone density inputs through cli_api and the fast path agree loosely
// on the endpoint. If the fast density-freeze (fast_trajectory.rs :616) were not fixed, the zone
// would have NO effect on the fast path and this would diverge. The CLI uses real zone RH for
// sound speed; the fast path retains its documented dry-sound fallback.
// ------------------------------------------------------------------------------------------------

/// Run the plain fast_integrate path with an optional AtmoSock; return endpoint velocity and y.
fn fast_endpoint(
    inputs: &BallisticInputs,
    atmo: &AtmosphericConditions,
    segments: Option<Vec<(f64, f64, f64, f64)>>,
) -> (f64, f64) {
    // Mirror monte_carlo.rs base-density derivation (standard-mode atmo_params).
    let (resolved_temp_c, resolved_pressure_hpa) =
        ballistics_engine::atmosphere::resolve_station_conditions(
            atmo.temperature,
            atmo.pressure,
            atmo.altitude,
        );
    let (air_density, _) = ballistics_engine::atmosphere::calculate_atmosphere(
        atmo.altitude,
        Some(resolved_temp_c),
        Some(resolved_pressure_hpa),
        atmo.humidity,
    );
    let base_ratio = air_density / 1.225;

    let v = inputs.muzzle_velocity;
    let a = inputs.muzzle_angle;
    let initial_state = [0.0, inputs.sight_height, 0.0, v * a.cos(), v * a.sin(), 0.0];
    let params = FastIntegrationParams {
        initial_state,
        t_span: (0.0, 30.0),
        horiz: TARGET_M,
        vert: 0.0,
        atmo_params: (atmo.altitude, resolved_temp_c, resolved_pressure_hpa, base_ratio),
        atmo_sock: segments.map(AtmoSock::new),
    };
    let sol = fast_integrate(inputs, &WindSock::new(vec![]), params);
    assert!(sol.success, "fast_integrate should succeed");
    let n = sol.t.len();
    let vx = sol.y[3][n - 1];
    let vy = sol.y[4][n - 1];
    let vz = sol.y[5][n - 1];
    ((vx * vx + vy * vy + vz * vz).sqrt(), sol.y[1][n - 1])
}

#[test]
fn atmo_sock_fast_path_responds_to_zone() {
    // The density-freeze fix means the plain fast path must now let a downrange zone change drag.
    let inputs = ref_308_175gr(0.008);
    let atmo = station_atmo();

    let hot_near = vec![
        (45.0, 1000.0, 50.0, TARGET_M * 0.5),
        (20.0, 1000.0, 50.0, TARGET_M * 5.0),
    ];

    let (v_no_zone, _) = fast_endpoint(&inputs, &atmo, None);
    let (v_hot, _) = fast_endpoint(&inputs, &atmo, Some(hot_near));

    // A hot near zone (thin air, less drag) must raise the fast-path retained velocity. If :616
    // still froze density, the zone would be inert here.
    assert!(
        v_hot > v_no_zone + 0.5,
        "fast path must respond to the hot zone (frozen-density regression): no-zone {v_no_zone:.2} \
         m/s vs hot {v_hot:.2} m/s"
    );
}

#[test]
fn atmo_sock_fast_and_cli_agree_on_zoned_endpoint() {
    let inputs = ref_308_175gr(0.008);
    let atmo = station_atmo();

    let hot_near = vec![
        (45.0, 1000.0, 50.0, TARGET_M * 0.5),
        (20.0, 1000.0, 50.0, TARGET_M * 5.0),
    ];

    let (v_fast, _) = fast_endpoint(&inputs, &atmo, Some(hot_near.clone()));

    let cli = solve_with(&inputs, &atmo, Some(hot_near), TARGET_M);
    let v_cli = impact_velocity(&cli);

    // Two independent integrators (fast fixed-step RK4 vs cli_api adaptive RK45) with the same
    // per-zone density base should agree within a few percent despite their documented sound-speed
    // humidity difference.
    let rel = (v_fast - v_cli).abs() / v_cli;
    assert!(
        rel < 0.05,
        "fast vs cli zoned endpoint velocity disagree by {:.1}% (fast {v_fast:.2}, cli {v_cli:.2})",
        rel * 100.0
    );
}

// ------------------------------------------------------------------------------------------------
// Reference-shot report: .308 175 gr, 1000 yd — drop with (a) no AtmoSock, (b) flat AtmoSock ==
// station, (c) hot-near / cold-far AtmoSock. Prints the table and asserts (a)==(b), (c) differs.
// ------------------------------------------------------------------------------------------------

#[test]
fn reference_shot_drop_table() {
    let inputs = ref_308_175gr(0.008);
    let atmo = station_atmo();

    let a = solve_with(&inputs, &atmo, None, TARGET_M);
    let b = solve_with(
        &inputs,
        &atmo,
        Some(vec![(20.0, 1000.0, 50.0, TARGET_M * 5.0)]),
        TARGET_M,
    );
    let c = solve_with(
        &inputs,
        &atmo,
        Some(vec![
            (45.0, 1000.0, 50.0, TARGET_M * 0.5), // hot near
            (-20.0, 1000.0, 50.0, TARGET_M * 5.0), // cold far
        ]),
        TARGET_M,
    );

    let drop_a = drop_at(&a, TARGET_M);
    let drop_b = drop_at(&b, TARGET_M);
    let drop_c = drop_at(&c, TARGET_M);
    let to_in = |m: f64| m / IN_TO_M;

    println!("MBA-1137 reference shot (.308 175gr G7 0.243, 1000 yd, station 20C/1000hPa/50%RH):");
    println!(
        "  (a) no AtmoSock          : drop {:+.4} m ({:+.2} in), v_impact {:.2} m/s, TOF {:.4} s",
        drop_a,
        to_in(drop_a),
        impact_velocity(&a),
        a.time_of_flight
    );
    println!(
        "  (b) flat AtmoSock=station: drop {:+.4} m ({:+.2} in), v_impact {:.2} m/s, TOF {:.4} s",
        drop_b,
        to_in(drop_b),
        impact_velocity(&b),
        b.time_of_flight
    );
    println!(
        "  (c) hot-near/cold-far    : drop {:+.4} m ({:+.2} in), v_impact {:.2} m/s, TOF {:.4} s",
        drop_c,
        to_in(drop_c),
        impact_velocity(&c),
        c.time_of_flight
    );

    // (a) == (b): a flat zone == station changes nothing.
    assert_eq!(
        drop_a.to_bits(),
        drop_b.to_bits(),
        "(a) no AtmoSock and (b) flat AtmoSock must be bit-identical"
    );
    // (c) differs sensibly from (a).
    assert!(
        (drop_c - drop_a).abs() > 1e-3,
        "(c) zoned atmosphere should change the drop vs (a): {drop_a:.4} m vs {drop_c:.4} m"
    );
}
