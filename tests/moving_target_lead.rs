// MBA-1287 acceptance gates for the moving-target lead calculation.
use ballistics_engine::{
    calculate_lead, AtmosphericConditions, BallisticInputs, DragModel, LeadError,
    TrajectorySolver, WindConditions,
};

fn base() -> BallisticInputs {
    let mut i = BallisticInputs::default();
    i.muzzle_velocity = 800.0;
    i.bc_value = 0.5;
    i.bc_type = DragModel::G7;
    i.bullet_mass = 0.0109;
    i.bullet_diameter = 0.00782;
    i.bullet_length = 0.0309;
    i.sight_height = 0.05;
    i.use_rk4 = true;
    i
}

/// TOF at exactly x, interpolated from an independent solve of the same inputs.
fn independent_tof(inputs: BallisticInputs, wind: WindConditions, x: f64, max_range: f64) -> f64 {
    let mut s = TrajectorySolver::new(inputs, wind, AtmosphericConditions::default());
    s.set_max_range(max_range);
    let r = s.solve().expect("solve");
    let pts = &r.points;
    for w in pts.windows(2) {
        if w[1].position.x >= x {
            let dx = w[1].position.x - w[0].position.x;
            let t = if dx.abs() < 1e-12 { 0.0 } else { (x - w[0].position.x) / dx };
            return w[0].time + t * (w[1].time - w[0].time);
        }
    }
    panic!("never reached {x}");
}

#[test]
fn perpendicular_lead_matches_analytic_identity() {
    // Ticket: lead_mil == v*TOF/R*1000 within 1e-6, TOF matches the solve endpoint.
    let (v, r) = (5.0, 500.0);
    let sol = calculate_lead(
        base(),
        WindConditions::default(),
        AtmosphericConditions::default(),
        v,
        90.0,
        r,
    )
    .expect("perpendicular solve");
    assert_eq!(sol.iterations, 0, "pure crossing must not iterate");
    assert!((sol.corrected_range_m - r).abs() < 1e-12);
    let expected_mil = v * sol.time_of_flight_s / r * 1000.0;
    assert!(
        (sol.lead_mil - expected_mil).abs() < 1e-6,
        "lead_mil {} vs analytic {expected_mil}",
        sol.lead_mil
    );
    // TOF cross-checked against an INDEPENDENT solve at the same range.
    let tof = independent_tof(base(), WindConditions::default(), r, r + 20.0);
    assert!(
        (sol.time_of_flight_s - tof).abs() < 1e-9,
        "TOF {} vs independent {tof}",
        sol.time_of_flight_s
    );
}

#[test]
fn outbound_45_converges_under_10_iterations_with_small_residual() {
    let sol = calculate_lead(
        base(),
        WindConditions::default(),
        AtmosphericConditions::default(),
        15.0,
        45.0,
        600.0,
    )
    .expect("outbound");
    assert!(sol.iterations < 10, "iterations {}", sol.iterations);
    // Residual: re-applying the fixed-point map moves the range by < 0.1 m.
    let tof_at_corrected =
        independent_tof(base(), WindConditions::default(), sol.corrected_range_m, 800.0);
    let reapplied = 600.0 + 15.0 * 45f64.to_radians().cos() * tof_at_corrected;
    assert!(
        (reapplied - sol.corrected_range_m).abs() < 0.1,
        "residual {} m",
        (reapplied - sol.corrected_range_m).abs()
    );
    assert!(sol.corrected_range_m > 600.0);

    // The reported solution must be evaluated AT the corrected intercept range,
    // not the original 600 m — this is the whole point of the range correction.
    assert!(
        (sol.time_of_flight_s - tof_at_corrected).abs() < 1e-9,
        "sol.time_of_flight_s {} vs independent TOF at corrected range {tof_at_corrected}",
        sol.time_of_flight_s
    );
    // Discriminating-power precondition: TOF at the corrected range must actually
    // differ from TOF at the original range, or the assertion above can't tell
    // them apart. If this ever fails, the scenario needs a faster/more radial
    // target, not a looser bound.
    let tof_at_original = independent_tof(base(), WindConditions::default(), 600.0, 800.0);
    assert!(
        (tof_at_corrected - tof_at_original).abs() > 1e-4,
        "test precondition: corrected vs original TOF must be distinguishable ({tof_at_corrected} vs {tof_at_original})"
    );

    // Pure recomputation pin: lead_m must equal v_lateral * reported TOF.
    let v_lateral = 15.0 * 45f64.to_radians().sin();
    assert!(
        (sol.lead_m - v_lateral * sol.time_of_flight_s).abs() < 1e-9,
        "lead_m {} vs v_lateral*TOF {}",
        sol.lead_m,
        v_lateral * sol.time_of_flight_s
    );
}

#[test]
fn inbound_target_shortens_corrected_range() {
    let sol = calculate_lead(
        base(),
        WindConditions::default(),
        AtmosphericConditions::default(),
        15.0,
        135.0, // inbound-crossing: cos(135°) < 0
        600.0,
    )
    .expect("inbound");
    assert!(sol.corrected_range_m < 600.0);
    assert!(sol.lead_m > 0.0, "sin(135°) > 0: still leads right");

    // The reported solution must be evaluated AT the corrected (shortened) intercept
    // range, not the original 600 m.
    let tof_at_corrected =
        independent_tof(base(), WindConditions::default(), sol.corrected_range_m, 700.0);
    assert!(
        (sol.time_of_flight_s - tof_at_corrected).abs() < 1e-9,
        "sol.time_of_flight_s {} vs independent TOF at corrected range {tof_at_corrected}",
        sol.time_of_flight_s
    );
}

#[test]
fn over_closing_target_returns_typed_error() {
    let r = calculate_lead(
        base(),
        WindConditions::default(),
        AtmosphericConditions::default(),
        2000.0,
        180.0,
        600.0,
    );
    match r {
        Err(LeadError::TargetOvertakesShooter { .. })
        | Err(LeadError::Convergence { .. })
        | Err(LeadError::BeyondSolvedSpan { .. }) => {}
        other => panic!("expected typed error, got {other:?}"),
    }
}

#[test]
fn zero_speed_gives_zero_lead() {
    let sol = calculate_lead(
        base(),
        WindConditions::default(),
        AtmosphericConditions::default(),
        0.0,
        90.0,
        500.0,
    )
    .expect("zero speed");
    assert_eq!(sol.lead_m, 0.0);
    assert_eq!(sol.lead_mil, 0.0);
    assert_eq!(sol.lead_moa, 0.0);
}

#[test]
fn moa_over_mil_is_exactly_3_438() {
    let sol = calculate_lead(
        base(),
        WindConditions::default(),
        AtmosphericConditions::default(),
        5.0,
        90.0,
        500.0,
    )
    .expect("solve");
    assert!((sol.lead_moa / sol.lead_mil - 3.438).abs() < 1e-12);
}

#[test]
fn tof_is_wind_aware_so_lead_follows() {
    // Strong headwind (blowing FROM downrange toward the shooter) slows the bullet:
    // TOF grows, so the required lead grows with it.
    let calm = calculate_lead(
        base(),
        WindConditions::default(),
        AtmosphericConditions::default(),
        5.0,
        90.0,
        600.0,
    )
    .expect("calm");
    let headwind = WindConditions {
        speed: 15.0,
        direction: 0.0, // engine convention post-0.19.0: 0 = headwind
        ..Default::default()
    };
    let windy = calculate_lead(
        base(),
        headwind,
        AtmosphericConditions::default(),
        5.0,
        90.0,
        600.0,
    )
    .expect("windy");
    assert!(
        windy.time_of_flight_s > calm.time_of_flight_s,
        "headwind must increase TOF: {} vs {}",
        windy.time_of_flight_s,
        calm.time_of_flight_s
    );
    assert!(windy.lead_mil > calm.lead_mil, "longer TOF must increase lead");
}
