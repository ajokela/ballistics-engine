//! MBA-1415: input conditioning must reach every public entry point into integration.
//!
//! `TrajectorySolver::new` is not the only way into the solver. `fast_integrate` and
//! `fast_integrate_with_segments` are `pub` and are called directly by the Python binding, so
//! conditioning that lives only in the solver constructor is invisible to those consumers: the
//! caller sets a field, gets no error, and the value is silently ignored.
//!
//! That failure shape has already cost this project twice — MBA-1296 (a dropped field that
//! zeroed Coriolis in production) and the MBA-1356 review catch (`cd_scale` dropped by the
//! segmented fast path). Both were found by review rather than by tests, which is precisely
//! what this file exists to change.
//!
//! The BC reference standard is the conditioning step used as the probe here because it is the
//! one with a non-default value today. It is a stand-in for the class: if a future step is added
//! to `normalize_for_solve` and one entry point stops calling it, these tests fail.

use ballistics_engine::atmosphere::{calculate_atmosphere, resolve_station_conditions};
use ballistics_engine::fast_trajectory::{
    fast_integrate, fast_integrate_with_segments, FastIntegrationParams,
};
use ballistics_engine::wind::WindSock;
use ballistics_engine::{
    AtmosphericConditions, BallisticInputs, BcReferenceStandard, DragModel, TrajectorySolver,
    WindConditions,
};

const TARGET_M: f64 = 457.2; // 500 yd

/// Parity tolerance, m/s. Chosen from measurement, not taste: at this range the two integrator
/// families differ by 0.14 m/s on identical inputs, while skipping the BC reference conversion
/// moves impact velocity by 4.15 m/s. 1.0 sits an order of magnitude above the former and well
/// below the latter, so the test tolerates integrator differences and nothing else. The
/// companion test below asserts that separation still holds, so this number cannot silently
/// become too loose to catch the bug it exists for.
const PARITY_TOL_MPS: f64 = 1.0;

fn asm_inputs() -> BallisticInputs {
    BallisticInputs {
        muzzle_velocity: 800.0,
        bullet_mass: 175.0 * 0.00006479891,
        bullet_diameter: 0.308 * 0.0254,
        bc_value: 0.243,
        bc_type: DragModel::G7,
        // The conditioning step under test: this BC is declared against Army Standard Metro and
        // must be converted to ICAO before any retardation math runs, on EVERY path.
        bc_reference_standard: BcReferenceStandard::ArmyStandardMetro,
        ..Default::default()
    }
}

/// Terminal velocity via the fast single-shot path, which external bindings call directly.
fn fast_impact_velocity(i: &BallisticInputs, segmented: bool) -> f64 {
    let (resolved_temp_c, resolved_pressure_hpa) =
        resolve_station_conditions(i.temperature, i.pressure, i.altitude);
    let (air_density, _) = calculate_atmosphere(
        i.altitude,
        Some(resolved_temp_c),
        Some(resolved_pressure_hpa),
        i.humidity_percent(),
    );
    let v = i.muzzle_velocity;
    let a = i.muzzle_angle;
    let params = FastIntegrationParams {
        initial_state: [0.0, i.muzzle_height, 0.0, v * a.cos(), v * a.sin(), 0.0],
        t_span: (0.0, 30.0),
        horiz: TARGET_M,
        vert: i.muzzle_height + i.sight_height,
        atmo_params: (
            i.altitude,
            resolved_temp_c,
            resolved_pressure_hpa,
            air_density / 1.225,
        ),
        atmo_sock: None,
    };
    let sol = if segmented {
        fast_integrate_with_segments(i, vec![], params)
    } else {
        fast_integrate(i, &WindSock::new(vec![]), params)
    };
    assert!(sol.success, "fast path should succeed");
    let n = sol.t.len();
    (sol.y[3][n - 1].powi(2) + sol.y[4][n - 1].powi(2) + sol.y[5][n - 1].powi(2)).sqrt()
}

/// Terminal velocity via the constructor path (CLI, WASM, FFI).
fn solver_impact_velocity(i: &BallisticInputs) -> f64 {
    let mut s = TrajectorySolver::new(
        i.clone(),
        WindConditions::default(),
        AtmosphericConditions::default(),
    );
    s.set_max_range(TARGET_M);
    let r = s.solve().expect("solve");
    r.points.last().expect("points").velocity_magnitude
}

/// Applying the conditioning twice must equal applying it once. Monte Carlo builds a solver per
/// sample and the fast entry points may receive inputs that already passed through one, so a
/// step that scales a value in place without disarming itself would double-apply in silence.
#[test]
fn normalization_is_idempotent() {
    let mut once = asm_inputs();
    once.normalize_for_solve();

    let mut twice = asm_inputs();
    twice.normalize_for_solve();
    twice.normalize_for_solve();

    assert_eq!(
        once.bc_value, twice.bc_value,
        "the BC was scaled a second time: {} then {}",
        once.bc_value, twice.bc_value
    );
    assert_eq!(once.muzzle_velocity, twice.muzzle_velocity);
    assert_eq!(once.caliber_inches, twice.caliber_inches);
    assert_eq!(once.weight_grains, twice.weight_grains);
}

/// The conversion must actually change the BC, or every other assertion here is vacuous.
#[test]
fn the_probe_step_is_not_a_no_op() {
    let mut normalized = asm_inputs();
    normalized.normalize_for_solve();
    assert!(
        (normalized.bc_value - 0.243).abs() > 1e-6,
        "an ASM-referenced BC must be converted; it stayed at {}",
        normalized.bc_value
    );
    assert!(
        matches!(
            normalized.bc_reference_standard,
            BcReferenceStandard::Icao
        ),
        "the step must disarm itself so a second application is a no-op"
    );
}

/// The whole point: identical inputs with a non-default conditioning step must produce the same
/// physics through the constructor path and through the fast path bindings use.
#[test]
fn every_entry_point_applies_the_same_conditioning() {
    let inputs = asm_inputs();

    let via_solver = solver_impact_velocity(&inputs);
    let via_fast = fast_impact_velocity(&inputs, false);
    let via_fast_segments = fast_impact_velocity(&inputs, true);

    // These are different integrators, so they are not expected to agree bit-for-bit; the
    // failure this guards against is one path never converting the BC at all, which moves the
    // impact velocity far more than integrator differences do.
    let tol = PARITY_TOL_MPS;
    assert!(
        (via_solver - via_fast).abs() < tol,
        "fast_integrate disagrees with the solver: {via_solver} vs {via_fast} m/s — \
         the fast path is probably not applying the BC reference conversion"
    );
    assert!(
        (via_solver - via_fast_segments).abs() < tol,
        "fast_integrate_with_segments disagrees with the solver: {via_solver} vs \
         {via_fast_segments} m/s"
    );
}

/// Guard the guard. A parity test is only worth its runtime if the bug it targets would actually
/// breach its tolerance — and the first draft of this file used a tolerance nearly twice the
/// effect size, so it would have passed with the conversion dropped entirely. This pins the
/// separation: skipping the conversion must move impact velocity by more than the tolerance the
/// parity test allows.
#[test]
fn the_parity_tolerance_is_tight_enough_to_catch_a_skipped_conversion() {
    let converted = solver_impact_velocity(&asm_inputs());

    let mut unconverted = asm_inputs();
    unconverted.bc_reference_standard = BcReferenceStandard::Icao;
    let unconverted = solver_impact_velocity(&unconverted);

    assert!(
        (converted - unconverted).abs() > PARITY_TOL_MPS,
        "skipping the conversion moves impact velocity by only {:.2} m/s, which is inside the \
         parity tolerance — the parity test would not catch a dropped conditioning step",
        (converted - unconverted).abs()
    );
}

