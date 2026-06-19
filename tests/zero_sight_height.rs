// MBA-951: a zero must target the LINE-OF-SIGHT height (sight_height) at the zero distance, not
// the bore line. That is the CLI's convention in every zero call ("target height = LOS height").
// The WASM standalone zero handler was passing 0.0 — a bore-line zero that ignored sight height.
// wasm.rs is cfg(target_arch="wasm32") and can't be compiled here, so this validates the
// canonical calculate_zero_angle_with_conditions the handler calls.
use ballistics_engine::{
    calculate_zero_angle_with_conditions, AtmosphericConditions, BallisticInputs, DragModel,
    WindConditions,
};

fn zero_inputs() -> BallisticInputs {
    let mut i = BallisticInputs::default();
    i.muzzle_velocity = 800.0; // m/s
    i.bc_value = 0.475;
    i.bc_type = DragModel::G1;
    i.bullet_mass = 168.0 * 0.00006479891; // kg
    i.bullet_diameter = 0.308 * 0.0254; // m
    i.sight_height = 0.05; // 50 mm above the bore
    i
}

#[test]
fn zero_targets_line_of_sight_not_bore() {
    let i = zero_inputs();
    let dist = 91.44; // 100 yd
    let sight = calculate_zero_angle_with_conditions(
        i.clone(),
        dist,
        i.sight_height, // LOS height (the fix)
        WindConditions::default(),
        AtmosphericConditions::default(),
    )
    .unwrap();
    let bore = calculate_zero_angle_with_conditions(
        i.clone(),
        dist,
        0.0, // the old WASM bug: bore-line zero
        WindConditions::default(),
        AtmosphericConditions::default(),
    )
    .unwrap();
    // Reaching the LOS height (above the bore) at the target needs a higher launch angle than
    // returning to the bore line, and the gap is the sight-height angle (~0.5 mrad here).
    assert!(
        sight > bore,
        "sight-line zero angle ({sight}) should exceed bore-line ({bore})"
    );
    assert!(
        (sight - bore) > 1e-4,
        "sight vs bore zero should differ measurably (got {})",
        sight - bore
    );
}

#[test]
fn zero_angle_invariant_to_muzzle_height() {
    // Heights above ground cancel: the zero angle is the same whether the muzzle is at ground
    // level or elevated, as long as the LOS target tracks it — confirming the handler is right to
    // ignore --muzzle-height/--target-height for the angle.
    let mk = |mh: f64| {
        let mut i = zero_inputs();
        i.muzzle_height = mh;
        calculate_zero_angle_with_conditions(
            i.clone(),
            91.44,
            mh + i.sight_height,
            WindConditions::default(),
            AtmosphericConditions::default(),
        )
        .unwrap()
    };
    assert!(
        (mk(0.0) - mk(1.5)).abs() < 1e-7,
        "zero angle should be invariant to muzzle height: {} vs {}",
        mk(0.0),
        mk(1.5)
    );
}
