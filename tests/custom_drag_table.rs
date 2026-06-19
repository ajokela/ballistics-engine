// MBA-940 regression guard: a user-supplied custom_drag_table must be HONORED. Previously it was
// plumbed onto BallisticInputs but never read — every solver used the G-model and silently
// ignored the user's curve. This drives the cli_api TrajectorySolver with a flat, very-low-drag
// custom curve and asserts it retains materially more velocity than the G7 model over the same
// shot (i.e. the custom Cd is actually used).
use ballistics_engine::drag::DragTable;
use ballistics_engine::{
    AtmosphericConditions, BallisticInputs, DragModel, TrajectorySolver, WindConditions,
};

fn base() -> BallisticInputs {
    let mut i = BallisticInputs::default();
    i.muzzle_velocity = 800.0; // m/s
    i.bc_value = 0.5;
    i.bc_type = DragModel::G7;
    i.bullet_mass = 168.0 * 0.00006479891; // kg
    i.bullet_diameter = 0.308 * 0.0254; // m
    i.bullet_length = 1.215 * 0.0254;
    i.caliber_inches = 0.308;
    i.weight_grains = 168.0;
    i.temperature = 15.0;
    i.pressure = 1013.25;
    i
}

fn velocity_at(inputs: BallisticInputs, range_m: f64) -> f64 {
    let mut s = TrajectorySolver::new(
        inputs,
        WindConditions::default(),
        AtmosphericConditions::default(),
    );
    s.set_max_range(range_m * 1.2);
    let r = s.solve().expect("solve");
    for p in &r.points {
        if p.position[0] >= range_m {
            return p.velocity_magnitude;
        }
    }
    r.points.last().unwrap().velocity_magnitude
}

#[test]
fn custom_drag_table_is_honored() {
    let range = 300.0;
    let g7_vel = velocity_at(base(), range);

    // Flat, very-low-drag custom curve (Cd ~0.05 across the Mach range) -> minimal velocity loss.
    let mut custom = base();
    custom.custom_drag_table = Some(DragTable::new(
        vec![0.0, 1.0, 2.0, 5.0],
        vec![0.05, 0.05, 0.05, 0.05],
    ));
    let custom_vel = velocity_at(custom, range);

    assert!(
        custom_vel > g7_vel + 50.0,
        "custom low-drag table should retain materially more velocity than G7 \
         (custom={custom_vel:.1} m/s, g7={g7_vel:.1} m/s) — table was likely ignored"
    );
}
