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

// Regression guard: a custom drag table's Cd is the projectile's ACTUAL drag coefficient, so
// the retardation denominator is the SECTIONAL DENSITY (lb/in²) derived from the bullet's
// mass/diameter — NOT bc_value (Cd_own / SD == Cd_ref / BC). Previously all solver paths
// divided the curve's Cd by bc_value, so custom-table trajectories wrongly scaled with
// whatever BC happened to be set.
#[test]
fn custom_drag_table_is_bc_invariant() {
    let table = DragTable::new(
        vec![0.0, 0.8, 1.0, 1.2, 2.0, 5.0],
        vec![0.20, 0.22, 0.40, 0.38, 0.30, 0.25],
    );
    let range = 500.0;

    let run = |bc: f64, with_table: bool| -> (f64, f64) {
        let mut i = base();
        i.bc_value = bc;
        if with_table {
            i.custom_drag_table = Some(table.clone());
        }
        let mut s = TrajectorySolver::new(
            i,
            WindConditions::default(),
            AtmosphericConditions::default(),
        );
        s.set_max_range(range * 1.2);
        let r = s.solve().expect("solve");
        for p in &r.points {
            if p.position[0] >= range {
                return (p.position[1], p.velocity_magnitude);
            }
        }
        let p = r.points.last().unwrap();
        (p.position[1], p.velocity_magnitude)
    };

    let (drop_a, vel_a) = run(0.505, true);
    let (drop_b, vel_b) = run(0.264, true);
    let (drop_c, vel_c) = run(0.150, true);

    assert!(
        (drop_a - drop_b).abs() < 1e-9 && (drop_b - drop_c).abs() < 1e-9,
        "trajectory with a custom drag table must not depend on bc_value \
         (drop: bc=0.505 -> {drop_a}, bc=0.264 -> {drop_b}, bc=0.150 -> {drop_c})"
    );
    assert!(
        (vel_a - vel_b).abs() < 1e-9 && (vel_b - vel_c).abs() < 1e-9,
        "velocity with a custom drag table must not depend on bc_value \
         (bc=0.505 -> {vel_a}, bc=0.264 -> {vel_b}, bc=0.150 -> {vel_c})"
    );

    // And the table run must actually differ from the G7 (no-table) run.
    let (drop_g7, vel_g7) = run(0.505, false);
    assert!(
        (drop_a - drop_g7).abs() > 1e-3 || (vel_a - vel_g7).abs() > 1.0,
        "custom-table run should differ from the no-table G7 run \
         (drop {drop_a} vs {drop_g7}, vel {vel_a} vs {vel_g7})"
    );
}

// MBA-1285: feeding the shipped G7 deck back in as a *custom* table for a projectile whose BC
// is set equal to its own sectional density (form factor i == 1: BC = SD / i) should reproduce
// the built-in G7 model's trajectory closely — Cd_own / SD == Cd_ref / BC when the two curves
// and denominators coincide (see BallisticInputs::custom_drag_denominator). The custom-table
// path skips the transonic shape correction the G-model path applies, so this is a
// self-consistency check, not a bit-for-bit equality guarantee.
#[test]
fn g7_deck_matches_g7_model_for_unit_form_factor() {
    use ballistics_engine::drag::DragTable;
    let g7 = DragTable::from_csv_str(include_str!("../data/g7.csv")).unwrap();

    let mut model = base();
    model.bc_type = ballistics_engine::DragModel::G7;

    // Choose mass/diameter (already set by base()) so sectional_density_lb_in2 == bc_value
    // (i == 1).
    let mut deck = model.clone();
    let sd = model
        .sectional_density_lb_in2()
        .expect("base() sets mass+diameter");
    model.bc_value = sd; // BC == SD => i == 1
    deck.bc_value = sd;
    deck.custom_drag_table = Some(g7);

    let v_model = velocity_at(model, 900.0);
    let v_deck = velocity_at(deck, 900.0);
    let rel = (v_model - v_deck).abs() / v_model;
    assert!(
        rel < 0.02,
        "G7 deck vs G7 model diverged {:.3}% at 900 m (model={v_model:.3} m/s, deck={v_deck:.3} m/s)",
        rel * 100.0
    );
}

// Integration-level guard on top of the unit tests in src/drag.rs: from_csv_str must reject
// malformed decks via the *file-parsing* path (not just try_new's already-parsed-values path),
// covering a negative Cd, a descending Mach axis, and an unparseable row.
#[test]
fn from_csv_str_rejects_bad_decks() {
    use ballistics_engine::drag::DragTable;
    assert!(DragTable::from_csv_str("0.5,0.2\n1.0,-0.1\n").is_err()); // negative Cd
    assert!(DragTable::from_csv_str("1.0,0.3\n0.5,0.3\n").is_err()); // descending Mach
    assert!(DragTable::from_csv_str("0.5,0.3\nbroken\n").is_err()); // malformed row
}

// The SD denominator itself: 168 gr / 7000 / 0.308^2 in² ≈ 0.25305 lb/in², derived from
// either the imperial mirror fields or the SI mass/diameter fallback; None when both are
// unusable (that case falls back to bc_value inside the solvers, with a warning).
#[test]
fn sectional_density_denominator_sources() {
    let i = base();
    let sd = i.sectional_density_lb_in2().expect("SD from imperial fields");
    assert!((sd - 168.0 / 7000.0 / (0.308 * 0.308)).abs() < 1e-12);

    // SI-only caller: imperial mirrors zeroed, SI kg/meters populated.
    let mut si = base();
    si.caliber_inches = 0.0;
    si.weight_grains = 0.0;
    let sd_si = si.sectional_density_lb_in2().expect("SD from SI fallback");
    assert!((sd_si - sd).abs() < 1e-9);

    // Degenerate inputs -> None (solvers fall back to bc_value instead of panicking).
    let mut none = base();
    none.caliber_inches = 0.0;
    none.weight_grains = 0.0;
    none.bullet_mass = 0.0;
    none.bullet_diameter = 0.0;
    assert!(none.sectional_density_lb_in2().is_none());
    assert_eq!(none.custom_drag_denominator(0.42), 0.42);
}
