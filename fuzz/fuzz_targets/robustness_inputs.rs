#![no_main]
use arbitrary::Unstructured;
use ballistics_engine::{TrajectorySolver, WindConditions};
use ballistics_engine_fuzz::domain;
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    // Exercise the solve-boundary validation gate with NaN, infinities, non-positive values, and
    // finite extremes. Returning Err is valid; every successful result must remain finite/sane.
    let Ok(inputs) = domain::hostile_inputs(&mut u) else {
        return;
    };
    // The air the draw flies in: its own wild air, a plausible station carrying its wild
    // humidity, a band edge, a hot saturated corner, or the default air (see
    // `domain::hostile_air`). `TrajectorySolver` takes the humidity from the atmosphere, so a
    // fixed default atmosphere would discard every drawn humidity before the validator.
    let Ok(air) = domain::hostile_air(&mut u, &inputs) else {
        return;
    };
    let max_range = inputs.target_distance;
    let mut solver = TrajectorySolver::new(inputs, WindConditions::default(), air);
    solver.set_max_range(max_range);
    // `solve()` returning Err is a VALID outcome (input cleanly rejected).
    // A panic, or NaN/Inf/absurd values inside an Ok result, is the bug.
    if let Ok(result) = solver.solve() {
        domain::assert_finite_sane(&result);
    }
});
