#![no_main]
use libfuzzer_sys::fuzz_target;
use arbitrary::Unstructured;
use ballistics_engine::{AtmosphericConditions, TrajectorySolver, WindConditions};
use ballistics_engine_fuzz::domain;

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let Ok(inputs) = domain::hostile_inputs(&mut u) else { return };
    let solver = TrajectorySolver::new(inputs, WindConditions::default(), AtmosphericConditions::default());
    // `solve()` returning Err is a VALID outcome (input cleanly rejected).
    // A panic, or NaN/Inf/absurd values inside an Ok result, is the bug.
    if let Ok(result) = solver.solve() {
        domain::assert_finite_sane(&result);
    }
});
