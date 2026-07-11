#![no_main]
use libfuzzer_sys::fuzz_target;
use arbitrary::Unstructured;
use ballistics_engine::{AtmosphericConditions, TrajectorySolver, WindConditions};
use ballistics_engine_fuzz::domain;

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    // Feed FINITE, physically-plausible inputs (Finding #1: the engine lacks a
    // finite/physical input-validation gate and returns Ok(NaN) for non-finite
    // inputs — logged separately). Asserting on valid inputs makes this a
    // soundness fuzzer: any plausible input must yield finite-or-Err, never
    // Ok(NaN)/Inf/absurd. Widen back to domain::hostile_inputs once the engine
    // gains an input-validation gate.
    let Ok(inputs) = domain::valid_inputs(&mut u) else { return };
    let solver = TrajectorySolver::new(inputs, WindConditions::default(), AtmosphericConditions::default());
    // `solve()` returning Err is a VALID outcome (input cleanly rejected).
    // A panic, or NaN/Inf/absurd values inside an Ok result, is the bug.
    if let Ok(result) = solver.solve() {
        domain::assert_finite_sane(&result);
    }
});
