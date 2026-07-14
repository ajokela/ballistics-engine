//! MBA-1293: frozen robustness-inputs artifact, decoded exactly as the harness does.
//!
//! Before the fix, this input produced `Ok(max_range = -50.588...)`: every field passes
//! the positivity/finiteness gate, but the drag time constant is far below the minimum
//! integration step, so the accepted minimum-dt RK45 step multiplied speed 13x and
//! reversed it. The solve must now reject the divergence with `Err`.

use arbitrary::Unstructured;
use ballistics_engine::{AtmosphericConditions, TrajectorySolver, WindConditions};
use ballistics_engine_fuzz::domain;

#[test]
fn repro_mba_1293_is_rejected() {
    let data = include_bytes!("crash-mba-1293.bin");
    let mut u = Unstructured::new(data);
    let inputs = domain::hostile_inputs(&mut u).expect("inputs");
    let max_range = inputs.target_distance;
    let mut solver = TrajectorySolver::new(
        inputs,
        WindConditions::default(),
        AtmosphericConditions::default(),
    );
    solver.set_max_range(max_range);

    match solver.solve() {
        Err(_) => {} // clean rejection: the outcome the harness treats as valid
        Ok(result) => panic!(
            "crash-mba-1293.bin must be rejected, got Ok with max_range {}",
            result.max_range
        ),
    }
}
