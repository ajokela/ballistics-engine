#![no_main]

use arbitrary::Unstructured;
use ballistics_engine::solve_json::SolveErrorEnvelopeV1;
use ballistics_engine::solve_v1;
use ballistics_engine_fuzz::solve_json_v1::{assert_finite_success, bounded_hostile_request};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let Ok(request) = bounded_hostile_request(&mut u) else {
        return;
    };

    // Do not catch panics: a panic, libFuzzer timeout, or OOM is a finding. The generator bounds
    // integration and allocation multipliers so those signals represent service defects rather
    // than an intentionally unbounded request.
    match solve_v1(request) {
        Ok(success) => assert_finite_success(&success),
        Err(error) => {
            let encoded = serde_json::to_vec(&error).expect("service error must serialize");
            let decoded: SolveErrorEnvelopeV1 =
                serde_json::from_slice(&encoded).expect("serialized service error must decode");
            assert_eq!(decoded, error, "service error changed on JSON round trip");
        }
    }
});
