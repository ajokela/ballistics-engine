#![no_main]

use ballistics_engine::solve_json::{decode_solve_request_v1, SolveErrorEnvelopeV1};
use libfuzzer_sys::fuzz_target;

const MAX_PROCESS_INPUT_BYTES: usize = 1024 * 1024;

fuzz_target!(|data: &[u8]| {
    // The process transport rejects anything larger than this. Keeping the same explicit ceiling
    // prevents a corpus file from turning parser fuzzing into an allocation benchmark.
    if data.len() > MAX_PROCESS_INPUT_BYTES {
        return;
    }

    let Ok(text) = std::str::from_utf8(data) else {
        // Invalid UTF-8 is a clean parser rejection at the process boundary. The target still sees
        // every byte string; only UTF-8 strings can reach the public string-based DTO decoder.
        return;
    };

    match decode_solve_request_v1(text) {
        Ok(request) => {
            let encoded = serde_json::to_vec(&request).expect("decoded request must serialize");
            let encoded = std::str::from_utf8(&encoded).expect("JSON serialization is UTF-8");
            let decoded =
                decode_solve_request_v1(encoded).expect("serialized request must decode again");
            assert_eq!(decoded, request, "request changed on JSON round trip");
        }
        Err(error) => {
            let encoded = serde_json::to_vec(&error).expect("parser error must serialize");
            let decoded: SolveErrorEnvelopeV1 =
                serde_json::from_slice(&encoded).expect("serialized parser error must decode");
            assert_eq!(decoded, error, "parser error changed on JSON round trip");
        }
    }
});
