#![no_main]

use arbitrary::Unstructured;
use ballistics_engine::solve_json::ResolvedWindV1;
use ballistics_engine::{decode_solve_request_v1, solve_v1, ResolvedSolveRequestV1, SolveRequestV1};
use ballistics_engine_fuzz::solve_json_v1::valid_request;
use libfuzzer_sys::fuzz_target;

/// Clear the two resolved-echo fields that `From<&ResolvedSolveRequestV1> for SolveRequestV1`
/// (`src/request_roundtrip.rs`) deliberately does NOT carry onto a rebuilt request, so a
/// second-generation `resolved_request` can be compared against the first without tripping on
/// an intentional difference.
///
/// Both fields are mode markers whose sibling value in the SAME struct is already the
/// post-conversion result: `atmosphere.pressure_pa` is already reduced from QNH to absolute
/// station pressure when `pressure_reference` was `"qnh"`, and every wind
/// `direction_from_rad` is already re-referenced from compass to shooter-relative when
/// `wind_reference` was `"compass"`. Echoing the mode back onto the rebuilt request alongside
/// an already-converted value would apply the conversion a second time on the next solve, so
/// `From` always emits the omitted-field default (absolute pressure, shooter-relative wind)
/// instead of the original mode. That makes `None` the CORRECT second-generation value even
/// when the first generation read `Some(..)` here -- see `request_roundtrip.rs`'s module doc
/// and its `roundtrip_does_not_double_reduce_a_qnh_pressure` /
/// `roundtrip_does_not_double_reference_a_compass_wind` tests, which pin exactly this
/// behavior. Do NOT "fix" this normalization by removing it: a whole-object comparison
/// without it fails on roughly half of all generated requests (the domain generator sets
/// `pressure_reference`/`wind_reference` to `Some(..)` about a quarter of the time each), not
/// because the round-trip is lossy, but because the exemption below is doing its job.
fn normalize(source: &ResolvedSolveRequestV1) -> ResolvedSolveRequestV1 {
    let mut normalized = source.clone();
    normalized.atmosphere.pressure_reference = None;
    match &mut normalized.wind {
        ResolvedWindV1::Constant(wind) => wind.wind_reference = None,
        ResolvedWindV1::Segmented(wind) => wind.wind_reference = None,
    }
    normalized
}

// Resolution must be idempotent through a round-trip for any request the service actually
// solves: resolve(to_request(resolve(r))) == resolve(r), modulo the two intentionally-dropped
// mode markers `normalize` strips above. A mismatch means `resolved_request` is lossy -- some
// input the first solve used never made it into the resolved echo -- so rebuilding a request
// from it and re-solving silently drops or changes that input. That is exactly the failure
// mode the 0.33.0 perturbation kernel cannot tolerate: it rebuilds a request from a resolved
// one to run a counterfactual, so a lossy round-trip would misattribute the counterfactual's
// effect to whatever the caller happened to be perturbing.
//
// A `resolved_request`-only comparison is not sufficient on its own: Task 2's review found a
// Critical bug (a silently dropped windage-zero convergence bias, see
// `roundtrip_preserves_the_windage_zero_bias` in `src/request_roundtrip.rs`) that left
// `resolved_request` byte-identical across the round-trip while the solved trajectory
// diverged. `summary` and `samples` are compared too so a lossy round-trip is caught even when
// it never touches a DTO field at all.
fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let Ok(request) = valid_request(&mut u) else {
        return;
    };
    // Go through the same public JSON boundary a real caller uses (encode, then the
    // structural + range-validating decoder) rather than handing the generated struct to
    // solve_v1 directly.
    let Ok(json) = serde_json::to_string(&request) else {
        return;
    };
    let Ok(decoded) = decode_solve_request_v1(&json) else {
        return;
    };
    let Ok(first) = solve_v1(decoded) else {
        return;
    };

    let rebuilt: SolveRequestV1 = (&first.resolved_request).into();
    let Ok(second) = solve_v1(rebuilt) else {
        panic!("a request rebuilt from a resolved request failed to solve");
    };

    assert_eq!(
        normalize(&first.resolved_request),
        normalize(&second.resolved_request),
        "resolved request changed after a round-trip (pressure_reference/wind_reference excluded -- see normalize())"
    );
    assert_eq!(
        first.summary, second.summary,
        "summary changed after a round-trip"
    );
    assert_eq!(
        first.samples, second.samples,
        "samples changed after a round-trip"
    );
});
