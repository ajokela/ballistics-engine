#!/usr/bin/env bash
# Build every fuzz target and run each for a few seconds. Used as a local/CI gate.
set -euo pipefail
# macOS aarch64 needs the deferred-symbol link workaround for the engine's cdylib under
# sancov-without-sanitizer (see the plan's Global Constraints). Not needed on Linux, where
# shared-object links tolerate the deferred symbols.
if [[ "$(uname)" == "Darwin" ]]; then
  export RUSTFLAGS="${RUSTFLAGS:-} -Clink-arg=-Wl,-undefined,dynamic_lookup"
fi
cd "$(dirname "$0")/../fuzz"
TARGETS=(robustness_inputs robustness_ffi robustness_montecarlo invariant_monotonic \
         invariant_symmetry differential_prev analytic_vacuum solver_zero \
         solve_json_parser solve_v1_hostile roundtrip_resolved)
# robustness_ffi fuzzes the current engine's C ABI, which lives behind the `ffi` feature. Every
# other target must be built WITHOUT it: differential_prev also links the 0.21.5 reference engine,
# whose ffi #[no_mangle] symbols would otherwise collide with the current engine's copy.
feature_args() { if [[ "$1" == "robustness_ffi" ]]; then echo "--features ffi"; fi; }
for t in "${TARGETS[@]}"; do
  echo "==> building $t"
  # shellcheck disable=SC2046 # intentional word-splitting of the optional --features flag
  cargo +nightly fuzz build --sanitizer none $(feature_args "$t") "$t"
done
for t in "${TARGETS[@]}"; do
  echo "==> smoke-running $t"
  # shellcheck disable=SC2046 # intentional word-splitting of the optional --features flag
  cargo +nightly fuzz run --sanitizer none $(feature_args "$t") "$t" -- -runs=20000 -max_total_time=20 -timeout=5
done
echo "all fuzz targets built and smoke-ran clean"
