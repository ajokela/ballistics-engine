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
         solve_json_parser solve_v1_hostile)
for t in "${TARGETS[@]}"; do
  echo "==> building $t"
  cargo +nightly fuzz build --sanitizer none "$t"
done
for t in "${TARGETS[@]}"; do
  echo "==> smoke-running $t"
  cargo +nightly fuzz run --sanitizer none "$t" -- -runs=20000 -max_total_time=20 -timeout=5
done
echo "all fuzz targets built and smoke-ran clean"
