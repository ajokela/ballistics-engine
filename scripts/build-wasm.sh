#!/usr/bin/env bash
# The single entry point for every WASM build (MBA: WASM command gating).
#
# WHY THIS EXISTS: the browser terminal's commands sit behind per-command cargo features
# (wasm-zero, wasm-monte-carlo, ...) so embedders can drop what they never call. Those
# features cannot live in `default` -- every wasm32 build must pass --no-default-features
# (the default pdf/online/cli/bridge/ffi set does not belong in, and partly does not compile
# for, the browser module), and cargo has no per-target default feature set, so a default-on
# gate would be stripped by the very flag every affected build already passes.
#
# That leaves one failure mode: forget `--features wasm-terminal` and you get a module that
# builds, deploys, and passes a version check while answering "Unknown command" to everything
# except `trajectory`. This script closes it from both ends:
#
#   1. The DEFAULT preset is `full`. Forgetting an argument gives you the complete terminal,
#      never a stripped one. Stripping is something you must ask for by name.
#   2. Every build is VERIFIED against the command set the preset promised, by inspecting the
#      emitted .wasm. A mismatch in either direction is a hard failure before anything ships.
#
# Usage:
#   scripts/build-wasm.sh                                  # full terminal, --target web
#   scripts/build-wasm.sh --preset slim                    # Calculator + trajectory only
#   scripts/build-wasm.sh --features wasm-zero,wasm-lead   # à la carte
#   scripts/build-wasm.sh --target nodejs --out-dir /tmp/x
#
# Any environment the build needs is passed through, e.g.
#   CARGO_PROFILE_RELEASE_OPT_LEVEL=z scripts/build-wasm.sh --out-dir .../ballistics-sh-site/wasm
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

# Every gateable command, as `feature:help-header`. The help header is what the verifier looks
# for: `show_help` emits each command's section from a chunk gated on that same feature, so the
# string is present in the binary if and only if the feature is on.
ALL_COMMANDS=(
  "wasm-zero:Zero Command:"
  "wasm-lead:Lead Command:"
  "wasm-monte-carlo:Monte Carlo Command:"
  "wasm-truing:True Velocity Command:"
  "wasm-truing:True Wind Command:"
  "wasm-estimate-bc:Estimate BC Command:"
  "wasm-bc-convert:BC Convert Command:"
  "wasm-drag-curve:Drag Curve Command:"
  "wasm-reticle:Reticle Command:"
  "wasm-powder:Powder Command:"
  "wasm-recoil:Recoil Command:"
  "wasm-power-factor:Power Factor Command:"
)

preset=""
features=""
target="web"
out_dir=""
extra=()

while [ $# -gt 0 ]; do
  case "$1" in
    --preset)   preset="${2:?--preset needs a value (full|slim)}"; shift 2 ;;
    --features) features="${2:?--features needs a comma-separated list}"; shift 2 ;;
    --target)   target="${2:?--target needs a value}"; shift 2 ;;
    --out-dir)  out_dir="${2:?--out-dir needs a path}"; shift 2 ;;
    --)         shift; extra=("$@"); break ;;
    -h|--help)  sed -n '2,30p' "${BASH_SOURCE[0]}"; exit 0 ;;
    *)          echo "error: unknown argument '$1' (see --help)" >&2; exit 2 ;;
  esac
done

if [ -n "$preset" ] && [ -n "$features" ]; then
  echo "error: pass --preset or --features, not both" >&2; exit 2
fi

# DEFAULT = full. See the header: a missing argument must never yield a stripped module.
if [ -z "$preset" ] && [ -z "$features" ]; then
  preset="full"
fi

# `expect_mode` is resolved from the PRESET, deliberately not from `$features`. If it were
# derived from the feature string, then losing that string -- the exact bug this script exists
# to catch -- would also lower the expectation, and the verifier would cheerfully agree that a
# stripped module is what `--preset full` asked for. `full` promises all twelve commands no
# matter how the feature list was computed, so a broken mapping fails loudly.
case "${preset:-custom}" in
  full)   features="wasm-terminal"; expect_mode="all" ;;
  slim)   features="";              expect_mode="none" ;;
  custom)                           expect_mode="from-features" ;;
  *)      echo "error: unknown preset '$preset' (expected full|slim)" >&2; exit 2 ;;
esac

# Resolve the command set this feature list is expected to produce, so the verifier below can
# compare against the artifact rather than against our intentions.
expected_headers=()
unexpected_headers=()
for entry in "${ALL_COMMANDS[@]}"; do
  feat="${entry%%:*}"; header="${entry#*:}"
  on=false
  case "$expect_mode" in
    all)  on=true ;;
    none) on=false ;;
    from-features)
      case ",$features," in
        *",wasm-terminal,"*) on=true ;;
        *",$feat,"*)         on=true ;;
      esac ;;
  esac
  if $on; then expected_headers+=("$header"); else unexpected_headers+=("$header"); fi
done

command -v wasm-pack >/dev/null 2>&1 || {
  echo "error: wasm-pack not found on PATH (https://rustwasm.github.io/wasm-pack/installer/)" >&2
  exit 1
}

built_dir="$out_dir"
if [ -z "$built_dir" ]; then built_dir="$repo_root/pkg"; fi

echo "==> wasm build: preset=${preset:-custom} target=$target features=${features:-<none>}"
# --no-default-features is mandatory (see header). The bare `--` is load-bearing: wasm-pack
# forwards only post-`--` arguments to cargo, so `--features` placed before it is consumed as
# an invalid wasm-pack flag.
pack_args=(build --target "$target" --no-default-features --out-dir "$built_dir")
cargo_args=()
[ -n "$features" ] && cargo_args+=(--features "$features")
[ ${#extra[@]} -gt 0 ] && cargo_args+=("${extra[@]}")
if [ ${#cargo_args[@]} -gt 0 ]; then
  wasm-pack "${pack_args[@]}" -- "${cargo_args[@]}"
else
  wasm-pack "${pack_args[@]}"
fi

# --- verification ---------------------------------------------------------------------
# Assert the artifact carries exactly the commands the preset promised. `show_help` emits each
# command's section from a chunk gated on that command's feature, so the presence of its header
# string in the .wasm is a direct readout of the cfg that was actually compiled -- not of what
# we meant to pass.
wasm_file="$built_dir/ballistics_engine_bg.wasm"
[ -f "$wasm_file" ] || { echo "error: no .wasm at $wasm_file" >&2; exit 1; }

# Written as explicit `if` blocks, not `grep && { ... }`: under `set -e` a failing `&&` list
# exits the script, which would turn "this command is correctly absent" into a silent abort.
fail=0
has () { grep -qa -- "$2" "$1"; }

if [ ${#expected_headers[@]} -gt 0 ]; then
  for header in "${expected_headers[@]}"; do
    if ! has "$wasm_file" "$header"; then
      echo "MISSING: '$header' (its feature should be ON)" >&2; fail=1
    fi
  done
fi
if [ ${#unexpected_headers[@]} -gt 0 ]; then
  for header in "${unexpected_headers[@]}"; do
    if has "$wasm_file" "$header"; then
      echo "PRESENT: '$header' (its feature should be OFF)" >&2; fail=1
    fi
  done
fi
# `trajectory` and the Calculator API are never gated -- if these ever vanish, the gating has
# eaten something it must not.
if ! has "$wasm_file" "Trajectory Command:"; then
  echo "MISSING: 'Trajectory Command:' (never gated)" >&2; fail=1
fi
# Scan every emitted binding file rather than one name: wasm-pack puts the bindings in
# ballistics_engine.js for the web/nodejs targets but in ballistics_engine_bg.js for
# --target bundler, and the .d.ts is generated for all of them.
# No pipe here on purpose: `grep -q` exits on the first match, so `cat | grep -q` gets
# SIGPIPE and, under `set -o pipefail`, reports failure on a SUCCESSFUL match.
if ! grep -qa -- "calculateTrajectory" "$built_dir"/*.js "$built_dir"/*.d.ts 2>/dev/null; then
  echo "MISSING: Calculator.calculateTrajectory in the JS bindings (never gated)" >&2; fail=1
fi

if [ "$fail" -ne 0 ]; then
  echo "error: the built module's command set does not match preset '${preset:-custom}'" >&2
  echo "       features passed: ${features:-<none>}" >&2
  exit 1
fi

size=$(wc -c < "$wasm_file" | tr -d ' ')
gz=$(gzip -9 -c "$wasm_file" | wc -c | tr -d ' ')
echo "==> verified: ${#expected_headers[@]} gateable command section(s) present, ${#unexpected_headers[@]} correctly absent"
echo "==> $wasm_file: $size bytes raw, $gz gzipped"
# Trimming changes the module's IMPORT list, not just its size: dropping wasm-monte-carlo
# removes the last user of `rand`, so the slim .wasm stops importing crypto.getRandomValues
# and the generated glue stops supplying it. A stale .wasm against fresh glue then fails to
# instantiate with a LinkError naming __wbg_getRandomValues. Say so on every build -- this
# already cost an embedder a debugging session.
echo "==> deploy ALL files in $built_dir together; the .wasm and its JS glue are a matched pair"
