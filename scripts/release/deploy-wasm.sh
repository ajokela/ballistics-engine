#!/usr/bin/env bash
# WASM x4 + both firebase sites + the hardcoded badge. Live sites are in firebase
# project ballistics-rs (NOT ballistics-buddy-lfx33 - the engine CLAUDE.md is stale
# on this). --no-default-features is REQUIRED (pdf/online don't build for wasm32).
# The build goes through scripts/build-wasm.sh, which defaults to the FULL terminal and
# then VERIFIES the emitted .wasm actually carries all twelve gateable commands. The
# browser terminal's non-trajectory commands are individually gated so embedders can drop
# them, and those gates cannot live in `default` (a default-on gate would be stripped by
# the --no-default-features every wasm build must pass). Building wasm-pack directly here
# would put a silent trajectory-only deploy one forgotten flag away.
set -euo pipefail
V="${1:?usage: deploy-wasm.sh VERSION}"
ENGINE="$HOME/projects/ballistics-engine"; SITES="$HOME/projects/ballistics.rs"
( cd "$ENGINE" && git describe --tags --exact-match 2>/dev/null | grep -q "v$V" ) \
  || { echo "engine checkout is not at tag v$V - refuse to build wasm from the wrong rev"; exit 1; }
( cd "$ENGINE" && ./scripts/build-wasm.sh --preset full --target web --out-dir "/tmp/wasm-$V" )
[ "$(python3 -c "import json;print(json.load(open('/tmp/wasm-$V/package.json'))['version'])")" = "$V" ] || { echo "pkg version mismatch"; exit 1; }
for d in wasm ballistics-sh-site/wasm ballistics-rs-site/wasm ballistics-rs-site/sh/wasm; do
  cp "/tmp/wasm-$V"/ballistics_engine_bg.wasm "/tmp/wasm-$V"/ballistics_engine_bg.wasm.d.ts \
     "/tmp/wasm-$V"/ballistics_engine.d.ts "/tmp/wasm-$V"/ballistics_engine.js \
     "/tmp/wasm-$V"/package.json "$SITES/$d/"
done
sed -i '' -E "s/Ballistics Engine v[0-9]+\.[0-9]+\.[0-9]+/Ballistics Engine v$V/" "$SITES/ballistics-rs-site/sh/index.html"
( cd "$SITES" && firebase deploy --only hosting --config firebase.ballistics-rs.json --project ballistics-rs )
sleep 5
LIVE=$(curl -s https://ballistics.rs/sh/ | grep -o "Ballistics Engine v[0-9.]*" | head -1)
[ "$LIVE" = "Ballistics Engine v$V" ] || { echo "LIVE BADGE MISMATCH: $LIVE"; exit 1; }
echo "OK: live terminal serves v$V"
