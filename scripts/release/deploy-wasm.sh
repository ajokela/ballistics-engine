#!/usr/bin/env bash
# WASM x4 + both firebase sites + the hardcoded badge. Live sites are in firebase
# project ballistics-rs (NOT ballistics-buddy-lfx33 - the engine CLAUDE.md is stale
# on this). --no-default-features is REQUIRED (pdf/online don't build for wasm32).
# --features wasm-terminal is REQUIRED too: the browser terminal's non-trajectory
# commands (zero, monte-carlo, reticle, ...) are individually gated so embedders can
# drop them, and those gates are NOT in `default` -- a default-on gate would be
# stripped by the --no-default-features above, i.e. by the very build it configures.
# Without it the live terminal silently degrades to trajectory + version only.
# Note the `--`: wasm-pack forwards only post-`--` args to cargo, so `--features`
# before it is consumed as an (invalid) wasm-pack flag and the build fails.
set -euo pipefail
V="${1:?usage: deploy-wasm.sh VERSION}"
ENGINE="$HOME/projects/ballistics-engine"; SITES="$HOME/projects/ballistics.rs"
( cd "$ENGINE" && git describe --tags --exact-match 2>/dev/null | grep -q "v$V" ) \
  || { echo "engine checkout is not at tag v$V - refuse to build wasm from the wrong rev"; exit 1; }
( cd "$ENGINE" && wasm-pack build --target web --no-default-features --out-dir "/tmp/wasm-$V" \
    -- --features wasm-terminal )
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
