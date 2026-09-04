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
# Overridable so the release can be cut from whichever checkout is actually AT the tag --
# an agent worktree, a detached checkout, a fresh clone -- rather than only from the one
# under $HOME. The tag gate below is what guarantees correctness; the path never did.
ENGINE="${ENGINE:-$HOME/projects/ballistics-engine}"; SITES="${SITES:-$HOME/projects/ballistics.rs}"
# grep -qx, not -q: an unanchored match accepts a NEIGHBOURING tag, so a request for
# 0.36.3 would pass on a checkout sitting at v0.36.30.
( cd "$ENGINE" && git describe --tags --exact-match 2>/dev/null | grep -qx "v$V" ) \
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

# Commit what was just published.
#
# This script deploys four copies of the wasm bundle and seds the terminal
# badge, then exited — so the sites repo was left dirty after every release and
# the drift was only ever noticed later. The badge in git sat at v0.22.9 for
# eleven releases while production served something newer, and the same gap
# reopened at 0.35.0.
#
# Staged by path, never `git add -A`: that repo accumulates unrelated
# work-in-progress and a release script must not sweep it up. Not pushed either
# — the deploy has already happened, so nothing is blocked by leaving the push
# to a human who has read the diff.
(
  cd "$SITES" || exit 0
  git rev-parse --git-dir >/dev/null 2>&1 || exit 0
  git add -- wasm ballistics-sh-site/wasm ballistics-rs-site/wasm \
              ballistics-rs-site/sh/wasm ballistics-rs-site/sh/index.html \
              .firebase 2>/dev/null || true
  if git diff --cached --quiet; then
    echo "sites repo: nothing to commit (already recorded)"
  else
    git commit -q -m "Publish engine $V to the browser terminals

Deployed by scripts/release/deploy-wasm.sh during the $V release: four copies of
the wasm bundle and the terminal badge. Recorded here so the repo matches what
ballistics.rs actually serves.

The .wasm and its JS glue are a matched pair and move together; a new .wasm
beside old glue fails at runtime, not at build time."
    echo "sites repo: committed the $V bundle."
    echo "            review with 'git -C $SITES show', then: git -C $SITES push origin main"
  fi
  if [ -n "$(git status --porcelain)" ]; then
    echo "sites repo: other uncommitted changes remain (left alone):"
    git status --short | sed "s/^/            /"
  fi
)
