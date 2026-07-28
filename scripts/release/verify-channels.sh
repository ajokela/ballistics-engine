#!/usr/bin/env bash
# The release IS NOT DONE until this passes. Every check uses a cache-proof
# endpoint (registry front-page APIs lag; the sparse index and /simple/ do not).
set -euo pipefail
V="${1:?usage: verify-channels.sh VERSION}"
fail=0
chk() { local name="$1" got="$2" want="$3"
  if [ "$got" = "$want" ]; then echo "ok   $name: $got"; else echo "FAIL $name: got '$got' want '$want'"; fail=1; fi; }

chk "crates.io sparse index" "$(curl -s https://index.crates.io/ba/ll/ballistics-engine | tail -1 | python3 -c 'import json,sys;print(json.load(sys.stdin)["vers"])')" "$V"
chk "PyPI /simple/ has wheels" "$(curl -s https://pypi.org/simple/ballistics-engine/ | grep -c "$V" | awk '{print ($1>0)?"yes":"no"}')" "yes"
chk "RubyGems" "$(curl -s https://rubygems.org/api/v1/gems/ballistics-engine.json | python3 -c 'import json,sys;print(json.load(sys.stdin)["version"])')" "$V"
chk "GH release assets (29 = 13 bins + 13 sha + 3 provenance)" "$(gh release view "v$V" --repo ajokela/ballistics-engine --json assets --jq '.assets|length')" "29"
chk "GCS objects" "$(gsutil ls "gs://ballistics-releases/$V/" 2>/dev/null | wc -l | tr -d ' ')" "29"
chk "live terminal badge" "$(curl -s https://ballistics.rs/sh/ | grep -o "Ballistics Engine v[0-9.]*" | head -1)" "Ballistics Engine v$V"
chk "docs" "$(curl -sL https://docs-ballistics-rs.web.app/generated-docs/ballistics_engine/index.html | grep -oc "$V" | awk '{print ($1>0)?"yes":"no"}')" "yes"
# wasm byte-parity: live vs local copy
LOCAL="$HOME/projects/ballistics.rs/ballistics-rs-site/sh/wasm/ballistics_engine_bg.wasm"
if [ -f "$LOCAL" ]; then
  chk "live wasm bytes == local" "$(curl -s https://ballistics.rs/sh/wasm/ballistics_engine_bg.wasm | shasum -a 256 | cut -d' ' -f1)" "$(shasum -a 256 "$LOCAL" | cut -d' ' -f1)"
fi
exit $fail
