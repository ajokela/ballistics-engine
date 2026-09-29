#!/usr/bin/env bash
# The release IS NOT DONE until this passes. Every check uses a cache-proof
# endpoint (registry front-page APIs lag; the sparse index and the JSON form of
# PyPI's /simple/ do not -- its HTML form DOES, see the PyPI check below).
set -euo pipefail
V="${1:?usage: verify-channels.sh VERSION}"
fail=0
chk() { local name="$1" got="$2" want="$3"
  if [ "$got" = "$want" ]; then echo "ok   $name: $got"; else echo "FAIL $name: got '$got' want '$want'"; fail=1; fi; }

chk "crates.io sparse index" "$(curl -s https://index.crates.io/ba/ll/ballistics-engine | tail -1 | python3 -c 'import json,sys;print(json.load(sys.stdin)["vers"])')" "$V"
# PyPI: the JSON form of /simple/ (PEP 691), which is what pip itself requests. The
# HTML form of the SAME URL is cached separately at the CDN edge, and for 0.44.0 an
# edge kept serving a page that ended at 0.43.0 for more than an hour after the
# upload -- while `pip install ballistics-engine==0.44.0` worked -- so this check
# failed a release that had shipped. Filenames are compared as fixed strings,
# anchored at both ends: a bare grep for the version once matched ANY sha256 hex
# digest (the dots are regex wildcards), passing for versions never published.
pypi_wheels() {
  curl -s -H 'Accept: application/vnd.pypi.simple.v1+json' https://pypi.org/simple/ballistics-engine/ \
    | python3 -c '
import json, sys
prefix = "ballistics_engine-" + sys.argv[1] + "-"
try:
    files = json.load(sys.stdin)["files"]
except (ValueError, KeyError):
    print("unreadable")  # e.g. the CDN answered HTML despite the Accept header
    sys.exit()
print("yes" if any(f["filename"].startswith(prefix) and f["filename"].endswith(".whl") for f in files) else "no")
' "$V"
}
chk "PyPI simple index (JSON, what pip reads) has wheels" "$(pypi_wheels)" "yes"
chk "RubyGems" "$(curl -s https://rubygems.org/api/v1/gems/ballistics-engine.json | python3 -c 'import json,sys;print(json.load(sys.stdin)["version"])')" "$V"
# npm was invisible here until MBA-1434, and it is exactly the channel that went missing:
# ten releases between 0.25.0 and 0.36.3 never reached npm and nothing said so. Reads the
# registry packument (registry.npmjs.org), not the npmjs.com package page, which lags it.
# `dist-tags.latest` rather than "the version exists somewhere": an out-of-order re-run
# publishes under the `backfill` dist-tag on purpose, so a version can be present on the
# registry without being the release. Only the release under test should be moving `latest`.
chk "npm" "$(curl -s https://registry.npmjs.org/ballistics-engine | python3 -c 'import json,sys;print(json.load(sys.stdin)["dist-tags"]["latest"])')" "$V"
chk "GH release assets (38 = 16 bins + 16 sha + 6 provenance)" "$(gh release view "v$V" --repo ajokela/ballistics-engine --json assets --jq '.assets|length')" "38"
chk "GCS objects" "$(gsutil ls "gs://ballistics-releases/$V/" 2>/dev/null | wc -l | tr -d ' ')" "38"
chk "live terminal badge" "$(curl -s https://ballistics.rs/sh/ | grep -o "Ballistics Engine v[0-9.]*" | head -1)" "Ballistics Engine v$V"
chk "docs" "$(curl -sL https://docs-ballistics-rs.web.app/generated-docs/ballistics_engine/index.html | grep -coF "$V" | awk '{print ($1>0)?"yes":"no"}')" "yes"
# wasm byte-parity: live vs local copy
LOCAL="$HOME/projects/ballistics.rs/ballistics-rs-site/sh/wasm/ballistics_engine_bg.wasm"
if [ -f "$LOCAL" ]; then
  chk "live wasm bytes == local" "$(curl -s https://ballistics.rs/sh/wasm/ballistics_engine_bg.wasm | shasum -a 256 | cut -d' ' -f1)" "$(shasum -a 256 "$LOCAL" | cut -d' ' -f1)"
fi
exit $fail
