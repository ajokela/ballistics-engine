#!/usr/bin/env bash
# Extract one version's CHANGELOG section as release-notes markdown (stdout).
set -euo pipefail
V="${1:?usage: extract-notes.sh VERSION}"
cd "$(git rev-parse --show-toplevel)"
python3 - "$V" <<'PY'
import re, sys
v = sys.argv[1]
s = open('CHANGELOG.md').read()
m = re.search(rf'^## \[{re.escape(v)}\][^\n]*\n(.*?)(?=^## \[|\Z)', s, re.M | re.S)
if not m:
    sys.exit(f"no CHANGELOG section for {v}")
print(m.group(1).strip())
print("\n---\nFull changelog: https://github.com/ajokela/ballistics-engine/blob/main/CHANGELOG.md")
PY
