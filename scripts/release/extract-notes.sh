#!/usr/bin/env bash
# Extract one version's CHANGELOG section as release-notes markdown (stdout).
set -euo pipefail
V="${1:?usage: extract-notes.sh VERSION}"
# Resolve the repo from this script's own location, not the caller's cwd:
# assemble.sh invokes us after cd-ing into the artifact directory, which is not
# a git repository (bit the 0.33.1 release).
cd "$(git -C "$(cd "$(dirname "$0")" && pwd)" rev-parse --show-toplevel)"
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
