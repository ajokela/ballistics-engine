#!/usr/bin/env bash
# Bump the engine to a new version: Cargo.toml, BOTH lockfiles, and the CHANGELOG cut.
#
# The fuzz lockfile is the step everyone forgets: it pins TWO copies of
# ballistics-engine (the local path crate AND the 0.21.5 reference engine that
# differential_prev links). Only the local one may move — `cargo update -p
# ballistics-engine@<old>` disambiguates. Forgetting this fails CI's --locked
# fuzz steps, which is exactly how dependabot's printpdf PR could never go green.
set -euo pipefail
NEW="${1:?usage: bump.sh NEW_VERSION (e.g. 0.30.1)}"
cd "$(git rev-parse --show-toplevel)"

OLD=$(grep -m1 '^version = ' Cargo.toml | cut -d'"' -f2)
[ "$OLD" != "$NEW" ] || { echo "already at $NEW"; exit 1; }
DATE=$(date +%Y-%m-%d)

sed -i '' "s/^version = \"$OLD\"/version = \"$NEW\"/" Cargo.toml
cargo update -p "ballistics-engine@$OLD" --precise "$NEW" 2>/dev/null || cargo update -p ballistics-engine
cargo update --manifest-path fuzz/Cargo.toml -p "ballistics-engine@$OLD" --precise "$NEW" 2>/dev/null \
  || cargo update --manifest-path fuzz/Cargo.toml -p "ballistics-engine@$OLD"

python3 - "$NEW" "$DATE" <<'PY'
import sys
new, date = sys.argv[1], sys.argv[2]
s = open('CHANGELOG.md').read()
if '## [Unreleased]' not in s:
    sys.exit("CHANGELOG has no [Unreleased] section - nothing to cut")
s = s.replace('## [Unreleased]', f'## [{new}] - {date}', 1)
open('CHANGELOG.md','w').write(s)
PY

echo "== pins =="
grep -A1 '^name = "ballistics-engine"$' Cargo.lock fuzz/Cargo.lock | grep version
echo "Bumped $OLD -> $NEW. Now: commit, push, WAIT FOR CI GREEN, then: git tag v$NEW <sha> && git push origin v$NEW"
