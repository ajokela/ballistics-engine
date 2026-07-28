#!/usr/bin/env bash
# Assemble the GitHub release + GCS upload from a directory of built artifacts.
# Refuses to run unless all 13 binaries are present and every checksum verifies.
set -euo pipefail
V="${1:?usage: assemble.sh VERSION [DIR]}"; DIR="${2:-$HOME/release-$V}"
cd "$DIR"
EXPECTED=(macos-aarch64 macos-x86_64 linux-x86_64 linux-aarch64 linux-riscv64 linux-mips64el
          windows-x86_64.exe freebsd-x86_64 openbsd-x86_64 netbsd-x86_64
          freebsd-aarch64 openbsd-aarch64 netbsd-aarch64)
missing=0
for p in "${EXPECTED[@]}"; do
  [ -f "ballistics-$V-$p" ] || { echo "MISSING: ballistics-$V-$p"; missing=1; }
done
[ "$missing" = 0 ] || exit 1
# Every producer (hosted CI, BSD VMs, fleet scripts) emits the two-space
# "<hash>  <bare-filename>" form, so this verifies rather than erroring out on
# format drift. If this ever fails with "no properly formatted lines", a producer
# regressed - fix the producer, not this check.
shasum -c ./*.sha256
scripts_dir=$(git -C "$(dirname "$0")" rev-parse --show-toplevel 2>/dev/null || echo "$HOME/projects/ballistics-engine")
"$scripts_dir/scripts/release/extract-notes.sh" "$V" > "/tmp/relnotes-$V.md"
gh release create "v$V" --repo ajokela/ballistics-engine --title "v$V" \
  --notes-file "/tmp/relnotes-$V.md" ballistics-"$V"-* \
  || gh release upload "v$V" --repo ajokela/ballistics-engine --clobber ballistics-"$V"-*
gsutil -m cp ballistics-"$V"-* "gs://ballistics-releases/$V/"
echo "GH assets: $(gh release view "v$V" --repo ajokela/ballistics-engine --json assets --jq '.assets|length')"
echo "GCS objects: $(gsutil ls "gs://ballistics-releases/$V/" | wc -l | tr -d ' ')"
