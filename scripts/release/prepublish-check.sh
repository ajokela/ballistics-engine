#!/usr/bin/env bash
# Gate `cargo publish` on the two things that actually went wrong in 0.33.1.
#
# WHAT HAPPENED. 0.33.1 shipped seven Applied Ballistics geometry files —
# ground-truth transcriptions, sampling frames, holdout manifests — derived from
# AB's proprietary bullet-card drawings. They were never committed to git and are
# not referenced by any code path. They were sitting untracked in the working
# directory, and the publish went out with `--allow-dirty`.
#
# WHY A FILE-LIST CHECK ALONE WOULD NOT HAVE HELPED. Cargo already refuses this.
# Verified against v0.33.2 with an untracked file planted under data/validation:
#
#   $ cargo package --list
#   error: 1 files in the working directory contain changes that were not yet
#   committed into git: data/validation/ab_experiment_probe.json
#   to proceed despite this and include the uncommitted changes, pass the
#   `--allow-dirty` flag
#
# The untracked file did not appear in the listing — cargo bailed instead. So the
# packaging default was never the hole; overriding it was. The primary check here
# is therefore "the working tree is pristine", and the `exclude` list in
# Cargo.toml is a second line rather than the first.
#
# The manifest comparison still earns its place: `exclude` only guards paths
# somebody thought of, and a file that is COMMITTED but unwanted sails past both
# the dirty check and the excludes. Diffing the real package contents against a
# committed expectation catches that, and forces any legitimate change to the
# shipped file set to be a visible, deliberate edit.
#
# Usage:  scripts/release/prepublish-check.sh [expected-version]
# Run it immediately before `cargo publish --locked`. Exits non-zero on any
# finding; prints every finding rather than stopping at the first.

set -uo pipefail
cd "$(git rev-parse --show-toplevel)" || exit 2

MANIFEST="scripts/release/packaged-files.txt"
EXPECTED_VERSION="${1:-}"
fail=0
note() { printf '  %s\n' "$*"; }
bad()  { printf 'FAIL: %s\n' "$*"; fail=1; }

echo "== pre-publish check =="

# 1. Pristine working tree. This is the one that would have stopped 0.33.1.
#    --porcelain reports modified AND untracked, which is what we want: the AB
#    files were untracked, so a modified-only check would have missed them.
dirty=$(git status --porcelain)
if [ -n "$dirty" ]; then
  bad "working tree is not pristine — publishing now risks shipping local files"
  printf '%s\n' "$dirty" | sed 's/^/        /'
  note "if any of these are research or licensed data, they must not be packaged;"
  note "commit, remove, or ignore them. Do NOT reach for --allow-dirty."
else
  note "working tree pristine (no modified, no untracked)"
fi

# 2. HEAD should be a release tag. Publishing from an untagged commit is how you
#    ship something that no tag, release, or changelog entry describes.
tag=$(git describe --exact-match --tags HEAD 2>/dev/null)
if [ -z "$tag" ]; then
  bad "HEAD is not at a tag ($(git rev-parse --short HEAD)); publish from the release tag"
else
  note "HEAD is at tag $tag"
fi

# 3. Cargo.toml version agrees with the tag (and with the caller's expectation).
cargo_ver=$(sed -n '/^\[package\]/,/^\[/p' Cargo.toml | sed -n 's/^version *= *"\([^"]*\)".*/\1/p' | head -1)
note "Cargo.toml version: ${cargo_ver:-<unreadable>}"
if [ -n "$tag" ] && [ "${tag#v}" != "$cargo_ver" ]; then
  bad "tag $tag disagrees with Cargo.toml version $cargo_ver"
fi
if [ -n "$EXPECTED_VERSION" ] && [ "${EXPECTED_VERSION#v}" != "$cargo_ver" ]; then
  bad "requested version $EXPECTED_VERSION disagrees with Cargo.toml version $cargo_ver"
fi

# 4. Packaged file set matches the committed expectation.
#    Deliberately NOT passing --allow-dirty: if the tree is dirty this call fails,
#    which is the correct outcome and is already reported by check 1.
listing=$(cargo package --list --locked 2>/tmp/prepublish-list.err)
if [ -z "$listing" ]; then
  bad "cargo package --list produced nothing; see /tmp/prepublish-list.err"
  head -5 /tmp/prepublish-list.err 2>/dev/null | sed 's/^/        /'
elif [ ! -f "$MANIFEST" ]; then
  bad "$MANIFEST is missing; regenerate with: cargo package --list --locked | sort > $MANIFEST"
else
  added=$(comm -23 <(printf '%s\n' "$listing" | sort) <(sort "$MANIFEST"))
  removed=$(comm -13 <(printf '%s\n' "$listing" | sort) <(sort "$MANIFEST"))
  if [ -n "$added" ]; then
    bad "files would be packaged that are not in $MANIFEST:"
    printf '%s\n' "$added" | sed 's/^/        + /'
    note "if these are intended, add them to the manifest in the same commit."
  fi
  if [ -n "$removed" ]; then
    bad "files in $MANIFEST would NOT be packaged:"
    printf '%s\n' "$removed" | sed 's/^/        - /'
  fi
  [ -z "$added$removed" ] && \
    note "packaged file set matches $MANIFEST ($(printf '%s\n' "$listing" | wc -l | tr -d ' ') files)"
fi

echo
if [ "$fail" -ne 0 ]; then
  echo "pre-publish check FAILED — do not publish."
  exit 1
fi
echo "pre-publish check passed."
