#!/usr/bin/env bash
# Build the exact directory that is published to npm as `ballistics-engine` (MBA-1434).
#
# BUILD ONLY. This never publishes. `.github/workflows/publish-npm.yml` calls it and then runs
# `npm publish` with an OIDC trusted-publisher identity; run it by hand to see precisely what a
# release would ship, without shipping anything:
#
#   scripts/release/npm-package.sh 0.36.3 /tmp/npm-0.36.3   # prints the tarball's file list
#   npm pack --dry-run /tmp/npm-0.36.3                      # the same list with sizes + shasum
#
# WHICH BUILD IS THE PUBLISHED ONE
#
# `--target web`, not `--target bundler`. Every version on npm from 0.25.0 onward is the WEB
# build: it is the same artifact `deploy-wasm.sh` ships to ballistics.rs, and it was published
# by hand out of that script's /tmp output directory (the runbook said so in as many words).
# Confirmed against the registry rather than assumed -- ballistics-engine 0.36.3's package.json
# carries the web target's `sideEffects: ["./snippets/*"]` and its three-file `files` list, not
# the bundler target's `ballistics_engine_bg.js`.
#
# That distinction is load-bearing. The web build's entry point takes an explicit
# `await init()`; the bundler build's is instantiated by the bundler and does not. Publishing
# the bundler output under the same package name would break every existing consumer's import
# in what looks like a patch release, so this script pins the target that is already on the
# registry -- NOT the one `scripts/build-npm.sh` calls "the package meant for npm publish".
# That header was written in MBA-1321 before the package existed and was never reconciled with
# what the first real publish actually did.
#
# Everything else here is `scripts/build-wasm.sh` plus `scripts/build-npm-postprocess.mjs`,
# which is what makes this safe to automate: the wasm build verifies that the emitted module
# really carries all twelve gateable terminal commands (a forgotten `--features` otherwise
# ships a trajectory-only module that installs and imports fine), and the post-process fixes
# npm metadata wasm-pack gets wrong or omits.
set -euo pipefail

V="${1:?usage: npm-package.sh VERSION [OUT_DIR]}"
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT="${2:-$repo_root/pkg-npm}"
cd "$repo_root"

PKG_NAME="ballistics-engine"

# Build the wrong revision and every check below still passes -- they all read the artifact,
# and the artifact would faithfully describe whatever was compiled. So gate on the source
# first. (The workflow additionally ties VERSION to the git tag; here we only know Cargo.toml.)
CARGO_VERSION=$(grep -m1 '^version = ' Cargo.toml | cut -d'"' -f2)
[ "$V" = "$CARGO_VERSION" ] || {
  echo "error: asked for $V but Cargo.toml at this checkout says $CARGO_VERSION" >&2; exit 1; }

command -v node >/dev/null 2>&1 || { echo "error: node not on PATH" >&2; exit 1; }
# For the tarball listing at the end. `npm pack --dry-run` writes nothing and contacts no
# registry, so needing npm here does not make this script a publish path.
command -v npm >/dev/null 2>&1 || { echo "error: npm not on PATH" >&2; exit 1; }

# Absolute, because build-wasm.sh cds to the repo root and a relative --out-dir would land
# somewhere the caller did not name.
case "$OUT" in /*) ;; *) OUT="$PWD/$OUT" ;; esac
rm -rf "$OUT"

echo "==> building the wasm-pack web target -> $OUT"
./scripts/build-wasm.sh --preset full --target web --out-dir "$OUT"

echo "==> post-processing package.json"
node scripts/build-npm-postprocess.mjs "$OUT/package.json" "$PKG_NAME"

# --- verification --------------------------------------------------------------------------
# build-wasm.sh has already proved the MODULE is complete. What is left to prove is that the
# PACKAGE around it is the one we mean to publish: right name, right version, right entry
# points, and both licences actually inside the tarball.
echo "==> verifying the package"
fail=0
say_fail () { echo "FAIL: $*" >&2; fail=1; }

# A here-string, not `< <(...)`: node's output has no trailing newline, and `read` hitting EOF
# without one returns non-zero, which under `set -e` aborts the script right here -- after a
# clean build, before printing a single check, with no message and a bare exit 1. Observed
# while writing this. `<<<` appends the newline that keeps `read` happy.
read -r got_name got_version got_main got_types <<< "$(
  node -e '
    const p = JSON.parse(require("fs").readFileSync(process.argv[1], "utf8"));
    process.stdout.write([p.name, p.version, p.main, p.types].map(v => v ?? "<unset>").join(" "));
  ' "$OUT/package.json"
)"
[ "$got_name" = "$PKG_NAME" ]                 || say_fail "package name is '$got_name', want '$PKG_NAME'"
[ "$got_version" = "$V" ]                     || say_fail "package version is '$got_version', want '$V'"
# The web target's entry points. If wasm-pack ever emits `module` instead of `main`, or renames
# the glue, consumers' `import ... from "ballistics-engine"` resolves differently -- catch that
# here rather than in someone's build.
[ "$got_main" = "ballistics_engine.js" ]      || say_fail "main is '$got_main', want 'ballistics_engine.js'"
[ "$got_types" = "ballistics_engine.d.ts" ]   || say_fail "types is '$got_types', want 'ballistics_engine.d.ts'"

for f in ballistics_engine_bg.wasm ballistics_engine.js ballistics_engine.d.ts \
         package.json README.md LICENSE LICENSE-APACHE; do
  [ -f "$OUT/$f" ] || say_fail "missing $f in $OUT"
done

# npm auto-includes README, package.json and `LICENSE` whatever `files` says, but NOT
# `LICENSE-APACHE` -- checked against the published 0.36.3 tarball, which ships LICENSE and
# omits LICENSE-APACHE even though wasm-pack had copied both into the directory. The crate is
# dual MIT/Apache-2.0, so the post-process adds it to `files`; assert that it stuck. This one
# names the cause, which is why it is worth keeping alongside the tarball check below: if it
# fires, the post-process is what changed.
node -e '
  const p = JSON.parse(require("fs").readFileSync(process.argv[1], "utf8"));
  if (!Array.isArray(p.files) || !p.files.includes("LICENSE-APACHE")) {
    console.error("FAIL: package.json files[] does not list LICENSE-APACHE");
    process.exit(1);
  }
' "$OUT/package.json" || fail=1

# ...and then ask npm what it would ACTUALLY pack, rather than inferring it from `files[]`.
# Everything above reads the directory and package.json; only npm knows the result of its own
# auto-include rules, `files[]`, and any .npmignore. `--dry-run` builds no tarball and touches
# no registry. Seven entries is the whole package: six that npm already shipped in 0.36.3, plus
# the LICENSE-APACHE that has been silently dropped from every release so far.
echo "==> what npm would actually pack"
packed=$(npm pack --dry-run --json "$OUT" | node -e '
  let s = ""; process.stdin.on("data", d => s += d).on("end", () => {
    const [tarball] = JSON.parse(s);
    process.stdout.write(tarball.files.map(f => f.path).sort().join("\n"));
  });
')
printf '%s\n' "$packed" | sed 's/^/      /'
for f in LICENSE LICENSE-APACHE README.md package.json \
         ballistics_engine_bg.wasm ballistics_engine.js ballistics_engine.d.ts; do
  printf '%s\n' "$packed" | grep -qxF "$f" || say_fail "npm would not pack $f"
done

[ "$fail" -eq 0 ] || { echo "error: $OUT is not publishable as described above" >&2; exit 1; }

echo "==> ok: $OUT is $PKG_NAME@$V (wasm-pack web target)"
echo "    inspect the tarball with: npm pack --dry-run $OUT"
