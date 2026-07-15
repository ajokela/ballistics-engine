#!/usr/bin/env bash
# Build npm-publishable package(s) from the WASM (wasm-bindgen) build (MBA-1321).
#
# PREP ONLY: this script builds and post-processes the package(s); it does NOT publish.
# Publishing needs an npm scope decided by the crate owner and an npm account with rights to
# it — neither exists yet as of this writing. See this repo's README.md ("WASM / npm Package"
# section) for the manual publish command, and README-npm.md's placeholder scope name.
#
# Produces two directories at the repo root, both gitignored and rebuilt from scratch on every
# run (safe to delete and re-run):
#
#   pkg/      --target bundler   The package meant for `npm publish`. Consumed via a native
#             `.wasm` ES import by bundlers that understand it (webpack with
#             experiments.asyncWebAssembly, Vite, Rollup + @rollup/plugin-wasm, Parcel).
#
#   pkg-web/  --target web       No-bundler build for direct `<script type="module">` browser
#             use or manual Node usage without a bundler. This is the SAME --target already
#             used to build ballistics.sh/ballistics.rs's WASM (see this repo's CLAUDE.md).
#             Documented and built here for completeness, but not published under the primary
#             package name in this PREP pass — see README.md.
#
# Why two directories instead of one dual-target npm package: `wasm-pack build --help` (checked
# against 0.13.1) takes exactly one `--target` per invocation and has no "build both" or "dual
# publish" mode. Stitching a bundler build and a web build into one package.json via manual
# "exports" conditions is possible in principle, but wasm-bindgen's bundler- and web-target
# output use different module wiring (the bundler target imports the .wasm file directly and
# expects the bundler to instantiate it; the web target ships an explicit async `init()` that
# fetches it) and wasm-pack itself doesn't generate or test that wiring for you. A single
# bundler-target package as the npm-published artifact — with the web build documented
# separately for the no-bundler case — is the ecosystem-standard shape for wasm-bindgen crates
# published to npm, and is what this script does.
#
# Both builds use --no-default-features: the crate's default `pdf`/`online` features pull in
# printpdf / ureq+ring, which do not compile for wasm32-unknown-unknown (see CLAUDE.md's
# "Updating the WASM Module").
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"
cd "$repo_root"

bundler_dir="pkg"
web_dir="pkg-web"

command -v wasm-pack >/dev/null 2>&1 || {
  echo "error: wasm-pack not found on PATH (https://rustwasm.github.io/wasm-pack/installer/)" >&2
  exit 1
}
command -v node >/dev/null 2>&1 || {
  echo "error: node not found on PATH (needed to post-process package.json, and to run npm pack/publish)" >&2
  exit 1
}

echo "==> building bundler target -> ${bundler_dir}/"
wasm-pack build --target bundler --no-default-features --out-dir "$bundler_dir"

echo "==> building web target -> ${web_dir}/"
wasm-pack build --target web --no-default-features --out-dir "$web_dir"

echo "==> post-processing ${bundler_dir}/package.json"
node "$script_dir/build-npm-postprocess.mjs" "$repo_root/$bundler_dir/package.json" "@SCOPE/ballistics-engine"

echo "==> post-processing ${web_dir}/package.json"
node "$script_dir/build-npm-postprocess.mjs" "$repo_root/$web_dir/package.json" "@SCOPE/ballistics-engine-web"

echo "==> installing README-npm.md as the package README"
cp "$repo_root/README-npm.md" "$repo_root/$bundler_dir/README.md"
cp "$repo_root/README-npm.md" "$repo_root/$web_dir/README.md"

cat <<EOF

==> done.

  ${bundler_dir}/      (bundler target -- the package meant for "npm publish")
  ${web_dir}/  (web target -- no-bundler / direct-browser build; documented, not
             published under the primary package name)

TODO before publishing: ${bundler_dir}/package.json's "name" is the placeholder
"@SCOPE/ballistics-engine" -- replace SCOPE with the maintainer's real npm org/user
scope. See README.md's "WASM / npm Package" section for the publish command.

Verify the tarball contents before publishing:
  cd ${bundler_dir} && npm pack --dry-run
EOF
