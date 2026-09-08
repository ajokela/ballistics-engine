#!/usr/bin/env node
// Post-processes a wasm-pack-generated package.json for npm packaging (MBA-1321).
// Invoked by build-npm.sh; not meant to be run standalone (it takes no defaults),
// though it's safe to — it only rewrites the given package.json in place.
//
// wasm-pack already reads most of this metadata from Cargo.toml (description, license,
// repository, keywords, homepage), but this script sets it explicitly rather than trusting
// that indefinitely — a future wasm-pack version or an edit to Cargo.toml's [package] table
// could silently change or drop a field with no build failure to catch it. It also fixes one
// gap in wasm-pack's own output:
//   - LICENSE-APACHE is copied into the output directory (wasm-pack copies every LICENSE*
//     file next to Cargo.toml) but is never added to package.json's "files" array, so
//     `npm pack`/`npm publish` silently drops it even though the crate is dual
//     MIT/Apache-2.0 licensed and both LICENSE files are sitting right there. Verified by
//     `npm pack --dry-run` before/after this fix during MBA-1321, and again against the
//     published 0.36.3 tarball, which ships LICENSE (npm always includes that name) and
//     omits LICENSE-APACHE (npm does not).
import { readFileSync, writeFileSync } from 'node:fs';

const [, , packageJsonPath, name] = process.argv;
if (!packageJsonPath || !name) {
  console.error('usage: build-npm-postprocess.mjs <path/to/package.json> <package-name>');
  process.exit(1);
}

const pkg = JSON.parse(readFileSync(packageJsonPath, 'utf8'));

pkg.name = name;
pkg.description =
  'High-performance ballistics trajectory engine with professional physics (WASM build)';
pkg.license = 'MIT OR Apache-2.0'; // must match Cargo.toml's [package].license
// Canonical `git+https://...git` form, which is what npm normalizes the shorter spellings
// to anyway. Written out in full because provenance attestation matches this field against the
// repository the publish actually ran in, and a mismatch fails the publish rather than quietly
// dropping the attestation.
pkg.repository = { type: 'git', url: 'git+https://github.com/ajokela/ballistics-engine.git' };
pkg.homepage = 'https://ballistics.rs/';
pkg.keywords = ['ballistics', 'trajectory', 'physics', 'simulation', 'wasm', 'webassembly'];
// A no-op for the unscoped `ballistics-engine` name that is actually published, and the whole
// story if this output is ever republished under a scope — scoped packages default to private,
// so a plain `npm publish` fails without it. Kept as the belt to `--access public`'s braces.
pkg.publishConfig = { access: 'public' };

pkg.files = Array.isArray(pkg.files) ? pkg.files : [];
if (!pkg.files.includes('LICENSE-APACHE')) {
  pkg.files.push('LICENSE-APACHE');
}

writeFileSync(packageJsonPath, JSON.stringify(pkg, null, 2) + '\n');
console.log(`  wrote ${packageJsonPath} (name=${name})`);
