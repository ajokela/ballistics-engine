#!/usr/bin/env node
// Post-processes a wasm-pack-generated package.json for npm packaging (MBA-1321).
// Invoked by build-npm.sh; not meant to be run standalone (it takes no defaults),
// though it's safe to — it only rewrites the given package.json in place.
//
// wasm-pack already reads most of this metadata from Cargo.toml (description, license,
// repository, keywords, homepage), but this script sets it explicitly rather than trusting
// that indefinitely — a future wasm-pack version or an edit to Cargo.toml's [package] table
// could silently change or drop a field with no build failure to catch it. It also fixes two
// gaps in wasm-pack's own output:
//   - the package name is unscoped ("ballistics-engine"); real publishing needs an npm scope,
//     which the tooling can't know, so an unmissable placeholder is substituted instead (see
//     README.md's "WASM / npm Package" section for the manual step that replaces it).
//   - LICENSE-APACHE is copied into the output directory (wasm-pack copies every LICENSE*
//     file next to Cargo.toml) but is never added to package.json's "files" array, so
//     `npm pack`/`npm publish` silently drops it even though the crate is dual
//     MIT/Apache-2.0 licensed and both LICENSE files are sitting right there. Verified by
//     `npm pack --dry-run` before/after this fix during MBA-1321.
import { readFileSync, writeFileSync } from 'node:fs';

const [, , packageJsonPath, name] = process.argv;
if (!packageJsonPath || !name) {
  console.error('usage: build-npm-postprocess.mjs <path/to/package.json> <package-name>');
  process.exit(1);
}

const pkg = JSON.parse(readFileSync(packageJsonPath, 'utf8'));

// TODO(owner): "@SCOPE" is a placeholder. Replace it with the real npm org/user scope in
// pkg/package.json before publishing — see README.md's "WASM / npm Package" section.
pkg.name = name;
pkg.description =
  'High-performance ballistics trajectory engine with professional physics (WASM build)';
pkg.license = 'MIT OR Apache-2.0'; // must match Cargo.toml's [package].license
pkg.repository = { type: 'git', url: 'https://github.com/ajokela/ballistics-engine' };
pkg.homepage = 'https://ballistics.rs/';
pkg.keywords = ['ballistics', 'trajectory', 'physics', 'simulation', 'wasm', 'webassembly'];
// Scoped packages default to private on npm; --access public is also documented as an
// explicit manual publish flag in README.md, but setting this too means a plain
// `npm publish` (no flag) works the first time as well.
pkg.publishConfig = { access: 'public' };

pkg.files = Array.isArray(pkg.files) ? pkg.files : [];
if (!pkg.files.includes('LICENSE-APACHE')) {
  pkg.files.push('LICENSE-APACHE');
}

writeFileSync(packageJsonPath, JSON.stringify(pkg, null, 2) + '\n');
console.log(`  wrote ${packageJsonPath} (name=${name})`);
