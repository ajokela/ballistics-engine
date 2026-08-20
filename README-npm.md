# ballistics-engine (WASM)

WebAssembly build of [`ballistics-engine`](https://github.com/ajokela/ballistics-engine), a
high-performance ballistics trajectory engine (RK4 integration, wind/Coriolis/spin-drift/Magnus
effects, custom drag tables, Monte Carlo, and more). This package exposes the same CLI-style
command surface as the native Rust binary, plus a small object-oriented `Calculator` class, to
JavaScript/TypeScript.

Built with `wasm-pack` and `--no-default-features` — the native crate's default `pdf`/`online`
features pull in dependencies that don't compile for `wasm32-unknown-unknown`, so the PDF dope-card
export and the online BC-estimation API are not available from WASM.

If you only use the `Calculator` API, you can rebuild this module without the terminal commands
you never call and cut it roughly in half — see **Size** under [Caveats](#caveats).

> This package is built from `pkg/` (the `wasm-pack --target bundler` output) via
> `scripts/build-npm.sh` in the source repo. The same script also produces a `pkg-web/` build
> (`--target web`) for use without a bundler — see "Browser without a bundler" below.

## Install

```bash
npm install @SCOPE/ballistics-engine
```

`@SCOPE` is a placeholder — see the source repo's `README.md` ("WASM / npm Package" section) for
the real published name once one exists.

## Quick start (bundler: webpack, Vite, Rollup, Parcel)

This package's `main` entry imports its `.wasm` file as a native ES module, which is how
`wasm-pack --target bundler` output is meant to be consumed. It works out of the box with Vite and
Rollup (`@rollup/plugin-wasm`), and with webpack once `experiments.asyncWebAssembly` (or
`experiments.syncWebAssembly`) is enabled — check your bundler's WASM docs if the import fails.

```js
import { WasmBallistics } from '@SCOPE/ballistics-engine';

const calc = new WasmBallistics();

// The command surface mirrors the native CLI (see CLI_USAGE.md in the source repo for the full
// flag reference): .308 Winchester, 168gr @ 2700 fps, zeroed at 200 yd, table out to 500 yd.
const table = calc.runCommand(
  'trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 500 --auto-zero 200',
);
console.log(table);
```

### Custom drag tables (`loadDragTable`)

Supply a measured or manufacturer-published Mach:Cd drag curve (Hornady CDM data, a Lapua/Doppler
deck, or your own) instead of a G1/G7 model + BC. Once loaded it's applied automatically to every
`trajectory`, `zero`, `lead`, and `monte-carlo` run — no extra flag needed.

```js
const csv = 'mach,cd\n0.5,0.220\n0.8,0.230\n1.0,0.520\n1.2,0.480\n1.5,0.400\n2.0,0.330\n2.5,0.300\n';
calc.loadDragTable(new TextEncoder().encode(csv));
calc.hasDragTable(); // true

calc.runCommand('trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 500');
```

`loadDragTable` takes raw bytes because WASM has no filesystem access — fetch the CSV yourself
(`fetch()` in the browser, `fs.readFileSync` in Node) and pass the bytes in. The CSV format is
documented in the source repo's `CLI_USAGE.md` ("Custom Drag Tables"); a matching
`loadBc5dTable(bytes)` / `hasBc5dTable()` pair exists for BC5D correction tables.

## Browser without a bundler

For a plain `<script type="module">` page (no build step), use the `pkg-web/` build instead —
produced by the same `scripts/build-npm.sh`, and the same build already deployed at
[ballistics.sh](https://ballistics.sh). It ships an explicit async `init()` you call once before
constructing `WasmBallistics`:

```html
<script type="module">
  import init, { WasmBallistics } from './ballistics_engine.js';

  await init(); // fetches and instantiates ballistics_engine_bg.wasm relative to this file
  const calc = new WasmBallistics();
  console.log(calc.runCommand('trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 500'));
</script>
```

Serve `ballistics_engine.js` and `ballistics_engine_bg.wasm` from the same directory, and make sure
your host serves `.wasm` with `Content-Type: application/wasm` (all major static hosts and CDNs do
this by default).

This `pkg-web/` build is not published to npm as part of this package in the current release —
it's built locally alongside `pkg/` for direct use or self-hosting. If you need it from npm, either
vendor the files from `pkg-web/` yourself or publish it as a second package.

### Node.js without a bundler

Plain Node `import`/`require` cannot load this package's bundler-target `.wasm` import directly.
Use the `pkg-web/` build instead, passing the file bytes explicitly (Node's `fetch()` does not
support `file://` URLs):

```js
import { readFileSync } from 'node:fs';
import init, { WasmBallistics } from './pkg-web/ballistics_engine.js';

const wasmBytes = readFileSync(new URL('./pkg-web/ballistics_engine_bg.wasm', import.meta.url));
await init({ module_or_path: wasmBytes });

const calc = new WasmBallistics();
console.log(calc.runCommand('trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 500'));
```

If you need a CommonJS (`require()`) Node build, generate one yourself:
`scripts/build-wasm.sh --target nodejs` (add `--preset slim` for a trajectory-only module —
see **Size** under Caveats).

## Caveats

- **Size**: this package ships the full command surface — 918 KB raw, 345 KB gzipped, 274 KB
  brotli (measured on 0.33.2; decimal KB). It is not code-split or lazily loaded, so the whole engine loads
  up front.

  **Most of that is the terminal, and you can drop it.** Every command except `trajectory` is
  behind its own cargo feature, while the `Calculator` API is never gated. Rebuilding with only
  what you call gets a trajectory-only module down to **483 KB raw / 191 KB gzipped** — 44% off
  the wire — with byte-identical trajectory output:

  ```bash
  # Calculator + trajectory only
  scripts/build-wasm.sh --preset slim --target bundler

  # ...or keep a subset
  scripts/build-wasm.sh --features wasm-zero,wasm-lead --target bundler
  ```

  The script verifies the module it produced actually matches the set you asked for.

  The per-command features and their measured savings are tabulated under
  "Trimming the WASM module" in the [source repo's
  README](https://github.com/ajokela/ballistics-engine#trimming-the-wasm-module). The biggest
  single win is `wasm-monte-carlo` (46 KB gzipped, and it takes the `--wez` sweep with it).
- **Single-threaded**: no SIMD/threads assumptions; no `SharedArrayBuffer` or
  cross-origin-isolation (COOP/COEP) headers required.
- **No filesystem/network**: table loaders (`loadDragTable`, `loadBc5dTable`) and any file-based
  CLI flags (e.g. native `--drag-table <FILE>`) need the host to fetch bytes and hand them in; see
  `loadDragTable` above.
- **`pdf`/`online` features are unavailable**: this build excludes them (see above), so
  PDF dope-card export and the online BC-estimation API are not part of the WASM surface.
- Full API surface (including the `Calculator` builder class) is documented in the bundled
  `.d.ts`; the full CLI flag reference `runCommand` accepts is documented in the source repo's
  `CLI_USAGE.md`.

## License

MIT OR Apache-2.0 — see `LICENSE` and `LICENSE-APACHE` in this package.

## Source

https://github.com/ajokela/ballistics-engine
