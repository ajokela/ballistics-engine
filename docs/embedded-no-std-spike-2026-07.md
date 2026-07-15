# MBA-1324: `no_std` + `alloc` Feasibility Spike (Embedded Targets)

**Date:** 2026-07-15
**Type:** Timeboxed feasibility spike (investigation only — no `src/` changes)
**Scope:** Could the ballistics-engine solver *core* run on embedded, device-class hardware
(rangefinders, dope computers) under `no_std` + `alloc`?

## Framing: this is NOT about the existing WASM deployment

A field tester already ships our WASM build in an app **served from** an ESP32 (the ESP32 is a
tiny HTTP server handing the page + `.wasm` to a phone or laptop that connects to it). That is a
completely different question from the one this spike answers:

- **What already works (WASM-served-from-ESP32):** the WASM binary executes in the *connecting
  browser's* JS engine — a phone or laptop with a full V8/SpiderMonkey/JSC runtime, gigabytes of
  RAM, and hardware double-precision floats. The ESP32 only serves bytes over Wi-Fi; it never
  runs a single instruction of our solver.
- **What this spike investigates:** compiling the solver itself to run natively on the
  microcontroller's own CPU — no browser, no OS, no filesystem, a few hundred KB of RAM, and
  (on the Cortex-M4-class part evaluated here) a single-precision-only FPU. This is the
  `no_std` + `alloc` question, and it is unrelated to whether WASM currently runs anywhere.

Nothing about the existing ESP32-hosting setup de-risks or informs this spike; they are
orthogonal deployment models that happen to share the word "ESP32."

## Executive summary

The dependency chain and the engine's own code are both **more no_std-friendly than expected,
though less cleanly isolated than an earlier draft of this report claimed** — corrected below and
in §2/§2a after a review pass caught the overclaim. Every core numeric dependency (nalgebra,
ndarray, rand, rand_distr, serde, serde_json, thiserror) already ships an `alloc`-based no_std
configuration, and in an *isolated probe crate* that explicitly sets `default-features = false`
on each of them, all seven compile cleanly together under `#![no_std]` + `extern crate alloc` on
`thumbv7em-none-eabihf` (a real Cortex-M4F target) — see §2a. **That isolation does not currently
hold inside `ballistics-engine` itself.** This crate's own `Cargo.toml` never sets
`default-features = false` on nalgebra/ndarray/serde (or rand/thiserror), so
`cargo check --target thumbv7em-none-eabihf --no-default-features` run against the real crate
fails today on *two* independent chains, not one: the un-gated `clap` chain, **and** the un-gated
numeric/serialization chain — `matrixmultiply` (pulled in by both nalgebra and ndarray) and
`serde_core` (pulled in by serde) each independently emit `error[E0463]: can't find crate for
std` in that exact run. A third, separate chain also fails: the unconditional `ndarray-npy`
dependency pulls in `zip` → `flate2`/`crc32fast` → `zlib-rs`, which likewise can't find `std`.
See §2/§2a for the corrected transcript and root cause. Grepping the actual solve path
(`cli_api.rs`, `trajectory_integration.rs`, `drag.rs`, `atmosphere.rs`, `wind.rs`, `wind_shear.rs`,
`moving_target.rs`, and the rest of the physics modules — spin drift, precession/nutation,
pitch damping, stability, transonic drag, form factor, Reynolds, aerodynamic jump) turned up
**zero** `std::fs`/`std::io`/`std::thread` usage and only a small, enumerable list of std-only
touchpoints (synchronization primitives — corrected count in §3, more than the "two" an earlier
draft claimed — one trait impl, a handful of `eprintln!` warning sites, and two `HashMap`s). None
of the CPU-heavy numerical hot path is std-coupled.

**Recommendation: qualified GO** — getting `cargo check --target thumbv7em-none-eabihf` green
for a "core" feature slice looks like roughly a one-week mechanical effort, not a rewrite. But
this spike only answers the *compiles* question. It does **not** answer whether the result is
fast/small/robust enough to be useful on real rangefinder-class hardware — see "Open questions"
and "Phase 3" below. Do not commit to shipping an embedded target until that hardware validation
happens.

---

## 1. Dependency audit

Investigated via `cargo tree`, reading each dependency's `Cargo.toml` `[features]` section
directly from the local registry checkout, and (for the "big six" numeric/serialization crates)
a standalone probe crate compiled against the real target — see §2.

| Crate | Version (locked) | no_std status | Notes |
|---|---|---|---|
| **nalgebra** | 0.35.0 | ✅ no_std-capable | `alloc` + `libm` features exist explicitly (`libm = ["simba/libm", ...]`); default features are `["std", "macros"]`, so `default-features = false` is required. **Compile-verified** on thumbv7em. |
| **ndarray** | 0.17.2 | ✅ no_std-capable | `#![cfg_attr(not(feature = "std"), no_std)]` in its own `lib.rs`; default is `["std"]` only. **Compile-verified.** Only used in `src/drag.rs`, and only by the peripheral filesystem-based G6/G8 table loader (see §3) — not on the hot path. |
| **rand** | 0.10.2 | ✅ no_std-capable (partial) | `alloc` feature exists; `std_rng` (ChaCha) does **not** require `std`. Default features (`std`, `std_rng`, `sys_rng`, `thread_rng`) do require std for OS entropy / TLS. A no_std build keeps `std_rng` for the algorithm but loses `sys_rng`/`thread_rng` — must be seeded manually or via a platform `getrandom` backend. **Compile-verified.** |
| **rand_distr** | 0.6.0 | ✅ no_std-capable | `alloc` feature exists (default is `std` = `alloc` + `rand/std`). **Compile-verified.** Used only in `cli_api.rs` for Monte Carlo. |
| **serde** (+derive) | 1.0.228 | ✅ no_std-capable | `alloc` feature exists (`alloc = ["serde_core/alloc"]`); `derive` is independent of `std`/`alloc`. **Compile-verified.** |
| **serde_json** | 1.0.149 | ✅ no_std-capable | Ships an `alloc` feature. **Compile-verified.** |
| **thiserror** | 2.0.18 | ✅ no_std by default | `std` is an opt-in feature (`default = ["std"]`, `std = []`); with it off, the derive targets `core::error::Error`, stable since Rust 1.81 (we're on 1.91.1). **Compile-verified.** Used only in `trajectory_observation.rs`. |
| **clap** / **clap_complete** | 4.5.60 / 4.5.66 | ❌ std-only | No no_std support (terminal detection, env var reads). **Not feature-gated in `Cargo.toml`** despite being used *only* by `src/main.rs` — see the important caveat in §1a below. |
| **csv** | 1.4.0 | ❌ std-only | No `[features]` section at all in its manifest (always requires `std::io`). Used only by `src/drag.rs`'s peripheral fs-based `load_drag_table()` (G6/G8 fallback) and by `src/main.rs`. The engine's *own* CSV parsing for the embedded G1/G7 tables (`DragTable::from_csv_str`, `parse_embedded_drag_table`) does **not** use this crate — it hand-parses `&str` with `str::lines()`, which is core-safe. |
| **dirs** | 6.0.0 | ❌ std-only, concept doesn't apply | OS home/config-dir resolution has no embedded equivalent. Used in `bc_table_download.rs` (already `online`-gated), `pdf_dope_card.rs` (bin-only, see §1a), and `main.rs`. |
| **strsim** | 0.11.1 | (irrelevant) | Only used in `main.rs` for CLI "did you mean" suggestions. Not a library concern either way. |
| **ureq** | 3.3.0 | ❌ std-only (networking) | Already correctly `optional = true` behind the `online` feature (`ureq = { ..., optional = true }`, `online = ["dep:ureq"]`). Not a blocker — already opt-out. |
| **printpdf** | 0.7.0 | ❌ std-only (files/fonts) | Already correctly `optional = true` behind the `pdf` feature. Not a blocker — already opt-out. Its only caller, `pdf_dope_card.rs`, isn't even part of the library (see §1a). |
| **ndarray-npy** | 0.10.0 | ❌ std-only | Wraps `zip`/file I/O internally; no std/no_std toggle. Unconditional dependency, but used only inside `drag.rs`'s peripheral `.npy`-file loader (same function as the `csv` fallback above). **Confirmed second, independent thumbv7em build failure today**, not just a theoretical concern: `cargo tree -i zlib-rs` shows `ndarray-npy → zip 6.0.0 → flate2 1.1.9 → zlib-rs 0.6.3`, and both `zlib_rs` and `crc32fast` emit `error[E0463]: can't find crate for std` in the exact `cargo check --target thumbv7em-none-eabihf --no-default-features` run in §2. This crate was missing from the Phase 0 gating list in an earlier draft — corrected in §4. |
| **rayon** | 1.11.0 | ❌ std-only (OS threads) | Unconditional dependency in `Cargo.toml` — **but is not referenced anywhere in `src/` today** (`grep -rn "rayon" src/` returns nothing). This is dead weight independent of the no_std question; worth dropping regardless. |
| **tikv-jemallocator** / **mimalloc** | — | N/A | `optional = true`, target-gated (`cfg(not(target_env = "msvc"))`), used only for allocator benchmarking. Irrelevant to embedded — a real device build would supply its own `#[global_allocator]` (e.g. a static-arena allocator), which is downstream-application responsibility, not this crate's. |
| **wasm-bindgen** / **web-sys** / **js-sys** / **wasm-bindgen-futures** / **serde-wasm-bindgen** / **getrandom** (wasm target) | — | N/A | Already gated under `[target.'cfg(target_arch = "wasm32")'.dependencies]`; irrelevant to a `thumbv7em` build. |

**What "Compile-verified" means in the table above:** verified in the isolated `/tmp/nostd-dep-probe`
crate (§2a) with `default-features = false` set explicitly on each crate — **not** verified as
part of `ballistics-engine`'s own `cargo check --target thumbv7em-none-eabihf
--no-default-features`, which fails today for two of these seven (`nalgebra`/`ndarray`'s shared
`matrixmultiply` dependency and `serde`'s `serde_core` dependency) because this crate's own
`Cargo.toml` doesn't yet apply that same `default-features = false` treatment. See the correction
in §2/§2a.

### 1a. "clap/ureq/printpdf are bin-or-gated" — confirmed, with one important caveat

`ureq` and `printpdf` **are** properly gated: both are `optional = true` in `[dependencies]` and
only pulled in via the `online`/`pdf` Cargo features respectively, which are part of the
crate's `default` feature set today but can be turned off with `--no-default-features`
(confirmed: the repo's own WASM build already does this).

`clap`, `clap_complete`, `csv`, `dirs`, and `strsim`, however, are **plain, non-optional
dependencies** in `Cargo.toml`. I confirmed by `grep` that `clap`/`clap_complete`/`strsim` are
used *exclusively* by `src/main.rs` (the CLI binary), and that `csv`/`dirs` are used by
`src/main.rs` plus a small number of peripheral, non-hot-path library functions. But because
Cargo resolves dependencies at the *package* level, not per-target-based-on-actual-usage,
**`cargo check --lib --no-default-features` still pulls in the entire `clap` → `anstyle` →
`is_terminal_polyfill` chain**, and that chain is what actually broke the `thumbv7em` build (see
§2). `pdf_dope_card.rs`, notably, isn't even declared as a library module (`mod pdf_dope_card;`
appears only in `main.rs`, not `pub mod pdf_dope_card;` in `lib.rs`) — so despite `printpdf`
being a Cargo-level "lib" dependency, its only consumer is bin-only.

**Implication for any real no_std work:** Cargo.toml needs the same `optional = true` +
feature-gate treatment applied to `clap`, `clap_complete`, `csv`, `dirs`, and `strsim` that
`ureq`/`printpdf` already have, before a `--lib`-only, `--no-default-features` build can even
attempt to run without the CLI's std-only dependency chain in tow. This is Phase 0 in the plan
below — mechanical, not risky, but currently missing.

---

## 2. `thumbv7em-none-eabihf` compile attempt

```
$ rustup target add thumbv7em-none-eabihf
info: downloading component 'rust-std' for 'thumbv7em-none-eabihf'
info: installing component 'rust-std' for 'thumbv7em-none-eabihf'

$ cargo check --target thumbv7em-none-eabihf --no-default-features
   ...
error[E0463]: can't find crate for `std`
  = note: `std` is required by `is_terminal_polyfill` because it does not declare `#![no_std]`
error: cannot resolve a prelude import
error[E0463]: can't find crate for `std`
   --> .../is_terminal_polyfill-1.70.2/src/lib.rs:43:5
    |
 43 |     std::fs::File,
error[E0463]: can't find crate for `std`
   --> .../memchr-2.8.0/src/lib.rs:198:1
    |
198 | extern crate std;
error[E0463]: can't find crate for `std`
  = note: `std` is required by `anstyle` because it does not declare `#![no_std]`
error[E0463]: can't find crate for `std`
  = note: `std` is required by `matrixmultiply` because it does not declare `#![no_std]`
error[E0463]: can't find crate for `std`
  = note: `std` is required by `serde_core` because it does not declare `#![no_std]`
error[E0463]: can't find crate for `std`
  = note: `std` is required by `zlib_rs` because it does not declare `#![no_std]`
error[E0463]: can't find crate for `std`
  = note: `std` is required by `crc32fast` because it does not declare `#![no_std]`
... (cascading errors: `derive`, `concat!` etc. "not found" once `core`'s prelude
    substitutes for `std`'s prelude — these are downstream symptoms of missing-std roots, not
    independent bugs; the full run produces 25 top-level `E0463`s naming 11 distinct crates)
```

Same result with `--lib` only (`cargo check --lib --target thumbv7em-none-eabihf
--no-default-features`).

**Correction (flagged by review, 2026-07-15):** an earlier draft of this section claimed "every
single compile error in this run traces back to the `clap` chain... none of them originate in...
the numeric dependency chain." **That claim is false**, and is disprovable by grepping this exact
transcript: `matrixmultiply` (a direct dependency of *both* nalgebra and ndarray — confirmed via
`cargo tree -i matrixmultiply`) and `serde_core` (a direct dependency of serde — confirmed via
`cargo tree -i serde_core`) each independently emit `error[E0463]: can't find crate for std` in
this run, squarely inside the "numeric/serialization dependency chain" §2a calls
"compile-verified." A third, independent chain also fails: `ndarray-npy` (an unconditional
dependency, flagged in §1's table) pulls in `zip → flate2 → zlib-rs` and `crc32fast`, both of
which likewise emit `E0463` (confirmed via `cargo tree -i zlib-rs` / `cargo tree -i crc32fast`,
and reproduced deterministically across two independent runs of the command above).

**Root cause is the same shape of bug, but at two different levels — that's the important
distinction, not "only clap":**
- The `clap`/`ndarray-npy` chains fail because `clap`, `clap_complete`, `csv`, `dirs`, `strsim`,
  and `ndarray-npy` are plain, non-optional dependencies with **no feature gate at all** — Cargo
  always resolves their std-requiring defaults, and `ndarray-npy` doesn't even expose a
  `[features]` section to opt out of even if we wanted to.
- The `matrixmultiply`/`serde_core` failures happen for a *different* reason: nalgebra, ndarray,
  and serde **are** no_std-capable — §2a proves that, in isolation — but `ballistics-engine`'s
  own `Cargo.toml` never sets `default-features = false` on them, so Cargo resolves their
  *default* (std-requiring) feature sets when building this crate. The isolated
  `/tmp/nostd-dep-probe` crate in §2a explicitly sets `default-features = false` on all seven
  numeric/serialization crates and therefore doesn't hit this. The narrower §2a claim ("the
  numeric/serialization crates are themselves no_std-capable") still holds; the broader claim
  this section originally made ("therefore this crate's build of them is clean today") does not.

Corrected summary: **two** root causes, not one — (a) un-gated CLI-only dependencies (`clap`,
`clap_complete`, `csv`, `dirs`, `strsim`, `ndarray-npy`), and (b) numeric/serialization
dependencies resolving to their default (std) feature sets because this crate's `Cargo.toml`
hasn't yet opted them into their `alloc`/no_std configurations. Both are Cargo.toml-level,
mechanical fixes — neither implicates the engine's own physics code — but Phase 0 (§4) needs to
address both; `ndarray-npy`'s gating was missing from that list in an earlier draft and is added
below.

### 2a. Isolating the real question: does the numeric/serialization chain compile no_std?

To separate "blocked by the CLI's dependency chain" from "blocked by the math dependencies," I
built a standalone probe crate (outside this worktree, in `/tmp/nostd-dep-probe/`, touching
nothing in `ballistics-engine`) pulling in the seven crates our physics code actually needs,
each configured with `default-features = false` plus their published `alloc`/no_std feature
flags, and exercising: `nalgebra::Vector3` arithmetic, `ndarray::Array1`, `rand::rngs::StdRng`
seeded manually, `rand_distr::Normal`, `#[derive(Serialize, Deserialize)]` via `serde`, and
`#[derive(thiserror::Error)]`:

```toml
[dependencies]
nalgebra   = { version = "=0.35.0", default-features = false, features = ["alloc", "libm"] }
ndarray    = { version = "=0.17.2", default-features = false }
rand       = { version = "=0.10.2", default-features = false, features = ["alloc", "std_rng"] }
rand_distr = { version = "=0.6.0",  default-features = false, features = ["alloc"] }
serde      = { version = "=1.0.228", default-features = false, features = ["alloc", "derive"] }
serde_json = { version = "1.0", default-features = false, features = ["alloc"] }
thiserror  = { version = "=2.0.18", default-features = false }
```

```
$ cargo check --target thumbv7em-none-eabihf
   Compiling ... (39 packages)
    Checking nalgebra v0.35.0
    Checking nostd-dep-probe v0.0.0 (/private/tmp/nostd-dep-probe)
warning: struct `Foo` is never constructed
    Finished `dev` profile [unoptimized + debuginfo] target(s) in 0.06s
```

**Clean compile, zero errors**, on the real target, with the real (locked) versions this repo
uses — **when each crate is explicitly configured with `default-features = false` plus its
documented `alloc`/no_std feature flags**, exactly as the `Cargo.toml` snippet above shows. This
is strong evidence that the numeric/serialization crates are *themselves* no_std-capable. It is
**not** evidence that `ballistics-engine`'s current build of them is clean — see the correction
in §2: because this crate's own `Cargo.toml` doesn't yet set `default-features = false` on
nalgebra/ndarray/serde (etc.), `matrixmultiply` and `serde_core` fail in the real crate today
with the exact same `E0463` class shown here, and `ndarray-npy` (which has no `[features]` to
even turn off) fails for the same "no gate at all" reason as `clap`. The gap between "this
dependency graph *can* compile no_std" (proven here) and "this crate's dependency graph *does*
compile no_std today" (false, per §2) is precisely the `default-features = false` +
feature-list work scoped into Phase 0 below — mechanical, not a discovery risk, but real work
that an earlier draft of this report undercounted by treating it as already done.

---

## 3. Inventory of our own std surface (core solve path + adjacent physics modules)

Grepped every `pub mod` under `lib.rs` that participates in the numeric solve path — `cli_api.rs`,
`trajectory_integration.rs`, `drag.rs`, `drag_tables.rs`, `atmosphere.rs`, `wind.rs`,
`wind_shear.rs`, `moving_target.rs`, `derivatives.rs`, `fast_trajectory.rs`,
`trajectory_sampling.rs`, `angle_calculations.rs`, `spin_drift.rs`, `spin_drift_advanced.rs`,
`spin_decay.rs`, `precession_nutation.rs`, `pitch_damping.rs`, `aerodynamic_jump.rs`,
`transonic_drag.rs`, `reynolds.rs`, `form_factor.rs`, `stability.rs`, `stability_advanced.rs`,
`monte_carlo.rs`, `bc_estimation.rs`, `cluster_bc.rs`, `constants.rs`, `drag_model.rs`,
`trajectory_observation.rs`, `bc_table.rs`, `bc_table_5d.rs` — for `std::fs`, `std::io`,
`std::time`, `std::thread`, `std::env`, `std::process`, `std::net`, `std::sync`, `HashMap`,
`println!`/`eprintln!`, and any other `std::` prefix, then manually classified each hit as
**core-blocking** (touches the numeric path itself) or **peripheral** (test-only, CLI/feature-
gated, or an auxiliary loader not reachable from the default `TrajectorySolver` path).

| Pattern | Files | Live-code hits | Classification |
|---|---|---|---|
| `std::fs`, `std::io`, `std::time`, `std::thread`, `std::env`, `std::process`, `std::net` | `drag.rs`, `bc_table.rs`, `bc_table_5d.rs` | `drag.rs`: 1 fs-based loader (`load_drag_table`, G6/G8 only — G1/G7 are baked in via `include_str!` and never touch `std::fs`); 1 `std::time::Instant` (inside a `#[test]` fn, not compiled for the target). `bc_table.rs`: `File`/`BufReader` in `load()`, called only from `main.rs` (CLI). `bc_table_5d.rs`: `File`/`BufReader` in a disk `load()` + `std::fs::read_dir` for auto-discovery, called only from `main.rs` and (notably) **not** from `wasm.rs`, which already routes through a byte-based `from_reader`/`from_bytes` constructor with an explicit code comment: *"there is no `std::fs`... the host JS/Node layer fetches the `.bin`."* | **Peripheral.** None of this is reachable from `cli_api.rs`'s core solve routines — confirmed by grepping for `bc_table::`/`bc_table_5d::` callers, which are only `main.rs` and `wasm.rs`. `wasm.rs` already demonstrates the exact pattern an embedded port would reuse (bytes in, no fs). |
| `std::sync::Once` | `cli_api.rs` (1 call site: `custom_drag_denominator`'s one-time stderr warning guard) | 1 | **Core-blocking, but tiny.** Needs a no_std-safe once-guard (hand-rolled `AtomicBool`/`critical-section`, or drop the warning). |
| `std::sync::LazyLock` | `drag.rs` (4 live statics: `G1_DRAG_TABLE`, `G7_DRAG_TABLE`, `G6_DRAG_TABLE`, `G8_DRAG_TABLE`) + `drag_tables.rs` (8 more: `G1`/`G7`/`G2`/`G5`/`G6`/`G8`/`GI`/`GS`, all `#[allow(dead_code)]` and, per `grep -rn "drag_tables::" src/`, never referenced anywhere in `src/`) | **12** — an earlier draft of this table claimed "2 (one each)", which undercounted `drag.rs` alone by 4x and missed `drag_tables.rs` entirely | **Core-blocking, but tiny, for the 4 live sites in `drag.rs`.** Needs a no_std lazy-init substitute (e.g. a `critical-section`-backed lazy cell). **The 8 in `drag_tables.rs` are dead code end-to-end** — the same category as the unused `rayon` dependency flagged in §1 — so the cheaper Phase-0 move is to delete `drag_tables.rs` outright rather than port 8 statics nothing calls. Either way they'd break a `#![no_std]` build today: `#[allow(dead_code)]` only silences the lint, it doesn't stop the compiler from type-checking and compiling the item. |
| `std::error::Error` | `cli_api.rs` (`impl Error for BallisticsError {}`) | 1 | **Core-blocking, trivial fix.** `core::error::Error` is stable since Rust 1.81 (repo's toolchain: 1.91.1) — this is a one-line swap. |
| `std::f64::consts::*`, `std::cmp::Ordering`, `std::fmt` | `cli_api.rs` (9+ sites), `monte_carlo.rs` (2 sites) | ~11 | **Core-blocking, but mechanical.** `core::f64::consts`, `core::cmp::Ordering`, `core::fmt` are byte-for-byte equivalent; this is a find/replace, no logic change. |
| `println!` | `wind.rs`, `stability.rs`, `stability_advanced.rs`, `trajectory_integration.rs` | **0** (all sites are inside `#[cfg(test)] mod tests { ... }`, confirmed by comparing line numbers against each file's `mod tests` boundary) | **Not blocking.** Test code isn't compiled into the library artifact for a downstream target unless someone runs `cargo test` there, which is not the deployment model. |
| `eprintln!` | `cli_api.rs` (1, same `Once`-guarded warning above), `trajectory_integration.rs` (6: RK45 adaptive-stepper diagnostics — non-finite tolerance, invalid `max_step`, non-finite minimum-step trial, max-iteration-limit warnings) | 7 | **Core-blocking, narrow, all diagnostic-only.** Every live site is a *warning on a degenerate/edge-case numeric condition* (bad input, solver not converging) — none are on the per-step numeric happy path. Needs redirecting through a `cfg(feature = "std")`-gated shim or a `log`/`defmt`-style facade (standard embedded pattern); the physics logic itself is untouched. |
| `HashMap` | `bc_table_5d.rs` (`tables: HashMap<i32, Bc5dTable>`, a caliber-keyed cache for the 5D BC correction feature), `trajectory_integration.rs` (one function explicitly commented `"Convert to Python-friendly format"` / `"legacy ... test function"`, returning `Vec<HashMap<String, f64>>`) | 2 sites | **Peripheral/legacy.** The `bc_table_5d.rs` cache belongs to the same fs-adjacent, main.rs/wasm.rs-only subsystem as the row above (not reachable from the core `TrajectorySolver` path). The `trajectory_integration.rs` `HashMap` is an explicitly-labeled legacy binding-compatibility shim, not the production solve path (`fast_integrate_with_segments` and friends are the real entry points). Both are mechanically swappable to `alloc::collections::BTreeMap` if ever needed no_std. |
| `atmosphere.rs`, `wind.rs`, `moving_target.rs` | — | **0** | **Fully clean already.** No `std::fs/io/time/thread`, no `HashMap`, no live `println!`/`eprintln!`. |
| `rayon` (thread-based parallelism) | — | **0** references anywhere in `src/` | **N/A** — unused dependency, not a blocker (see §1). |
| `ffi.rs` (C ABI for iOS/Android) | `src/ffi.rs` | 97 `std::` references (`std::alloc`, `std::mem`, `std::ptr`, `std::slice`, `std::os::raw::c_char`) | **Out of scope for this spike.** All of these have `core`/`alloc` equivalents in principle (this is exactly the shape of a C-ABI FFI layer), but it's a separate, sizeable surface deserving its own review — not part of "core solve path" as scoped by the ticket. Flagged here only so it isn't mistaken for a clean result; a future embedded port would need to decide whether it reuses this layer or defines a leaner one. |

**Bottom line:** across roughly 24,000 lines spanning the entire physics/solver module set
(cli_api.rs's solver internals down through spin drift, precession/nutation, stability, drag,
atmosphere, wind, and moving-target lead), the code that would need to change to compile under
`#![no_std]` is: 5 live synchronization-primitive call sites (1 `Once` + 4 `LazyLock`, all in
`cli_api.rs`/`drag.rs`), 8 more `LazyLock` statics in the entirely-dead-code `drag_tables.rs`
(cheaper to delete than port), one trait impl, about a dozen mechanically-portable
`std::`→`core::` path swaps, seven diagnostic `eprintln!`s, and two `HashMap`s in
already-peripheral subsystems. Everything else — the actual physics — is already written in
ordinary, allocation-only Rust. (An earlier draft undercounted the `LazyLock` sites as "two"; the
corrected count doesn't change the mechanical, no-discovery-risk character of the list — see
§4 — but is recorded here for accuracy, per review.)

---

## 4. Phased plan (if pursued)

**Phase 0 — Cargo.toml surgery (0.5–1 day).**
Gate `clap`, `clap_complete`, `csv`, `dirs`, `strsim`, **and `ndarray-npy`** behind a new
opt-in-by-default `cli` feature, mirroring the existing `online`/`pdf` pattern (`optional = true`
+ `cli = ["dep:clap", "dep:clap_complete", "dep:csv", "dep:dirs", "dep:strsim",
"dep:ndarray-npy"]`). **`ndarray-npy` was missing from this list in an earlier draft**, even
though it's the source of a second, independent thumbv7em failure today (§1, §2:
`ndarray-npy` → `zip` → `flate2`/`crc32fast` → `zlib-rs`, all std-only with no feature toggle of
their own to fall back on). Introduce a `std` feature (on by default) that forwards to each
numeric dependency's own `std` feature (`nalgebra/std`, `ndarray/std`, `rand/std`,
`rand_distr/std`, `serde/std`, `serde_json/std`, `thiserror/std`), with those dependencies
switched to `default-features = false` plus explicit `alloc`/`libm` feature lists **directly in
this crate's own `Cargo.toml`** so the "no `std` feature" configuration resolves to the no_std
variant validated in isolation in §2a — this is the specific step that's currently missing and
that makes `matrixmultiply`/`serde_core` fail in the real crate today (§2). Verify
`cargo check --lib --no-default-features` (i.e. no `std`, no `cli`) now produces a *clean
dependency graph* — no longer pulling in `clap`, `ndarray-npy`'s zip chain, or any numeric
crate's default std-requiring feature set.

**Phase 1 — fix the enumerated blockers (2–4 days).**
Add `#![cfg_attr(not(feature = "std"), no_std)]` + `extern crate alloc;` to `lib.rs`. Apply the
fixes cataloged in §3, all of which are now fully known in advance (no discovery risk left, per
the corrected counts there): swap `std::sync::Once` (1 site) and `std::sync::LazyLock` (4 live
sites) for a no_std-safe lazy-init in `cli_api.rs`/`drag.rs`, and delete `drag_tables.rs`'s 8
dead-code `LazyLock` statics outright rather than port them; swap
`impl std::error::Error` → `impl core::error::Error`; mechanical
`std::f64`/`std::cmp`/`std::fmt` → `core::` path swaps; feature-gate the ~7 `eprintln!`
diagnostics behind `cfg(feature = "std")` (or route through a minimal logging trait); bake G6/G8
drag tables into the binary via `include_str!` the same way G1/G7 already are, and
feature-gate the fs-based `load_drag_table`/`find_drag_tables_dir` behind `std`; feature-gate
`bc_table.rs`/`bc_table_5d.rs`'s disk-`load()`/`read_dir` entry points behind `std` (their
byte-based constructors already work no_std, per the existing `wasm.rs` precedent).

**Phase 2 — prove it (1–2 days).**
`cargo check --lib --target thumbv7em-none-eabihf --no-default-features --features
<new-core-embedded-alias>` until clean; pin this as a CI job so a future PR can't silently
reintroduce a std dependency into the core physics modules.

**Phase 3 — on-target validation (unscoped/exploratory, separate follow-up).**
This is the part that actually decides whether to ship it, and it's explicitly *not* answered
by "it compiles":
- Flash to real Cortex-M4F/M7 hardware and measure actual single-shot trajectory-solve latency,
  peak RAM, and flash footprint for a representative rangefinder/dope-computer use case.
- Resolve the f64-on-single-precision-FPU question (see Open questions below) with real numbers,
  not assumptions.
- Decide the fate of Monte Carlo (needs a real entropy source on-device — a product decision, not
  just a technical one), the `solve_json`/`solve_v1` JSON I/O layer (String/serde_json-heavy;
  possibly worth a leaner byte-oriented ABI for a memory-constrained MCU, akin to `ffi.rs`), and
  whether `ffi.rs` is reused or a narrower embedded-specific ABI is written.
- Run the existing fuzz-crate's known numeric-robustness gaps (see repo memory: unvalidated
  `Ok(NaN)` paths) against the embedded slice specifically — a NaN or panic that merely returns a
  bad number on desktop is a `panic = "abort"` hard reset on a rangefinder with no OS to catch it.

**Effort estimate:** Phases 0–2 (getting `cargo check --target thumbv7em-none-eabihf` green for
a real "core" feature slice) is roughly **one focused engineer-week**, given how narrow and
already-enumerated the blocker list turned out to be. Phase 3 (hardware bring-up and the
performance/robustness verdict) is open-ended exploratory work and should be scoped as its own
ticket with its own hardware budget before any commitment to ship.

---

## 5. Open questions / risks not resolved by this spike

1. **f64 on a single-precision FPU.** `thumbv7em-none-eabihf`'s "hf" hard-float ABI covers
   single-precision (Cortex-M4F's VFPv4-sp-d16 has no double-precision unit); this engine's
   physics is `f64` throughout. `simba`'s `libm` feature (validated in §2a) supplies software
   implementations for transcendentals (sin/cos/sqrt/etc.), but ordinary `f64` add/multiply/divide
   also fall back to compiler-generated software floating point on this class of part — not just
   the transcendental calls. That's a real, unmeasured latency cost this spike did not quantify.
   Cortex-M7 parts with a double-precision FPU option (e.g. STM32H7) would avoid it; a
   single-precision-only rangefinder MCU would not. This must be measured on real hardware
   (Phase 3), not assumed either way.
2. **Table storage / RAM budget.** The embedded G1/G7 tables are baked in as CSV text
   (`include_str!`) and parsed into heap `Vec<f64>` on first use via `LazyLock`. That's fine on
   desktop/WASM; on a device with a few hundred KB of RAM, runtime string-parsing plus heap
   allocation for lookup tables may be worth replacing with genuinely `const`/build-time-baked
   binary tables. Out of scope for "does it compile," but likely necessary before "is it good."
3. **Global allocator.** A `no_std` + `alloc` library requires the *consuming application* to
   register a `#[global_allocator]` (e.g. a static-arena allocator crate). That's downstream
   responsibility, not ours, but worth stating explicitly since it's a hard prerequisite.
4. **Monte Carlo's entropy source.** `rand`'s `std_rng` works no_std with a manually-supplied
   seed (validated in §2a), but loses `sys_rng`/`thread_rng`. A real embedded consumer needs an
   actual entropy source (ADC noise, RTC jitter, a hardware TRNG peripheral) wired through
   `getrandom`'s custom-backend hook or a directly-supplied seed — a per-target integration
   detail, and possibly a reason to simply exclude Monte Carlo from a v1 embedded feature slice.
5. **Panic semantics.** No code in this crate defines a `#[panic_handler]` (correctly — that's an
   application concern), but a `.unwrap()` or index-out-of-bounds that today unwinds cleanly on
   desktop typically means `panic = "abort"` and a hard device reset on a microcontroller with no
   OS underneath. This raises the stakes on the crate's existing numeric-robustness gaps (see
   repo memory on the fuzz crate's open `Ok(NaN)` finding) — worth a dedicated robustness pass
   before shipping an embedded target, independent of the no_std mechanics themselves.
6. **`ffi.rs` (iOS/Android C ABI)** is a structurally similar but separately-scoped 97-reference
   std surface that this spike deliberately did not audit in depth (it's not "core solve path"
   per the ticket). A real embedded port needs to decide whether it reuses this layer (after its
   own no_std pass) or defines a narrower ABI purpose-built for a microcontroller.

---

## 6. GO / NO-GO

**Qualified GO — for Phases 0–2 only.**

The compile-time feasibility question this spike was scoped to answer comes back positive, with
one correction from a review pass: the numeric/serialization dependency chain is
compile-verified no_std-capable **in isolation** (§2a) on a real Cortex-M4F target — it does not
yet compile that way *inside this crate*, because `Cargo.toml` hasn't been given the same
`default-features = false` treatment the probe crate got (§2). Our own code's std surface on the
actual solve path is still small, fully enumerated, and mechanical to fix (no discovery risk
remains — every blocker in §3, including the corrected `LazyLock` count, is a known, bounded
change), and Phase 0's scope now explicitly includes `ndarray-npy` alongside the CLI-only
crates. This is not a "rewrite the physics" problem; it's a "restructure some Cargo features and
swap about two dozen call sites" problem, still estimated at about one engineer-week — the
correction changes what Phase 0 needs to gate, not the size of Phase 1's mechanical work.

**This is explicitly not a GO to ship an embedded target.** Whether the result is fast enough
(f64-on-single-FPU performance), small enough (RAM for table storage), and safe enough
(numeric-robustness under abort-on-panic) on real rangefinder/dope-computer hardware is
completely unanswered by a `cargo check` passing, and none of those questions can be resolved
without Phase 3 hardware bring-up. Recommend filing Phases 0–2 as a bounded, low-risk follow-up
ticket, and treating Phase 3 as a separate, hardware-budgeted spike gated on Phases 0–2 landing
first.

---

## Commands run (reference)

```
git worktree add .claude/worktrees/mba-1324 -b mba-1324-nostd-spike main
cargo test                                                      # baseline: all green
cargo clippy --all-targets                                      # baseline: 0 warnings
cargo check --target wasm32-unknown-unknown --no-default-features   # baseline: succeeds (pre-existing warnings only)
rustup target add thumbv7em-none-eabihf
cargo check --target thumbv7em-none-eabihf --no-default-features        # fails: clap AND numeric/serde chains want std
cargo check --lib --target thumbv7em-none-eabihf --no-default-features  # same failures, same root causes
# standalone probe crate at /tmp/nostd-dep-probe (outside the worktree, no repo files touched):
cargo check --target thumbv7em-none-eabihf   # nalgebra+ndarray+rand+rand_distr+serde+serde_json+thiserror: clean
                                              # (only with default-features=false set explicitly — see §2a)

# post-review correction (fix pass, 2026-07-15) — isolating exactly which crates emit E0463
# in the REAL crate's --no-default-features build, not just the clap chain an earlier draft
# blamed exclusively:
cargo tree -i matrixmultiply      # confirms: nalgebra AND ndarray both pull it in
cargo tree -i serde_core          # confirms: serde pulls it in
cargo tree -i zlib-rs             # confirms: ndarray-npy -> zip -> flate2 -> zlib-rs
cargo tree -i crc32fast           # confirms: ndarray-npy -> zip -> flate2 -> crc32fast
cargo check --target thumbv7em-none-eabihf --no-default-features 2>&1 | grep -B2 "E0463"
  # confirms matrixmultiply, serde_core, zlib_rs, crc32fast, anstyle, anstyle_query, clap_lex,
  # is_terminal_polyfill, strsim, ucd_trie all independently emit E0463 in this one run
  # (25 top-level E0463s naming 11 distinct crates); reproduced deterministically twice
grep -rn "drag_tables::" src/         # confirms drag_tables.rs's 8 LazyLock statics are dead code
```

Full command transcripts and outputs are in the accompanying execution log,
`.superpowers/mba-1324-report.md`, in this worktree.
