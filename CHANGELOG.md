# Changelog

All notable changes to the ballistics-engine project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.20.0] - 2026-06-21

Dependency modernization — a coordinated upgrade of the major dependencies. No new
features, CLI flags, or physics; numerical trajectory output is unchanged. The one
observable difference is in Monte Carlo: under the new RNG the per-run sample
*sequence* differs, but the statistical results (mean, std-dev, CEP, confidence
ellipse) are unchanged.

### Changed
- **`ndarray` 0.15 → 0.17** and **`ndarray-npy` 0.8 → 0.10** (upgraded together; the
  drag-table `.npy` loader is unaffected).
- **`rand` 0.8 → 0.9** and **`rand_distr` 0.4 → 0.5**. Monte Carlo uses an unseeded
  RNG, so this only changes the sample sequence, not the statistical behavior.
- **`ureq` 2.12 → 3.1** (online feature). Migrated to the 3.x API; the public error
  mapping (`ServerError(code, body)`, timeout vs. network) and the unbounded
  BC-table read are preserved.
- **`getrandom` 0.2 → 0.3** with the `wasm_js` backend feature (WASM builds).
- CI: `actions/upload-artifact` 6 → 7, `actions/download-artifact` 7 → 8,
  `cross-platform-actions/action` 0.30 → 0.32.

## [0.19.0] - 2026-06-20

Adds downrange-segmented wind to the CLI and WASM, and **corrects the wind-direction
sign on the `cli_api`/`TrajectorySolver` path** (a breaking output change — see below).

### Added
- **Downrange-segmented wind** — wind can now vary along the trajectory (e.g. a muzzle
  reading plus downrange sensor stations).
  - CLI: repeatable `--wind-segment SPEED:ANGLE:UNTIL_DISTANCE` on the `trajectory`
    command. SPEED and UNTIL_DISTANCE follow `--units` (mph & yards imperial, m/s &
    meters metric); ANGLE is degrees in the wind-FROM convention (same as
    `--wind-direction`). Each segment applies from the previous boundary out to its
    UNTIL_DISTANCE; wind is a step function with **zero wind beyond the last segment**
    (a coverage warning is printed when the segments don't reach `--max-range`).
    Segments **override** scalar `--wind-speed`/`--wind-direction`, and are **not
    compatible with `--enable-wind-shear`** (rejected with an error).
  - WASM: `runCommand` accepts the same repeatable `--wind-segment` flag, and the
    `Calculator` gains `addWindSegment(speedMph, directionDeg, untilYards)` and
    `clearWindSegments()`.
  - API: `TrajectorySolver::set_wind_segments(Vec<WindSegment>)` performs a per-step
    downrange lookup via `WindSock`, preserving all solver physics (Coriolis, spin
    drift, aerodynamic jump, zero-finding, full `TrajectoryResult` output). A solve
    with no segments is numerically identical to before. New public
    `wind::parse_wind_segment_str(&str, imperial)`.

### Fixed
- **CLI `--twist-right` can now select left-hand twist.** It was a presence-only flag
  (always true; `--twist-right false` errored), so left-hand twist — and the
  aerodynamic-jump direction reversal it produces — was unreachable from the CLI. It now
  takes an optional value: `--twist-right` / `--twist-right true` = right-hand (default),
  `--twist-right false` = left-hand (the bare flag still works, so existing usage is
  unaffected).

### Changed (BREAKING)
- **Wind direction sign corrected on the `cli_api`/`TrajectorySolver` path.** The scalar
  wind vector was built with the opposite sign to the standard wind-clock convention
  (and to the engine's own `WindSock`, the Monte-Carlo path, and the Python/Ruby
  bindings): `--wind-direction 0` produced a *tailwind* and `90°` drifted the bullet the
  wrong way. It is now correct: **0° = headwind, 90° = from the right (drifts left),
  180° = tailwind, 270° = from the left** — matching `WindSock`/the bindings.
  **Trajectories computed with a non-zero wind direction will change** (head/tail effect
  and crosswind drift direction flip); no change when wind speed is 0. This affects
  **every consumer of the `TrajectorySolver`/cli_api path**: the CLI `trajectory`,
  `wind-card`, `range-table`, and `come-ups` subcommands; the WASM `runCommand` and
  `Calculator`; the C FFI (`FFIWindConditions.direction`, for iOS/Android callers); and
  any direct `TrajectorySolver` user. The Monte-Carlo path and the Python/Ruby bindings
  already used the correct convention and are unchanged. The aerodynamic-jump crosswind sign
  on this path was aligned to match (right twist + wind-from-right → impact up, drift
  left), bringing `cli_api` into agreement with the fast-integrate path. The wind-shear
  path was already sign-aligned with the scalar path and flips with it.
- The `--wind-direction` help text was corrected (previously the inaccurate
  "0=North, 90=East").

## [0.18.3] - 2026-06-19

Extends aerodynamic jump (MBA-959) to the fast-integrate path. Default-off; no change to existing trajectories.

### Added
- **Aerodynamic jump in `fast_trajectory::fast_integrate`** — the fast fixed-step kernel (used by the `ballistics-engine-py` `fast_integrate` binding and the Monte-Carlo path) now applies aerodynamic jump when `enable_aerodynamic_jump` is set, by perturbing the prebuilt launch velocity vertically. New public `fast_trajectory::aerodynamic_jump_launch_offset_rad(inputs, atmo_params)` returns the Litz vertical offset (radians) from the engine's Miller Sg; crosswind is taken from `wind_speed`/`wind_angle` (BallisticInputs convention: 0 = headwind, +90° = from the right). Previously only the `cli_api::TrajectorySolver` (Euler/RK4/RK45) path applied AJ; the fast path ignored the flag.

## [0.18.2] - 2026-06-19

Adds an opt-in aerodynamic-jump effect (MBA-959). Default-off and additive — no change to existing trajectories.

### Added
- **Aerodynamic jump** (`enable_aerodynamic_jump`, default off): the vertical point-of-impact shift a crosswind imparts to a spin-stabilized bullet, applied as a muzzle launch-angle perturbation in the solver. Uses Bryan Litz's estimator `Y = 0.01·Sg − 0.0024·L + 0.032` (MOA per mph of crosswind), fed by the engine's own Miller stability factor (Sg) and bullet length in calibers. Exposed on `TrajectoryResult.aerodynamic_jump`, via the CLI `--enable-aerodynamic-jump` flag, and as the public `aerodynamic_jump::litz_crosswind_jump_moa()`. The jump is purely vertical; direction follows Litz (right twist: crosswind from the right → impact up). The zero is found on the bare bore, so the jump appears as an additive POI shift rather than being absorbed by the zero. Most accurate near Sg ≈ 1.75 (Litz regression validity).

### Notes
- The legacy heuristic `aerodynamic_jump::calculate_aerodynamic_jump` is retained for backward compatibility but is superseded by the Litz estimator and is no longer used by the solver.

## [0.18.1] - 2026-06-18

Two follow-up fixes after 0.18.0. Neither changes trajectory output.

### Fixed
- **BCCR drag-correction loader** now verifies the stored CRC32 (a standard IEEE/zlib CRC over the data section) instead of reading and discarding it — corrupt-but-in-range files are rejected with `ChecksumMismatch` instead of loading silently. Gated on a non-zero stored checksum so older checksum-less files still load (MBA-953).
- **Precession/nutation angular diagnostics** — corrected four dimensional errors: the precession and nutation frequencies are now rad/s via the standard epicyclic form `φ = (Iₓp/2Iy)[1 ± √(1 − 1/Sg)]` (with the dimensionless Miller stability factor), and the yaw no longer random-walks. Diagnostic-only — these feed only the opt-in `max_yaw_angle` / `max_precession_angle` outputs (`--enable-precession-nutation`); no trajectory change (MBA-941).

## [0.18.0] - 2026-06-18

A cross-solver consistency and physics-correction pass (JIRA MBA-938..958). 17 commits since v0.17.0; build and full test suite (210+) green. Several changes correct real bugs that **shift numerical results** on the Coriolis, transonic/subsonic drag, spin-drift-at-altitude, and WASM-zero paths — see Changed and validate against a reference before deploying downstream.

### Fixed
- **Coriolis deflection direction** — corrected to the physical `-2 Ω×v` with the right lateral ω sign; a Northern-hemisphere shot now drifts right (was left), validated vs py_ballisticcalc. The fast/Monte-Carlo path (used by the Python/Ruby bindings) applied no Coriolis at all and now matches the default solver to ~1%.
- **Transonic drag resolution** — the engine silently used a coarse 21-point G1/G7 fallback (a runtime loader looked for a directory that never exists at runtime); the high-resolution 79/84-point tables are now baked in via `include_str!`. Transonic drop error vs py_ballisticcalc fell from +43" to +3" at 1300 yd.
- **Custom drag tables** (`custom_drag_table`) are now honored by all three solvers (were plumbed onto the inputs but never read — the G-model was always used).
- **Named bullet shapes** (boat/round/flat in `bullet_model`) now drive the transonic shape on every solver (the integrate_trajectory path ignored the name and used a caliber/weight heuristic).
- **`use_form_factor`** is honored consistently across all three solvers (was derivatives-only; no-op by default).
- **WASM standalone zero** targets the line of sight (sight height), not the bore line — it was off by the sight-height angle (~2 MOA at 100 yd / 2" sight).
- **PDF output** is gated behind the `pdf` feature so `--no-default-features` builds the binary.
- **Sampled `--full` CSV** reports drop/drift in inches (Imperial) / meters (Metric), matching the table and API (was yards).
- Removed a dead WASM `--wind-dir-std` setter; documented the wind-direction-spread approximation.

### Changed
> **FLAGGED for downstream:** these correct real bugs but **shift numerical output**. Re-validate dope and any cached/expected outputs.
- **Coriolis** — direction fix plus the previously-absent Monte-Carlo/binding application; affects any Coriolis-enabled trajectory.
- **Transonic/subsonic drag** — the high-resolution G1/G7 tables shift drag below ~Mach 1.3; the low-velocity **Reynolds correction was removed** (it ran on only one of three solvers, py_ballisticcalc does not model it, and the default solver already matches pbc subsonically), shifting subsonic (<200 m/s) shots on the binding path.
- **Spin drift at altitude** — the Miller stability density correction is now the canonical linear `(T/T0)(P0/P)` (was `sqrt`), raising spin drift at altitude (sea level unchanged). The same correction now also feeds the Magnus / yaw-of-repose Sg (off by default).
- **WASM zero** shifts by the sight-height angle (~2 MOA at 100 yd) — now a proper sight-line zero.

### Added
- `transonic_drag::resolve_projectile_shape(...)` — shared bullet-shape resolution (name first, then heuristic) used by all three solvers.
- `TrajectoryParams.ground_threshold` — the integrate_trajectory ground plane is configurable (was hardcoded `-1000.0`; default preserved).

### Performance (output-identical)
- Velocity-BC segments are estimated **once** at solver setup on the integrate_trajectory path (the Python-binding hot path) instead of being rebuilt — allocating a string + a vector — on every derivative evaluation. Byte-identical output.

## [0.17.0] - 2026-06-17

An autonomous adversarial hardening pass plus a follow-up altitude/temperature fix. 82 fixes across 78 commits since v0.16.0; build and full test suite green. Most changes are correctness/robustness fixes that preserve output, but several correct real bugs that **shift numerical results** — see Changed below and validate against a reference before deploying downstream.

### Fixed
- **RK45 time-of-flight desync** — the default solver advanced `time` by the *next* step's dt, corrupting time_of_flight and per-point/sample times.
- **Arden-Buck vapor pressure** — the linear factor sat outside `exp()`, over-estimating saturation pressure (~7x) and corrupting humidity-dependent air density.
- **Gyroscopic spin-drift twist sign** — applied twice (canceling), so left-twist barrels drifted the wrong way (right-twist/default unchanged).
- **API twist-rate unit** — sent in meters where inches/turn was expected (~33x off); now sent verbatim as inches/turn across API/FFI/WASM.
- **Wind-layer extrapolation** — above the top custom layer the profile extrapolated garbage instead of clamping; WindSock cursor no longer skips segments on large query jumps.
- **Transonic drag** applied exactly once everywhere (was double-applied in `derivatives.rs`, missing in `fast_trajectory.rs`, shape-misclassified in `cli_api.rs`); **Reynolds drag** applied once and capped at 5x (was double-applied / uncapped ~16x spike).
- **Euler solver** now uses the shared acceleration kernel — previously ignored the ballistic coefficient (diverging up to ~2.3x) and omitted Magnus/Coriolis.
- **Projectile shape (default solver)** received caliber in meters instead of inches, under-applying transonic drag for G1/G2/G5/G6/G8/GI/GS.
- **Monte Carlo** — fixed long-range truncation and a 100x humidity mis-scale (fraction vs percent); BC-correction table no longer silently zeros the BC when queried with the other drag model.
- **Drag model selection** — G2/G5/GI/GS now explicitly alias to G1 (documented).
- **Output** — corrected dope-card Drop MIL sign, `--full` table downrange/drift column swap, dope-card footer weekday (was off by 4 days), stability min-twist units label, and estimate-bc verification range.
- **WASM** — auto-zero solves to line-of-sight (not bore); wind direction consumed as radians; Monte Carlo honors `--drag-model`; full-trajectory JSON keys fixed; wasm32 target now compiles; help text synced to the parser and dead zero flags dropped.
- **Edge cases / divide-by-zero / NaN guards** — sectional-density and BC zero divisors, aero-jump zero/negative twist and `mass_kg<=0`, moist-air speed of sound at zero pressure, air density at non-positive pressure, zero-angle solver reporting non-convergence as success.
- **Robustness** — BC5D/BCCR table loaders use checked-arithmetic dimension bounds and NaN-safe lookups; FFI stops leaking a `CString` per version call and caps oversized/negative Monte Carlo sim counts; fixed panics on PDF UTF-8 header truncation, `percentile(p>1)` OOB, and Euler zero-relative-velocity.
- **Input validation** — CLI rejects non-finite numeric arguments; `target_distance` must be finite and positive (was `> 0`, letting `+inf` through); only `-inf` maps to the ignore-ground sentinel.

### Changed
> **FLAGGED for downstream:** these correct real bugs but **shift numerical output** — altitude, air-density, gravity, and transonic results all move. Re-validate dope and any cached/expected outputs.
- **Air density unified** across all three integrators — `cli_api` no longer double-counts altitude (`exp(-alt/8000)` atop station pressure) and now applies humidity; Monte Carlo no longer reprojects shooter-altitude density to sea level. Shifts *every* trajectory (~0.32% at sea level, larger at altitude).
- **Altitude now drives both station pressure and temperature** — with a default (15 °C / sea-level) input, `--altitude` derives ICAO station pressure *and* applies the −6.5 °C/km lapse rate; an explicit pressure or temperature stays authoritative. Previously altitude-only ran too warm/thin. Validated against py_ballisticcalc to ~0.04% (0–3000 m).
- **Gravity** is the SI-canonical `9.80665` in the default RK4/RK45 kernel (was `9.81`).
- **Transonic single-application** and **Reynolds once + capped** change near-Mach-1 and Monte Carlo drag.
- **Spin twist sign** and **aero-jump Sg** (kg→grains, was ~1000x off) affect opt-in spin/jump diagnostics; left-twist spin-drift direction changes.
- **Default `bullet_length` 4.0→4.5x caliber** (now unified across library/FFI/WASM) — shifts opt-in spin/Magnus/Miller diagnostics only, not the base trajectory.
- **Zero solver returns `Err` on non-convergence** (was a best-effort `Ok(angle)`) — downstream callers must handle `Err`.
- **WASM default drag model G7→G1** to match the G1-scale default BC.
- **Monte Carlo range/velocity means exclude short-falling draws** when no explicit target is given (keeps CEP measurable and FFI arrays equal-length); pass an explicit target for unbiased per-target stats.

### Added
- `atmosphere::resolve_station_pressure(pressure_hpa, altitude_m) -> Option<f64>` — an explicit station pressure stays authoritative; a default pressure with real altitude derives the ICAO station pressure.
- `atmosphere::resolve_station_temperature(temperature_c, altitude_m) -> Option<f64>` — mirrors the above for temperature (returns `None` at the 15 °C default so the ICAO lapse rate applies).

### Performance (output-identical)
- Build `BallisticInputs` once per integration instead of per derivative eval (eliminates per-step String alloc and Vec/DragTable clones on the Monte Carlo hot path).
- Binary-search the drag Mach axis; precompute WindSock segment vectors; hoist RK45 air-density/wind and per-step pitch-damping/drag-coefficient allocations; skip BC-segment re-sort when already sorted; lowercase the bullet model once.

## [0.13.30] - 2026-01-26

### Added
- **Terms of Service Acceptance** - When using `--online` for the first time, users must accept the Terms of Service from https://ballistics.rs/terms.txt
  - TOS is fetched from the server and displayed in the terminal
  - User must type 'y' or 'yes' to accept
  - Acceptance is stored in `~/.ballistics/tos.json`
  - Subsequent runs skip the prompt unless TOS version changes

### Dependencies
- Added `dirs` crate for cross-platform home directory detection

## [0.13.29] - 2026-01-26

### Added
- **Expanded Test Coverage** - Added 36 new unit tests to modules with lighter coverage:
  - `trajectory_integration`: RK4/RK45 consistency, ground impact detection, target distance, wind handling (+6 tests)
  - `fast_trajectory`: Solution interpolation edge cases, BC segment boundaries, event arrays (+6 tests)
  - `stability_advanced`: Bullet type parameters, atmospheric corrections, dynamic stability (+11 tests)
  - `spin_drift_advanced`: Drift direction, edge cases, transonic correction, yaw of repose (+13 tests)

### Changed
- Total test count increased from 156 to 192 tests

## [0.13.28] - 2026-01-26

### Fixed
- **Wind Direction CSV Override** - Fixed default check comparing against 90.0 instead of 0.0 (the actual CLI default), which prevented location CSV wind direction overrides from working

## [0.13.27] - 2026-01-26

### Fixed
- **Location CSV Overrides** - Fixed bug where `final_humidity` and `final_wind_direction` were set but not used, causing location CSV overrides for humidity and wind direction to be ignored
- **Compiler Warnings** - Eliminated all unused variable warnings in main.rs

### Changed
- Location CSV now properly overrides humidity (`HUMIDITY` column) and wind direction (`WIND_DIR` column)

**Related Tickets:**
- MBA-600: Fix unused variable warnings in ballistics-engine CLI

## [0.13.26] - 2026-01-26

### Added
- **Longitude Parameter** - `--longitude` for weather zone features (degrees, -180 to 180)
- **Shot Direction Parameter** - `--shot-direction` for azimuth (0=North, 90=East)

### Notes
- Negative values require equals format: `--longitude=-115.2`
- These parameters are required when using `--enable-weather-zones` or `--enable-3d-weather`

**Related Tickets:**
- MBA-601: Add weather integration control parameter for --online CLI mode

## [0.13.25] - 2026-01-26

### Added
- **Weather Zone Control** - `--enable-weather-zones` flag for automatic weather zone generation
- **3D Weather Control** - `--enable-3d-weather` flag for altitude-dependent atmospheric corrections
- **Wind Shear Model Selection** - `--wind-shear-model` parameter (none, logarithmic, power_law, ekman_spiral)
- **Weather Zone Interpolation** - `--weather-zone-interpolation` parameter (linear, cubic, step)

### Notes
- All weather parameters are passed through to Flask API when using `--online` mode
- Weather features require latitude, longitude, and shot-direction to be specified

**Related Tickets:**
- MBA-601: Add weather integration control parameter for --online CLI mode

## [0.13.24] - 2026-01-26

### Fixed
- **Online Mode Ground Impact** - Fixed `--ignore-ground-impact` flag not being passed to Flask API when using `--online` mode, causing trajectories to be truncated when bullet dropped below -100m threshold

### Added
- `ground_threshold` field in API client TrajectoryRequest

**Related Tickets:**
- Online solver returning partial results investigation

## [0.13.23] - 2026-01-25

### Added
- **BC Truing** - `--bc-adjustment` parameter for multiplying BC by a correction factor (e.g., 0.85)
- **Velocity Truing** - `--velocity-adjustment` parameter for adding offset to base velocity from chronograph data
- **CSV Profile Import** - `--profile` and `--profile-row` to load gun configurations from CSV files
- **CSV Location Import** - `--location` and `--site` to load shooting location environmental data from CSV files
- CSV parser handles Glenn's format (header with `#` prefix, comma-separated values)

### Changed
- `--velocity` and `--bc` are now optional when using `--profile` (loaded from CSV)
- Integration tests now use main `ballistics` binary instead of `ballistics-cli`

**Related Tickets:**
- MBA-591: CLI: Add --bc-adjustment parameter for BC truing
- MBA-592: CLI: Add --velocity-adjustment parameter
- MBA-593: CLI: Add CSV profile import (--profile)
- MBA-594: CLI: Add location presets (--location)

## [0.5.0] - 2025-01-31

### Added
- **Advanced Spin Drift Module** (`spin_drift_advanced.rs`) - Enhanced gyroscopic drift calculations with velocity-dependent effects (MBA-138)
- **Advanced Stability Module** (`stability_advanced.rs`) - Comprehensive stability analysis including transonic behavior (MBA-139)
- Both modules ported from private ballistics_rust repository as part of strategic reconciliation (MBA-144)

### Fixed
- **Critical Wind Vector Calculation** - Fixed swapped sin/cos in `wind.rs` `calc_vec()` function that caused 90° rotation of wind effects (MBA-146)
- **Wind Shear Bugs** - Multiple critical fixes in `wind_shear.rs` for proper altitude-dependent wind modeling (MBA-145)

### Verified
- **Transonic Drag Parity** - Verified parity with private repository implementation, all public functions available (MBA-140)

### Context
This release represents significant progress in the strategic reconciliation between the private ballistics_rust and public ballistics-engine repositories. The strategy established in MBA-144 v0.30.0 makes ballistics-engine the canonical physics source, with improvements flowing from private → public → wrapped in private.

**Module Reconciliation Progress:**
- 7 modules wrapped in ballistics_rust using ballistics-engine
- 2,164 lines of duplicate code eliminated (25.5% reduction)
- ballistics_rust now depends on ballistics-engine as authoritative source

**Related Tickets:**
- MBA-144: Strategic reconciliation with ballistics-engine as canonical source
- MBA-138: Port spin_drift_advanced.rs
- MBA-139: Port stability_advanced.rs
- MBA-140: Verify transonic_drag.rs parity
- MBA-142: Align Rust dependencies (nalgebra 0.33, rayon 1.10)
- MBA-145: Fix critical wind_shear.rs bugs
- MBA-146: Wind vector calculation bug
- MBA-147: Reynolds.rs missing FlowRegime enum

## [0.4.3] - 2025-01-20

### Fixed
- Fixed metric unit calculations
- Updated WASM to use RK45 integration method

## [0.4.0] - 2025-01-15

### Added
- RK45 integration method as default (improved accuracy)
- Muzzle height and target height parameters for realistic trajectory simulation

### Changed
- Major physics improvements across trajectory calculations

## [0.3.0] - 2024-12-01

### Added
- Initial public release
- Core trajectory physics
- RK4 and Euler integration methods
- Monte Carlo simulation support
- WASM support for web deployment

[0.5.0]: https://github.com/ajokela/ballistics-engine/compare/v0.4.3...v0.5.0
[0.4.3]: https://github.com/ajokela/ballistics-engine/compare/v0.4.0...v0.4.3
[0.4.0]: https://github.com/ajokela/ballistics-engine/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/ajokela/ballistics-engine/releases/tag/v0.3.0
