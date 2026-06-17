# Changelog

All notable changes to the ballistics-engine project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
