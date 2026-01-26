# Changelog

All notable changes to the ballistics-engine project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
