# McCoy Coordinate Migration — Downstream Notes (ballistics_rust)

As of the commit "Migrate engine to McCoy coordinate convention", the engine's
internal coordinate system changed from:

- **OLD:** `X = lateral`, `Y = vertical`, `Z = downrange`  *(left-handed)*

to the McCoy exterior-ballistics convention:

- **NEW (McCoy):** `X = downrange`, `Y = vertical`, `Z = lateral (right)`  *(right-handed)*

## What did NOT change

All **formatted output** is byte-for-byte identical:

- CLI `table` / `json` / `csv` output (columns/keys keep their meaning: `drift`
  is still lateral, range is still downrange, the JSON `x` field is still
  lateral, etc.).
- FFI debug logging (`lateral=`, `downrange=` labels).
- `windage_meters` / `wind_drift` / `drop` semantic fields.

The formatters map at the output boundary, so consumers of *formatted* output
need no changes.

## What DID change — raw vectors are now McCoy-ordered

Any code that reads **raw position/velocity vectors or raw state arrays** now
receives McCoy ordering (`[downrange, vertical, lateral]`). Affected public
surfaces:

| API | Before (per element) | After (McCoy) |
|-----|----------------------|---------------|
| FFI `FFITrajectoryPoint.position_x/_y/_z` | lateral, vert, downrange | **downrange, vert, lateral** |
| FFI `impact_positions_x/_y/_z` | lateral, vert, downrange | **downrange, vert, lateral** |
| `TrajectoryPoint.position` (`.x/.y/.z`) | lateral, vert, downrange | **downrange, vert, lateral** |
| `trajectory_integration::solve_trajectory_rust` dict keys `x`,`z`,`vx`,`vz` | lateral, downrange | **downrange, lateral** |
| `fast_trajectory` state arrays / `FastSolution.y[0]`,`y[2]` | lateral, downrange | **downrange, lateral** |
| `monte_carlo::solve_trajectory_for_monte_carlo` (internal state) | — | now McCoy |
| `derivatives::compute_derivatives` state in/out | lateral, downrange | **downrange, lateral** |

`Coriolis` note: migrating left-handed → right-handed reflects the frame, so the
cross-product term flipped sign internally (`-2Ω×v` → `+2Ω×v`). Output deflection
is unchanged; only relevant if you reimplement Coriolis against engine vectors.

## ballistics_rust changes to make (when bumping the engine dependency)

`ballistics_rust` currently pins `ballistics-engine` at `v0.13.38` (pre-McCoy),
so **nothing breaks until the tag is bumped**. When updating to the McCoy
version, apply these:

### 1. `ballistics_rust/src/lib.rs` — local integrator (~line 1576)
It builds `pos = Vector3::new(initial_state[0], initial_state[1], initial_state[2])`
and calls `ballistics_engine::derivatives::compute_derivatives`. After the bump,
`compute_derivatives` treats `state[0]` as downrange. Ensure the `initial_state`
this integrator receives is McCoy-ordered (downrange in `[0]`, velocity downrange
in `[3]`), and that `positions_x/_z` / `velocities_x/_z` exports are interpreted
as McCoy by callers. The existing `pos.x > initial_state[0] + 2000.0` range guard
is already correct for downrange-in-x.

### 2. Python: `ballistics/solvers/trajectory_helpers_refactored.py`
This reads raw state from `solution.y`:
- `final_x = final_state[0]` is now **downrange** (was lateral)
- `final_z = final_state[2]` is now **lateral** (was downrange)

Swap the downstream uses:
- `wind_drift_m = final_x`  →  `wind_drift_m = final_z`  *(lateral)*
- `apex_distance = final_z` and `solution.y[2][i]`  →  use `final_x` / `solution.y[0][i]`  *(downrange)*
- `expected_y = -sight_height_m + final_z * tan(angle)`  →  use `final_x` (downrange) in the LOS term

### 3. Python: `ballistics/solvers/integrator.py`
Debug logs and angle math read `state[2]` as downrange and `vz_downrange`
(index 5). If this path consumes engine/Rust state (vs. a pure-Python scipy
integrator — verify), update `state[2]` → `state[0]` for downrange and use
velocity index `[3]` for downrange velocity.

### 4. Verify, don't assume
`ballistics_rust` and the Python layer have their own integrators and tests.
Confirm which result-consuming paths are fed by the **Rust/engine** solver
(now McCoy) vs. a **pure-Python** solver (unchanged) before swapping indices,
then run that repo's test suite.

## How the engine change was verified

- CLI output: golden-master diff across formats (table/json/csv) × units
  (imperial/metric) × effects (spin drift, Coriolis @ 2 azimuths, Magnus, wind
  @ 2 directions, zero) — **byte-identical**.
- FFI: `tests/test_ffi.c` output — **byte-identical** (exercises wind shear +
  sampling).
- Secondary solver (`fast_trajectory`/`monte_carlo`): direct HEAD-vs-branch
  comparison of `solve_trajectory_for_monte_carlo` with Coriolis+Magnus+spin
  enabled — **byte-identical**.
- Full `cargo test` suite: 221 passing.

---

# SI Unit Unification — Downstream Notes (ballistics_rust)

A follow-up release makes `bullet_mass`/`bullet_diameter`/`bullet_length`
**SI-canonical (kg, meters, meters)** across the *entire* engine, including the
secondary solver path that `ballistics_rust` drives. Previously that path read
those fields as **grains/inches** (a hidden imperial sub-convention). The
imperial mirror fields `caliber_inches` (inches) and `weight_grains` (grains)
remain and must be kept in sync.

## Breaking changes for ballistics_rust

| API (engine) | Old contract | New contract |
|---|---|---|
| `derivatives::compute_derivatives` (driven directly in `ballistics_rust/src/lib.rs`) | `inputs.bullet_mass` grains, `bullet_diameter`/`bullet_length` inches | **kg, meters, meters**; also set `caliber_inches`/`weight_grains` (used by stability/Magnus/drag-shape helpers) |
| `fast_trajectory::fast_integrate` | `inputs.bullet_mass` grains | **kg** |
| `monte_carlo::solve_trajectory_for_monte_carlo` | inputs in **yards / fps / grains / feet / inches / km·h / degrees** | **meters / m·s / kg / meters(alt) / meters(sight) / m·s(wind) / radians(wind & azimuth)** |

### What to change
1. Wherever `ballistics_rust` builds a `BallisticInputs` for these calls, populate
   `bullet_mass` in kg, `bullet_diameter`/`bullet_length` in meters, and also set
   `caliber_inches`/`weight_grains` (the engine no longer derives them on this path).
2. The Monte Carlo wrapper (`ballistics_rust/src/monte_carlo.rs`) historically fed
   `solve_trajectory_for_monte_carlo` legacy imperial values — convert those to SI
   before the call (target_distance m, muzzle_velocity m/s, mass kg, altitude m,
   sight_height m, wind_speed m/s, wind_angle radians).

## Bug fixes in this release (output changes even for correct callers)

- **Enhanced spin drift**: a `×15.432358` (grams→grains) factor wrongly inflated
  the grains mass, making spin decay ~3.93× too slow. Spin-drift magnitudes change
  (grows with range).
- **cluster_bc** (`use_cluster_bc=true`): caliber was passed in meters where
  inches were expected, misclassifying every bullet. BC degradation changes.
- **Monte Carlo trajectory**: `atmo_params` were packed in the wrong order and
  `base_ratio` was set from humidity, scrambling base density to ~417 kg/m³
  (~340× drag) — the Monte Carlo bullet previously collapsed at ~33 m. Now
  `solve_trajectory_for_monte_carlo` agrees with the cli_api solver within ~0.3%.
  (`atmo_params` itself is gone since 0.50.0; see the next section.)
- **Reynolds correction** (subsonic tail, mach<1 & <200 m/s): diameter was passed
  in inches where meters were expected.

---

# Humidity Consistency (0.50.0) — Downstream Notes (bindings, rockets-engine, embedders)

0.50.0 flies every solver family on the same humid air, through one per-step lookup,
`atmosphere::air_at`. The kernels used to take a bare `(f64, f64, f64, f64)` atmosphere tuple
with no humidity in it and flew a 0 % sound speed whatever the station's humidity. Their
atmosphere is now a `KernelAtmosphere`, which carries the humidity, so the density and the
sound speed cannot disagree about it.

## Breaking changes (Rust source)

| 0.49.0 | 0.50.0 |
|---|---|
| `TrajectoryParams.atmos_params: (alt, T, P, ratio)` | `TrajectoryParams.atmosphere: KernelAtmosphere` |
| `FastIntegrationParams.atmo_params: (alt, T, P, ratio)` | `FastIntegrationParams.atmosphere: KernelAtmosphere` |
| `compute_derivatives(.., atmos_params: (f64, f64, f64, f64), ..)` | `compute_derivatives(.., atmosphere: &KernelAtmosphere, ..)` |
| `try_solve_trajectory_rust` / `solve_trajectory_rust` atmosphere tuple | `atmosphere: KernelAtmosphere` |
| `aerodynamic_jump_launch_offset_rad(inputs, (f64, f64, f64, f64))` | `aerodynamic_jump_launch_offset_rad(inputs, &KernelAtmosphere)` |
| `try_integrate_trajectory` -> `Result<_, WindSegmentError>` | `Result<_, KernelInputError>` (a wind segment or an atmosphere zone) |
| `BallisticInputs.humidity`, a 0-1 FRACTION, and `humidity_percent()` | `BallisticInputs.humidity_percent`, a PERCENT (default 50); `humidity_percent()` is deprecated for one release and returns the field (0.49.0's clamped into 0-100; 0.50.0 clamps nowhere) |
| `atmosphere::get_local_atmosphere` | deprecated (dry by contract); use `air_at(&KernelAtmosphere::Station(StationAir::new(..)?), None, 0.0, altitude_m)`, or `StationAir::with_density_ratio` for a caller-held density ratio |
| `wez::compute_wez(..24 arguments..)`, which always flew 50 % | `wez::compute_wez_with_humidity(.., humidity: RelativeHumidity)`; the 24-argument `compute_wez` is deprecated for one release and flies `RelativeHumidity::DEFAULT` |
| an exhaustive `match` on `wind_scenarios::WindScenarioError` | two new variants: `InvalidHumidity(AtmosphereError)` (a load's `humidity_pct` beyond the sensor band, NaN or infinite, at the field `humidity_pct`; 0.49.0 reported a NaN as `InvalidLoad { field: "humidity", .. }` and flew 150 as 100 %) and `ZeroFailed { message }` (a zero that cannot be solved, not a `SolveFailed` of a scenario named `<zero>`) |
| `wasm::Calculator::set_humidity(..) -> Calculator` (wasm32) | `-> Result<Calculator, JsValue>`: it refuses a humidity beyond the sensor band itself |
| a legacy tuple `compute_derivatives` flew: a negative or infinite direct-mode value, a NaN, -inf or >= 2 ratio, a negative pressure, a temperature at or below absolute zero or above water's critical point, air that cannot exist at the stated humidity | refused by `from_legacy_params` with a named `AtmosphereError` (see below); `(0, c, 0, 0)` is still an exact vacuum |
| `tolerance::ToleranceAxisV1 { .. }` literal | names `near_at_physical_boundary` and `far_at_physical_boundary` too |

`compute_derivatives` and `fast_integrate*` now measure a point's height from the
`KernelAtmosphere` station's altitude, not from `inputs.altitude`, so `pos.y` is a height above
the station. A caller that passes an absolute altitude in `pos.y` (with `inputs.altitude` at 0)
subtracts the station's altitude, `StationAir::altitude_m()`, from it; then a dry flight's
acceleration is 0.49.0's to the bit. `air_at` itself takes the absolute altitude:

```rust
// rockets-engine's shape: an absolute altitude in the state, the pad as the station.
let pos = Vector3::new(x, altitude_m - station.altitude_m(), z);
let kernel = KernelAtmosphere::Station(station);
let accel = compute_derivatives(pos, vel, &inputs, wind, &kernel, bc, None, t, None);
let air_here = air_at(&kernel, None, x, altitude_m);
```

## What to change

```rust
use ballistics_engine::{KernelAtmosphere, RelativeHumidity, StationAir};

// 0.49.0
//   atmos_params: (altitude_m, temperature_c, pressure_hpa, density_ratio),
// 0.50.0: the station itself, with its humidity (density and sound speed both follow it).
let humidity = RelativeHumidity::try_percent(humidity_percent)?;
let atmosphere = KernelAtmosphere::Station(StationAir::new(
    altitude_m,
    temperature_c,
    pressure_hpa,
    humidity,
)?);
// A caller that must keep its own density ratio (a legacy 4-vector) converts the tuple, at
// the humidity the ratio was built at; direct mode's (density, sound speed, 0, 0) becomes
// KernelAtmosphere::uniform(density, sound_speed) the same way (a density of 0 is an exact
// vacuum: gravity and no drag).
let legacy = KernelAtmosphere::from_legacy_params(
    (altitude_m, temperature_c, pressure_hpa, density_ratio),
    humidity,
)?;
```

- `BallisticInputs.humidity = rh / 100.0` becomes `humidity_percent = rh`: delete the `/ 100`.
  A fraction from a wire converts with `RelativeHumidity::try_fraction(f)?.percent()`.
- Every constructor validates: a humidity more than 5 points outside 0-100 %
  (`HUMIDITY_SENSOR_TOLERANCE_PERCENT`), NaN, or air that cannot exist is an
  `AtmosphereError`; a reading inside the band is stored as exactly 100 or 0.
- The crate root now exports `RelativeHumidity` (with `StationAir`, `KernelAtmosphere`,
  `LocalAir`, `air_at`, `AtmosphereError`, `AtmosphereErrorReason`, `DEFAULT_HUMIDITY_PERCENT`
  and `HUMIDITY_SENSOR_TOLERANCE_PERCENT`). `InputAxis` has a `RelativeHumidity` variant, so
  glob-importing both `ballistics_engine::*` and `InputAxis::*` makes the name ambiguous; name
  one of them explicitly.

## Legacy tuples that no longer convert

`from_legacy_params` validates what 0.49.0's kernels flew unchecked.
`tests/kernel_atmosphere_legacy_tuples.rs` lists every tuple class 0.49.0's
`compute_derivatives` flew, with what 0.50.0 does with it:

| 0.49.0 tuple | 0.49.0 flew it as | 0.50.0 |
|---|---|---|
| `(0, c, 0, 0)` or `(-0, c, 0, 0)` | an exact vacuum | the same vacuum, `KernelAtmosphere::uniform(0.0, c)` |
| `(rho < 0, c, 0, 0)` | a forward drag | `UniformDensityOutOfRange` on `density_kg_m3` |
| `(rho, inf, 0, 0)` | Mach 0 | `NotFinite` on `speed_of_sound_mps` |
| `(alt, t, p, NaN or -inf)` | ratio 1.0 | `NotFinite` on `density_ratio` |
| `(alt, t, p, ratio >= 2)` | that ratio (`fast_integrate` refused it) | `DensityRatioOutOfRange` on `density_ratio` |
| `(alt, t, p < 0, ratio)` | its ratio, ignoring the pressure | `NonPositivePressure` on `pressure_hpa` |
| `(alt, t <= -273.15, p, ratio)` | a finite flight | `AtOrBelowAbsoluteZero` on `temperature_c` |
| `(alt, t > 373.946, p, ratio)` | pure steam's sound speed | `ImpossibleAir` on `temperature_c`, at any humidity |
| any station whose humidity needs more vapour pressure than `p` | dry air (the kernel never saw the humidity) | `ImpossibleAir` on `humidity_percent` |

A finite ratio of 0 or less still means "not supplied" and flies 1.0. Tuples 0.49.0 flew to a
NaN acceleration are each a named refusal now.

## Numeric effect

Two changes move numbers at a humidity above 0; with the station and every zone at 0 % nothing
moves, bit for bit.

- **Kernel runs.** Their sound speed rises to the station's own moist value (+0.133 % at
  15 °C / 50 %, +0.787 % at 35 °C / 90 %), and inside an `AtmoSock` zone to the zone's own, so
  Mach falls and the drag coefficient moves along the drag curve. The station density is
  unchanged.
- **Every path, `TrajectorySolver` included: the air aloft.** Away from the station's height
  the air now holds the station's water-vapour mole fraction for its sound speed, its density
  and its (virtual-temperature) pressure lapse. Flat fire moves by under 0.1 mm at 1000 m; a
  30° lob moves by tenths of a metre to metres.

CHANGELOG.md's 0.50.0 section has the measured shifts of both.

To reproduce a 0.49.0 kernel result exactly (for re-baselining, not for new work), pass the
same humid ratio at `RelativeHumidity::DRY`:
`KernelAtmosphere::from_legacy_params((alt, t, p, ratio), RelativeHumidity::DRY)` flies
0.49.0's half-humid air (the humid density with a dry sound speed). It was checked bit for bit
against the 0.49.0 build on `try_integrate_trajectory` and `fast_integrate`, for a flight with
no `AtmoSock` zones or only dry ones. A zone with a humidity above 0 now flies its own moist
sound speed and lapse, at its own height too, and no API reproduces 0.49.0 there; migrate a
zoned caller with `StationAir::new` at the station's humidity, so the station and its zones fly
one model.
