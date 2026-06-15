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
