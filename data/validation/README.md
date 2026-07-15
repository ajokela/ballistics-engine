# Golden-reference physics validation dataset (MBA-1318)

A versioned, provenance-tagged set of reference trajectories with per-source tolerance bands.
Each `*.json` file in this directory is one **case**: a complete set of solver inputs plus a list
of expected observables (drop / drift / time-of-flight / velocity at specific downrange stations)
with the tolerance and the *reason for that tolerance* attached to every single expectation.

The cases are replayed through the real `TrajectorySolver` by `tests/golden_physics.rs`, gated
behind the non-default `validation` cargo feature. The goal is that a physics regression shows up
as a red CI check with a rich, aligned delta table — not as a surprised email from a field tester.

```bash
# run the whole harness (all cases must be green at their tolerances)
cargo test --features validation --test golden_physics

# print the delta table for EVERY expectation (pass or fail) to see the head-room
GOLDEN_DUMP=1 cargo test --features validation --test golden_physics -- --nocapture
```

The harness adds **zero** dependencies: `serde` / `serde_json` are already core dependencies.

---

## Case-file schema

One case per file. Fields:

```jsonc
{
  "id": "string",                 // unique, matches the file name (sans .json)
  "description": "string",        // human summary of what the case exercises
  "source": {
    "kind": "analytic" | "cross-implementation" | "published",
    "citation": "string",         // exact provenance (formula / tool+version / publication)
    "retrieved": "YYYY-MM-DD",    // optional; when the reference was produced/obtained
    "generator_version": "string" // optional; e.g. the pinned generator + tool version
  },
  "inputs": { /* see below */ },
  "expectations": [
    {
      "observable": "drop_m" | "drift_m" | "tof_s" | "velocity_mps",
      "range_m": 900.0,           // downrange distance (McCoy X) at which to evaluate
      "value": -11.23038,         // the reference value
      "tol_abs": 0.01,            // optional absolute tolerance (meters / seconds / m·s⁻¹)
      "tol_rel": 0.008,           // optional relative tolerance (fraction of |value|)
      "tolerance_justification": "why this band, physically"
    }
  ],
  "known_discrepancy": false,     // optional; see "Known discrepancies" below
  "tracking_note": "string"       // required iff known_discrepancy is true
}
```

### Observables and the reference frame

The solver integrates in the McCoy frame: **X = downrange, Y = vertical (up positive), Z = lateral
(right positive)**. The origin is the muzzle. Observables are obtained by interpolating the solved
trajectory points at the requested `range_m`:

| observable     | meaning                                                                     |
|----------------|-----------------------------------------------------------------------------|
| `drop_m`       | signed vertical position (McCoy Y). **Up is positive**, so a real bullet drop below the bore line is *negative*. For a level (0°) shot this is the classic drop. |
| `drift_m`      | lateral position (McCoy Z). Right (for a right-hand twist / from-the-left wind) is positive. |
| `tof_s`        | time of flight to `range_m`.                                                |
| `velocity_mps` | speed (velocity magnitude) at `range_m`.                                    |

### Tolerance semantics

A single expectation passes when

```
|actual − value| ≤ tol_abs + tol_rel·|value|
```

(the NumPy `allclose` convention — `tol_abs` is the floor near zero, `tol_rel` scales with
magnitude). **At least one** of `tol_abs` / `tol_rel` must be present, and every expectation must
carry a non-empty `tolerance_justification`; the harness asserts both.

### Inputs

Inputs are **engine-native** (SI, engine conventions). Anything omitted takes the documented engine
default. Fields consumed by the harness (`build_inputs` / `build_wind` / `build_atmosphere` in
`tests/golden_physics.rs`):

| field                     | unit / meaning                                                        |
|---------------------------|-----------------------------------------------------------------------|
| `muzzle_velocity_mps`     | m/s                                                                    |
| `bullet_mass_kg`          | kg                                                                     |
| `bullet_diameter_m`       | m (imperial mirrors `caliber_inches` / `weight_grains` are derived)   |
| `bullet_length_m`         | m (matters only for the spin-drift / stability path)                  |
| `bc_value`, `bc_type`     | ballistic coefficient and drag model (`"G1"`, `"G7"`, …)              |
| `muzzle_angle_rad`        | launch elevation (0 = level bore)                                     |
| `shooting_angle_rad`      | uphill (+) / downhill (−) incline of the shot frame                   |
| `twist_rate_in`           | barrel twist, inches/turn (0 = no spin drift)                         |
| `is_twist_right`          | right-hand twist when true                                            |
| `sight_height_m`, `muzzle_height_m` | m                                                           |
| `altitude_m`, `temperature_c`, `pressure_hpa`, `humidity_frac` | station atmosphere (humidity is a 0–1 fraction) |
| `wind_speed_mps`, `wind_angle_rad`, `vertical_wind_mps` | wind (`wind_angle` 0 = headwind, π/2 = from the right, π = tail, 3π/2 = from the left) |
| `use_rk4`, `use_adaptive_rk45` | integrator selection (default RK45 adaptive — what ships)        |
| `use_enhanced_spin_drift`, `enable_advanced_effects`, `enable_magnus`, `enable_coriolis` | effect flags (default off) |
| `latitude_deg`, `shot_azimuth_rad` | for Coriolis                                             |
| `zero_drag`               | **harness directive**: install an all-zero custom drag table so gravity is the only force (exact vacuum trajectory). Used by the analytic anchors. |
| `ground_threshold_m`      | optional override; defaults to a very negative value so the shot always reaches the furthest observation station. |

---

## Source kinds and provenance

### `analytic` — exact truths (unimpeachable anchors)

Closed-form references with (near) machine-precision tolerances. They pin the integrator + gravity
frame independently of any external tool. With `zero_drag` the only force is constant gravity
(`g = 9.80665 m/s²`, `src/constants.rs`), and RK4/RK45 integrate the resulting quadratic **exactly**;
the only residual is the harness's linear interpolation of a parabola between stored points
(~0.1 mm — see the head-room table below). Cases:

* `analytic_vacuum_flat_drop` — level vacuum parabola: drop / tof / speed.
* `analytic_vacuum_angled_apex_range` — 35° launch: apex height, return-to-launch range, both tofs.
* `analytic_vacuum_uphill` / `analytic_vacuum_downhill` — shot-frame gravity decomposition on a ±20° incline.

### `cross-implementation` — py-ballisticcalc

Reference trajectories from **py-ballisticcalc 2.2.10**, an independent open-source point-mass
solver (its default **RK4 fixed-step** engine, standard G1/G7 drag tables, Miller stability + Litz
spin drift). Generated by `generators/gen_pybc_refs.py` (committed, reproducible). The two frames
are made to coincide exactly by firing flat (0° barrel elevation, no zeroing) with sight height 0
and look angle 0, so the line of sight is the horizontal bore line; then py-bc `height` == engine
`drop_m`, `windage` == `drift_m`, `velocity` == `velocity_mps`, `time` == `tof_s`. Wind is mapped
analytically (`pybc direction_from = engine wind_angle + 180°`, derived from the two wind-vector
definitions). Humidity is pinned to 0 on both sides.

Agreement is bounded by **solver-class differences** — adaptive RK45 (engine) vs fixed-step RK4
(py-bc), and the two solvers' distinct G1/G7 table resolutions/interpolants — not by machine
precision. Tolerances are set at roughly **15–25× the observed cross-solver delta**: loose enough to
absorb those legitimate differences, tight enough to fail on a field-noticeable regression (a ~1%
long-range drop shift ≈ ~1 MOA). They are **not** reverse-fitted to the pass/fail boundary; each
`tolerance_justification` states the mechanism. Cases span the regime matrix:

| case                        | drag | regime               | feature exercised            |
|-----------------------------|------|----------------------|------------------------------|
| `xref_g1_supersonic_308`    | G1   | supersonic           | drop/tof/vel; lateral ≈ 0    |
| `xref_g7_supersonic_308`    | G7   | supersonic           | drop/tof/vel                 |
| `xref_g7_65mm_supersonic`   | G7   | supersonic (6.5 mm)  | second G7 caliber            |
| `xref_g7_transonic_reach`   | G7   | through transonic    | integrates the Mach 1.2→1.0 drag rise |
| `xref_g1_crosswind`         | G1   | supersonic           | crosswind drift, from-right  |
| `xref_g7_crosswind`         | G7   | supersonic           | crosswind drift, from-left   |
| `xref_g7_spin_drift_on`     | G7   | supersonic           | Litz+Miller spin drift ON    |
| `xref_g7_spin_drift_off`    | G7   | supersonic           | spin-drift flag gated OFF (lateral ≈ 0 despite a twist) |

To regenerate:

```bash
python3 -m venv venv && venv/bin/pip install "py-ballisticcalc==2.2.10"
venv/bin/python data/validation/generators/gen_pybc_refs.py
```

### `published` — **TODO (not yet shipped)**

The ticket's ambition was 2–3 externally *published* reference points (published G7 drop tables,
manufacturer trajectory tables, or a textbook worked example) with real citations. This was **not
met for the initial dataset**, and no `published`-kind case is shipped, deliberately — provenance
was not fabricated. Attempts made (2026-07-15): Sierra's online .308 175 gr MatchKing chart (behind
an age-gated SPA — no table served), R. McCoy *Modern Exterior Ballistics* on the Internet Archive
(full-text fetch truncated to the front matter), a Gorilla Ammunition manufacturer data sheet (a
marketing PDF exposing only muzzle velocity, no trajectory table), and several web searches (only
approximate MOA dope, never an exactly-citable table with full stated conditions).

**Slot spec for a future published case.** Prefer `velocity_mps` (and `tof_s`) as the observables:
retained velocity/time at range are essentially independent of the zero / sight height / small
launch angle, so a published *zeroed* table still maps cleanly onto this flat-fire harness without
replicating the zeroing convention. Record the source's muzzle velocity, BC (model + value), and
atmosphere verbatim in `source.citation`, set `kind: "published"`, and justify each tolerance by the
source's stated precision plus the solver-class gap. `drop_m` from a published table requires also
reproducing that table's zero range and sight height and is discouraged for a first published case.

---

## Known discrepancies

If a case genuinely disagrees with the engine and investigation shows the *engine* is the outlier
(the reference is sound and the tolerance is justified), keep the case, set `known_discrepancy: true`
with a `tracking_note`, and open a ticket. The harness then reports that case in a separate
informational section and **excludes it from pass/fail**, so a real, tracked discrepancy stays
visible without wedging CI. There are **no** known discrepancies in the initial dataset — every
analytic and cross-implementation expectation passes well inside its band.

## Observed head-room (initial dataset, engine 0.24.1)

Every expectation currently sits at `delta/tol ≤ 0.07`. Representative worst cases:

| case                      | observable @ range  | |delta|        | tol     | delta/tol |
|---------------------------|---------------------|----------------|---------|-----------|
| analytic (all)            | drop / tof / vel    | ≤ 0.11 mm      | ≥ 3 mm  | ≤ 0.04    |
| xref supersonic           | drop @ 900 m        | ~4 mm / 11 m   | ~0.1 m  | ~0.04     |
| xref supersonic           | velocity @ 900 m    | ~0.19 m/s      | ~2.6 m/s| ~0.07     |
| xref crosswind            | drift @ 800 m       | ~2.4 mm / 3.6 m| ~0.09 m | ~0.03     |
| xref spin-drift ON        | drift @ 900 m       | ~0.16 mm / 0.31 m | ~0.016 m | ~0.01  |
| xref transonic reach      | drop @ 1400 m       | ~1.6 cm / 42 m | ~0.54 m | ~0.03     |

The engine and py-ballisticcalc agree to roughly **0.03–0.05%** across the whole matrix — an
independent confirmation of the engine's supersonic/transonic drag, wind, and spin-drift models.
