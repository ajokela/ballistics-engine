# Wind field modeling in the `wind` module (planned, additive)

Status: decided 2026-10-24; implemented 2026-10-25 per
[WIND_FIELD_IMPLEMENTATION.md](WIND_FIELD_IMPLEMENTATION.md), unreleased. See "As built" below.
Companion: [winds-engine plan](https://github.com/ajokela/winds-engine/blob/main/docs/winds-engine-plan.md)
(GRIB acquisition and radiosonde validation live there; this note covers only
what lands in this crate).

## Decision

Wind **modeling** is upgraded inside this crate's existing `wind` module —
**additively**. `WindSock`/`WindSegment` and every existing API stay exactly
as they are; nothing is replaced and no API expectations change. The
acquisition/validation pipeline (GRIB2/HRRR/GFS ingestion, sounding
datasets) lives in the separate `winds-engine` crate, scoped to rocket
modeling, so this crate's release train is not coupled to data-pipeline
iteration.

## What lands here (one minor release)

All pure additions to the `wind` module:

- **`WindField` trait** — `wind(time_s, position_m) -> Vec3`, stateless and
  deterministic. Same sampling pattern downstream consumers already rely on
  for RK4 substeps, snapshots and seed reproduction (rockets-engine's
  frozen-field turbulence established this pattern).
- **`ConstantField`** — today's uniform vector as a field.
- **`ProfileField`** — piecewise-linear vertical profile
  `Vec<(altitude_m, speed_mps, direction_deg)>`.
- **Analytic profiles** — `LogProfile` (roughness boundary layer) and
  `EkmanProfile` (geostrophic wind, latitude, roughness), consistent with
  this crate's atmosphere tables.
- **`WindSock` adapter** — so the existing segmented wind satisfies
  `WindField`; old and new interop with no consumer changes.

## What does not land here

- GRIB2/HRRR/GFS acquisition and decoding, HTTP subset downloads, offline
  caches, radiosonde validation datasets — winds-engine's job.
- Diagnostic mass-consistent terrain flow — deferred; revisit only with a
  concrete artillery (moving-batteries) requirement. If ever built, its
  home is this crate, not winds-engine.
- Any change to existing wind APIs.

## Determinism policy

Unchanged: no network, no wall-clock, no unseeded randomness in the crate;
everything reproducible from inputs.

## Validation

- Unit tests: profile interpolation, analytic-profile sanity (log law
  reproduces its reference point; Ekman endpoints and turning sense; see "As
  built" on the log-layer limit), adapter equivalence with `WindSock` on the
  segment ladders.
- Cross-check: winds-engine's radiosonde datasets compare `ProfileField`
  providers against sounding archives on a handful of station-days.

## Milestone

Kernel `WindField` minor release, then rockets-engine
`with_wind_field`/`--wind-profile` end to end against two sounding days.

## As built

- `EkmanProfile` follows the implementation brief: geostrophic speed, latitude,
  a constant eddy viscosity K (not a roughness length) and the geostrophic
  direction. The constant-K solution is linear near the ground, not
  logarithmic, so it does not reduce to the log layer; `LogProfile` is the
  surface-layer model.
- Its turning sense is Holton's: in the northern hemisphere the surface wind
  is backed about 45 degrees from the geostrophic wind and veers with height.
  The staged winds-engine version rotated the other way; a test pins both
  hemispheres here.
- `LogProfile` and `EkmanProfile` take y as height above the ground.
  `ProfileField` uses its levels' datum. `HeightOffset` is the one datum
  conversion: its offset is the height of the caller's y = 0 in the field's
  datum. rockets-engine's y is altitude above sea level (its CLI's ground sits
  at `--ground-altitude-m`, sea level by default), so a sounding profile in
  metres above sea level needs no conversion there, while a log or Ekman
  profile takes minus the ground elevation. This crate's `TrajectorySolver`
  measures y from the level `muzzle_height` is measured from and takes the
  ground to run along the line of sight, as its ground check and built-in
  shear do, so a log or Ekman profile there takes the bore's real height
  above the ground minus `muzzle_height` on any shot. Sea level stays level:
  a sea-level sounding takes the atmosphere's `altitude` on a level shot and
  the solver's own air altitude (`altitude + x·sin θ + y·cos θ`) on an
  inclined one. Convert in one place only.
- winds-engine's `validate()` takes its own staged `ProfileField`, so the
  radiosonde cross-check above reaches this crate's `ProfileField` only once
  winds-engine moves to these types, which changes winds-engine's API.
