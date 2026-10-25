# Implementation brief: additive wind field modeling in the `wind` module

Status: ready to implement. Tracks [issue #90](https://github.com/ajokela/ballistics-engine/issues/90).
Design record: [WIND_FIELD_MODELING.md](WIND_FIELD_MODELING.md). Data pipeline
(GRIB, radiosonde sets) lives in [winds-engine](https://github.com/ajokela/winds-engine)
and is **not** part of this task.

## Task in one paragraph

Add a stateless, deterministic `WindField` trait to this crate's existing
`wind` module, plus four implementations — `ConstantField`, `ProfileField`
(piecewise-linear vertical profile), `LogProfile`, `EkmanProfile` — and an
adapter so the existing `WindSock` satisfies the trait. **Everything is
additive**: `WindSegment`/`WindSock` and every existing API stay exactly as
they are. No API is renamed, removed, or re-conventioned.

## Non-negotiable constraints

- **Additive only.** Existing wind APIs (`WindSock`, `WindSegment`,
  `wind_vector`, `WindSegmentField`, etc.) must not change signature or
  behaviour. Downstream (rockets-engine, moving-batteries) must compile
  unchanged against the result.
- **Determinism.** No network, no wall clock, no unseeded randomness, no
  `std::time`. Same inputs → same outputs, bit-exact.
- **no_std-compatible shape.** `alloc::vec::Vec` and `nalgebra::Vector3` are
  fine (both already used by the module). No new dependencies; if a new
  dependency seems necessary, stop and flag it before proceeding.
- **Wind conventions must match the existing module** (see
  `src/wind.rs`): world-frame McCoy coordinates (x downrange, y up, z
  lateral), wind-FROM directions with 0° = headwind (from +x, vector
  pointing −x), 90° = from the right (−z). Speeds in the *new* API are
  m/s; `WindSegment` keeps its internal km/h (adapter converts, factor 3.6).
- **Commits as the human developer** (Alex Jokela <alex@camulus.com>), no
  AI/assistant attribution anywhere.

## API to add (in `src/wind.rs` or a new `src/wind/` submodule — your call,
keep the module's public path stable)

```rust
/// A deterministic, stateless wind field: same time and position always
/// give the same wind, so integrator substeps, snapshots and seeded
/// reproduction are exact.
pub trait WindField {
    /// World-frame air velocity (m/s) at a time and position.
    fn wind(&self, time_s: f64, position_m: Vector3<f64>) -> Vector3<f64>;
}

/// Uniform wind everywhere.
pub struct ConstantField(pub Vector3<f64>);

/// One level of a vertical profile.
pub struct ProfileLevel {
    pub altitude_m: f64,
    pub speed_mps: f64,          // nonnegative
    pub direction_from_deg: f64, // wind-FROM, 0° = headwind, matching WindSegment
}

/// Piecewise-linear vertical profile: speed interpolates linearly in
/// altitude, direction along the shortest circular arc; nearest level holds
/// outside the range.
pub struct ProfileField { /* sorted Vec<ProfileLevel> */ }
impl ProfileField {
    pub fn new(levels: Vec<ProfileLevel>) -> Result<Self>; // validate: finite, speeds >= 0, strictly ascending altitudes after sorting
    pub fn levels(&self) -> &[ProfileLevel];
    pub fn wind_at(&self, altitude_m: f64) -> Vector3<f64>;
}

/// Log-law boundary layer: u(z) = u_ref · ln(z/z0) / ln(z_ref/z0), z > z0;
/// zero below z0. Constant direction with height.
pub struct LogProfile { /* reference_speed_mps, reference_height_m, roughness_m, direction_from_deg */ }

/// Ekman spiral (constant eddy viscosity K; Holton's classic solution):
/// u = G(1 − e^(−γz)·cos(γz)), v = sign(f)·G·e^(−γz)·sin(γz),
/// γ = sqrt(|f|/2K), f = 2Ω·sin(latitude). Zero at the surface, approaches
/// the geostrophic wind aloft, rotated to `direction_from_deg`.
pub struct EkmanProfile { /* geostrophic_mps, latitude_deg, eddy_viscosity_m2s, direction_from_deg */ }
```

- `impl WindField for ConstantField`, `ProfileField`, `LogProfile`,
  `EkmanProfile`.
- **`WindSock` adapter**: `impl WindField for WindSock` (or a
  `SegmentedField(WindSock)` newtype if a blanket impl would be surprising)
  returning `vector_for_range_stateless(position_m.x)` — time-independent,
  exactly the existing behaviour.
- Constructors validate and return the crate's `Result`/`Error` types
  (follow the module's existing error style; see `WindSock::try_new`).

## Reference implementation to adapt

winds-engine stages byte-for-byte candidates (MIT, same author):
- `src/field.rs` — `Wind`/`WindField`/`ConstantField`/`ProfileField` with
  the circular direction interpolation already worked out.
- `src/analytic.rs` — `LogProfile`/`EkmanProfile` with endpoint behaviour
  pinned by tests.
- `tests/profile.rs` — 12 tests worth porting (direction round-trips,
  profile interpolation incl. the 350°→10° short-arc case, log-law
  reference/roughness points, Ekman endpoints, validation errors).

Differences to apply while promoting: use `nalgebra::Vector3<f64>` instead
of the staged plain struct (this crate's convention), match this crate's
`Error` type, and use the module's wind-FROM convention note (the staged
`0° = from +X downrange` is the same convention under a different name —
verify with a test that `0°` gives a headwind vector, matching
`wind_vector(speed, 0, 0)`).

## Tests required

- All ported winds-engine tests, adjusted for `Vector3`.
- **Adapter equivalence**: a `WindSock` over several zones returns, via
  `WindField::wind`, exactly `vector_for_range_stateless(x)` at assorted
  positions (including across a zone boundary and at negative x).
- Direction-convention test: `ProfileField` level at 0° produces the same
  vector as `wind_vector(speed, 0.0, 0.0)`; 90° the same as
  `wind_vector(speed, FRAC_PI_2, 0.0)`.
- `WindSegment`-style validation errors (empty levels, unsorted duplicates,
  negative speed, non-finite values).
- No behaviour change: the full existing test suite must pass untouched.

## Docs and release

- Extend `src/wind.rs` module docs with the new types (one short section).
- Add a line to the crate changelog/release notes file if one exists.
- Minor version bump (additive). Publication goes through the full release
  runbook — see the `ballistics-engine-release` skill; do **not** run
  `cargo publish` ad hoc.
- After release: rockets-engine will consume this via
  `BallisticForces::with_wind_field` (separate task, not this one).

## Explicitly out of scope

- GRIB2/HRRR/GFS acquisition, HTTP, file I/O (winds-engine's job).
- Radiosonde datasets and validation (winds-engine's job; it will compare
  its soundings against these profile types).
- Terrain flow / mass-consistent models (deferred entirely).
- Turbulence/gusts (rockets-engine already owns a frozen-field model).
- Any change to existing wind APIs or to the trajectory kernel.
