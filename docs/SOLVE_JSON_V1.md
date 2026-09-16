# Solve JSON v1 contract

Solve JSON is the binding-neutral request and result contract for one deterministic trajectory
solve. It is separate from `BallisticInputs`, command-line profile JSON, FFI structs, and any
language binding. Those implementation APIs may change without changing this contract.

The Rust DTOs and the envelope-producing decoder are public in
`ballistics_engine::solve_json`. The transport-free Rust service is
`ballistics_engine::solve_v1`; it resolves defaults, validates physical and cross-field rules,
runs the engine, and constructs either a success value or a structured error envelope. Input and
output transport remain outside this contract.

## Process transport

The additive CLI transport reads one v1 request from standard input and writes one compact JSON
envelope followed by a newline to standard output:

```text
ballistics solve-json < request.json > response.json
```

Input is limited to 1 MiB (1,048,576 bytes); the exact limit is accepted. The command does not
read profiles, files, or the network. All dimensional fields remain explicit SI even when the
global `--units` option is supplied. Once command-line parsing has selected `solve-json`, stdout is
reserved exclusively for the protocol envelope; handled failures use the same v1 error shape as
the library service.

| Exit status | Meaning |
| --- | --- |
| `0` | Success envelope. |
| `1` | Standard-I/O or internal failure. |
| `2` | Malformed JSON, schema, shape, or semantic validation failure. |
| `3` | Resource limit or engine solve failure. |

Malformed JSON includes one-based `line` and `column` coordinates. A contained panic becomes a
generic `internal_error`; panic payloads and backtraces are never placed in the JSON envelope.
Failures before successful command selection (for example, an unknown command-line option) remain
ordinary command-line errors rather than protocol responses.

## Versioning and compatibility

Every request and envelope contains the integer:

```json
"schema_version": 1
```

Version 1 accepts only the value `1`; any other integer returns
`unsupported_schema_version`. Version dispatch happens before v1 field validation, so a v2
document containing fields unknown to v1 still receives `unsupported_schema_version`, not
`unknown_field`. The public Rust `SchemaVersionV1` invariant type always serializes as integer `1`
and deserializes only that integer.

All objects reject unknown fields. A producer must therefore emit only fields documented for v1,
and a consumer must not silently reinterpret a v1 field. Removing a field, changing a unit or
sign, renaming an enum value, or changing a field's meaning requires a new schema version — these
are the invariants v1 guarantees and they hold absolutely.

v1 does, however, grow by ADDITION under three strict rules, so that a feature need not force a
whole new schema:

1. A new REQUEST field must be optional and default to the exact pre-existing behavior when
   omitted. Every request valid before the field remains valid and produces identical results
   (`zero_poi_up_m`, `zero_poi_right_m`, `sight_offset_lateral_m`, `drops_reference`,
   `wind_reference` were added this way).
2. A new RESPONSE field must be omitted whenever its feature is inactive, so a response that does
   not exercise the feature is byte-identical to one produced before the field existed
   (`reticle_hold` appears only when the request carries a `reticle` block).
3. An existing REQUEST enum may gain a new accepted value when it names newly supported behavior.
   The new value is opt-in and may be echoed in the resolved response; requests using older values
   keep identical behavior and output. Producers that must work with older engines must not send
   the new value (`G2`, `G5`, `GI`, `GS`, and `RA4` drag models were added this way in MBA-1442).

Consumers must therefore tolerate response fields they do not recognize (ignore, do not reject) if
they want to parse across engine versions; a consumer that pins a response DTO with
`deny_unknown_fields` is pinning an engine version, not the v1 contract. **One documented
exception to rule 2:** `summary.equivalent_horizontal_range_m` (MBA-1395) appears on any inclined
shot (`shooting_angle_rad != 0`) without an explicit opt-in, because an incline-corrected shoot-to
range is meaningful exactly when a look angle is set. A strict consumer that must reject unknown
fields will see this on inclined-shot responses; flat-fire responses are unaffected.

Object member order, indentation, insignificant whitespace, human-readable diagnostic messages,
the engine version string, and the exact decimal spelling of floating-point values are not wire
contracts. Checked regression fixtures use
`abs(actual - expected) <= 1e-10 + 1e-9 * max(abs(actual), abs(expected))` to detect accidental
numeric drift. An intentional physics correction may update those expected values with its
regression evidence without requiring v2, provided the documented field meaning is unchanged.

Representative request, success, error, resource-limit, and early-termination documents live in
[`tests/fixtures/solve_json_v1`](../tests/fixtures/solve_json_v1/). Tests compare parsed object
shapes and tolerant numeric values rather than serialized bytes. A minimal Lattice process
consumer is checked in as
[`examples/solve_json_v1_lattice.lat`](../examples/solve_json_v1_lattice.lat); it intentionally
uses only a fixed trusted fixture until Lattice provides its shell-free argv/stdin process API.

Fields documented with defaults may be omitted inside their section. The request DTO preserves
that omission: serializing a decoded request does not insert defaults, and an explicitly supplied
value remains present even when it equals the documented default. This lets the solve service
distinguish caller intent while it performs one centralized resolution step. The eight top-level
request sections themselves are required so the shape remains explicit.

Omission means the member is not present. Explicit JSON `null` is invalid for every request field,
including fields that may be omitted; the envelope-producing decoder reports `invalid_value` at
that field's path rather than treating null as omission.

Success responses use a distinct `ResolvedSolveRequestV1` representation for
`resolved_request`. Every documented default is a required, concrete value there. Values that are
semantically inapplicable may remain absent: projectile length and latitude, for example, do not
acquire invented values. In particular, a successful response always materializes
`temperature_k`, `pressure_pa`, and the effective `shot.muzzle_angle_rad`, even when their input
members were omitted.

## Units and conventions

All dimensional input and result fields use SI units, made explicit by their suffix:

- `_kg` is kilograms, `_m` is metres, `_mps` is metres per second, and `_j` is joules.
- `_s` is seconds, `_pa` is pascals, `_k` is kelvin, and `_rad` is radians.
- `relative_humidity` is a dimensionless fraction from `0` through `1`.
- Wind angles are wind-from directions relative to the shot axis: `0` is a headwind and positive
  `pi / 2` is wind from the shooter's right.
- `drop_m` is positive below the line of sight.
- `windage_m` is positive to the shooter's right.
- `twist_direction` is viewed from the breech toward the muzzle.

Numbers must be finite. Physical ranges and service-level cross-field constraints are checked by
the solve service and reported as `invalid_value` or `conflicting_fields`; parsing a DTO does not
itself run a trajectory. The request decoder does enforce structural conflicts explicitly defined
by this contract, including the mutually exclusive effects below.

## Request

A representative request is:

```json
{
  "schema_version": 1,
  "projectile": {
    "mass_kg": 0.01134,
    "diameter_m": 0.00671,
    "length_m": 0.031,
    "drag_model": "G7",
    "ballistic_coefficient": 0.243
  },
  "rifle": {
    "muzzle_velocity_mps": 823.0,
    "sight_height_m": 0.0381,
    "muzzle_height_m": 0.0,
    "twist_rate_m_per_turn": 0.2032,
    "twist_direction": "right"
  },
  "shot": {
    "max_range_m": 1000.0,
    "zero_distance_m": 100.0,
    "aim_azimuth_rad": 0.0,
    "shot_azimuth_rad": 0.0,
    "shooting_angle_rad": 0.0,
    "cant_angle_rad": 0.0,
    "target_height_m": 0.0381,
    "ground_threshold_m": -100.0
  },
  "atmosphere": {
    "altitude_m": 250.0,
    "temperature_k": 288.15,
    "pressure_pa": 101325.0,
    "relative_humidity": 0.5,
    "latitude_rad": 0.7853981633974483
  },
  "wind": {
    "speed_mps": 4.4704,
    "direction_from_rad": 1.5707963267948966,
    "vertical_speed_mps": 0.0
  },
  "solver": {
    "method": "rk45",
    "time_step_s": 0.001
  },
  "effects": {
    "magnus": false,
    "coriolis": true,
    "enhanced_spin_drift": true
  },
  "sampling": {
    "interval_m": 10.0
  }
}
```

### `projectile`

| Field | Required | Meaning |
| --- | --- | --- |
| `mass_kg` | yes | Projectile mass. |
| `diameter_m` | yes | Projectile diameter. |
| `length_m` | no | Projectile length; required by effects that need geometry. |
| `drag_model` | yes | One of `G1`, `G2`, `G5`, `G6`, `G7`, `G8`, `GI`, `GS`, or `RA4`. |
| `ballistic_coefficient` | yes | BC for the selected reference drag model. |

All nine built-in reference drag models are backed by distinct tables and are accepted by v1.
The enum spellings are exact and case-sensitive. Custom drag files and tables remain outside
this wire format; see [Deliberate v1 exclusions](#deliberate-v1-exclusions).

### `rifle`

| Field | Required | Default | Meaning |
| --- | --- | --- | --- |
| `muzzle_velocity_mps` | yes | — | Projectile speed at the muzzle. |
| `sight_height_m` | no | `0.05` | Sight height above the bore. |
| `muzzle_height_m` | no | `0` | Bore height above the ground reference. |
| `twist_rate_m_per_turn` | no | `0.3048` | Rifling travel per full turn. Omitting it does **not** mean "no twist" and does not disable anything that reads it: the service substitutes `0.3048` m — exactly 1:12 inches — raises a `default_applied` assumption notice, and solves against that barrel. See [Omitting the twist rate](#omitting-the-twist-rate). |
| `twist_direction` | no | `right` | `left` or `right`. |
| `sight_offset_lateral_m` | no | `0` | Lateral sight-to-bore mount offset (MBA-1396): positive = sight RIGHT of bore. The trajectory starts that far left of the sight line; with `zero_distance_m` the windage zero converges it onto the sight line at the zero range. Must be finite and smaller than 0.5 m in magnitude. Echoed in `resolved_request.rifle.sight_offset_lateral_m` when supplied; omitting it is byte-identical to requests that predate it, with no assumption notice for its absence. |

#### Omitting the twist rate

`twist_rate_m_per_turn` defaults like any other optional field, which invites reading the absence
as "do not model spin". It is not. The resolved request carries `0.3048` m either way, and
nothing downstream can tell an assumed barrel from a stated one, so every model that reads the
twist runs — against 1:12.

There is no configuration in which omitting it is free. `summary.stability_factor` reads the
resolved twist on **every** solve — no flag opts into it — so an omitted field publishes the Sg of
a 1:12 barrel. On a .308 175 gr with no effects enabled at all, Sg is `1.6650154161603787` with the
twist omitted, bit-identically `1.6650154161603787` with 1:12 stated, and `3.7462846863608514` with
1:8 stated: 2.25x, across the line a shooter reads Sg to decide. That case raises
`stability_factor_assumed_twist_rate` at `$.rifle.twist_rate_m_per_turn`, since there is no flag to
attach it to.

What an omitted twist does *not* move, as long as no spin-driven effect is enabled, is the
trajectory: `drop_m` and `windage_m` are bit-identical for every twist rate. The trajectory becomes
twist-dependent as soon as one of these is in the request:

- `effects.enhanced_spin_drift` — the Litz drift scales with the muzzle Sg, which is a function of
  the twist. On a .308 175 gr at 800 m, `summary.spin_drift_m` is 0.183 m with the twist omitted,
  the identical 0.183 m with 1:12 stated, and 0.316 m with 1:8 stated: 73% more drift from the
  barrel alone, and the omitted request is indistinguishable from the 1:12 one.
- `effects.magnus` — the side force and yaw of repose are driven by the spin rate. Its absolute
  contribution is small on a flat-fire shot (0.19 mm of drop at 800 m for the assumed 1:12) but it
  is entirely twist-bound: a stated 1:6 makes the same contribution 0.75 mm, four times as large.
- `effects.aerodynamic_jump` — covered in [`effects.aerodynamic_jump`](#effectsaerodynamic_jump)
  below, which also reads `projectile.length_m`.

Enabling `magnus` or `enhanced_spin_drift` without stating the twist therefore raises
`spin_effect_assumed_twist_rate` at the enabled flag's path, and `aerodynamic_jump` raises
`aerodynamic_jump_assumed_geometry` at its own — all in addition to the `default_applied`
assumption notice for the default itself, which says a default was applied without saying that
anything now depends on it.

**These codes are not mutually exclusive.** Each names the consumer it belongs to, so one omitted
`twist_rate_m_per_turn` raises one warning per consumer that ran. A request enabling both
`aerodynamic_jump` and `magnus` without a stated twist gets three warnings from the single
omission — `stability_factor_assumed_twist_rate` at `$.rifle.twist_rate_m_per_turn`,
`spin_effect_assumed_twist_rate` at `$.effects.magnus`, and `aerodynamic_jump_assumed_geometry` at
`$.effects.aerodynamic_jump` — each with its own message naming that same missing field. Match on
the code you care about rather than assuming at most one is present.

Twist *direction* defaults separately, to `right`, and flips the sign of the spin drift rather than
its magnitude. It carries its own `default_applied` notice.

### `shot`

| Field | Required | Default | Meaning |
| --- | --- | --- | --- |
| `max_range_m` | yes | — | Requested downrange termination distance. |
| `zero_distance_m` | no | absent | Solve the muzzle elevation for this zero distance. |
| `muzzle_angle_rad` | no | absent | Supply muzzle elevation directly. |
| `aim_azimuth_rad` | no | `0` | Small horizontal aim offset in the sight frame. |
| `shot_azimuth_rad` | no | `0` | Compass bearing used for Earth-rotation effects; `0` is north. |
| `shooting_angle_rad` | no | `0` | Uphill/downhill line-of-sight angle. |
| `cant_angle_rad` | no | `0` | Clockwise rifle cant is positive from the shooter's view. |
| `target_height_m` | no | the line of sight when a zero is solved, else `0` | World-vertical target height above the ground reference, used as the height the elevation search converges on. When the search runs (`zero_distance_m` present, `muzzle_angle_rad` absent) and this field is omitted, it defaults to the LINE OF SIGHT — `rifle.muzzle_height_m + rifle.sight_height_m` — so the solved trajectory crosses the line of sight at the zero distance. Otherwise it defaults to `0`. See the zeroing notes below. |
| `ground_threshold_m` | no | `-100` | Stop after the projectile falls below this height. |
| `zero_poi_up_m` | no | `0` | Deliberate vertical POI offset AT the zero range (MBA-1359, Kestrel "zero height"): positive = deliberately zeroed to impact HIGH by this much at `zero_distance_m`. Must be finite and smaller than 1 m in magnitude. |
| `zero_poi_right_m` | no | `0` | Deliberate horizontal POI offset AT the zero range (MBA-1359, Kestrel "zero offset"): positive = impacts RIGHT. Same bounds as `zero_poi_up_m`. |
| `drops_reference` | no | `"los"` | Which plane sample `drop_m` values are referenced to (MBA-1403). `"los"` = perpendicular to the line of sight (the historical behavior); `"target"` = vertical in the target plane: `drop_m` divided by `cos(shooting_angle_rad)` (JBM's "target plane" reference). |

`zero_distance_m` solves the muzzle elevation; `muzzle_angle_rad` supplies it directly. A request
may supply either, both, or neither. With only `zero_distance_m` present, the service solves for
the muzzle angle that hits it. With `muzzle_angle_rad` present — alone, or together with
`zero_distance_m` (0.33.0 decision-support: this is exactly what rebuilding a request from a
previous `resolved_request` produces, via `From<&ResolvedSolveRequestV1> for SolveRequestV1`) —
the elevation search does not run and the supplied angle is used directly. `zero_distance_m`, when
also present in that case, no longer re-derives the elevation, but is not otherwise inert: it is
still validated, still widens the required wind coverage, still gates
`summary.equivalent_horizontal_range_m` (below), and still drives the windage-convergence bias
described under `sight_offset_lateral_m` and `zero_poi_right_m`. The service emits a
`zero_distance_elevation_not_resolved` warning whenever both fields are supplied together, naming
exactly this. When neither field is present, the service uses a zero muzzle angle and records that
assumption in the response. In `resolved_request`, `muzzle_angle_rad` is always the effective angle
used by the engine. If the caller supplied a zero distance, the resolved shot contains both the
original `zero_distance_m` intent and the muzzle angle used for it, whether that angle was solved
or supplied directly.

`zero_poi_up_m` and `zero_poi_right_m` describe an angular zero-state bias (offset divided by the
zero distance). `zero_poi_up_m` is applied only by the elevation search itself, so it has no effect
whenever `muzzle_angle_rad` is supplied directly, regardless of whether `zero_distance_m` is also
present. `zero_poi_right_m` instead shares the windage-convergence bias `sight_offset_lateral_m`
uses: it has no effect only when `zero_distance_m` is entirely absent (including when
`muzzle_angle_rad` is supplied alone), and still applies whenever `zero_distance_m` is present,
even alongside an explicit `muzzle_angle_rad`. They are echoed in
`resolved_request.shot.zero_poi_up_m` / `zero_poi_right_m` when supplied; the resolved
`muzzle_angle_rad` separately and always reports the biased effective angle regardless of whether
the bias fields themselves were supplied. Omitting both fields is byte-identical to requests that
predate them, and no assumption notice is emitted for their absence.

The zero search uses the request's resolved projectile, atmosphere, wind (including downrange
segments), and integration method, and the request's effects apart from the one carve-out named
below. It follows the engine's level-rifle convention by solving with zero cant; the requested
`cant_angle_rad` is applied only to the subsequent trajectory. `target_height_m` remains an
absolute world-vertical height above the local ground datum, as named above; inclined zeroing
projects the shot-frame trajectory back into that world frame.

**What the search converges on.** `target_height_m` is the height the search drives the bullet to
at `zero_distance_m`. Omitting it while the search runs defaults it to the line of sight of a
level rifle — `rifle.muzzle_height_m + rifle.sight_height_m` — so for a level shot "zero at 100 m"
means what a shooter means by it: the trajectory crosses the line of sight there. A supplied value
always wins, `0.0` included, which is how a caller asks for a zero against the ground datum
instead. The service emits a `default_applied` assumption at `$.shot.target_height_m` naming the
height it used whenever it applies this default, because the defaulted value normally moves every
elevation number in the response. It need not: a request whose `sight_height_m` and
`muzzle_height_m` are both `0` resolves the default to `0`, and then only the notice differs. The default is gated on the search actually running: with
`muzzle_angle_rad` supplied there is no zero to frame, and `target_height_m` then still defaults
to `0` (where it feeds only the `drops_reference: "target"` sampler datum, which references
`max_range_m` rather than `zero_distance_m`).

That default is a LEVEL line of sight, and this field is a world-vertical height, so it does not
describe the sight line of an inclined shot. A request that pairs `zero_distance_m` with a nonzero
`shooting_angle_rad` and omits `target_height_m` is not zeroing to its own line of sight: uphill
it fails to converge, and downhill it converges on the world height the default names rather than
on the sight line. Supply the height explicitly for an inclined zero — the sight line at the zero
distance, projected into the world frame, is
`zero_distance_m * sin(shooting_angle_rad) + (muzzle_height_m + sight_height_m) * cos(shooting_angle_rad)`.

**`effects.aerodynamic_jump` is excluded from the search**, deliberately. A rifle is zeroed in
calm air; letting a crosswind-driven jump term into the zero trials would bake the wind of the
zeroing session into the stored elevation and into every solve made from it. The search therefore
runs its trials with the jump off (`zero_trial_height_at`, MBA-959) and the jump stays what it is
meant to be — an additive fire-time launch-angle perturbation. The visible consequence, which is
physics and not a defect: with `aerodynamic_jump` enabled and a crosswind, the solved elevation is
unchanged, so the trajectory sits off its own stated zero AT the zero distance by the jump. That
offset is reported as `summary.aerodynamic_jump_moa`. Wind and Coriolis are not carved out; both
reach the search and move the solved elevation.

`drops_reference` is an output-mode toggle only: it rescales each sample's `drop_m` and changes
nothing else — not the solved trajectory, not `windage_m`, not the `summary` block, and not
zeroing (which keeps its own `target_height_m` semantics). With `"target"` and
`|shooting_angle_rad| >= 90 degrees` the transform is undefined and the request fails with a
solve error. It is echoed in `resolved_request.shot.drops_reference` when supplied. Omitting the
field is byte-identical to requests that predate it, and no assumption notice is emitted for its
absence; explicitly supplying `"los"` solves identically to omitting the field (`"los"` is the
behavioral default), but is no longer byte-identical at the envelope level, since the echo itself
then appears where an omitted-field response leaves it absent.

### `atmosphere`

| Field | Required | Default | Meaning |
| --- | --- | --- | --- |
| `altitude_m` | no | `0` | Station altitude. |
| `temperature_k` | no | ICAO at `altitude_m` | Authoritative station temperature when present. |
| `pressure_pa` | no | ICAO at `altitude_m` | Station or sea-level pressure when present; see `pressure_reference`. |
| `pressure_reference` | no | `"absolute"` | Whether `pressure_pa` is absolute station pressure or a QNH altimeter setting. |
| `relative_humidity` | no | `0.5` | Relative-humidity fraction. |
| `latitude_rad` | no | absent | Geodetic latitude; needed when Coriolis is enabled. |

An empty atmosphere object selects ICAO standard conditions at the resolved altitude (sea level
when `altitude_m` is also omitted). Enabling Coriolis without a latitude is an `invalid_value`
error rather than a silently chosen latitude.

Explicit temperature and pressure values are authoritative station conditions, including values
equal to `288.15 K` and `101325 Pa` at nonzero altitude. The solve service must preserve that
explicit intent and must not apply the legacy CLI rule that treats those exact values as omitted
standard-atmosphere sentinels. MBA-1302 must bypass that CLI sentinel inference: omitted fields
select ICAO-at-altitude resolution, while present fields are passed as authoritative values. The
resolved values are recorded as explicit numbers in a successful response's
`resolved_request`; the resolved altitude and relative humidity are concrete there as well.

#### `pressure_reference` (MBA-1397)

`pressure_pa` can mean two different physical quantities, and the caller must say which:

- `"absolute"` (the default, and the only meaning before this field existed): `pressure_pa` is
  already the absolute station pressure at `altitude_m`. Used as-is.
- `"qnh"`: `pressure_pa` is a sea-level-corrected altimeter setting (a weather-report barometer
  / METAR-style QNH reading). It is reduced to station pressure at `altitude_m` via the ICAO
  inverse-barometric formula (`station = QNH * (1 - 0.0065*h/288.15)^5.25588`) before use, and
  the resolved `pressure_pa` in a successful response's `resolved_request` is the REDUCED
  station pressure, not the raw QNH the caller sent. The reduction is recorded as an
  `assumptions` entry with code `qnh_reduced_to_station_pressure`.

The mode itself is separately echoed in `resolved_request.atmosphere.pressure_reference`
whenever the request supplies one, including an explicit `"absolute"` — which therefore solves
identically to omission (`"absolute"` is the behavioral default) but is no longer byte-identical
at the envelope level, since the echo itself then appears where an omitted-field response leaves
it absent (the same separate-mode-echo convention `wind_reference` below uses).

`pressure_reference` has no effect when `pressure_pa` is omitted: an omitted pressure always
resolves to the ICAO standard station pressure at `altitude_m`, which is mathematically the same
result as reducing a QNH of exactly `101325 Pa` (the ICAO sea-level standard).

### `wind`

An empty object means still air. Constant wind uses:

- `speed_mps` and `direction_from_rad`, which must be supplied together;
- optional `vertical_speed_mps`, defaulting to zero when constant wind is selected.

Downrange wind uses `segments` instead:

```json
{
  "segments": [
    {
      "until_distance_m": 300.0,
      "speed_mps": 2.0,
      "direction_from_rad": 0.0,
      "vertical_speed_mps": 0.0
    },
    {
      "until_distance_m": 1000.0,
      "speed_mps": 5.0,
      "direction_from_rad": 1.5707963267948966
    }
  ]
}
```

Segment boundaries must increase strictly. `segments` conflicts with all three constant-wind
fields. A partial constant wind or overlapping segment boundaries are reported as
`conflicting_fields` or `invalid_value` by the service. Segments may end before the requested
range; the engine uses still air beyond the final boundary and the response includes a
`partial_wind_coverage` warning. Coverage is checked through the farther of `max_range_m` and an
applicable `zero_distance_m`, because the zero trial uses the same segmented wind.

Input presence is preserved for wind too: an omitted `segments` member is distinct from an
explicit array, and an omitted segment `vertical_speed_mps` remains absent until resolution.
Resolved wind is exactly one of two object shapes: a constant object with concrete `speed_mps`,
`direction_from_rad`, and `vertical_speed_mps`, or a segmented object whose segments each have a
concrete vertical speed. Still air resolves to the constant shape with all three values set to
zero.

Optional `wind_reference` (MBA-1368) selects the frame every wind direction in the request is
entered in: omitted means shooter-relative wind-FROM radians, byte-identical to requests that
predate the field, with no assumption notice for its absence; `"compass"` means earth-fixed
bearings (0 = north) — the constant `direction_from_rad` AND every segment's — which the service
re-references against the shot azimuth at resolve time as `bearing - shot.shot_azimuth_rad`,
normalized to `[0, 2π)`. The RESOLVED wind echo always reports the converted shooter-relative
direction (the same fold-into-the-resolved-value convention QNH pressure uses); the mode itself is
separately echoed in `resolved_request.wind.wind_reference` whenever the request supplies one,
including an explicit `"shooter"` — which therefore solves identically to omission (`"shooter"` is
the behavioral default) but is no longer byte-identical at the envelope level, since the echo
itself then appears where an omitted-field response leaves it absent. `"compass"` requires an
explicit `shot.shot_azimuth_rad` — omitting it is a `conflicting_fields` error at
`$.wind.wind_reference`, never a silent treat-as-shooter-relative. A wind FROM north
(`direction_from_rad: 0`) on a shot fired due north (`shot_azimuth_rad: 0`) is a pure headwind.

### `solver`, `effects`, and `sampling`

| Section and field | Default | Meaning |
| --- | --- | --- |
| `solver.method` | `rk45` | `rk45`, `rk4`, or `euler`. |
| `solver.time_step_s` | `0.001` | Fixed step for RK4 and Euler; accepted but ignored by adaptive RK45. |
| `effects.magnus` | `false` | Enable the engine's Magnus-force model. |
| `effects.coriolis` | `false` | Enable Earth-rotation deflection. |
| `effects.enhanced_spin_drift` | `false` | Enable enhanced spin-drift modeling. |
| `effects.wind_shear_model` | omitted (`"none"`) | Altitude-dependent wind shear: `none`, `logarithmic`, `power_law`, or `ekman_spiral` (alias `ekman`). |
| `effects.aerodynamic_jump` | omitted (`false`) | Enable crosswind aerodynamic (gyroscopic) jump as a muzzle launch-angle perturbation. |
| `sampling.interval_m` | `10` | Regular downrange result interval. |

Supplying `solver.time_step_s` with `rk45` is valid, but RK45 owns its adaptive step size. The
resolved request retains the supplied value and the response includes an
`rk45_time_step_ignored` warning.

A v1 success response contains at most 10,000 trajectory samples. Exactly 10,000 is valid;
10,001 is not. Sampling is evaluated against the trajectory's actual reached range, so an early
ground or time termination can keep a fine-grid response within the limit. If the completed
trajectory would produce more than 10,000 observations, the service fails with `resource_limit` at
`$.sampling.interval_m` before allocating or serializing the response. The service must not
truncate or thin the requested sample sequence to fit the limit.

Effects remain opt-in. The service may require projectile length, twist data, latitude, or other
documented prerequisites when a corresponding effect is enabled.

Magnus, enhanced spin drift and aerodynamic jump are experimental engine models. Enabling any of
them produces an `experimental_effect` warning at the corresponding request path.

`effects.magnus` and `effects.enhanced_spin_drift` cannot both be true in v1. The engine's legacy
solver silently suppresses Magnus in that combination; the request decoder instead reports
`conflicting_fields` so the resolved request never misstates which physics ran.

Both of them read `rifle.twist_rate_m_per_turn`, and omitting it does not disable them — it solves
them against the assumed 1:12 barrel and raises `spin_effect_assumed_twist_rate`. See
[Omitting the twist rate](#omitting-the-twist-rate).

### `effects.aerodynamic_jump`

Crosswind aerodynamic jump: the fixed angular departure a spinning projectile takes as it leaves
the constrained bore, applied as an initial launch-angle offset rather than as a downrange force.
Omitted or `false` is byte-identical to every response from before the field existed. An explicit
`false` is echoed in `resolved_request.effects` while an omitted field is not, so a round-tripped
request says exactly what the original said.

The model is Bryan Litz's regression, `Y = 0.01*Sg - 0.0024*L + 0.032` MOA per mph of crosswind,
fed by the engine's own Miller stability factor. It is a fit that is best near `Sg` ~ 1.75, not a
first-principles derivation, which is why enabling it raises `experimental_effect`.

Three properties a caller has to know, because none of them is visible in the trajectory:

- **It is vertical.** The jump perturbs elevation only. Windage is unchanged to within rounding,
  and the effect is independent of `magnus`, `coriolis` and `enhanced_spin_drift` — there is no
  suppression rule between them.
- **It is computed from the crosswind AT THE MUZZLE.** A shot with no crosswind at the muzzle is
  unaffected however much wind it meets downrange, including a segmented-wind request whose first
  segment is calm. That case is a present, exactly-zero `summary.aerodynamic_jump_moa`, which is a
  different response from the effect never running (the field is absent then).
- **It reads the barrel and the bullet.** The jump scales with stability and bullet length, so it
  consumes `rifle.twist_rate_m_per_turn`, `rifle.twist_direction` and `projectile.length_m`.
  Omitting them does **not** disable the correction and does not zero it: the solve proceeds
  against the assumed 1:12 twist and the mass/diameter length estimate and returns a confident
  number computed for a rifle the request never described. A .308 175 gr at 800 m in a 10 mph
  full-value crosswind moves 10.79 cm with a stated 1:10 twist and 9.08 cm with the twist omitted.
  Enabling the flag without either field therefore raises `aerodynamic_jump_assumed_geometry` at
  `$.effects.aerodynamic_jump`, in addition to the ordinary assumption notices for the defaults
  themselves. Jump is not the only consumer of the barrel — `effects.magnus` and
  `effects.enhanced_spin_drift` read it under the separate `spin_effect_assumed_twist_rate` code,
  and `summary.stability_factor` under `stability_factor_assumed_twist_rate` on every solve. The
  codes are distinct so each names its own consumer; they are not mutually exclusive, so enabling
  jump alongside `magnus` with the twist omitted returns all three at once. See
  [Omitting the twist rate](#omitting-the-twist-rate).

The applied jump is reported as `summary.aerodynamic_jump_moa`, in MOA, positive up. It is present
only when the effect ran, because a launch-angle offset is folded into every drop in the table and
is otherwise unobservable: a caller comparing two solves cannot tell a jump that applied from one
that silently came out at zero.

### `effects.wind_shear_model`

Additive and optional: omitting it is byte-identical to requests from before it existed, echo
included, and solves exactly as an explicit `"none"` does. Wind shear scales the request's wind by
a boundary-layer profile keyed off the projectile's height above the ground — the operative wind
is a floor, so the profile only ever increases it. Below the 10 m meteorological reference height
the ratio is exactly 1.0, which means an ordinary flat-fire solve is unaffected by any model; the
shear becomes visible on lofted, high-angle, and ELR trajectories that climb well clear of it.

The height driving the profile is height above the *muzzle* (plus an assumed 1.5 m muzzle height),
not above sea level. `atmosphere.altitude_m` is deliberately not an input to it: boundary-layer
shear is relative to the local ground, and altitude's effect on the shot is air density, which
`atmosphere.altitude_m` already carries.

| Value | Meaning |
| --- | --- |
| omitted / `none` | No altitude dependence. The request's wind applies at every height. |
| `logarithmic` | `ln(z / z0) / ln(z_ref / z0)`, with `z0 = 0.03 m` (short grass) and `z_ref = 10 m`. |
| `power_law` | `(z / z_ref)^(1/7)`, the neutral-stability power law. |
| `ekman_spiral` (alias `ekman`) | Accepted for parity with the CLI's `--wind-shear-model`, but this solve path's evaluator has no near-ground profile for it and leaves the wind at the operative value. Emits a `wind_shear_model_not_modeled` warning; it is never silently inert. |

An unrecognized model is an `invalid_value` error at `$.effects.wind_shear_model` listing the
accepted spellings — never a silent fall back to `"none"`, which would return unsheared numbers
indistinguishable from sheared ones. The engine's `custom_layers` model is deliberately not
exposed: it needs a caller-supplied layer table this contract has no field for.

The model is echoed at `resolved_request.effects.wind_shear_model` whenever the request supplies
one, including an explicit `"none"` — which therefore solves identically to omission but is no
longer byte-identical at the envelope level, since the echo then appears where an omitted-field
response leaves it absent. Aliases are canonicalized in the echo: `"ekman"` in, `"ekman_spiral"`
out.

Wind shear cannot be combined with `wind.segments`: downrange segments plus altitude shear is not
a defined model, and the solver's segment lookup takes precedence over its shear branch, so the
pair would drop the shear while the resolved request still named a model. The combination is a
`conflicting_fields` error at `$.effects.wind_shear_model`, matching the CLI's and WASM front
end's refusal of `--wind-segment` with `--enable-wind-shear`.

### Optional `corrections` block

| Field | Required | Meaning |
| --- | --- | --- |
| `bc5d_table_path` | no | Local filesystem path to a caliber-specific BC5D correction table (`bc5d_<caliber>.bin`). |

Additive and optional: omitting the block is byte-identical to requests from before it existed.
When `bc5d_table_path` is supplied, the service loads the named table (the exact dual-CRC binary
format the CLI's `--bc-table-dir` consumes — the caller, typically a mobile app, downloads the
table itself and hands the engine a path), verifies the stored CRC32, generates the same
velocity-keyed BC segment ladder the CLI generates for this load, and folds the table's muzzle
correction into the scalar BC as the schedule's interior-gap fallback. Both the zero-distance
elevation search and the trajectory integrate with that schedule, so a zeroed solve cannot be
mis-zeroed against a different BC than the flight uses. A table that carries no meaningful
correction for the load (every sampled cell ≈ 1.0) leaves the request's constant BC in place.

The BC5D tables carry G1 and G7 planes only; a request naming any other reference drag model is
solved with the lookup typed as G1 and receives a `bc5d_drag_model_coerced` warning.

The table must be for the shot's caliber. Its header declares one, and the service compares that
against `projectile.diameter_m` by the same 3-digit key that names the file — `round(caliber ×
1000)`, so a `308` table covers `[0.3075, 0.3085)` in — and refuses a mismatch as an
`invalid_value` naming both values (`table is for 0.224, shot is 0.308`). This is a refusal, not a
downgrade: a wrong-caliber table never fails a lookup (values clamp to its edge bins), so applying
one would silently bias every row, and answering with an unannotated uncorrected solve would be
indistinguishable from a corrected one. Use the bridge's `bc5d.info` (`caliber_key`) to pre-check a
downloaded table.

A missing or unreadable file is an `io_error`; a file that is not a valid BC5D table (bad magic,
unsupported version, impossible dimensions, or a CRC mismatch) is an `invalid_value`, both at
`$.corrections.bc5d_table_path`. On engine builds without filesystem access (wasm32), supplying
the field is an `invalid_value` at the same path, never a silent no-op.

The block is echoed at `resolved_request.corrections` when supplied. Segments are always
regenerated from the request's published `ballistic_coefficient`, so re-solving a resolved
request re-applies the same correction rather than compounding it.

## Success envelope

`status` is the literal `"ok"`. `engine_version` identifies the engine implementation that
produced the result. `resolved_request` is a `ResolvedSolveRequestV1`, not a replay of the
presence-aware input DTO: all rifle, shot, atmosphere, wind, solver, effects, and sampling
defaults are materialized so the calculation is reproducible. `assumptions` and `warnings` are
arrays of objects with a stable `code`, a human-readable `message`, and an optional request
`path`.

The service emits notices in deterministic request-field order. Stable v1 assumption codes are
`default_applied` for literal defaults, `icao_standard_temperature` and
`icao_standard_pressure` for omitted station values resolved from the requested altitude,
`qnh_reduced_to_station_pressure` for an explicit `pressure_reference: "qnh"` pressure reduced
to station pressure (MBA-1397), and `estimated_projectile_length` when the engine needs
inferred projectile geometry. Stable v1 warning codes are `partial_wind_coverage`,
`experimental_effect`, `rk45_time_step_ignored`, `zero_distance_elevation_not_resolved` (an
explicit `muzzle_angle_rad` supplied together with `zero_distance_m`, so the elevation search did
not run — see `shot.muzzle_angle_rad` above), `bc5d_drag_model_coerced` (a
`corrections.bc5d_table_path` request whose drag model is outside the table's G1/G7 planes — see
the optional `corrections` block above), `wind_shear_model_not_modeled` (an accepted
`effects.wind_shear_model` this solve path has no profile for, so the wind is left unchanged — see
`effects.wind_shear_model` above), `aerodynamic_jump_assumed_geometry` (`effects.aerodynamic_jump`
enabled without `rifle.twist_rate_m_per_turn` or `projectile.length_m`, so the jump was computed
from an assumed barrel rather than disabled — see `effects.aerodynamic_jump` above), and
`spin_effect_assumed_twist_rate` (`effects.magnus` or `effects.enhanced_spin_drift` enabled without
`rifle.twist_rate_m_per_turn`, so the spin-driven model was computed from the assumed 1:12 twist
rather than disabled), and `stability_factor_assumed_twist_rate` (`rifle.twist_rate_m_per_turn`
omitted on a solve that reported `summary.stability_factor`, so the Sg published is the Sg of the
assumed 1:12 barrel — this one is not gated on any flag, because Sg is computed on every solve; see
[Omitting the twist rate](#omitting-the-twist-rate) above for both). More than one of these can be
present at once: they are keyed to the consumer, not to the missing field. Messages are descriptive
text rather than a compatibility surface.

```json
{
  "schema_version": 1,
  "engine_version": "0.24.1",
  "status": "ok",
  "resolved_request": { "...": "the complete v1 request" },
  "assumptions": [],
  "warnings": [],
  "summary": {
    "actual_range_m": 1000.0,
    "maximum_height_m": 3.2,
    "time_of_flight_s": 1.6,
    "terminal_speed_mps": 360.0,
    "terminal_energy_j": 734.8,
    "stability_factor": 1.5,
    "spin_drift_m": 0.12,
    "termination": "max_range"
  },
  "samples": [
    {
      "distance_m": 1000.0,
      "time_s": 1.6,
      "speed_mps": 360.0,
      "energy_j": 734.8,
      "drop_m": 8.1,
      "windage_m": 0.43,
      "mach": 1.06,
      "flags": ["transonic", "terminal"]
    }
  ]
}
```

Summary fields have fixed evaluation frames:

- `maximum_height_m` is the greatest world-vertical projectile height above the same local
  ground/reference datum used by `rifle.muzzle_height_m`. It is not the inclined shot frame's Y
  coordinate and is not height above the line of sight.
- `stability_factor` is the dimensionless muzzle gyroscopic stability factor Sg, evaluated after
  resolving projectile geometry, muzzle velocity, twist, and the station atmosphere. It is absent
  only when the service cannot calculate Sg from the resolved inputs. The twist it uses is the
  *resolved* one, so a request that omitted `rifle.twist_rate_m_per_turn` reports the Sg of a 1:12
  barrel — and, because that is reported rather than merely assumed, raises
  `stability_factor_assumed_twist_rate`. It is the only twist-dependent output warned on every
  solve rather than behind an `effects` flag. See
  [Omitting the twist rate](#omitting-the-twist-rate).
- `spin_drift_m` is the signed gyroscopic spin-drift contribution at the terminal sample, positive
  to the shooter's right. It excludes wind drift and is absent when enhanced spin drift is disabled
  or cannot be calculated.
- `aerodynamic_jump_moa` (MBA-959) is the vertical crosswind aerodynamic jump actually applied at
  the muzzle, in MOA, positive up. Present only when `effects.aerodynamic_jump` was enabled and the
  solver produced components; absent otherwise, so every response that predates the field is
  byte-identical. A present `0.0` means the effect ran against no muzzle crosswind — a different
  fact from an absent field, which means it never ran. It is reported because a launch-angle offset
  is folded into every drop rather than appearing as its own column, and is otherwise invisible.
- `equivalent_horizontal_range_m` (MBA-1395) is the flat-fire range whose angular elevation
  correction — against the same solved zero — matches the inclined solution's at the terminal
  range: the BDC "shoot-to" range (SIG AMR / Leica EHR / Gunwerks style). It is defined by
  angular match over a flat re-solve, not the rifleman's-rule cosine approximation. Present only
  when `shot.shooting_angle_rad` is nonzero, `zero_distance_m` is present (whether or not it
  re-derived the elevation — see above), and the inverse is well-defined (terminal range past the
  zero range with a positive correction); absent otherwise, so flat solves and pre-existing
  requests serialize byte-identically.

### Optional `reticle` block (MBA-1361)

A request may carry an optional top-level `reticle` object. When it does — and only then —
the success envelope gains a top-level `reticle_hold` object. Requests without it, and every
response that predates the field, serialize byte-identically.

The raw block is also echoed verbatim onto `resolved_request.reticle` when supplied (0.33.0
decision-support), completing the resolved request as a full description of the solve. Omitting
`reticle` keeps `resolved_request` byte-identical to requests that predate the field.

```json
"reticle": {
  "range_m": 600.0,
  "magnification": 5.0,
  "description": {
    "name": "MyScope MIL",
    "focal_plane": "sfp",
    "reference_magnification": 10.0,
    "marks": [
      {"down_mil": 0.0, "right_mil": 0.0, "kind": "center"},
      {"down_mil": 2.0, "right_mil": 0.0, "kind": "hash", "label": "600"}
    ]
  }
}
```

The envelope is strict (`range_m`, `magnification` and `description` are all required, and no
other key is accepted). `description` is the shared reticle schema — the same JSON
`ballistics reticle generate -o json` emits — and is deliberately permissive about extra keys
inside it, so a front end's render metadata round-trips.

`reticle_hold` reports, in milliradians from the optical center with `down_mil` positive BELOW
center and `right_mil` positive to the shooter's RIGHT:

```json
"reticle_hold": {
  "range_m": 600.0, "magnification": 5.0,
  "down_mil": 3.18, "right_mil": 0.74, "mark_scale": 2.0,
  "nearest_mark_index": 1, "nearest_mark_label": "600",
  "nearest_mark_distance_mil": 0.22, "off_reticle": false
}
```

The angular values are read from the response's OWN `samples` (linearly interpolated at
`range_m`, milliradian small-angle definition), so the hold and the sample rows can never
describe different trajectories. `mark_scale` is `reference_magnification / magnification` for a
second-focal-plane reticle and exactly `1.0` for first focal plane; the hold coordinates are
always true angular, and only the marks are rescaled. A `range_m` outside the sampled
trajectory returns a structured `invalid_value` error at `$.reticle.range_m` — the service never
extrapolates a hold off a trajectory that did not get there.

`summary.termination` is one of `max_range`, `ground_threshold`, `time_limit`, or
`velocity_floor`. The solve service must populate it from explicit termination metadata returned
by the engine. It must not infer a reason heuristically from the last distance, height, or speed.

Regular interval sampling never drops the terminal observation. `samples` includes the terminal
sample exactly once even when its distance is not an interval boundary, and that sample carries
the `terminal` flag. Sample flags are `transonic`, `subsonic`, `terminal`, and
`ground_threshold`. `transonic` denotes the inclusive Mach 0.8–1.2 band and `subsonic` denotes
Mach below 1.0, so both flags intentionally appear from Mach 0.8 through values just below 1.0.

Rust service implementations must call `SolveSuccessV1::validate_for_serialization` immediately
before encoding a success envelope. The corresponding public ceiling is
`MAX_SOLVE_JSON_SAMPLES_V1` (`10_000`). This explicit check lets the service return the structured
`resource_limit` error envelope. The DTO serializer also rejects an oversized `samples` array as a
fail-closed backstop, but that serializer error is intentionally not a substitute for the protocol
error envelope.

## Error envelope

`status` is the literal `"error"`. Structurally valid JSON uses `path` and null line/column fields;
malformed JSON uses one-based `line` and `column` and a null path. These forms are mutually
exclusive: a path cannot appear with a source location, line and column must appear together, and
neither source coordinate may be zero. Errors unrelated to a particular input location use null
for all three fields. Because `serde_json` reports column zero for some end-of-file errors, the
decoder normalizes that parser-only value to column one before creating the envelope.

```json
{
  "schema_version": 1,
  "status": "error",
  "error": {
    "code": "unknown_field",
    "message": "unknown field `balistic_coefficient`",
    "path": "$.projectile.balistic_coefficient",
    "line": null,
    "column": null
  }
}
```

The v1 error codes are:

| Code | Meaning |
| --- | --- |
| `invalid_json` | The input is not a JSON document; line and column identify the parser error. |
| `unsupported_schema_version` | `schema_version` is not `1`. |
| `unknown_field` | An object contains a field not defined by v1. |
| `missing_field` | A required field or top-level section is absent. |
| `invalid_value` | A value has the wrong type, range, enum, or physical validity. |
| `conflicting_fields` | Individually valid fields cannot be used together. |
| `resource_limit` | The request or requested result exceeds a documented service limit. |
| `solve_failed` | The engine could not complete a valid trajectory solve. |
| `io_error` | The transport could not read input or write output. |
| `internal_error` | An unexpected implementation failure was contained. |

Human-readable messages are diagnostic and are not a compatibility surface. Consumers branch on
`code` and may use `path`, `line`, and `column` to highlight the input.

`solve_v1` is deterministic and transport-free: it performs no filesystem or network access,
does not load profiles, and does not write to stdout or stderr.

## Deliberate v1 exclusions

V1 does not expose custom drag files or tables, caller-authored velocity/Mach-dependent BC
schedules, powder temperature curves, atmosphere zones, cluster-BC degradation, wind shear,
pitch damping, precession/nutation, or angular diagnostics. (The optional
[`corrections`](#optional-corrections-block) block is the one deliberate carve-out: a
velocity-keyed BC schedule generated by the engine itself from a verified BC5D table, not an
arbitrary caller-supplied schedule.) It also does not expose network access or filesystem access
beyond reading the one named correction table. Those features require a later schema or a
separately versioned extension after their input semantics and result provenance are stable.
