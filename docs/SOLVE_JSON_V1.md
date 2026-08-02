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

v1 does, however, grow by ADDITION under two strict rules, so that a feature need not force a
whole new schema:

1. A new REQUEST field must be optional and default to the exact pre-existing behavior when
   omitted. Every request valid before the field remains valid and produces identical results
   (`zero_poi_up_m`, `zero_poi_right_m`, `sight_offset_lateral_m`, `drops_reference`,
   `wind_reference` were added this way).
2. A new RESPONSE field must be omitted whenever its feature is inactive, so a response that does
   not exercise the feature is byte-identical to one produced before the field existed
   (`reticle_hold` appears only when the request carries a `reticle` block).

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
    "target_height_m": 0.0,
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
| `drag_model` | yes | One of `G1`, `G6`, `G7`, or `G8`. |
| `ballistic_coefficient` | yes | BC for the selected reference drag model. |

Only reference drag models backed by distinct tables in the engine are part of v1. The engine's
legacy `G2`, `G5`, `GI`, and `GS` enum values currently fall back to the G1 table, so v1 rejects
them rather than silently claiming a distinct model.

### `rifle`

| Field | Required | Default | Meaning |
| --- | --- | --- | --- |
| `muzzle_velocity_mps` | yes | — | Projectile speed at the muzzle. |
| `sight_height_m` | no | `0.05` | Sight height above the bore. |
| `muzzle_height_m` | no | `0` | Bore height above the ground reference. |
| `twist_rate_m_per_turn` | no | `0.3048` | Rifling travel per full turn. |
| `twist_direction` | no | `right` | `left` or `right`. |
| `sight_offset_lateral_m` | no | `0` | Lateral sight-to-bore mount offset (MBA-1396): positive = sight RIGHT of bore. The trajectory starts that far left of the sight line; with `zero_distance_m` the windage zero converges it onto the sight line at the zero range. Must be finite and smaller than 0.5 m in magnitude. Echoed in `resolved_request.rifle.sight_offset_lateral_m` when supplied; omitting it is byte-identical to requests that predate it, with no assumption notice for its absence. |

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
| `target_height_m` | no | `0` | Target height above the ground reference for zeroing. |
| `ground_threshold_m` | no | `-100` | Stop after the projectile falls below this height. |
| `zero_poi_up_m` | no | `0` | Deliberate vertical POI offset AT the zero range (MBA-1359, Kestrel "zero height"): positive = deliberately zeroed to impact HIGH by this much at `zero_distance_m`. Must be finite and smaller than 1 m in magnitude. |
| `zero_poi_right_m` | no | `0` | Deliberate horizontal POI offset AT the zero range (MBA-1359, Kestrel "zero offset"): positive = impacts RIGHT. Same bounds as `zero_poi_up_m`. |
| `drops_reference` | no | `"los"` | Which plane sample `drop_m` values are referenced to (MBA-1403). `"los"` = perpendicular to the line of sight (the historical behavior); `"target"` = vertical in the target plane: `drop_m` divided by `cos(shooting_angle_rad)` (JBM's "target plane" reference). |

`zero_distance_m` solves the muzzle elevation; `muzzle_angle_rad` supplies it directly. A request
may supply either, both, or neither. With only `zero_distance_m` present, the service solves for
the muzzle angle that hits it. With `muzzle_angle_rad` present — alone, or together with
`zero_distance_m` (0.33.0 decision-support: this is exactly what rebuilding a request from a
previous `resolved_request` produces, via `From<&ResolvedSolveRequestV1> for SolveRequestV1`) —
it is used directly and no zero search runs; `zero_distance_m`, if also present, is then retained
only as the caller's original zeroing intent and has no effect on the solve. When neither is
present, the service uses a zero muzzle angle and records that assumption in the response. In
`resolved_request`, `muzzle_angle_rad` is always the effective angle used by the engine. If the
caller supplied a zero distance, the resolved shot contains both the original `zero_distance_m`
intent and the muzzle angle used for it, whether that angle was solved or supplied directly.

`zero_poi_up_m` and `zero_poi_right_m` describe an angular zero-state bias (offset divided by the
zero distance, applied to the solved elevation and azimuth after the zero search converges). They
are meaningful only together with `zero_distance_m`; with `muzzle_angle_rad` supplied directly (or
neither field present) no zero is solved and they have no effect. They are echoed in
`resolved_request.shot.zero_poi_up_m` / `zero_poi_right_m` when supplied; the resolved
`muzzle_angle_rad` separately and always reports the biased effective angle regardless of whether
the bias fields themselves were supplied. Omitting both fields is byte-identical to requests that
predate them, and no assumption notice is emitted for their absence.

The zero search uses the request's resolved projectile, atmosphere, wind (including downrange
segments), effects, and integration method. It follows the engine's level-rifle convention by
solving with zero cant; the requested `cant_angle_rad` is applied only to the subsequent trajectory.
`target_height_m` remains an absolute world-vertical height above the local ground datum, as named
above; inclined zeroing projects the shot-frame trajectory back into that world frame.

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

Magnus and enhanced spin drift are experimental engine models. Enabling either produces an
`experimental_effect` warning at the corresponding request path.

`effects.magnus` and `effects.enhanced_spin_drift` cannot both be true in v1. The engine's legacy
solver silently suppresses Magnus in that combination; the request decoder instead reports
`conflicting_fields` so the resolved request never misstates which physics ran.

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
`experimental_effect`, and `rk45_time_step_ignored`. Messages are descriptive text rather than
a compatibility surface.

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
  only when the service cannot calculate Sg from the resolved inputs.
- `spin_drift_m` is the signed gyroscopic spin-drift contribution at the terminal sample, positive
  to the shooter's right. It excludes wind drift and is absent when enhanced spin drift is disabled
  or cannot be calculated.
- `equivalent_horizontal_range_m` (MBA-1395) is the flat-fire range whose angular elevation
  correction — against the same solved zero — matches the inclined solution's at the terminal
  range: the BDC "shoot-to" range (SIG AMR / Leica EHR / Gunwerks style). It is defined by
  angular match over a flat re-solve, not the rifleman's-rule cosine approximation. Present only
  when `shot.shooting_angle_rad` is nonzero, a `zero_distance_m` was solved, and the inverse is
  well-defined (terminal range past the zero range with a positive correction); absent otherwise,
  so flat solves and pre-existing requests serialize byte-identically.

### Optional `reticle` block (MBA-1361)

A request may carry an optional top-level `reticle` object. When it does — and only then —
the success envelope gains a top-level `reticle_hold` object. Requests without it, and every
response that predates the field, serialize byte-identically.

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

V1 does not expose custom drag files or tables, velocity/Mach-dependent BC schedules, powder
temperature curves, atmosphere zones, cluster-BC degradation, wind shear, pitch damping,
precession/nutation, or angular diagnostics. It also does not expose arbitrary filesystem or
network access. Those features require a later schema or a separately versioned extension after
their input semantics and result provenance are stable.
