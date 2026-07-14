# Solve JSON v1 contract

Solve JSON is the binding-neutral request and result contract for one deterministic trajectory
solve. It is separate from `BallisticInputs`, command-line profile JSON, FFI structs, and any
language binding. Those implementation APIs may change without changing this contract.

The Rust DTOs and the envelope-producing decoder are public in
`ballistics_engine::solve_json`. Transport and conversion into engine inputs are outside this
contract.

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
sign, renaming an enum value, or changing a field's meaning requires a new schema version.

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

`zero_distance_m` and `muzzle_angle_rad` conflict in an input request. A request supplies at most
one. When neither is present, the service uses a zero muzzle angle and records that assumption in
the response. In `resolved_request`, `muzzle_angle_rad` is always the effective angle used by the
engine. If the caller requested a zero distance, the resolved shot contains both the original
`zero_distance_m` intent and the muzzle angle calculated for it.

### `atmosphere`

| Field | Required | Default | Meaning |
| --- | --- | --- | --- |
| `altitude_m` | no | `0` | Station altitude. |
| `temperature_k` | no | ICAO at `altitude_m` | Authoritative station temperature when present. |
| `pressure_pa` | no | ICAO at `altitude_m` | Authoritative station pressure when present. |
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

Segment boundaries must increase strictly and cover the requested range. `segments` conflicts
with all three constant-wind fields. A partial constant wind, overlapping segment boundaries, or
an uncovered range is reported as `conflicting_fields` or `invalid_value` by the service.

Input presence is preserved for wind too: an omitted `segments` member is distinct from an
explicit array, and an omitted segment `vertical_speed_mps` remains absent until resolution.
Resolved wind is exactly one of two object shapes: a constant object with concrete `speed_mps`,
`direction_from_rad`, and `vertical_speed_mps`, or a segmented object whose segments each have a
concrete vertical speed. Still air resolves to the constant shape with all three values set to
zero.

### `solver`, `effects`, and `sampling`

| Section and field | Default | Meaning |
| --- | --- | --- |
| `solver.method` | `rk45` | `rk45`, `rk4`, or `euler`. |
| `solver.time_step_s` | `0.001` | Fixed step for RK4 and Euler; accepted but ignored by adaptive RK45. |
| `effects.magnus` | `false` | Enable the engine's Magnus-force model. |
| `effects.coriolis` | `false` | Enable Earth-rotation deflection. |
| `effects.enhanced_spin_drift` | `false` | Enable enhanced spin-drift modeling. |
| `sampling.interval_m` | `10` | Regular downrange result interval. |

Effects remain opt-in. The service may require projectile length, twist data, latitude, or other
documented prerequisites when a corresponding effect is enabled.

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

`summary.termination` is one of `max_range`, `ground_threshold`, `time_limit`, or
`velocity_floor`. The solve service must populate it from explicit termination metadata returned
by the engine. It must not infer a reason heuristically from the last distance, height, or speed.

Regular interval sampling never drops the terminal observation. `samples` includes the terminal
sample exactly once even when its distance is not an interval boundary, and that sample carries
the `terminal` flag. Sample flags are `transonic`, `subsonic`, `terminal`, and
`ground_threshold`.

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

## Deliberate v1 exclusions

V1 does not expose custom drag files or tables, velocity/Mach-dependent BC schedules, powder
temperature curves, atmosphere zones, cluster-BC degradation, wind shear, pitch damping,
precession/nutation, or angular diagnostics. It also does not expose arbitrary filesystem or
network access. Those features require a later schema or a separately versioned extension after
their input semantics and result provenance are stable.
