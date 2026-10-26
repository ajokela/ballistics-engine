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

v1 does, however, grow by ADDITION under four strict rules, so that a feature need not force a
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
4. A combination of individually valid request fields that was refused as `conflicting_fields`
   may become accepted when it names newly supported behavior. Every request valid before keeps
   identical behavior and output. Producers that must work with older engines must not send the
   combination, since those engines still refuse it (`shot.zero_distance_m` beside
   `shot.muzzle_angle_rad` was accepted this way in 0.33.0, and a `reference`-kind
   `projectile.drag_table` beside `projectile.bc_segments` after 0.46.0, a pair that engines up
   to 0.46.0 refuse at `$.projectile.bc_segments`).

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
    "temperature_k": 291.15,
    "pressure_pa": 98500.0,
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

The `atmosphere` block here is a measured station reading: `pressure_pa` is the absolute pressure
at the 250 m station, a little above the ICAO standard there (98357.65 Pa), not a
sea-level-corrected value. Send `pressure_pa` and `temperature_k` only when you have a reading.
When you do not, omit them: the engine resolves the ICAO standard at `altitude_m` and says so in
`assumptions` (`icao_standard_pressure`, `icao_standard_temperature`). Do not send the sea-level
standards `101325 Pa` and `288.15 K` as placeholders. A present value is authoritative, so
`101325 Pa` at 250 m describes a day of unusually high pressure (about 1044 hPa at sea level), and
with a temperature also present `altitude_m` no longer sets the air density. Higher up, above
about 645 m, the placeholder also draws a `station_pressure_implausible_for_altitude` warning
(see [`atmosphere`](#atmosphere)); at 250 m it does not, because 101325 Pa is a possible reading
there. An altimeter setting (QNH), such as a METAR figure or a weather meter's baro reading, goes
in with `"pressure_reference": "qnh"`; see [`pressure_reference`](#pressure_reference-mba-1397),
and [ATMOSPHERE_DENSITY_ALTITUDE.md](ATMOSPHERE_DENSITY_ALTITUDE.md#request) for other pressure
sources.

### `projectile`

| Field | Required | Meaning |
| --- | --- | --- |
| `mass_kg` | yes | Projectile mass. |
| `diameter_m` | yes | Projectile diameter. |
| `length_m` | no | Projectile length; required by effects that need geometry. |
| `drag_model` | yes | One of `G1`, `G2`, `G5`, `G6`, `G7`, `G8`, `GI`, `GS`, or `RA4`. Does not drive drag while `drag_table` is present. `bc_segments` band BCs are referenced to it, or to a `drag_table` of kind `reference` when one is sent. |
| `ballistic_coefficient` | yes | BC for the selected reference drag model, or for `drag_table` when its `kind` is `reference`. Still required beside `bc_segments`, but then it flies only in an interior gap between bands, and blended into the bands within up to 25 fps either side of a gap edge. |
| `drag_table` | no | A caller-supplied drag curve replacing the built-in model; see below. |
| `bc_segments` | no | A velocity-banded BC ladder, the wire form of the CLI's `--bc-segment`; see below. |

All nine built-in reference drag models are backed by distinct tables and are accepted by v1.
The enum spellings are exact and case-sensitive.

### `projectile.drag_table` (MBA-1597)

A Cd-vs-Mach curve sent inline, the wire form of the CLI's `--drag-table`. Use it for a drag law
the engine does not ship — an airgun pellet law such as GA2, say — or for a projectile's own
measured curve.

```json
"drag_table": {
  "kind": "reference",
  "name": "Example (invented points)",
  "points": [
    {"mach": 0.0, "cd": 0.210},
    {"mach": 0.5, "cd": 0.200},
    {"mach": 0.9, "cd": 0.350}
  ]
}
```

| Field | Required | Meaning |
| --- | --- | --- |
| `kind` | yes | `reference` or `projectile`, below. There is no default. |
| `name` | no | A label for the curve, such as `"GA (ChairGun Elite)"`: 1 to 64 characters (Unicode scalar values), not blank, one line (no control characters and no U+2028/U+2029 separators). Echoed verbatim and never read by the solver. |
| `points` | yes | 2 to 4096 `{mach, cd}` objects. Mach finite, at least 0 and strictly ascending; Cd finite and above 0. |

`kind` decides what `ballistic_coefficient` means, and there is no default because the two
readings differ by the projectile's form factor — a wrong guess gives confident, wrong drops:

- **`reference`** — a standard curve that `ballistic_coefficient` was measured against (a GA2 BC
  with the GA2 curve). The pair flies exactly as a built-in model flies with its own BC: sending
  the engine's own G1 table as a `reference` curve with a G1 BC gives the same trajectory as
  `"drag_model": "G1"` to rounding. The BC still matters, so the error budget and the other
  kernels that perturb it move the result as they do for a built-in model. A
  [`bc_segments`](#projectilebc_segments-mba-1598) ladder measured against the same curve can be
  sent beside it; see **With a BC ladder** below.
- **`projectile`** — this projectile's own measured drag coefficient (Doppler radar, for
  instance). `ballistic_coefficient` is not used, exactly as with `--drag-table`; the engine
  divides the curve's Cd by the sectional density from `mass_kg` and `diameter_m`. Since the BC
  is unused, perturbing it moves nothing, and a `bc_segments` ladder is refused beside it.

Outside the table the nearest end value is held, so give a curve covering every Mach the flight
reaches. With either kind, `drag_model` is still required and echoed but has no effect on drag.

`name` is for people reading a request or a response later. `drag_model` cannot say which curve
flew, and two programs can ship different curves under the same name, so say it here, source
included. It changes nothing in the result, and `explain` treats two tables that differ only in
name as the same curve. When it is omitted, nothing about the table's echo changes. Engines up to
0.47.1 refuse it as `unknown_field` at `$.projectile.drag_table.name`, so a producer that must
also work with older engines must not send it.

Errors name the exact entry: `$.projectile.drag_table.kind` for a missing or unknown kind,
`$.projectile.drag_table.name` for a name that is not a string, is blank, holds a control
character or line break, or runs past 64 characters, `$.projectile.drag_table.points` for a wrong point count,
and
`$.projectile.drag_table.points[3].mach` (or `.cd`) for a bad value. A table cannot be combined
with [`corrections.bc5d_table_path`](#optional-corrections-block): BC5D corrects a G1/G7 BC and a
table replaces the curve, so the pair is `conflicting_fields` at
`$.corrections.bc5d_table_path`. A `projectile`-kind table also cannot be combined with
[`projectile.bc_segments`](#projectilebc_segments-mba-1598): its curve is the projectile's own
Cd, with no BC for the bands to change, so that pair is `conflicting_fields` at
`$.projectile.bc_segments`. A ladder that breaks one of its own rules, including a one-band
ladder whose BC does not match the scalar, is refused at that band first, by the decoder, beside
either kind.

**With a BC ladder.** A `reference` table and a `bc_segments` ladder may be sent together. Use
the pair when BCs are published per velocity band against a curve the engine does not ship,
such as an airgun pellet law. The table replaces `drag_model`'s curve for the bands as it does
for the scalar, and the bands are selected, gapped and blended exactly as on a built-in model:
at each airspeed the engine divides the table's Cd by the band BC that flies there, or by the
scalar `ballistic_coefficient` in an interior gap. Sending the engine's own G7 table as a
`reference` curve beside a ladder gives the same samples as `"drag_model": "G7"` with the same
ladder, under every solver method. The band BCs must have been measured against this table's
curve. The engine cannot check that, and bands quoted against G1 but flown over a pellet law
give confident, wrong drops. 0.46.0 refuses the pair as `conflicting_fields` at
`$.projectile.bc_segments`, and earlier engines refuse `bc_segments` itself as `unknown_field`,
so a producer that must also work with older engines must not send it
([versioning rule 4](#versioning-and-compatibility)).

The table is echoed verbatim at `resolved_request.projectile.drag_table`, beside any ladder at
`resolved_request.projectile.bc_segments`, so a resolved request re-solves on the same curve and
bands. The engine ships no airgun curves: the caller supplies the table.

### `projectile.bc_segments` (MBA-1598)

A velocity-banded BC ladder sent inline, the wire form of the CLI's `--bc-segment`. Use it to fly
a published multi-BC schedule, where a maker quotes one BC per velocity band, or BCs you measured
at different speeds.

```json
"projectile": {
  "mass_kg": 0.01089,
  "diameter_m": 0.00782,
  "drag_model": "G1",
  "ballistic_coefficient": 0.462,
  "bc_segments": [
    {"velocity_min_mps": 640.0, "velocity_max_mps": 1200.0, "ballistic_coefficient": 0.462},
    {"velocity_min_mps": 490.0, "velocity_max_mps": 640.0, "ballistic_coefficient": 0.447},
    {"velocity_min_mps": 0.0, "velocity_max_mps": 490.0, "ballistic_coefficient": 0.424}
  ]
}
```

The member is an array of 1 to 64 bands (`MIN_BC_SEGMENTS_V1` and `MAX_BC_SEGMENTS_V1` in
`ballistics_engine::solve_json`). Each band is an object with exactly these three members, all
required numbers:

| Field | Meaning |
| --- | --- |
| `velocity_min_mps` | The band's lower edge, in m/s, inclusive. Finite and at least 0. |
| `velocity_max_mps` | The band's upper edge, in m/s, exclusive. Finite and above `velocity_min_mps`. There is no upper bound. |
| `ballistic_coefficient` | The BC that flies inside the band. Same units as the scalar `ballistic_coefficient` (lb/in²), referenced to the same curve (`drag_model`, or a `reference`-kind `drag_table`), and held to the same bounds, `0.0001` to `100`. |

**Order.** List the bands from the fastest to the slowest. `velocity_min_mps` must strictly
descend, and each band must end at or below the start of the band before it:
`bc_segments[i].velocity_max_mps <= bc_segments[i-1].velocity_min_mps`. Bands may touch (equal
edges) or leave a gap between them, but they may not overlap. The order is only a validation
rule: the engine's lookup does not depend on it, and the rule leaves exactly one legal order for
each set of bands.

**Which BC flies.** Before the edge blending below, a band covers the speeds from its
`velocity_min_mps` up to, but not including, its `velocity_max_mps`, so a touching edge belongs to
the faster band. At any speed:

- inside a band, that band's BC flies;
- in an interior gap between two bands, the scalar `ballistic_coefficient` flies;
- below the slowest band, or at or above the top of the fastest, the nearest edge band applies.
  The ladder is clamped at its ends, not padded with the scalar.

`ballistic_coefficient` therefore stays required and validated, but under a gapless ladder (the
example above) it never flies. A one-band ladder flies that band's BC at every speed, so its BC
must match `ballistic_coefficient` (see **Errors** below), and its two edges change nothing as
long as they are valid. A fastest band that ends below the muzzle velocity still
applies from the muzzle down. Band BCs are referenced to the same curve as the scalar:
`drag_model`'s, so a G7 ladder needs `"drag_model": "G7"`, or a `reference`-kind `drag_table`'s
when one is sent (see [`projectile.drag_table`](#projectiledrag_table-mba-1597), **With a BC
ladder**).

**Floor-style tables.** Many calculators take one speed per BC row, meaning "this BC from this
speed up to the next faster row", and show a single row with a speed of 0. Such a table maps onto
bands directly. List the rows fastest first, and give each band the row's speed as
`velocity_min_mps` (0 for the last row) and the next faster row's speed as `velocity_max_mps`. The
fastest band's `velocity_max_mps` only has to be valid, because the ladder is clamped above it. A
value at least 200 fps (60.96 m/s) above its `velocity_min_mps` keeps that band from narrowing the
blend at the edge below it (see **Edges are blended**). The CLI applies the same floor rule to
the BC tables in saved profiles, ending the fastest band just above 5000 fps (1524 m/s). Send a
one-row table as one band from 0 to 1524 m/s with the row's BC as `ballistic_coefficient` too, or
simply send the scalar alone: both fly the same.

**Airspeed, not ground speed.** Bands are selected by airspeed, the bullet's speed relative to the
air mass. That is the same speed the drag is computed from, so in a wind the band edges are
crossed at a different distance than in still air.

**Edges are blended.** The BC does not jump at an edge. The engine blends every interior
band-to-band edge, and both edges of every interior gap, with a smoothstep (`3t^2 - 2t^3`). The
blend window is centred on the edge and is as wide as the smaller of two values: 50 fps
(15.24 m/s), or a quarter of the narrower of the two spans that meet there (a band or a gap). So
the BC starts to change at most 25 fps (7.62 m/s) either side of an edge, and closer to it where a
band or gap is narrow. The two outer edges of the ladder are clamps, with nothing on the far side
to blend toward, so the BC stays flat across them. A one-band ladder is never blended. This is the
engine's own lookup (`velocity_segment_bc` in `src/bc_estimation.rs`), the same one
`--bc-segment` and `card.*` use.

**Errors.** The decoder checks every rule above before anything is solved, and the service checks
the value rules again for a Rust caller that builds the request without the decoder. Each error
names the exact entry:

- `$.projectile.bc_segments`: `invalid_value` for `null`, a non-array, an empty array, or more than
  64 bands. To send no ladder, omit the member; `[]` is not read as "no ladder".
- `$.projectile.bc_segments[2]`: `invalid_value` for a band that is not an object.
- `$.projectile.bc_segments[2].<member>`: `missing_field` for a missing member, `unknown_field` for
  a member not listed above, and `invalid_value` for a non-number.
- `$.projectile.bc_segments[2].velocity_min_mps`: `invalid_value` when it is negative or not
  finite, or when it is not below the previous band's. The order message names the previous band
  and says to list the bands from the highest speed to the lowest, or to reverse the array.
- `$.projectile.bc_segments[2].velocity_max_mps`: `invalid_value` when it is not above
  `velocity_min_mps` or not finite, or when it reaches above the start of the previous band. The
  overlap message names that band.
- `$.projectile.bc_segments[2].ballistic_coefficient`: `invalid_value` outside `0.0001` to `100`,
  in the same form as the scalar's message.
- `$.projectile.bc_segments[0].ballistic_coefficient`: `conflicting_fields` when the ladder has
  exactly one band and its BC differs from the scalar `ballistic_coefficient` by more than one
  part in a million of the larger (`ONE_BAND_BC_RELATIVE_TOLERANCE_V1`). A one-band ladder flies
  its band at every speed, so the two would be two different BCs for one projectile. The
  tolerance forgives rounding nobody chose, such as a BC held as a 32-bit float or converted
  between units and back, and within it the band's BC flies; 0.172 against 0.1724 is still
  refused. The message names both values. Send the same value in both, or omit the ladder.

The values are checked band by band in array order, and the first fault is reported. A band
nested inside the one before it starts above that band's start, so the order rule reports it, at
its `velocity_min_mps`. The one-band check runs after the band's own value rules.

**Conflicts.** A ladder cannot be combined with either of these. The service refuses both as
`conflicting_fields` at `$.projectile.bc_segments`, before any correction table is opened:

- a [`projectile.drag_table`](#projectiledrag_table-mba-1597) of kind `projectile`: that curve is
  the projectile's own Cd, divided by its sectional density, so no BC reaches the drag and the
  bands would change nothing. A `reference`-kind table is accepted, and the band BCs are then
  measured against its curve instead of `drag_model`'s.
- [`corrections.bc5d_table_path`](#optional-corrections-block): the correction table builds its own
  BC ladder from `ballistic_coefficient`. Send one ladder or the other.

A request that carries a ladder, a drag table of either kind and a correction table together
fails the older drag-table-with-correction-table check first, so it is `conflicting_fields` at
`$.corrections.bc5d_table_path` instead. All of these conflicts are checked after the ladder's own
rules: a band that breaks one, including a one-band BC that does not match the scalar, is refused
at that band first, by the decoder.

**Echo.** The ladder is echoed verbatim at `resolved_request.projectile.bc_segments`, so a request
rebuilt from a resolved request re-solves on the same ladder. The zero search and the trajectory
both fly it. Omitting the member changes nothing, and the key never appears in a response.

**Known limitation: BC sensitivity.** The CLI's `error-budget` and `tolerance` perturb the scalar
`ballistic_coefficient`, and no input axis carries the band BCs of a ladder with two or more
bands. A one-band ladder is the exception: its BC must equal the scalar, so the BC axis moves the
band with it, and the slope is exactly that of the same request without the ladder. Under a
gapless ladder of two or more bands the scalar never flies, so its slope is exactly 0, and `error-budget` reports a declared BC sigma as
contributing nothing rather than listing BC under `unavailable_sources`. This is the same blind
spot as a `kind: "projectile"` drag table. With an interior gap, the scalar moves the result only
through the gap and the band speeds within up to 25 fps of a gap edge, where it is blended into
the bands. Beside a `reference`-kind drag table all of this holds: the slopes are exactly those
of the built-in model of the same curve with the same ladder. The one difference is in the last
digits: there a one-band ladder matches the same request without a ladder only to rounding,
because a reference table without a ladder keeps its older scaled path. `explain` is not
affected: its `ProjectileDrag` swap carries the whole ladder, and any drag table, with
`drag_model` (MBA-1603).

**Platforms.** The ladder is plain request data and no file is read, so it works wherever
`solve` runs, including through the WASM module's `bridgeCall`, where `corrections.bc5d_table_path`
is refused for want of a filesystem.

**Compared with `card.*`.** The bridge's `card.*` commands take a BC ladder too, under different
rules, so a ladder stored for one is not always valid for the other:

| | `solve`: `projectile.bc_segments` | `card.*`: `bc_segments` |
| --- | --- | --- |
| Band members | `velocity_min_mps`, `velocity_max_mps`, `ballistic_coefficient` | `velocity_min`, `velocity_max`, `bc` |
| Speed unit | m/s, always | the request's `velocity_unit` (fps imperial, m/s metric by default) |
| Order | fastest first, no overlaps | any order; no overlaps (touching edges and gaps are fine) |
| `[]` | refused | read as no ladder |
| With a correction table | `conflicting_fields` | the ladder wins and the table is ignored |
| Band BC | `0.0001` to `100` | `0.0001` to `100` |
| One band beside the scalar BC | must match it (within one part in a million) | not compared; the band flies |
| Band count | 1 to 64 | up to 64 |

Both are converted to the engine's fps ladder and fly through the same lookup. The CLI's
`--bc-segment` follows `card.*` in units, order and precedence: its speeds follow `--units`, any
order is accepted, and it overrides the correction tables. Unlike `card.*`, it still accepts
overlapping bands, more than 64 bands and any positive band BC.

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

#### Sight-in (zero-day) conditions

A request carries one atmosphere, one wind and one set of effects, and the zero search above runs
in them. A rifle is usually sighted in on another day, often somewhere else: a different
temperature, pressure, altitude or humidity. The CLI models that with `--zero-temperature`,
`--zero-pressure`, `--zero-pressure-type`, `--zero-humidity`, `--zero-altitude`,
`--zero-density-altitude` and `--zero-velocity`: it solves the zero in the sight-in air, in still
air, and then flies the trajectory in the shot-day air. solve-json has no second atmosphere. Send
two `solve` requests instead:

1. **The zero request** sends the same projectile and rifle, `shot.zero_distance_m` (and
   `shot.zero_poi_up_m` for a deliberately high zero), the sight-in atmosphere and still air
   (`"wind": {}`). Leave `pressure_pa` out when no pressure was measured at the sight-in, so the
   engine uses the ICAO pressure for that `altitude_m`. Keep it level and uncanted, as a zeroing
   range usually is: omit `shooting_angle_rad` and `cant_angle_rad`. If the muzzle velocity on the
   sight-in day differed, for example at another powder temperature, send that velocity here.
   Coriolis reaches the search, so with `effects.coriolis` on, give this request the sight-in
   `atmosphere.latitude_rad` and `shot.shot_azimuth_rad`. Read the solved angle from
   `result.resolved_request.shot.muzzle_angle_rad`.
2. **The shot request** sends the shot-day atmosphere, wind, geometry and effects, with that angle
   in `shot.muzzle_angle_rad`. Leave `zero_distance_m` out unless the request relies on the
   windage-convergence terms (`rifle.sight_offset_lateral_m`, `shot.zero_poi_right_m`), which need
   it. With both fields present, the `zero_distance_elevation_not_resolved` warning described
   above is expected.

The pinned angle carries the zero shift: the shot request flies the barrel angle that hit the zero
in the sight-in air. This reproduces the CLI exactly. Engine 0.48.0 output, for a .308 (11.34 g,
7.82 mm, G7 0.243, 792.48 m/s, sight height 0.0381 m) zeroed at 91.44 m at sea level, 288.15 K
and 101325 Pa, then shot at 3500 m, 303.15 K and 65760 Pa, with 50% humidity on both days and
all effects off:

| | Zero angle | Drop at 914.4 m |
| --- | --- | --- |
| Two `solve` requests as above | 0.00116729736328125 rad (0.066881212353°) | 8.1301 m |
| CLI `trajectory --units metric --auto-zero 91.44 --zero-altitude 0 --zero-temperature 15 --zero-pressure 1013.25 --zero-humidity 50` with the same inputs | 0.066881212353° | 8.13 m |
| One request, zero solved in the shot-day air | — | 8.1427 m |

Here the sight-in conditions move the solution by 0.047 MOA. One request with `zero_distance_m`
is the right model only when the rifle was zeroed in the air it is being shot in.

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
when `altitude_m` is also omitted). The same object, sent to the bridge's
`atmosphere.density_altitude` command, reports the pressure altitude, density altitude and air
density of the air it resolves to — see [ATMOSPHERE_DENSITY_ALTITUDE.md](ATMOSPHERE_DENSITY_ALTITUDE.md). Enabling Coriolis without a latitude is an `invalid_value`
error rather than a silently chosen latitude.

Explicit temperature and pressure values are authoritative station conditions, including values
equal to `288.15 K` and `101325 Pa` at nonzero altitude. The solve service must preserve that
explicit intent and must not apply the legacy CLI rule that treats those exact values as omitted
standard-atmosphere sentinels. MBA-1302 must bypass that CLI sentinel inference: omitted fields
select ICAO-at-altitude resolution, while present fields are passed as authoritative values. The
resolved values are recorded as explicit numbers in a successful response's
`resolved_request`; the resolved altitude and relative humidity are concrete there as well.

An absolute `pressure_pa` that cannot be a station reading at `altitude_m` is still flown as
sent, and the response carries a `station_pressure_implausible_for_altitude` warning at
`$.atmosphere.pressure_pa`. It is raised when `pressure_pa` is outside 0.85 to 1.08 times the
ICAO standard pressure at the resolved `altitude_m` (0 m when omitted). The recorded sea-level
pressure extremes, 870 hPa and 1089.1 hPa, are 0.859 and 1.075 times standard, so weather does
not raise it. Above the band the likely cause is a sea-level placeholder, or an altimeter setting
sent without `"pressure_reference": "qnh"`; below it, an `altitude_m` that is not the station
elevation, or a pressure that is not in pascals. The message names the causes for its side. In
engine output, 101325 Pa at Cusco's 3310.128 m is 1.503 times the standard there (67402.29 Pa)
and is warned. For a .308 (11.34 g, G7 0.243, 792.48 m/s, zeroed at 91.44 m, 16 °C, 52% RH) it
gives 10.811 m of drop at 914.4 m, against 8.368 m from a real reading of 68830 Pa, which is not
warned. The same 68830 Pa with `altitude_m` omitted is 0.679 times the 101325 Pa standard at
0 m and is warned, though with a temperature present the altitude only sets how the air changes
along the bullet's path, and the drop at 914.4 m differs from the 3310.128 m request by 0.3 µm.
The band catches gross mismatches only: 101325 Pa stays inside it up to about 645 m. The warning
is never raised for a `"qnh"` pressure or an omitted one, and it changes no number in the
response.

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

The height driving the profile is height above the *ground*, not above sea level: the projectile's
height above the ground reference `rifle.muzzle_height_m` is measured from. When
`rifle.muzzle_height_m` is omitted or `0` the bore's real height is unknown, so the bore is taken
to be 1.5 m above the ground. On an inclined shot the ground is taken to run along the line of
sight, as it does between a shooter and a target both standing on a slope and as
`shot.ground_threshold_m` already assumes, so the height is measured perpendicular to the line of
fire. That differs from `summary.maximum_height_m`, which reports a world-vertical height.
`atmosphere.altitude_m` is deliberately not an input to it: boundary-layer shear is relative to
the local ground, and altitude's effect on the shot is air density, which `atmosphere.altitude_m`
already carries.

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
from an assumed barrel rather than disabled — see `effects.aerodynamic_jump` above),
`spin_effect_assumed_twist_rate` (`effects.magnus` or `effects.enhanced_spin_drift` enabled without
`rifle.twist_rate_m_per_turn`, so the spin-driven model was computed from the assumed 1:12 twist
rather than disabled), `stability_factor_assumed_twist_rate` (`rifle.twist_rate_m_per_turn`
omitted on a solve that reported `summary.stability_factor`, so the Sg published is the Sg of the
assumed 1:12 barrel — this one is not gated on any flag, because Sg is computed on every solve; see
[Omitting the twist rate](#omitting-the-twist-rate) above for both), and
`station_pressure_implausible_for_altitude` (an absolute `atmosphere.pressure_pa` outside 0.85 to
1.08 times the ICAO standard pressure at `altitude_m`, flown as sent — a sea-level placeholder,
an altimeter setting sent as station pressure, or a wrong `altitude_m`; see
[`atmosphere`](#atmosphere) above). More than one of these can be
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
No current solver emits `velocity_floor`; the value is reserved.

A projectile that loses all of its downrange speed, such as one at the bottom of the accepted
`ballistic_coefficient` range, a `bc_segments` band that low, or a long fall toward a very low
`ground_threshold_m`, ends like any other early termination (MBA-1599). The response is
`status: ok`, `termination` is `ground_threshold` (or `time_limit` if it is still falling at
100 s), and `actual_range_m` is the distance at which downrange motion stopped, with the terminal
sample at the impact state. Every integration step counts against the engine's 250,000-point
limit, including the steps spent at the stall distance, so with a fine `euler` or `rk4`
`time_step_s` a stall can reach that limit before the ground and return `solve_failed` with
the message `trajectory point limit of 250000 exceeded`. A trajectory that moves back toward
the muzzle before it terminates (a headwind carrying a stalled projectile back, a launch past
vertical, or an explicit `euler`/`rk4` time step too coarse for the drag) cannot be reported
against downrange distance and returns `solve_failed`.

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

## Deriving holds from samples

`samples` say where the bullet is, not where to aim. Turning a sample into a hold takes a sign
and one division, and the engine's own hold outputs (`reticle_hold`, the bridge's `card.*`
commands, `ballistics lead`) do exactly that arithmetic, so an app that repeats it on the same
drop gets their numbers.

**Signs.** `drop_m` is positive below the line of sight, so a positive `drop_m` is a hold or a
dial up, and a negative one, with the bullet above the line of sight, is a hold down.
`windage_m` is positive to the shooter's right, so a positive `windage_m` is a hold left and a
negative one a hold right.

**Angles.** The engine converts with the small-angle ratio, never an arctangent:

```text
elevation_mil = drop_m    / distance_m * 1000    (positive: hold up)
windage_mil   = windage_m / distance_m * 1000    (positive: hold left)
```

`reticle_hold` is this arithmetic at `reticle.range_m`: its `down_mil` is `elevation_mil` and its
`right_mil` is `windage_mil`, unchanged. Both are the reticle point to put on the target,
measured from the optical center: a bullet that drops is held over with a point below center,
and a bullet blown right is held left with a point to the right of center, so neither sign
needs flipping. The `card.*` commands keep the same signs in their angular columns: a positive
`drop_adj` is a hold up, and `wind_adj` has the sign of `windage_m`.

**MOA.** The engine's MOA figures do not all use one constant:

| Output | MOA per radian |
| --- | --- |
| `card.come_ups`, `card.range_table`, `card.wind` and `card.pdf` columns in `"moa"`, the PDF's Lead column included | 3438 |
| The same columns in `"smoa"` or `"iphy"` | 3600 (one inch per 100 yards) |
| `ballistics lead` and the `lead_from_tof` formula below | 3438 |
| `reticle.catalog`'s MOA ladders, and MOA clicks in `true.dial_plan` | 3438 |
| `summary.aerodynamic_jump_moa`, which the solver divides by this to get its launch-angle offset | 10800/π = 3437.7467707849 |
| `true.fit` and `true.dsf` drops sent and returned with `drop_unit` `"moa"` | 10800/π |

3438 is the engine's printed-table constant. An MOA column that is meant to agree with
`card.range_table` converts with it, `moa = drop_m / distance_m * 3438`, which is
`elevation_mil * 3.438`. The exact 10800/π gives values 0.0074% smaller, which is enough to
change the last printed digit on some rows.

**Lead.** The bridge has no lead command. `ballistics lead` and the `card.pdf` Lead column
compute lead with `lead_from_tof` (`src/moving_target.rs`), and an app can apply it to a
sample:

```text
lead_m   = target_speed_mps * sin(angle) * time_s
lead_mil = lead_m / distance_m * 1000
lead_moa = lead_m / distance_m * 3438
```

`angle` is the direction the target travels relative to the line of sight: 0° directly away,
90° crossing left to right, 180° directly toward the shooter, 270° crossing right to left. A
positive lead is a hold to the right and a negative one a hold to the left, so either way the
hold is ahead of the target. The ratio is the same small-angle one as above, and it is zero
under 1 m. Lead is target motion only: the wind is already in `windage_m`, and the lead adds to
it on the same axis, so the whole horizontal hold, positive left as above, is
`windage_mil - lead_mil`.

`ballistics lead` also moves the range for a target that is not crossing square. It iterates
the intercept range `R = R0 + target_speed_mps * cos(angle) * TOF(R)` until a step changes it by
less than 0.1 m, and evaluates the lead there. `lead_from_tof` on the sample at `R0` does not,
so the two agree only for a crossing target (90° or 270°). An app that needs the correction can
iterate the same way, interpolating `time_s` between samples.

### Worked example

The representative request above, with `effects` all `false` (the `card.*` commands and
`ballistics lead` model neither Coriolis nor spin drift), has this sample at 300 m:

```json
{"distance_m": 300.0, "time_s": 0.4100179737381975, "speed_mps": 650.9376231518361, "energy_j": 2402.491204959965, "drop_m": 0.4583304990888585, "windage_m": -0.20339285743668, "mach": 1.8998340035015642, "flags": []}
```

- **Elevation:** `0.4583304990888585 / 300 * 1000` = 1.527768 mil, held up, or 5.252468 MOA.
  The same request with a `reticle` block at `range_m` 300 returns `reticle_hold.down_mil`
  1.527768330296195.
- **Windage:** `-0.20339285743668 / 300 * 1000` = −0.677976 mil. The wind from the right
  carries the bullet left, so the hold is 0.677976 mil right, or 2.330882 MOA; `reticle_hold`
  gives `right_mil` −0.6779761914556001.
- **Against `card.range_table`:** the same load and air as a metric card,

  ```json
  {"units": "metric", "muzzle_velocity": 823.0, "ballistic_coefficient": 0.243, "mass": 11.34,
   "diameter": 6.71, "drag_model": "g7", "sight_height": 38.1, "zero_distance": 100.0,
   "temperature": 18.0, "pressure": 985.0, "humidity": 50.0, "altitude": 250.0,
   "wind_speed": 4.4704, "wind_direction_deg": 90.0, "start": 100.0, "end": 500.0,
   "step": 100.0, "adjustment_unit": "moa"}
  ```

  prints `wind_adj` −2.330882146224353 at 300 m, which is `windage_m / 300 * 3438` from this
  sample to the last digit. Its `drop_adj` is 5.252455654337169, because its own solve puts the
  drop at 458.329 mm where this sample has 458.330 mm. This sample's drop with the exact
  constant would give 5.252081.
- **Lead:** a target crossing left to right (`angle` 90°) at 3 m/s has
  `lead_m = 3 * 0.4100179737381975` = 1.230054 m, so `lead_mil` = 4.100180 and `lead_moa` =
  14.096418, held right. With the wind, the horizontal hold is `-0.677976 - 4.100180` =
  −4.778156 mil, which is 4.778156 mil right.
- **Against `ballistics lead`:** with the same load, air and target,

  ```text
  ballistics lead -u metric -v 823 -b 0.243 --drag-model g7 -m 11.34 -d 6.71 \
    --sight-height 38.1 --temperature 18 --pressure 985 --humidity 50 --altitude 250 \
    --wind-speed 4.4704 --wind-direction 90 --target-speed 3 --target-angle 90 \
    --start 300 --end 300 -o json
  ```

  reports `lead_mil` 4.100183254173634 and `lead_moa` 14.096430027848953, over a `tof_s` of
  0.41001832541736344. It flies its own trajectory from a level bore, with no zero, so its time
  of flight differs from the sample's 0.4100179737381975 s in the seventh significant digit,
  and the lead with it.

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

V1 does not read drag files (a drag curve is sent inline as
[`projectile.drag_table`](#projectiledrag_table-mba-1597)) or Mach-keyed BC tables (a
velocity-keyed BC ladder is sent inline as
[`projectile.bc_segments`](#projectilebc_segments-mba-1598)). It has
no field for powder temperature curves, atmosphere zones, cluster-BC degradation, custom wind-shear
layer tables (the built-in profiles are
[`effects.wind_shear_model`](#effectswind_shear_model)), pitch damping, precession/nutation, or
angular diagnostics. (The optional [`corrections`](#optional-corrections-block) block is the one
place the engine builds a BC ladder itself, from a verified BC5D table, instead of reading one
from the request.) It also does not expose network access or filesystem access beyond reading the
one named correction table. Those features need a later schema, a separately versioned
extension, or an optional addition under the rules in
[Versioning and compatibility](#versioning-and-compatibility), once their input semantics and
result provenance are stable.
