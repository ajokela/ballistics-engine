# Ballistics Engine CLI Tool

Comprehensive command-line interface for professional ballistics trajectory calculations with advanced drag modeling and automatic zeroing.

## Installation

```bash
# Build from source
cargo build --release

# Binary location
./target/release/ballistics
```

## Unit Systems

The CLI supports two unit systems, selectable with the `--units` flag (default: Imperial)

### Imperial Units (Default)
- Velocity: feet per second (fps)
- Mass: grains
- Distance: yards
- Diameter: inches
- Sight Height: inches
- Bore Height: feet
- Temperature: Fahrenheit
- Pressure: inHg

### Metric Units
- Velocity: meters per second (m/s)
- Mass: grams
- Distance: meters
- Diameter: millimeters (mm)
- Sight Height: millimeters (mm)
- Bore Height: meters
- Temperature: Celsius
- Pressure: hPa

## Turret Adjustment Units

Separate from `--units` (imperial/metric, above), every sweep-table command takes an
`--adjustment-unit` flag selecting how angular dial/hold columns (Drop, Wind, Lead, the
mover Ring) are displayed. Five values are accepted (MBA-724, MBA-1355, MBA-1410):

| Value | Meaning | Conversion |
|---|---|---|
| `mil` (default) | Milliradians | `(drop_yd / range_yd) × 1000` |
| `moa` | True Minutes of Angle | `(drop_yd / range_yd) × 3438` — the CLI's locked printed-table dial constant, deliberately not the geometrically exact 3437.7467 (MBA-724) |
| `smoa` | Shooter's MOA (exactly 1 inch per 100 yards) | `(drop_yd / range_yd) × 3600` |
| `iphy` | Inches per hundred yards | numerically identical to `smoa` — same conversion, different header text |
| `clicks` | Whole turret clicks | see below — requires a click graduation, not a fixed factor |

Which columns the unit reaches depends on the command. `come-ups`, `wind-card`,
`range-table`, and `compare` always print dial columns, so it always applies there. On
`trajectory` it reaches only two places — the PDF dope card's Drop/Wind/Lead columns
(`-o pdf`) and the mover Ring column (`--target-speed`) — so a `trajectory` run with
neither leaves the printed table byte-identical. Since 0.28.2 that case warns on stderr
rather than passing silently; the browser terminal warns in the table output (it has no
PDF card and no stderr).

### Independent elevation vs. windage units (MBA-1410)

`--adjustment-unit` sets the **elevation** axis (Drop, and the mover Ring column, which
is treated as elevation-like). Commands with a separate horizontal column — `range-table`
and `compare` (Wind), plus `trajectory`'s PDF dope card (Wind/Lead) — also accept an
independent:

- **`--windage-unit <mil|moa|smoa|iphy|clicks>`** — falls back to `--adjustment-unit` when
  omitted, so e.g. `--adjustment-unit mil --windage-unit moa` prints mil elevation holds
  next to moa windage holds in the same run. `clicks` is only accepted here when
  `--adjustment-unit` is **also** `clicks` — turret click output needs the elevation
  graduation resolved first (see below); an explicit `--windage-unit clicks` paired with a
  non-clicks `--adjustment-unit` is rejected:

  ```
  error: --windage-unit clicks requires --adjustment-unit clicks (turret click output needs the resolved elevation graduation)
  ```

Commands with only ONE angular column resolve it through whichever axis the column
actually is, with no separate `--windage-unit` flag to set: `wind-card`'s drift and
`lead`'s hold are both windage-type (horizontal), so `--adjustment-unit` on those two
commands drives the windage axis directly. `come-ups` is the mirror case — a pure
elevation table, so `--adjustment-unit` there is always the elevation axis (see the
`--windage-click-value` note below).

### Scope tracking correction factors (`--elevation-cf`, `--windage-cf`, `tall-target`) — MBA-1358

Real turrets rarely move exactly their nominal click value. A **tall-target test**
(Litz, *Accuracy and Precision for Long Range Shooting*; the same feature AB, Strelok
Pro, Shooter, TRASOL, ColdBore and ChairGun expose) measures the error. The engine
stores the published tall-target ratio — **`CF = actual measured travel / dialed
travel`** (0.95 = the scope under-tracks by 5%) — and **divides** every dial-unit
output by it exactly once, at the shared conversion boundary. The direction is the
physics: to land a true angular correction `N` on a scope whose dialed unit delivers
`CF` true units, the number set on the dial must be `N / CF`, so an under-tracking
scope (CF < 1) needs **more** dial.

- **`--elevation-cf <FACTOR>`** — corrects elevation dial values: `come-ups`,
  `range-table`/`compare` Drop columns, the PDF dope card's Drop column, and the
  `zero` command's MOA/mrad outputs (native and browser terminal).
- **`--windage-cf <FACTOR>`** — corrects windage/lead dial values: `wind-card` drift,
  `range-table`/`compare` Wind columns, `lead` (mover lead is a dialed quantity), the
  mover Ring column/fields (`--target-speed`), and the PDF Wind/Lead columns.

**Worked example**: you dial 10 mil in the tall-target test and measure only 9.5 mil of
actual group travel (34.2 in at 100 yd) — `CF = 9.5 / 10 = 0.95`. With
`--elevation-cf 0.95`, a solution needing a true 10 mil correction displays as
`10 / 0.95 = 10.53` mil: dialing 10.53 indicated on that scope moves the impact the
true 10 mil.

The factors apply to **dial units only** (mil/MOA/SMOA/IPHY and whole clicks — the
angle is corrected before click quantization): raw drop/drift **inches and degree
bore-angle echoes are never touched**. Both must be **strictly between 0.5 and 1.5** —
anything outside that band means the tall-target test went wrong, and both the flags
and the stored profile fields (`elevation_cf` / `windage_cf`, validated on load,
carried forward on re-save) reject it with a hard error naming the field. `1.0` (the
default) is byte-identical to not passing the flag. CLI flags override the profile
fields. (Convention note: Kestrel's Scope Cal *multiplies* its factor into the solution
because it stores the reciprocal, dialed/actual; this engine divides because it stores
the published actual/dialed — same physics, opposite bookkeeping.)

**Deriving the factor** — `ballistics tall-target` computes `CF = actual measured travel
/ dialed travel`, pure arithmetic over the same locked angular factors every dial column
uses (no trajectory solve):

```bash
# Dialed 10 mil; the group moved 36.0 inches at 100 yd (exactly 10 mil): CF = 1.0000
ballistics tall-target --dialed 10 --measured 36 --range 100

# The group moved 34.2 inches (9.5 mil): CF = 0.9500 — enter as --elevation-cf 0.95
ballistics tall-target --dialed 10 --measured 34.2 --range 100 --unit mil
```

**Truing converts in the opposite direction** (`true-velocity --elevation-cf`): dialed
observations (mil/moa drops; single-observation mode is always mil) are **multiplied**
by the CF before the fit consumes them — a dialed reading of `D` scope units delivers
`D × CF` true angular units — so scope tracking error is not baked into the trued
MV/BC; the report's dial-unit values are then shown back in scope units (÷CF). Linear
`in` drops are tape measurements and never scale. (The uncertainty-aware mode converts
its inputs the same way but reports in true, scope-corrected units — its report is a
versioned analysis document, and it says so when a CF is active.)

### `clicks`: whole-click output

`--adjustment-unit clicks` rounds the angular adjustment to the nearest **whole turret
click** instead of printing an angle — ties round away from zero, sign is preserved, and
ranges under 1 yard/meter are defined as zero adjustment (same short-range guard as
every other adjustment unit). It needs a **click graduation** — the angular size of one
click on your turret — which has no default and must come from one of two places:

- **`--elevation-click-value <SIZE><UNIT>`** / **`--windage-click-value <SIZE><UNIT>`** —
  CLI flags, e.g. `--elevation-click-value 0.25moa` or `--elevation-click-value 0.1mil`.
  The suffix is mandatory and selects the graduation's own base unit — `mil`, `moa`,
  `smoa`, or `iphy` (`iphy` is accepted as an alias for `smoa`, the identical unit); the
  magnitude must be a positive, finite number. `wind-card` and `lead` have only one
  (windage) column, so they take just `--windage-click-value` — there is no
  `--elevation-click-value` flag on those two.
- **A saved profile's `elevation_click` / `windage_click` fields** — see [Importing
  profiles (.a7p)](#importing-profiles-a7p) below for saved-profile basics; set the click
  fields with `profile save --elevation-click <SIZE><UNIT> --windage-click <SIZE><UNIT>`,
  validated with the same parser at save time so a profile can never store an invalid
  graduation. `compare` has no single `--profile` (it compares multiple loads, each
  optionally from its own saved profile), so its click graduations come only from the
  explicit CLI flags, never a per-load profile.

**Resolution order** (checked once, eagerly, before any calculation): an explicit CLI
flag beats the saved profile's field for that axis. **Windage falls back to elevation**
when neither `--windage-click-value` nor the profile's `windage_click` is set — most
turrets share one graduation between the two knobs (on `wind-card`/`lead`, which have no
elevation column of their own, the single click value still falls back from
`--windage-click-value` to the profile's `windage_click`, then its `elevation_click`).
**Elevation must resolve from at least one source** on the commands that have an
elevation column — clicks output has nowhere else to get a graduation from — so a run
with neither an elevation flag nor a profile elevation click fails fast:

```
error: --adjustment-unit clicks requires a turret elevation graduation: pass --elevation-click-value <SIZE><UNIT> (e.g. 0.25moa or 0.1mil), or save one on the profile with `profile save --elevation-click`
```

Where clicks resolves, the header/column suffix follows the same `(mil)`/`(moa)`
convention as every other unit — e.g. the mover Ring column reads `Ring(clicks)`, and
come-ups' Drop column reads `Drop (CLICKS)` — and the values print as whole integers
instead of a decimal angle. Drop/Ring use the **elevation** click graduation; Wind/Lead
(on `range-table`, `compare`, and the PDF dope card), plus `wind-card`'s drift and
`lead`'s hold, use the **windage** one. `come-ups` has no windage column, so its
`--windage-click-value` flag can never affect its output — since 0.29.0 it warns on
stderr rather than silently accepting it:

```
warning: --windage-click-value has no effect on come-ups (it has no windage column); did you mean --elevation-click-value?
```

The default `mil`/`moa`/`smoa`/`iphy` output is completely unaffected by any of this —
`clicks` (and `--windage-unit`) are strictly additive.

## Commands

### Trajectory Calculation

Calculate ballistic trajectories with advanced physics modeling:

```bash
# Basic trajectory (Imperial - default)
./ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308

# With automatic zeroing at 200 yards
./ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --auto-zero 200

# Metric units
./ballistics trajectory --units metric -v 823 -b 0.475 -m 10.9 -d 7.82

# Full example with environmental conditions
./ballistics trajectory \
  -v 2700          # Velocity (fps)
  -b 0.475         # Ballistic coefficient
  -m 168           # Mass (grains)
  -d 0.308         # Diameter (inches)
  --drag-model g7  # G7 drag model
  --auto-zero 200  # Zero at 200 yards
  --max-range 1000 # Max range (yards)
  --wind-speed 10  # Wind (mph)
  --wind-direction 90 # Wind from right
  --temperature 59 # Temp (°F)
  --pressure 29.92 # Pressure (inHg)
  --humidity 50    # Humidity (%)
  --altitude 5000  # Altitude (feet)
  --full          # Show all points
```

#### Advanced BC Options

```bash
# Enable velocity-based BC segmentation
./ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 \
  --use-bc-segments \
  --auto-zero 600
```

`--use-bc-segments`' automatic characteristic-based estimator (used when no explicit
`--bc-segment`/`--bc-table-dir` schedule is supplied) is **G1/G7 only** — a wider drag
model warns and is treated as G1 for the estimate, same as `true-velocity`/`plan-truing`.

### Custom Drag Tables

Supply a measured or manufacturer-published drag curve — Hornady CDM data, a Lapua/Doppler-radar-derived deck, or your own — instead of relying on a G1/G7 reference curve plus a single BC value. Available via `--drag-table <FILE>` on the `trajectory`, `zero`, and `monte-carlo` subcommands.

**CSV format:** two columns, `mach,cd`, one point per line.
- A single leading header row (e.g. `mach,cd`) is tolerated and skipped.
- Blank lines and lines starting with `#` are ignored.
- Mach must be strictly ascending, with at least 2 data points.
- Cd must be finite and greater than 0.
- **Mach-keyed only.** Velocity-keyed decks (e.g. raw Doppler output in fps/m/s) must be converted by you first: `mach = velocity / speed_of_sound` at the conditions the velocity was measured under.

**`.drg` format:** `--drag-table` also accepts a file with a `.drg` extension (case-insensitive), a small vendor text format used to distribute Doppler-radar-measured drag curves (e.g. Lapua's free QuickTARGET Unlimited downloads). It tolerates a leading name/description line, tab/comma/semicolon-separated fields, and either `(mach, cd)` or `(cd, mach)` column order (detected automatically); any error names the file, the `.drg` format, and the offending line number. Dispatch is purely by file extension — a `.drg`-suffixed file is never parsed as CSV, and every other extension (including plain `.csv` or none) always goes through the CSV format above, unchanged. The engine ships no vendor drag curves of its own and cannot download them for you (licensing): obtain a `.drg` or CSV file from your manufacturer/measurement source yourself and point `--drag-table` at the local copy.

**Worked example:**

```bash
cat > deck.csv <<'EOF'
mach,cd
0.5,0.220
0.8,0.230
1.0,0.520
1.2,0.480
1.5,0.400
2.0,0.330
2.5,0.300
EOF

./ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --drag-table deck.csv --max-range 500
```

`-b`/`--bc` may still be supplied but its value is ignored once `--drag-table` is set. (On `trajectory` it is optional and defaults to 0.5; on `zero` and `monte-carlo` it remains a required argument, though it is likewise ignored for drag when a table is active.) — the deck supplies Cd directly. `-m/--mass` and `-d/--diameter` remain **required** (grains/inches under imperial, grams/mm under `--units metric`): the engine divides the deck's Cd by the projectile's sectional density (derived from mass and diameter) in place of the usual BC-based retardation denominator.

**Precedence:** a custom drag table completely replaces the G1/G7 model and any BC. It also takes precedence over `--use-bc-segments` / `--bc-segment`; if both are supplied, the drag table wins and a warning is printed:

```
Warning: --drag-table and BC segments were both provided; the drag table takes precedence and BC segments are ignored.
```

**Out-of-range policy:** Mach numbers outside the table's measured domain **hold the nearest tabulated Cd** rather than extrapolating. On `trajectory`, if the shot's Mach range (muzzle to impact) extends beyond the table's domain, a coverage warning is printed:

```
Warning: shot Mach range [1.47, 2.42] extends beyond the drag table domain [0.80, 1.20]; the nearest tabulated Cd is held outside that range (approximate).
```

**Monte Carlo caveat:** when a custom drag table is active, `--bc-std` dispersion is a no-op — the table fixes Cd directly, so perturbing the (ignored) BC value has no effect on drag. Velocity, angle, and wind dispersion still vary normally.

**WASM:** the browser/Node build has no filesystem, so `--drag-table <FILE>` isn't a CLI flag there. Instead, call `wasm.loadDragTable(bytes)` with the raw bytes of the same `mach,cd` CSV (parsed by the identical `DragTable::from_csv_str`) before running a command; `wasm.hasDragTable()` reports whether one is loaded. With no file extension available to dispatch on, `loadDragTable` tries CSV first exactly as before, and only on CSV failure — if the text looks like a `.drg` deck — retries it through the same `.drg` parser the native CLI uses; if neither format parses, the error names both. Once loaded, the table is applied automatically to every `trajectory`, `zero`, `lead`, and `monte-carlo` run — including `lead`, which has no native `--drag-table` flag of its own — until a new table replaces it. See `loadDragTable`'s doc comment in `src/wasm.rs` for the exact contract.

#### Whole-Curve Drag Scale (`--cd-scale`)

For after-the-fact truing of a custom deck against chronograph/observed data, `--cd-scale <FACTOR>` multiplies every interpolated Cd by a constant factor — `Cd_used = table.interpolate(mach) * cd_scale` — the same mechanism a Hornady 4DOF "AFF" (Aerodynamic Fudge Factor) or an Applied Ballistics "CDF" (Cd Factor) applies to true a whole measured curve. This is distinct from `--bc-adjustment`, which trues a scalar BC on the G1/G7 path; a custom drag table has no BC to adjust, so `--cd-scale` is its equivalent. `1.0` is neutral (byte-identical to omitting the flag); the typical truing range is `0.90`-`1.10`. Works identically with both the CSV and `.drg` deck formats — it scales the interpolated Cd after either loader, so the source format doesn't matter.

Available on `trajectory`, `zero`, and `monte-carlo` — the same trio as `--drag-table` — and **requires** it: supplying `--cd-scale` without `--drag-table` fails before any solve, naming `--bc-adjustment` as the G1/G7 alternative:

```
error: --cd-scale requires --drag-table (for G1/G7 use --bc-adjustment instead)
```

A value far outside the typical truing range (outside `[0.5, 2.0]`) is still accepted — the engine's own gate is only finite and `> 0` — but warns once on stderr:

```
warning: --cd-scale 3 is far outside the typical truing range (0.90-1.10)
```

**WASM:** pass `--cd-scale <FACTOR>` as a terminal argument to `trajectory`, `zero`, or `monte-carlo` alongside a table loaded via `loadDragTable`; the pairing requirement and range warning are identical (the pairing failure surfaces as a rejected promise/`Err` instead of a process exit, and the range warning is prepended to the table-style output rather than printed to a separate stderr stream). `lead` also accepts `--cd-scale` in the WASM terminal (MBA-1411) — since it already applies a loaded table unconditionally (see the drag-table note above), a table trued via `--cd-scale` elsewhere needs the same scale here or its truing is lost; native `lead` has no `--drag-table`/`--cd-scale` of its own, so this is WASM-only, like the drag-table application itself.

#### BC Reference Standard (`--bc-reference`)

A published ballistic coefficient is not a pure property of the bullet — it also encodes an assumed reference air density. This engine's retardation math is calibrated to the **ICAO Standard Atmosphere** (the basis for `CD_TO_RETARD`), which most modern published BCs already assume. Some vendors — notably Sierra, Hornady, and Barnes for many bullets — instead publish BCs referenced to the older, denser **Army Standard Metro** atmosphere. Feeding an Army-Standard-Metro BC into this engine as if it were ICAO under-predicts drag by about 1.8%.

`--bc-reference <icao|army-standard-metro>` declares which one `--bc` (and any BC segments) is referenced to. `icao` is the default and is a complete no-op — omitting the flag, or passing `--bc-reference icao` explicitly, is byte-identical to every run before this flag existed. `army-standard-metro` multiplies the declared BC by a fixed ratio derived from the two reference densities (`ASM_DENSITY_LB_FT3 / ICAO_DENSITY_LB_FT3 ≈ 0.98237`) exactly once, at input normalization, before any retardation math runs — so the conversion applies uniformly to a scalar `--bc`, `--bc-segment` entries, Monte Carlo BC sampling, and a `zero`/`trajectory` auto-zero solve, with no risk of it being applied twice or missed on one code path.

```bash
# Sierra's published G1 BC for this bullet (0.475) is Army-Standard-Metro-referenced.
# The engine internally solves with 0.475 * 0.98237 ≈ 0.4666.
./ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --bc-reference army-standard-metro
```

Available on `trajectory`, `zero`, `monte-carlo`, and `profile save` (stored in the saved profile so it round-trips with the profile; profiles saved before this flag existed omit the key and load as `icao`).

**No effect together with `--drag-table`:** a custom drag deck supplies the projectile's actual Cd and divides by sectional density, not by a BC value, so there is no BC to reference-convert. Combining `--bc-reference army-standard-metro` with `--drag-table` is accepted (not an error) but warns once:

```
warning: --bc-reference army-standard-metro has no effect together with a custom drag table (--drag-table): the deck's Cd is divided by sectional density, not a BC value, so no BC-reference conversion applies
```

**`estimate-bc` outputs are always ICAO-referenced.** The fitted BC that `estimate-bc` (and any truing built on it) reports is searched directly against this engine's own ICAO-calibrated retardation math — it is never itself declared Army-Standard-Metro internally. Do not pass a fitted value back in with `--bc-reference army-standard-metro` on a later `trajectory`/`monte-carlo`/`zero`/`profile save` — that would double-convert it (dividing out the correction an ICAO-referenced value never needed).

**WASM:** pass `--bc-reference <icao|army-standard-metro>` as a terminal argument to `trajectory`, `zero`, or `monte-carlo`; behavior matches native exactly (same one-time normalization site), except the inert-with-`--drag-table` warning has no visible stderr in the browser terminal, so it is prepended to the table-style output text instead (same convention as the `--cd-scale` range warning above).

**FFI:** `FFIBallisticInputs` is an ABI-frozen `repr(C)` struct for iOS/Android consumers and gains no new field. Instead, call the standalone `ballistics_bc_for_reference_standard(bc, reference_standard)` (`reference_standard`: `0` = ICAO, `1` = Army Standard Metro) once on a raw BC value and write the result into `FFIBallisticInputs.bc_value` before calling any existing `ballistics_calculate_trajectory*`/`ballistics_calculate_zero_angle*`/`ballistics_monte_carlo*` export unchanged. This is a pure addition to the C ABI — existing callers that never call the new function are unaffected and need no recompile.

#### Pressure Reference: Absolute vs. QNH (`--pressure-type`)

A barometric pressure reading can mean two different physical quantities, and mixing them up is a classic source of vertical error: **station pressure** (the actual air pressure at the shooter's altitude) versus **QNH** (a sea-level-corrected altimeter setting — the kind of "barometer" value read off a weather report or METAR, which has been artificially inflated to what the pressure WOULD be at sea level). Kestrel, AB Mobile, Shooter, JBM, and Hornady 4DOF all let the user declare which one they are entering; this engine's `--pressure` previously only supported station pressure, so a QNH value entered at a real altitude silently over-stated air density and under-stated drop.

`--pressure-type <absolute|qnh>` declares which one `--pressure` is. `absolute` is the default and is a complete no-op — omitting the flag, or passing `--pressure-type absolute` explicitly, is byte-identical to every run before this flag existed. `qnh` reduces the declared pressure to station pressure at `--altitude` via the ICAO inverse-barometric formula, the exact inverse of the standard-atmosphere formula this engine already uses for altitude-derived pressure:

```
station = QNH * (1 - 0.0065 * h_geopotential_m / 288.15) ^ 5.25588
```

`h_geopotential_m` is the geopotential height ISA is defined against, not the geometric
altitude you type in; the engine converts before applying the formula. The difference is
small (1500 m geometric is 1499.65 m geopotential, worth about 0.04 hPa) but it is why
substituting your `--altitude` straight into the expression above will not reproduce the
engine's number exactly.

```bash
# A weather report gives 30.31 inHg (1026.5 hPa) QNH at a 5,000 ft (1524 m) firing position.
# Entered as absolute (the historical default) this over-states density; entered as qnh it is
# correctly reduced to ~25.22 inHg (~854.1 hPa) station pressure before the solve runs.
./ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 \
  --altitude 5000 --pressure 30.31 --pressure-type qnh
```

An omitted `--pressure` is unaffected by `--pressure-type`: the sea-level standard default (29.92 inHg / 1013.25 hPa) reduces to precisely the same ICAO-standard station pressure at any altitude under either mode, so there is nothing to reduce.

Available on **every subcommand that accepts `--pressure`** (MBA-1416): `trajectory` (`--pressure-type`, plus an independent `--zero-pressure-type` for the `--auto-zero` zero-day override — it defaults to the shot-day `--pressure-type` when not given, matching the existing "omitting all `--zero-*` flags reuses the shot-day values" contract), `zero`, `estimate-bc`, `true-velocity`, `plan-truing`, `mpbr`, `come-ups`, `lead`, `wind-card`, `stability`, `range-table`, `compare`, and `profile save`. On the calculator subcommands the mode reduces `--pressure` against that command's own `--altitude`, so `qnh` means the same thing everywhere; omitting it is exactly `absolute`, leaving every existing invocation byte-identical. **`profile save` records the mode, and since MBA-1417 `--saved-profile` applies it — but only to the pressure it describes.** The profile's atmosphere (temperature, pressure, humidity, altitude, pressure mode, density altitude) now loads as a SET. Precedence per value is: an explicit CLI flag beats the profile, which beats the built-in standard. The stored MODE, however, is honored only when the profile's stored pressure VALUE is the one in use — supply `--pressure` yourself, or a `--location` CSV pressure, and the stored mode is ignored because it describes a different reading. That pairing is the point: applying a mode saved on another day to a pressure you just typed was a measured 77 fps error at 300 yd. A stored `density_altitude` is stricter still, since it supersedes altitude and pressure outright: it applies only when the run supplies no pressure or altitude of its own. Profiles saved before these flags existed omit the keys and load as `absolute` with no density altitude, unchanged.

**Which pressure the mode describes (MBA-1421).** A `--pressure-type` you type on the command
line describes **whatever pressure that run actually uses**, including one supplied by a
`--location` CSV's `PRESSURE` column. There is exactly one pressure in the run for it to
describe, and the mode is your present-tense declaration about it:

```bash
# 24.90 inHg in the CSV, read as a sea-level altimeter setting and reduced to station pressure
ballistics trajectory ... --location sites.csv --site TestSite --pressure-type qnh
```

That is the opposite of how a profile-STORED mode behaves, and deliberately so. A stored mode
describes the value the profile stored alongside it; applying it to a different pressure is the
silent-error case 0.29.0 fixed. A typed mode has no other value it could be describing.

As of 0.30.1, `--location` REQUIRES `--site` (the half-pair was silently ignored before),
and an unreadable file or a missing site row is a hard error instead of a warn-and-continue
on the default atmosphere. The reverse stays legal: a standalone `--site` (like a standalone
`--profile-row`) is a PDF-label-only input with a pre-existing job. Note the row selector is **`--site <NAME>`**, not `--location-name` (which only overrides the
location string printed on a PDF header). `--location` without `--site` currently loads nothing,
silently — tracked as MBA-1425.

**Not yet available** on `estimate-bc`, `true-velocity`, `plan-truing`, `mpbr`, `come-ups`, `lead`, `wind-card`, `stability`, `range-table`, or `compare` — these subcommands still treat `--pressure` as absolute station pressure only. This is a tracked follow-up, not an oversight.

**WASM:** pass `--pressure-type <absolute|qnh>` (and `--zero-pressure-type` for the auto-zero override) as a terminal argument to `trajectory`; behavior matches native exactly (same reduction formula, same Authoritative-constructor bypass of the legacy default-sentinel heuristic).

**Monte Carlo / FFI:** the engine's `BallisticInputs.pressure` and `FFIAtmosphericConditions.pressure` fields have always meant absolute station pressure and continue to. `FFIAtmosphericConditions` is an ABI-frozen `repr(C)` struct with no new field. A QNH-aware FFI caller should call the standalone `ballistics_reduce_qnh_pressure(qnh_hpa, altitude_m)` once and write the result into `FFIAtmosphericConditions.pressure` before calling any existing `ballistics_calculate_trajectory*`/`ballistics_calculate_zero_angle*`/`ballistics_monte_carlo*` export unchanged — a pure addition to the C ABI requiring no recompile for existing callers.

**solve-json v1:** `atmosphere.pressure_reference` (`"absolute"` or `"qnh"`) mirrors `--pressure-type`; see [docs/SOLVE_JSON_V1.md](docs/SOLVE_JSON_V1.md#pressure_reference-mba-1397). A QNH reduction is recorded in the response's `assumptions` with code `qnh_reduced_to_station_pressure`.

#### Density Altitude as a Direct Input (`--density-altitude`)

Shooter, Nosler, AB Analytics, Ballistic AE, and TRASOL all let a shooter carry conditions in the field as a single **density altitude** value instead of separately dialing in altitude, pressure, and temperature — it is how many US shooters actually record a day's conditions. `--density-altitude` accepts that same single value (feet imperial / meters metric, like `--altitude`).

**This is NOT the direct-density bypass.** `--density-altitude` back-solves an ISA-equivalent station altitude/temperature/pressure and feeds them through the exact same altitude-lapse pipeline every other atmosphere input uses — Mach number, the lapse rate, and downrange-segmented atmosphere (`--wind-segment`'s atmosphere analogue) all behave normally as the shot climbs or falls. Feeding density directly (bypassing that pipeline) would freeze the speed of sound for the whole flight, which is why that path is never used here.

The inverse is the exact algebraic inverse of the published NWS/FAA pressure-altitude model this engine already uses to REPORT density altitude on the PDF dope card:

```
pressure_alt_ft = 145366.45 * (1 - (station_hpa / 1013.25)^0.190284)
isa_temp_f      = 59.0 - pressure_alt_ft / 1000.0 * 3.57
density_alt_ft  = pressure_alt_ft + (120*5/9) * (station_temp_f - isa_temp_f)
```

solved for `station_hpa` given a target density altitude and a station temperature. This model performs no geopotential correction (unlike the ICAO/QNH formula above) — the altitude it returns is treated as GEOMETRIC altitude, the same convention as `--altitude` everywhere else in this engine.

**Precedence** (shares this rule with `--pressure-type` above — read together):
- `--density-altitude` **supersedes `--altitude` and `--pressure`/`--pressure-type` entirely**. A `--location`/CSV `DA`/`DENSITY_ALTITUDE` column is used only when the flag itself is omitted (same CLI-wins-over-CSV rule as every other environmental field).
- The resulting temperature defaults to the **ISA temperature at that density altitude** — algebraically this makes the back-solved station altitude equal the density altitude exactly (see the round-trip test `density_altitude_default_temp_pressure_alt_equals_density_alt` in `src/atmosphere.rs`).
- An **explicit `--temperature` wins outright** over that default, because real consumers need a real air temperature (powder-temperature sensitivity and the moist speed of sound both read the resolved temperature). The station pressure is then re-solved so the SAME density altitude is still reproduced at the real temperature — density is honored either way; only the implied pressure/altitude differ.
- A `--pressure`/`--pressure-type` given ALONGSIDE `--density-altitude` is ignored outright (a note is printed to stderr, or folded into the table-only warning block in WASM, per the same convention `--bc-reference`'s inert-with-`--drag-table` warning uses).

```bash
# A shooter's Kestrel reads 3000 ft density altitude on a warm day; they also know the real
# air temperature (95 F / 35 C) for correct powder-temperature sensitivity.
./ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 \
  --density-altitude 3000 --temperature 95
```

**BC reference standard (MBA-1365) interaction:** independent. `--bc-reference` converts the BC value itself (`inputs.bc_value`); `--density-altitude` resolves the atmosphere (`AtmosphericConditions`/`BallisticInputs.altitude`/`.temperature`/`.pressure`). Neither reads the other's fields, so any combination of the two flags composes exactly as if applied separately.

Available on `trajectory` (CLI flag + `--location`/CSV `DA`/`DENSITY_ALTITUDE` columns) and `profile save --density-altitude` (stored for round-trip parity, alongside the pre-existing `temperature`/`pressure`/`altitude`/`humidity` profile fields — like those, `trajectory --saved-profile` does not read any of them back today; only `--location`/CSV feeds shot-day atmosphere there). **Not yet available** on `zero` or any of the subcommands `--pressure-type` also excludes (`estimate-bc`, `true-velocity`, `plan-truing`, `mpbr`, `come-ups`, `lead`, `wind-card`, `stability`, `range-table`, `compare`) — tracked alongside that same follow-up.

**WASM:** pass `--density-altitude <value>` as a terminal argument to `trajectory`; behavior matches native exactly (same back-solve, same Authoritative-constructor bypass of the legacy default-sentinel heuristic — including when `--density-altitude` is combined with `--pressure-type`).

**FFI:** `FFIAtmosphericConditions` is the same ABI-frozen `repr(C)` struct as above, with no new field. A density-altitude-aware FFI caller calls three standalone conversions once — `ballistics_density_altitude_altitude_m`, `ballistics_density_altitude_temperature_c`, `ballistics_density_altitude_pressure_hpa` (all `(density_altitude_m, explicit_temperature_c)`, where `explicit_temperature_c = NAN` means the ISA-at-DA default) — and writes the three results into `FFIAtmosphericConditions.altitude`/`.temperature`/`.pressure` before calling any existing `ballistics_calculate_trajectory*`/`ballistics_calculate_zero_angle*`/`ballistics_monte_carlo*` export unchanged. A pure addition to the C ABI; no recompile is required for existing callers (only for a caller that wants to opt into density-altitude support).

#### BC5D Correction Tables

BC5D (5-Dimensional BC Correction) tables provide ML-derived, velocity-dependent BC corrections for specific calibers. These tables capture how BC changes throughout the flight envelope based on weight, BC, muzzle velocity, current velocity, and drag model.

**Auto-Download Mode (Requires `--online` feature):**

Tables are automatically downloaded from the server and cached locally:

```bash
# Auto-download tables (downloads on first use)
./ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --bc-table-auto

# Force re-download cached tables
./ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --bc-table-auto --bc-table-refresh

# Use custom server URL
./ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 \
  --bc-table-auto \
  --bc-table-url https://your-server.com/bc5d
```

**Local Directory Mode:**

```bash
# Use predownloaded tables from a local directory
./ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 \
  --bc-table-dir ./bc_tables/
```

**Available Calibers:** .224, .243, .264, .277, .284, .308, .338

**Cache Locations:**
- macOS: `~/Library/Caches/ballistics-engine/bc5d/`
- Linux: `~/.cache/ballistics-engine/bc5d/`
- Windows: `%LOCALAPPDATA%\ballistics-engine\cache\bc5d\`

When a caliber isn't available, you'll see a helpful message:
```
Warning: No BC5D table available for caliber 0.375 (9.5mm)
         Available calibers: .224, .243, .264, .277, .284, .308, .338
         Continuing without BC5D correction table.
```

Both correction-table lookups — the single-file `--bc-table` correction and the
caliber-specific `--bc-table-dir`/BC5D tables above — are **G1/G7 only**: a wider drag
model warns and is treated as G1 for the lookup, same as `true-velocity`/`plan-truing`.

#### BC and Velocity Truing

Adjust BC and velocity based on real-world chrono data and field observations:

```bash
# BC truing - multiply stated BC by adjustment factor (e.g., 0.85 = 85%)
./ballistics trajectory -v 2822 -b 0.270 -m 140 -d 0.264 \
  --bc-adjustment 0.85

# Velocity truing - add offset to base velocity from chronograph data
./ballistics trajectory -v 2822 -b 0.270 -m 140 -d 0.264 \
  --velocity-adjustment 53   # Adds 53 fps to base velocity

# Combined truing
./ballistics trajectory -v 2822 -b 0.270 -m 140 -d 0.264 \
  --bc-adjustment 0.85 \
  --velocity-adjustment 53
# Result: velocity=2875 fps, BC=0.2295
```

#### CSV Profile and Location Support

Load gun profiles and shooting locations from CSV files for batch processing:

**Gun Profile CSV Format** (`gun_profiles.csv`):
```csv
#RIFLE_NAME,VELOCITY,BC,BC_TYPE,BULLET_WEIGHT,CALIBER,ZERO_TEMP,ZERO_ALT,VELOCITY_ADJ,BC_ADJ
AR22,1115,0.138,G1,40,0.22,32,1370,1,1.0
R700_65CM,2822,0.270,G7,140,0.264,57,1806,53,0.85
```

**Location CSV Format** (`locations.csv`):
```csv
LOCATION_NAME,ALTITUDE,PRESSURE,TARGET_TEMP
KF_LR,2506,27.29,32
Home_Range,500,29.92,70
```

A `DA` or `DENSITY_ALTITUDE` column is also recognized (MBA-1366) and, when present, supersedes that row's `ALTITUDE`/`PRESSURE` exactly like the `--density-altitude` flag — see [Density Altitude](#density-altitude-as-a-direct-input---density-altitude) for the full precedence.

Note which flag wins here: an explicit `--density-altitude` overrides the CSV `DA` column (CLI beats CSV, as everywhere else), but `--altitude` and `--pressure` do **not** — density altitude supersedes altitude and pressure whichever source it came from, because it is the stronger specification. A run combining a CSV `DA` row with `--altitude`/`--pressure` prints a note saying so.

**Usage:**
```bash
# Load from profile CSV
./ballistics trajectory \
  --profile gun_profiles.csv \
  --profile-row R700_65CM \
  -m 140 -d 0.264 \
  --max-range 1000

# Load profile + location
./ballistics trajectory \
  --profile gun_profiles.csv --profile-row R700_65CM \
  --location locations.csv --site KF_LR \
  -m 140 -d 0.264

# CLI args override CSV values
./ballistics trajectory \
  --profile gun_profiles.csv --profile-row R700_65CM \
  --velocity 2900 \   # Overrides CSV velocity
  -m 140 -d 0.264
```

### Importing profiles (.a7p)

Import an ArcherBC2 `.a7p` profile into the local profile store:

```bash
# Preview what would be imported (nothing is written)
ballistics profile import my-rifle.a7p --dry-run

# Import under the file's own profile name (sanitized)
ballistics profile import my-rifle.a7p

# Import under a chosen name; fail hard on checksum mismatch
ballistics profile import my-rifle.a7p --name match-338 --strict
```

The import prints a full mapping report: every field it imported (source
value, converted value, destination), every field it could NOT map (for
example powder temperature sensitivity, scope click offsets, and the
device's range-card list), and any warnings (checksum mismatch). Imported
profiles are stored in metric units and can be recalled by name with
`--profile <name>` (on `mpbr`, `come-ups`, `lead`, `wind-card`, `stability`,
and `range-table`) or `--saved-profile <name>` on `trajectory`, which
reserves `--profile` for CSV gun-profile files.

**Multi-BC and CUSTOM drag curves (MBA-1323 Phase 2).** A `.a7p` file with
more than one G1/G7 coefficient row imports ALL rows as a velocity-banded BC
schedule (`bc_segments` in the saved profile JSON), not just the fastest
row — the scalar `bc` field is still set to the fastest row's value for
tools that only understand one BC. A file with `bc_type CUSTOM` (a full
Mach/Cd drag curve) imports as `drag_model: "CUSTOM"` with the curve stored
in `drag_curve`; because no single coefficient applies to a full curve, the
saved profile's `bc` field is set to `0.0` — an intentionally invalid
sentinel, physically inert once the curve is in use, that makes any
consumer which still expects a scalar BC fail loudly instead of silently
solving under an assumed G1 model.

These two fields are consumed automatically by `trajectory --saved-profile`,
`come-ups --profile`, and `lead --profile`: a saved profile's `bc_segments`
feeds the same velocity-keyed BC schedule as `--bc-segment`/`--bc-table-dir`
(only when the run did not already supply one of those), and `drag_curve`
feeds the same custom drag table as `--drag-table` (only when the run did
not already supply `--drag-table`). `mpbr`, `wind-card`, `stability`, and
`range-table` do not yet consume these two fields from a saved profile —
a profile with multi-BC/CUSTOM data still works there, but only via its
scalar `bc`/`drag_model` fallback. `profile show` prints a summary line
(row/point count and range) for whichever of the two is present.

**DSF (drop-scale-factor) table (MBA-1357).** A profile can also carry `dsf_points`, a
Mach-keyed table of drop corrections accumulated by the `dsf` command (see
[DSF Truing](#dsf-drop-scale-factor-truing) below) — `None`/absent for a profile with no
DSF calibration. `profile show` renders it, one line per point:
```
DSF table (2 points, Mach 0.65-0.95):
  Mach 0.65  DSF 1.0820
  Mach 0.95  DSF 1.0310
```
`profile save NAME ... --clear-dsf` removes an existing table; without `--clear-dsf`,
re-saving a profile (e.g. to tweak an unrelated field) carries its DSF table forward
unchanged — `profile save` has no flags of its own to express DSF points, so silently
dropping them on every unrelated edit would be hostile. `trajectory --saved-profile` and
`come-ups --profile` automatically apply a profile's DSF table to the solved drop and
print a table-output-only note; see [DSF Truing](#dsf-drop-scale-factor-truing).

**Turret click graduations (MBA-1355).** `profile save` also accepts
`--elevation-click <SIZE><UNIT>` and `--windage-click <SIZE><UNIT>` (e.g. `0.1mil` or
`0.25moa`), stored as the profile's `elevation_click`/`windage_click` fields and shown by
`profile show`. Both are validated with the same parser `--elevation-click-value` uses,
at save time, so a saved profile can never contain a graduation `resolve_click_values`
would later reject. They are angular graduations, not linear measurements, so unlike most
profile fields they are **not** rescaled when a profile is loaded under the other
`--units` system. See [Turret Adjustment Units](#turret-adjustment-units) above for how
they're resolved against the `--elevation-click-value`/`--windage-click-value` CLI flags.

**Reading a v2 profile with an older tool: one-way forward-incompatibility.**
`bc_segments` and `drag_curve` are additive JSON keys (unknown-key-tolerant,
default-on-absence), so a profile saved by this version always deserializes
cleanly in an older `ballistics` build, an un-updated WASM build, or a
binding that hasn't been regenerated against this schema. For a `CUSTOM`
(full drag-curve) profile that's safe: the older reader still sees
`drag_model: "CUSTOM"` and the inert `bc: 0.0` sentinel, so it refuses to
solve rather than guessing (`bc_value must be finite and greater than
zero`). For a **multi-BC** profile it is *not* safe: the older reader has no
way to know `bc_segments` exists, silently ignores it, and solves using only
the scalar `bc` (the fastest row) for the entire trajectory. This produces a
materially different, unwarned answer whenever the slower bands matter — for
example, a profile whose bands span a wide velocity range showed a ~639 m/s
impact velocity under the scalar-only fallback versus ~411 m/s with the full
schedule applied, for the identical saved profile. There is no equivalent
sentinel guard available for `bc_segments` (a real, plausible BC value is
required there for single-BC tools to keep working at all), so this
particular skew — new profile, old reader, multi-BC case — degrades
silently by design. Do not exchange saved profile JSON across
`ballistics`/WASM/binding versions that straddle this feature (MBA-1323
Phase 2) unless the older side only ever reads single-BC, non-`CUSTOM`
profiles.

**Muzzle velocity and downrange chronographs (MBA-1377).** The `.a7p` format has no field
for the screen distance a device's muzzle velocity was measured at, so a profile's stored
velocity is always taken as a true muzzle velocity, unadjusted. If your chronograph reads
downrange (most do — 10-15 ft, or 25 m, is typical) run `true-velocity` with
`--chrono-velocity`/`--chrono-distance` first to back-solve the corrected muzzle velocity,
then save (or re-import) the profile with that corrected value rather than the raw reading.

The `.a7p` wire format is implemented independently for interoperability;
no third-party code is bundled.

### Sight geometry and zero state

Two related but deliberately distinct input families describe how the sight sits and how
the rifle was zeroed. They must not be conflated:

* **Zero POI offset** (`--zero-poi-up` / `--zero-poi-right`, MBA-1359) — an *angular zero
  state*: "my rifle is deliberately zeroed 0.1 in high / 0.2 in left AT the zero range"
  (sub-click zero error, suppressor on/off shift). Kestrel ZH/ZO semantics.
* **Lateral sight offset** (`--sight-offset`, MBA-1396) — *physical mount geometry*: "my
  optic sits 0.5 in right of the bore axis" (offset mounts, side-mounted clip-ons). The
  windage counterpart of `--sight-height`.

Both may be set together; their effects are additive and distinct — the POI offset shows
up AT the zero range, while a converged mount offset contributes nothing there (see the
convergence rule below).

#### Deliberate Zero POI Offset (`--zero-poi-up`, `--zero-poi-right`)

Values are **inches** under imperial units and **centimeters** under metric, signed:

* `--zero-poi-up 0.1` — the rifle is zeroed to impact 0.1 in **HIGH** at the zero range
  (negative = low).
* `--zero-poi-right -0.2` — zeroed to impact 0.2 in **LEFT** at the zero range
  (positive = right).

The offset is a *linear distance at the zero range*, converted once to the equivalent
angular bias (`offset / zero distance`) and applied to the solved zero elevation and
azimuth after the zero search converges — the standard Kestrel/AB "zero height / zero
offset" model, not a re-zero. It therefore requires a zero solve to act on: pair it with
`--auto-zero` (or a saved profile's zero) on `trajectory`, and it is accepted anywhere a
zero is solved (`mpbr`, `come-ups`, `wind-card`, `range-table`, `compare`).

```bash
# Zeroed deliberately 0.1" high and 0.2" left at 100 yd:
ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 \
  --auto-zero 100 --zero-poi-up 0.1 --zero-poi-right -0.2 --max-range 600
```

The whole solution shifts by the equivalent angle: ~0.1 in high at 100 yd becomes ~0.2 in
at 200 yd, and the drift column carries the constant angular windage bias. Saved profiles
store the offsets as `zero_poi_up_m` / `zero_poi_right_m` (always meters, regardless of
the profile's display units); explicit CLI flags override the stored values.

`.a7p` import: the ArcherBC2 format stores its zeroing state (`zero_x`/`zero_y`) as raw
device **click counts** without the device's click size. Pass
`profile import <file> --zero-click 0.1mil` (or `0.25moa`, etc.) to convert them into the
profile's zero POI offset at the file's zero distance; without `--zero-click` they remain
reported as unmapped, exactly as before.

#### Lateral Sight Offset (`--sight-offset`)

Values are **inches** under imperial units and **millimeters** under metric (the same
units as `--sight-height` — this is its windage counterpart), signed:

* `--sight-offset 0.5` — the sight sits 0.5 in **RIGHT** of the bore axis.
* `--sight-offset -20` (metric) — the sight sits 20 mm **LEFT** of the bore axis.

This is physical mount geometry, present from the instant the bullet leaves the bore: the
bullet starts `offset` LEFT of the sight line (for a sight right of bore). What happens
next depends on whether a zero is solved:

* **With a zero** (`--auto-zero`, or the zero every table command solves): the windage
  zero converges the trajectory onto the sight line — an azimuth term of
  `offset / zero distance` — so drift reads `-offset` at the muzzle, `0` at the zero
  range, and continues right past it. This matches what AB Mobile, JBM, and Shooter
  compute; it is **not** a constant parallel offset.
* **Without a zero** (an explicit `--angle` flight, or `monte-carlo`, which solves no
  zero): only the constant physical displacement applies — drift starts at `-offset` and
  stays uncorrected. The flag's help text states this explicitly.

```bash
# Optic mounted 0.5" right of bore, zeroed at 100 yd:
ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 \
  --auto-zero 100 --sight-offset 0.5 --max-range 600
```

Composes with `--cant` (the cant term swings the bore's sight-HEIGHT offset laterally;
this flag adds the fixed mount displacement on top) and with the zero POI offsets above
(angular zero state vs physical geometry — both azimuth terms share one convergence
point in the zero solve). Saved profiles store it as `sight_offset_lateral_m` (always
meters); the explicit CLI flag overrides the stored value.

#### Named Zero Sets and Per-Load Offsets (`profile zero-set`, `--zero-set`)

Saved profiles can carry **multiple named zero conditions** (MBA-1360 — the Lapua
Sight-In POI / ATrag zero-zone / Strelok multi-zero feature class): alternate zero
distances and constant per-load dial corrections, stored under the profile's
`zero_sets` key and selected per run.

```bash
# Store a 200 yd suppressed zero that also needs 0.3 mil LESS elevation and
# 0.1 mil MORE right windage than the master zero:
ballistics profile zero-set add mycreedmoor --name suppressed \
  --zero-distance 200 --poi-up -0.3 --poi-right 0.1 --notes "suppressed load"

ballistics profile zero-set list mycreedmoor
ballistics profile zero-set remove mycreedmoor --name suppressed

# Use it:
ballistics come-ups --profile mycreedmoor --zero-set suppressed \
  --zero-distance 200 --start 200 --end 800 --step 100
```

Semantics:

* `--poi-up` / `--poi-right` are **dial corrections in MILs**, ADDED to elevation /
  windage adjustments when the set is selected: positive = dial UP / RIGHT more. A
  load that impacts 0.25 mil HIGH relative to the master zero therefore stores
  `--poi-up -0.25`. This is a *constant angular* correction (unlike the
  linear-at-zero-range `--zero-poi-up` above, which biases the solved zero itself);
  the two compose. In non-MIL display units the correction is rescaled through the
  same locked factor table as every other dial conversion, and with a tracking CF
  (MBA-1358) the order is fixed: `dial = (true need + set correction) / CF`.
* A set's `zero_distance` (display units, converted with the profile) feeds the
  auto-zero exactly as the profile's own zero would; **explicit CLI zero flags always
  win** over the set.
* `--zero-set NAME` exists on `trajectory`, `come-ups`, `wind-card`, `range-table`,
  `dsf`, and `plan-truing`. The dial corrections reach *total-correction* dial outputs
  (come-ups, wind-card drift, range-table Drop/Wind, the PDF dope card's Drop/Wind);
  component holds (mover Ring, lead) deliberately stay uncorrected — they compose on
  top of a wind dial that already carries the correction. On `dsf`/`plan-truing` only
  the set's zero distance participates (a DSF observation is a physical impact, not a
  dial value). `compare`, `lead`, `mpbr`, and `stability` take no `--zero-set` (no
  single profile, no total dial output, or no fixed zero to select).
* An unknown name is a **hard error listing the available sets**; nothing is ever
  silently ignored.

Data sources that feed sets automatically:

* **Profile CSVs** (`trajectory --profile file.csv --profile-row R1`): the
  long-allowlisted `V_OFFSET_MIL` / `H_OFFSET_MIL` columns now form an ephemeral zero
  set named after the row (dial-correction convention, as above) — select it with
  `--zero-set R1`. Without `--zero-set`, the columns remain inert exactly as before.
* **`.a7p` import with `--zero-click`**: in addition to the MBA-1359 engine-field
  mapping described above (unchanged), the file's `zero_x`/`zero_y` click state is
  recorded as a set named `a7p-zero` in dial-correction convention (the negated
  angular POI offset). Select it only on a profile whose `zero_poi_up_m`/
  `zero_poi_right_m` fields have been cleared — the engine fields already shift the
  whole solution, so applying both double-counts the same click state.

Compatibility: profiles carrying `zero_sets` load and solve identically on older
binaries **as long as no alternate set is requested** (the master zero is untouched;
`--zero-set` is an unknown flag there and fails loudly). Re-saving a profile with an
older binary drops the stored sets — the same documented one-way skew as
`bc_segments`.

### Canted Shooting

Model a rifle that is zeroed level but *fired* with the scope/receiver rotated about the
line of sight — the classic "canted rifle" error. Available via `--cant <DEGREES>`
(alias `--cant-angle`) on the `trajectory` and `monte-carlo` subcommands. Default `0` =
level (bit-identical to a solve without the flag).

**Not available on `zero`.** Zeroing always solves the un-canted trajectory — cant is
applied only at fire time. This models "zero the rifle level, then shoot it canted," not
"the rifle was canted while zeroing," which would (mostly) cancel out.

**Sign convention:** positive degrees = clockwise cant as seen from behind the rifle (the
top of the scope tips to the right). For a rifle with an upward zero elevation
correction — the normal case — that rotates the correction partly into windage, so point
of impact moves **right and low** relative to the un-canted zero.

**Error model:** cant rotates the sight-frame aim offsets (elevation and windage) about
the line of sight, and swings the bore's sight-height offset laterally with it. For small
cant angles the combined lateral error at a given range is approximately

```
lateral_error ≈ (D − sight_height) · sin(cant)
```

where `D` is the height the zero's elevation correction adds at that range (i.e. how much
higher the zeroed, un-canted trajectory sits than a flat 0°-elevation shot would). Since
`D` grows with range (and with how much the load has dropped), **cant error grows with
range** — a small cant is barely noticeable at zero range and increasingly costly beyond
it. This is validated in `tests/canted_fire.rs` to within 5% (windage) / 10% (elevation)
of the analytic prediction at 300 m and 600 m.

**Which elevation rotates — matching the field rule of thumb.** `--cant` rotates the
elevation that is *in the gun* (`muzzle_angle`, i.e. whatever the zero/dial put there).
Two scenarios that are easy to conflate:

- **Zeroed at a near distance, never re-dialed** (e.g. `--auto-zero 100`, then read
  windage at long range): only the small 100 yd zero elevation rotates, so the lateral
  error is nearly **constant in mils** across range (≈ `zero_elevation · sin(cant)` —
  roughly 0.1 mil for 5° of cant on a typical rifle zero).
- **Dialed (or held over) for the engagement range, then canted** — the realistic
  long-range case: model it by zeroing at the engagement distance
  (`--auto-zero 1000 --cant 5`). The *entire* come-up rotates, reproducing the classic
  field rule `lateral ≈ come_up · sin(cant)`. Example: a .224 77 gr at 2650 fps (G7 BC 0.372) dialed
  from 100 yd to 1000 yd carries ≈ 9.7 mil of total launch elevation; at 5° of cant the
  engine puts the shot 0.86 mil right — `9.7 · sin(5°) ≈ 0.85 mil`. The two rules agree;
  they just answer different questions.

**Worked example** (build first with `cargo build`; `--sample-interval` is always meters
regardless of `--units`, so `91.44` below is exactly 100 yd):

```bash
# Level (no cant)
./target/debug/ballistics trajectory -v 2700 -m 168 -d 0.308 --bc 0.5 \
  --auto-zero 100 --max-range 600 \
  --sample-trajectory --sample-interval 91.44 -o csv --full

# Same load, 10 degrees of clockwise cant
./target/debug/ballistics trajectory -v 2700 -m 168 -d 0.308 --bc 0.5 \
  --auto-zero 100 --max-range 600 --cant 10 \
  --sample-trajectory --sample-interval 91.44 -o csv --full
```

Level output — `drift_in` is exactly zero at every range:

```
distance_yd,drop_in,drift_in,velocity_fps,energy_ft-lb,time_s
0.00,2.00,0.00,2700.00,2718.96,0.0000
100.00,0.00,0.00,2519.55,2367.68,0.1150
200.00,3.48,0.00,2346.19,2053.06,0.2384
300.00,13.29,0.00,2179.61,1771.88,0.3711
400.00,30.43,0.00,2019.83,1521.62,0.5141
500.00,56.11,0.00,1867.27,1300.44,0.6686
```

10-degree canted output — the right-and-low effect grows with range, so `drift_in` and
`drop_in` climb steadily above the level case with distance. Right at the muzzle it's the
opposite: the bore itself swings toward the cant pivot before the zero's elevation
correction has any range to leak into windage, so `drop_in` briefly runs *lower* than
level (1.97 vs 2.00 at 0 yd) before flipping higher by 100 yd:

```
distance_yd,drop_in,drift_in,velocity_fps,energy_ft-lb,time_s
0.00,1.97,-0.35,2700.00,2718.96,0.0000
100.00,0.04,0.43,2519.55,2367.68,0.1150
200.00,3.59,1.22,2346.19,2053.06,0.2384
300.00,13.46,2.00,2179.61,1771.88,0.3711
400.00,30.67,2.78,2019.83,1521.62,0.5141
500.00,56.42,3.56,1867.27,1300.44,0.6686
```

The small offset already present at the muzzle (`drift_in = -0.35` at 0 yd) is the bore
itself swinging laterally below the canted sight; windage then climbs through positive
(rightward) values as the zero's elevation correction leaks into windage with range.

**Monte Carlo caveat:** cant is a *systematic* aim bias, not a dispersion source. Because
`monte-carlo` reports statistics as deviations about its own (canted) mean, `--cant` shifts
the whole cloud together and has almost no effect on the reported spread — expect the
dispersion numbers to look essentially the same as a level run. Use `trajectory --cant` to
see the point-of-impact shift itself.

#### Drops Reference: LOS vs Target Plane (`--drops-reference`) — MBA-1403

**The convention.** The engine solves inclined fire in the shot frame (gravity is rotated
by `--shooting-angle`), so a sampled drop value is natively measured **perpendicular to
the line of sight** — that is the `los` reference, the default, and the only convention
every output used before this flag existed. JBM's "target plane" checkbox instead reports
drop **vertically in the plane of the target**, which matters when holding over a steeply
inclined shot against a vertical target face: the two differ by exactly
`1 / cos(shooting angle)` (about 15% at 30°). `--drops-reference target` applies that
transform to the sampled drop values — and, because the mil/MOA conversions derive from
the same drop, to those as well. Wind drift is a lateral quantity and is never touched.

```bash
# LOS-referenced (default; byte-identical with or without the flag)
ballistics trajectory -v 2700 -m 175 -d 0.308 --bc 0.243 --drag-model g7 \
  --shooting-angle 30 --auto-zero 100 --max-range 600 \
  --sample-trajectory --sample-interval 100

# Target-plane referenced: sampled drops are divided by cos(30) and the drop
# column is relabeled "Drop (target)" (CSV: drop_target_in)
ballistics trajectory -v 2700 -m 175 -d 0.308 --bc 0.243 --drag-model g7 \
  --shooting-angle 30 --auto-zero 100 --max-range 600 \
  --sample-trajectory --sample-interval 100 --drops-reference target
```

Details worth knowing:

- The table and CSV drop columns are labeled with the active reference (`Drop (target)`
  / `drop_target_in`), so a card never leaves ambiguity about its convention.
- This is an **output-mode toggle only**: the solved trajectory, drift, velocity, energy,
  time, zeroing, and `--auto-zero` semantics are all unchanged.
- At `--shooting-angle 0` the two references coincide exactly.
- `target` mode is rejected at shooting angles of 90° or beyond, where the transform is
  undefined.
- The solve-json v1 surface accepts the same toggle as `$.shot.drops_reference`
  (`"los"`/`"target"`); the WASM browser terminal accepts `--drops-reference` with the
  same semantics.

#### Equivalent Horizontal Range (BDC shoot-to) — MBA-1395

When a look angle is set (`--shooting-angle` nonzero) **and** a zero was solved
(`--auto-zero` or a profile zero), the trajectory table prints one extra summary line:

```
Equivalent horizontal range: 469 yd (shoot-to for BDC)
```

This is the range SIG BDX (AMR), Leica (EHR), and Gunwerks BR2 report: the **flat-fire
range whose angular elevation correction — against the same zero — matches the inclined
solution's** at the solved range. A shooter with a fixed BDC turret or reticle
calibrated on flat ground dials or holds for the printed range instead of the true
slant range. The match is **angular**, computed by an inverse lookup over one extra
flat solve (McDonald's Sierra inclined-fire treatment / Litz), *not* the rifleman's-rule
`range × cos(angle)` approximation — matching linear drop that way reproduces
rifleman's-rule-class error at long range.

The line is deliberately omitted where the inverse is ill-defined: flat shots (output
byte-identical to before this feature), targets at or inside the zero range, and
non-positive corrections (bullet at/above the line of sight). solve-json v1 carries the
same value as `summary.equivalent_horizontal_range_m` (absent under the same
conditions); the WASM terminal prints the same summary line after the trajectory table.

### Moving-Target Lead

Calculate the hold needed to hit a target moving at a constant ground speed across (or
along) the line of sight, swept over a range of distances. Available via the `lead`
subcommand.

```bash
./ballistics lead -v 2700 -m 168 -d 0.308 -b 0.5 --target-speed 3
```

Besides the usual load/atmosphere/wind arguments shared with `trajectory` (`-v -b -m -d
--drag-model --sight-height --temperature --pressure --humidity --altitude --wind-speed
--wind-direction`, plus `--profile`), `lead` also accepts `trajectory`'s powder-temperature
flags — `--use-powder-sensitivity`, `--powder-temp-sensitivity`, `--powder-temp`, and
`--powder-temp-curve` — plumbed identically (MBA-1325: both commands build the same
`BallisticInputs` and resolve the correction in the same place,
`cli_api::TrajectorySolver::new`), so a `lead` run can reproduce a powder-corrected muzzle
velocity without a separate `trajectory` call first. See the [Trajectory
Command](#trajectory-command) parameters table for their units/defaults; omitting all four is
unchanged from before (`-v` used verbatim).

`lead` then adds its own moving-target arguments:

- **`--target-speed <SPEED>`** (required) — target ground speed, mph under imperial units,
  m/s under metric.
- **`--target-angle <DEGREES>`** — direction of target *travel* relative to the line of
  sight (default `90`):
  - `0` = directly away (outbound)
  - `90` = crossing left-to-right (full broadside)
  - `180` = directly toward (inbound)
  - `270` = crossing right-to-left
- **`--target-length <LENGTH>`** — target body length (inches imperial, mm metric). When
  given, an extra `Bodies` column reports lead as a multiple of the target's length
  (`lead ÷ target_length`) — a common visual hold reference ("hold one body-length ahead").
- **`--start` / `--end` / `--step`** — range sweep in yards (imperial) or meters (metric),
  like the other sweep tables; defaults `100`/`600`/`100`.
- **`--adjustment-unit <mil|moa|smoa|iphy|clicks>`** — angular unit for the `Lead` column
  (default `mil`); the column is windage-type (a horizontal hold), so `clicks` resolves
  through `--windage-click-value` (MBA-1410) — there is no `--elevation-click-value` flag
  here. See [Turret Adjustment Units](#turret-adjustment-units) for the full unit list.
- **`--windage-click-value <SIZE><UNIT>`** — turret click graduation for
  `--adjustment-unit clicks`, e.g. `0.25moa` or `0.1mil` (falls back to the saved
  profile's `windage_click`, then its `elevation_click`, when omitted).
- **`-o, --output <table|json|csv|pdf>`** — output format (`pdf` renders the same as `table`
  on this subcommand).

**Locked conventions:**

- **Positive lead = hold in the target's direction of travel.** For a 90° (left-to-right)
  crosser, positive lead is a hold to the right.
- **`lead_mil = (lead/range)·1000`; `lead_moa = (lead/range)·3438`** — the same dial
  convention used by every other hold table in this CLI (MBA-724): MOA is exactly 3.438×
  MIL, not the geometrically exact 3437.7467×.
- **Lead is pure target motion, additive to your wind-corrected hold.** Time of flight
  comes from the engine's wind-aware trajectory solve, but wind deflection itself is *not*
  folded into the lead number — it stays in the separate wind column of your dope; add the
  two holds together.
- **Time of flight is wind-aware** — the underlying solve accounts for wind drag effects on
  TOF even though the lead figure itself reports only target-motion offset.
- **Non-perpendicular motion (angle ≠ 90°/270°) shifts the intercept range.** An outbound
  or inbound target has moved farther or closer by the time the bullet arrives, so
  `calculate_lead` fixed-point iterates `R = R₀ + v_radial·TOF(R)` until the correction is
  below 0.1 m, and reports TOF/lead at that corrected range — the table's `Intercept`
  column shows the range actually used, which differs from the requested `Range` for any
  non-perpendicular angle. A target closing faster than the geometry allows (or one whose
  corrected range runs past the solved trajectory span) produces a typed error printed
  inline in place of that row's data instead of a bogus number.

**Worked example** (build first with `cargo build`):

```bash
./target/debug/ballistics lead -v 2700 -m 168 -d 0.308 -b 0.5 --target-speed 3
```

```
Moving-Target Lead Table (target speed: 3.0 mph, angle: 90°, MIL)
Positive lead = hold in the direction of target travel.

┌──────────┬──────────┬──────────┬──────────┬──────────┐
│Range (yd)│TOF (s)   │Lead ( yd)│Lead (MIL)│Intercept │
├──────────┼──────────┼──────────┼──────────┼──────────┤
│      100 │    0.115 │     0.17 │    1.687 │    100.0 │
│      200 │    0.238 │     0.35 │    1.748 │    200.0 │
│      300 │    0.371 │     0.54 │    1.814 │    300.0 │
│      400 │    0.514 │     0.75 │    1.885 │    400.0 │
│      500 │    0.669 │     0.98 │    1.961 │    500.0 │
│      600 │    0.836 │     1.23 │    2.043 │    600.0 │
└──────────┴──────────┴──────────┴──────────┴──────────┘
```

For a pure 90° crosser the `Intercept` column always equals `Range` (no radial motion to
correct for). The MIL figure climbs slowly with range because the bullet slows down, so
each added yard of range costs more time of flight — and more lateral target travel per
yard — than the last.

`-o json` (trimmed to the first two rows):

```json
{
  "adjustment_unit": "MIL",
  "distance_unit": "yd",
  "rows": [
    {
      "intercept_range": 100.0,
      "iterations": 0,
      "lead": 0.16871387472963065,
      "lead_mil": 1.6871387472963066,
      "lead_moa": 5.800383013204702,
      "range": 100.0,
      "tof_s": 0.11503218731565724
    },
    {
      "intercept_range": 200.0,
      "iterations": 0,
      "lead": 0.34969310252646174,
      "lead_mil": 1.7484655126323088,
      "lead_moa": 6.011224432429878,
      "range": 200.0,
      "tof_s": 0.2384271153589512
    }
  ],
  "target_angle": 90.0,
  "target_speed": 3.0,
  "target_speed_unit": "mph",
  "units": "imperial"
}
```

`lead_moa / lead_mil` is 3.438 on every row, confirming the dial convention above.
`iterations: 0` on a perpendicular crosser — there's no radial motion to fixed-point
iterate on.

Non-perpendicular motion with body-length holds:

```bash
./target/debug/ballistics lead -v 2700 -m 168 -d 0.308 -b 0.5 --target-speed 15 \
  --target-angle 45 --target-length 40 --start 200 --end 400 --step 100
```

```
Moving-Target Lead Table (target speed: 15.0 mph, angle: 45°, MIL)
Positive lead = hold in the direction of target travel.

┌──────────┬──────────┬──────────┬──────────┬──────────┬──────────┐
│Range (yd)│TOF (s)   │Lead ( yd)│Lead (MIL)│Intercept │Bodies    │
├──────────┼──────────┼──────────┼──────────┼──────────┼──────────┤
│      200 │    0.240 │     1.24 │    6.184 │    201.2 │     1.12 │
│      300 │    0.374 │     1.94 │    6.419 │    301.9 │     1.74 │
│      400 │    0.518 │     2.69 │    6.671 │    402.7 │     2.42 │
└──────────┴──────────┴──────────┴──────────┴──────────┴──────────┘
```

At 45° the target has an outbound (receding) component, so `Intercept` runs a yard or two
past the requested `Range` — the iteration converges on the slightly longer range the
bullet actually has to cover by the time it arrives. `Bodies` reports lead as a multiple of
the 40" target length (e.g. "hold 1.1 body-lengths ahead" at 200 yd).

**Library API:** for programmatic use, `ballistics_engine::calculate_lead(inputs, wind,
atmo, target_speed_mps, target_angle_deg, range_m) -> Result<LeadSolution, LeadError>` runs
the same wind-aware solve and intercept-range iteration directly, without going through the
CLI. `LeadSolution` carries `time_of_flight_s`, `lead_m`, `lead_mil`, `lead_moa`,
`corrected_range_m`, and `iterations`; `LeadError` is a typed enum covering invalid input,
an over-closing (`TargetOvertakesShooter`) target, iteration `Convergence` failure, a
corrected range that runs `BeyondSolvedSpan`, and an underlying trajectory-solve failure
(`Solver`).

### Mover Ring (`--target-speed`)

A field-tested alternative to `lead` for engaging movers (MBA-1325): instead of computing a
directional hold, `trajectory --target-speed <SPEED>` derives a **ring radius** —
`target_speed × time-of-flight-to-that-point` — around your hold point at every printed/
exported trajectory point. Watch the target through your optic; fire the instant it crosses
into the ring. Because ring size only needs time of flight (which the trajectory solve already
produced) and target speed, it falls out of an already-solved trajectory as pure
post-processing — no second command, no re-entered ballistic data, and no assumed crossing
angle (unlike `lead --target-angle`).

```bash
./ballistics trajectory -v 2700 -m 168 -d 0.308 -b 0.5 --target-speed 3 --full
```

- **`--target-speed <SPEED>`** — mph under imperial units, m/s under metric (same convention
  and same `0`–`300` accepted range as `lead --target-speed`; out-of-range values are
  rejected, not clamped). `0` (the default) leaves every output format byte-identical to a
  run without the flag. This is the same flag that drives the PDF dope card's `Lead` column
  (see [PDF Dope Card Format](#pdf-dope-card-format)) — setting it turns on both at once.
- **`--adjustment-unit <mil|moa|smoa|iphy|clicks>`** — angular unit for the ring
  **table** column only (default `mil`; the flag trajectory already exposes for the PDF
  dope card). With `moa` the column reads `Ring(moa)` with `ring_moa = ring_mil × 3.438`
  — the CLI's locked printed-table dial convention (MBA-724, deliberately not the
  exact-angle 3437.7467/1000), so Ring keeps the same MIL/MOA ratio as every other hold
  column; `smoa`/`iphy` share that ratio too (see [Turret Adjustment
  Units](#turret-adjustment-units)). With `clicks` the column reads `Ring(clicks)` and
  rounds to whole turret clicks against the resolved **elevation** graduation (the Ring
  isn't cleanly an elevation- or windage-axis hold, so it reuses the same graduation as
  the dope card's Drop column) — requires `--elevation-click-value` or a saved profile's
  `elevation_click`. CSV keeps `ring_mil` and JSON keeps `mover_ring_m`/`mover_ring_mil`
  regardless — their names are the unit contract.

**Table** (`--full -o table`) gains a `Ring(mil)` column (`Ring(moa)`/`Ring(smoa)`/
`Ring(iphy)`/`Ring(clicks)` under the matching `--adjustment-unit`). The muzzle point
prints `-` (no flight time has elapsed, so the ring has no defined angular size there
yet):

```
Trajectory Points:
┌──────────┬──────────┬──────────┬──────────┬──────────┬──────────┐
│ Time (s) │  X (yd)  │  Y (yd)  │ Vel(fps) │Energy(ft-lb)│ Ring(mil)│
├──────────┼──────────┼──────────┼──────────┼──────────┼──────────┤
│    0.000 │     0.00 │     1.67 │  2700.00 │  2718.96 │        - │
│    0.152 │   130.61 │     1.55 │  2465.74 │  2267.62 │     1.71 │
│    0.302 │   248.75 │     1.21 │  2264.16 │  1912.00 │     1.78 │
│    0.452 │   357.43 │     0.66 │  2086.99 │  1624.49 │     1.85 │
│    0.582 │   444.86 │     0.03 │  1950.48 │  1418.91 │     1.92 │
└──────────┴──────────┴──────────┴──────────┴──────────┴──────────┘
```

**JSON** (`--full -o json`) adds two per-point fields, present only when `--target-speed > 0`:

```json
{
  "time": 0.17187200000000002,
  "x": 0.0,
  "y": 1.513531761659063,
  "z": 146.9565529362786,
  "velocity": 2437.286215359467,
  "energy": 2215.5854770501155,
  "mover_ring_m": 0.23050097664000005,
  "mover_ring_mil": 1.715329655579473
}
```

- **`mover_ring_m`** — linear ring radius in **meters**, always, regardless of `--units` — the
  unit is in the field name so it can't be silently misread as something else (the legacy
  `trajectory[].x`/`y`/`z` fields don't have that luxury, being bare numbers; see the `legend`
  block documented under [JSON Format](#json-format), MBA-1315, for how those are labeled
  instead). Present at every point once the flag is on, including the muzzle (`0.0`, since
  `target_speed × 0 = 0`).
- **`mover_ring_mil`** — `mover_ring_m / downrange_m × 1000`. **Omitted** at the muzzle
  (`downrange = 0`, the ratio is undefined), not emitted as `0` or `null`.

**CSV** (`--full -o csv`) gains a trailing `ring_mil` column — the header carries the unit; the
muzzle row's field is empty rather than `0`:

```
time_s,x_yd,y_yd,z_yd,velocity_fps,energy_ft-lb,ring_mil
0.0000,0.00,1.67,0.00,2700.00,2718.96,
0.0010,0.00,1.67,0.90,2698.34,2715.63,1.630
0.0028,0.00,1.67,2.52,2695.37,2709.63,1.631
```

**Reading the ring correctly:**

- **It's a worst-case bound, not an exact intercept solution — it assumes the hold point sits
  on the mover's track.** If the mover is heading straight for your hold point, firing the
  instant it enters the ring puts the bullet there exactly when the mover arrives (both cover
  their respective distance — bullet's flight, mover's ring radius — in the same time). If your
  hold point is only *near* the mover's actual line of travel rather than on it, entering the
  ring no longer means "arrives at the hold point in one ToF"; treat the ring as "no later than
  now," not a promise of an exact hit. Use `lead --target-angle` instead when you know the
  crossing angle and want the angle-aware exact hold — the ring trades that precision for not
  needing to know the angle at all.
- **The mil value is only near-constant across a short stage — it grows slowly with range, not
  a single number for the whole sheet.** In the table above it climbs from 1.71 to 1.92 mil
  over 130–445 yd: the bullet decelerates, so later range gates cost more time of flight per
  yard, and the mover's linear ring radius grows faster than downrange distance does. Read the
  ring size for the range you're actually engaging at.
- **Doesn't fold in wind.** Like `lead`'s hold, the ring is pure target-motion bookkeeping —
  it stays separate from the wind column/dial on the same dope.

### Terminal Chart (`--plot`)

Render four stacked inline terminal charts after the normal `trajectory` output (MBA-1320,
extended by MBA-1394): drop vs. range, lateral drift vs. range, velocity vs. range, then
energy vs. range, in that order. Pure Rust, zero new dependencies — no terminal-graphics
crate, no terminal-size detection, no ANSI colors.

```bash
# Bare --plot: default Unicode braille-dot renderer
./ballistics trajectory -v 2700 -m 168 -d 0.308 -b 0.5 --wind-speed 10 --wind-direction 90 --plot

# --plot ascii: '*'-per-cell fallback for terminals/fonts without braille glyph coverage
./ballistics trajectory -v 2700 -m 168 -d 0.308 -b 0.5 --wind-speed 10 --wind-direction 90 --plot ascii
```

- **`--plot`** (bare) — the default renderer: each terminal character cell packs a 2 (wide)
  x 4 (tall) grid of dots into one Unicode Braille Patterns glyph (`U+2800`–`U+28FF`),
  giving roughly 4x the vertical and 2x the horizontal resolution of one character cell.
- **`--plot ascii`** — a `'*'`-per-cell fallback for terminals/fonts without full braille
  glyph coverage. Same dot-addressed layout and axis scaling; only the dot canvas itself
  changes between the two styles — the frame (`┌─┐│└┘`) and axis-range text
  stay ordinary Unicode box-drawing either way, since those glyphs have near-universal font
  support (unlike the braille block, which is why it needs a fallback at all).
- Only affects `-o table` (the default output format). `-o json`/`-o csv`/`-o pdf` are
  completely unaffected, so scripts parsing those formats never see chart text, and
  omitting `--plot` leaves every output format byte-identical to a pre-MBA-1320 run.
- Fixed 72x12-cell canvas per chart. There's no terminal-size detection — that would need a
  dependency this feature deliberately doesn't take on.
- Deliberately monochrome: no ANSI color/SGR codes anywhere in the renderer. That sidesteps
  `NO_COLOR` (<https://no-color.org/>) entirely — there's nothing to suppress — and keeps
  output byte-identical whether the terminal honors color, redirects to a file, or is a
  dumb pipe.
- All four charts plot the SAME per-point data the `--full` "Trajectory Points:" table
  prints (`result.points`, the raw un-decimated integration output — `--plot` works
  without `--full` too, it just doesn't print the table itself). Drop is the table's `Y`
  column (vertical position), lateral drift is the table's downrange-paired `Z` column
  (not printed by the table by default) — both in the same range unit (yd/m) the rest of
  the table uses, **not** inches. This is a different, deliberate convention from
  `--sample-trajectory`'s sight-line-relative `drop_m`/`wind_drift_m` (see
  [Trajectory Sampling for Analysis](#trajectory-sampling-for-analysis)); don't conflate
  the two.
- Velocity and energy panels (MBA-1394) plot the same points' velocity magnitude and
  kinetic energy, in the same units the summary box / `--full` CSV / `-o json` trajectory
  columns already use for this command: fps/ft-lb under imperial units, m/s/J under
  metric.

Example (`--plot ascii`, 10 mph 90° crosswind):

```
Drop vs Range:
┌ drop (yd) — y:[0.00, 1.67] ────────────────────────────────────────────┐
│** ** *** ** *** *** *                                                  │
│                      *** *** *                                         │
│                               *** ***                                  │
│                                      ** **                             │
│                                           **** *                       │
│                                                 ***                    │
│                                                    *** *               │
│                                                         ***            │
│                                                            ***         │
│                                                               ***      │
│                                                                  ** *  │
│                                                                      **│
└ x:[0.00, 448.79] ──────────────────────────────────────────────────────┘

Lateral Drift vs Range:
┌ drift (yd) — y:[-0.44, 0.00] ──────────────────────────────────────────┐
│** ** *** ** *** ***                                                    │
│                     **** *** *                                         │
│                               *** **                                   │
│                                     *** **                             │
│                                           ****                         │
│                                                ****                    │
│                                                    ***                 │
│                                                        ****            │
│                                                            ***         │
│                                                               ***      │
│                                                                  ** *  │
│                                                                      **│
└ x:[0.00, 448.79] ──────────────────────────────────────────────────────┘

Velocity vs Range:
┌ velocity (fps) — y:[1944.46, 2700.00] ─────────────────────────────────┐
│** **                                                                   │
│      *** *                                                             │
│           * ***                                                        │
│                 *** **                                                 │
│                       ** **                                            │
│                            * ****                                      │
│                                   *****                                │
│                                         *****                          │
│                                              * *****                   │
│                                                     ** ***             │
│                                                           *******      │
│                                                                  ** ***│
└ x:[0.00, 448.79] ──────────────────────────────────────────────────────┘

Energy vs Range:
┌ energy (ft-lb) — y:[1410.18, 2718.96] ─────────────────────────────────┐
│** *                                                                    │
│    * ***                                                               │
│          ** **                                                         │
│               * ***                                                    │
│                     ****                                               │
│                          *** **                                        │
│                                ** ***                                  │
│                                      ** ***                            │
│                                            *** ***                     │
│                                                   **** **              │
│                                                          *******       │
│                                                                 *** ***│
└ x:[0.00, 448.79] ──────────────────────────────────────────────────────┘
```

### Wind Card

Generate a wind-drift dope card: deflection at a sweep of ranges, one column per wind
speed. Available via the `wind-card` subcommand.

```bash
./ballistics wind-card -v 2700 -m 168 -d 0.308 -b 0.5 --zero-distance 100
```

Besides the usual load/atmosphere arguments shared with `trajectory` (`-v -b -m -d
--drag-model --sight-height --temperature --pressure --humidity --altitude`, plus
`--profile`), `wind-card` adds:

- **`--zero-distance <DIST>`** (required) — zero distance, yards imperial / meters metric.
- **`--wind-speeds <CSV>`** — comma-separated wind speeds, one column per value (default
  `5,10,15,20`, mph imperial / m/s metric).
- **`--wind-angle <DEG>`** — a single wind angle in degrees, wind-FROM convention (same as
  `--wind-direction` on `trajectory`): **`0` = headwind, `90` = from the right (full
  value), `180` = tailwind, `270` = from the left.** Mutually exclusive with
  `--wind-angles`.
- **`--wind-angles <CSV>`** — comma-separated wind angles; emits one *complete* card per
  angle (e.g. `--wind-angles 30,60,90`). Mutually exclusive with `--wind-angle`.
- **`--start` / `--end` / `--step`** — range sweep, like the other sweep tables; defaults
  `100`/`1000`/`100`.
- **`--adjustment-unit <mil|moa|smoa|iphy|clicks>`** — angular unit for the drift columns
  (default `mil`); the card is windage-type (crosswind drift), so `clicks` resolves
  through `--windage-click-value` (MBA-1410) — there is no `--elevation-click-value` flag
  here. See [Turret Adjustment Units](#turret-adjustment-units).
- **`--windage-click-value <SIZE><UNIT>`** — turret click graduation for
  `--adjustment-unit clicks`, e.g. `0.25moa` or `0.1mil` (falls back to the saved
  profile's `windage_click`, then its `elevation_click`, when omitted).
- **`-o, --output <table|json|csv|pdf>`** — output format (`pdf` renders the same as
  `table` on this subcommand).

**Default (no flags) is the classic full-value 90° card, unchanged.** With neither
`--wind-angle` nor `--wind-angles`, the card is computed at a fixed 90° (full-value
crosswind from the right) and the output is byte-identical to the pre-oblique-angle CLI:
table title says "full-value crosswind" with no angle suffix, and JSON carries
`"crosswind": "full-value (90°)"` instead of a `wind_angle` key.

**Sign convention:** drift values are signed the same way as `drift_in`/`x_yd` elsewhere in
this CLI — positive = rightward (dial right), negative = leftward (dial left). Because the
default card is wind *from the right*, its values run negative; `--wind-angle 270` (wind
from the left) produces the exact mirror image — equal magnitude, opposite sign. `0°`
(headwind) and `180°` (tailwind) both drift to `0.0` (no crosswind component).

**Each angle's cells are a real solve, not a scaled copy.** `wind-card` doesn't take the
90° column and multiply by `sin`/`cos` of the requested angle — every `(speed, angle)`
pair runs its own full trajectory solve. The result tracks `sin(angle)` of the full-value
column closely (within ~1% for the pairs checked) but is not an exact match, because drag
and time-of-flight differ slightly between a purely crosswind solve and an oblique one.

**CSV/JSON self-identify when non-default.** A card produced by `--wind-angle` or
`--wind-angles` prints a `# wind_angle=<DEG>` comment line above its CSV header (one per
card, blank-line separated for multi-angle output) and carries a `"wind_angle"` key in its
JSON object. The legacy no-flag card has neither — it keeps the original `"crosswind"` key
and no CSV comment line, so existing parsers of the default card don't need to change.
Requesting more than one angle (`--wind-angles`) changes the JSON shape from a single
object to an array of one object per angle, in the order given.

**Worked example** (build first with `cargo build`):

```bash
# Default: full-value 90° crosswind card
./target/debug/ballistics wind-card -v 2700 -m 168 -d 0.308 -b 0.5 --zero-distance 100 --end 300

# Same load, a single oblique 45° card
./target/debug/ballistics wind-card -v 2700 -m 168 -d 0.308 -b 0.5 --zero-distance 100 --end 300 --wind-angle 45
```

Default card:

```
Wind Card (zero: 100 yd, MIL, full-value crosswind)
┌──────────┬──────────┬──────────┬──────────┬──────────┐
│Range (yd)│       5 mph │      10 mph │      15 mph │      20 mph │
├──────────┼──────────┼──────────┼──────────┼──────────┤
│      100 │     -0.1 │     -0.2 │     -0.3 │     -0.4 │
│      200 │     -0.2 │     -0.4 │     -0.6 │     -0.8 │
│      300 │     -0.3 │     -0.6 │     -0.9 │     -1.2 │
└──────────┴──────────┴──────────┴──────────┴──────────┘
```

45° card — smaller-magnitude drift than the full-value 90° card, at every speed and range,
because only part of the wind is perpendicular to the shot:

```
Wind Card (zero: 100 yd, MIL) — wind angle 45° (wind-FROM: 0=head, 90=right, 180=tail, 270=left)
┌──────────┬──────────┬──────────┬──────────┬──────────┐
│Range (yd)│       5 mph │      10 mph │      15 mph │      20 mph │
├──────────┼──────────┼──────────┼──────────┼──────────┤
│      100 │     -0.1 │     -0.1 │     -0.2 │     -0.3 │
│      200 │     -0.1 │     -0.3 │     -0.4 │     -0.6 │
│      300 │     -0.2 │     -0.4 │     -0.7 │     -0.9 │
└──────────┴──────────┴──────────┴──────────┴──────────┘
```

### Wind Direction Entry (degrees, clock positions) — MBA-1367

Every `--wind-direction` flag (`trajectory`, `monte-carlo`, `come-ups`, `lead`,
`range-table`, `compare`, `profile save`) accepts the same three forms:

* **Bare number — degrees** (unchanged, always): wind-FROM convention, `0` = headwind,
  `90` = from the right, `180` = tailwind, `270` = from the left. Any value that parsed
  before still parses identically.
* **`<H>oc`** — a 1-12 o'clock shooter-relative clock position: `3oc` = 90°, `6oc` =
  180°, `9oc` = 270°, `12oc` = headwind = 0°.
* **`<H>h<MM>`** — clock position with minutes: `10h30` = 315° (minutes count 0.5° each:
  `(H % 12) × 30 + MM × 0.5`).
* **`<H>:<MM>`** — the colon spelling of the same thing (`10:30` = 315°), legal on
  **standalone flags only**.

The markers are mandatory: a bare `3` is 3 degrees, never 3 o'clock — that ambiguity is
why unmarked clock numbers are not accepted. Inside `--wind-segment`, whose grammar is
colon-delimited (`SPEED:ANGLE:DIST[:VERTICAL]`), the ANGLE field takes the **colon-free**
marked forms only: `--wind-segment 10:3oc:400` is a 10 mph wind from 3 o'clock out to
400 yd, while `10:30:400` keeps its historical meaning (10 mph from 30 degrees).
Malformed clock tokens (hour outside 1-12, minutes outside 0-59) are hard errors naming
the rule.

```bash
# Identical solves, three spellings:
ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --wind-speed 10 --wind-direction 90
ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --wind-speed 10 --wind-direction 3oc
ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --wind-speed 10 --wind-segment 10:3oc:1000
```

Clock positions are **shooter-relative by definition** (12 o'clock = dead ahead), so they
cannot be combined with the earth-fixed compass mode described below.

**Sentinel fix (behavior change, deliberate):** `trajectory` used to treat
`--wind-direction 0` as "not set" and let a `--location` CSV's `WIND_DIR` column replace
it. Explicit presence now decides: an explicit `--wind-direction 0` (or `12oc`, which
maps to 0°) wins over the CSV; omitting the flag still inherits the CSV value exactly as
before.

#### Earth-Fixed Compass Bearings (`--wind-ref compass`) — MBA-1368

Orthogonal to the SYNTAX above is the REFERENCE FRAME the direction is entered in
(`--wind-ref {shooter|compass}` on `trajectory` and `monte-carlo`):

* **`shooter`** (default) — wind-FROM angles relative to the line of fire, exactly as
  always. Byte-identical to omitting the flag.
* **`compass`** — every wind direction the run consumes (`--wind-direction`, a
  `--location` CSV's `WIND_DIR`, and each `--wind-segment` angle) is an **absolute
  earth-fixed bearing** (0 = north, 90 = east), the Vortex Wind Bearing Capture / Lapua
  Ballistics convention. It is re-referenced against the shot azimuth once at the input
  boundary — `relative = bearing − shot azimuth`, wind-FROM on both sides, normalized to
  [0°, 360°) — so the solver's physics never sees a bearing.

Compass mode **requires `--shot-direction`** (the shot's own compass bearing; a hard
error names the flag when missing — silently treating bearings as shooter-relative is
exactly the bug this mode exists to prevent). An explicit `--shot-direction 0` is a real
bearing (due north), not "unset". Clock positions are rejected in compass mode: a
"3 o'clock bearing" is meaningless, because clock positions are shooter-relative by
definition.

```bash
# Wind FROM the north-east (bearing 45) on a shot fired due east (azimuth 90):
# relative = 45 - 90 = -45 -> 315 = wind from the shooter's front-left.
ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --wind-speed 10 \
  --wind-ref compass --shot-direction 90 --wind-direction 45

# Sanity anchor: wind FROM north, shot due north = pure headwind (relative 0).
ballistics trajectory ... --wind-ref compass --shot-direction 0 --wind-direction 0
```

On `monte-carlo` the base wind bearing is converted **before** any dispersion sampling,
so `--wind-direction-std` disperses around the converted shooter-relative direction
(`monte-carlo` gains `--shot-direction` for exactly this; it models no Coriolis). The
solve-json v1 surface takes the same mode as the optional `wind.wind_reference` field
(`"compass"` requires an explicit `shot.shot_azimuth_rad`; see docs/SOLVE_JSON_V1.md).
Commands without a shot-azimuth concept (`come-ups`, `lead`, `wind-card`,
`range-table`, `compare`, `dsf`) stay shooter-relative only, as do the FFI bindings
(convert bearings yourself as `bearing − shot azimuth` before the call).

### Vertical Wind

Model wind with a vertical component — a thermal updraft, a downdraft off the lee side of
a ridge, orographic lift along a slope — as a direct shift in point of impact, on top of
the horizontal drift `--wind-speed`/`--wind-direction` already model. Available via
`--wind-vertical <SPEED>` on both the `trajectory` and `monte-carlo` subcommands, and as an
optional 4th field on `--wind-segment`.

**Sign convention:** positive = updraft, negative = downdraft. An updraft raises point of
impact; a downdraft lowers it. Units follow `--units` the same way `--wind-speed` does: mph
under imperial, m/s under metric. Default `0` (no vertical wind, bit-identical to a solve
without the flag).

```bash
--wind-vertical 5     # metric: 5 m/s updraft, raises POI
--wind-vertical -8    # imperial: 8 mph downdraft, lowers POI
```

**Downrange segments — the optional 4th field.** `--wind-segment` accepts an optional 4th
colon-separated field: `SPEED:ANGLE:UNTIL_DISTANCE[:VERTICAL]`, letting each segment carry
its own vertical wind:

```bash
ballistics trajectory -v 2700 -m 168 -d 0.308 --bc 0.5 --max-range 1000 \
  --wind-segment 8:90:300:5 \
  --wind-segment 12:90:1000:10
```

- **VERTICAL is always m/s, positive = updraft — regardless of `--units`.** This is a
  deliberate asymmetry with SPEED (and UNTIL_DISTANCE), which *do* follow `--units`: SPEED
  matches the display system so a wind-meter reading in mph can be pasted directly, but
  there's no comparably universal display convention for vertical wind, so it's pinned to
  the engine's native m/s to keep the field unambiguous no matter what `--units` is active.
- Omitting the 4th field is unchanged from before — that segment's vertical wind is `0.0`.
  Every existing 3-field `--wind-segment SPEED:ANGLE:UNTIL_DISTANCE` string keeps working
  exactly as it did.

**Shear pass-through.** `--enable-wind-shear`'s boundary-layer models (logarithmic /
power-law / Ekman spiral) scale **horizontal** wind speed with altitude only. Vertical wind
passes through unscaled wherever shear is layered on top of it — a 5 m/s updraft is a
5 m/s updraft at every altitude the shot climbs through, shear on or off. (`--wind-segment`
itself is not compatible with `--enable-wind-shear` — see below — so in practice this rule
matters for the scalar `--wind-speed`/`--wind-vertical` + `--enable-wind-shear`
combination.)

**Precedence.** `--wind-segment`, when given, overrides the scalar wind entirely —
including `--wind-vertical`. Once segments are set, each segment's own 4th field is the
only source of vertical wind; the scalar `--wind-vertical` value is ignored (a note is
printed to stderr, the same override behavior `--wind-speed`/`--wind-direction` already have
against `--wind-segment`).

**Drag-symmetry fact.** Aerodynamic drag doesn't distinguish "up" from "sideways" — a
projectile flying into a 5 m/s updraft is deflected vertically by essentially the same
mechanism that deflects it laterally in a 5 m/s crosswind. `tests/vertical_wind.rs` locks
this as a regression gate: a 5 m/s updraft's vertical deflection and a 5 m/s crosswind's
lateral deflection must agree within 5% at 300 m and 600 m. The two solves in fact agree
far more tightly than the bound requires — about 0.001% at 600 m in the current build.

**Monte Carlo caveat:** `--wind-vertical` on `monte-carlo` sets the *base* (mean) vertical
wind shared by every simulated shot — it is a systematic input, not a dispersion source.
There is no `--wind-vertical-std`; each sampled wind draw carries the same vertical
component through un-dispersed, while horizontal wind speed/direction still vary per
`--wind-std`/`--wind-direction-std`. Expect `--wind-vertical` to shift the whole impact
cloud's mean point of impact without changing its reported spread — the same pattern as
`--cant` (see the Canted Shooting Monte Carlo caveat above).

**Worked example** (build first with `cargo build`):

```bash
# Calm baseline
./target/debug/ballistics trajectory --units metric -v 823 -m 10.9 -d 7.82 --bc 0.5 \
  --drag-model g7 --auto-zero 500 --max-range 600 \
  --sample-trajectory --sample-interval 100 -o csv --full

# Same load, 5 m/s updraft
./target/debug/ballistics trajectory --units metric -v 823 -m 10.9 -d 7.82 --bc 0.5 \
  --drag-model g7 --auto-zero 500 --max-range 600 --wind-vertical 5 \
  --sample-trajectory --sample-interval 100 -o csv --full
```

Calm baseline:

```
distance_m,drop_m,drift_m,velocity_m/s,energy_J,time_s
0.00,0.05,0.00,823.00,3691.44,0.0000
100.00,-0.30,0.00,792.49,3422.84,0.1238
200.00,-0.49,0.00,762.61,3169.54,0.2525
300.00,-0.51,0.00,733.35,2931.05,0.3862
400.00,-0.36,0.00,704.75,2706.89,0.5253
500.00,0.00,0.00,676.82,2496.58,0.6701
600.00,0.57,0.00,649.57,2299.59,0.8209
```

5 m/s updraft — `drop_m` runs smaller at every range past the zero, and by 600 m point of
impact sits about 8 cm higher than the calm run (`0.49` vs `0.57`):

```
distance_m,drop_m,drift_m,velocity_m/s,energy_J,time_s
0.00,0.05,0.00,823.00,3691.44,0.0000
100.00,-0.25,0.00,792.49,3422.85,0.1238
200.00,-0.41,0.00,762.61,3169.55,0.2525
300.00,-0.44,0.00,733.35,2931.06,0.3862
400.00,-0.30,0.00,704.75,2706.89,0.5253
500.00,-0.00,0.00,676.82,2496.56,0.6701
600.00,0.49,0.00,649.57,2299.55,0.8209
```

### Zero Calculation

Calculate sight adjustments for specific distances:

```bash
# Calculate zero for 200 yards
./ballistics zero -v 2700 -b 0.475 -m 168 -d 0.308 --target-distance 200

# With custom sight height
./ballistics zero -v 2700 -b 0.475 -m 168 -d 0.308 \
  --target-distance 300 \
  --sight-height 0.055  # 2.2 inches in yards

# Metric
./ballistics zero --units metric -v 823 -b 0.475 -m 10.9 -d 7.82 \
  --target-distance 200  # 200 meters
```

Output provides:
- Zero angle in degrees
- MOA adjustment
- Mrad adjustment
- Maximum ordinate

#### Solving Range From a Stored Angle (`--from-angle`)

Hornady and Kestrel 4DOF treat the zero angle as the *portable* quantity: capture it once
(from this command's own `zero_angle_degrees` output, or from `trajectory`'s auto-zero echo
below), then recover the zero range(s) it implies later — on a different day, under
different conditions — without re-solving from a target distance. `--from-angle <DEGREES>`
(MBA-1402) does this by running the trajectory at the given bore angle and reporting BOTH
line-of-sight crossings it finds.

**This is *not* the exact inverse of the default `--target-distance` mode.** A rifle sighted
above the bore generally crosses a level line of sight TWICE on a rising shot: once
ascending, close to the muzzle (the near zero), and once descending past the apex (the far
zero). Both are equally real zero ranges for the *same* bore angle — this is the classic
25/300-yard battle-zero relationship: a rifle zeroed at 25 yd is, for that same bore angle,
also zeroed again around 300 yd. `--target-distance`'s forward solve returns whichever root
matches the distance you asked it to hit — the near root for a short/flat zero like 25 or 50
yd, the far root for a conventional 100+ yd zero — so `--from-angle` cannot know in advance
which one you mean and reports both instead.

**The handoff only reproduces the original zero under the original conditions.** A stored bore
angle is portable, but the RANGE it implies is not a property of the angle alone: it depends on
the wind, the bore height, and the zero-day atmosphere in force when you feed it back. Recover a
range under a different crosswind, a different scope height, or a materially different density
altitude and you will get a different — and correct — answer for those conditions, not the
number you originally zeroed at. That is the point of the feature, but it surprises people who
expect a round trip to be exact. Re-supply the same `--wind-*`, `--sight-height`, and zero-day
atmosphere flags if you want the original range back.

**`max_ordinate` is conditional** (MBA-1419): present in JSON/CSV — and printed in the
table — only when the diagnostics trajectory actually turned over inside its solve envelope.
When the solve is still climbing at its limit (a steep launch under a slick drag model), the
table says `still climbing at limit` and the JSON key / CSV row are omitted rather than
carrying the truncation height, which is not a max ordinate. Before 0.30.1 the JSON and CSV
emitted that truncation height unconditionally.

**`primary_crossing`** (MBA-1419) names which of the two roots the single `zero_range` value is:
`"far"` when a far crossing was found (the conventional "sighted in at D yards" meaning), and
`"near"` when it was not, which happens when the far crossing lies beyond the search envelope.
Both roots are also reported individually as `near_zero_range` / `far_zero_range` whenever they
exist, so a consumer never has to infer which answer it received.

**`--drag-model`** (0.30.1): the zero command previously had no drag-model selection — a G7
BC silently ran against the G1 reference while every sibling command (including the browser
terminal's own zero) could choose. Default `g1` preserves the old behavior byte-for-byte.

**Browser terminal**: the terminal's `zero` previously accepted NO atmosphere flags at all
(both its solve branches ran on the default atmosphere). As of 0.30.1 it takes
`--temperature`/`--pressure`/`--pressure-type`/`--humidity`/`--altitude`/`--density-altitude`,
feeding both the target-distance and `--from-angle` solves, with the density-altitude
supersede notice riding the output text. The terminal's `true-velocity`, `estimate-bc`, and
`lead` also gained `--pressure-type`, and `trajectory` gained `--zero-density-altitude`.

`--from-angle` and `--target-distance` are mutually exclusive — supply exactly one:

```bash
# Solve the angle for a 200-yard zero...
./ballistics zero -v 2700 -b 0.475 -m 168 -d 0.308 --target-distance 200
# ... then, later, recover the range(s) that same stored angle (0.0997°) implies. One angle
# generally gives TWO zeros: this one reports a near crossing at 37.3 yd (the bullet rising
# through the line of sight) and a far crossing at 199.9 yd (descending back through it).
# The 200 yd zero you asked for is the far one — the same relationship that makes a 25 yd
# zero and a ~300 yd zero the same bore angle:
./ballistics zero -v 2700 -b 0.475 -m 168 -d 0.308 --from-angle 0.0997
```

`--output table` prints a "Near Zero"/"Far Zero" line for each crossing, or "not within the
solved range" when that crossing doesn't occur (e.g. a short/shallow angle whose far crossing
lies beyond the solve envelope). `--output json`/`--output csv` add `near_zero_range` and
`far_zero_range` values (display units), each present only when that crossing was found, plus
a `zero_range` value kept for continuity (the far crossing when present, else the near one) —
alongside the `zero_angle_degrees`/`zero_angle_moa`/`zero_angle_mrad` echo of the angle you
supplied. `sight_adjustment_moa`, `max_ordinate`, and `point_blank_range` are computed the
same way as the forward mode, referenced against that same `zero_range` value. Round-trip
precision — recovering the ORIGINAL distance as one of the two reported crossings — follows
from the forward solver's own convergence granularity: well under 1% at typical zero
distances, with the shortest/near-tangent zeros (very flat trajectories, e.g. a 100 yd zero
for a high-BC load) at the loose end of that range.

### Load Comparison (`compare`)

Run several loads through identical conditions and see them side by side. Each load is
zeroed independently at the shared `--zero-distance`, then solved twice (a no-wind pass
for pure drop, a wind pass for drift), exactly like `range-table`:

```bash
# Two loads by inline spec: NAME:DRAG:BC:MASS:VELOCITY[:DIAMETER]
ballistics compare \
  --load "175 SMK:g7:0.243:175:2650" \
  --load "168 ELD-M:g7:0.523:168:2700" \
  --zero-distance 100 --end 800 --step 100

# Mix inline specs with saved profiles; MOA adjustments; machine output
ballistics compare --load "Factory:g1:0.475:168:2700" --profile my-match-load \
  --zero-distance 100 --adjustment-unit moa -o json
```

Load-spec fields follow the session `--units`: `MASS` is grains (imperial) or grams
(metric), `VELOCITY` fps or m/s, and the optional `DIAMETER` inches or mm (defaulting to
.308 in / 7.82 mm). `DRAG` is any of `g1`/`g2`/`g5`/`g6`/`g7`/`g8`/`gi`/`gs`/`ra4`, and
`NAME` may not contain `:`. Between 2 and 8 loads are accepted, from `--load` and/or
`--profile` in any combination. A saved
profile's velocity-BC segments and custom Cd(Mach) drag curve (e.g. from an `.a7p`
import) ARE consumed here — they drive both the load's zeroing and its trajectory runs,
and such loads are tagged `[BC segments]` / `[custom drag curve]` in the table legend
(and flagged in JSON). Inline `--load` specs use the scalar BC.

The table shows per-load drop (`--adjustment-unit <mil|moa|smoa|iphy|clicks>`, the
elevation axis) and drift (`--windage-unit`, falling back to `--adjustment-unit` when
omitted — MBA-1410; see [Turret Adjustment Units](#turret-adjustment-units)), and velocity
at each range. `clicks` requires `--elevation-click-value`/`--windage-click-value` — since
`compare` has no single `--profile`, these always come from the explicit flags, never a
per-load profile's own graduation. `-o json` adds linear drop/drift, energy, time of
flight, each load's zero angle, and per-row deltas against load #1 (`delta_drop`,
`delta_drift`, `delta_velocity`, `delta_energy` — zero for the baseline itself); `-o csv`
emits one column group per load (names sanitized for CSV). PDF output is not supported
for this command.

### Powder Temperature Velocity (`powder`)

Resolve the powder-temperature-adjusted muzzle velocity without running a trajectory.
The physics is the exact resolution the `trajectory` and `lead` solvers apply
internally (one shared implementation). One flag difference: `powder` always applies
the linear model, while `trajectory`/`lead` only apply it when you also pass
`--use-powder-sensitivity` (a measured curve applies there unconditionally) — carry
that flag along or your trajectory will fly the nominal velocity:

```bash
# Linear model: 2800 fps load (measured at the default 70 °F reference), 40 °F day
ballistics powder -v 2800 --temperature 40
```

```
Powder Temperature Velocity
===========================
  Model:              linear, 1.00 fps/°F
  Reference temp:     70.0 °F
  Nominal velocity:   2800.0 fps

  Shot temp:          40.0 °F
  Resolved velocity:  2770.0 fps  (-30.0)
```

```bash
# Measured curve (overrides the linear model): interpolate at a 55 °F powder temp
ballistics powder --powder-temp-curve "40:2620,70:2700,100:2760" --powder-temp 55

# Velocity ladder across a temperature range, with muzzle energy
ballistics powder -v 2700 -m 168 --sweep 20:110:30
```

```
Powder Temperature Velocity
===========================
  Model:              linear, 1.00 fps/°F
  Reference temp:     70.0 °F
  Nominal velocity:   2700.0 fps

   Temp (°F)  Velocity (fps)   Shift (fps)   Energy (ft·lb)
        20.0          2650.0         -50.0             2619
        50.0          2680.0         -20.0             2679
        80.0          2710.0          10.0             2739
       110.0          2740.0          40.0             2800
```

Flags follow the session `--units` (fps + °F imperial; m/s + °C metric) and carry the
same meanings as on `trajectory`/`lead`: `--powder-temp-sensitivity` defaults to
1.0 fps/°F (0.54864 m/s/°C); `--powder-temp` is the linear model's *reference*
temperature — the temperature the stated `-v` velocity was measured at, defaulting to
70 °F — or, with `--powder-temp-curve`, the powder temperature the curve is
interpolated at (defaulting to `--temperature`, i.e. powder at air temperature; the
curve is clamped at its endpoints, never extrapolated). With a curve the sweep
temperatures are powder temperatures, and `-v` is optional — the curve supplies the
velocity, `-v` only anchors the reported shift. `-m/--mass` (grains/grams) adds muzzle
energy (ft·lb / J). Output: table (default), `-o json`, or `-o csv`; PDF is not
supported for this command.

### Free Recoil (`recoil`)

Free recoil energy, velocity, and impulse from SAAMI's own momentum-balance formula
("Gun Recoil - Technical", Rev. 7/9/2018, freely downloadable from saami.org): the
firearm's momentum is equal and opposite to the ejecta (bullet) plus propellant gas
momentum leaving the muzzle, with the propellant gas *mass* equated to the powder
charge weight (SAAMI's own simplifying assumption — the gas itself is hard to weigh).

**Gas-velocity convention.** The propellant gas leaves the muzzle faster than the
bullet. This command defaults to **SAAMI's own type-keyed multiplier** (`V_gas = f *
V_muzzle`), sourced from the document's own 1929-era British ballistic testing:

| `--firearm-type` | `f` |
|---|---|
| `rifle` (default) | 1.75 |
| `pistol` | 1.50 |
| `shotgun-average` | 1.50 |
| `shotgun-long` | 1.25 |

This is preferred over the popular fixed "~4700 fps for smokeless powder" rule of
thumb quoted by some reloading references because it scales with muzzle velocity (a
fixed constant over- or under-states gas momentum for very slow or very fast loads)
and is keyed to firearm class the way SAAMI itself publishes it. Both escape hatches
are still exposed: `--gas-velocity-factor <F>` overrides the SAAMI factor with your
own ratio, and `--gas-velocity <VEL>` pins an absolute gas velocity (fps/m·s⁻¹)
independent of muzzle velocity — e.g. to reproduce the fixed-constant convention or
plug in a chronographed value. Precedence: `--gas-velocity` > `--gas-velocity-factor`
> `--firearm-type`.

```bash
# .308 Win: 168 gr bullet, 43 gr charge, 2700 fps, 8.5 lb rifle
ballistics recoil -b 168 -c 43 -v 2700 -f 8.5
```

```
Free Recoil (SAAMI Momentum Balance)
=====================================
  Firearm weight:      8.50 lb
  Bullet weight:       168.0 gr
  Charge weight:       43.0 gr
  Muzzle velocity:     2700.0 fps
  Gas velocity model:  saami-rifle  (4725.0 fps)

  Recoil velocity:     11.04 fps
  Recoil energy:       16.09 ft-lb
  Recoil impulse:      2.916 lb-s
```

Flags follow the session `--units`: bullet/charge weight in grains (imperial) or
grams (metric); firearm weight in **pounds** (imperial) or **kilograms** (metric) —
note this differs from bullet/charge weight's grains/grams, matching how a firearm is
actually weighed. `-b/--bullet-weight`, `-c/--charge-weight`, `-v/--velocity`, and
`-f/--firearm-weight` are all required (no profile integration — this is a standalone
calculator). Output: table (default), `-o json`, or `-o csv`; PDF is not supported.

### Reference Drag Curve (`drag-curve`)

Prints a built-in reference drag function as `(Mach, Cd)` data, so a chart or an audit can use
the engine's own numbers instead of a re-vendored copy that drifts from them as the tables are
refined.

```bash
ballistics drag-curve --drag-model g7            # table (default)
ballistics drag-curve --drag-model g7 -o csv     # mach,cd
ballistics drag-curve --drag-model gs  -o json   # with provenance and domain
```

All nine models are available: `g1`, `g2`, `g5`, `g6`, `g7`, `g8`, `gi`, `gs`, `ra4`. Output is
`table` (default), `csv`, or `json`; PDF is not supported.

Points are emitted **verbatim from the table** — no resampling and no interpolation onto a
uniform grid — so what you plot is exactly what the solver interpolates.

**The tables do not share a Mach domain.** Most run to Mach 5; `gs` and `ra4` stop at Mach 4.
Read the emitted points (or `mach_min`/`mach_max` in the JSON form) rather than assuming a range.

These are the public-domain Aberdeen/BRL reference functions as tabulated in McCoy, *Modern
Exterior Ballistics*, plus the British RA 1929 function for `ra4`. Each `data/*.csv` in the
source tree carries its own provenance header.

**This is the reference curve for the model's standard projectile**, not the effective Cd a
particular bullet flew. For that — form-factor scaled, with segmented-BC band steps — use
`trajectory --full --with-drag-coefficient`. A chart typically wants both: the reference curve as
a backdrop and the load's own curve over it.

Library callers can reach the same data via `ballistics_engine::drag::reference_drag_table`.

**Browser terminal:** `drag-curve` works there too, with the same flags. Both surfaces print the
same shared formatter (`drag::format_reference_drag_curve`), so their table, CSV and JSON output
are byte-identical by construction rather than by discipline.

### Power Factor (`power-factor`)

Power factor (PF) — `bullet weight × velocity / 1000` — gates ammunition into scoring
categories across essentially every American practical/action shooting sport. This
command reports PF plus a per-organization pass/fail against a small, cheap-to-update
data table of published rulebook minimums:

| Organization | Class | Min PF | Velocity limit | Source |
|---|---|---|---|---|
| USPSA | Minor | 125 | — | USPSA Competition Rules, March 2026, Sec. 5.6 |
| USPSA | Major | 165 | — | USPSA Competition Rules, March 2026, Sec. 5.6 |
| IDPA | BUG | 95 | — | 2024 IDPA Rulebook, Sec. 8.3.4.5 |
| IDPA | Stock Revolver / CCP | 105 | — | 2024 IDPA Rulebook, Sec. 8.3.4.3 |
| IDPA | SSP / ESP / CO | 125 | — | 2024 IDPA Rulebook, Sec. 8.3.4.1 |
| IDPA | PCC | 135 | — | 2024 IDPA Rulebook, Sec. 8.3.4.6 |
| IDPA | Enhanced Revolver | 155 | — | 2024 IDPA Rulebook, Sec. 8.3.4.4 |
| IDPA | CDP | 165 | — | 2024 IDPA Rulebook, Sec. 8.3.4.2 |
| SASS | Smokeless, Pistol/Revolver | 60 | 400-1000 fps | SASS CAS Shooter's Handbook, Version 28, Jan 2026 |
| SASS | Smokeless, Rifle | 60 | 400-1400 fps | SASS CAS Shooter's Handbook, Version 28, Jan 2026 |

**Truncation.** USPSA's and IDPA's own rulebooks both say the scored PF *truncates*
decimals rather than rounding (USPSA Sec. 5.6 Rule 36: "a result of 124.9999 is not
125"; IDPA Sec. 8.3.4.7: "ignore numbers to the right of the decimal"). This command
reports both `Power factor (raw)` (unrounded) and `Power factor (scored)`
(`floor(weight × velocity / 1000)`, applied uniformly to every organization including
SASS, whose own handbook doesn't spell out a truncation rule) — pass/fail is always
judged against the scored value.

**Boundary semantics.** Every rulebook cited above uses "meets or exceeds" / "at
least" / "not less than" language, never a strict "greater than" — so a load landing
*exactly on* a threshold **passes**. Every threshold here is `>=`, every velocity cap
is inclusive (`<=`/`>=`), never a strict inequality.

```bash
# 9mm major/minor check across every organization
ballistics power-factor -w 147 -v 900
```

```
Power Factor
============
  Weight:              147.0 gr
  Velocity:            900.0 fps
  Power factor (raw):  132.30
  Power factor (scored): 132

  Org       Class                         Min PF    Pass      Velocity limit
  USPSA     Minor                            125    PASS                   -
  USPSA     Major                            165    FAIL                   -
  IDPA      BUG                               95    PASS                   -
  IDPA      Stock Revolver / CCP             105    PASS                   -
  IDPA      SSP / ESP / CO                   125    PASS                   -
  IDPA      PCC                              135    FAIL                   -
  IDPA      Enhanced Revolver                155    FAIL                   -
  IDPA      CDP                              165    FAIL                   -
  SASS      Smokeless - Pistol-Revolver       60    PASS        400-1000 fps
  SASS      Smokeless - Rifle                 60    PASS        400-1400 fps
```

```bash
# Filter to one organization
ballistics power-factor -w 147 -v 900 --organization sass
```

Flags follow the session `--units` (weight in grains/grams, velocity in fps/m·s⁻¹),
but power factor is intrinsically a grains/fps quantity by every rulebook's own
definition, so metric inputs are converted internally before the PF arithmetic runs
— the echoed weight/velocity stay in your chosen display units. `--organization
<uspsa|idpa|sass>` filters the table to one organization (case-insensitive); omit it
to see all three. Output: table (default), `-o json`, or `-o csv` — the JSON/CSV
`thresholds` rows carry `pf_pass`, `velocity_pass` (`null` when the organization has
no velocity constraint), and the combined `pass`. PDF is not supported.

### Monte Carlo Simulation

Statistical analysis with parameter variations:

```bash
# Basic Monte Carlo
./ballistics monte-carlo -v 2700 -b 0.475 -m 168 -d 0.308 -n 1000

# With variations and target distance
./ballistics monte-carlo \
  -v 2700         # Base velocity (fps)
  -b 0.475        # Base BC
  -m 168          # Mass (grains)
  -d 0.308        # Diameter (inches)
  -n 1000         # Simulations
  --velocity-std 10    # Velocity std dev
  --angle-std 0.5      # Angle std dev
  --bc-std 0.01        # BC std dev
  --wind-std 2         # Wind std dev
  --target-distance 600 # For hit probability
```

### Confidence-Controlled Sampling (`--adaptive`) — MBA-1352

Every `monte-carlo --target-distance` run reports a hit probability, but a bare point estimate
does not say how much sampling backs it -- a 20% estimate from 30 trials and a 20% estimate
from 30,000 trials look identical unless the report also states an interval. This section
covers two ways to get one: an **additive interval on the existing fixed-`n` run**, and an
**opt-in `--adaptive` mode** that keeps sampling until the interval is as tight as you asked
for.

**Flags:**

| Flag | Default | Meaning |
| --- | --- | --- |
| `--confidence <90\|95\|99>` | `95` | Confidence level of the reported interval. Applies to the fixed-count Wilson line/JSON key below AND to `--adaptive`'s stopping rule and interval. Any other value is a usage error naming `--confidence` and listing the three choices. |
| `--adaptive` | off | Switch from a fixed `--num-sims` to confidence-controlled sampling (below). Incompatible with `--wez`, which stays fixed-count. |
| `--target-ci-half-width <p>` | `0.02` | `--adaptive` only. Stop once the interval's half-width is at or below this many **probability units** (`0.02` = the hit probability is known to within about ±2 percentage points). Must be finite and greater than zero. Requires `--adaptive`. |
| `--min-samples <n>` | `1000` | `--adaptive` only. Never stop before this many trials, even if the interval already looks tight. Matches `--num-sims`'s own legacy default. Requires `--adaptive`. |
| `--max-samples <n>` | `100000` | `--adaptive` only. Never run more than this many trials. Hitting it is an honest, successful run (exit `0`), not an error. Requires `--adaptive`. |
| `--mc-batch-size <n>` | `500` | `--adaptive` only. How many trials to run between convergence checks. Requires `--adaptive`. |
| `--seed <u64>` | none (fresh randomness each run) | Seeds the Monte Carlo RNG for reproducible output, on **both** the default fixed-count path and `--adaptive`. |

`--adaptive` ignores `--num-sims` (`--min-samples` is the floor instead -- the two modes have
different stopping semantics) and `--wind-direction-std` (the adaptive driver disperses no wind
direction, matching `run_monte_carlo_with_wind`); both are silently unused rather than errors,
since they still apply to every *other* `monte-carlo` mode on the same command line.

**Fixed-count mode gains one line.** No existing number or field changes shape -- `--seed` is
new too, added only so this example reproduces exactly if you run it yourself; the run behaves
identically to an unseeded one, just with fixed randomness -- but the hit-probability line now
has a companion:

```bash
./ballistics monte-carlo -v 2700 -b 0.475 -m 168 -d 0.308 -n 1000 \
  --velocity-std 10 --angle-std 0.5 --bc-std 0.01 --wind-std 2 \
  --target-distance 600 --seed 42
```
```
╔════════════════════════════════════════╗
║      MONTE CARLO RESULTS               ║
║      1000 simulations                   ║
╠════════════════════════════════════════╣
║ Mean Range:          480.36 m          ║
║ Std Dev Range:       311.66 m          ║
║ Mean Impact Vel:     564.32 m/s        ║
║ Std Dev Velocity:    141.30 m/s        ║
║ CEP (arrivals):        3.73 m          ║
║ Target Shortfall:      61.9 %          ║
║ Hit Probability:        0.2 %          ║
╚════════════════════════════════════════╝
Hit probability 95% CI: [0.001, 0.007] (Wilson, n=1000)
```

This is the textbook point precisely because the estimate is small: "0.2%" alone reads as
precise, but the 95% Wilson interval says the true probability could plausibly be anywhere from
about 0.05% to 0.73% at only 1000 trials -- a better-than-13x range the bare percentage hides.
`-o json`/`-o full` (same command, `-o json` added) gain the identical fact as an additive key,
`hit_probability_ci`, alongside the unchanged existing keys:

```json
{
  "mean_range": 480.35763209156374,
  "std_range": 311.65807211958315,
  "mean_impact_velocity": 564.3227501383241,
  "std_impact_velocity": 141.29683010644374,
  "cep": 3.7330017194509684,
  "target_shortfall_fraction": 0.619,
  "hit_probability": 0.002,
  "hit_probability_ci": {
    "method": "wilson_fixed_n",
    "confidence_percent": 95,
    "low": 0.000548643588120728,
    "high": 0.007262807863492177,
    "samples": 1000
  }
}
```

`hit_probability_ci` is `null` (present, not absent, matching every other optional field this
command's JSON already carries) whenever `hit_probability` itself is -- i.e. whenever
`--target-distance` was not given.

**`--adaptive` decides its own sample size.** Instead of guessing `--num-sims`, tell it how
precise you need the answer and it samples in batches (`--mc-batch-size`) until the interval
reaches `--target-ci-half-width` or `--max-samples` is exhausted:

```bash
./ballistics monte-carlo -v 2700 -b 0.475 -m 168 -d 0.308 \
  --velocity-std 5 --angle-std 0.1 --bc-std 0.005 --wind-std 1.5 --target-distance 400 \
  --adaptive --confidence 90 --target-ci-half-width 0.03 \
  --min-samples 200 --max-samples 20000 --seed 42
```
```
Adaptive Monte Carlo (confidence-controlled sampling, MBA-1352)
  method:            anytime_beta_binomial_mixture_cs_v1
  confidence level:  90%
  hit probability:   0.1912  (90% CI: [0.1645, 0.2200])
  samples:           2500  (attempts: 2500, arrivals: 1725)
  stop reason:       target_half_width_met
  hit radius:        0.3000 m
  target distance:   365.7600 m
  impact velocity:   mean 560.851 m/s, std 31.466 m/s
  drop at target:    mean 0.3291 m, std 0.4407 m
  drift at target:   mean -0.0043 m, std 0.3236 m
```

It ran 2500 of the 20,000 permitted trials -- exactly what the requested 90% / ±0.03 precision
needed on this load, no more -- and stopped for a data-dependent reason
(`target_half_width_met`), not because a fixed count ran out. `samples`, `attempts`, and
`arrivals` are three different counts, in that order: `attempts` is every trial drawn,
`samples` (`= attempts` here) is every trial that produced an outcome and is the `n` behind the
hit probability and its interval, and `arrivals` (`1725`, well under `samples`) is only the
trials that reached the target plane -- the population behind the six velocity/drop/drift
statistics below `stop reason`. `--output full` -- `json` is accepted as an alias for the same
value, so `-o json` works too -- prints the identical report as its versioned wire form,
`AdaptiveMcReportV1`, verbatim:

```json
{
  "schema_version": 1,
  "method": "anytime_beta_binomial_mixture_cs_v1",
  "assumptions": [
    "Sampling uncertainty only: intervals cover Monte Carlo sampling error, not model error in the trajectory solver or its inputs.",
    "Anytime-valid stopping: the beta-binomial mixture confidence sequence keeps its coverage guarantee despite stopping the moment the target half-width is met.",
    "Input dispersions are the independent normal distributions declared in MonteCarloParams; correlations between inputs are not modeled.",
    "Continuous statistics are streaming Welford moments over trials that reached the target plane, reported with sample (n-1) standard deviations; hit probability's denominator includes all trials."
  ],
  "confidence_percent": 90,
  "hit_probability": 0.1912,
  "ci_low": 0.16448890323739063,
  "ci_high": 0.2199595698086259,
  "samples": 2500,
  "attempts": 2500,
  "arrivals": 1725,
  "stop_reason": "target_half_width_met",
  "hit_radius_m": 0.3,
  "target_distance_m": 365.76,
  "mean_impact_velocity_mps": 560.8513010242867,
  "std_impact_velocity_mps": 31.46555841729584,
  "mean_drop_at_target_m": 0.3290765670086509,
  "std_drop_at_target_m": 0.44069648856048765,
  "mean_wind_drift_at_target_m": -0.004288891866863569,
  "std_wind_drift_at_target_m": 0.32356953067182526
}
```

`-o statistics` (the CSV format) is not available under `--adaptive`: the adaptive driver keeps
no per-trial ranges or velocities to tabulate (that constant-memory property, unbounded trial
count at the cost of the CSV/full trial history, is the point of it), so this is a named usage
error rather than a silently-empty or fabricated CSV.

**Honest limits, stated in the payload itself.** Every `AdaptiveMcReportV1` carries its own
`assumptions` array so a consumer rendering just the numbers still has the caveats attached.
The first two, verbatim:

> Sampling uncertainty only: intervals cover Monte Carlo sampling error, not model error in the
> trajectory solver or its inputs.
>
> Anytime-valid stopping: the beta-binomial mixture confidence sequence keeps its coverage
> guarantee despite stopping the moment the target half-width is met.

The first is the same caveat every interval in this file carries: it bounds sampling noise, not
whether the underlying physics model or its inputs are themselves correct. The second is
specific to `--adaptive` -- it is *not* safe to run the fixed-count Wilson interval above
repeatedly and stop the instant it "looks" tight; that inflates the true error rate with every
peek. `--adaptive` is safe to stop early because it uses a different interval (an anytime-valid
confidence sequence) built for exactly that use, at the cost of a strictly wider interval than
Wilson would report at the same `n` -- paying that width is what buys the early stop.

**Units rule.** `--target-ci-half-width`, and every low/high/`ci_low`/`ci_high` bound this
section's flags and JSON keys report, are in **probability units** (`0`-`1`), the same space
`hit_probability` itself is in -- never a percentage, even though the fixed-count box above
prints `hit_probability` itself as a percentage for the existing `Hit Probability:` line.

### WEZ (Weapon Employment Zone) Sweep

> Also available in the WASM terminal (ballistics.sh): `monte-carlo --wez` with the
> WEZ flags (`--target-size`, `--wind-call-error`, `--wez-start/-end/-step`, `-o`)
> and output matching the native CLI. The terminal's `monte-carlo` does not expose
> the base-wind or hold flags (`--wind-speed`, `--wind-direction`, `--wind-vertical`,
> `--cant`, `--target-distance`, `--target-radius`), which stay at their defaults;
> it does accept `--drag-model`, which the native command lacks. A sweep runs
> num-sims full solves per range step in the browser — prefer `-n 300` for
> interactive use.

`monte-carlo --wez` answers a different question than the base command above. Instead of a
single summary at one `--target-distance`, it sweeps a range of distances and reports **hit
probability on a fixed target size at each range, holding a single zero** — the classic
"point-blank range" question: *how far out can I engage this target size without holding over
or dialing elevation?*

This matters because it means ballistic **drop below your line of sight counts as a miss
source**, exactly like it would for a real shot fired with that one zero. That is different from
the base `monte-carlo --target-distance` command, whose hit probability is measured against
*that run's own* point of aim — i.e. it implicitly assumes you re-dial correctly for every
range. A WEZ sweep does not assume that; it uses the elevation you pass with `-a/--angle` (from
`ballistics zero`, typically) for every step of the sweep.

```bash
# Zero for 300 yards first
./ballistics zero -v 2700 -b 0.475 -m 168 -d 0.308 --target-distance 300
# -> Zero Angle: 0.1432°

# Sweep 200-1000 yd in 100 yd steps against an 18"x30" target (e.g. IPSC/steel silhouette),
# with a 3 mph wind-call error on top of a 1 mph physical wind-speed uncertainty
./ballistics monte-carlo -v 2700 -b 0.475 -m 168 -d 0.308 -a 0.1432 \
  --wez --target-size 18x30 \
  --wind-call-error 3 --wind-std 1

# WEZ sweep: 1000 sims/step, wind call 3.00 mph + wind std 1.00 mph (quadrature) = 3.16 mph effective
# ┌────────────┬──────────┬───────────────┬───────────┬───────────┬───────────┐
# │ Range ( yd) │  P(hit)  │ Dominant      │ Wind call │  MV SD    │ Other/grp │
# ├────────────┼──────────┼───────────────┼───────────┼───────────┼───────────┤
# │      200.0 │    61.1% │ other         │      0.0% │      0.0% │    100.0% │
# │      300.0 │    37.5% │ other         │      0.0% │      0.0% │    100.0% │
# │      400.0 │    18.7% │ other         │      0.0% │      0.0% │    100.0% │
# │      500.0 │     9.0% │ other         │      0.0% │      0.0% │    100.0% │
# │      600.0 │     3.2% │ n/a           │      0.0% │      0.0% │      0.0% │
# │      700.0 │     0.4% │ n/a           │      0.0% │      0.0% │      0.0% │
# │      800.0 │     0.0% │ n/a           │      0.0% │      0.0% │      0.0% │
# │      900.0 │     0.0% │ n/a           │      0.0% │      0.0% │      0.0% │
# │     1000.0 │     0.0% │ n/a           │      0.0% │      0.0% │      0.0% │
# └────────────┴──────────┴───────────────┴───────────┴───────────┴───────────┘
```

`Dominant` and the three share columns read `n/a` instead of a value whenever the variance
attribution isn't meaningful for that row, for any of three different reasons that the table
does not distinguish between:

- **The range is past where the (undispersed) trajectory can reach at all** — a nearly flat
  shot from a normal bore height eventually crosses the ground plane, e.g. the `600`-`1000` yd
  rows above.
- **The shared attribution kernel hit a structural refusal** on one of the shot's error sources
  for that particular range — an input combination the underlying sensitivity solve can't
  evaluate.
- **The drag configuration has no representation in the kernel's wire format at all** — a
  loaded `--drag-table`, or a `--drag-model` of `G2`, `G5`, `GI`, `GS`, or `RA4` (only `G1`,
  `G6`, `G7`, and `G8` are supported). This is the case most likely to catch you off guard:
  every row of a `--drag-table` sweep, or a sweep using one of those unsupported drag models,
  reads `n/a` on every range — including ranges the bullet clears easily. It is NOT a claim
  that your bullet fails to reach the target.

`P(hit)` is unaffected in all three cases — it always comes from the fully-dispersed Monte Carlo
run directly, never from the attribution kernel (and will simply read 0% under the first case,
since every sample there is a miss).

**Wind call vs. ballistic wind.** `--wind-call-error` is the shooter's own uncertainty in
*reading* the wind (e.g. "I think it's 8-10 mph, call it 9") — a human estimation error. It is
distinct from `--wind-std`, which models physical gust-to-gust variability in the wind itself.
Both perturb the same physical channel (the wind speed fed to the trajectory solve), so as
independent random errors they combine **in quadrature**, not by simple addition:

```
effective_wind_std = sqrt(wind_std^2 + wind_call_error^2)
```

`--wind-call-error` defaults to `0.0` (perfect wind call — only `--wind-std` and the other
dispersion sources contribute), matching the base command's existing behavior when `--wez` is
not set.

**Target size** (`--target-size WIDTHxHEIGHT` or `--target-size RADIUS`):
- `18x30` — a rectangle, 18" wide x 30" tall (cm under `--units metric`), centered on the line
  of sight.
- `12` — a single number falls back to a circular hit radius (same units, same "distance from
  point of aim" semantics as the base command's `--target-radius`), just expressed in
  target-size units (inches/cm) instead of range units (yards/meters).

**Variance attribution.** Each row's `Dominant` column and the three share columns (`Wind
call`, `MV SD`, `Other/grp`) report which uncertainty source contributes the most to the miss
variance at that range, and their approximate shares (they sum to ~100%). `Other/grp` bundles
mechanical/ammo group dispersion (`--angle-std`, the derived azimuth spread, `--bc-std`) with
the *ballistic* (non-call) share of wind uncertainty (`--wind-std`, `--wind-direction-std`).

This is computed **analytically** from a central-difference derivative of the impact point with
respect to each independent source's own input, rather than by re-running the full Monte Carlo
sample set once per bucket with that source zeroed out. Each derivative is taken at a small step
size intrinsic to that input — independent of the sigma you configured for it — and only
afterward scaled by the source's declared standard deviation (`--wind-std`, `--angle-std`, etc.)
to turn it into a displacement, so the differencing step and the uncertainty being attributed are
deliberately two different numbers. A full decomposed re-run would multiply the sweep's cost by
the number of buckets; this instead costs a handful of extra deterministic trajectory solves per
range and keeps the default sweep's runtime in the same ballpark as the base `monte-carlo`
command (a few seconds for the default 9-step sweep). This is a first-order (linear-response)
error budget, not an exact decomposition of the Monte Carlo variance.

**Sweep range**: `--wez-start`, `--wez-end` (inclusive), `--wez-step` (all in yards for
imperial, meters for metric; defaults `200`/`1000`/`100`).

**Output**: `-o summary` (default) prints the table above; `-o full` prints JSON with
unit-labeled fields (`range_m`, `p_hit`, `dominant_error_source`, `wind_call_share`,
`mv_sd_share`, `other_share`, plus the resolved `target_size`, `wind_speed_std_mps`,
`wind_call_error_mps`, and `combined_wind_speed_std_mps`); `-o statistics` prints one CSV row
per range step.

**Performance**: `--num-sims` (`-n`, default 1000) is respected per range step, same as the
base command — a WEZ sweep is `--num-sims` trajectories times the number of sweep steps, plus a
handful of cheap deterministic solves per step for attribution.

### BC Estimation

Estimate ballistic coefficient from observed data. Supports both the **G1 and G7** drag
models and two fit bases — a **drop** curve or a downrange **velocity** curve (the latter is
immune to zero / sight-height / launch-angle error). A row is printed for each drag model ×
data basis you supply.

```bash
# Legacy two-point drop input (G1 + G7 by default)
./ballistics estimate-bc \
  -v 2700 -m 168 -d 0.308 \
  --distance1 100 --drop1 0.0 \
  --distance2 200 --drop2 0.023

# n-point drop series, G7 only
./ballistics estimate-bc -v 2650 -m 77 -d 0.224 \
  --data "300,29.0;500,89.9;700,204.6" --drag-model g7

# All four variants: G1/G7 x drop/velocity
./ballistics estimate-bc -v 2650 -m 77 -d 0.224 \
  --data "300,29.0;500,89.9;700,204.6" \
  --velocity-data "300,1980;500,1560;700,1240" \
  --drag-model both
```

Options: `--data "dist,drop;..."` (yd,in / m,mm), `--velocity-data "dist,vel;..."`
(yd,fps / m,m/s), `--drag-model g1|g7|both` (default `both`), `-o table|json|csv`.

**Dope-card (zeroed) data — use `--zero-range` and match the atmosphere.** A dope card's
drops are measured below your line of sight from a rifle **zeroed** at some range (so the
drop is ~0 at the zero and grows downrange). Pass `--zero-range` so the fit matches that
frame, and give the conditions the card was made at — BC only means something relative to
air density:

```bash
./ballistics estimate-bc -v 2650 -m 77 -d 0.224 \
  --data "100,0;300,14.2;500,61.4;700,162.4;900,343.0;1100,643.0" \
  --zero-range 100 --sight-height 2.0 \
  --temperature 59 --pressure 29.92 --altitude 0 --drag-model g7
```

Without `--zero-range`, drop is treated as **bore-referenced** (flat-fire drop below the
extended bore) — correct only for a bore-drop table, not a dope card; the tool warns if your
data looks zeroed. A fit that can't determine a value from the data (too few/short-range
points, or wrong zero/atmosphere) is flagged **UNRELIABLE** rather than returning a bogus
number. Atmosphere flags: `--temperature` (°F/°C), `--pressure` (inHg/hPa), `--humidity`
(%), `--altitude` (ft/m); `--zero-range` (yd/m), `--sight-height` (in/mm).

### True Velocity Calculation

Find the effective muzzle velocity that produces a measured drop at a known range. This helps "true" your ballistic system by identifying discrepancies between chronograph readings and real-world ballistic performance.

```bash
# Basic true velocity calculation (offline)
./ballistics true-velocity \
  --measured-drop 5.1    # Measured drop in MILs
  --range 600            # Range where drop was measured (yards)
  --bc 0.27              # Ballistic coefficient
  --drag-model g7        # G7 drag model
  --mass 140             # Bullet mass (grains)
  --diameter 0.264       # Bullet diameter (inches)
  --offline              # Use local calculation

# With chronograph velocity for comparison
./ballistics true-velocity \
  --measured-drop 5.1 --range 600 \
  --bc 0.27 --drag-model g7 \
  --mass 140 --diameter 0.264 \
  --chrono-velocity 2822 \  # Compare against chrono reading
  --offline

# Chronograph reading taken downrange, not at the muzzle (MBA-1377): most
# chronographs read 10-15 ft (or 25 m) downrange, so the raw reading is
# slightly low. --chrono-distance back-solves the true muzzle velocity from
# the measured reading using the SAME --bc/--drag-model/atmosphere, via
# secant iteration on the real forward drag model (McCoy; JBM's published
# "Distance to Chronograph" calculator).
./ballistics true-velocity \
  --measured-drop 5.1 --range 600 \
  --bc 0.27 --drag-model g7 \
  --mass 140 --diameter 0.264 \
  --chrono-velocity 2822 --chrono-distance 15 \  # measured 15 ft from the muzzle
  --offline

# With BC5D tables for improved accuracy
./ballistics true-velocity \
  --measured-drop 5.1 --range 600 \
  --bc 0.27 --drag-model g7 \
  --mass 140 --diameter 0.264 \
  --bc-table-auto \        # Auto-download BC5D tables
  --offline

# Using online API (with fallback)
./ballistics true-velocity \
  --measured-drop 5.1 --range 600 \
  --bc 0.27 --drag-model g7 \
  --mass 140 --diameter 0.264 \
  --offline-fallback       # Try API, fall back to local if fails

# Metric units
./ballistics true-velocity --units metric \
  --measured-drop 5.1 --range 549 \  # 549 meters ≈ 600 yards
  --bc 0.27 --drag-model g7 \
  --mass 9.07 --diameter 6.71 \      # grams, mm
  --offline
```

#### True Velocity Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| --measured-drop | Measured drop at `--range` (MIL by default; follows `--drop-unit` in multi-observation mode) | Required |
| --range | Range where drop was measured | Required |
| --observed | Additional observed impact `RANGE:DROP` (repeatable) — enables joint MV+BC calibration | None |
| --drop-unit | Drop unit for `--measured-drop`/`--observed` in multi-observation mode (`mil`/`moa`/`in`) | mil |
| --bc | Ballistic coefficient (starting value; fitted when observations allow) | Required |
| --drag-model | Drag model (G1/G7) | g1 |
| --mass | Bullet mass | Required |
| --diameter | Bullet diameter | Required |
| --chrono-velocity | Chronograph velocity for comparison | None |
| --chrono-distance | Distance `--chrono-velocity` was measured at, downrange of the muzzle (ft/m). Nonzero back-solves the true muzzle velocity via the real forward drag model (secant iteration); zero/absent is an exact no-op. Valid range 1-98 ft / 0.3-30 m (100 ft is out of range); requires `--chrono-velocity` | None (no-op) |
| --zero-range | Zero range | 100 yd/m |
| --sight-height | Sight height above bore | 2.0 in/50mm |
| --bullet-length | Bullet length (for BC5D lookup) | Auto-calculated |
| --offline | Force offline mode (local calculation) | false |
| --offline-fallback | Fall back to local if API fails | false |
| --bc-table-dir | Directory with BC5D tables | None |
| --bc-table-auto | Auto-download BC5D tables | false |

#### Output

The command outputs:
- **Effective Velocity**: The calculated muzzle velocity that produces the measured drop
- **Velocity Adjustment**: Difference from chrono velocity (if provided)
- **Adjustment Percent**: Percentage adjustment from chrono
- **Confidence**: High/Medium/Low based on convergence quality
- **Iterations**: Number of iterations to converge
- **Final Error**: Remaining error in MILs

Example output:
```
True Velocity Results
═════════════════════
Effective Velocity:    2740 fps
Chrono Velocity:       2822 fps
Velocity Adjustment:   -82 fps (-2.91%)
Confidence:            high
Iterations:            12
Final Error:           0.001 MIL
Calculated Drop:       5.10 MIL
```

### Wind-Call Truing (`true-wind`) — MBA-1392

`true-velocity` trues the VERTICAL axis. `true-wind` trues the other one: given where a
group actually landed left or right of the aim point, it back-solves the constant
crosswind that reproduces that miss through the real forward model — the number you
compare against the wind you *called*, so wind calls get calibrated instead of guessed at.

A horizontal miss is not purely wind, so the command insists on being able to separate it:

- **`--twist-rate` is required.** Gyroscopic spin drift is a lateral effect of the same
  order as a light wind at long range (a 1:11" .308 drifts ~3.5 in right at 700 yd), so
  without a twist rate the fit would report spin drift as wind. It is modelled and
  therefore subtracted, always.
- **`--latitude` + `--shot-direction` are optional but travel together.** Supply both and
  Coriolis is modelled and subtracted too; supply neither and it stays absorbed in the
  solved wind, which the report says out loud. Supplying only one is a hard error — a
  half-specified Earth frame is a mistake, not a partially usable input.

Anything the model has no data for is listed under *NOT subtracted*, so a contaminated
number is never presented as pure wind.

#### Signs (all of them, in one place)

| Quantity | Positive means |
|----------|----------------|
| `--miss` RIGHT value | the group landed **right** of the aim point |
| solved crosswind | wind **FROM the shooter's LEFT** (9 o'clock), which pushes impacts **right** |
| `--called-wind` | same convention as the solved value |
| `--twist-right` | right-hand twist (drifts right); `--twist-right false` = left-hand |
| `--shot-direction` | compass bearing fired *along*, 0 = North, 90 = East |

The solved wind's sign follows the deflection it produces, and it is a **full-value**
crosswind (a 90° wind), not a half-value component of some other bearing. Internally that
is the engine's wind-FROM convention (`0` = headwind, `π/2` = from the right, `3π/2` = from
the left) — the convention established by the **0.19.0 wind-direction sign fix**, which
flipped `0` from tailwind to headwind.

#### Units

`--miss RANGE:RIGHT[:SIGMA]`: **RANGE follows `--units`** (yards imperial, meters metric),
while **RIGHT and SIGMA are LINEAR INCHES in both unit systems** — a tape measurement off
the target, exactly like `true-velocity --drop-unit in`. Wind speeds (`--called-wind` and
the solved value) follow `--units` like `--wind-speed` does: mph imperial, m/s metric.

Because `--miss` is a linear measurement and not a dial reading, the MBA-1358 scope
tracking correction factors do **not** apply to it — only *dialed* observations are
CF-converted. This command deliberately exposes no `--elevation-cf`/`--windage-cf`.

#### Worked example

```bash
# Two groups, both right of aim: 14.0 in at 500 yd and 29.5 in at 700 yd.
# You called a 9 mph left-to-right wind.
./ballistics true-wind \
  --miss 500:14.0 --miss 700:29.5 \
  -v 2700 -b 0.475 -m 168 -d 0.308 --drag-model g7 \
  --twist-rate 11 \
  --called-wind 9
```

```
=== EFFECTIVE WIND TRUING (from observed horizontal miss) ===

  Range (yd)     Miss (in)   Spin/Cor (in)        Wind (mph)  Resid (in)
  ----------------------------------------------------------------------
       500.0        +14.00           +1.75            +12.68      -0.000
       700.0        +29.50           +3.49            +13.23      +0.000
  ----------------------------------------------------------------------

  Effective crosswind:   +12.96 mph  (mean of 2 observations)
  Called wind:            +9.00 mph   ->  wind-call correction factor 1.44
    (multiply your wind calls by 1.44 to match what actually hit)
  Effects subtracted:  spin drift
  NOT subtracted (absorbed into the solved wind): Coriolis (supply --latitude and --shot-direction to subtract it)
```

Read that as: the wind was really ~13 mph, you called 9, so you habitually under-call by a
factor of ~1.44. The *Spin/Cor* column is the part of each miss the model attributed to
spin drift (and Coriolis, when an Earth frame is supplied) rather than to wind.

```bash
# Subtract Coriolis too, and weight the observation by its measurement error.
./ballistics true-wind \
  --miss 600:16.5:0.75 \
  -v 2700 -b 0.475 -m 168 -d 0.308 --drag-model g7 --twist-rate 11 \
  --latitude 45 --shot-direction 90 \
  -o json
```

With a sigma on every `--miss`, the combined answer is inverse-variance weighted and
carries its own 1σ; with a sigma on none, it is the plain mean. Mixing the two is a hard
error rather than a silent blend of weighting schemes.

#### `true-wind` Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--miss` | Observed horizontal miss `RANGE:RIGHT[:SIGMA]`, repeatable. RANGE follows `--units`; RIGHT/SIGMA are inches in both systems; RIGHT is signed (positive = right of aim) | Required |
| `-v`, `--velocity` | **Known** muzzle velocity (true it up with `true-velocity` first if unsure) | Required |
| `-b`, `--bc` | Ballistic coefficient | Required |
| `--drag-model` | Drag model (G1/G7) | g1 |
| `-m`, `--mass` | Bullet mass | Required |
| `-d`, `--diameter` | Bullet diameter | Required |
| `--twist-rate` | Barrel twist (in/turn imperial, mm/turn metric) | **Required** |
| `--twist-right` | Right-hand twist (`--twist-right false` = left-hand) | true |
| `--latitude` | Firing latitude, degrees. Needs `--shot-direction` | None |
| `--shot-direction` | Shot azimuth, degrees (0 = N, 90 = E). Needs `--latitude` | None |
| `--called-wind` | The wind you called (mph/m·s⁻¹) — adds the wind-call correction factor | None |
| `--zero-distance` | Zero distance | 100 yd/m |
| `--sight-height` | Sight height above bore | 2 in / 50 mm |
| `--temperature`, `--pressure`, `--pressure-type`, `--humidity`, `--altitude` | Atmosphere (same meanings as elsewhere) | 59 °F / 29.92 inHg / 50% / 0 |
| `--units` | Unit system | imperial |
| `-o`, `--output` | `table`, `json`, or `csv` | table |

#### Reported quantities

- **Effective crosswind** — signed, mph or m/s. Per observation and combined.
- **Spin/Cor** — the model's zero-wind lateral: what was attributed to spin drift, plus
  Coriolis when an Earth frame was supplied.
- **Wind-call correction factor** — `solved ÷ called`. Greater than 1 means you under-call
  the wind; less than 1 that you over-call it; **negative means you called the wrong side**.
- **Validity note** — an observation whose impact moves less than 0.25 in per mph of
  crosswind is flagged as weakly identified (the same idea as `true-velocity`'s
  MV-calibration window: short-range groups barely constrain the wind).
- An observed miss no crosswind within ±100 mph can produce is **rejected** with a
  diagnostic naming the solvable band, rather than clamped into a plausible-looking answer.

The command is entirely offline — there is no API for it and no flag that could reach one.
It is available on the WASM terminal with identical output (both surfaces render through
one shared formatter).

### Reticle Hold Points (`reticle`) — MBA-1361

Most of this tool answers "what do I dial?". `reticle` answers the other question: *where
in my reticle do I put the target?* It takes an angular firing solution and reports the
hold point plus the nearest mark, FFP/SFP aware.

```bash
# Build a description (this one is pure geometry, no solve involved)
ballistics reticle generate tree --rows 6 --row-spacing 1.0 --spread-step 0.5 \
  --name "MyScope MIL" --focal-plane sfp --reference-mag 12 -o json > myscope.json

# Ask where a solution lands in it, at the magnification you're actually on
ballistics reticle hold --reticle-json myscope.json --mag 6 --drop-mil 4.2 --wind-mil 1.1

# Or let it solve the trajectory for you
ballistics reticle hold --reticle-json myscope.json --mag 6 --range 600 \
  -v 2700 -b 0.475 -m 168 -d 0.308 --zero-distance 100 --wind-speed 10 --wind-direction 90
```

#### Angular conventions (all of them, in one place)

Every angle is a **milliradian**, measured from the optical center:

| Quantity | Positive means |
|----------|----------------|
| `down_mil` / `--drop-mil` | **below** center — a holdover point |
| `right_mil` / `--wind-mil` | to the shooter's **right** of center |

The mapping from a firing solution to a hold is deliberately direct: if the bullet falls
`d` mil below the line of sight you put a point `d` mil *below* center on the target, and
if the wind pushes it `w` mil right you put a point `w` mil *right* of center on the
target. The command's real work is the mark search, not the coordinates.

#### Focal plane

Published optics-manual math, nothing more:

- **FFP** — the reticle is magnified with the image, so a mark subtends the same angle at
  every magnification. `--mag` is validated but changes nothing.
- **SFP** — the reticle is a fixed size at the eyepiece while the image scales, so a
  mark's **true** subtension is `nominal × reference-mag ÷ magnification`. A 2 mil mark on
  a reticle calibrated at 10x covers **4 mil** of target at 5x, and 1 mil at 20x.

The hold point is a property of the trajectory, so it is always **true angular**; only the
MARKS are rescaled, and the nearest-mark distance is measured after that rescaling. A
non-physical `--mag` (zero or negative) is rejected on both planes, and an SFP description
with no usable `reference_magnification` is rejected too.

`off_reticle` fires when the hold falls outside the marks' bounding box grown by **20 % of
that box's span, per axis**. An axis on which every mark shares a coordinate (a pure BDC
ladder has no windage marks) has zero margin, so any deviation there reads as off-reticle
— which is honest: such a reticle offers nothing to hold on in that direction.

#### Generators

| Layout | Shape |
|--------|-------|
| `mil-grid --spacing --extent` | a plain mil-hash **cross** — marks along both stadia, *not* a filled 2-D grid |
| `tree --rows --row-spacing --spread-step` | a generic widening holdover tree: row *N* carries windage dots at ±1..*N* steps |
| `bdc --drop RANGE:DROP_MIL ...` | a labeled ladder from **already-solved** drops |

`bdc` deliberately runs no solve of its own — produce the drops with `come-ups` or
`trajectory` first, so the ladder's provenance (load, atmosphere, zero) stays yours.
All three accept `--name`, `--focal-plane`, `--reference-mag` and `-o table|json`.
`-o json` emits exactly what `reticle hold --reticle-json` and
`profile save --reticle-json` consume.

#### Intellectual-property exclusions

Horus grid reticles and Time-of-Flight Wind Dots are actively patented, and Horus
monetizes app integration through its own licensed app. So, by design and permanently
unless licensed: **no TREMOR-family or Horus grid layouts**, **no wind-dot calibration**
(wind enters only as an angular deflection you already solved), and **no vendor reticle
catalog**. Manufacturer subtension sheets are published facts and a legally viable catalog
source, but curating one is a separate per-vendor IP-reviewed data project.

#### Saved profiles

```bash
ballistics profile save myrifle -v 2700 -b 0.475 -m 168 -d 0.308 --reticle-json myscope.json
ballistics reticle hold --profile myrifle --mag 6 --drop-mil 4.2
```

The attachment is carried forward by an unrelated re-save (like the DSF table) and dropped
by `--clear-reticle`. Mark coordinates are angular, so `--units` conversion leaves them
alone — for the same reason turret click graduations are left alone.

#### Other surfaces

- **solve-json v1** gains an optional request block
  `reticle: {range_m, magnification, description}` and, only when it is present, a
  `reticle_hold` object in the response. Omitting it leaves every existing response
  byte-identical. The hold is read from the response's own samples, so the two can never
  describe different trajectories; a range outside the sampled trajectory is a structured
  error, never an extrapolation.
- **WASM terminal**: `reticle generate` (all three layouts) and `reticle hold`, rendering
  through the *same* formatter as native, so identical inputs produce identical bytes.
  Two deliberate differences: `--reticle-json` there takes the description as **inline
  JSON text** (there is no filesystem in the browser), and `--range`/`--profile` are
  native-only — run `trajectory` first and pass the drop as `--drop-mil`.
- **FFI**: a new appended `FFIReticleHold` struct and
  `ballistics_hold_point_in_reticle(...)`. Nothing existing moved, so no consumer needs a
  recompile. `marks_len` is bounds-checked against `MAX_FFI_RETICLE_MARKS` before a single
  element is read.

### Robust Hold Corridors (`hold-corridor`) — MBA-1349

You usually have a handful of concrete plausible wind calls — a low call and a high call,
or two different downrange patterns — not a distribution you actually believe. A nominal
trajectory hides that ambiguity, and Monte Carlo asks for a probability model you do not
have. `hold-corridor` takes the middle path: it solves each *named* scenario exactly and
reports all of them, plus the corridor they span.

```bash
cat > calls.json <<'JSON'
{"version": 1, "nominal": "low",
 "scenarios": [
   {"name": "low",    "segments": ["4:90:1000"]},
   {"name": "high",   "segments": ["14:90:1000"]},
   {"name": "switch", "segments": ["10:90:400", "8:270:1000"]}]}
JSON

ballistics hold-corridor --scenarios calls.json --ranges 200,400,600 \
  --target rect:12x18 -v 2700 -b 0.475 -m 168 -d 0.308 --zero-distance 100
```

Segments use exactly the `SPEED:ANGLE:UNTIL[:VERTICAL]` tokens `--wind-segment` takes, and
`--units` applies to them identically (mph/yards imperial, m/s and meters metric). One
grammar, one parser.

#### What it reports, per range

- **Every scenario's hold**, in milliradians: elevation positive = hold UP, windage
  positive = hold RIGHT.
- **The corridor** — the min and max over the scenarios, on each axis.
- **The minimax hold** — the Chebyshev center of the scenario holds under the target's
  metric, and the **worst-case miss** from it, naming the scenario that produces it.
- **What the nominal would have cost**, when the set names one: the same worst-case
  measure computed from holding that scenario instead.
- **`fits_target`** — whether that one hold keeps every scenario on the target.
  **Boundary contact counts as a fit.**

#### The two metrics (both documented, because they differ)

| `--target` | Metric | Minimax hold | Worst case |
|---|---|---|---|
| `rect:WIDTHxHEIGHT` or omitted | per-axis (L∞) | the **midpoint of the extremes on each axis independently** | the larger per-axis half-span |
| `circle:DIAMETER` | Euclidean (L2) | the center of the **minimum enclosing circle** of the scenario holds | that circle's radius |

The rectangle's axes are independent, so the per-axis midpoint is optimal for each and
therefore for the rectangle. For the circle the objective is genuinely two-dimensional, so
the answer is the 1-center of the point set — computed exhaustively over the candidate
circles (one point, a pair's diameter, a triple's circumcircle), which with ≤8 scenarios is
under a hundred candidates and, unlike the usual randomized construction, **deterministic**.
With no `--target` there is no shape to judge against, so the per-axis metric is used and
no fit is reported. Target dimensions are inches (imperial) or centimeters (metric).

#### Guarantees

- **The zero is solved once, in calm air**, and reused by every scenario. A rifle has one
  zero; re-zeroing per scenario would fold each hypothesis into its own datum and collapse
  the very corridor this exists to show.
- **Caps are enforced before any solving**: ≤ 8 scenarios and ≤ 64 ranges, with structured
  errors. A malformed segment, an unknown `nominal`, a duplicate scenario name, a duplicate
  range and an unsupported `version` are all rejected the same way, before a single
  trajectory runs.
- **Reordering scenarios changes nothing.** They are sorted by name internally before
  anything is solved, so the corridor, the minimax hold and the named worst-case scenario
  are all independent of the order in the file.

#### Non-goals (deliberate, permanent)

No probability is assigned to any scenario, no continuous guarantee is claimed between
them, and the finite set is never folded into a standard deviation. **A three-scenario
corridor is not a confidence interval**, and the output says so. Statistical dispersion
remains `monte-carlo --wez`'s job.

`-o json` emits the versioned `RobustHoldReportV1` (SI throughout). Native-only this train;
the core and its formatter are already shared-ready for the WASM follow-up.

### Reticle/BDC Inverse Solvers — MBA-1362

Three read-only solvers over an existing load. All three share **one** drop-vs-range
helper — a single solved, finely sampled angular-drop curve with a forward lookup and a
monotone bisection inverse — so none of them can disagree with the others (or with
`reticle hold --range`) about what "the drop at 500 yards" is. They are CLI-only this
train; WASM parity is a tracked follow-up.

#### `mark-to-range` — what range is this mark good for?

The inverse of a come-up table (Nightforce / Nikon Spot On / Swarovski / TRACT).

```bash
ballistics mark-to-range --mark 2 --mark 4 --mark 6 --mark 8 \
  --focal-plane sfp --reference-mag 12 --mag 8 \
  -v 2700 -b 0.475 -m 168 -d 0.308 --zero-distance 100
```

Each subtension is converted to TRUE angular for the focal plane and magnification in use
(FFP: unchanged; SFP: `nominal × reference-mag ÷ mag`), then the range where angular drop
equals it is found by bisection past the far zero — the only interval on which drop-vs-range
is monotone. Each mark reports its range, time of flight, remaining velocity and energy.

**Marks that do not map to a range are reported, never dropped**, so the table always has
one row per mark and the reticle's usable limits are visible:

- `not_a_holdover` — a subtension at or above the optical center. Angular drop is exactly
  zero at the far zero and only grows past it, so these correspond to ranges *inside* the
  zero, where the relationship is not single-valued.
- `beyond_max_range` — the load's drop never grows that far within `--max-range`, and the
  row says what the deepest reachable hold actually is.

`--reticle-json <FILE>` takes the subtensions, focal plane and reference magnification from
a description instead of `--mark` flags. Marks on the vertical stadium are all used
(including any above center, which come back as `not_a_holdover`); the optical center and
off-axis windage marks are excluded, and their count is printed.

#### `bdc-match` — what magnification makes this BDC reticle fit?

The Zeiss Rapid-Z question, and **second focal plane only** — an FFP reticle's subtensions
do not depend on magnification, so no magnification fits better than another, and the
command says so rather than returning a meaningless number.

```bash
ballistics bdc-match --mark-range 2:300 --mark-range 4:500 --mark-range 6:700 \
  --reference-mag 12 -v 2700 -b 0.475 -m 168 -d 0.308 --zero-distance 100
```

Each `--mark-range MIL:RANGE` is a mark's etched subtension and the range it is meant to
hit; at least two are required. The apparent subtension of mark *i* at magnification *M* is
`s_i × reference_mag ÷ M`, so substituting `u = reference_mag ÷ M` makes the residual
linear in `u` and the best fit is the ordinary least-squares slope through the origin,
`u* = Σ(s_i·t_i) ÷ Σ(s_i²)`. It is closed form — exact, deterministic, and independent of
any starting guess. The report gives the fitted magnification, the per-mark residuals, and
a **warning above `--residual-warn`** (default 0.2 mil) saying in plain words that the
reticle does not fit the load at any magnification and the answer is the least-bad
compromise, not a calibration.

#### `optimal-zero` — one zero for a whole stage

GeoBallistics HDZ class. Unlike `mpbr`, which sizes one vital zone around one zero, this
balances a specific list of ranges against each other.

```bash
ballistics optimal-zero --target 200:10 --target 400:12 --target 600:18 \
  -v 2700 -b 0.475 -m 168 -d 0.308
```

`--target RANGE[:HEIGHT]` is repeatable, 2 to 16 entries; HEIGHT is the target's full
vertical size (inches imperial / cm metric) and falls back to `--vital` when omitted. The
search minimizes, over zero distance, the LARGEST hold any target needs — a golden-section
search over a bracketed zero range with fixed constants, so the answer is reproducible
bit for bit. The report gives the optimal zero, every target's hold, the signed miss a
dead-center hold would produce (positive = the bullet lands low), and whether each target
is inside its vital zone.

Targets are sorted by range internally, so the answer does not depend on the order they
were typed in. A correct min-max solution equalizes the extremes: the nearest and furthest
targets end up needing the same magnitude of hold in opposite directions.

### DSF (Drop-Scale-Factor) Truing

`dsf` is the second stage of Applied Ballistics' two-stage truing workflow (MBA-1357).
Once `true-velocity` has fixed the supersonic muzzle velocity/BC against a chronograph
and a supersonic (Mach > 1.2) observed drop, drop discrepancies that grow through the
transonic region and into the subsonic regime are no longer fixable by a single MV
correction — the residual is a slowly-varying function of Mach. `dsf` records a handful
of *observed drop / predicted drop* ratios at specific (Mach <= 1.2) ranges and stores
them, keyed by Mach, on a saved profile; `trajectory --saved-profile` and
`come-ups --profile` then scale predicted drop by the nearest table entry automatically.

Unlike every other command, `dsf` takes **no ballistic parameters of its own** — it
solves the named saved profile's own trajectory (same physics `trajectory
--saved-profile NAME` would fly with no other flags) and derives everything else from
`--range` and `--observed-drop`.

```bash
# Record an observed 5.1 mil drop at 900 yards on a profile already MV-trued.
./ballistics dsf --saved-profile my-rifle --range 900 --observed-drop 5.1mil

# Same, drop given in MOA or linear inches instead.
./ballistics dsf --saved-profile my-rifle --range 900 --observed-drop 17.4moa
./ballistics dsf --saved-profile my-rifle --range 1000 --observed-drop 42.0in
```

#### DSF Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--saved-profile` | Saved profile to solve, calibrate, and write the point back to | Required |
| `--range` | Range at which the drop was observed (yards imperial, meters metric) | Required |
| `--observed-drop` | Observed drop, value and unit with NO separator: `mil`, `moa`, or `in` (e.g. `5.1mil`, `17.4moa`, `42.0in`) | Required |

#### How a point is derived

1. Solves the saved profile's trajectory and finds the predicted drop and Mach number at
   `--range` (Mach = velocity ÷ the solver's station speed of sound — the same divisor
   `apply_dsf` uses later, so a point derived here lands back on the identical Mach when
   applied).
2. `dsf = observed / predicted`, expressed in the same unit `--observed-drop` used.
3. **Staging gates:**
   - An observation whose target-range Mach is **> 1.2** is rejected outright — that's
     `true-velocity`'s territory, not DSF's:
     ```
     error: observation is supersonic (Mach 1.35); calibrate muzzle velocity first (true-velocity), then collect DSF points at Mach <= 1.2
     ```
   - A warning (not fatal) when the observation range is beyond 90% of the trajectory's
     solved max range — the solution's reliability degrades past this point:
     ```
     warning: observation at 950 yd is beyond 90% of the solved range; solution reliability degrades past this point
     ```
   - A warning (not fatal) when the table's highest-Mach point still sits below Mach 0.9
     — nothing in the table covers the 0.9-1.2 transonic band:
     ```
     warning: no DSF point in the transonic band (Mach 1.2-0.9); transonic drops remain uncorrected
     ```
   - A note (stdout, always) naming the DSF *validity window* — at or beyond 90% of the
     trajectory's downward Mach 0.9 crossing — using the same solve `dsf` already
     performed above:
     ```
     note: DSF window: at or beyond 620.4 yd (90% of the Mach 0.9 distance)
     ```
     or, if the trajectory never goes subsonic within the solved range:
     ```
     note: no subsonic window inside the solved range
     ```
4. The point is added to the profile's table (up to 6 distinct points): a new point
   within 0.05 Mach of an existing one **supersedes** it (reported on stdout); a 7th
   distinct point is rejected, naming the 6-point cap and `--clear-dsf` to make room.
5. The profile is saved and the resulting table is printed.

#### Auto-apply

`trajectory --saved-profile NAME` and `come-ups --profile NAME` automatically apply a
profile's DSF table to every solved point's drop (byte-identical velocity/energy/time of
flight/windage — this is a **drop-only** correction) and, for **table output only**,
print:
```
DSF table active (2 points, Mach 0.65-0.95)
```
JSON and CSV output carry the corrected drop numbers too, but get no equivalent text or
extra top-level field — the note above is purely a human-facing display detail.

> The WASM terminal (ballistics.sh) supports DSF truing on `trajectory` too, but has no
> saved-profile storage to carry a table between calls — pass it per call instead with
> one or more repeatable `--dsf-point MACH:DSF` flags (e.g. `--dsf-point 0.65:1.082
> --dsf-point 0.95:1.031`), up to 6. Validation and the auto-apply/note behavior are
> identical to the native CLI's saved-profile path above (MBA-1411). One difference from the
> native `dsf` verb: per-call points are NOT merged — the native verb supersedes points within
> 0.05 Mach of each other, but `--dsf-point` passes your list through as-is, so keep the Mach
> keys distinct.

### MCP Server (`mcp`)

`ballistics mcp` runs a [Model Context Protocol](https://modelcontextprotocol.io/) server over
the stdio transport, so an MCP-aware AI assistant can drive the engine directly instead of
shelling out to the CLI. It speaks newline-delimited JSON-RPC 2.0: one JSON-RPC message per
line on standard input, one JSON-RPC message per line on standard output, exactly as the MCP
stdio transport specifies. Nothing else is ever written to stdout.

```bash
# Run directly (an MCP client normally launches this for you as a subprocess)
./ballistics mcp
```

#### Claude Desktop configuration

Add an entry to Claude Desktop's `claude_desktop_config.json` pointing at the built binary:

```json
{
  "mcpServers": {
    "ballistics-engine": {
      "command": "/absolute/path/to/ballistics",
      "args": ["mcp"]
    }
  }
}
```

Restart Claude Desktop after editing the config. The `command` must be an absolute path; the
binary is discovered from a build (`cargo build --release`, binary at
`target/release/ballistics`) or from an installed release.

#### Implemented methods

| Method | Behavior |
| --- | --- |
| `initialize` | Returns `serverInfo` (`name: "ballistics-engine"`, this crate's version), `capabilities: {"tools": {}}`, and echoes back whatever `protocolVersion` the client requested (falling back to a fixed recent version when the client omits it). |
| `notifications/initialized` | Accepted as a no-op, per the MCP lifecycle. |
| `tools/list` | Returns the two tools below with their JSON Schema `inputSchema`. |
| `tools/call` | Invokes `solve` or `engine_info`; see below. |
| `ping` | Returns an empty result. |

Any other method returns JSON-RPC error `-32601` (Method not found). Malformed JSON returns
`-32700` (Parse error); a structurally invalid JSON-RPC message (not an object, wrong
`jsonrpc` version, missing `method`) returns `-32600` (Invalid Request). A single JSON-RPC
message is capped at 1 MiB; an oversized line is rejected with `-32700` without buffering it in
full, and the session keeps running. Malformed input never terminates the server — only closing
stdin (EOF) does, which exits `0`.

#### Tools

**`solve`** — arguments *are* a [solve-json v1](docs/SOLVE_JSON_V1.md) request object
(`schema_version`, `projectile`, `rifle`, `shot`, `atmosphere`, `wind`, `solver`, `effects`,
`sampling`; explicit SI units throughout). The tool result's text content is the solve-json v1
response JSON — either a success envelope with `resolved_request`/`summary`/`samples`, or a
solve-json v1 error envelope. Arguments that are not even a structurally valid solve-json v1
request (unknown or missing fields, wrong JSON types, an unsupported `schema_version`) are
instead rejected as a JSON-RPC `-32602` (Invalid params) protocol error, with the solve-json v1
error object attached as `error.data`. A well-formed request the engine cannot solve (an
out-of-range value, a resource limit, a genuine solve failure) is reported as a normal
`tools/call` result with `isError: true` instead — see the doc comment in
[`src/mcp_command.rs`](src/mcp_command.rs) for the full rationale behind this split.

**`engine_info`** — no arguments. Returns this crate's version, the drag models solve-json v1
accepts (`G1`, `G6`, `G7`, `G8`), and the crate feature flags this binary was compiled with.

Only these two tools are exposed in this pass; other CLI subcommands (Monte Carlo, BC
estimation, true velocity, profile import, and so on) are not wrapped as MCP tools yet.
#### Joint MV + BC Calibration (multiple observed impacts)

A single observed drop cannot separate muzzle velocity from ballistic
coefficient — mid-range (fully supersonic) drops are dominated by time of flight,
which is set by muzzle velocity, while BC only bites once the bullet bleeds into
the transonic region. Supply **two or more** observed impacts with the repeatable
`--observed RANGE:DROP` flag and `true-velocity` fits **both** muzzle velocity and
BC jointly against the real forward model (a full trajectory solve per candidate),
using damped Gauss-Newton (Levenberg-Marquardt).

The primary `--range`/`--measured-drop` pair is the first observation; each
`--observed` adds another. `RANGE` follows `--units` (yd imperial / m metric) and
`DROP` follows `--drop-unit` (`mil` default, or `moa` / `in`). This mode is always
computed locally.

```bash
# .308 168gr — three drops spanning supersonic -> near-transonic
./ballistics true-velocity \
  --range 300 --measured-drop 1.30 \   # first observation (yd : mil)
  --observed 600:4.40 \                # second observation
  --observed 900:9.00 \                # third observation
  --bc 0.45 --drag-model g1 \          # starting BC (fitted from here)
  --mass 168 --diameter 0.308

# MOA drops instead of MILs
./ballistics true-velocity \
  --range 300 --measured-drop 4.47 \
  --observed 600:15.1 --observed 900:30.9 \
  --drop-unit moa \
  --bc 0.45 --mass 168 --diameter 0.308
```

Example output:
```
=== VELOCITY + BC TRUING (multi-observation) ===

  Fitted muzzle velocity:    2676.2 fps
  Fitted BC:                 0.4813  (input 0.4500)

  Range (yd)  Observed (mil)  Predicted (mil)   Resid (mil)
  --------------------------------------------------------
       300.0           1.300           1.274        -0.026
       600.0           4.400           4.415        +0.015
       900.0           9.000           8.997        -0.003
  --------------------------------------------------------
  RMS residual: 0.018 mil   |   iterations: 4

  Joint MV+BC fit, excellent: RMS residual 0.018 mil, conditioning 148
  Diagnostics: BC sensitivity ratio 0.3013, conditioning 148
  MV-calibration window: 656.7-729.7 yd (90-100% of the Mach 1.2 distance)
  for optimal observation ranges run: ballistics plan-truing
```

**MV-calibration window.** The finally fitted load is re-solved (independent of
the observation set) to find where it crosses downward through Mach 1.2 — the
90-100% span of that distance is the range band where a drop residual most
cleanly identifies muzzle velocity. Table output only. If the trajectory never
crosses downward through Mach 1.2 within a generous fixed envelope, a note
prints instead, and its text depends on *why* there is no crossing. A load
still supersonic at the end of the envelope prints:
```
note: no MV window: trajectory is supersonic through 3109.4 yd; MV is identifiable at any range
```
A load that launches below Mach 1.2 (e.g. a subsonic/suppressed build) and so
never crosses downward at all prints a different note instead:
```
note: no MV window: trajectory never reaches Mach 1.2; calibrate muzzle velocity with a chronograph, then collect DSF points
```
Any observation outside the window gets a per-observation warning on stderr
(regardless of `-o`):
```
warning: observation at 300.0 is outside the MV-calibration window (656.7-729.7); MV fits from this range are weakly identified
```
See [`plan-truing`](#design-an-identifiable-truing-experiment-plan-truing) below
for choosing observation ranges up front instead of diagnosing them after the fact.

**Identifiability / honest refusal.** Before fitting BC, the command measures how
strongly the observation set constrains it: a *BC sensitivity ratio* (relative
influence of a fractional BC change vs a fractional MV change on the predicted
drops) and a *condition number* (collinearity of the MV and BC effects). If the
observations are all short and closely spaced — so a BC change is indistinguishable
from an MV change — the command **refuses the joint fit**, fits muzzle velocity
only, holds BC at the input value, and prints the reason. It never reports a
BC pulled out of ill-conditioned data as if it were precise. To constrain BC, add
a longer-range / transonic observation.

```
=== VELOCITY + BC TRUING (multi-observation) ===

  Fitted muzzle velocity:    2675.6 fps
  BC:                        0.4500  (held; not fitted)
  ...
  MV-only fit, excellent: RMS residual 0.005 mil (BC held fixed)
  Note: observations do not constrain BC (BC sensitivity ratio 0.09 < 0.20
        threshold); BC held at input 0.450. Add a longer-range / transonic
        observation to fit BC.
```

`-o json` reports the fitted values, per-observation residuals, RMS, the
identifiability diagnostics, and a self-describing `legend` with unit-labelled
field names (`range_yd`, `observed_drop_mil`, `predicted_drop_mil`,
`residual_mil`, `rms_residual_mil`, `fitted_muzzle_velocity` + `velocity_unit`).
It also carries the MV-calibration window as two additive fields,
`mv_window_start_m` / `mv_window_end_m` (meters, `null` when there is no
window) — JSON/CSV never get the note text above, only these numbers.

> With zero `--observed` flags the command behaves exactly as the classic
> single-observation velocity truing described above.

> The WASM terminal (ballistics.sh) supports the `true-velocity` command —
> single- and multi-observation — with output matching the native CLI. The
> gaps: the online/BC5D flags (`--bc-table-dir`, `--bc-table-auto`,
> `--bc-table-url`, `--offline-fallback`, `--api-url`, `--api-timeout`) and
> `--bullet-length` are not exposed there (`--offline` is accepted as a no-op;
> the terminal always calculates locally).

#### Design an identifiable truing experiment (`plan-truing`)

`plan-truing` answers the question *before* impacts are collected: which of the
target distances actually available at this facility best separate muzzle
velocity from a scalar BC? It uses the same trajectory solver and central
finite-difference perturbations as `true-velocity`, and it never fabricates an
observation or changes a saved profile.

```bash
# Explicit discrete target stations
./ballistics plan-truing \
  -v 2700 -b 0.475 --drag-model g1 -m 168 -d 0.308 \
  --candidate-ranges 200,300,400,500,600,700,800,900 \
  --observation-count 3 \
  --minimum-separation 100 \
  --measurement-resolution 0.03 \
  --drop-unit mil

# Or explicitly discretize an interval; this expands to 200,300,...,1000
./ballistics plan-truing \
  --profile 308-match \
  --range-grid 200:1000:100 \
  --observation-count 3 \
  --minimum-separation 150 \
  --measurement-resolution 0.05 \
  --output json
```

Ranges and minimum separation follow the global `--units`. The candidate list
and `--range-grid START:END:STEP` conflict: the interval is never silently
discretized at an engine-chosen step. Every selected range is one of the
declared candidates, the requested count is exact, and all selected pairs honor
the minimum separation. For small candidate sets the optimizer exhaustively
checks every feasible combination; large sets use deterministic greedy
construction plus one-for-one exchange. Input order does not change the answer.

`--measurement-resolution` means the independent **1σ standard deviation of one
impact reading** in `--drop-unit` (`mil`, `moa`, or `in`). The report repeats this
assumption and shows how the weak-axis fractional uncertainty scales at half and
twice that sigma. The discrete design is also re-optimized at those two
resolutions; if either assumption changes the selected station set, the report
warns that the recommendation is resolution-sensitive instead of implying it is
robust. It also reports:

- the chosen stations and predicted nominal drop;
- fractional MV and BC sensitivity at each station;
- each selected station's leave-one-out information gain;
- BC sensitivity ratio, condition number, singular values, log determinant, and
  weak-axis fractional 1σ uncertainty;
- eligible-but-unselected and rejected/unreachable candidates with reasons;
- an explicit `mv_only` recommendation when the available ranges do not separate
  MV from BC.

The finite information-gain score is `0.5 log det(I + F)` in fractional MV/BC
coordinates. Here `I` is a disclosed identity reference information matrix
(unit reference covariance for fractional parameter changes), used only to make
the experiment-design score finite. It is not a prior injected into subsequent
truing; the unregularized singular values, determinant, and condition number are
reported alongside it.

Saved profiles are resolved by the native CLI; explicit load/atmosphere flags
override their scalar values. V1 deliberately supports one scalar G1/G7 BC.
The nominal design point must lie inside the joint truing bounds: 1000–5000 fps
after unit conversion and scalar BC 0.05–2.0.
Profiles with velocity-banded BC values or a custom Mach/Cd curve are rejected:
varying a single BC in either model would be physically meaningless. Supporting
a fitted scale for a complete drag deck is a different parameter model.

#### Uncertainty-aware MV + BC truing

The legacy joint fit is intentionally conservative: it reports a point estimate
only when an identifiability gate passes, otherwise it holds BC fixed. For an
explicit probabilistic analysis, declare the measurement errors:

```bash
./ballistics true-velocity \
  --range 500 --measured-drop 3.179 \
  --observed 600:4.349:0.03 \
  --observed 900:8.891:0.02 \
  --observation-sigma 0.03 \
  --bc 0.45 --drag-model g1 --mass 168 --diameter 0.308 \
  --predict-range 1000 \
  --prediction-sigma 0.03 \
  --output json
```

`--observation-sigma SIGMA` supplies the known absolute 1σ error for the primary
`--range`/`--measured-drop` pair and the default for every additional impact.
An optional third field in `--observed RANGE:DROP:SIGMA` overrides it for that
reading. Drop and sigma use `--drop-unit`; range uses `--units`. All sigmas must
be positive and finite. This is weighting, not shot dispersion inferred from the
residuals: a reading with half the sigma receives four times the least-squares
weight.

Optional independent normal priors are explicit:

```text
--mv-prior MEAN:SIGMA     # fps or m/s according to --units
--bc-prior MEAN:SIGMA     # dimensionless
```

No prior is inferred from the input BC, chronograph velocity, a saved profile, or
the residual RMS. `--chrono-velocity` remains a display-only comparison — it never
feeds the joint fit, and neither does `--chrono-distance` (MBA-1377): the distance
correction only adjusts the value being *compared*, before that comparison happens.
The MAP
minimizes the weighted observation chi-square plus only the priors shown in the
request. Because observation sigmas are declared as absolute known standard
deviations, the local covariance is `(Jᵀ W J + prior precision)⁻¹`; it is **not**
multiplied by residual variance.

The scalar BC and every prior mean must lie within the constrained joint-fit
bounds (MV 1000–5000 fps after conversion; BC 0.05–2.0). These bounds apply only
to the opt-in uncertainty fit; the legacy point path retains its historical
input contract.

The v1 report contains:

- the weighted joint MAP, per-observation sigma/residual/standardized residual;
- MV and BC 95% local-Gaussian intervals, physical covariance, and correlation;
- data/prior chi-square, effective fitted-parameter count, effective residual
  degrees of freedom, and reduced chi-square when defined;
- warnings for weak BC sensitivity, collinearity, prior domination, interval/bound
  overlap, extrapolation, or an unavailable approximation;
- a latent model-drop interval at every repeatable `--predict-range`;
- when `--prediction-sigma` is supplied, a distinct, wider future-observation
  interval that adds that declared reading error in quadrature.

The approximation is Laplace/Gauss-Newton around a constrained MAP, not MCMC and
not a claim that a multimodal posterior is Gaussian. If the optimizer is not at a
verified stationary point, the information matrix is singular/non-finite, or the
MAP lies on a parameter bound, the MAP and residuals remain available but the
covariance/bands are replaced by a structured failure. Short, clustered ranges
therefore produce broad/prior-dominated BC uncertainty or an explicit inability
to approximate it, never a precise-looking fixed BC.

The diagnostics state how the MAP was verified. `scaled_gradient` means the
Gauss-Newton half-gradient met its scaled tolerance. If the production solver's
fine numerical texture disagrees with that broad finite-difference stencil,
`objective_mesh` means a deterministic eight-direction poll of the actual
penalized chi-square found no improvement larger than `1e-8` down to a scaled
`1e-7` mesh; its radius, largest observed improvement, and evaluation count are
reported rather than hiding the alternate convergence criterion.

This uncertainty surface currently runs in the native CLI and library. With no
uncertainty flag, `true-velocity` takes the existing path unchanged: point
estimates, table/CSV text, and JSON schema remain compatible.

### Solution Diff Attribution (`explain`) — MBA-1345

`explain`, `error-budget`, and `tolerance` (below) are a small analysis train built on the same
shared perturbation kernel. Like `solve-json`, these take a
[solve-json v1](docs/SOLVE_JSON_V1.md) request file as their primary input rather than
`-v`/`-b`/`-m`/`-d` flags -- a plain JSON document, explicit SI units throughout, that fully
describes one firing solution. You author one by hand (or generate it from another tool);
`ballistics solve-json < request.json` will solve one standalone and echo back the resolved
values if you want to sanity-check it first.

`explain` answers a question none of this crate's other commands do: not "what does this load
do," but *why do two solved requests disagree*. It resolves `--a` and `--b` independently and
attributes the difference in the observed impact, at each requested range, to a symmetric
counterfactual swap of each of seven input groups -- drag, muzzle velocity, zero/sight geometry,
atmosphere, wind, shot geometry, effects. Whatever nonlinear interaction between groups the seven
do not explain is reported once per range as an explicit **interaction remainder**, and is never
distributed across groups.

```bash
cat > a.json <<'JSON'
{"schema_version": 1,
 "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7", "ballistic_coefficient": 0.243},
 "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
 "shot": {"max_range_m": 900.0, "zero_distance_m": 100.0},
 "atmosphere": {}, "wind": {"speed_mps": 3.0, "direction_from_rad": 1.5707963267948966},
 "solver": {}, "effects": {}, "sampling": {"interval_m": 25.0}}
JSON
cat > b.json <<'JSON'
{"schema_version": 1,
 "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7", "ballistic_coefficient": 0.243},
 "rifle": {"muzzle_velocity_mps": 870.0, "sight_height_m": 0.05},
 "shot": {"max_range_m": 900.0, "zero_distance_m": 100.0},
 "atmosphere": {}, "wind": {"speed_mps": 6.0, "direction_from_rad": 1.5707963267948966},
 "solver": {}, "effects": {}, "sampling": {"interval_m": 25.0}}
JSON
ballistics explain --a a.json --b b.json --ranges 600
```

```
explain -- method: symmetric_group_counterfactual

range 600.000 m
  total                  drop -0.375886 m   windage -0.540216 m   time -0.057616 s   velocity  +35.9830 m/s
  ProjectileDrag         drop +0.000000 m   windage +0.000000 m   time +0.000000 s   velocity   +0.0000 m/s
  MuzzleVelocity         drop -0.375901 m   windage +0.082054 m   time -0.057619 s   velocity  +35.9817 m/s
  ZeroSightGeometry      drop +0.000000 m   windage +0.000000 m   time +0.000000 s   velocity   +0.0000 m/s
  Atmosphere             drop +0.000000 m   windage +0.000000 m   time +0.000000 s   velocity   +0.0000 m/s
  Wind                   drop +0.000015 m   windage -0.622271 m   time +0.000003 s   velocity   +0.0013 m/s
  ShotGeometry           drop +0.000000 m   windage +0.000000 m   time +0.000000 s   velocity   +0.0000 m/s
  Effects                drop +0.000000 m   windage +0.000000 m   time +0.000000 s   velocity   +0.0000 m/s
  interaction remainder  drop +0.000000 m   windage +0.000000 m   time +0.000000 s   velocity   +0.0000 m/s   (unexplained by any single group)
```

`total` is the raw observed difference at that range (`b` minus `a`); each group line is that
group's own symmetric counterfactual contribution -- the mean of swapping just that group's
inputs in each direction, so the split does not depend on which file is `--a` and which is `--b`.
Here `b.json` differs from `a.json` in both muzzle velocity (+47 m/s) and wind call (+3 m/s), and
the table attributes the drop change almost entirely to `MuzzleVelocity` and the windage change to
a mix of `MuzzleVelocity` (a faster bullet spends less time exposed to the crosswind) and `Wind`
itself -- with a near-zero interaction remainder, meaning these two particular changes combine
almost additively at this range. `-o json` carries the identical numbers as a versioned `SolutionDiffReportV1`
(`schema_version`, `method`, `assumptions`, `skipped_axes`, `rows[].total` /
`.contributions[]` / `.interaction_remainder`); a table run also lists any `skipped_axes` (an
input that could not be swapped on one or both requests at all, e.g. a wind axis absent under
segmented wind) with the group and reason.

**Units.** `--ranges` is always meters, solve-json v1's own SI unit, and is never affected by
`--units` -- `explain` has no `--target`/`--sigma`/`--domain`-style flag for the global `--units`
flag to apply to, so it is simply inert here.

**Honest limits**, carried in the report's own `assumptions` field, not only here: group
contributions are symmetric counterfactuals, so the split never depends on replacement order or
which file is `--a`/`--b`; nonlinear interaction between groups is reported ONLY as the explicit
interaction remainder and is never distributed across groups -- for genuinely correlated inputs no
unique causal attribution exists, and this command does not manufacture one; and an axis that
could not be swapped for one or both requests is excluded from BOTH directions of its group (see
`skipped_axes`), with any real effect it would have had folded into the interaction remainder
rather than silently dropped or partially attributed.

### Per-Input Error Budget (`error-budget`) — MBA-1347

A shooter usually has time to improve exactly one input before a shot: rezero, get a better wind
call, chronograph another string. `error-budget` takes a solve-json v1 request (like `explain`,
above) and a set of declared one-sigma input uncertainties (`--sigma AXIS=VALUE`, repeatable,
kebab-case axis names, e.g. `wind-speed`, `muzzle-velocity-mps`, `cant` -- an unrecognized name is
rejected with the full accepted list), propagates each one to impact covariance at every
requested range through the same central-difference kernel `explain`/`tolerance` share, and ranks
the sources by their own share of impact variance. Every declared source is kept individually and
ranked -- **never** collapsed into an "other" bucket the way `monte-carlo --wez`'s attribution is
-- so the report ends with a concrete answer: which one input is worth improving here.

```bash
ballistics error-budget --request a.json --ranges 600 --units imperial \
  --target rect:10x15 \
  --sigma wind-speed=1.5 --sigma muzzle-velocity-mps=10.0
```

```
error-budget -- method: central_difference_first_order_propagation_gl20_panelled_pm6sigma

range 600.000 m
  impact spread: sigma_drop 0.087680 m   sigma_windage 0.325048 m   covariance -0.001101 m^2
  95% ellipse: semi-major 0.795650 m   semi-minor 0.214436 m   area 0.536007 m^2
  p_hit: 29.5%
  ranked sources (share of impact variance):
     1. wind-speed                     93.1%   gain if perfected: p_hit +67.523%
     2. muzzle-velocity-mps             6.9%   gain if perfected: p_hit +0.909%
  priority: WindSpeed dominates at 600 m (93.1% of impact variance). Measuring it better is the highest-value single improvement here.
```

A 1.5 m/s wind-call sigma and a 10 m/s muzzle-velocity sigma leave this 10x15 in target at only a
29.5% hit probability, almost entirely (93.1%) because of the wind call -- and the `gain if
perfected` column turns that into a decision: perfecting the wind call alone would recover 67.5
points of hit probability, while perfecting muzzle velocity recovers under 1 point. That gain
column, computed from `p_hit_gain_if_perfect`, appears only when `--target` is supplied (omit it
to get the ranked sources and impact spread without a hit probability). A source the kernel
cannot evaluate for this request -- a categorical axis, a wind axis absent under segmented wind,
one outside its own physical domain -- is named in a separate `unavailable sources` section with
its reason, and **never** silently dropped from the ranking: a dropped source would read as
"contributes no uncertainty," a different fact from "could not be measured." `-o json` carries
the identical numbers as a versioned `ErrorBudgetReportV1`, including each source's finite-
difference `scheme` (the table names it only when it is not the ordinary central difference) and
`unavailable_sources[].code`/`.reason`.

**Units.** Two conventions are in play in the example above: `--target` follows `--units`
(inches under imperial, cm under metric -- matching `tolerance --target` below), while `--ranges`
and every `--sigma` value are always in solve-json v1's own SI unit (meters, m/s, radians, ...),
regardless of `--units`.

**Honest limits**, carried in the report's own `assumptions` field: declared sources are treated
as independent, with no correlation between them modelled; propagation is first-order
(local-linear) about the nominal solution, using each axis's own small default differencing step
regardless of the declared sigma -- a large declared sigma is a linear extrapolation from a slope
measured over a much smaller window, not exact for large or non-Gaussian input uncertainty; the
95% ellipse and the hit probability both assume an approximately Gaussian impact distribution; a
source's derivative may be one-sided rather than central at a physical domain boundary (still air
for wind speed is the routine case, not an exotic one); and hit probability reflects only the
declared input uncertainty, never model error, over a target always centred on the nominal
(zero-mean) impact point.

### Tolerance Envelopes (`tolerance`) — MBA-1350

`tolerance` answers the complementary question to `error-budget`: not "how much uncertainty does
this input contribute," but "how wrong may this one input be before the shot leaves the target."
It takes a solve-json v1 request, resolves it, and bisects each requested `--axis` outward --
toward its own configured `--domain`'s lower bound and, independently, its upper bound -- until
the impact crosses the target boundary. `--domain AXIS=LO:HI` is required once per `--axis`, in
that axis's own physical unit (no default domain is ever invented, even for a categorical axis,
whose bounds then simply go unused).

```bash
ballistics tolerance --request a.json --range 600 --units imperial --target rect:10x15 \
  --axis wind-speed --domain wind-speed=0:12 \
  --axis muzzle-velocity-mps --domain muzzle-velocity-mps=750:900
```

```
tolerance -- method: one_variable_deterministic_bisection
range 600.000 m

axis wind-speed (nominal 3.000000 m/s, target margin 0.1270 m)
  near (toward the domain's lower bound): bound at 2.413496 m/s -- crosses the Right edge
  far  (toward the domain's upper bound): bound at 3.586503 m/s -- crosses the Left edge

axis muzzle-velocity-mps (nominal 823.000000 m/s, target margin 0.1270 m)
  near (toward the domain's lower bound): bound at 802.154006 m/s -- crosses the Bottom edge
  far  (toward the domain's upper bound): bound at 845.896312 m/s -- crosses the Top edge
```

On this 10x15 in target, the wind call may drift only to 2.41-3.59 m/s (well inside the
configured 0-12 m/s search domain) before the shot leaves the target sideways, while muzzle
velocity has much more room, 802-846 m/s, before it leaves top or bottom. A `--domain` searched to
its edge with no crossing reports `no bound within the configured domain` (the impact never
leaves the target there) rather than extrapolating past it; an axis that provably cannot move the
impact at this fixed range at all (`target-distance`, for example, which only matters when the
range itself changes) reports `no measurable effect on the impact` instead -- a different fact
from merely being unbounded, so the two are never rendered with the same text. An axis the kernel
refuses outright (a categorical toggle, an axis absent under segmented wind, ...) appears in a
separate `unavailable axes` section instead of a fabricated bound.

**`target-distance`'s domain must not reach below `--range`.** The `no measurable effect on the
impact` claim above holds only while `--domain target-distance=LO:HI`'s own `LO` stays at or
above `--range`: bisecting `target-distance` (`shot.max_range_m`) DOWNWARD past the range you are
observing at asks the solver for an impact beyond where the trajectory was even computed, which
fails outright (a genuine, whole-command-aborting error, not a graceful `unavailable axes` entry).
`target-distance` is also not an answer to "how far off could my rangefinder reading be" in the
first place -- it changes how far the trajectory is COMPUTED, never the muzzle angle, so it has no
effect on the impact at a fixed range at all. `zero-distance` is the axis that actually answers
that question: it re-runs the elevation search for a different assumed zero and observes the
impact at your true, fixed `--range`, exactly what an inaccurate rangefinder does to a real shot.

**Units.** The same mixed convention as `error-budget` above: `--target` follows `--units`
(inches under imperial, cm under metric), while `--range` and every `--domain` bound are always
in solve-json v1's own SI unit (meters, m/s, radians, ...), regardless of `--units`.

`-o json` carries the identical numbers as a versioned `ToleranceReportV1` (`schema_version`,
`method`, `assumptions`, `range_m`, `axes[]` with each axis's `nominal`, `near_bound`/`far_bound`,
`near_limiting_boundary`/`far_limiting_boundary`, `unbounded_in_domain`, and the
`near_has_no_effect`/`far_has_no_effect` flags) and `unavailable_axes[].code`/`.reason`.

**Honest limits**, carried in the report's own `assumptions` field: each bound holds ONE input at
its limit while every other input stays at its nominal value -- bounds from different axes may
**not** be assumed to hold simultaneously, since two inputs each at their own individual limit
will generally miss even though neither alone would; no probability distribution is assumed or
implied by any bound here, these are deterministic limits of a one-variable search, not confidence
intervals or a probability of hit; and a bound is reported only when found strictly within the
axis's own configured search domain, never extrapolated beyond it.

### Constrained Dial & Hold Planning (`dial-plan`) — MBA-1348

Turns an already-known TRUE angular correction into ranked, EXECUTABLE plans against a real
optic: whole clicks with a direction and (when the turret declares `clicks_per_revolution`) a
revolution count, a reticle hold in continuous TRUE angular mil, or a hybrid split (dial what
the turret can reach, hold the true angular remainder). All three strategies -- `dial_all`,
`hold_all`, `hybrid` -- are always evaluated and ranked feasible-first. An optic that cannot
fully execute a strategy is never silently clamped: the plan is still returned, `feasible:
false`, with the limiting mechanism named (`TravelExceeded`, `HoldBoundExceeded`,
`NoTravelData`, `NoHoldBoundData`, ...) -- an infeasibility analysis is a normal, successful
run (exit 0), never an error.

The optic comes from `--profile NAME` (its saved turret mechanics, hold bounds, and
`elevation_cf`/`windage_cf`) or from inline flags mirroring `profile save`'s own set
(`--elevation-click`, `--travel-up`/`--travel-down`, `--windage-travel-left`/
`--windage-travel-right`, `--hold-up`/`--hold-down`/`--hold-left`/`--hold-right`,
`--turret-elev`/`--turret-wind`, `--clicks-per-rev`, `--zero-stop`) -- the two sources are
mutually exclusive.

**The come-ups → dial-plan pairing.** `come-ups` computes a correction from a load and a
zero; `dial-plan` turns an already-known correction into an executable plan. Run them in
sequence:

```bash
ballistics come-ups -v 2700 -b 0.475 -m 175 -d 0.308 --drag-model g7 \
  --zero-distance 100 --start 100 --end 600 --step 100
```
```

Come-Up Table (zero: 100 yd, MIL)
┌──────────┬──────────┬──────────┬──────────┬──────────┬──────────┐
│Range (yd)|Drop (MIL)|Come-Up   │ Vel (fps)│Energy    │ Time (s) │
├──────────┼──────────┼──────────┼──────────┼──────────┼──────────┤
│      100 │    0.000 │    —     │     2604 │     2634 │    0.113 │
│      200 │    0.434 │    0.434 │     2509 │     2446 │    0.231 │
│      300 │    1.090 │    0.656 │     2417 │     2269 │    0.352 │
│      400 │    1.832 │    0.742 │     2326 │     2102 │    0.479 │
│      500 │    2.634 │    0.802 │     2237 │     1945 │    0.610 │
│      600 │    3.490 │    0.856 │     2151 │     1797 │    0.747 │
└──────────┴──────────┴──────────┴──────────┴──────────┴──────────┘
```

At 600 yards the "Drop (MIL)" column is the total elevation correction from zero, 3.490 mil --
feed that straight into `dial-plan`:

```bash
ballistics dial-plan --elevation 3.49mil --range 600 \
  --elevation-click 0.1mil --travel-up 30mil --travel-down 5mil
```
```
dial-plan -- method: dial_space_quantization_v1
range: 548.640 m

#1 strategy: dial_all [FEASIBLE]
  elevation: dial UP 35 clicks; hold 0.000 mil
  windage: dial RIGHT 0 clicks; hold 0.000 mil
  residual: elevation -0.010 mil, windage 0.000 mil -> 0.005 m at 548.640 m

#2 strategy: hybrid [INFEASIBLE]
  elevation: dial UP 35 clicks; hold -0.010 mil
  windage: dial RIGHT 0 clicks; hold 0.000 mil
  residual: elevation 0.000 mil, windage 0.000 mil -> 0.000 m at 548.640 m
  limit: elevation NoHoldBoundData -- needed -0.010 mil, no data declared for this axis

#3 strategy: hold_all [INFEASIBLE]
  elevation: dial UP 0 clicks; hold 3.490 mil
  windage: dial RIGHT 0 clicks; hold 0.000 mil
  residual: elevation 0.000 mil, windage 0.000 mil -> 0.000 m at 548.640 m
  limit: elevation NoHoldBoundData -- needed 3.490 mil, no data declared for this axis

assumptions:
  - Linear miss at range uses the small-angle approximation (mil / 1000 * range); it is not exact at extreme angles.
  - Elevation and windage are planned independently; no cant-induced coupling between axes is modeled.
  - Reticle holds are assumed continuous and unquantized, unlike turret clicks.
  - Travel limits and turret state are trusted exactly as declared in the optic profile, not sensed or independently verified.
  - MOA-graduated clicks convert to milliradians using the locked printed-table constant 3438, not the exact geometric 3437.7467.
```

`dial_all` wins: 35 clicks up on a 0.1 mil turret with 30 mil of declared up-travel leaves only
a -0.010 mil residual (5 mm of linear miss at 600 yards) -- negligible, and feasible outright.
`hybrid`/`hold_all` are marked infeasible here only because this optic declared no
`--hold-up`/`--hold-down`/`--hold-left`/`--hold-right`; add all four (a real reticle's usable
extent) to make those two strategies evaluable as well.

**Units.** `--range` follows `--units` (yards imperial / meters metric) like every other
calculator command, but the report's own `range_m` field -- and everything the table above
prints -- is always meters, regardless of `--units`: the example requested 600 yards and the
report says `range: 548.640 m`. `--elevation`/`--windage` (and every inline optic travel/hold/
turret-state flag) are TRUE angular with a mandatory unit suffix (mil/moa/smoa/iphy); `+` means
up/right on both axes.

`-o csv`/`-o pdf` are rejected by name (`dial-plan has no CSV/PDF form; use -o table or -o
json`); `-o json` is `DialPlanReportV1` pretty-printed verbatim (`schema_version`, `method`,
`assumptions`, `range_m`, `plans[]`).

**Honest limits**, carried in the report's own `assumptions` field (quoted above verbatim): the
linear-miss figure is a small-angle approximation, not exact at extreme angles; elevation and
windage are planned independently with no cant-induced coupling modeled; reticle holds are
continuous and unquantized, unlike turret clicks; declared travel/turret state is trusted
exactly as given, never sensed or independently verified; and MOA-graduated clicks convert on
this crate's locked printed-table constant (3438), not the exact geometric 3437.7467.

### Adaptive Range Cards (`adaptive-card`) — MBA-1351

The smallest range card that PROVABLY reconstructs the trajectory within a stated elevation/
windage error budget: greedy worst-point insertion over the solved trajectory's own sample
grid, starting from both range ends plus any `--anchor`s and adding whichever audited point is
furthest outside budget until none is (or the row cap binds, or the remaining error is an
irreducible click-quantization floor). A SEPARATE dense pass then MEASURES the finished card
and reports its true worst-case error and the grid it was verified against.

**What this buys you, and what it does not.** Measured against the most generous possible
fixed-step card, adaptive placement does **not** reliably produce fewer rows than a
well-chosen uniform step on a smooth trajectory -- the engine's own
`fixed_step_comparison_is_measured_not_assumed` test (`src/card.rs`) found adaptive placement
lost 10, tied 5, and won 5 across five domains and four budgets. What it buys instead: a
MEASURED worst-case error bound (not an assumed one), every range you asked for as an
`--anchor` always present in the card, and no step size to guess.

```bash
ballistics adaptive-card -v 2700 -b 0.475 -m 175 -d 0.308 --drag-model g7 \
  --zero-distance 100 --start 300 --end 700
```
```

Adaptive Range Card (zero: 100 yd) -- method: greedy_worst_point_insertion_on_holdcurve_grid_v1
┌───────┬─────────┬─────────┬─────────┬─────────┬─────────┬─────────┬───────┐
│Range  │Drop     │Drop     │Wind     │Wind     │Vel      │Energy   │ToF    │
│(yd)   │(in)     │(MIL)    │(in)     │(MIL)    │(fps)    │(ft-lb)   │(s)    │
├───────┼─────────┼─────────┼─────────┼─────────┼─────────┼─────────┼───────┤
│   300 │    11.8 │   1.090 │     0.0 │    0.00 │    2417 │    2269 │  0.35 │
│   494 │    46.0 │   2.584 │     0.0 │    0.00 │    2243 │    1954 │  0.60 │
│   700 │   110.9 │   4.400 │     0.0 │    0.00 │    2066 │    1659 │  0.89 │
└───────┴─────────┴─────────┴─────────┴─────────┴─────────┴─────────┴───────┘

budget met: yes
rows: 3 of 25 max
worst error: elevation 0.0286 MIL, windage 0.0000 MIL @ 395.000 yd
verification grid: 1.0000 yd step
```

With no `--elevation-budget`/`--windage-budget` given (as here), the default is half the
profile's own click size -- the honest floor documented below -- when `--profile` supplies
one, or 0.1 true mil otherwise. Three rows reconstruct 300-700 yards within that 0.1 mil
default: interpolating linearly between any two adjacent rows above is never more than 0.0286
mil off the real trajectory, worst case, at 395 yards -- and that worst point, and the grid it
was measured against, are stated rather than assumed. (The `(ft-lb)` header cell overflowing
its column by one character is a pre-existing cosmetic quirk shared with `range-table`'s
identical imperial header, not something new here.)

Table/CSV/PDF (`-o table|csv|pdf`) all carry the same footer block shown above (`budget met` /
`rows` / `worst error @ range` / `verification grid`) -- `-o csv`/`-o pdf` write it to STDERR
so their own data channel (the CSV rows; the PDF file) stays exactly that and nothing else.
`-o json` carries no separate footer, because every fact in it is already a field on the
report:

```bash
ballistics adaptive-card -v 2700 -b 0.475 -m 175 -d 0.308 --drag-model g7 \
  --zero-distance 100 --start 300 --end 400 \
  --elevation-budget 0.05mil --windage-budget 0.05mil -o json
```
```
{
  "schema_version": 1,
  "method": "greedy_worst_point_insertion_on_holdcurve_grid_v1",
  "assumptions": [
    "Verification is limited to the hold curve's declared sample grid (verification_grid_step_m) together with the card's own rows; no claim is made about ranges between those audited points.",
    "The reader of this card interpolates linearly between adjacent rows.",
    "Errors are measured in printed-value space -- the same zero-set bias, tracking-correction division and click quantization the printed rows carry -- so a constant zero-set bias cancels out of the interpolation error and the tracking correction factor is already inside the numbers being compared.",
    "Rows quantized to an optic's clicks carry an irreducible error of up to half a click at the rows themselves, which no number of added rows can remove.",
    "A budget tighter than that half-click floor is reported as budget_met: false; the requested tolerance is never silently relaxed."
  ],
  "rows": [
    {
      "range": 274.32,
      "drop_linear": 0.29906997487840825,
      "drop_adj": 1.0902230055351714,
      "come_up": null,
      "wind_linear": 0.0,
      "wind_adj": 0.0,
      "velocity": 736.5692110050417,
      "energy": 3076.1171845457698,
      "time": 0.35235712510991274,
      "lead_adj": null,
      "wind_columns": []
    },
    {
      "range": 365.76,
      "drop_linear": 0.670060323809346,
      "drop_adj": 1.8319672020159286,
      "come_up": null,
      "wind_linear": 0.0,
      "wind_adj": 0.0,
      "velocity": 708.9565577781669,
      "energy": 2849.8040653051507,
      "time": 0.4788992677109733,
      "lead_adj": null,
      "wind_columns": []
    }
  ],
  "budget_met": true,
  "rows_capped": false,
  "worst_elevation_error": 0.008452678351470011,
  "worst_windage_error": 0.0,
  "worst_error_range_m": 318.2112,
  "verification_grid_step_m": 0.9144
}
```

**Units.** `--start`/`--end`/`--anchor` follow `--units` for input (yards imperial / meters
metric), exactly like `range-table`, and the Table/CSV/PDF Range column is likewise in
`--units`. `-o json`, however, is `AdaptiveCardReportV1` pretty-printed **verbatim**:
`range`/`drop_linear`/`wind_linear`/`worst_error_range_m`/`verification_grid_step_m` stay in
the engine's own native METRES regardless of `--units` -- notice the example above requested
`--start 300` (yards) and the JSON's first row range reads `274.32` (= 300 x 0.9144, not
re-converted for display). `drop_adj`/`wind_adj` are the one exception either way: they are
already expressed in the `--adjustment` unit (`mil` or `moa`) in every output format, and are
never re-converted. `--elevation-budget`/`--windage-budget` are TRUE angular with a mandatory
unit suffix (mil/moa/smoa/iphy), converted into whichever `--adjustment` unit is selected --
the same conversion rule the default (profile-click-derived) budget uses, so an explicit
budget and the default can never silently disagree about which space they live in.

**Table/CSV Range formatting.** `HoldCurve`'s native grid is spaced every 0.9144 m, so under
`--units metric` almost every algorithm-chosen interior row -- and any `--anchor` that doesn't
happen to land on a whole display unit, in either unit system -- is genuinely fractional. The
Table and CSV Range column render it as a bare integer only when within floating-point noise
of one, and with exactly one decimal place otherwise (the identical rule the PDF path's
`format_range` already used) -- **not** rounded to the nearest whole unit, which would silently
disagree with the exact value `-o json` and `-o pdf` both carry for that same row. For example,
requesting `--anchor 412.5` on a 200-500 m card:

```bash
ballistics --units metric adaptive-card -v 800 -b 0.223 -m 10.9 -d 7.82 --drag-model g7 \
  --sight-height 45 --zero-distance 100 --start 200 --end 500 --anchor 412.5
```
```

Adaptive Range Card (zero: 100 m) -- method: greedy_worst_point_insertion_on_holdcurve_grid_v1
┌───────┬─────────┬─────────┬─────────┬─────────┬─────────┬─────────┬───────┐
│Range  │Drop     │Drop     │Wind     │Wind     │Vel      │Energy   │ToF    │
│( m)   │(mm)     │(MIL)    │(mm)     │(MIL)    │(m/s)    │(   J)   │(s)    │
├───────┼─────────┼─────────┼─────────┼─────────┼─────────┼─────────┼───────┤
│   200 │   137.6 │   0.688 │     0.0 │    0.00 │     670 │    2449 │  0.27 │
│ 307.2 │   529.8 │   1.724 │     0.0 │    0.00 │     606 │    2003 │  0.44 │
│ 412.5 │  1214.2 │   2.943 │     0.0 │    0.00 │     546 │    1628 │  0.62 │
│   500 │  2057.7 │   4.115 │     0.0 │    0.00 │     499 │    1357 │  0.79 │
└───────┴─────────┴─────────┴─────────┴─────────┴─────────┴─────────┴───────┘

budget met: yes
rows: 4 of 25 max
worst error: elevation 0.0274 MIL, windage 0.0000 MIL @ 251.460 m
verification grid: 0.9144 m step
```

The requested anchor prints as `412.5`, exactly what `-o json` reports for that row, not `412`
-- and the algorithm-chosen row at the curve's own native grid point 307.2384 m prints `307.2`,
the same one-decimal rounding `-o pdf`'s Range column already applied before this fix.

The optic (for click quantization) and its tracking correction factors come only from
`--profile NAME` here -- unlike `dial-plan`, there is no inline-flag optic path for
`adaptive-card`, since quantization is optional for a range card: the exact, unrounded angle
is a perfectly good answer with no optic declared at all.

**Honest limits**, carried in the report's own `assumptions` field (quoted above verbatim):
verification is limited to the hold curve's declared sample grid together with the card's own
rows -- no claim is made about ranges between those audited points; the reader of the card is
assumed to interpolate linearly between adjacent rows; errors are measured in printed-value
space (the same zero-set bias, tracking-correction division, and click quantization the
printed rows themselves carry); rows quantized to an optic's clicks carry an irreducible error
of up to half a click at the rows themselves, which no number of added rows can remove; and a
budget tighter than that half-click floor is reported `budget_met: false` -- the requested
tolerance is never silently relaxed.

## Output Formats

### Units by Output Surface

The same solved trajectory is expressed differently depending on which output you
ask for — when comparing numbers across surfaces, check which one you are reading:

| Surface | Vertical | Units |
|---|---|---|
| Native `-o json` (`--full` points) | world-frame `x`/`y`/`z`, `y` = height above ground | yd (imperial) / m (metric), per its legend |
| Native table / CSV drop & drift | below the line of sight | inches (imperial) / meters (metric) |
| `solve-json` v1 | below the line of sight (`drop_m`) | always SI (meters) |
| WASM terminal `-o json` | below the line of sight, unit in the key name | `drop_inches`/`drift_inches` (imperial), `drop_cm`/`drift_cm` (metric), per its `units` legend |

Two related cross-system notes:

- The **default bore height** is the round number of each system: 60 in (imperial)
  vs 1500 mm (metric) — 1.524 m vs 1.5 m. The 2.4 cm difference is visible in
  `max_height` and in high-precision imperial-vs-metric parity checks; pass
  `--bore-height` explicitly when comparing systems digit-for-digit.
- The **mover ring** renders as an angular hold (`--adjustment-unit`: mil, moa, smoa,
  iphy, or clicks — see [Turret Adjustment Units](#turret-adjustment-units)) in the
  human table, while CSV keeps `ring_mil` and JSON keeps `mover_ring_m` /
  `mover_ring_mil` regardless — machine columns carry their unit in the name and
  stay stable.


All commands support four output formats via `-o`:

### Table Format (default)
```bash
./ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 -o table
```

### JSON Format
```bash
./ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 -o json > trajectory.json
```

**Units and axes (MBA-1315).** The legacy `trajectory[]` points are bare `x`/`y`/`z` numbers
with no unit suffix in the field name, and their axis order surprises tooling written against
a "x=lateral, y=up, z=depth" 3D convention — in this output `x` is lateral and `z` is
downrange. A field tester's tooling once misread `x`/`z` (in **yards**, not feet) as a depth
axis in feet and misdiagnosed a bug. The document carries a `legend` block (appended after
`trajectory`; every pre-existing key/value is unchanged) that states both explicitly:

```json
{
  "units": "imperial",
  "max_range": 25.0,
  "max_height": 1.6666666666666667,
  "time_of_flight": 0.02803675408939502,
  "impact_velocity": 2651.8043404763275,
  "impact_energy": 2622.759081598047,
  "stability_coefficient": 2.000083126080045,
  "spin_drift": null,
  "trajectory": [
    {
      "time": 0.0,
      "x": 0.0,
      "y": 1.6666666666666667,
      "z": 0.0,
      "velocity": 2700.0,
      "energy": 2718.96097071048
    }
  ],
  "legend": {
    "units": {
      "distance": "yd",
      "velocity": "fps",
      "energy": "ft-lb"
    },
    "axes": {
      "x": "lateral offset from the muzzle's initial aiming direction; positive means to the shooter's right (e.g. a crosswind FROM the left, --wind-direction 270, drifts x positive; FROM the right, --wind-direction 90, drifts x negative). Zero at the muzzle.",
      "y": "height above the ground in the world frame; positive means up. Starts near bore height and reaches 0 at ground impact. This is NOT height above the line of sight (compare solve-json v1's LOS-relative drop_m).",
      "z": "downrange distance from the muzzle; zero at the muzzle, always increasing."
    }
  }
}
```

- **`legend.units`** — concrete abbreviation for each numeric quantity family: `distance`
  (`trajectory[].x`/`y`/`z`, `max_range`, `max_height`, `spin_drift`), `velocity`
  (`trajectory[].velocity`, `impact_velocity`), and `energy` (`trajectory[].energy`,
  `impact_energy`). These are `yd`/`fps`/`ft-lb` under `--units imperial` (the default) and
  `m`/`m/s`/`J` under `--units metric` — the same mapping as the top-level `units` field, spelled
  out numerically instead of left for tooling to infer. Angle fields (`max_yaw_angle`,
  `max_precession_angle`, `final_pitch_angle`, `final_yaw_angle`, present only with
  `--enable-precession`) are always radians; time fields are always seconds; neither varies
  with `--units`, so neither is covered by `legend.units`.
- **`legend.axes`** — verified empirically (not assumed) from controlled solves: a pure
  crosswind run and a no-wind run, comparing which point component moved, in which sign
  direction, against the table output. `x` is **lateral** (positive = shooter's right; wind
  FROM the left, `--wind-direction 270`, drifts it positive), `y` is **height above the
  ground** (positive = up; this is a world-frame height, NOT height above the line of sight),
  and `z` is **downrange distance** from the muzzle (always increasing). Note the reversal
  from some 3D conventions: `x` is lateral and `z` is downrange here, not the other way round.

**Mover ring** (`mover_ring_m`/`mover_ring_mil`, only with `--target-speed > 0`) and the
`min_pitch_damping`/`transonic_mach`/`max_yaw_angle`/`max_precession_angle`/
`final_pitch_angle`/`final_yaw_angle` diagnostics (only with `--enable-pitch-damping` /
`--enable-precession`) carry their unit in the field name already and are covered above and
in the [Mover Ring](#mover-ring---target-speed) section, not restated in `legend`.

**`drag_coefficient`** (MBA-1423): present on each `trajectory[]` point only when
`--with-drag-coefficient` was supplied alongside `--full`. Dimensionless, so it does not vary
with `--units`.

This is the **projectile's own** drag coefficient, not the reference table's. A G-model Cd
describes the *standard* projectile, so reporting it directly would hand back the same G7
curve for every bullet sharing a drag model. It is scaled by the form factor to describe this
bullet, from the identity that both descriptions must produce the same retardation:

```
Cd_own / SD == Cd_ref / BC    =>    Cd_own == Cd_ref × SD / BC
```

Consequences worth knowing before you chart it:

- A **velocity- or Mach-segmented BC** (`--bc-segment`) shows its band steps here, because the
  BC in that denominator is the one resolved for that point. This is the part of a real load's
  drag curve you cannot reconstruct from a published BC.
- A **custom drag table** (`--drag-table`) already supplies the projectile's actual Cd and is
  divided by sectional density rather than a BC, so the scale factor collapses to exactly 1
  and the curve passes through unchanged (including `--cd-scale`).
- It covers **only the speeds the bullet actually flew** — it is not a 0–3 Mach sweep. If you
  need the reference curves themselves across their full domain, that is a separate concern.
- Cd is paired with the **station** speed of sound, the same one the document reports Mach
  against, so plotting Cd against Mach is self-consistent. Using each step's local speed of
  sound would attribute a Cd to a Mach this document never shows.
- **Absent** (not `null`) on a point whose sectional density is unknown, since the
  projectile's own Cd is undefined without both mass and diameter.

JSON only. CSV's sampled-interval branch (`--sample-interval`) reports a different point type
that carries no Cd, so a column there would be populated in one branch and empty in the other.
Without the flag, output is byte-identical to before this field existed.

**Browser terminal:** `--with-drag-coefficient` works there too on `-o json` (MBA-1427). The
value comes straight off the same solver-annotated points — it is not re-derived in the WASM
layer, which would reintroduce exactly the drift this feature exists to prevent. One pre-existing
surface difference to know: the browser JSON always includes the per-point array, while native
gates it on `--full`.

**`zero_angle_degrees`** (MBA-1402): present only when auto-zero ran (`--auto-zero`, or a
saved profile's `auto_zero`/`zero_distance`) — the same bore angle `TrajectoryConfig.angle`
was set to for this run, echoed in degrees so it can be captured and fed straight into
`zero --from-angle` later (see [Solving Range From a Stored Angle](#solving-range-from-a-stored-angle---from-angle)).
Absent (not `null`) on a bare `--angle` run. Always degrees regardless of `--units`, matching
the `zero` command's own `zero_angle_degrees` convention. `-o csv`'s non-`--full` summary adds
a matching `zero_angle_degrees,<value>,degrees` row under the same condition; the table
summary adds a `Zero Angle: <value>°` line. Every surface is additive: with auto-zero absent,
output is byte-identical to before this field existed.

**Note:** `--auto-zero` itself is not a new flag. This means every existing `--auto-zero`
invocation gains a new table line / JSON key / CSV row / WASM banner term as of this
release, on a flag that predates it — anything already parsing `trajectory --auto-zero`
output (scripts, the WASM banner text, downstream bindings) should be updated to expect it.

**WASM `-o json`** (browser/`ballistics.sh` interface) uses a different, already
self-describing shape — `range_yards`/`drop_inches`/`drift_inches` (or the `_meters`/`_cm`
metric equivalents) instead of bare `x`/`y`/`z` — but until MBA-1315 it still left the sign
convention of `drop`/`drift` unstated. It carries the same `legend` block, adapted to its own
field names:

```json
{
  "trajectory": [ { "range_yards": 0.0, "drop_inches": 0.0, "drift_inches": 0.0, "velocity_fps": 2700.0, "energy_ftlb": 2718.96, "time_seconds": 0.0 } ],
  "summary": { "max_range_yards": 25.0, "max_height_inches": 60.0, "time_of_flight_seconds": 0.028, "impact_velocity_fps": 2651.8 },
  "legend": {
    "units": { "range": "yd", "drop": "in", "drift": "in", "velocity": "fps", "energy": "ft-lb" },
    "axes": {
      "range": "downrange distance from the muzzle; zero at the muzzle, always increasing.",
      "drop": "vertical miss from the line of sight; positive means the bullet is below the line of sight (has fallen below the aim point). Not the same reference as the native CLI legacy JSON's world-frame `y`.",
      "drift": "lateral miss from the line of sight; positive means to the shooter's right (e.g. a crosswind FROM the left, --wind-direction 270, drifts it positive; FROM the right, --wind-direction 90, drifts it negative). Same sign and source as the native CLI legacy JSON's `x`."
    }
  }
}
```

`drift` is read off the same underlying value as the native CLI's `x` (both are `position.z`
in the engine's internal frame), so the two agree in sign; `drop` is LOS-relative (subtracts
the line-of-sight height) where the native CLI's `y` is an absolute world-frame height, so the
two do **not** share a reference despite both being "vertical".

### CSV Format
```bash
./ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 -o csv > trajectory.csv
```

### PDF Dope Card Format
Generate a printable dope card with two-column layout, color-coded values, and alternating row stripes for field readability:

```bash
./ballistics trajectory -v 2550 -b 0.236 -m 175 -d 0.308 --drag-model g7 \
  --auto-zero 100 --max-range 1000 \
  --wind-speed 5 --wind-direction 90 \
  --temperature 55 --pressure 27.32 --altitude 2500 \
  --sample-trajectory --sample-interval 9.144 \
  --ignore-ground-impact \
  --target-speed 4 \
  --powder "IMR4320" --bullet-name "SMK" \
  --location-name "General" \
  --profile-row "R700_308" \
  -o pdf --output-file dope_card.pdf
```

**PDF-specific options:**
| Parameter | Description |
|-----------|-------------|
| `--output-file` | Output file path (required for PDF) |
| `--adjustment-unit` | Angular unit for the Drop column and the mover Ring (elevation axis): `mil` (default), `moa`, `smoa`, `iphy`, or `clicks` — see [Turret Adjustment Units](#turret-adjustment-units) |
| `--windage-unit` | Independent angular unit for the Wind/Lead columns (MBA-1410); falls back to `--adjustment-unit` when omitted — see [Turret Adjustment Units](#turret-adjustment-units) |
| `--elevation-click-value` / `--windage-click-value` | Turret click graduations for `--adjustment-unit clicks`, e.g. `0.1mil` or `0.25moa` — see [Turret Adjustment Units](#turret-adjustment-units) |
| `--target-speed` | Target speed for the Lead column (mph imperial / m/s metric — MBA-1325 fixed this to follow `--units` like every other speed flag; previously always read as mph). Also enables the [Mover Ring](#mover-ring---target-speed) column/fields on every output format |
| `--powder` | Powder type (shown in footer) |
| `--bullet-name` | Bullet name (shown in footer) |
| `--location-name` | Location name (shown in header) |
| `--profile-row` | Rifle name (shown in header) |
| `--font-scale` / `--font-preset` | Data-table font size |
| `--bold-data` | Bold font for data cells |

**PDF features:**
- Two-column table layout with Range (yd) and Drop/Wind/Lead in **MIL, MOA, SMOA, IPHY,
  or whole turret clicks** (via `--adjustment-unit`; drop uses the elevation click
  graduation, Wind/Lead the windage one)
- Color coding: Black=Range, Red=Drop, Green=Wind, Blue=Lead
- Alternating row stripes for easy tracking in field conditions
- Header with rifle, location, density altitude, atmospheric data
- Footer with timestamp, load data, BC, and velocity

## Parameters Reference

### Trajectory Command

| Parameter | Description | Default | Imperial | Metric |
|-----------|-------------|---------|----------|--------|
| -v, --velocity | Muzzle velocity | Required | fps | m/s |
| -a, --angle | Launch angle | 0.0° | degrees | degrees |
| -b, --bc | Ballistic coefficient | Required | - | - |
| --bc-reference | BC reference standard (icao/army-standard-metro); see [BC Reference Standard](#bc-reference-standard---bc-reference). Also on `zero`, `monte-carlo`, `profile save` | icao | - | - |
| -m, --mass | Projectile mass | Required | grains | grams |
| -d, --diameter | Projectile diameter | Required | inches | mm |
| --drag-model | Drag model (g1/g2/g5/g6/g7/g8/gi/gs/ra4) | g1 | - | - |
| --auto-zero | Auto-zero distance | None | yards | meters |
| --zero-velocity | Zero-day muzzle velocity (auto-zero only); overrides both powder models | shot-day velocity | fps | m/s |
| --zero-temperature | Zero-day air temperature (auto-zero only); also resolves linear powder velocity unless `--zero-velocity` is set | shot-day temperature | °F | °C |
| --zero-pressure | Zero-day barometric pressure (auto-zero only) | shot-day pressure | inHg | hPa |
| --zero-pressure-type | Whether `--zero-pressure` is absolute or QNH (auto-zero only); see [Pressure Reference](#pressure-reference-absolute-vs-qnh---pressure-type) | shot-day `--pressure-type` | - | - |
| --zero-humidity | Zero-day relative humidity (auto-zero only) | shot-day humidity | percent | percent |
| --zero-altitude | Zero-day altitude (auto-zero only) | shot-day altitude | feet | meters |
| --zero-powder-temp | Zero-day powder temp for the curve lookup (auto-zero only); otherwise uses explicit --zero-temperature, or inherits shot-day --powder-temp when zero temperature is unchanged. --zero-velocity still wins | zero air / inherited shot powder | °F | °C |
| --powder-temp-curve | Measured `TEMP:VEL,...` powder-temp→velocity table (interpolated at the powder temp, clamped; overrides --powder-temp-sensitivity) | none | °F & fps | °C & m/s |
| --powder-temp | With a curve: powder temp the curve is looked up at (default --temperature). With the linear model: reference temp (default 70/21) | --temperature (curve) / 70°F (linear) | °F | °C |
| --sight-height | Sight height above bore | 0.05 | yards | meters |
| --bore-height | Bore height above ground | 5 | feet | meters |
| --ignore-ground-impact | Disable ground impact detection | false | - | - |
| --max-range | Maximum range | 1000 | yards | meters |
| --time-step | Integration time step — RK4/Euler only (the adaptive RK45 default steps adaptively and ignores this) | 0.001 | seconds | seconds |
| --wind-speed | Wind speed | 0 | mph | m/s |
| --wind-direction | Wind direction (0=headwind, 90=from right, 180=tailwind, 270=from left) | 0° | degrees | degrees |
| --wind-vertical | Vertical wind; positive = updraft (raises POI), negative = downdraft. Also on `monte-carlo` | 0 | mph | m/s |
| --wind-segment | Downrange wind segment `SPEED:ANGLE:UNTIL_DISTANCE[:VERTICAL]` (repeatable); VERTICAL is always m/s regardless of `--units` | — | mph & yd | m/s & m |
| --temperature | Temperature | 59 | °F | °C |
| --pressure | Barometric pressure | 29.92 | inHg | hPa |
| --pressure-type | Whether `--pressure` is absolute station pressure or a QNH altimeter setting; see [Pressure Reference](#pressure-reference-absolute-vs-qnh---pressure-type). Also on `zero`, `profile save` | absolute | - | - |
| --humidity | Relative humidity | 50 | % | % |
| --altitude | Altitude | 0 | feet | meters |
| --density-altitude | Density altitude; supersedes `--altitude`/`--pressure`/`--pressure-type` entirely, see [Density Altitude](#density-altitude-as-a-direct-input---density-altitude). Also on `profile save` | — (unset) | feet | meters |
| --use-bc-segments | Enable BC segmentation | false | - | - |
| --bc-segment | Manual velocity-keyed BC segment `VMIN:VMAX:BC` (repeatable) | — | fps | m/s |
| --print-bc-segments | Print the BC5D-generated segment ladder as ready-to-paste `--bc-segment` arguments (requires `--bc-table-dir`) | false | fps | m/s |
| --full | Show all trajectory points | false | - | - |
| --enable-magnus | Enable Magnus effect | false | - | - |
| --enable-coriolis | Enable Coriolis effect | false | - | - |
| --enable-spin-drift | Enable empirical Litz spin drift | false | - | - |
| --twist-rate | Barrel twist rate | 12 | inches/turn | inches/turn |
| --twist-right | Right-hand twist | false | - | - |
| --latitude | Latitude for Coriolis/weather | None | degrees | degrees |
| --longitude | Longitude for weather zones | None | degrees | degrees |
| --shot-direction | Shot azimuth (0=N, 90=E) | None | degrees | degrees |
| --shooting-angle | Incline angle (up/down) | 0 | degrees | degrees |
| --cant | Rifle cant about the line of sight (alias `--cant-angle`); positive = clockwise, POI right and low. Also on `monte-carlo`, not `zero` | 0 | degrees | degrees |
| --enable-wind-shear | Wind shear with altitude | false | - | - |
| --sample-trajectory | Sample at regular intervals | false | - | - |
| --sample-interval | Sampling interval (always meters, not unit-system dependent) | 10 | meters | meters |
| --enable-pitch-damping | Transonic stability analysis | false | - | - |
| --enable-precession | Angular motion physics | false | - | - |
| --use-rk4-fixed | Use fixed-step RK4 instead of adaptive RK45 | false | - | - |
| --target-speed | Moving-target speed, 0–300 (see [Mover Ring](#mover-ring---target-speed)); also drives the PDF dope card's Lead column. `0` disables both | 0 | mph | m/s |
| --plot | Inline terminal chart: drop, lateral drift, velocity, then energy vs. range (see [Terminal Chart](#terminal-chart---plot)); bare = braille, `--plot ascii` = ASCII fallback. `-o table` only | off | - | - |

### Manual BC Segments (`--bc-segment`)

A bullet's effective BC changes with its **velocity** (it degrades as the bullet slows,
sharpest through transonic). `--bc-segment VMIN:VMAX:BC` (repeatable) lets you supply your
own velocity-keyed BC ladder — the given BC applies while the bullet's current speed is in
`[VMIN, VMAX)`:

```bash
# BC 0.243 above 1800 fps, 0.228 from 1500-1800, 0.205 from 1200-1500
ballistics trajectory -v 2600 -b 0.243 -m 175 -d 0.308 --drag-model g7 --max-range 1000 \
  --bc-segment 1800:4000:0.243 \
  --bc-segment 1500:1800:0.228 \
  --bc-segment 1200:1500:0.205
```

- **VMIN/VMAX** follow `--units` (fps imperial, m/s metric); **BC** is dimensionless.
- Segments are keyed to **velocity**, not distance — this is orthogonal to `--wind-segment`
  (which is distance-keyed). You can combine both; each applies on its own axis.
- Passing any `--bc-segment` implies `--use-bc-segments` and **overrides** `--bc-table` and
  `--bc-table-dir` (manual pairs are highest priority). An interior gap between segments falls
  back to the manually adjusted base `--bc`; outside the global coverage, the nearest segment
  is used.
- To run BC5D-equivalent corrections on a device that cannot hold the tables (e.g. the
  WASM CLI), run once with `--bc-table-dir ... --use-bc-segments --print-bc-segments`:
  the generated ladder prints as ready-to-paste `--bc-segment` lines (velocities in the
  active `--units`). Pasting the full ladder reproduces the table trajectory to well
  under 1%. Note `--bullet-length` is informational for BC5D: the v2 table axes are
  drag type x weight x BC x muzzle velocity x current velocity — length is not a lookup
  dimension.

### Downrange Wind Segments (`--wind-segment`)

Real wind varies along the bullet's path. `--wind-segment SPEED:ANGLE:UNTIL_DISTANCE[:VERTICAL]`
(repeatable) lets you describe wind that changes with downrange distance — for example a
muzzle reading plus downrange sensor stations. The optional 4th field adds a per-segment
vertical wind component — see [Vertical Wind](#vertical-wind) for the sign convention, its
m/s-regardless-of-`--units` field, and the shear/precedence rules:

```bash
# 8 mph at the muzzle, 12 mph past 300 yd, 18 mph past 600 yd (all from the right)
ballistics trajectory -v 2600 -b 0.243 -m 175 -d 0.308 --max-range 1000 \
  --wind-segment 8:90:300 \
  --wind-segment 12:90:600 \
  --wind-segment 18:90:1000
```

- **SPEED** and **UNTIL_DISTANCE** follow `--units` (mph & yards imperial, m/s & meters
  metric). **ANGLE** is degrees in the wind-FROM convention, same as `--wind-direction`
  (0 = headwind, 90 = from the right, 180 = tailwind, 270 = from the left).
- Each segment applies from the previous boundary out to its `UNTIL_DISTANCE`. The wind
  is a **step function** — there is no interpolation between segments.
- **Wind is zero beyond the last segment.** If your segments don't reach `--max-range`,
  a coverage warning is printed; extend the last segment past the target to avoid it.
- `--wind-segment` **overrides** `--wind-speed`/`--wind-direction` (a note is printed if
  both are given), and is **not compatible with `--enable-wind-shear`**.

### Online Mode Parameters (--online)

When using `--online`, calculations are routed through the Flask API for ML-enhanced predictions:

| Parameter | Description | Default |
|-----------|-------------|---------|
| --online | Route through Flask API | false |
| --api-url | API endpoint URL | https://api.ballistics.7.62x51mm.sh |
| --api-timeout | Request timeout (seconds) | 10 |
| --offline-fallback | Fall back to local if API fails | false |
| --compare | Compare local vs API results | false |
| --enable-weather-zones | Enable weather zone generation | false |
| --enable-3d-weather | Enable altitude weather corrections | false |
| --wind-shear-model | Wind shear model (none/logarithmic/power_law/ekman_spiral) | logarithmic |
| --weather-zone-interpolation | Zone interpolation (linear/cubic/step) | linear |

**Note:** Weather features require `--latitude`, `--longitude`, and `--shot-direction`. Negative values need equals format: `--longitude=-115.2`

`--online` (and `--compare`, which builds the same remote request) routes through the
Flask API, a separate HTTP service that is **G1/G7 only** — a wider drag model warns and
is treated as G1 for that remote request, same as `true-velocity`/`plan-truing`.

### BC5D Table Parameters

BC5D (5-Dimensional BC Correction) tables provide ML-derived corrections for improved accuracy:

| Parameter | Description | Default |
|-----------|-------------|---------|
| --bc-table-dir | Directory with BC5D table files | None |
| --bc-table-auto | Auto-download BC5D tables (online feature) | false |
| --bc-table-url | Base URL for BC5D downloads (online feature) | https://ballistics.tools/downloads/bc5d |
| --bc-table-refresh | Force re-download even if cached (online feature) | false |

**Note:** `--bc-table-auto`, `--bc-table-url`, and `--bc-table-refresh` require the `online` feature. Use `--bc-table-dir` for fully offline operation with pre-downloaded tables.

### True Velocity Command

| Parameter | Description | Default | Imperial | Metric |
|-----------|-------------|---------|----------|--------|
| --measured-drop | Measured drop in MILs | Required | MIL | MIL |
| --range | Range where drop was measured | Required | yards | meters |
| -b, --bc | Ballistic coefficient | Required | - | - |
| --drag-model | Drag model (G1/G7) | g1 | - | - |
| -m, --mass | Bullet mass | Required | grains | grams |
| -d, --diameter | Bullet diameter | Required | inches | mm |
| --chrono-velocity | Chronograph velocity for comparison | None | fps | m/s |
| --chrono-distance | Screen distance from the muzzle for `--chrono-velocity`; nonzero back-solves true muzzle velocity (requires `--chrono-velocity`; range 1-98 / 0.3-30) | None (no-op) | feet | meters |
| --zero-range | Zero range | 100 | yards | meters |
| --sight-height | Sight height above bore | 2.0 | inches | mm |
| --bullet-length | Bullet length (for BC5D lookup) | auto | inches | mm |
| --temperature | Temperature | 59 | °F | °C |
| --pressure | Barometric pressure | 29.92 | inHg | hPa |
| --humidity | Relative humidity | 50 | % | % |
| --altitude | Altitude | 0 | feet | meters |
| --offline | Force offline mode | false | - | - |
| --offline-fallback | Fall back to local if API fails | false | - | - |
| --bc-table-dir | Directory with BC5D tables | None | - | - |
| --bc-table-auto | Auto-download BC5D tables | false | - | - |

**Note:** The true-velocity command works in both online and offline modes. Use `--offline` for fully local calculation, or omit for API-based calculation (requires `online` feature).


## Practical Examples

### Hunting Zero at 200 Yards
```bash
# Calculate zero
./ballistics zero -v 2650 -b 0.460 -m 180 -d 0.308 --target-distance 200

# Verify with trajectory
./ballistics trajectory -v 2650 -b 0.460 -m 180 -d 0.308 \
  --auto-zero 200 --max-range 400 --full
```

### Long Range Precision
```bash
./ballistics trajectory \
  -v 2850 -b 0.690 -m 230 -d 0.338 \
  --drag-model g7 \
  --use-bc-segments \
  --auto-zero 100 \
  --max-range 1500 \
  --wind-speed 10 \
  --wind-direction 270 \
  --altitude 5000 \
  --full
```

### Load Development Comparison
```bash
# Load 1: Higher velocity
./ballistics monte-carlo -v 2750 -b 0.475 -m 168 -d 0.308 \
  -n 1000 --velocity-std 15 --target-distance 600

# Load 2: More consistent
./ballistics monte-carlo -v 2680 -b 0.475 -m 168 -d 0.308 \
  -n 1000 --velocity-std 8 --target-distance 600
```

### Varmint Trajectory
```bash
./ballistics trajectory \
  -v 3200 -b 0.242 -m 55 -d 0.224 \
  --auto-zero 200 \
  --max-range 500
```

### Wind Shear and Atmospheric Effects
```bash
# Enable wind shear for altitude-dependent wind
./ballistics trajectory \
  -v 2700 -b 0.475 -m 168 -d 0.308 \
  --wind-speed 10 \
  --wind-direction 90 \
  --enable-wind-shear \
  --altitude 5000 \
  --max-range 1000
```

### Trajectory Sampling for Analysis
```bash
# Sample trajectory at 25-yard intervals
./ballistics trajectory \
  -v 2700 -b 0.475 -m 168 -d 0.308 \
  --sample-trajectory \
  --sample-interval 25 \
  --max-range 1000 -o json > sampled_trajectory.json
```

### Transonic Stability Analysis
```bash
# Enable pitch damping for transonic stability warnings
./ballistics trajectory \
  -v 3000 -b 0.475 -m 168 -d 0.308 \
  --enable-pitch-damping \
  --max-range 2000
```

### Precession and Nutation Physics
```bash
# Enable angular motion modeling
./ballistics trajectory \
  -v 2700 -b 0.475 -m 168 -d 0.308 \
  --twist-rate 10 \
  --enable-precession \
  --max-range 1000
```

### Advanced Physics - Magnus and Spin Drift
```bash
# Enable Magnus effect and spin drift for precision calculation
./ballistics trajectory \
  -v 2700 -b 0.475 -m 168 -d 0.308 \
  --twist-rate 10 \
  --twist-right \
  --enable-magnus \
  --enable-spin-drift \
  --wind-speed 10 \
  --wind-direction 90 \
  --max-range 1000

# Left-hand twist barrel (omit --twist-right)
./ballistics trajectory \
  -v 2700 -b 0.475 -m 168 -d 0.308 \
  --twist-rate 12 \
  --enable-magnus \
  --enable-spin-drift \
  --max-range 1000
```

### Coriolis Effect for Extreme Long Range
```bash
# Northern hemisphere shot, eastward
./ballistics trajectory \
  -v 2850 -b 0.690 -m 230 -d 0.338 \
  --enable-coriolis \
  --latitude 45 \
  --shot-direction 90 \
  --max-range 2000

# Complete advanced physics
./ballistics trajectory \
  -v 3000 -b 0.750 -m 250 -d 0.338 \
  --drag-model g7 \
  --twist-rate 8.5 \
  --twist-right \
  --enable-magnus \
  --enable-coriolis \
  --enable-spin-drift \
  --latitude 38.5 \
  --shooting-angle 45 \
  --wind-speed 15 \
  --wind-direction 270 \
  --altitude 6000 \
  --temperature 25 \
  --pressure 25.5 \
  --humidity 30 \
  --max-range 3000
```

**Shot direction matters (Eötvös effect, fixed in 0.21.0):** with `--enable-coriolis`
and `--latitude`, the `--shot-direction` bearing (0=N, 90=E, 180=S, 270=W) changes the
vertical correction. An **east** shot is lifted (`+2Ω·cos(latitude)·v_east`) and prints
slightly higher; a **west** shot is depressed and prints lower; north/south sit in
between. The horizontal (left/right) Coriolis drift is essentially direction-independent
in the northern hemisphere (always to the right). Prior to 0.21.0 `--shot-direction` was
ignored by the local solver, so east and west gave identical output.

### Online Mode with ML Enhancements
```bash
# Basic online calculation
./ballistics trajectory \
  -v 2700 -b 0.475 -m 168 -d 0.308 \
  --auto-zero 100 --max-range 1000 \
  --online

# Online with weather zones and 3D weather
./ballistics trajectory \
  -v 2850 -b 0.690 -m 230 -d 0.338 \
  --drag-model g7 \
  --auto-zero 100 --max-range 2000 \
  --latitude 36.6 --longitude=-115.2 --shot-direction 90 \
  --enable-weather-zones \
  --enable-3d-weather \
  --wind-shear-model logarithmic \
  --online

# Compare local vs API results
./ballistics trajectory \
  -v 2700 -b 0.475 -m 168 -d 0.308 \
  --auto-zero 100 --max-range 1000 \
  --online --compare
```

### Extreme Weather Conditions
```bash
# Cold, low pressure, high humidity (poor conditions)
./ballistics trajectory \
  -v 2700 -b 0.475 -m 168 -d 0.308 \
  --auto-zero 100 \
  --temperature -10 \  # Very cold
  --pressure 28.50 \   # Low pressure storm
  --humidity 95 \      # Near saturation
  --altitude 7000 \    # High altitude
  --max-range 500

# Hot, dry, high pressure (good conditions) 
./ballistics trajectory \
  -v 2700 -b 0.475 -m 168 -d 0.308 \
  --auto-zero 100 \
  --temperature 95 \   # Hot day
  --pressure 30.50 \   # High pressure
  --humidity 10 \      # Very dry
  --altitude 0 \       # Sea level
  --max-range 500
```

## Advanced Features

### Drag Models
`--drag-model` (`-d` on some commands) accepts the full standard-projectile family, every
one backed by its own real Mach-indexed reference table (MBA-1386) — none of them fall
back to another curve:

- **G1**: Standard flat-base projectile (most common)
- **G2**: Aberdeen J projectile
- **G5**: Short 9° boat-tail
- **G6**: Flat-base secant-ogive
- **G7**: Long 7.5° boat-tail (better for long range)
- **G8**: Flat-base 10° secant-ogive
- **GI**: Blunt/flat-nose, flat base
- **GS**: Round-nose sphere
- **RA4**: British RA 1929 reference drag function (McCoy, *Modern Exterior Ballistics*)
- Full drag tables with Mach-indexed coefficients
- Transonic corrections applied automatically
- Standard drag tables are used without an automatic Reynolds multiplier; a low-Re helper remains available through the Rust API only
- An unrecognized `--drag-model` string (a typo, or a family the library doesn't know) still prints a `warning:` to stderr and falls back to G1 for that run.
- `true-velocity` and `plan-truing`'s forward model is deliberately **G1/G7 only** — an
  unsupported family on those two commands warns and coerces to G1, same as before.
- **GL** (a lower-drag long-range family) is explicitly out of scope: its only public
  source is velocity-domain data, which doesn't fit this engine's Mach-indexed table
  format. Not planned unless a Mach-indexed source turns up.

### BC Modeling
- **BC Segmentation**: Velocity-dependent BC based on bullet type
- **Form Factor**: Additional corrections for bullet shape
- Automatic bullet type identification from parameters

### Physics Engine
- **Integration Methods**:
  - RK45 (Dormand-Prince adaptive) - default for best accuracy
  - RK4 (Runge-Kutta 4th order fixed-step) - available with `--use-rk4-fixed` flag
- Full 3D trajectory integration with six-state modeling
- Magnus effect for spin drift
- Coriolis effect (with latitude input)
- Variable atmospheric conditions
- **Wind Shear**: Altitude-dependent wind profiles
  - Power law model
  - Logarithmic model
  - Exponential decay model
- **Trajectory Sampling**: Regular interval data collection
- **Transonic Effects**:
  - Automatic drag corrections in transonic regime
  - Pitch damping analysis for stability
  - Wave drag modeling
- **Angular Motion**:
  - Precession physics
  - Nutation modeling
  - Gyroscopic stability calculations
- Ground impact detection

#### Advanced Physics Notes
- **Spin Drift**: Requires `--enable-magnus` or `--enable-coriolis` plus `--enable-spin-drift`
- **Magnus Effect**: Side force from spinning projectile, requires `--twist-rate` specification
- **Coriolis Effect**: Earth rotation effects, requires `--latitude` and `--shooting-angle`
- **Twist Direction**: Use `--twist-right` for right-hand twist, omit for left-hand twist
- **Wind Shear**: Models wind speed increase with altitude, affects long-range shots
- **Trajectory Sampling**: Use with JSON/CSV output for detailed analysis
- **Pitch Damping**: Warns about transonic instability (Mach 0.8-1.2)
- **Precession/Nutation**: Models angular motion of spinning projectiles
- **Integration Method**: RK45 adaptive is default (most accurate), RK4 fixed-step available for speed
- Both Magnus and spin drift work together to model the complete gyroscopic effects

### Atmospheric Modeling
- **Temperature Effects**: Affects air density and speed of sound
- **Pressure Effects**: Direct impact on air density (drag)
- **Humidity Effects**: 
  - Humid air is less dense (reduces drag)
  - Increases speed of sound slightly
  - Uses Arden Buck equations for vapor pressure
- **Altitude Effects**: Automatic pressure/density reduction with elevation
- **ICAO Standard Atmosphere**: Full implementation up to 84km
- **CIPM Formula**: Precise air density calculations with humidity

## Online Reverse Solvers

These subcommands call the hosted service at `https://api.ballistics.7.62x51mm.sh` (override
with `--api-url`) instead of computing locally. They require a CLI access token — create one at
`ballisticsinsight.com/account` and save it with `login`. Every other subcommand is fully local
and needs no token or network access.

```bash
# Save your CLI access token (or run `login` with no --token to be prompted)
ballistics login --token bpat_XXXXXXXX
# ...and remove it later
ballistics logout
```

The token is read from the `BALLISTICS_API_TOKEN` environment variable if set, otherwise from
`~/.ballistics/credentials.toml` (written with `0600` permissions).

```bash
# Powder/charge recommendations for a load
ballistics recommend-powder --cartridge "308 Winchester" --bullet-weight 175 \
    --desired-velocity 2600            # optional: --barrel-length 24 --velocity-tolerance 25

# Optimal barrel twist rate
ballistics recommend-twist --caliber 0.308 --weight 175 --overall-length 1.24 --nose-length 0.65

# Recommended cartridge overall length
ballistics recommend-col --cartridge "308 Winchester"   # optional: --bullet-weight 175 --bullet-type SST

# Estimate a bullet's BC from its dimensions
ballistics calibrate-bc --diameter 0.308 --length 1.24 --weight 175   # optional: --bullet-type "boat tail"
```

If the service asks you to sign in, the CLI prints an actionable hint. These subcommands are
available whenever the `online` feature is enabled (on by default).

## Notes

- Default units are Imperial (fps, grains, yards)
- All internal calculations use SI units for precision
- BC values are dimensionless (same for G1 and G7)
- Wind direction: 0° = headwind, 90° = from right, 180° = tailwind, 270° = from left
- Trajectory stops at ground impact or max range
- Sight height default is 1.8 inches (0.05 yards) above bore
- Bore height default is 5 feet (1.5 meters) above ground - adjust for shooting position (e.g., 2ft prone, 4ft sitting, 5ft standing)
