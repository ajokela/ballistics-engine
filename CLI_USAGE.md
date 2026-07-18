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

### Custom Drag Tables

Supply a measured or manufacturer-published drag curve — Hornady CDM data, a Lapua/Doppler-radar-derived deck, or your own — instead of relying on a G1/G7 reference curve plus a single BC value. Available via `--drag-table <FILE>` on the `trajectory`, `zero`, and `monte-carlo` subcommands.

**CSV format:** two columns, `mach,cd`, one point per line.
- A single leading header row (e.g. `mach,cd`) is tolerated and skipped.
- Blank lines and lines starting with `#` are ignored.
- Mach must be strictly ascending, with at least 2 data points.
- Cd must be finite and greater than 0.
- **Mach-keyed only.** Velocity-keyed decks (e.g. raw Doppler output in fps/m/s) must be converted by you first: `mach = velocity / speed_of_sound` at the conditions the velocity was measured under.

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

**WASM:** the browser/Node build has no filesystem, so `--drag-table <FILE>` isn't a CLI flag there. Instead, call `wasm.loadDragTable(bytes)` with the raw bytes of the same `mach,cd` CSV (parsed by the identical `DragTable::from_csv_str`) before running a command; `wasm.hasDragTable()` reports whether one is loaded. Once loaded, the table is applied automatically to every `trajectory`, `zero`, `lead`, and `monte-carlo` run — including `lead`, which has no native `--drag-table` flag of its own — until a new table replaces it. See `loadDragTable`'s doc comment in `src/wasm.rs` for the exact contract.

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

The `.a7p` wire format is implemented independently for interoperability;
no third-party code is bundled.

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
- **`--adjustment-unit <mil|moa>`** — angular unit for the `Lead (MIL/MOA)` column
  (default `mil`).
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
- **`--adjustment-unit <mil|moa>`** — angular unit for the ring **table** column only
  (default `mil`; the flag trajectory already exposes for the PDF dope card). With `moa`
  the column reads `Ring(moa)` with `ring_moa = ring_mil × 3.438` — the CLI's locked
  printed-table dial convention (MBA-724, deliberately not the exact-angle 3437.7467/1000),
  so Ring keeps the same MIL/MOA ratio as every other hold column. CSV keeps `ring_mil` and
  JSON keeps `mover_ring_m`/`mover_ring_mil` regardless — their names are the unit contract.

**Table** (`--full -o table`) gains a `Ring(mil)` column (`Ring(moa)` under
`--adjustment-unit moa`). The muzzle point prints `-` (no flight time has elapsed, so the
ring has no defined angular size there yet):

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

Render two stacked inline terminal charts after the normal `trajectory` output (MBA-1320):
drop vs. range, then lateral drift vs. range. Pure Rust, zero new dependencies — no
terminal-graphics crate, no terminal-size detection, no ANSI colors.

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
- Both charts plot the SAME per-point data the `--full` "Trajectory Points:" table prints
  (`result.points`, the raw un-decimated integration output — `--plot` works without
  `--full` too, it just doesn't print the table itself). Drop is the table's `Y` column
  (vertical position), lateral drift is the table's downrange-paired `Z` column (not
  printed by the table by default) — both in the same range unit (yd/m) the rest of the
  table uses, **not** inches. This is a different, deliberate convention from
  `--sample-trajectory`'s sight-line-relative `drop_m`/`wind_drift_m` (see
  [Trajectory Sampling for Analysis](#trajectory-sampling-for-analysis)); don't conflate
  the two.

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
- **`--adjustment-unit <mil|moa>`** — angular unit for the drift columns (default `mil`).
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
.308 in / 7.82 mm). `DRAG` is `g1` or `g7`, and `NAME` may not contain `:`. Between 2 and
8 loads are accepted, from `--load` and/or `--profile` in any combination. A saved
profile's velocity-BC segments and custom Cd(Mach) drag curve (e.g. from an `.a7p`
import) ARE consumed here — they drive both the load's zeroing and its trajectory runs,
and such loads are tagged `[BC segments]` / `[custom drag curve]` in the table legend
(and flagged in JSON). Inline `--load` specs use the scalar BC.

The table shows per-load drop, drift (both in the `--adjustment-unit`), and velocity at
each range. `-o json` adds linear drop/drift, energy, time of flight, each load's zero
angle, and per-row deltas against load #1 (`delta_drop`, `delta_drift`, `delta_velocity`,
`delta_energy` — zero for the baseline itself); `-o csv` emits one column group per load
(names sanitized for CSV). PDF output is not supported for this command.

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

Past a range the (undispersed) trajectory can no longer physically reach — because a nearly
flat shot from a normal bore height eventually crosses the ground plane — `Dominant` and the
share columns read `n/a`: the variance attribution isn't meaningful there, though `P(hit)` is
still correct (and will simply read 0%, since every sample is a miss).

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

This is computed **analytically** from a linearized (one-standard-deviation finite-difference)
sensitivity of the impact point to each independent source, rather than by re-running the full
Monte Carlo sample set once per bucket with that source zeroed out. A full decomposed re-run
would multiply the sweep's cost by the number of buckets; the linearized estimate instead costs
a handful of extra deterministic trajectory solves per range and keeps the default sweep's
runtime in the same ballpark as the base `monte-carlo` command (a few seconds for the default
9-step sweep). At the magnitude of a single sigma the ballistic response is close enough to
linear for this to be a reasonable first-order error budget — not an exact decomposition.

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
```

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

> With zero `--observed` flags the command behaves exactly as the classic
> single-observation velocity truing described above.

> The WASM terminal (ballistics.sh) supports the `true-velocity` command —
> single- and multi-observation — with output matching the native CLI. The
> gaps: the online/BC5D flags (`--bc-table-dir`, `--bc-table-auto`,
> `--bc-table-url`, `--offline-fallback`, `--api-url`, `--api-timeout`) and
> `--bullet-length` are not exposed there (`--offline` is accepted as a no-op;
> the terminal always calculates locally).

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
- The **mover ring** renders as an angular hold (`--adjustment-unit`, mil or MOA)
  in the human table, while CSV keeps `ring_mil` and JSON keeps `mover_ring_m` /
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
| `--adjustment-unit` | Angular unit for Drop/Wind/Lead columns: `mil` (default) or `moa` |
| `--target-speed` | Target speed for the Lead column (mph imperial / m/s metric — MBA-1325 fixed this to follow `--units` like every other speed flag; previously always read as mph). Also enables the [Mover Ring](#mover-ring---target-speed) column/fields on every output format |
| `--powder` | Powder type (shown in footer) |
| `--bullet-name` | Bullet name (shown in footer) |
| `--location-name` | Location name (shown in header) |
| `--profile-row` | Rifle name (shown in header) |
| `--font-scale` / `--font-preset` | Data-table font size |
| `--bold-data` | Bold font for data cells |

**PDF features:**
- Two-column table layout with Range (yd) and Drop/Wind/Lead in **MIL or MOA** (via `--adjustment-unit`)
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
| -m, --mass | Projectile mass | Required | grains | grams |
| -d, --diameter | Projectile diameter | Required | inches | mm |
| --drag-model | Drag model (g1/g7) | g1 | - | - |
| --auto-zero | Auto-zero distance | None | yards | meters |
| --zero-velocity | Zero-day muzzle velocity (auto-zero only); overrides both powder models | shot-day velocity | fps | m/s |
| --zero-temperature | Zero-day air temperature (auto-zero only); also resolves linear powder velocity unless `--zero-velocity` is set | shot-day temperature | °F | °C |
| --zero-pressure | Zero-day barometric pressure (auto-zero only) | shot-day pressure | inHg | hPa |
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
| --humidity | Relative humidity | 50 | % | % |
| --altitude | Altitude | 0 | feet | meters |
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
| --plot | Inline terminal chart, drop then lateral drift vs. range (see [Terminal Chart](#terminal-chart---plot)); bare = braille, `--plot ascii` = ASCII fallback. `-o table` only | off | - | - |

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
- **G1**: Standard projectile (most common)
- **G7**: Boat-tail bullets (better for long range)
- Full drag tables with Mach-indexed coefficients
- Transonic corrections applied automatically
- Standard drag tables are used without an automatic Reynolds multiplier; a low-Re helper remains available through the Rust API only

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

## Notes

- Default units are Imperial (fps, grains, yards)
- All internal calculations use SI units for precision
- BC values are dimensionless (same for G1 and G7)
- Wind direction: 0° = headwind, 90° = from right, 180° = tailwind, 270° = from left
- Trajectory stops at ground impact or max range
- Sight height default is 1.8 inches (0.05 yards) above bore
- Bore height default is 5 feet (1.5 meters) above ground - adjust for shooting position (e.g., 2ft prone, 4ft sitting, 5ft standing)
