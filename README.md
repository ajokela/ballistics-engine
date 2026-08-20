# Ballistics Engine

A high-performance ballistics trajectory calculation engine with comprehensive physics modeling, automatic zeroing, and statistical analysis capabilities.

**Project Website:** [https://ballistics.rs/](https://ballistics.rs/)

## Features

- **Full 3D Trajectory Integration** - Six-state ballistic modeling with adaptive RK45 and fixed-step RK4 integration methods
- **Advanced Drag Models** - Full standard-projectile family (G1, G2, G5, G6, G7, G8, GI, GS, and the British RA4 reference function), each backed by its own real Mach-indexed table with automatic transonic corrections, plus user-supplied custom Cd(Mach) drag tables (`--drag-table`, used as-is with endpoint hold outside their measured domain, no transonic correction applied — see [CLI_USAGE.md](CLI_USAGE.md#custom-drag-tables); `bc_value` is ignored while a custom table is active), with an optional `--cd-scale <FACTOR>` whole-curve truing multiplier (Hornady AFF / AB CDF style; `1.0` = neutral, typical range 0.90-1.10; requires `--drag-table` — see [CLI_USAGE.md](CLI_USAGE.md#whole-curve-drag-scale---cd-scale)). `true-velocity`/`plan-truing`'s forward model is deliberately G1/G7 only. GL is out of scope — its only public source is velocity-domain data, which doesn't fit this engine's Mach-indexed table format.
- **Automatic Zeroing** - Calculate sight adjustments and apply zero angles automatically; `trajectory` echoes the solved bore angle (degrees, additive across table/JSON/CSV) whenever auto-zero ran, and `zero --from-angle <DEGREES>` solves the zero RANGE(S) a previously solved/stored bore angle produces — a bore angle generally implies two zeros (the classic 25/300-yard relationship), so both are reported rather than one being silently picked — so the angle can be captured once and reused later, independent of the day it was solved — see [CLI_USAGE.md](CLI_USAGE.md#solving-range-from-a-stored-angle---from-angle)
- **Canted-Rifle Modeling** - Model a rifle zeroed level but fired canted (`--cant <DEGREES>`, alias `--cant-angle`, on `trajectory`/`monte-carlo`); clockwise cant shifts point of impact right and low downrange for a rifle with an upward zero correction — see [CLI_USAGE.md](CLI_USAGE.md#canted-shooting)
- **Deliberate Zero POI Offset** - Record that the rifle is deliberately zeroed off — e.g. 0.1 in high / 0.2 in left at the zero range (Kestrel ZH/ZO semantics) — and shift the whole solution by the equivalent angular bias (`--zero-poi-up`/`--zero-poi-right`, inches imperial / cm metric, on `trajectory` and every zero-solving subcommand; solve-json `shot.zero_poi_up_m`/`zero_poi_right_m`; saved-profile fields; `.a7p` zero click counts convertible via `profile import --zero-click`) — see [CLI_USAGE.md](CLI_USAGE.md#sight-geometry-and-zero-state)
- **Lateral Sight Offset** - Model offset-mounted optics (`--sight-offset`, inches imperial / mm metric, positive = sight right of bore): the bullet starts that far left of the sight line and the windage zero converges it onto the sight line at the zero range — physical mount geometry, distinct from and additive with the zero POI offset (solve-json `rifle.sight_offset_lateral_m`; saved-profile field; also on `monte-carlo`/WEZ mirroring `--cant`) — see [CLI_USAGE.md](CLI_USAGE.md#sight-geometry-and-zero-state)
- **Multiple Named Zeroes / Per-Load Offsets in Profiles** - Store alternate zero conditions on a saved profile (Lapua Sight-In POI / ATrag zero-zone class): each named set carries an optional zero distance plus constant per-load dial corrections in mils (`profile zero-set add|remove|list`; `--zero-set NAME` on `trajectory`/`come-ups`/`wind-card`/`range-table`/`dsf`/`plan-truing`; corrections are added to total-correction dial outputs before the tracking-CF division; profile-CSV `V_OFFSET_MIL`/`H_OFFSET_MIL` columns and `.a7p --zero-click` import feed sets automatically; unknown names fail loudly listing the available sets) — see [CLI_USAGE.md](CLI_USAGE.md#named-zero-sets-and-per-load-offsets-profile-zero-set---zero-set)
- **Drops Reference: LOS vs Target Plane** - Reference sampled drops to the line of sight (default, unchanged) or to the vertical target plane on steep inclined shots (`trajectory --drops-reference target`: drop ÷ cos(shooting angle), JBM's "target plane" checkbox; sampled table/CSV column relabeled `Drop (target)`/`drop_target_in`; solve-json `shot.drops_reference`) — see [CLI_USAGE.md](CLI_USAGE.md#drops-reference-los-vs-target-plane---drops-reference--mba-1403)
- **Scope Tracking Correction Factors** - Compensate turrets that don't track their nominal click value: derive `CF = actual/dialed` from a tall-target test (`tall-target` subcommand, pure arithmetic), then `--elevation-cf`/`--windage-cf` (or validated saved-profile fields) DIVIDE every dial-unit output once at the shared conversion boundary — come-ups, range/compare/wind cards, mover lead/Ring, PDF dope card, zero MOA/mrad — never raw inches (an under-tracking scope, CF < 1, needs more dial); `true-velocity` multiplies dialed observations by the CF (scope-dial to true angular) so scope error is not baked into trued MV/BC — see [CLI_USAGE.md](CLI_USAGE.md#scope-tracking-correction-factors---elevation-cf---windage-cf-tall-target--mba-1358)
- **Equivalent Horizontal Range (BDC shoot-to)** - Inclined zeroed shots print the flat range whose angular correction against the same zero matches the inclined solution (SIG AMR / Leica EHR / Gunwerks style, angular-match inversion over one flat re-solve — not the rifleman's-rule cosine), so fixed BDC turrets/reticles can dial as if flat (`trajectory` summary line; solve-json `summary.equivalent_horizontal_range_m`; public `TrajectorySolver::equivalent_horizontal_range`) — see [CLI_USAGE.md](CLI_USAGE.md#equivalent-horizontal-range-bdc-shoot-to--mba-1395)
- **Reticle Hold Points** - Place a firing solution where you actually read it: a point in your own reticle, FFP/SFP aware with second-focal-plane subtensions rescaled for the magnification in use (`reticle hold`, native + browser terminal + a new appended C ABI export; `reticle generate mil-grid|tree|bdc` builds descriptions in one shared serde schema that saved profiles, `--reticle-json` and solve-json's optional `reticle` block all speak). Horus/TREMOR grid layouts, wind-dot calibration and any vendor catalog are deliberately excluded — see [CLI_USAGE.md](CLI_USAGE.md#reticle-hold-points-reticle--mba-1361)
- **Reticle/BDC Inverse Solvers** - Three read-only solvers over an existing load, sharing one drop-vs-range root find: `mark-to-range` maps each reticle subtension to the range where it lands (Nightforce / Nikon Spot On / Swarovski / TRACT), reporting marks the load cannot reach rather than dropping them; `bdc-match` fits the magnification that makes an SFP BDC reticle match the load (Zeiss Rapid-Z; closed-form least squares, with a residual warning when nothing fits); `optimal-zero` min-max searches the one zero that minimizes the largest hold a whole target list needs and reports whether a dead-center hold keeps each inside its vital zone (GeoBallistics HDZ). CLI-only this train — see [CLI_USAGE.md](CLI_USAGE.md#reticlebdc-inverse-solvers--mba-1362)
- **Robust Hold Corridors** - Solve a bounded set of NAMED segmented-wind scenarios at once (`hold-corridor --scenarios set.json --ranges 200,400,600 [--target rect:WxH|circle:D]`) and get, at every range, each scenario's hold, the min/max corridor they span, the minimax (Chebyshev-center) hold, the worst-case miss from it, and whether one hold keeps every scenario inside the target — with a versioned `RobustHoldReportV1` JSON form. Caps (≤8 scenarios, ≤64 ranges) and malformed segments are structured errors raised *before* any solving, and reordering scenarios cannot change the answer. No probabilities are assigned anywhere: the corridor is the span of the hypotheses you supplied, not a confidence interval — see [CLI_USAGE.md](CLI_USAGE.md#robust-hold-corridors-hold-corridor--mba-1349)
- **Moving-Target Lead** - Wind-aware hold tables for targets moving at a constant speed/angle, with iterative intercept-range correction for non-perpendicular motion (`lead` subcommand; public `ballistics_engine::calculate_lead` API) — see [CLI_USAGE.md](CLI_USAGE.md#moving-target-lead)
- **Mover Ring** - Field-tested alternative for engaging movers: a per-point ring radius (`target_speed × time-of-flight`) falls out of an already-solved trajectory with no second command or re-entered ballistic data (`trajectory --target-speed`, additive across table/JSON/CSV output); `lead` also gained `trajectory`'s powder-temperature flags for muzzle-velocity parity between the two — see [CLI_USAGE.md](CLI_USAGE.md#mover-ring---target-speed)
- **Side-by-Side Load Comparison** - Compare 2-8 loads at identical conditions with per-load independent zeroing (`compare --load "NAME:DRAG:BC:MASS:VELOCITY[:DIAMETER]"`, mixable with saved profiles); JSON/CSV output carries per-row deltas against the first load — see [CLI_USAGE.md](CLI_USAGE.md#load-comparison-compare)
- **Powder Temperature Command** - Resolve the temperature-adjusted muzzle velocity standalone, without a trajectory solve (`powder` subcommand): linear fps-per-degree model or a measured temperature→velocity curve, optional `--sweep` velocity ladder and muzzle energy; shares the solvers' exact resolution code (public `resolve_powder_adjusted_velocity` API) — see [CLI_USAGE.md](CLI_USAGE.md#powder-temperature-velocity-powder)
- **Unit Conversion** - Seamless switching between Imperial (default) and Metric units
- **BC Segmentation** - Velocity-dependent ballistic coefficient modeling with automatic estimation
- **Atmospheric Modeling** - Temperature, pressure, humidity, and altitude effects with ICAO standard atmosphere; also accepts a single **density altitude** reading (`trajectory --density-altitude`, feet imperial / meters metric) as a direct alternative to entering altitude/pressure/temperature separately — back-solves an ISA-equivalent atmosphere (preserving Mach/lapse-rate/segmented-atmosphere behavior, not a density-only shortcut) and supersedes `--altitude`/`--pressure`/`--pressure-type` entirely, with an explicit `--temperature` still honored for correct powder-temperature sensitivity — see [CLI_USAGE.md](CLI_USAGE.md#density-altitude-as-a-direct-input---density-altitude)
- **Clock-Position Wind Entry** - Enter wind direction as the dominant field convention: marked clock positions (`--wind-direction 3oc`, `10h30`, or `10:30`; 12 o'clock = headwind, minutes count 0.5°) alongside plain degrees on every wind-direction flag and the WASM terminal; inside `--wind-segment` the colon-free forms apply (`10:3oc:400`) while `10:30:400` keeps its numeric SPEED:ANGLE:DIST meaning; bare numbers stay degrees everywhere — see [CLI_USAGE.md](CLI_USAGE.md#wind-direction-entry-degrees-clock-positions--mba-1367)
- **Earth-Fixed Compass Wind Bearings** - Store wind as absolute compass bearings and let the solver re-reference them against the shot azimuth (`--wind-ref compass` + `--shot-direction` on `trajectory`/`monte-carlo`; covers the single direction, location-CSV WIND_DIR, and every `--wind-segment` angle; Monte Carlo converts before dispersion sampling; solve-json `wind.wind_reference`; WASM builder `setWindReference`/`setShotDirection`; wind FROM north on a shot due north = pure headwind, pinned) — see [CLI_USAGE.md](CLI_USAGE.md#earth-fixed-compass-bearings---wind-ref-compass--mba-1368)
- **Wind Effects** - 3D wind calculations with altitude-dependent wind shear modeling, **downrange-segmented wind** (`--wind-segment SPEED:ANGLE:DIST[:VERTICAL]`, repeatable — model wind that varies along the path, e.g. muzzle plus downrange sensor readings), and **vertical wind** (`--wind-vertical <SPEED>` on `trajectory`/`monte-carlo`, or the segment's optional 4th field; positive = updraft, raises point of impact) — see [CLI_USAGE.md](CLI_USAGE.md#vertical-wind)
- **Oblique Wind-Drift Cards** - Wind dope cards at any wind-FROM angle, not just full-value 90° crosswind (`wind-card --wind-angle <DEG>` or `--wind-angles <CSV>` for one card per angle); each cell is a real trajectory solve, default (no flags) unchanged from the classic full-value 90° card — see [CLI_USAGE.md](CLI_USAGE.md#wind-card)
- **Monte Carlo Simulations** - Statistical analysis with parameter uncertainties
- **BC Estimation** - Estimate ballistic coefficients from trajectory data
- **Advanced Physics**:
  - **Spin Effects**: Magnus effect and empirical Litz spin drift
  - **Earth Effects**: Coriolis effect with latitude-dependent calculations
  - **Angular Motion**: Gyroscopic precession and nutation physics
  - **Transonic Analysis**: Pitch damping coefficients and stability warnings
  - **Trajectory Sampling**: Regular interval data collection for analysis
  - **Form Factor Corrections**: Bullet-specific drag adjustments
- **Multiple Output Formats** - JSON, CSV, formatted tables, and printable PDF dope cards
- **Terminal Chart** - Inline drop, drift, velocity, and energy vs. range charts right in the terminal (`trajectory --plot`, Unicode braille-dot canvas by default, `--plot ascii` fallback); pure Rust, zero new dependencies, no ANSI colors — see [CLI_USAGE.md](CLI_USAGE.md#terminal-chart---plot)
- **Profile import**: `ballistics profile import file.a7p` — imports ArcherBC2 `.a7p` profiles (rifle, bullet, atmosphere, zero) with a full mapping report; `--dry-run` previews without saving
- **Online Reverse Solvers** - Optional `login` + `recommend-powder`/`recommend-twist`/`recommend-col`/`calibrate-bc` subcommands query the hosted service for load, twist, cartridge-overall-length, and BC suggestions using a CLI access token saved from your ballisticsinsight.com account (`BALLISTICS_API_TOKEN` env var, or `~/.ballistics/credentials.toml`); all local subcommands work offline with no token — see [CLI_USAGE.md](CLI_USAGE.md#online-reverse-solvers)
- **Solution Diff Attribution** — explain why two resolved solutions differ, attributed by input group with an explicit interaction remainder (`explain`) — see [CLI_USAGE.md](CLI_USAGE.md#solution-diff-attribution-explain--mba-1345)
- **Per-Input Error Budget** — rank which input is worth measuring better by its share of impact uncertainty, with the hit-probability gain if it were perfected (`error-budget`) — see [CLI_USAGE.md](CLI_USAGE.md#per-input-error-budget-error-budget--mba-1347)
- **Tolerance Envelopes** — how wrong one input may be before the shot leaves the target (`tolerance`) — see [CLI_USAGE.md](CLI_USAGE.md#tolerance-envelopes-tolerance--mba-1350)
- **Constrained Dial & Hold Planning** — rank whole-click dial, reticle-hold, and hybrid execution plans for a TRUE angular correction against a real optic's turret mechanics, travel, and hold bounds, with infeasibility naming the limiting mechanism rather than a silent clamp (`dial-plan`) — see [CLI_USAGE.md](CLI_USAGE.md#constrained-dial--hold-planning-dial-plan--mba-1348)
- **Adaptive Range Cards** — a range card that provably reconstructs the trajectory within a stated elevation/windage error budget, with click rounding from a saved optic and a footer stating the measured worst-case error and the grid it was verified against — a MEASURED error bound and always-present anchors, not a claim of fewer rows than a well-chosen fixed step (`adaptive-card`) — see [CLI_USAGE.md](CLI_USAGE.md#adaptive-range-cards-adaptive-card--mba-1351)
- **Confidence-Controlled Monte Carlo** — every *sampled* hit-probability estimate (the fixed-count run and the opt-in `--adaptive` mode — not the separate `--wez` sweep, which reports a bare point estimate by design) states its sample count, method, confidence level, and interval: an additive Wilson companion line/JSON key on the existing fixed-count run, and `--adaptive` sampling in batches until an anytime-valid confidence sequence meets a requested half-width instead of guessing `--num-sims` (`monte-carlo --adaptive`, MBA-1352) — see [CLI_USAGE.md](CLI_USAGE.md#confidence-controlled-sampling---adaptive--mba-1352)

## Installation

### From crates.io

```bash
cargo install ballistics-engine
```

### From Source

```bash
git clone https://github.com/ajokela/ballistics-engine.git
cd ballistics-engine
cargo build --release
```

The binary will be at: `target/release/ballistics`

### Feature Flags

| Feature | Default | Description |
|---------|---------|-------------|
| `online` | ✅ Yes | HTTP client for API integration (`--online` flag) |

To build without network capabilities:
```bash
cargo build --release --no-default-features
```

## Quick Start

### Basic Trajectory (Imperial Units - Default)

```bash
# .308 Winchester, 168gr bullet at 2700 fps
./ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 1000

# With automatic zeroing at 200 yards
./ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --auto-zero 200 --max-range 500
```

### Metric Units

```bash
# Same bullet in metric units
./ballistics trajectory --units metric -v 823 -b 0.475 -m 10.9 -d 7.82 --max-range 1000
```

## Unit Systems

The engine supports two unit systems, selectable with the `--units` flag:

### Imperial (Default)
- **Velocity**: feet per second (fps)
- **Mass**: grains
- **Distance**: yards
- **Diameter**: inches
- **Temperature**: Fahrenheit
- **Pressure**: inHg
- **Wind**: mph

### Metric
- **Velocity**: meters per second (m/s)
- **Mass**: grams
- **Distance**: meters
- **Diameter**: millimeters
- **Temperature**: Celsius
- **Pressure**: hPa (millibars)
- **Wind**: m/s

## Commands

### Trajectory Calculation

Calculate ballistic trajectory with environmental conditions:

```bash
# Imperial units (default)
./ballistics trajectory \
  -v 2700          # Velocity (fps)
  -b 0.475         # Ballistic coefficient
  -m 168           # Mass (grains)
  -d 0.308         # Diameter (inches)
  --drag-model g7  # G7 drag model
  --angle 0        # Launch angle (degrees)
  --max-range 1000 # Maximum range (yards)
  --wind-speed 10  # Wind speed (mph)
  --wind-direction 90  # Wind from right (degrees)
  --temperature 59     # Temperature (Fahrenheit)
  --pressure 29.92    # Pressure (inHg)
  --humidity 50       # Relative humidity (%)
  --altitude 0        # Altitude (feet)
  --full             # Show all trajectory points
```

#### Auto-Zero Feature

Automatically calculate and apply the zero angle for a specific distance:

```bash
# Zero at 200 yards and show trajectory to 500 yards
./ballistics trajectory \
  -v 2700 -b 0.475 -m 168 -d 0.308 \
  --auto-zero 200 \  # Automatically zero at 200 yards
  --max-range 500 \
  --full

# Custom sight height for auto-zero
./ballistics trajectory \
  -v 2700 -b 0.475 -m 168 -d 0.308 \
  --auto-zero 100 \
  --sight-height 0.055  # 2.2 inches in yards
```

#### Zero-Day Conditions (zero shift)

A rifle's zero is a fixed barrel angle set on the day you sighted in. If you later shoot
in different weather — or with a different muzzle velocity (e.g. a cold vs. warm powder
temperature) — the point of impact shifts. By default `--auto-zero` solves the zero angle
using the same conditions you pass for the shot, which assumes you zeroed in today's
conditions. The `--zero-*` flags let you decouple the two: the zero **angle** is solved
under the conditions the rifle was actually zeroed in, while the trajectory itself runs
under the current shot-day conditions.

```bash
# Zeroed on a cold morning (28 F) at 2600 fps; shooting this afternoon at 85 F / 2700 fps.
# The zero angle is solved for the cold/slow load, then the warm/fast trajectory is
# computed against it — so the dope correctly shows the point of impact drifting high.
./ballistics trajectory \
  -v 2700 -b 0.19 -m 77 -d 0.224 --drag-model g7 \
  --temperature 85 --pressure 29.92 \
  --auto-zero 100 --max-range 1000 --full \
  --zero-velocity 2600 \
  --zero-temperature 28
```

Available overrides (each independently optional; any omitted flag falls back to the
shot-day value, so leaving them all off reproduces the previous behavior exactly):

| Flag | Meaning | Units (imperial / metric) |
|------|---------|---------------------------|
| `--zero-velocity` | Muzzle velocity on the zeroing day | fps / m·s⁻¹ |
| `--zero-temperature` | Air temperature on the zeroing day | °F / °C |
| `--zero-pressure` | Barometric pressure on the zeroing day | inHg / hPa |
| `--zero-humidity` | Relative humidity on the zeroing day | percent |
| `--zero-altitude` | Altitude on the zeroing day | feet / meters |

#### Powder Temperature

Propellant temperature changes muzzle velocity. Two models are available:

**Linear** — a constant sensitivity (fps or m/s per degree) applied relative to the
temperature the load was chronographed at:

```bash
./ballistics trajectory -v 2700 -b 0.19 -m 77 -d 0.224 --drag-model g7 \
  --temperature 85 --use-powder-sensitivity \
  --powder-temp-sensitivity 1.2 --powder-temp 70   # +1.2 fps per F above 70 F
```

**Measured curve (non-linear)** — real powders aren't perfectly linear (temperature-
stable powders flatten; others steepen when hot). If you've chronographed the load at
several temperatures, pass the points directly and the muzzle velocity is interpolated
at the powder temperature (clamped at the endpoints — no extrapolation). This
**overrides** `--powder-temp-sensitivity` when supplied:

```bash
./ballistics trajectory -v 2700 -b 0.19 -m 77 -d 0.224 --drag-model g7 \
  --temperature 85 \
  --powder-temp-curve "40:2620,70:2700,100:2760"   # TEMP:VELOCITY points
```

**Powder temperature vs air temperature.** The curve maps *powder* temperature to
velocity, while `--temperature` drives air *density*. These are decoupled: the curve is
looked up at `--powder-temp` when given, otherwise at `--temperature` (powder assumed at
air temperature). So a load left in a hot chamber or a cold pocket:

```bash
# 85 F air (density), but the powder is at 60 F (velocity from the curve at 60 F)
./ballistics trajectory ... --temperature 85 --powder-temp 60 \
  --powder-temp-curve "40:2620,70:2700,100:2760"
```

Both powder models compose with `--auto-zero`, symmetrically. For the linear model,
`--zero-temperature` resolves zero-day velocity relative to the reference `--powder-temp`.
For a curve, `--zero-powder-temp` overrides the powder lookup; otherwise an explicit
`--zero-temperature` is used, or the shot-day `--powder-temp` is inherited when no zero-day
temperature was supplied. Zero-day atmosphere flags still drive air density independently.
An explicit `--zero-velocity` takes precedence over either powder model.

#### Bore Height and Ground Impact

Control bore height above ground and ground impact detection:

```bash
# Set bore height for prone shooting position (2 feet)
./ballistics trajectory \
  -v 2700 -b 0.475 -m 168 -d 0.308 \
  --auto-zero 100 \
  --bore-height 2  # 2 feet (imperial) or meters (metric)

# Disable ground impact detection for full trajectory to max range
./ballistics trajectory \
  -v 2700 -b 0.475 -m 168 -d 0.308 \
  --auto-zero 100 \
  --max-range 1000 \
  --ignore-ground-impact
```

Bore height defaults: 5 feet (imperial) / 1.5 meters (metric) - standing position.

#### Advanced BC Modeling

Enable velocity-dependent BC modeling for more accurate long-range predictions:

```bash
# Enable BC segmentation (velocity-based BC changes)
./ballistics trajectory \
  -v 2700 -b 0.475 -m 168 -d 0.308 \
  --use-bc-segments \
  --auto-zero 600 \
  --max-range 1000
```

#### Advanced Physics - Magnus and Spin Drift

Enable advanced gyroscopic and aerodynamic effects:

```bash
# Magnus effect and spin drift calculation
./ballistics trajectory \
  -v 2700 -b 0.475 -m 168 -d 0.308 \
  --twist-rate 10      # 1:10" barrel twist
  --twist-right        # Right-hand twist
  --enable-magnus      # Enable Magnus effect
  --enable-spin-drift  # Enable empirical Litz spin drift
  --wind-speed 10 \
  --wind-direction 90 \
  --max-range 1000

# Coriolis effect for extreme long range
./ballistics trajectory \
  -v 3000 -b 0.750 -m 250 -d 0.338 \
  --enable-coriolis \
  --latitude 45        # Shooting latitude
  --shooting-angle 90  # Azimuth (0=N, 90=E)
  --max-range 2000
```

### Zero Calculation

Calculate the sight adjustment needed to zero at a specific distance:

```bash
# Calculate zero for 200 yards
./ballistics zero \
  -v 2700 -b 0.475 -m 168 -d 0.308 \
  --target-distance 200

# With custom sight height (default is 0.05 yards / 1.8 inches)
./ballistics zero \
  -v 2700 -b 0.475 -m 168 -d 0.308 \
  --target-distance 300 \
  --sight-height 0.055  # 2.2 inches

# Metric example
./ballistics zero --units metric \
  -v 823 -b 0.475 -m 10.9 -d 7.82 \
  --target-distance 200  # 200 meters
```

Output includes:
- Zero angle in degrees
- Adjustment in MOA (Minutes of Angle)
- Adjustment in mrad (milliradians)
- Maximum ordinate (highest point of trajectory)

### Monte Carlo Simulation

Run statistical analysis with parameter variations:

```bash
./ballistics monte-carlo \
  -v 2700         # Base velocity (fps)
  -b 0.475        # Base BC
  -m 168          # Mass (grains)
  -d 0.308        # Diameter (inches)
  -n 1000         # Number of simulations
  --velocity-std 10    # Velocity std dev (fps)
  --angle-std 0.5     # Angle std dev (degrees)
  --bc-std 0.01       # BC std dev
  --wind-std 2        # Wind speed std dev (mph)
  --wind-direction-std 5  # Wind direction std dev (degrees)
  --target-distance 300  # Target distance for hit probability
```

### BC Estimation

Estimate ballistic coefficient from observed trajectory data:

```bash
./ballistics estimate-bc \
  -v 2700 -m 168 -d 0.308 \
  --distance1 100 --drop1 0.0   # First data point
  --distance2 200 --drop2 0.023  # Second data point
```

### True Velocity (Velocity Truing)

Calculate the effective muzzle velocity that produces a measured drop at a known range. This helps "true" your ballistic system by identifying discrepancies between chronograph readings and real-world performance.

```bash
# Basic offline calculation
./ballistics true-velocity \
  --measured-drop 5.1 --range 600 \
  --bc 0.27 --drag-model g7 \
  --mass 140 --diameter 0.264 \
  --offline

# With chronograph comparison
./ballistics true-velocity \
  --measured-drop 5.1 --range 600 \
  --bc 0.27 --drag-model g7 \
  --mass 140 --diameter 0.264 \
  --chrono-velocity 2822 \
  --offline

# Chronograph measured downrange, not at the muzzle (MBA-1377): most
# chronographs read 10-15 ft (or 25 m) downrange, so --chrono-distance
# back-solves the true muzzle velocity from the raw reading
./ballistics true-velocity \
  --measured-drop 5.1 --range 600 \
  --bc 0.27 --drag-model g7 \
  --mass 140 --diameter 0.264 \
  --chrono-velocity 2822 --chrono-distance 15 \
  --offline

# With BC5D tables for improved accuracy
./ballistics true-velocity \
  --measured-drop 5.1 --range 600 \
  --bc 0.27 --drag-model g7 \
  --mass 140 --diameter 0.264 \
  --bc-table-auto --offline

# Joint MV + BC calibration from multiple observed impacts
./ballistics true-velocity \
  --range 300 --measured-drop 1.30 \
  --observed 600:4.40 --observed 900:9.00 \
  --bc 0.45 --drag-model g1 \
  --mass 168 --diameter 0.308
```

Use case: A shooter measures 5.1 MIL of drop at 600 yards. Their chronograph showed 2822 fps. The command calculates the effective velocity is actually ~2740 fps, suggesting a -82 fps adjustment for accurate ballistic predictions.

**Downrange chronograph correction (MBA-1377).** Most chronographs (and radar units) read some distance downrange rather than at the muzzle — 10-15 ft is typical for optical screens, 25 m for Lapua/JBM's reference distance — so the raw reading is a few fps low. `--chrono-distance` back-solves the true muzzle velocity from that reading via secant iteration on the same forward drag model (BC/`--drag-model`/atmosphere) the rest of the command uses; zero or absent is an exact no-op. It's a pure display-side correction (`--chrono-velocity` never feeds the drop-based solve either way) and validates its input to a sane 1-98 ft / 0.3-30 m band (100 ft is out of range) rather than silently extrapolating a bad distance into a bad velocity.

**Joint velocity + BC truing.** With two or more `--observed RANGE:DROP` impacts spanning supersonic to transonic ranges, `true-velocity` fits *both* muzzle velocity and ballistic coefficient against the real trajectory solver. When the observation set is too short/closely-spaced to separate the two, it refuses the joint fit, trues velocity only, and says so — no false-precision BC. See [CLI_USAGE.md](CLI_USAGE.md#joint-mv--bc-calibration-multiple-observed-impacts) for details.

**Plan the observations before shooting (MBA-1346).** `plan-truing` evaluates a
discrete set of ranges with the same forward model and finite-difference Jacobian
used by the fitter, then chooses an exact-size, minimum-separation-compliant design.
It reports information gain, singular values, conditioning, rejected/unreachable
candidates, and an explicit MV-only recommendation when the available facility
cannot identify BC:

```bash
./ballistics plan-truing \
  -v 2700 -b 0.475 --drag-model g1 -m 168 -d 0.308 \
  --candidate-ranges 200,300,400,500,600,700,800,900 \
  --observation-count 3 --minimum-separation 100 \
  --measurement-resolution 0.03 --drop-unit mil
```

`--measurement-resolution` is the assumed independent **one-standard-deviation**
impact-reading error, not a tolerance or extreme bound. A saved scalar G1/G7
profile may replace the explicit load flags (`--profile NAME`); velocity-banded BC
profiles and custom drag curves are rejected because they do not have one scalar BC
parameter to identify.

**Quantify what the observations actually learned (MBA-1353).** Add
`--observation-sigma` to opt into a weighted joint MV/BC MAP fit and local Gaussian
uncertainty report. A third `--observed RANGE:DROP:SIGMA` field overrides the default
for one reading. Optional priors are always visible and explicit; predictive output
separates uncertainty in the modeled drop from the wider interval for a future
reading:

```bash
./ballistics true-velocity \
  --range 500 --measured-drop 3.18 \
  --observed 600:4.35:0.03 --observed 900:8.89:0.02 \
  --observation-sigma 0.03 \
  --bc 0.45 --drag-model g1 --mass 168 --diameter 0.308 \
  --predict-range 1000 --prediction-sigma 0.03 --output json
```

The report includes MV/BC 95% intervals, covariance and correlation, chi-square,
effective degrees of freedom, prior-domination/weak-identification warnings, and
propagated drop bands. Declared sigmas are treated as absolute known errors, so the
covariance is not rescaled by residual RMS. With no uncertainty flags, the existing
point estimate and output schema are unchanged.

### Wind-Call Truing (`true-wind`)

`true-velocity` trues the vertical axis; `true-wind` (MBA-1392) trues the other one.
Give it where your groups actually landed left/right of the aim point and it back-solves
the crosswind that reproduces that miss through the real forward model, plus a wind-call
correction factor against the wind you *called*:

```bash
./ballistics true-wind \
  --miss 500:14.0 --miss 700:29.5 \
  -v 2700 -b 0.475 -m 168 -d 0.308 --drag-model g7 \
  --twist-rate 11 --called-wind 9
```

A horizontal miss is not purely wind, so the command separates it: `--twist-rate` is
**required** and gyroscopic spin drift is always modelled and subtracted (a 1:11" .308
drifts ~3.5 in right at 700 yd — read as wind, that alone is several mph of error), and
`--latitude` with `--shot-direction` adds Coriolis. Anything the model had no data for
stays absorbed in the solved wind and is named in the report, so a contaminated number is
never presented as pure wind. Signs are documented in one block: `--miss` positive =
impact **right** of aim, solved wind positive = wind **from the shooter's left**
(9 o'clock) pushing impacts right. `--miss` values are linear inches off the target, not
dial readings, so scope tracking correction factors deliberately do not apply.
See [CLI_USAGE.md](CLI_USAGE.md#wind-call-truing-true-wind--mba-1392).

### DSF (Drop-Scale-Factor) Truing

Second stage of the Applied Ballistics-style two-stage truing workflow (MBA-1357).
Once `true-velocity` has fixed the supersonic (Mach > 1.2) muzzle velocity/BC, the
drop discrepancies that grow through the transonic region and into subsonic flight are
no longer fixable by a single MV correction — the residual is a slowly-varying function
of Mach. `dsf` records observed-drop/predicted-drop ratios at specific Mach <= 1.2
ranges and keys them, one saved profile at a time, to a Mach-indexed table:

```bash
# Stage 1: true the muzzle velocity from a supersonic-range drop reading (as above).
./ballistics true-velocity \
  --measured-drop 3.2 --range 500 \
  --bc 0.475 --mass 168 --diameter 0.308 \
  --offline

# Stage 2: record a subsonic/transonic drop observation on the trued, saved profile.
./ballistics dsf --saved-profile my-rifle --range 900 --observed-drop 5.1mil
```

`dsf` takes no ballistic parameters of its own — it solves the named saved profile's
own trajectory (no CLI overrides) and derives everything else from `--range` and
`--observed-drop` (`mil`, `moa`, or `in`, no separator between number and unit). An
observation whose target-range Mach exceeds 1.2 is rejected outright, pointing back to
`true-velocity`. Up to 6 distinct Mach-keyed points accumulate per profile; a new point
within 0.05 Mach of an existing one supersedes it (announced on stdout); a 7th distinct
point is rejected, naming `--clear-dsf` to make room. `trajectory --saved-profile` and
`come-ups --profile` then auto-apply the table as a **drop-only** correction — velocity,
energy, and time of flight are byte-identical to the untrued solve — printing a
table-output-only note; JSON/CSV carry the corrected drop numbers with no equivalent
text. `profile save NAME ... --clear-dsf` removes an existing table. See
[CLI_USAGE.md](CLI_USAGE.md#dsf-drop-scale-factor-truing) for the full staging-gate
reference.

## Advanced Features

### Online Mode (API Integration)

The CLI can query a remote ballistics API server instead of calculating locally. This enables access to enhanced BC data, ML-augmented predictions, and doppler-derived drag curves.

> **Important:** The `--online` feature connects to a **proprietary cloud service** that is not covered by the MIT license. When using `--online`, trajectory parameters and your IP address are transmitted to our servers. See [ONLINE_SERVICE.md](ONLINE_SERVICE.md) for full terms, privacy policy, and data handling practices.

```bash
# Use online mode to query the API
./ballistics trajectory \
  -v 2700 -b 0.475 -m 168 -d 0.308 \
  --online \
  --max-range 1000

# Custom API endpoint
./ballistics trajectory \
  -v 2700 -b 0.475 -m 168 -d 0.308 \
  --online \
  --api-url https://your-api.example.com/v1/calculate \
  --max-range 1000
```

**Default API**: `https://api.ballistics.7.62x51mm.sh/v1/calculate`

Online mode benefits:
- **Enhanced BC data** - Access to doppler-derived ballistic coefficients
- **ML predictions** - Machine learning augmented trajectory calculations
- **BC segments** - Velocity-dependent BC modeling from measured data
- **Form factor corrections** - Bullet-specific drag adjustments

**Data transmitted when using --online:**
- All trajectory parameters (BC, mass, velocity, wind, atmospheric conditions, etc.)
- Your IP address and client version
- Request logs retained for 30 days, then deleted

To use only local calculations (no network, no data transmission):
```bash
cargo install ballistics-engine --no-default-features
```

### Integration Methods

The engine supports two numerical integration methods:

- **RK45 (Dormand-Prince Adaptive)** - Default method, provides best accuracy with adaptive step sizing
- **RK4 (Runge-Kutta 4th Order Fixed-Step)** - Available with `--use-rk4-fixed` flag for faster computation

### Wind Shear Modeling

Model altitude-dependent wind variations:

```bash
./ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 \
  --wind-speed 10 --wind-direction 90 \
  --enable-wind-shear \
  --max-range 1000
```

### Transonic Stability Analysis

Analyze projectile stability through the transonic regime:

```bash
./ballistics trajectory -v 3000 -b 0.475 -m 168 -d 0.308 \
  --enable-pitch-damping \
  --max-range 2000
```

Provides warnings about transonic instability and minimum pitch damping coefficients.

### Trajectory Sampling

Collect trajectory data at regular intervals for detailed analysis:

```bash
./ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 \
  --sample-trajectory \
  --sample-interval 25  # Sample every 25 meters
  --max-range 1000 -o json
```

### Angular Motion Physics

Model precession and nutation of spinning projectiles:

```bash
./ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 \
  --twist-rate 10 \
  --enable-precession \
  --max-range 1000
```

### Complete Advanced Physics Example

```bash
./ballistics trajectory \
  -v 2850 -b 0.690 -m 230 -d 0.338 \
  --drag-model g7 \
  --twist-rate 8.5 --twist-right \
  --enable-magnus \
  --enable-coriolis \
  --enable-spin-drift \
  --enable-wind-shear \
  --enable-pitch-damping \
  --enable-precession \
  --sample-trajectory \
  --latitude 38.5 \
  --shooting-angle 45 \
  --wind-speed 15 --wind-direction 270 \
  --altitude 6000 \
  --max-range 2000
```

## Physics Modeling

The ballistics engine implements comprehensive physics modeling for accurate trajectory prediction:

### Aerodynamic Effects
- **Drag Modeling** - Full standard-projectile family (G1, G2, G5, G6, G7, G8, GI, GS, RA4) plus custom Cd(Mach) curves, with transonic flow corrections
- **Form Factor** - Projectile efficiency corrections based on shape and design
- **Reynolds Number Effects** - Reynolds diagnostics and an opt-in helper for genuinely low-Re flow; standard drag tables are not multiplied by an extra correction

### Gyroscopic Effects  
- **Spin Drift** - Lateral deviation due to gyroscopic and Magnus effects
- **Precession** - Gyroscopic precession of spinning projectile
- **Nutation** - Oscillatory motion superimposed on precession
- **Spin Decay** - Reduction in spin rate over time due to aerodynamic damping
- **Pitch Damping** - Aerodynamic moments opposing angular motion

### Environmental Effects
- **Coriolis Effect** - Earth's rotation influence on long-range trajectories
- **Magnus Effect** - Force from spinning projectile in crossflow
- **Wind Shear** - Altitude-dependent wind variations
- **Atmospheric Stratification** - Density and sound speed variations with altitude

### Stability Modeling
- **Dynamic Stability** - Gyroscopic and aerodynamic stability calculations
- **Yaw of Repose** - Gravity/gyroscopic equilibrium yaw; crosswind yaw is a transient handled by aerodynamic jump
- **Limit Cycle Yaw** - Bounded oscillatory motion analysis

## Language Bindings

Official language bindings are maintained as separate projects:

- **Python**: [ballistics-engine-py](https://github.com/ajokela/ballistics-engine-py) - PyO3 bindings via maturin
- **Ruby**: [ballistics-engine-rb](https://github.com/ajokela/ballistics-engine-rb) - Magnus bindings via rb_sys

These bindings depend on the `ballistics-engine` crate published on [crates.io](https://crates.io/crates/ballistics-engine).

### WASM / npm Package

The engine also compiles to WebAssembly (`src/wasm.rs`, `wasm-bindgen`) and already powers
[ballistics.sh](https://ballistics.sh) and [ballistics.rs](https://ballistics.rs) in the browser.
It is not yet published to npm for third-party use — `scripts/build-npm.sh` builds and prepares a
publish-ready package; publishing itself is a manual step (see below).

```bash
scripts/build-npm.sh
```

This builds two `wasm-bindgen` targets, both with `--no-default-features` (the default
`pdf`/`online` features pull in `printpdf`/`ureq`+`ring`, which do not compile for
`wasm32-unknown-unknown`) plus
`--features wasm-terminal` (the browser terminal's command set — see
[Trimming the WASM module](#trimming-the-wasm-module)):

- **`pkg/`** — `--target bundler`, the package meant for `npm publish`. Consumed via a native
  `.wasm` ES import by bundlers that understand it (webpack with `experiments.asyncWebAssembly`,
  Vite, Rollup + `@rollup/plugin-wasm`, Parcel).
- **`pkg-web/`** — `--target web`, a no-bundler build for direct `<script type="module">` browser
  use or manual Node usage without a bundler — the same `--target` already used to build
  ballistics.sh/ballistics.rs's WASM. Documented and built for completeness; not published under
  the primary package name in this initial pass.

`wasm-pack` has no built-in dual-target/"publish both" mode, and stitching bundler- and web-target
output into one package.json via manual `exports` conditions isn't something `wasm-pack` generates
or tests for you — see the comment header of `scripts/build-npm.sh` for the full reasoning. A
single bundler-target package as the published npm artifact, with the web build documented
separately, is the ecosystem-standard shape for `wasm-bindgen` crates on npm.

### Trimming the WASM module

The published module carries two independent surfaces: the **`Calculator`** builder API
(`setBC`, `setDragModel`, `setWind`, `enableSpinDrift`, `enableCoriolis`, `calculateTrajectory`,
`getFullTrajectory`, …) and the **browser terminal** (`WasmBallistics.runCommand`) that powers
ballistics.sh. An app that only solves trajectories pays for the terminal's other twelve
commands, which is most of the binary.

Each non-trajectory command sits behind its own cargo feature, so you can select the subset you
actually call. `trajectory`, `version`, and the whole `Calculator` API are **never** gated —
`Calculator` composes a `trajectory` command line internally, so it keeps working with every
feature below turned off.

Build through `scripts/build-wasm.sh`, which is the one entry point every WASM build uses —
the ballistics.rs deploy and `build-npm.sh` included:

```bash
# Everything. Also what you get with no arguments at all — the default is deliberately the
# complete terminal, so a forgotten flag can never silently ship a stripped module.
scripts/build-wasm.sh

# Trajectory only — the Calculator API and nothing else
scripts/build-wasm.sh --preset slim

# À la carte
scripts/build-wasm.sh --features wasm-zero,wasm-lead

# --target and --out-dir pass through; so does the environment
CARGO_PROFILE_RELEASE_OPT_LEVEL=z scripts/build-wasm.sh --target nodejs --out-dir /tmp/pkg
```

After every build the script **verifies the artifact against the preset it was asked for** —
it reads the emitted `.wasm` and checks that exactly the promised commands are present, failing
the build otherwise. `--preset full` expects all twelve regardless of how the feature list was
computed, so a dropped flag is a hard error rather than a terminal that deploys cleanly and
then answers `Unknown command` to everything but `trajectory`.

If you invoke `wasm-pack` directly instead, note the bare `--`: it forwards only post-`--`
arguments to cargo, so `--features` placed before it is consumed as an (invalid) `wasm-pack`
flag.

Measured on 0.33.2, `--target web`, default release profile (`opt-level = 3`, LTO), against the
full build's 918 KB raw / 345 KB gzip (all sizes decimal KB):

| feature | command(s) removed | raw | gzip |
|---|---|---:|---:|
| `wasm-monte-carlo` | `monte-carlo`, including its `--wez` sweep | 115 KB | 46 KB |
| `wasm-truing` | `true-velocity`, `true-wind` | 92 KB | 31 KB |
| `wasm-bc-convert` | `bc-convert` | 65 KB | 22 KB |
| `wasm-reticle` | `reticle` | 55 KB | 21 KB |
| `wasm-lead` | `lead` | 21 KB | 7 KB |
| `wasm-powder` | `powder` | 17 KB | 6 KB |
| `wasm-estimate-bc` | `estimate-bc` | 17 KB | 7 KB |
| `wasm-zero` | `zero` | 15 KB | 4 KB |
| `wasm-recoil` | `recoil` | 12 KB | 3 KB |
| `wasm-power-factor` | `power-factor` | 11 KB | 4 KB |
| `wasm-drag-curve` | `drag-curve` | 7 KB | 3 KB |
| **all of the above** | **`Calculator` + `trajectory` only** | **434 KB** | **153 KB** |

Each row is that feature's marginal cost, measured by dropping it from the full set. The
commands share almost nothing, so the rows are close to additive: they sum to 427 KB raw /
153 KB gzip against a measured all-removed saving of 434 KB / 153 KB — pick any subset and the
rows add up. A trajectory-only module is **483,496 bytes raw, 191,421 gzip, 154,797 brotli**,
against 917,924 / 344,592 / 273,503 for the full build — 44% off the wire.

Splitting the help text into per-command chunks costs the full build about 3 KB raw / 1 KB
gzipped (35 `push_str` calls where there was one literal). That is the price of the table
above; every configuration that drops a command is far ahead.

**The `.wasm` and the JS glue are a matched pair — replace both together.** `wasm-bindgen`
generates the glue to match one specific module, and trimming genuinely changes the module's
import list: dropping `wasm-monte-carlo` removes the last user of `rand`, so the slim `.wasm`
no longer imports `crypto.getRandomValues` and the slim glue no longer supplies it. Ship a
stale full `.wasm` against new slim glue and instantiation fails outright:

```
LinkError: WebAssembly.Instance(): Import #4 "./ballistics_engine_bg.js"
  "__wbg_getRandomValues_..." function import requires a callable
```

The reverse pairing — slim `.wasm` with full glue — happens to load, because the extra import
simply goes unused. Do not rely on that: it is a coincidence of which imports differ today, not
a compatibility guarantee, and it will not hold for a different feature subset. Copy every file
`build-wasm.sh` emits, from the same run, and clear any bundler cache that may hold the old one.

Removing a command does not change any number the remaining ones produce: the full-terminal
build is byte-identical to an ungated build across every command, and `Calculator` output is
byte-identical between the full and trajectory-only builds. A command compiled out reports
`Unknown command`, and the `help` text lists only what is actually present.

Two things are *not* separable, because they are not separate to begin with:

- **`--wez`** is a flag on `monte-carlo`, not a command, so it leaves with `wasm-monte-carlo`.
- **`explain`, `error-budget`, `tolerance`, `dial-plan`, `adaptive-card`** (0.33.x
  decision-support) are native-CLI-only — they were never wired into the WASM dispatch, and
  dead-code elimination already keeps them out of the module. There is nothing to remove.

The script also post-processes each `package.json` (name, description, license, repository,
keywords, and the `files` list — including an `LICENSE-APACHE` entry `wasm-pack` itself omits even
though it copies the file) and installs `README-npm.md` as the package's `README.md`.

**Before publishing**, edit `pkg/package.json`'s `"name"` — it ships as the placeholder
`"@SCOPE/ballistics-engine"`. Replace `SCOPE` with the maintainer's real npm org/user scope (a
scope decision, plus an npm account with publish rights to it, are both needed and don't exist yet
as of this writing). Then:

```bash
scripts/build-npm.sh
cd pkg
npm pack --dry-run   # sanity-check the tarball contents first
npm publish --access public
```

(`--access public` is required the first time a scoped package is published, since scoped packages
default to private on free npm accounts; `pkg/package.json` also sets `publishConfig.access` to
`public` so a plain `npm publish` works too.)

## FFI Layer

The library includes a Foreign Function Interface (FFI) layer for integration with iOS, Android, and other platforms. The FFI provides C-compatible bindings for all major functionality.

<img src="ios.png" alt="iOS Integration Example" width="35%">

### FFI Features
- **C-Compatible Structures** - All data structures use C-compatible layouts
- **Safe Memory Management** - Proper handling of memory across language boundaries
- **iOS/Swift Integration** - Ready for use with Swift through bridging headers
- **Android/JNI Support** - Compatible with Java Native Interface
- **Monte Carlo Simulation** - Statistical analysis with parameter variations
- **Error Handling** - Graceful error propagation across FFI boundary

### Example FFI Usage (C/Swift)
```c
// Create input parameters
FFIBallisticInputs inputs = {
    .muzzle_velocity = 823.0,        // m/s
    .ballistic_coefficient = 0.475,
    .mass = 0.0109,                  // kg
    .diameter = 0.00782,             // meters
    .drag_model = 0,                 // G1
    .sight_height = 0.05,            // meters
    .temperature = 15.0,             // Celsius
    .altitude = 0.0
};

// Calculate trajectory. The final argument is the integration step in milliseconds
// (minimum 0.1 ms; smaller or non-finite values return NULL).
FFITrajectoryResult* result = ballistics_calculate_trajectory(&inputs, NULL, NULL, 1000.0, 0.1);

// NULL also reports invalid inputs or the 250,000-point resource ceiling.
// Increase the step, reduce the range, or use adaptive RK45 for an over-budget solve.
if (result != NULL) {
    printf("Max range: %.2f meters\n", result->max_range);
    ballistics_free_trajectory_result(result);
}
```

### Monte Carlo Simulation via FFI
```c
// Set up Monte Carlo parameters
FFIMonteCarloParams params = {
    .num_simulations = 1000,
    .velocity_std_dev = 10.0,       // m/s variation
    .angle_std_dev = 0.001,         // radian variation (elevation)
    .bc_std_dev = 0.01,             // BC variation
    .wind_speed_std_dev = 2.0,      // m/s wind variation
    .target_distance = 600.0,       // Target at 600m
    .azimuth_std_dev = 0.001        // radian variation (horizontal)
};

// Run simulation with an independent 0.1-radian wind-direction sigma.
// Use ballistics_monte_carlo(...) when no direction variation is desired.
FFIMonteCarloResults* results =
    ballistics_monte_carlo_with_direction_std_dev(&inputs, NULL, &params, 0.1);

// Use statistical results
printf("Mean range: %.2f m (σ=%.2f)\n", results->mean_range, results->std_dev_range);
printf("Hit probability at 600m: %.1f%%\n", results->hit_probability * 100);

// Access individual shots
for (int i = 0; i < results->num_results; i++) {
    printf("Shot %d: Range %.2f m, Impact velocity %.2f m/s\n", 
           i, results->ranges[i], results->impact_velocities[i]);
}

// Clean up
ballistics_free_monte_carlo_results(results);
```

## Output Formats

All commands support three output formats via the `-o` flag:

- **table** (default) - Formatted ASCII table for terminal display
- **json** - Complete data in JSON format for programmatic use
- **csv** - Comma-separated values for spreadsheet analysis

## Practical Examples

### Hunting Zero

Zero a hunting rifle at 200 yards with environmental conditions:

```bash
# Calculate zero angle
./ballistics zero \
  -v 2650 -b 0.460 -m 180 -d 0.308 \
  --target-distance 200

# Verify trajectory with auto-zero
./ballistics trajectory \
  -v 2650 -b 0.460 -m 180 -d 0.308 \
  --auto-zero 200 \
  --max-range 400 \
  --wind-speed 15 \
  --wind-direction 270 \
  --temperature 32 \
  --humidity 30 \
  --altitude 5000 \
  --full
```

### Long Range Shooting

Analyze trajectory for 1000-yard shot:

```bash
./ballistics trajectory \
  -v 2850 -b 0.690 -m 230 -d 0.338 \
  --drag-model g7 \
  --auto-zero 100 \
  --max-range 1100 \
  --wind-speed 10 \
  --wind-direction 45 \
  --full \
  -o json > trajectory.json
```

### Load Development

Compare different loads using Monte Carlo:

```bash
# Load 1: Higher velocity, more variation
./ballistics monte-carlo \
  -v 2750 -b 0.475 -m 168 -d 0.308 \
  -n 1000 \
  --velocity-std 15 \
  --target-distance 600

# Load 2: Lower velocity, more consistent
./ballistics monte-carlo \
  -v 2680 -b 0.475 -m 168 -d 0.308 \
  -n 1000 \
  --velocity-std 8 \
  --target-distance 600
```

## Advanced Features

### BC Segmentation

Velocity-dependent BC modeling accounts for how ballistic coefficient changes as the bullet slows down. Enable with `--use-bc-segments`:

- Automatically estimates BC segments based on bullet characteristics
- No external data required - uses caliber, weight, and BC
- Identifies bullet type (Match, Hunting, VLD, etc.) from parameters
- Applies physics-based BC degradation curves

Example:
```bash
./ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --use-bc-segments --max-range 1000
```

**Manual velocity-keyed BC segments** — supply your own `VMIN:VMAX:BC` pairs (repeatable,
velocities in `--units`) instead of the auto-estimated/table ones. Keyed to velocity, so it
composes with distance-keyed `--wind-segment`; implies `--use-bc-segments` and overrides
`--bc-table` and `--bc-table-dir`:
```bash
./ballistics trajectory -v 2600 -b 0.243 -m 175 -d 0.308 --drag-model g7 --max-range 1000 \
  --bc-segment 1800:4000:0.243 --bc-segment 1500:1800:0.228 --bc-segment 1200:1500:0.205
```

### BC5D Correction Tables

BC5D tables provide ML-derived, 5-dimensional BC corrections indexed by weight, BC, muzzle velocity, current velocity, and drag model. Tables are caliber-specific and capture the complete velocity-dependent behavior.

**Auto-Download Mode** (requires `online` feature):
```bash
# Downloads tables automatically on first use
./ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --bc-table-auto

# Force refresh cached tables
./ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --bc-table-auto --bc-table-refresh
```

**Offline Mode** with pre-downloaded tables:
```bash
./ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --bc-table-dir ./bc_tables/
```

**Available calibers:** .224, .243, .264, .277, .284, .308, .338

**Cache locations:**
- macOS: `~/Library/Caches/ballistics-engine/bc5d/`
- Linux: `~/.cache/ballistics-engine/bc5d/`
- Windows: `%LOCALAPPDATA%\ballistics-engine\cache\bc5d\`

Tables are approximately 1-1.5 MB each and include CRC32 validation to ensure data integrity.

### Advanced Physics Modeling

When enabled, the engine calculates:
- **Magnus Effect** - Side force from spinning projectiles
- **Spin Drift** - Lateral drift due to gyroscopic effects
- **Coriolis Effect** - Earth rotation effects (with latitude input)
- **Transonic Drag** - Enhanced drag modeling in transonic regime
- **Low-Reynolds Helper** - Opt-in viscous correction below the standard projectile-table regime

## Building from Source

### Requirements

- Rust 1.70 or later
- Cargo build system

### Build Commands

```bash
# Debug build
cargo build

# Release build (optimized)
cargo build --release

# Run tests
cargo test

# Build documentation
cargo doc --open
```

## Library Usage

Use as a Rust library in your own projects:

```rust
use ballistics_engine::{
    BallisticInputs, TrajectorySolver, 
    WindConditions, AtmosphericConditions
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let inputs = BallisticInputs {
        muzzle_velocity: 823.0,  // m/s
        launch_angle: 0.0,       // radians
        ballistic_coefficient: 0.475,
        mass: 0.0109,            // kg
        diameter: 0.00782,       // meters
        sight_height: 0.05,      // meters
        ..Default::default()
    };
    
    let wind = WindConditions {
        speed: 5.0,              // m/s
        direction: 1.5708,       // 90 degrees in radians
        ..Default::default()
    };
    
    let atmosphere = AtmosphericConditions {
        temperature: 15.0,       // Celsius
        pressure: 1013.25,       // hPa
        humidity: 50.0,          // %
        altitude: 0.0,           // meters
        ..Default::default()
    };
    
    let solver = TrajectorySolver::new(inputs, wind, atmosphere);
    let result = solver.solve()?;
    
    println!("Max range: {:.2} m", result.max_range);
    println!("Max height: {:.2} m", result.max_height);
    println!("Time of flight: {:.3} s", result.time_of_flight);
    
    Ok(())
}
```

## Performance

Optimized Rust implementation provides:
- Single trajectory (1000m): ~5ms
- Monte Carlo (1000 runs): ~500ms  
- BC estimation: ~50ms
- Zero calculation: ~10ms

## Common Ballistic Coefficients

| Caliber | Weight | BC (G1) | BC (G7) | Description |
|---------|--------|---------|---------|-------------|
| .223    | 55gr   | 0.250   | -       | FMJ |
| .223    | 77gr   | 0.362   | 0.182   | Match |
| .308    | 168gr  | 0.475   | 0.224   | Match |
| .308    | 175gr  | 0.505   | 0.253   | Match |
| .308    | 180gr  | 0.480   | -       | Hunting |
| .338    | 300gr  | 0.768   | 0.383   | Match |
| 6.5mm   | 140gr  | 0.620   | 0.310   | Match |
| .50     | 750gr  | 1.050   | 0.520   | Match |

## Troubleshooting

### Trajectory hits ground early
- Check if you're using `--auto-zero` or setting `--angle` manually
- Default angle is 0° (horizontal), which will hit ground quickly
- Use `--auto-zero <distance>` to automatically calculate proper angle

### Units confusion
- Default is Imperial (fps, grains, yards)
- Use `--units metric` for metric system
- All inputs must match the selected unit system

### Unexpected BC behavior
- G1 and G7 models have different BC values for same bullet
- G7 typically better for boat-tail bullets
- BC segmentation automatically applied based on bullet type

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Add tests for new features
4. Run `cargo test` and `cargo fmt`
5. Submit a pull request

## License

This project is licensed under the MIT License - see LICENSE file for details.

**Note:** The MIT license applies to the open source ballistics-engine library, CLI, and FFI bindings. The `--online` feature connects to a proprietary cloud service with separate terms. See [ONLINE_SERVICE.md](ONLINE_SERVICE.md) for details.

## Acknowledgments

- Ballistics physics based on Robert McCoy's "Modern Exterior Ballistics"
- Drag tables from military ballistics research
- BC segmentation algorithms from Bryan Litz's research
- Community contributions and testing

## Support

For issues, questions, or contributions:
- GitHub Issues: [github.com/ajokela/ballistics-engine/issues](https://github.com/ajokela/ballistics-engine/issues)
