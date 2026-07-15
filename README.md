# Ballistics Engine

A high-performance ballistics trajectory calculation engine with comprehensive physics modeling, automatic zeroing, and statistical analysis capabilities.

**Project Website:** [https://ballistics.rs/](https://ballistics.rs/)

## Features

- **Full 3D Trajectory Integration** - Six-state ballistic modeling with adaptive RK45 and fixed-step RK4 integration methods
- **Advanced Drag Models** - Support for G1 and G7 reference curves (with automatic transonic corrections) plus user-supplied custom Cd(Mach) drag tables (`--drag-table`, used as-is with endpoint hold outside their measured domain, no transonic correction applied — see [CLI_USAGE.md](CLI_USAGE.md#custom-drag-tables); `bc_value` is ignored while a custom table is active)
- **Automatic Zeroing** - Calculate sight adjustments and apply zero angles automatically
- **Canted-Rifle Modeling** - Model a rifle zeroed level but fired canted (`--cant <DEGREES>`, alias `--cant-angle`, on `trajectory`/`monte-carlo`); clockwise cant shifts point of impact right and low downrange for a rifle with an upward zero correction — see [CLI_USAGE.md](CLI_USAGE.md#canted-shooting)
- **Moving-Target Lead** - Wind-aware hold tables for targets moving at a constant speed/angle, with iterative intercept-range correction for non-perpendicular motion (`lead` subcommand; public `ballistics_engine::calculate_lead` API) — see [CLI_USAGE.md](CLI_USAGE.md#moving-target-lead)
- **Mover Ring** - Field-tested alternative for engaging movers: a per-point ring radius (`target_speed × time-of-flight`) falls out of an already-solved trajectory with no second command or re-entered ballistic data (`trajectory --target-speed`, additive across table/JSON/CSV output); `lead` also gained `trajectory`'s powder-temperature flags for muzzle-velocity parity between the two — see [CLI_USAGE.md](CLI_USAGE.md#mover-ring---target-speed)
- **Unit Conversion** - Seamless switching between Imperial (default) and Metric units
- **BC Segmentation** - Velocity-dependent ballistic coefficient modeling with automatic estimation
- **Atmospheric Modeling** - Temperature, pressure, humidity, and altitude effects with ICAO standard atmosphere
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

# With BC5D tables for improved accuracy
./ballistics true-velocity \
  --measured-drop 5.1 --range 600 \
  --bc 0.27 --drag-model g7 \
  --mass 140 --diameter 0.264 \
  --bc-table-auto --offline
```

Use case: A shooter measures 5.1 MIL of drop at 600 yards. Their chronograph showed 2822 fps. The command calculates the effective velocity is actually ~2740 fps, suggesting a -82 fps adjustment for accurate ballistic predictions.

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
- **Drag Modeling** - Multiple drag functions (G1-G8, JBM, custom curves) with transonic flow corrections
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
