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

### BC Estimation

Estimate ballistic coefficient from observed trajectory:

```bash
./ballistics estimate-bc \
  -v 2700 -m 168 -d 0.308 \
  --distance1 100 --drop1 0.0 \
  --distance2 200 --drop2 0.023
```

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
| --measured-drop | Measured drop in MILs (positive = below LOS) | Required |
| --range | Range where drop was measured | Required |
| --bc | Ballistic coefficient | Required |
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

## Output Formats

All commands support four output formats via `-o`:

### Table Format (default)
```bash
./ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 -o table
```

### JSON Format
```bash
./ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 -o json > trajectory.json
```

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
| `--target-speed` | Target speed in mph for lead calculation |
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
| --zero-velocity | Zero-day muzzle velocity (auto-zero only) | shot-day velocity | fps | m/s |
| --zero-temperature | Zero-day air temperature (auto-zero only) | shot-day temperature | °F | °C |
| --zero-pressure | Zero-day barometric pressure (auto-zero only) | shot-day pressure | inHg | hPa |
| --zero-humidity | Zero-day relative humidity (auto-zero only) | shot-day humidity | percent | percent |
| --zero-altitude | Zero-day altitude (auto-zero only) | shot-day altitude | feet | meters |
| --powder-temp-curve | Measured `TEMP:VEL,...` powder-temp→velocity table (interpolated at ambient temp, clamped; overrides --powder-temp-sensitivity) | none | °F & fps | °C & m/s |
| --sight-height | Sight height above bore | 0.05 | yards | meters |
| --bore-height | Bore height above ground | 5 | feet | meters |
| --ignore-ground-impact | Disable ground impact detection | false | - | - |
| --max-range | Maximum range | 1000 | yards | meters |
| --time-step | Integration time step — RK4/Euler only (the adaptive RK45 default steps adaptively and ignores this) | 0.001 | seconds | seconds |
| --wind-speed | Wind speed | 0 | mph | m/s |
| --wind-direction | Wind direction (0=headwind, 90=from right, 180=tailwind, 270=from left) | 0° | degrees | degrees |
| --wind-segment | Downrange wind segment `SPEED:ANGLE:UNTIL_DISTANCE` (repeatable) | — | mph & yd | m/s & m |
| --temperature | Temperature | 59 | °F | °C |
| --pressure | Barometric pressure | 29.92 | inHg | hPa |
| --humidity | Relative humidity | 50 | % | % |
| --altitude | Altitude | 0 | feet | meters |
| --use-bc-segments | Enable BC segmentation | false | - | - |
| --full | Show all trajectory points | false | - | - |
| --enable-magnus | Enable Magnus effect | false | - | - |
| --enable-coriolis | Enable Coriolis effect | false | - | - |
| --enable-spin-drift | Enable enhanced spin drift | false | - | - |
| --twist-rate | Barrel twist rate | 12 | inches/turn | inches/turn |
| --twist-right | Right-hand twist | false | - | - |
| --latitude | Latitude for Coriolis/weather | None | degrees | degrees |
| --longitude | Longitude for weather zones | None | degrees | degrees |
| --shot-direction | Shot azimuth (0=N, 90=E) | None | degrees | degrees |
| --shooting-angle | Incline angle (up/down) | 0 | degrees | degrees |
| --enable-wind-shear | Wind shear with altitude | false | - | - |
| --sample-trajectory | Sample at regular intervals | false | - | - |
| --sample-interval | Sampling interval (always meters, not unit-system dependent) | 10 | meters | meters |
| --enable-pitch-damping | Transonic stability analysis | false | - | - |
| --enable-precession | Angular motion physics | false | - | - |
| --use-rk4-fixed | Use fixed-step RK4 instead of adaptive RK45 | false | - | - |

### Downrange Wind Segments (`--wind-segment`)

Real wind varies along the bullet's path. `--wind-segment SPEED:ANGLE:UNTIL_DISTANCE`
(repeatable) lets you describe wind that changes with downrange distance — for example a
muzzle reading plus downrange sensor stations:

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
- Reynolds number corrections for low velocities

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