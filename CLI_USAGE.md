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
- Temperature: Fahrenheit
- Pressure: inHg

### Metric Units
- Velocity: meters per second (m/s)
- Mass: grams
- Distance: meters
- Diameter: millimeters
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

## Output Formats

All commands support three output formats via `-o`:

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
| --sight-height | Sight height above bore | 0.05 | yards | meters |
| --max-range | Maximum range | 1000 | yards | meters |
| --time-step | Integration time step | 0.001 | seconds | seconds |
| --wind-speed | Wind speed | 0 | mph | m/s |
| --wind-direction | Wind direction | 0° | degrees | degrees |
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
| --latitude | Latitude for Coriolis | None | degrees | degrees |
| --shooting-angle | Azimuth angle | 0 | degrees | degrees |
| --enable-wind-shear | Wind shear with altitude | false | - | - |
| --sample-trajectory | Sample at regular intervals | false | - | - |
| --sample-interval | Sampling interval | 10 | yards/meters | yards/meters |
| --enable-pitch-damping | Transonic stability analysis | false | - | - |
| --enable-precession | Angular motion physics | false | - | - |
| --use-euler | Use Euler integration | false | - | - |


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
  --shooting-angle 90 \
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
  - RK4 (Runge-Kutta 4th order) - default for accuracy
  - Euler method - available with `--use-euler` flag
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
- **Integration Method**: RK4 is default (more accurate), Euler available for speed
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