# Ballistics Engine CLI Tool

A command-line interface for ballistics trajectory calculations.

## Building

```bash
# Build the standalone CLI
rustc --edition 2021 \
  -L target/debug/deps \
  --extern clap=$(ls target/debug/deps/libclap-*.rlib | head -1) \
  --extern serde=$(ls target/debug/deps/libserde-*.rlib | head -1) \
  --extern serde_json=$(ls target/debug/deps/libserde_json-*.rlib | head -1) \
  src/bin/ballistics_cli.rs \
  -o ballistics-cli
```

## Usage

### Basic Commands

```bash
# Show help
./ballistics-cli --help

# Show version and info
./ballistics-cli info
```

### Trajectory Calculation

Calculate ballistic trajectories with various parameters:

```bash
# Basic trajectory with default parameters
./ballistics-cli trajectory -v 800 -a 45

# Full example with all parameters
./ballistics-cli trajectory \
  -v 800        # Initial velocity (m/s)
  -a 45         # Launch angle (degrees)
  -b 0.5        # Ballistic coefficient
  -m 0.02       # Mass (kg)
  -d 0.00762    # Diameter (meters)
  --max-time 10 # Maximum simulation time (seconds)
  --time-step 0.01 # Integration time step (seconds)
  -o table      # Output format (table/json/csv)
```

### Monte Carlo Simulation

Run statistical analysis with parameter variations:

```bash
# Basic Monte Carlo simulation
./ballistics-cli monte-carlo -v 800 -a 45 -n 1000

# With custom standard deviations
./ballistics-cli monte-carlo \
  -v 800        # Base velocity (m/s)
  -a 45         # Base angle (degrees)
  -n 1000       # Number of simulations
  --velocity-std 10   # Velocity std dev (m/s)
  --angle-std 1.0     # Angle std dev (degrees)
  --bc-std 0.05       # BC std dev
  -o summary    # Output format
```

#### Monte Carlo Output Formats

- **summary** (default): Formatted statistics table
- **full**: Complete JSON output with all statistics
- **statistics**: CSV format for data analysis

```bash
# Summary output
./ballistics-cli monte-carlo -v 700 -a 30 -n 500 -o summary

# JSON output
./ballistics-cli monte-carlo -v 700 -a 30 -n 500 -o full > monte_carlo.json

# CSV statistics
./ballistics-cli monte-carlo -v 700 -a 30 -n 500 -o statistics > stats.csv
```

### Output Formats

#### Table Format (default)
Displays results in a formatted ASCII table:
```bash
./ballistics-cli trajectory -v 800 -a 30 -o table
```

#### JSON Format
Outputs complete trajectory data as JSON:
```bash
./ballistics-cli trajectory -v 800 -a 30 -o json > trajectory.json
```

#### CSV Format
Outputs trajectory points as CSV for data analysis:
```bash
./ballistics-cli trajectory -v 800 -a 30 -o csv > trajectory.csv
```

## Examples

### Long Range Shot
```bash
./ballistics-cli trajectory -v 900 -a 35 -m 0.045 -b 0.7 --max-time 20
```

### High Angle Shot
```bash
./ballistics-cli trajectory -v 500 -a 70 -m 0.01 --time-step 0.001
```

### Monte Carlo Analysis
```bash
# Analyze precision at long range
./ballistics-cli monte-carlo -v 850 -a 35 -n 1000 --velocity-std 5 --angle-std 0.2

# High-precision system analysis
./ballistics-cli monte-carlo -v 900 -a 30 -n 2000 --velocity-std 2 --angle-std 0.1 --bc-std 0.01
```

### Export for Analysis
```bash
# Generate CSV data for plotting
./ballistics-cli trajectory -v 750 -a 45 --time-step 0.1 -o csv > data.csv

# Generate JSON for programmatic processing
./ballistics-cli trajectory -v 750 -a 45 -o json | jq '.max_range'
```

## Parameters

### Trajectory Command

| Parameter | Short | Description | Default | Unit |
|-----------|-------|-------------|---------|------|
| velocity | -v | Initial muzzle velocity | Required | m/s |
| angle | -a | Launch angle | 0.0 | degrees |
| bc | -b | Ballistic coefficient | 0.5 | - |
| mass | -m | Projectile mass | 0.01 | kg |
| diameter | -d | Projectile diameter | 0.00762 | meters |
| max-time | - | Maximum simulation time | 10.0 | seconds |
| time-step | - | Integration time step | 0.01 | seconds |
| output | -o | Output format (table/json/csv) | table | - |

### Monte Carlo Command

| Parameter | Short | Description | Default | Unit |
|-----------|-------|-------------|---------|------|
| velocity | -v | Base velocity | Required | m/s |
| angle | -a | Base launch angle | 0.0 | degrees |
| bc | -b | Base ballistic coefficient | 0.5 | - |
| mass | -m | Projectile mass | 0.01 | kg |
| diameter | -d | Projectile diameter | 0.00762 | meters |
| num-sims | -n | Number of simulations | 1000 | - |
| velocity-std | - | Velocity standard deviation | 5.0 | m/s |
| angle-std | - | Angle standard deviation | 0.5 | degrees |
| bc-std | - | BC standard deviation | 0.02 | - |
| output | -o | Output format (summary/full/statistics) | summary | - |

## Notes

- The current implementation uses a simplified drag model
- Atmospheric conditions are fixed at sea level standard atmosphere
- The trajectory stops when the projectile hits the ground (y < 0)
- All calculations use SI units internally