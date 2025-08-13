# Ballistics Engine

A high-performance ballistics trajectory calculation engine with comprehensive physics modeling and statistical analysis capabilities.

## Features

- **6DOF Trajectory Integration** - Full six degrees of freedom ballistic modeling
- **Advanced Drag Models** - Support for G1, G7, and custom drag curves
- **Atmospheric Modeling** - Temperature, pressure, humidity, and altitude effects
- **Wind Effects** - Crosswind, headwind/tailwind drift calculations
- **Monte Carlo Simulations** - Statistical analysis with parameter uncertainties
- **BC Estimation** - Estimate ballistic coefficients from trajectory data
- **Multiple Output Formats** - JSON, CSV, formatted tables
- **Command-Line Interface** - Full-featured CLI tool for all calculations

## Installation

### From Source

```bash
git clone https://github.com/yourusername/ballistics-engine.git
cd ballistics-engine
cargo build --release
```

### Building the CLI Tool

```bash
# Build the command-line interface
cargo build --release --bin ballistics-cli

# The binary will be at: target/release/ballistics-cli
```

## Quick Start

### Command-Line Interface

```bash
# Basic trajectory calculation
./ballistics-cli trajectory -v 850 -a 30 -b 0.5 -m 0.02

# Monte Carlo simulation
./ballistics-cli monte-carlo -v 850 -a 35 -n 1000 --velocity-std 5

# Show help
./ballistics-cli --help
```

### As a Rust Library

```rust
use ballistics_engine::{BallisticInputs, TrajectorySolver};

fn main() {
    let inputs = BallisticInputs {
        muzzle_velocity: 850.0,  // m/s
        launch_angle: 0.523599,  // 30 degrees in radians
        ballistic_coefficient: 0.5,
        mass: 0.02,              // kg
        diameter: 0.00762,       // meters
        ..Default::default()
    };
    
    let solver = TrajectorySolver::new(
        inputs, 
        Default::default(),  // Wind conditions
        Default::default()   // Atmosphere
    );
    
    let result = solver.solve().unwrap();
    println!("Max range: {:.2} m", result.max_range);
}
```

## CLI Commands

### Trajectory Calculation

Calculate single trajectory with environmental conditions:

```bash
./ballistics-cli trajectory \
  -v 850          # Velocity (m/s)
  -a 30           # Angle (degrees)
  -b 0.5          # Ballistic coefficient
  -m 0.02         # Mass (kg)
  -d 0.00762      # Diameter (meters)
  --wind-speed 10 # Wind speed (m/s)
  --wind-direction 90  # Wind direction (degrees)
  --temperature 20     # Temperature (Celsius)
  --pressure 1013     # Pressure (hPa)
  -o json         # Output format
```

### Monte Carlo Simulation

Run statistical analysis with parameter variations:

```bash
./ballistics-cli monte-carlo \
  -v 850          # Base velocity (m/s)
  -a 35           # Base angle (degrees)
  -n 1000         # Number of simulations
  --velocity-std 5    # Velocity std dev
  --angle-std 0.5     # Angle std dev
  -o summary      # Output format
```

### Output Formats

- **table** (default) - Formatted ASCII table
- **json** - Complete data in JSON format
- **csv** - Comma-separated values for analysis

## Examples

The `examples/` directory contains practical demonstrations:

| Example | Description |
|---------|-------------|
| `01_basic_trajectory.rs` | Simple trajectory calculation |
| `02_monte_carlo.rs` | Statistical simulation with uncertainties |
| `03_wind_effects.rs` | Wind drift and compensation |
| `04_bc_estimation.rs` | Estimate BC from observed data |
| `05_drag_comparison.rs` | Compare different projectiles |
| `06_data_export.rs` | Export data in multiple formats |

Run examples with:
```bash
cargo run --example 01_basic_trajectory
```

## Physics Model

### Trajectory Integration

The engine uses numerical integration with configurable time steps to solve the equations of motion:

- **Drag Force**: `F_d = 0.5 * ρ * C_d * A * v²`
- **Gravity**: Standard gravity (9.80665 m/s²)
- **Wind Effects**: 3D wind vector calculations
- **Atmospheric Model**: ISA atmosphere with corrections

### Drag Models

Supported drag models:
- **G1** - Standard projectile
- **G7** - Boat-tail bullets
- **Custom** - User-defined drag curves

### Statistical Analysis

Monte Carlo simulations provide:
- Mean and standard deviation
- CEP (Circular Error Probable)
- Hit probability calculations
- Distribution analysis

## API Documentation

### Core Types

```rust
pub struct BallisticInputs {
    pub muzzle_velocity: f64,      // m/s
    pub launch_angle: f64,          // radians
    pub ballistic_coefficient: f64,
    pub mass: f64,                  // kg
    pub diameter: f64,              // meters
    pub drag_model: DragModel,
    pub sight_height: f64,          // meters
}

pub struct TrajectoryResult {
    pub max_range: f64,
    pub max_height: f64,
    pub time_of_flight: f64,
    pub impact_velocity: f64,
    pub impact_energy: f64,
    pub points: Vec<TrajectoryPoint>,
}
```

### Key Functions

- `TrajectorySolver::solve()` - Calculate trajectory
- `run_monte_carlo()` - Run Monte Carlo simulation
- `calculate_zero_angle()` - Calculate sight zero
- `estimate_bc_from_trajectory()` - Estimate BC from data

## Performance

- Optimized Rust implementation
- Parallel Monte Carlo simulations
- Configurable time steps for accuracy/speed trade-off
- Efficient drag table interpolation

### Benchmarks

| Operation | Time |
|-----------|------|
| Single trajectory (1000m) | ~5ms |
| Monte Carlo (1000 runs) | ~500ms |
| BC estimation | ~50ms |

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

## Python Support

The engine includes Python bindings via PyO3:

```python
import ballistics_engine

# Use the Python API
result = ballistics_engine.calculate_trajectory(
    velocity=850,
    angle=30,
    bc=0.5,
    mass=0.02
)
```

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Development Setup

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Run tests: `cargo test`
5. Submit a pull request

## License

This project is dual-licensed under:
- MIT License
- Apache License 2.0

See LICENSE files for details.

## Acknowledgments

- Ballistics physics based on Robert McCoy's "Modern Exterior Ballistics"
- Drag tables from military ballistics research
- Community contributions and testing

## Support

For issues, questions, or contributions:
- GitHub Issues: [github.com/yourusername/ballistics-engine/issues](https://github.com/yourusername/ballistics-engine/issues)
- Documentation: [Full API Docs](https://docs.rs/ballistics-engine)

## Roadmap

- [ ] GPU acceleration for Monte Carlo
- [ ] Additional drag models (G2, G5, G6, G8)
- [ ] Coriolis effect calculations
- [ ] Spin drift modeling
- [ ] Web API interface
- [ ] Mobile app support

## Citation

If you use this engine in research, please cite:

```bibtex
@software{ballistics_engine,
  title = {Ballistics Engine: High-Performance Trajectory Calculator},
  author = {Ballistics Engine Contributors},
  year = {2025},
  url = {https://github.com/yourusername/ballistics-engine}
}
```