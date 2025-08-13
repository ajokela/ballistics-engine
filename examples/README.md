# Ballistics Engine Examples

This directory contains practical examples demonstrating various features and use cases of the ballistics engine.

## Examples Overview

### 01_basic_trajectory.rs
**Basic Trajectory Calculation**
- Simple ballistic trajectory computation
- Key trajectory parameters (range, height, time of flight)
- Step-by-step trajectory points
- Energy calculations

```bash
cargo run --example 01_basic_trajectory
```

### 02_monte_carlo.rs
**Monte Carlo Simulation**
- Statistical analysis with parameter uncertainties
- Normal distribution sampling
- Range and impact statistics
- Dispersion analysis (CEP, R95)
- Distribution histogram visualization

```bash
cargo run --example 02_monte_carlo
```

### 03_wind_effects.rs
**Wind Effects on Trajectories**
- Headwind and tailwind effects
- Crosswind drift calculations
- Wind correction tables
- MOA (Minutes of Angle) conversions
- Distance-based drift analysis

```bash
cargo run --example 03_wind_effects
```

### 04_bc_estimation.rs
**Ballistic Coefficient Estimation**
- BC estimation from observed trajectory data
- Iterative optimization methods
- Least squares fitting
- Error analysis and verification
- BC sensitivity analysis

```bash
cargo run --example 04_bc_estimation
```

### 05_drag_comparison.rs
**Projectile Comparison**
- Multiple caliber comparisons (.22 LR to .50 BMG)
- Drop tables at various distances
- Time of flight calculations
- Velocity retention analysis
- Energy comparisons
- Performance metrics (supersonic range, max range)

```bash
cargo run --example 05_drag_comparison
```

### 06_data_export.rs
**Data Export and Visualization**
- Export to multiple formats:
  - CSV for spreadsheets
  - JSON for web applications
  - Gnuplot data files
  - MATLAB/Octave scripts
  - Python numpy arrays
- Ballistic tables generation
- Automatic plot script creation
- Visualization templates

```bash
cargo run --example 06_data_export
```

## Running Examples

To run any example:

```bash
# Run a specific example
cargo run --example <example_name>

# Run with release optimizations
cargo run --release --example <example_name>

# Run all examples
for example in 01_basic_trajectory 02_monte_carlo 03_wind_effects 04_bc_estimation 05_drag_comparison 06_data_export; do
    echo "Running $example..."
    cargo run --example $example
    echo ""
done
```

## Key Concepts Demonstrated

### Physics Modeling
- Drag force calculations
- Atmospheric effects
- Gravity integration
- Wind drift modeling
- Energy dissipation

### Numerical Methods
- Euler integration
- Root finding (zeroing)
- Least squares optimization
- Monte Carlo sampling
- Statistical analysis

### Data Processing
- Trajectory interpolation
- Statistical calculations
- Error analysis
- Data formatting
- Visualization preparation

## Common Use Cases

### Long Range Shooting
- Use `03_wind_effects.rs` for wind compensation
- Use `04_bc_estimation.rs` to determine BC from field data
- Use `02_monte_carlo.rs` for hit probability analysis

### Ballistic Research
- Use `05_drag_comparison.rs` for projectile selection
- Use `06_data_export.rs` for detailed analysis
- Use `01_basic_trajectory.rs` for quick calculations

### System Development
- Use examples as templates for integration
- Modify parameters for specific projectiles
- Extend with custom drag models
- Add specialized output formats

## Output Files

Some examples generate output files:

- `trajectory.csv` - Trajectory data in CSV format
- `trajectory.json` - Trajectory data in JSON format
- `trajectory.dat` - Gnuplot-compatible data
- `trajectory.m` - MATLAB/Octave script
- `trajectory.py` - Python data arrays
- `ballistic_table.txt` - Formatted ballistic table
- `plot_trajectory.gnuplot` - Gnuplot visualization script
- `plot_trajectory_py.py` - Python matplotlib script

## Customization

Each example can be easily modified:

1. **Change projectile parameters**: Modify mass, diameter, BC values
2. **Adjust environmental conditions**: Temperature, pressure, altitude
3. **Modify output formats**: Add custom fields or formatting
4. **Extend calculations**: Add new physics models or corrections

## Performance Notes

- Examples use simplified physics models for clarity
- Production use should employ the full engine capabilities
- Time steps can be adjusted for accuracy vs. performance
- Monte Carlo simulations benefit from parallelization

## Contributing

To add a new example:

1. Create a new file: `examples/XX_feature_name.rs`
2. Include clear documentation in the code
3. Add the example to this README
4. Ensure it can run standalone
5. Include sample output in comments