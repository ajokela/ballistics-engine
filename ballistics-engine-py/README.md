# ballistics-engine

High-performance ballistics trajectory engine with professional physics modeling.

## Features

- Professional-grade trajectory calculations with multiple drag models (G1, G7, G8)
- Advanced physics including wind effects and atmospheric modeling
- Fast Rust implementation with Python bindings via PyO3
- Imperial units API (grains, fps, yards, inches) with automatic metric conversion

## Installation

```bash
pip install ballistics-engine
```

## Quick Start

```python
from ballistics_engine import BallisticInputs, TrajectorySolver, WindConditions, AtmosphericConditions, DragModel

# Create ballistic inputs (all imperial units)
inputs = BallisticInputs(
    bc=0.505,                        # G7 BC
    bullet_weight_grains=168,        # grains
    muzzle_velocity_fps=2650,        # feet per second
    bullet_diameter_inches=0.308,    # inches
    bullet_length_inches=1.24,       # inches
    sight_height_inches=1.5,         # inches above bore
    zero_distance_yards=100,         # yards
    twist_rate_inches=11.25,         # inches per turn
)
inputs.drag_model = DragModel.g7()  # Use G7 drag model

# Create wind conditions (optional)
wind = WindConditions(
    speed_mph=10,                    # mph
    direction_degrees=90,            # degrees (0=headwind, 90=from right)
)

# Create atmospheric conditions (optional)
atmosphere = AtmosphericConditions(
    temperature_f=59,                # Fahrenheit
    pressure_inhg=29.92,             # inHg
    humidity_percent=50,             # percent
    altitude_feet=0,                 # feet
)

# Solve trajectory
solver = TrajectorySolver(inputs, wind=wind, atmosphere=atmosphere)
result = solver.solve()

# Print results
print(f"Max range: {result.max_range_yards:.1f} yards")
print(f"Time of flight: {result.time_of_flight:.2f} seconds")
print(f"Impact velocity: {result.impact_velocity_fps:.1f} fps")
print(f"Impact energy: {result.impact_energy_ftlbs:.1f} ft-lbs")

# Iterate through trajectory points
for point in result.points:
    print(f"Time: {point.time:.2f}s, X: {point.x:.1f}yd, Y: {point.y:.3f}yd, V: {point.velocity_fps:.1f}fps")
```

## Units

The Python API uses **imperial units** for convenience:

- **Mass**: grains (gr)
- **Velocity**: feet per second (fps)
- **Distance**: yards (yd) and inches (in)
- **Pressure**: inches of mercury (inHg)
- **Temperature**: Fahrenheit (°F)
- **Wind speed**: miles per hour (mph)

All conversions to metric (used internally by the Rust engine) are handled automatically.

## API Reference

### `BallisticInputs`

Main input parameters for trajectory calculation.

**Parameters:**
- `bc` (float): Ballistic coefficient
- `bullet_weight_grains` (float): Bullet mass in grains
- `muzzle_velocity_fps` (float): Muzzle velocity in fps
- `bullet_diameter_inches` (float): Bullet diameter in inches
- `bullet_length_inches` (float): Bullet length in inches
- `sight_height_inches` (float): Sight height above bore in inches
- `zero_distance_yards` (float): Zero distance in yards
- `shooting_angle_degrees` (float): Uphill/downhill angle in degrees
- `twist_rate_inches` (float): Barrel twist rate (inches per turn)
- `is_right_twist` (bool): True for right-hand twist (default: True)

### `WindConditions`

Wind parameters.

**Parameters:**
- `speed_mph` (float): Wind speed in mph (default: 0)
- `direction_degrees` (float): Wind direction in degrees (0=headwind, 90=from right, default: 0)

### `AtmosphericConditions`

Atmospheric parameters.

**Parameters:**
- `temperature_f` (float): Temperature in Fahrenheit (default: 59)
- `pressure_inhg` (float): Barometric pressure in inHg (default: 29.92)
- `humidity_percent` (float): Relative humidity percentage (default: 50)
- `altitude_feet` (float): Altitude in feet (default: 0)

### `TrajectorySolver`

Trajectory calculation engine.

**Methods:**
- `__init__(inputs, wind=None, atmosphere=None)`: Create solver with inputs
- `solve()`: Calculate trajectory, returns `TrajectoryResult`

### `TrajectoryResult`

Trajectory calculation results.

**Properties:**
- `max_range_yards` (float): Maximum range in yards
- `max_height_yards` (float): Maximum height in yards
- `time_of_flight` (float): Total flight time in seconds
- `impact_velocity_fps` (float): Impact velocity in fps
- `impact_energy_ftlbs` (float): Impact energy in ft-lbs
- `points` (list[TrajectoryPoint]): List of trajectory points

### `TrajectoryPoint`

Individual point along trajectory.

**Properties:**
- `time` (float): Time in seconds
- `x` (float): Downrange distance in yards
- `y` (float): Vertical position in yards (relative to line of sight)
- `z` (float): Lateral position in yards
- `velocity_fps` (float): Velocity in fps
- `energy_ftlbs` (float): Kinetic energy in ft-lbs

### `DragModel`

Drag model selection.

**Static methods:**
- `DragModel.g1()`: G1 drag model (flat base)
- `DragModel.g7()`: G7 drag model (boat tail)
- `DragModel.g8()`: G8 drag model (boat tail with meplat)

## License

Dual licensed under MIT or Apache-2.0.

## Links

- **Homepage**: https://ballistics.rs/
- **Repository**: https://github.com/ajokela/ballistics-engine
- **Documentation**: https://docs.rs/ballistics-engine
