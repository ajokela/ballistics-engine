# Ballistics Engine - Ruby Bindings

High-performance ballistics calculations library for Ruby, powered by Rust.

## Features

- **4-DOF Trajectory Modeling** - Complete trajectory calculations with realistic physics
- **Multiple Drag Models** - G1, G7, and G8 ballistic coefficients
- **Wind Deflection** - Accurate wind drift calculations
- **Atmospheric Effects** - Temperature, pressure, humidity, and altitude compensation
- **Unit Conversion** - Automatic handling of imperial/metric conversions
- **High Performance** - Rust-based calculations for maximum speed

## Installation

Add this line to your application's Gemfile:

```ruby
gem 'ballistics-engine'
```

And then execute:

```bash
$ bundle install
```

Or install it yourself as:

```bash
$ gem install ballistics-engine
```

## Quick Start

```ruby
require 'ballistics_engine'

# Create ballistic inputs (168gr .308 Winchester)
inputs = BallisticsEngine::BallisticInputs.new(
  0.223,      # BC (G7)
  168.0,      # bullet weight (grains)
  2650.0,     # muzzle velocity (fps)
  0.308,      # bullet diameter (inches)
  1.2,        # bullet length (inches)
  1.5,        # sight height (inches)
  100.0,      # zero distance (yards)
  0.0,        # shooting angle (degrees)
  11.25,      # twist rate (inches)
  true        # right-hand twist
)

# Optional: Add wind conditions
wind = BallisticsEngine::WindConditions.new(
  10.0,       # speed (mph)
  90.0        # direction (degrees, 90 = from right)
)

# Optional: Add atmospheric conditions
atmosphere = BallisticsEngine::AtmosphericConditions.new(
  59.0,       # temperature (F)
  29.92,      # pressure (inHg)
  50.0,       # humidity (%)
  0.0         # altitude (feet)
)

# Create solver and calculate trajectory
solver = BallisticsEngine::TrajectorySolver.new(inputs, wind, atmosphere)
result = solver.solve

# Access results
puts "Max range: #{result.max_range_yards.round(1)} yards"
puts "Time of flight: #{result.time_of_flight.round(3)} seconds"
puts "Impact velocity: #{result.impact_velocity_fps.round(1)} fps"
puts "Impact energy: #{result.impact_energy_ftlbs.round(1)} ft-lbs"

# Iterate through trajectory points
result.points.each do |point|
  puts "Range: #{point.x.round(1)}yd, Drop: #{point.y.round(2)}yd, Velocity: #{point.velocity_fps.round(1)}fps"
end
```

## API Reference

### Classes

#### `BallisticInputs`

Ballistic calculation parameters.

**Constructor:**
```ruby
BallisticInputs.new(
  bc,                       # Ballistic coefficient
  bullet_weight_grains,     # Bullet weight in grains
  muzzle_velocity_fps,      # Muzzle velocity in fps
  bullet_diameter_inches,   # Bullet diameter in inches
  bullet_length_inches,     # Bullet length in inches
  sight_height_inches,      # Sight height in inches
  zero_distance_yards,      # Zero distance in yards
  shooting_angle_degrees,   # Shooting angle in degrees
  twist_rate_inches,        # Barrel twist rate in inches
  is_right_twist            # Right-hand twist? (boolean)
)
```

**Attributes:**
- All constructor parameters are accessible as read/write attributes

#### `WindConditions`

Wind parameters.

**Constructor:**
```ruby
WindConditions.new(
  speed_mph,               # Wind speed in mph
  direction_degrees        # Wind direction in degrees (0=headwind, 90=from right)
)
```

**Attributes:**
- `speed_mph` - Wind speed in mph (read/write)
- `direction_degrees` - Wind direction in degrees (read/write)

#### `AtmosphericConditions`

Atmospheric parameters.

**Constructor:**
```ruby
AtmosphericConditions.new(
  temperature_f,           # Temperature in Fahrenheit
  pressure_inhg,           # Pressure in inches of mercury
  humidity_percent,        # Relative humidity (0-100)
  altitude_feet            # Altitude in feet
)
```

**Attributes:**
- `temperature_f` - Temperature in Fahrenheit (read/write)
- `pressure_inhg` - Pressure in inHg (read/write)
- `humidity_percent` - Humidity percentage (read/write)
- `altitude_feet` - Altitude in feet (read/write)

#### `TrajectorySolver`

Trajectory calculation engine.

**Constructor:**
```ruby
TrajectorySolver.new(
  inputs,                  # BallisticInputs object
  wind = nil,              # WindConditions object (optional)
  atmosphere = nil         # AtmosphericConditions object (optional)
)
```

**Methods:**
- `solve()` - Calculate trajectory, returns `TrajectoryResult`

#### `TrajectoryResult`

Trajectory calculation results.

**Attributes:**
- `max_range_yards` - Maximum range in yards
- `max_height_yards` - Maximum height in yards
- `time_of_flight` - Time of flight in seconds
- `impact_velocity_fps` - Impact velocity in fps
- `impact_energy_ftlbs` - Impact energy in ft-lbs
- `points` - Array of `TrajectoryPoint` objects

#### `TrajectoryPoint`

Individual point along trajectory.

**Attributes:**
- `time` - Time in seconds
- `x` - Downrange distance in yards
- `y` - Vertical position in yards (relative to line of sight)
- `z` - Lateral drift in yards
- `velocity_fps` - Velocity in fps
- `energy_ftlbs` - Energy in ft-lbs

#### `DragModel`

Ballistic coefficient drag model.

**Class Methods:**
- `DragModel.g1()` - G1 drag model
- `DragModel.g7()` - G7 drag model
- `DragModel.g8()` - G8 drag model

## Development

After checking out the repo, run `bundle install` to install dependencies.

To build the native extension:

```bash
bundle exec rake compile
```

To run tests:

```bash
bundle exec rake test
```

## Contributing

Bug reports and pull requests are welcome on GitHub at https://github.com/ajokela/ballistics-engine.

## License

This project is licensed under either of:

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
- MIT license ([LICENSE](LICENSE) or http://opensource.org/licenses/MIT)

at your option.
