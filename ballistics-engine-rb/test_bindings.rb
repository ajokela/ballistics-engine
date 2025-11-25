#!/usr/bin/env ruby
# frozen_string_literal: true

require_relative 'lib/ballistics_engine'

puts "Testing Ballistics Engine Ruby Bindings"
puts "=" * 50

# Test 1: Basic trajectory calculation
puts "\n1. Basic trajectory calculation (168gr .308 Winchester)"
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

puts inputs.to_s

# Create solver without wind/atmosphere (uses defaults)
solver = BallisticsEngine::TrajectorySolver.new(inputs, nil, nil)
result = solver.solve

puts "\nResults (no wind):"
puts "  Max range: #{result.max_range_yards.round(1)} yards"
puts "  Max height: #{result.max_height_yards.round(2)} yards"
puts "  Time of flight: #{result.time_of_flight.round(3)} seconds"
puts "  Impact velocity: #{result.impact_velocity_fps.round(1)} fps"
puts "  Impact energy: #{result.impact_energy_ftlbs.round(1)} ft-lbs"
puts "  Number of points: #{result.points.length}"

# Test 2: Trajectory with wind
puts "\n2. Trajectory with wind (10 mph from right)"
wind = BallisticsEngine::WindConditions.new(10.0, 90.0)
solver_wind = BallisticsEngine::TrajectorySolver.new(inputs, wind, nil)
result_wind = solver_wind.solve

puts "\nResults (with wind):"
puts "  Max range: #{result_wind.max_range_yards.round(1)} yards"
puts "  Time of flight: #{result_wind.time_of_flight.round(3)} seconds"
puts "  Impact velocity: #{result_wind.impact_velocity_fps.round(1)} fps"

# Test 3: Custom atmospheric conditions
puts "\n3. Custom atmospheric conditions (high altitude, cold)"
atmosphere = BallisticsEngine::AtmosphericConditions.new(
  20.0,       # temperature (F)
  25.84,      # pressure (inHg) - ~5000ft elevation
  30.0,       # humidity (%)
  5000.0      # altitude (feet)
)
solver_altitude = BallisticsEngine::TrajectorySolver.new(inputs, nil, atmosphere)
result_altitude = solver_altitude.solve

puts "\nResults (high altitude):"
puts "  Max range: #{result_altitude.max_range_yards.round(1)} yards"
puts "  Impact velocity: #{result_altitude.impact_velocity_fps.round(1)} fps"

# Test 4: Trajectory point details
puts "\n4. Sample trajectory points (every 100 yards):"
puts "  Range(yd)  Drop(yd)   Drift(yd)  Velocity(fps)  Energy(ft-lbs)"
result_wind.points.each do |point|
  if point.x % 100.0 < 1.0  # Close to 100 yard intervals
    puts "  #{point.x.round(0).to_s.rjust(8)}  " \
         "#{point.y.round(2).to_s.rjust(8)}  " \
         "#{point.z.round(2).to_s.rjust(9)}  " \
         "#{point.velocity_fps.round(1).to_s.rjust(13)}  " \
         "#{point.energy_ftlbs.round(1).to_s.rjust(13)}"
  end
end

# Test 5: Drag models
puts "\n5. Testing drag models"
puts "  G1: #{BallisticsEngine::DragModel.g1.to_s}"
puts "  G7: #{BallisticsEngine::DragModel.g7.to_s}"
puts "  G8: #{BallisticsEngine::DragModel.g8.to_s}"

puts "\n✓ All tests passed!"
