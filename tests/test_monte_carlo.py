#!/usr/bin/env python3
import subprocess
import json
import sys

# Run Monte Carlo simulation
cmd = [
    "./target/debug/ballistics", "monte-carlo",
    "--velocity", "800",
    "--angle", "45",
    "--bc", "0.5",
    "--mass", "0.01",
    "--diameter", "0.008",
    "--num-sims", "50",
    "--velocity-std", "2",
    "--angle-std", "0.05",
    "--bc-std", "0.005",
    "--output", "json"
]

result = subprocess.run(cmd, capture_output=True, text=True)

if result.returncode != 0:
    print(f"Error running command: {result.stderr}")
    sys.exit(1)

try:
    data = json.loads(result.stdout)
    
    # Check if we have impact positions
    if "impact_positions" in data:
        positions = data["impact_positions"]
        print(f"Number of impact positions: {len(positions)}")
        
        # Extract x and z coordinates (horizontal spread)
        x_coords = [pos[0] for pos in positions]
        z_coords = [pos[2] for pos in positions]
        
        # Calculate spread statistics
        x_min, x_max = min(x_coords), max(x_coords)
        z_min, z_max = min(z_coords), max(z_coords)
        
        print(f"X spread: {x_min:.2f} to {x_max:.2f} m (range: {x_max - x_min:.2f} m)")
        print(f"Z spread: {z_min:.2f} to {z_max:.2f} m (range: {z_max - z_min:.2f} m)")
        
        # Simple ASCII visualization
        print("\nImpact pattern (top-down view, scaled):")
        
        # Scale to 40x20 grid
        grid_width = 40
        grid_height = 20
        
        if x_max > x_min and z_max > z_min:
            grid = [[' ' for _ in range(grid_width)] for _ in range(grid_height)]
            
            for x, z in zip(x_coords, z_coords):
                # Normalize to grid coordinates
                grid_x = int((x - x_min) / (x_max - x_min) * (grid_width - 1))
                grid_z = int((z - z_min) / (z_max - z_min) * (grid_height - 1))
                grid[grid_z][grid_x] = '*'
            
            # Print grid
            for row in grid:
                print(''.join(row))
        else:
            print("All impacts at same position - checking for vertical spread only")
            
    else:
        print("No impact positions in output")
        print("Available keys:", list(data.keys()))
        
except json.JSONDecodeError as e:
    print(f"Failed to parse JSON: {e}")
    print("Output:", result.stdout[:500])