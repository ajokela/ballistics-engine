/// Data Export Example
/// 
/// This example demonstrates how to export trajectory data in various formats
/// for analysis, visualization, and integration with other tools.

use std::fs::File;
use std::io::Write;
use std::f64::consts::PI;

fn main() -> std::io::Result<()> {
    println!("=== Trajectory Data Export Example ===\n");
    
    // Projectile parameters
    let velocity = 850.0;      // m/s
    let angle_deg = 35.0;      // degrees  
    let mass = 0.01166;        // kg (180 grains)
    let diameter = 0.00762;    // meters (.308)
    let bc = 0.45;
    
    println!("Generating trajectory data for:");
    println!("  Velocity: {} m/s", velocity);
    println!("  Angle: {}°", angle_deg);
    println!("  Mass: {} kg", mass);
    println!("  BC: {}", bc);
    println!();
    
    // Calculate full trajectory
    let trajectory = calculate_full_trajectory(velocity, angle_deg, mass, diameter, bc);
    
    // Export to different formats
    export_csv(&trajectory, "trajectory.csv")?;
    println!("✓ Exported to trajectory.csv");
    
    export_json(&trajectory, "trajectory.json")?;
    println!("✓ Exported to trajectory.json");
    
    export_gnuplot(&trajectory, "trajectory.dat")?;
    println!("✓ Exported to trajectory.dat (Gnuplot format)");
    
    export_matlab(&trajectory, "trajectory.m")?;
    println!("✓ Exported to trajectory.m (MATLAB/Octave format)");
    
    export_python(&trajectory, "trajectory.py")?;
    println!("✓ Exported to trajectory.py (Python format)");
    
    // Generate summary statistics
    println!("\nTrajectory Summary:");
    println!("  Total points: {}", trajectory.len());
    println!("  Time span: {:.2} seconds", trajectory.last().unwrap().time);
    println!("  Max range: {:.2} m", trajectory.last().unwrap().x);
    
    let max_height = trajectory.iter()
        .map(|p| p.y)
        .fold(f64::NEG_INFINITY, f64::max);
    println!("  Max height: {:.2} m", max_height);
    
    // Create a simplified ballistic table
    export_ballistic_table(&trajectory, "ballistic_table.txt")?;
    println!("\n✓ Exported ballistic table to ballistic_table.txt");
    
    // Generate plotting script
    create_gnuplot_script("plot_trajectory.gnuplot")?;
    println!("✓ Created Gnuplot script: plot_trajectory.gnuplot");
    println!("  Run: gnuplot plot_trajectory.gnuplot");
    
    create_python_plot_script("plot_trajectory_py.py")?;
    println!("✓ Created Python plot script: plot_trajectory_py.py");
    println!("  Run: python plot_trajectory_py.py");
    
    Ok(())
}

#[derive(Debug, Clone)]
struct TrajectoryPoint {
    time: f64,
    x: f64,
    y: f64,
    z: f64,
    vx: f64,
    vy: f64,
    vz: f64,
    velocity: f64,
    energy: f64,
    mach: f64,
}

fn calculate_full_trajectory(
    velocity: f64,
    angle_deg: f64,
    mass: f64,
    diameter: f64,
    bc: f64,
) -> Vec<TrajectoryPoint> {
    let angle_rad = angle_deg * PI / 180.0;
    let dt = 0.01;
    let mut trajectory = Vec::new();
    
    let mut t = 0.0;
    let mut x = 0.0;
    let mut y = 0.0;
    let mut z = 0.0;
    let mut vx = velocity * angle_rad.cos();
    let mut vy = velocity * angle_rad.sin();
    let mut vz = 0.0;
    
    let g = 9.80665;
    let air_density = 1.225;
    let area = PI * (diameter / 2.0).powi(2);
    let cd = 0.47 / bc;
    let speed_of_sound = 343.0;
    
    while y >= 0.0 && t < 100.0 {
        let v = (vx * vx + vy * vy + vz * vz).sqrt();
        let energy = 0.5 * mass * v * v;
        let mach = v / speed_of_sound;
        
        trajectory.push(TrajectoryPoint {
            time: t,
            x, y, z,
            vx, vy, vz,
            velocity: v,
            energy,
            mach,
        });
        
        let drag = 0.5 * air_density * cd * area * v;
        let ax = -drag * vx / mass;
        let ay = -drag * vy / mass - g;
        let az = -drag * vz / mass;
        
        vx += ax * dt;
        vy += ay * dt;
        vz += az * dt;
        
        x += vx * dt;
        y += vy * dt;
        z += vz * dt;
        t += dt;
    }
    
    trajectory
}

fn export_csv(trajectory: &[TrajectoryPoint], filename: &str) -> std::io::Result<()> {
    let mut file = File::create(filename)?;
    
    // Write header
    writeln!(file, "time,x,y,z,vx,vy,vz,velocity,energy,mach")?;
    
    // Write data
    for point in trajectory {
        writeln!(file, "{:.4},{:.3},{:.3},{:.3},{:.2},{:.2},{:.2},{:.2},{:.2},{:.3}",
            point.time, point.x, point.y, point.z,
            point.vx, point.vy, point.vz,
            point.velocity, point.energy, point.mach)?;
    }
    
    Ok(())
}

fn export_json(trajectory: &[TrajectoryPoint], filename: &str) -> std::io::Result<()> {
    let mut file = File::create(filename)?;
    
    writeln!(file, "{{")?;
    writeln!(file, "  \"trajectory\": [")?;
    
    for (i, point) in trajectory.iter().enumerate() {
        write!(file, "    {{")?;
        write!(file, "\"time\": {:.4}, ", point.time)?;
        write!(file, "\"position\": [{:.3}, {:.3}, {:.3}], ", point.x, point.y, point.z)?;
        write!(file, "\"velocity\": [{:.2}, {:.2}, {:.2}], ", point.vx, point.vy, point.vz)?;
        write!(file, "\"speed\": {:.2}, ", point.velocity)?;
        write!(file, "\"energy\": {:.2}, ", point.energy)?;
        write!(file, "\"mach\": {:.3}", point.mach)?;
        write!(file, "}}")?;
        
        if i < trajectory.len() - 1 {
            writeln!(file, ",")?;
        } else {
            writeln!(file)?;
        }
    }
    
    writeln!(file, "  ]")?;
    writeln!(file, "}}")?;
    
    Ok(())
}

fn export_gnuplot(trajectory: &[TrajectoryPoint], filename: &str) -> std::io::Result<()> {
    let mut file = File::create(filename)?;
    
    writeln!(file, "# Trajectory data for Gnuplot")?;
    writeln!(file, "# time x y z velocity energy mach")?;
    
    for point in trajectory {
        writeln!(file, "{:.4} {:.3} {:.3} {:.3} {:.2} {:.2} {:.3}",
            point.time, point.x, point.y, point.z,
            point.velocity, point.energy, point.mach)?;
    }
    
    Ok(())
}

fn export_matlab(trajectory: &[TrajectoryPoint], filename: &str) -> std::io::Result<()> {
    let mut file = File::create(filename)?;
    
    writeln!(file, "% Trajectory data for MATLAB/Octave")?;
    writeln!(file, "% Generated by ballistics engine")?;
    writeln!(file)?;
    
    // Export as MATLAB arrays
    write!(file, "t = [")?;
    for (i, p) in trajectory.iter().enumerate() {
        if i > 0 { write!(file, ", ")?; }
        write!(file, "{:.4}", p.time)?;
    }
    writeln!(file, "];")?;
    
    write!(file, "x = [")?;
    for (i, p) in trajectory.iter().enumerate() {
        if i > 0 { write!(file, ", ")?; }
        write!(file, "{:.3}", p.x)?;
    }
    writeln!(file, "];")?;
    
    write!(file, "y = [")?;
    for (i, p) in trajectory.iter().enumerate() {
        if i > 0 { write!(file, ", ")?; }
        write!(file, "{:.3}", p.y)?;
    }
    writeln!(file, "];")?;
    
    write!(file, "velocity = [")?;
    for (i, p) in trajectory.iter().enumerate() {
        if i > 0 { write!(file, ", ")?; }
        write!(file, "{:.2}", p.velocity)?;
    }
    writeln!(file, "];")?;
    
    writeln!(file)?;
    writeln!(file, "% Plot the trajectory")?;
    writeln!(file, "figure;")?;
    writeln!(file, "subplot(2,1,1);")?;
    writeln!(file, "plot(x, y);")?;
    writeln!(file, "xlabel('Distance (m)');")?;
    writeln!(file, "ylabel('Height (m)');")?;
    writeln!(file, "title('Ballistic Trajectory');")?;
    writeln!(file, "grid on;")?;
    writeln!(file)?;
    writeln!(file, "subplot(2,1,2);")?;
    writeln!(file, "plot(x, velocity);")?;
    writeln!(file, "xlabel('Distance (m)');")?;
    writeln!(file, "ylabel('Velocity (m/s)');")?;
    writeln!(file, "title('Velocity vs Distance');")?;
    writeln!(file, "grid on;")?;
    
    Ok(())
}

fn export_python(trajectory: &[TrajectoryPoint], filename: &str) -> std::io::Result<()> {
    let mut file = File::create(filename)?;
    
    writeln!(file, "#!/usr/bin/env python3")?;
    writeln!(file, "# Trajectory data for Python")?;
    writeln!(file, "# Generated by ballistics engine")?;
    writeln!(file)?;
    writeln!(file, "import numpy as np")?;
    writeln!(file)?;
    
    // Export as numpy arrays
    write!(file, "t = np.array([")?;
    for (i, p) in trajectory.iter().enumerate() {
        if i > 0 { write!(file, ", ")?; }
        write!(file, "{:.4}", p.time)?;
    }
    writeln!(file, "])")?;
    
    write!(file, "x = np.array([")?;
    for (i, p) in trajectory.iter().enumerate() {
        if i > 0 { write!(file, ", ")?; }
        write!(file, "{:.3}", p.x)?;
    }
    writeln!(file, "])")?;
    
    write!(file, "y = np.array([")?;
    for (i, p) in trajectory.iter().enumerate() {
        if i > 0 { write!(file, ", ")?; }
        write!(file, "{:.3}", p.y)?;
    }
    writeln!(file, "])")?;
    
    write!(file, "velocity = np.array([")?;
    for (i, p) in trajectory.iter().enumerate() {
        if i > 0 { write!(file, ", ")?; }
        write!(file, "{:.2}", p.velocity)?;
    }
    writeln!(file, "])")?;
    
    writeln!(file)?;
    writeln!(file, "# Access the data:")?;
    writeln!(file, "# print(f'Max range: {{x[-1]:.2f}} m')")?;
    writeln!(file, "# print(f'Max height: {{y.max():.2f}} m')")?;
    writeln!(file, "# print(f'Time of flight: {{t[-1]:.2f}} s')")?;
    
    Ok(())
}

fn export_ballistic_table(trajectory: &[TrajectoryPoint], filename: &str) -> std::io::Result<()> {
    let mut file = File::create(filename)?;
    
    writeln!(file, "BALLISTIC TABLE")?;
    writeln!(file, "===============")?;
    writeln!(file)?;
    writeln!(file, " Range | Drop  | Drift | Time  | Velocity | Energy | Mach")?;
    writeln!(file, "  (m)  |  (m)  |  (m)  | (s)   |  (m/s)   |  (J)   |  #")?;
    writeln!(file, "-------|-------|-------|-------|----------|--------|------")?;
    
    // Export at specific ranges
    let ranges = vec![100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0];
    
    for target_range in ranges {
        // Find the closest point to target range
        if let Some(point) = trajectory.iter().find(|p| p.x >= target_range) {
            writeln!(file, " {:5.0} | {:5.2} | {:5.2} | {:5.3} | {:8.1} | {:6.0} | {:.2}",
                point.x, -point.y, point.z, point.time,
                point.velocity, point.energy, point.mach)?;
        }
    }
    
    Ok(())
}

fn create_gnuplot_script(filename: &str) -> std::io::Result<()> {
    let mut file = File::create(filename)?;
    
    writeln!(file, "#!/usr/bin/gnuplot")?;
    writeln!(file, "# Gnuplot script for trajectory visualization")?;
    writeln!(file)?;
    writeln!(file, "set terminal png size 1200,800")?;
    writeln!(file, "set output 'trajectory_plot.png'")?;
    writeln!(file)?;
    writeln!(file, "set multiplot layout 2,2")?;
    writeln!(file)?;
    writeln!(file, "# Plot 1: Trajectory")?;
    writeln!(file, "set xlabel 'Distance (m)'")?;
    writeln!(file, "set ylabel 'Height (m)'")?;
    writeln!(file, "set title 'Ballistic Trajectory'")?;
    writeln!(file, "set grid")?;
    writeln!(file, "plot 'trajectory.dat' using 2:3 with lines title 'Trajectory' lw 2")?;
    writeln!(file)?;
    writeln!(file, "# Plot 2: Velocity")?;
    writeln!(file, "set xlabel 'Distance (m)'")?;
    writeln!(file, "set ylabel 'Velocity (m/s)'")?;
    writeln!(file, "set title 'Velocity vs Distance'")?;
    writeln!(file, "plot 'trajectory.dat' using 2:5 with lines title 'Velocity' lw 2")?;
    writeln!(file)?;
    writeln!(file, "# Plot 3: Energy")?;
    writeln!(file, "set xlabel 'Time (s)'")?;
    writeln!(file, "set ylabel 'Energy (J)'")?;
    writeln!(file, "set title 'Kinetic Energy vs Time'")?;
    writeln!(file, "plot 'trajectory.dat' using 1:6 with lines title 'Energy' lw 2")?;
    writeln!(file)?;
    writeln!(file, "# Plot 4: Mach Number")?;
    writeln!(file, "set xlabel 'Distance (m)'")?;
    writeln!(file, "set ylabel 'Mach Number'")?;
    writeln!(file, "set title 'Mach Number vs Distance'")?;
    writeln!(file, "plot 'trajectory.dat' using 2:7 with lines title 'Mach' lw 2")?;
    writeln!(file)?;
    writeln!(file, "unset multiplot")?;
    
    Ok(())
}

fn create_python_plot_script(filename: &str) -> std::io::Result<()> {
    let mut file = File::create(filename)?;
    
    writeln!(file, "#!/usr/bin/env python3")?;
    writeln!(file, "import matplotlib.pyplot as plt")?;
    writeln!(file, "import pandas as pd")?;
    writeln!(file)?;
    writeln!(file, "# Load trajectory data")?;
    writeln!(file, "data = pd.read_csv('trajectory.csv')")?;
    writeln!(file)?;
    writeln!(file, "# Create figure with subplots")?;
    writeln!(file, "fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 8))")?;
    writeln!(file)?;
    writeln!(file, "# Plot trajectory")?;
    writeln!(file, "ax1.plot(data['x'], data['y'], 'b-', linewidth=2)")?;
    writeln!(file, "ax1.set_xlabel('Distance (m)')")?;
    writeln!(file, "ax1.set_ylabel('Height (m)')")?;
    writeln!(file, "ax1.set_title('Ballistic Trajectory')")?;
    writeln!(file, "ax1.grid(True)")?;
    writeln!(file)?;
    writeln!(file, "# Plot velocity")?;
    writeln!(file, "ax2.plot(data['x'], data['velocity'], 'r-', linewidth=2)")?;
    writeln!(file, "ax2.set_xlabel('Distance (m)')")?;
    writeln!(file, "ax2.set_ylabel('Velocity (m/s)')")?;
    writeln!(file, "ax2.set_title('Velocity vs Distance')")?;
    writeln!(file, "ax2.grid(True)")?;
    writeln!(file)?;
    writeln!(file, "# Plot energy")?;
    writeln!(file, "ax3.plot(data['time'], data['energy'], 'g-', linewidth=2)")?;
    writeln!(file, "ax3.set_xlabel('Time (s)')")?;
    writeln!(file, "ax3.set_ylabel('Energy (J)')")?;
    writeln!(file, "ax3.set_title('Kinetic Energy vs Time')")?;
    writeln!(file, "ax3.grid(True)")?;
    writeln!(file)?;
    writeln!(file, "# Plot Mach number")?;
    writeln!(file, "ax4.plot(data['x'], data['mach'], 'm-', linewidth=2)")?;
    writeln!(file, "ax4.axhline(y=1.0, color='k', linestyle='--', alpha=0.5, label='Mach 1')")?;
    writeln!(file, "ax4.set_xlabel('Distance (m)')")?;
    writeln!(file, "ax4.set_ylabel('Mach Number')")?;
    writeln!(file, "ax4.set_title('Mach Number vs Distance')")?;
    writeln!(file, "ax4.grid(True)")?;
    writeln!(file, "ax4.legend()")?;
    writeln!(file)?;
    writeln!(file, "plt.tight_layout()")?;
    writeln!(file, "plt.savefig('trajectory_plot.png', dpi=150)")?;
    writeln!(file, "plt.show()")?;
    
    Ok(())
}