#include <stdio.h>
#include <math.h>

// FFI structure definitions (normally these would be in a header file)
typedef struct {
    double muzzle_velocity;
    double launch_angle;
    double ballistic_coefficient;
    double mass;
    double diameter;
    int drag_model;
    double sight_height;
    double target_distance;
    double temperature;
    double twist_rate;
    int is_twist_right;
    double shooting_angle;
    double altitude;
    double latitude;
    double azimuth_angle;
    int use_rk4;
    int enable_wind_shear;
    int enable_trajectory_sampling;
    double sample_interval;
    int enable_pitch_damping;
    int enable_precession_nutation;
    int enable_spin_drift;
    int enable_magnus;
    int enable_coriolis;
} FFIBallisticInputs;

typedef struct {
    double speed;
    double direction;
} FFIWindConditions;

typedef struct {
    double temperature;
    double pressure;
    double humidity;
    double altitude;
} FFIAtmosphericConditions;

typedef struct {
    double distance;
    double time;
    double velocity_mps;
    double energy_joules;
    double drop_meters;
    double windage_meters;
    double mach;
    double spin_rate_rps;
} FFITrajectorySample;

typedef struct {
    double max_range;
    double max_height;
    double time_of_flight;
    double impact_velocity;
    double impact_energy;
    void* points;
    int point_count;
    FFITrajectorySample* sampled_points;
    int sampled_point_count;
    double min_pitch_damping;
    double transonic_mach;
    double final_pitch_angle;
    double final_yaw_angle;
    double max_yaw_angle;
    double max_precession_angle;
} FFITrajectoryResult;

// External FFI functions
extern FFITrajectoryResult* ballistics_calculate_trajectory(
    const FFIBallisticInputs* inputs,
    const FFIWindConditions* wind,
    const FFIAtmosphericConditions* atmosphere,
    double max_range,
    double step_size
);

extern void ballistics_free_trajectory_result(FFITrajectoryResult* result);

int main() {
    printf("Testing Ballistics Engine FFI with Advanced Features\n");
    printf("====================================================\n\n");
    
    // Set up inputs for .308 Winchester, 168gr
    FFIBallisticInputs inputs = {
        .muzzle_velocity = 823.0,      // 2700 fps in m/s
        .launch_angle = 0.0872665,     // 5 degrees in radians
        .ballistic_coefficient = 0.475,
        .mass = 0.0109,                // 168 grains in kg
        .diameter = 0.00782,            // 0.308 inches in meters
        .drag_model = 0,                // G1
        .sight_height = 0.05,           // meters
        .target_distance = 200.0,       // meters
        .temperature = 15.0,            // Celsius
        .twist_rate = 10.0,             // 1:10 twist
        .is_twist_right = 1,            // Right-hand twist
        .shooting_angle = 0.0,
        .altitude = 0.0,
        .latitude = 45.0,               // 45 degrees north
        .azimuth_angle = 0.0,
        .use_rk4 = 1,                   // Use RK4 integration
        .enable_wind_shear = 1,         // Enable wind shear
        .enable_trajectory_sampling = 1, // Enable sampling
        .sample_interval = 50.0,        // Sample every 50 meters
        .enable_pitch_damping = 1,      // Enable pitch damping
        .enable_precession_nutation = 1,// Enable precession/nutation
        .enable_spin_drift = 1,         // Enable spin drift
        .enable_magnus = 1,             // Enable Magnus effect
        .enable_coriolis = 0            // Disable Coriolis for now
    };
    
    FFIWindConditions wind = {
        .speed = 5.0,                   // 5 m/s wind
        .direction = 1.5708             // 90 degrees (from right)
    };
    
    FFIAtmosphericConditions atmosphere = {
        .temperature = 15.0,
        .pressure = 1013.25,
        .humidity = 50.0,
        .altitude = 0.0
    };
    
    // Calculate trajectory
    FFITrajectoryResult* result = ballistics_calculate_trajectory(
        &inputs, &wind, &atmosphere, 1000.0, 1.0
    );
    
    if (result == NULL) {
        printf("Error: Failed to calculate trajectory\n");
        return 1;
    }
    
    printf("Basic Trajectory Results:\n");
    printf("-------------------------\n");
    printf("Max Range: %.2f m\n", result->max_range);
    printf("Max Height: %.2f m\n", result->max_height);
    printf("Time of Flight: %.3f s\n", result->time_of_flight);
    printf("Impact Velocity: %.2f m/s\n", result->impact_velocity);
    printf("Impact Energy: %.2f J\n", result->impact_energy);
    printf("Point Count: %d\n", result->point_count);
    printf("\n");
    
    // Check advanced features
    printf("Advanced Physics Results:\n");
    printf("------------------------\n");
    
    if (!isnan(result->min_pitch_damping)) {
        printf("Min Pitch Damping: %.3f\n", result->min_pitch_damping);
        if (result->min_pitch_damping > 0) {
            printf("  Warning: Transonic instability detected!\n");
        }
    }
    
    if (!isnan(result->transonic_mach)) {
        printf("Transonic Entry Mach: %.2f\n", result->transonic_mach);
    }
    
    if (!isnan(result->final_pitch_angle)) {
        printf("Final Pitch Angle: %.4f rad (%.2f°)\n", 
               result->final_pitch_angle, 
               result->final_pitch_angle * 180.0 / 3.14159265359);
    }
    
    if (!isnan(result->final_yaw_angle)) {
        printf("Final Yaw Angle: %.4f rad (%.2f°)\n", 
               result->final_yaw_angle,
               result->final_yaw_angle * 180.0 / 3.14159265359);
    }
    
    if (!isnan(result->max_yaw_angle)) {
        printf("Max Yaw Angle: %.4f rad (%.2f°)\n", 
               result->max_yaw_angle,
               result->max_yaw_angle * 180.0 / 3.14159265359);
    }
    
    if (!isnan(result->max_precession_angle)) {
        printf("Max Precession: %.4f rad\n", result->max_precession_angle);
    }
    
    printf("\n");
    
    // Check sampled points
    if (result->sampled_point_count > 0 && result->sampled_points != NULL) {
        printf("Trajectory Samples (%d points):\n", result->sampled_point_count);
        printf("Distance(m)  Time(s)  Velocity(m/s)  Drop(m)  Windage(m)\n");
        printf("---------------------------------------------------------\n");
        
        for (int i = 0; i < result->sampled_point_count && i < 10; i++) {
            FFITrajectorySample* sample = &result->sampled_points[i];
            printf("%10.1f  %7.3f  %12.1f  %7.3f  %9.3f\n",
                   sample->distance, sample->time, sample->velocity_mps,
                   sample->drop_meters, sample->windage_meters);
        }
        if (result->sampled_point_count > 10) {
            printf("... (%d more samples)\n", result->sampled_point_count - 10);
        }
    }
    
    // Clean up
    ballistics_free_trajectory_result(result);
    
    printf("\n✓ FFI test completed successfully!\n");
    return 0;
}