#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
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
    int use_adaptive_rk45;
    int enable_wind_shear;
    int enable_trajectory_sampling;
    double sample_interval;
    int enable_pitch_damping;
    int enable_precession_nutation;
    int enable_spin_drift;
    int enable_magnus;
    int enable_coriolis;
    double shot_azimuth;
    double cant_angle;      // appended (must match FFIBallisticInputs field order): rifle cant angle, radians
    double zero_poi_vertical;   // appended (must match FFIBallisticInputs field order): deliberate vertical POI offset at the zero range, meters, + = impacts high (MBA-1359)
    double zero_poi_horizontal; // appended (must match FFIBallisticInputs field order): deliberate horizontal POI offset at the zero range, meters, + = impacts right (MBA-1359)
    double sight_offset_lateral;// appended (must match FFIBallisticInputs field order): lateral sight-to-bore mount offset, meters, + = sight right of bore (MBA-1396)
} FFIBallisticInputs;

typedef struct {
    double speed;
    double direction;
    double vertical_speed;  // appended (must match FFIWindConditions field order): vertical wind m/s, + = updraft
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

// MBA-1361: appended reticle hold-point export. A NEW struct and a NEW function -- no
// pre-existing layout above is affected, so this mirror can be added without touching
// anything else in this file.
typedef struct {
    double down_mil;
    double right_mil;
    int nearest_mark;
    double nearest_mark_distance_mil;
    int off_reticle;
    double mark_scale;
} FFIReticleHold;

#define FFI_RETICLE_FIRST_FOCAL_PLANE 0
#define FFI_RETICLE_SECOND_FOCAL_PLANE 1
#define FFI_RETICLE_OK 0
#define FFI_RETICLE_ERR_INVALID_ARGUMENT (-1)
#define FFI_RETICLE_ERR_MAGNIFICATION (-2)

// External FFI functions
extern FFITrajectoryResult* ballistics_calculate_trajectory(
    const FFIBallisticInputs* inputs,
    const FFIWindConditions* wind,
    const FFIAtmosphericConditions* atmosphere,
    double max_range,
    double step_size
);

extern void ballistics_free_trajectory_result(FFITrajectoryResult* result);

// 0.49.0 station-presence exports: the same structs, plus a bit set saying which of the
// atmosphere's temperature / pressure the caller actually supplied.
#define BALLISTICS_STATION_TEMPERATURE_SUPPLIED 1u
#define BALLISTICS_STATION_PRESSURE_SUPPLIED 2u

typedef struct {
    int num_simulations;
    double velocity_std_dev;
    double angle_std_dev;
    double bc_std_dev;
    double wind_speed_std_dev;
    double target_distance;
    double base_wind_speed;
    double base_wind_direction;
    double azimuth_std_dev;
} FFIMonteCarloParams;

typedef struct {
    double* ranges;
    double* impact_velocities;
    double* impact_positions_x;
    double* impact_positions_y;
    double* impact_positions_z;
    int num_results;
    double mean_range;
    double std_dev_range;
    double mean_impact_velocity;
    double std_dev_impact_velocity;
    double hit_probability;
} FFIMonteCarloResults;

extern double ballistics_calculate_zero_angle(
    const FFIBallisticInputs* inputs,
    const FFIWindConditions* wind,
    const FFIAtmosphericConditions* atmosphere,
    double zero_distance
);

extern FFITrajectoryResult* ballistics_calculate_trajectory_with_station_presence(
    const FFIBallisticInputs* inputs,
    const FFIWindConditions* wind,
    const FFIAtmosphericConditions* atmosphere,
    double max_range,
    double step_size,
    const double* drag_mach,
    const double* drag_cd,
    int drag_table_len,
    double cd_scale,
    uint32_t station_presence
);

extern double ballistics_calculate_zero_angle_with_station_presence(
    const FFIBallisticInputs* inputs,
    const FFIWindConditions* wind,
    const FFIAtmosphericConditions* atmosphere,
    double zero_distance,
    const double* drag_mach,
    const double* drag_cd,
    int drag_table_len,
    double cd_scale,
    uint32_t station_presence
);

extern FFIMonteCarloResults* ballistics_monte_carlo(
    const FFIBallisticInputs* inputs,
    const FFIAtmosphericConditions* atmosphere,
    const FFIMonteCarloParams* params
);

extern FFIMonteCarloResults* ballistics_monte_carlo_with_station_presence(
    const FFIBallisticInputs* inputs,
    const FFIAtmosphericConditions* atmosphere,
    const FFIMonteCarloParams* params,
    double wind_direction_std_dev,
    uint32_t station_presence
);

extern void ballistics_free_monte_carlo_results(FFIMonteCarloResults* results);

extern int ballistics_hold_point_in_reticle(
    double drop_mil,
    double wind_mil,
    double magnification,
    const double* marks,       // flat [down_0, right_0, down_1, right_1, ...], nominal mil
    int marks_len,             // number of MARKS (half the array length)
    int focal_plane,
    double reference_magnification,
    FFIReticleHold* out
);

// MBA-1361 smoke test: FFP invariance, SFP scaling, and the marks_len bounds guard.
static int test_reticle_hold(void) {
    // center, 2 mil, 4 mil, and one windage dot at (2, 1).
    const double marks[8] = {0.0, 0.0, 2.0, 0.0, 4.0, 0.0, 2.0, 1.0};
    FFIReticleHold hold;
    int rc;

    printf("Reticle Hold Point (MBA-1361):\n");
    printf("------------------------------\n");

    rc = ballistics_hold_point_in_reticle(4.0, 0.0, 6.0, marks, 4,
                                          FFI_RETICLE_FIRST_FOCAL_PLANE, 0.0, &hold);
    if (rc != FFI_RETICLE_OK || hold.nearest_mark != 2 ||
        fabs(hold.nearest_mark_distance_mil) > 1e-12 || hold.mark_scale != 1.0) {
        printf("FATAL: FFP hold mismatch (rc=%d mark=%d dist=%.6f scale=%.6f)\n",
               rc, hold.nearest_mark, hold.nearest_mark_distance_mil, hold.mark_scale);
        return 1;
    }
    printf("FFP  4.00 mil drop -> mark #%d, %.3f mil away, scale %.2fx\n",
           hold.nearest_mark, hold.nearest_mark_distance_mil, hold.mark_scale);

    // SFP at half the reference magnification: the etched 2 mil mark reads 4 mil true.
    rc = ballistics_hold_point_in_reticle(4.0, 0.0, 5.0, marks, 4,
                                          FFI_RETICLE_SECOND_FOCAL_PLANE, 10.0, &hold);
    if (rc != FFI_RETICLE_OK || hold.nearest_mark != 1 ||
        fabs(hold.nearest_mark_distance_mil) > 1e-12 || hold.mark_scale != 2.0) {
        printf("FATAL: SFP hold mismatch (rc=%d mark=%d dist=%.6f scale=%.6f)\n",
               rc, hold.nearest_mark, hold.nearest_mark_distance_mil, hold.mark_scale);
        return 1;
    }
    printf("SFP  4.00 mil drop at 5x (ref 10x) -> mark #%d, %.3f mil away, scale %.2fx\n",
           hold.nearest_mark, hold.nearest_mark_distance_mil, hold.mark_scale);

    // The marks_len guard (MBA-1407 lesson) must reject before reading anything.
    if (ballistics_hold_point_in_reticle(1.0, 0.0, 10.0, marks, 0,
                                         FFI_RETICLE_FIRST_FOCAL_PLANE, 0.0, &hold)
        != FFI_RETICLE_ERR_INVALID_ARGUMENT) {
        printf("FATAL: marks_len 0 was not rejected\n");
        return 1;
    }
    if (ballistics_hold_point_in_reticle(1.0, 0.0, 10.0, marks, 1 << 30,
                                         FFI_RETICLE_FIRST_FOCAL_PLANE, 0.0, &hold)
        != FFI_RETICLE_ERR_INVALID_ARGUMENT) {
        printf("FATAL: oversized marks_len was not rejected\n");
        return 1;
    }
    if (ballistics_hold_point_in_reticle(1.0, 0.0, 0.0, marks, 4,
                                         FFI_RETICLE_FIRST_FOCAL_PLANE, 0.0, &hold)
        != FFI_RETICLE_ERR_MAGNIFICATION) {
        printf("FATAL: zero magnification was not rejected\n");
        return 1;
    }
    printf("Bounds and magnification guards: OK\n\n");
    return 0;
}

// 0.49.0: a pressure inside the old 1013.25 +/- 0.5 hPa band at 150 m. Flags 0 must equal the
// struct export (the legacy rule replaces it with ICAO air); the pressure bit must fly it.
static int test_station_presence(const FFIBallisticInputs* base) {
    FFIBallisticInputs inputs = *base;
    inputs.altitude = 150.0;
    inputs.temperature = 20.0;
    inputs.enable_wind_shear = 0;
    FFIAtmosphericConditions atmosphere = {
        .temperature = 20.0,
        .pressure = 1013.0,
        .humidity = 50.0,
        .altitude = 150.0
    };

    printf("Station presence (0.49.0):\n");
    printf("--------------------------\n");

    double legacy = ballistics_calculate_zero_angle(&inputs, NULL, &atmosphere, 300.0);
    double none = ballistics_calculate_zero_angle_with_station_presence(
        &inputs, NULL, &atmosphere, 300.0, NULL, NULL, 0, 1.0, 0u);
    double supplied = ballistics_calculate_zero_angle_with_station_presence(
        &inputs, NULL, &atmosphere, 300.0, NULL, NULL, 0, 1.0,
        BALLISTICS_STATION_PRESSURE_SUPPLIED);
    double unknown = ballistics_calculate_zero_angle_with_station_presence(
        &inputs, NULL, &atmosphere, 300.0, NULL, NULL, 0, 1.0, 4u);
    if (!isfinite(legacy) || none != legacy || !isfinite(supplied) || supplied == legacy ||
        !isnan(unknown)) {
        printf("FATAL: zero angle presence mismatch (legacy %.17g, flags 0 %.17g, flags 2 %.17g)\n",
               legacy, none, supplied);
        return 1;
    }
    printf("Zero angle: legacy %.9f rad, flags 0 equal, pressure supplied %.9f rad\n",
           legacy, supplied);

    FFITrajectoryResult* old_result = ballistics_calculate_trajectory(
        &inputs, NULL, &atmosphere, 600.0, 1.0);
    FFITrajectoryResult* zero_bits = ballistics_calculate_trajectory_with_station_presence(
        &inputs, NULL, &atmosphere, 600.0, 1.0, NULL, NULL, 0, 1.0, 0u);
    FFITrajectoryResult* typed = ballistics_calculate_trajectory_with_station_presence(
        &inputs, NULL, &atmosphere, 600.0, 1.0, NULL, NULL, 0, 1.0,
        BALLISTICS_STATION_PRESSURE_SUPPLIED);
    int bad = old_result == NULL || zero_bits == NULL || typed == NULL ||
              zero_bits->impact_velocity != old_result->impact_velocity ||
              typed->impact_velocity == old_result->impact_velocity;
    if (!bad) {
        printf("Trajectory: legacy %.6f m/s, flags 0 equal, pressure supplied %.6f m/s\n",
               old_result->impact_velocity, typed->impact_velocity);
    }
    ballistics_free_trajectory_result(old_result);
    ballistics_free_trajectory_result(zero_bits);
    ballistics_free_trajectory_result(typed);
    if (bad) {
        printf("FATAL: trajectory presence mismatch\n");
        return 1;
    }

    // No dispersion at all: every trial is the baseline, so the unseeded runs compare exactly.
    FFIMonteCarloParams params = {
        .num_simulations = 4,
        .velocity_std_dev = 0.0,
        .angle_std_dev = 0.0,
        .bc_std_dev = 0.0,
        .wind_speed_std_dev = 0.0,
        .target_distance = 300.0,
        .base_wind_speed = 0.0,
        .base_wind_direction = 0.0,
        .azimuth_std_dev = 0.0
    };
    FFIMonteCarloResults* mc_old = ballistics_monte_carlo(&inputs, &atmosphere, &params);
    FFIMonteCarloResults* mc_zero = ballistics_monte_carlo_with_station_presence(
        &inputs, &atmosphere, &params, 0.0, 0u);
    FFIMonteCarloResults* mc_typed = ballistics_monte_carlo_with_station_presence(
        &inputs, &atmosphere, &params, 0.0, BALLISTICS_STATION_PRESSURE_SUPPLIED);
    bad = mc_old == NULL || mc_zero == NULL || mc_typed == NULL ||
          mc_zero->mean_impact_velocity != mc_old->mean_impact_velocity ||
          mc_typed->mean_impact_velocity == mc_old->mean_impact_velocity;
    if (!bad) {
        printf("Monte Carlo: legacy %.6f m/s, flags 0 equal, pressure supplied %.6f m/s\n\n",
               mc_old->mean_impact_velocity, mc_typed->mean_impact_velocity);
    }
    ballistics_free_monte_carlo_results(mc_old);
    ballistics_free_monte_carlo_results(mc_zero);
    ballistics_free_monte_carlo_results(mc_typed);
    if (bad) {
        printf("FATAL: Monte Carlo presence mismatch\n");
        return 1;
    }
    return 0;
}

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
        .use_adaptive_rk45 = 0,         // (must match FFIBallisticInputs field order)
        .enable_wind_shear = 1,         // Enable wind shear
        .enable_trajectory_sampling = 1, // Enable sampling
        .sample_interval = 50.0,        // Sample every 50 meters
        .enable_pitch_damping = 1,      // Enable pitch damping
        .enable_precession_nutation = 1,// Enable precession/nutation
        .enable_spin_drift = 1,         // Enable spin drift
        .enable_magnus = 1,             // Enable Magnus effect
        .enable_coriolis = 0,           // Disable Coriolis for now
        .shot_azimuth = 0.0,            // North; appended ABI field
        .cant_angle = 0.0,              // No cant; appended ABI field
        .zero_poi_vertical = 0.0,       // No deliberate zero POI offset; appended ABI field (MBA-1359)
        .zero_poi_horizontal = 0.0,     // No deliberate zero POI offset; appended ABI field (MBA-1359)
        .sight_offset_lateral = 0.0     // Sight directly above the bore; appended ABI field (MBA-1396)
    };

    FFIWindConditions wind = {
        .speed = 5.0,                   // 5 m/s wind
        .direction = 1.5708,            // 90 degrees (from right)
        .vertical_speed = 0.0           // No updraft/downdraft; appended ABI field
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

    // Sanity gate: with struct mirrors out of sync with src/ffi.rs, the C caller and
    // the Rust callee disagree on field layout/struct size, which is undefined
    // behavior. That UB has previously manifested as a wildly wrong-but-still-"valid"
    // -looking max_height (e.g. 890 m collapsing to 73 m under -fstack-protector-all)
    // rather than a crash, letting a silently broken fixture "pass". Fail loudly
    // instead of printing garbage as if it were a real result.
    if (!isfinite(result->max_height) || result->max_height <= 0.0 || result->max_height >= 10000.0) {
        printf("FATAL: implausible Max Height %.2f m (expected 0 < h < 10000) -- likely FFI struct/ABI mismatch\n",
               result->max_height);
        ballistics_free_trajectory_result(result);
        exit(1);
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

    printf("\n");
    if (test_reticle_hold() != 0) {
        return 1;
    }
    if (test_station_presence(&inputs) != 0) {
        return 1;
    }

    printf("\n✓ FFI test completed successfully!\n");
    return 0;
}
