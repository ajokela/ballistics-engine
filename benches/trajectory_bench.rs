use ballistics_engine::*;
use criterion::{criterion_group, criterion_main, Criterion};
use std::hint::black_box;

fn bench_trajectory_g7(c: &mut Criterion) {
    c.bench_function("trajectory_g7_1000yd", |b| {
        b.iter(|| {
            let inputs = BallisticInputs {
                bc_value: 0.275,
                bc_type: DragModel::G7,
                bullet_mass: 0.01088,       // 168gr in kg
                muzzle_velocity: 822.96,    // 2700fps in m/s
                bullet_diameter: 0.00782,   // .308 in meters
                target_distance: 914.4,     // 1000yd in meters
                muzzle_angle: 0.005,        // small launch angle
                sight_height: 0.0508,       // 2 inches
                use_adaptive_rk45: true,
                use_rk4: true,
                ..BallisticInputs::default()
            };
            let wind = WindConditions { speed: 0.0, direction: 0.0, vertical_speed: 0.0 };
            let atmo = AtmosphericConditions {
                temperature: 15.0, pressure: 1013.25, humidity: 50.0, altitude: 0.0,
            };
            let solver = TrajectorySolver::new(black_box(inputs), wind, atmo);
            solver.solve()
        })
    });
}

fn bench_zero_finding(c: &mut Criterion) {
    c.bench_function("zero_100yd", |b| {
        b.iter(|| {
            let inputs = BallisticInputs {
                bc_value: 0.475,
                bc_type: DragModel::G1,
                muzzle_velocity: 822.96, // 2700fps in m/s
                sight_height: 0.0508,    // 2 inches
                ..BallisticInputs::default()
            };
            calculate_zero_angle(
                black_box(inputs),
                black_box(91.44), // 100yd in meters
                0.0, // target height
            )
        })
    });
}

criterion_group!(benches, bench_trajectory_g7, bench_zero_finding);
criterion_main!(benches);
