//! The wind-field API as a downstream crate sees it: public paths, trait objects, the
//! `WindSock` adapter and the shared direction convention.

use ballistics_engine::wind::{
    wind_vector, ConstantField, EkmanProfile, HeightOffset, LogProfile, ProfileField, ProfileLevel,
    WindField, WindFieldError, WindSegment, WindSock,
};
use ballistics_engine::{AtmosphericConditions, BallisticInputs, TrajectorySolver, WindConditions};
use nalgebra::Vector3;
use std::f64::consts::FRAC_PI_2;

fn level(altitude_m: f64, speed_mps: f64, direction_from_deg: f64) -> ProfileLevel {
    ProfileLevel {
        altitude_m,
        speed_mps,
        direction_from_deg,
    }
}

#[test]
fn every_field_type_works_behind_one_trait_object() {
    let sock = WindSock::try_new(vec![
        WindSegment::new(18.0, 90.0, 200.0),
        WindSegment::new(36.0, 0.0, 800.0),
    ])
    .unwrap();
    let fields: Vec<Box<dyn WindField>> = vec![
        Box::new(ConstantField(Vector3::new(-3.0, 0.0, 1.0))),
        Box::new(ProfileField::new(vec![level(0.0, 2.0, 90.0), level(500.0, 8.0, 120.0)]).unwrap()),
        Box::new(LogProfile::new(6.0, 10.0, 0.1, 270.0).unwrap()),
        Box::new(EkmanProfile::new(15.0, 52.0, 10.0, 240.0).unwrap()),
        Box::new(sock.clone()),
    ];
    let p = Vector3::new(300.0, 120.0, -2.0);
    for field in &fields {
        let w = field.wind(1.5, p);
        assert!(w.iter().all(|c| c.is_finite()));
        assert_eq!(field.wind(1.5, p), w, "a field must be deterministic");
    }
    assert_eq!(
        fields[4].wind(0.0, p),
        sock.vector_for_range_stateless(300.0)
    );
}

#[test]
fn profile_levels_follow_the_wind_vector_convention() {
    let at0 = ProfileField::new(vec![level(0.0, 7.0, 0.0)]).unwrap();
    assert_eq!(at0.wind(0.0, Vector3::zeros()), wind_vector(7.0, 0.0, 0.0));
    let at90 = ProfileField::new(vec![level(0.0, 7.0, 90.0)]).unwrap();
    assert_eq!(
        at90.wind(0.0, Vector3::zeros()),
        wind_vector(7.0, FRAC_PI_2, 0.0)
    );
}

#[test]
fn constructors_return_typed_errors() {
    assert_eq!(ProfileField::new(vec![]), Err(WindFieldError::EmptyProfile));
    let err = LogProfile::new(6.0, 10.0, 10.0, 0.0).unwrap_err();
    assert_eq!(
        err.to_string(),
        "roughness_m must be finite, greater than zero and below the reference height"
    );
    let err = EkmanProfile::new(15.0, 0.0, 10.0, 0.0).unwrap_err();
    assert_eq!(
        err.to_string(),
        "latitude_deg must be finite, nonzero and strictly between -90 and 90 degrees"
    );
}

#[test]
fn a_launch_frame_reads_a_sea_level_profile_through_the_site_elevation() {
    // A sounding-style profile in metres above sea level and a launch site 1609 m up.
    let msl =
        ProfileField::new(vec![level(1609.0, 2.0, 180.0), level(3609.0, 12.0, 270.0)]).unwrap();
    let launch_frame: Box<dyn WindField> =
        Box::new(HeightOffset::new(msl.clone(), 1609.0).unwrap());
    assert_eq!(
        launch_frame.wind(0.0, Vector3::zeros()),
        msl.wind_at(1609.0)
    );
    assert_eq!(
        launch_frame.wind(0.0, Vector3::new(0.0, 1000.0, 0.0)),
        msl.wind_at(2609.0)
    );
}

/// The height of the first trajectory point for a given `muzzle_height`.
fn launch_height_m(muzzle_height: f64) -> f64 {
    let inputs = BallisticInputs {
        muzzle_velocity: 800.0,
        muzzle_angle: 0.01,
        target_distance: 100.0,
        muzzle_height,
        ..BallisticInputs::default()
    };
    let mut solver = TrajectorySolver::new(
        inputs,
        WindConditions::default(),
        AtmosphericConditions::default(),
    );
    solver.set_max_range(50.0);
    solver.solve().unwrap().points[0].position.y
}

#[test]
fn a_solver_frame_offsets_a_log_profile_by_only_the_bore_height_it_leaves_out() {
    // On a level shot TrajectorySolver measures y from the level muzzle_height is measured from,
    // so the offset for a ground-relative profile is the real bore height minus muzzle_height.
    let bore_height_m = 1.5;
    let log = LogProfile::new(8.0, 10.0, 0.03, 90.0).unwrap();

    // muzzle_height is the real bore height: y is already height above the ground.
    let y = launch_height_m(bore_height_m);
    assert_eq!(y, bore_height_m);
    assert_eq!(
        log.wind(0.0, Vector3::new(0.0, y, 0.0)),
        log.wind_at(bore_height_m)
    );

    // muzzle_height left at its default 0: y = 0 is the bore, so the profile takes +1.5.
    assert_eq!(BallisticInputs::default().muzzle_height, 0.0);
    let y = launch_height_m(0.0);
    assert_eq!(y, 0.0);
    let from_bore = HeightOffset::new(log, bore_height_m - 0.0).unwrap();
    assert_eq!(
        from_bore.wind(0.0, Vector3::new(0.0, y, 0.0)),
        log.wind_at(bore_height_m)
    );
}
