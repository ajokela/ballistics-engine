//! The wind-field API as a downstream crate sees it: public paths, trait objects, the
//! `WindSock` adapter and the shared direction convention.

use ballistics_engine::wind::{
    wind_vector, ConstantField, EkmanProfile, LogProfile, ProfileField, ProfileLevel, WindField,
    WindFieldError, WindSegment, WindSock,
};
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
