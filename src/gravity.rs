//! Internal uniform-gravity model. Public entry points select standard Earth gravity.
//!
//! Gravity is per simulation, never mutable global state. This does not select an atmosphere,
//! planetary radius, rotation rate, or a new calibration for empirical aerodynamic models.

use nalgebra::Vector3;

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Gravity {
    magnitude_mps2: f64,
}

impl Gravity {
    pub(crate) const EARTH: Self = match Self::new(crate::constants::G_ACCEL_MPS2) {
        Some(gravity) => gravity,
        None => panic!("invalid standard Earth gravity"),
    };

    /// A finite, nonnegative acceleration magnitude in m/s²; zero is valid.
    pub(crate) const fn new(magnitude_mps2: f64) -> Option<Self> {
        if magnitude_mps2.is_finite() && magnitude_mps2 >= 0.0 {
            Some(Self { magnitude_mps2 })
        } else {
            None
        }
    }

    pub(crate) fn magnitude(self) -> f64 {
        self.magnitude_mps2
    }

    /// Project downward acceleration into a frame inclined above the horizontal.
    pub(crate) fn in_frame(self, inclination_rad: f64) -> Vector3<f64> {
        Vector3::new(
            -self.magnitude_mps2 * inclination_rad.sin(),
            -self.magnitude_mps2 * inclination_rad.cos(),
            0.0,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_invalid_magnitudes_and_accepts_zero() {
        for magnitude in [-1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(Gravity::new(magnitude).is_none());
        }
        assert_eq!(Gravity::new(0.0).unwrap().in_frame(0.7), Vector3::zeros());
    }

    #[test]
    fn projection_preserves_magnitude_and_earth_components() {
        for angle in [-1.2_f64, 0.0, 0.7] {
            let acceleration = Gravity::EARTH.in_frame(angle);
            assert_eq!(acceleration.x, -9.80665 * angle.sin());
            assert_eq!(acceleration.y, -9.80665 * angle.cos());
            assert!((Gravity::new(1.62).unwrap().in_frame(angle).norm() - 1.62).abs() < 1e-14);
        }
    }
}
