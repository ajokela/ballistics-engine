use crate::InternalBallisticInputs as BallisticInputs;

/// Calculate the gyroscopic stability coefficient (SG) for the bullet.
///
/// This function uses the Miller stability formula. An SG value greater than 1.5
/// is generally considered to indicate adequate stability.
///
/// # Arguments
/// * `inputs` - Ballistic input parameters
/// * `atmo_params` - Atmospheric parameters (altitude, temp_c, pressure_hpa, density_ratio)
///
/// # Returns
/// * Stability coefficient (dimensionless)
pub fn compute_stability_coefficient(
    inputs: &BallisticInputs,
    atmo_params: (f64, f64, f64, f64),
) -> f64 {
    // Check for required parameters
    if inputs.twist_rate == 0.0 || inputs.bullet_length == 0.0 || inputs.bullet_diameter == 0.0 {
        return 0.0;
    }

    // Pre-calculated constants for efficiency
    const MILLER_CONST: f64 = 30.0;
    const VEL_REF_FPS: f64 = 2800.0;
    const TEMP_REF_K: f64 = 288.15; // 15°C
    const PRESS_REF_HPA: f64 = 1013.25;

    // Calculate intermediate values
    // Convert twist rate from inches to meters for consistency
    let twist_rate_m = inputs.twist_rate.abs() * 0.0254; // inches to meters
    let twist_calibers = twist_rate_m / inputs.bullet_diameter;
    let length_calibers = inputs.bullet_length / inputs.bullet_diameter;

    // Convert units for Miller formula
    let mass_grains = inputs.bullet_mass / 0.00006479891; // kg to grains
    let diameter_inches = inputs.bullet_diameter / 0.0254; // meters to inches
    let velocity_fps = inputs.muzzle_velocity * 3.28084; // m/s to fps

    // Miller stability formula components
    let mass_term = MILLER_CONST * mass_grains;
    let geom_term = twist_calibers.powi(2)
        * diameter_inches.powi(3)
        * length_calibers
        * (1.0 + length_calibers.powi(2));

    if geom_term == 0.0 {
        return 0.0;
    }

    // Extract atmospheric parameters
    let (_, temp_c, current_press_hpa, _) = atmo_params;
    let temp_k = temp_c + 273.15;

    // Atmospheric density correction factor
    // Ratio of reference density to current density
    let density_correction = (temp_k / TEMP_REF_K) * (PRESS_REF_HPA / current_press_hpa);

    // Velocity correction factor
    let velocity_correction = (velocity_fps / VEL_REF_FPS).powf(1.0 / 3.0);

    // Final stability calculation

    (mass_term / geom_term) * velocity_correction * density_correction
}

/// Effective density (kg/m^3) used by the mass-based bullet-length estimate.
///
/// A jacketed lead-core bullet is not a solid cylinder of lead (11340 kg/m^3): the ogive/boat-tail
/// taper and the copper jacket lower the *effective* density that relates total mass to a simple
/// `length * frontal-area` cylinder. 7600 kg/m^3 was fit so common lead-core match bullets land at
/// their real length / caliber ratio (see reference-bullet unit tests below). MBA-1135.
pub const BULLET_LENGTH_RHO_EFF_KG_M3: f64 = 7600.0;

/// Lower / upper bound on the estimated length-to-diameter ratio (calibers). Keeps degenerate
/// mass/diameter combinations inside the range of real small-arms projectiles.
const MIN_LENGTH_CALIBERS: f64 = 2.5;
const MAX_LENGTH_CALIBERS: f64 = 6.5;

/// Gyroscopic-stability target used to synthesize a default twist when the shooter does not
/// supply one (MBA-1135).
///
/// Chosen = **2.0** (see [`default_twist_inches`] docs for the reference-bullet comparison that
/// motivated 2.0 over 1.5). A synthesized default should slightly *over*-stabilize rather than
/// under-stabilize, and 2.0 reproduces the real factory twists of the twist-sensitive high-BC
/// match bullets (.308/175 gr -> ~1:11", 6.5mm/140 gr -> ~1:8") most faithfully.
pub const DEFAULT_TWIST_SG_TARGET: f64 = 2.0;

const GRAINS_PER_KG: f64 = 1.0 / 0.00006479891;
const METERS_PER_INCH: f64 = 0.0254;
const MPS_TO_FPS: f64 = 3.28084;
const MILLER_VEL_REF_FPS: f64 = 2800.0;

/// Estimate a bullet's length (meters) from its diameter and mass using a constant-effective-density
/// cylinder model (MBA-1135).
///
/// Replaces the historical mass-blind `diameter * 4.5` fallback: two bullets of the same caliber but
/// very different weights get very different lengths, which is what the Miller stability formula and
/// the enhanced spin-drift / Magnus models actually need.
///
/// Model: `L = mass / (rho_eff * (pi/4) * d^2)` with `rho_eff = 7600 kg/m^3`, then the length/diameter
/// ratio is clamped to `[2.5, 6.5]` calibers so pathological inputs stay physical.
///
/// Returns `0.0` for non-physical inputs (`mass_kg <= 0` or `diameter_m <= 0`); callers that need a
/// non-zero length should fall back to their historical literal in that case.
pub fn estimate_bullet_length_m(diameter_m: f64, mass_kg: f64) -> f64 {
    if mass_kg <= 0.0 || diameter_m <= 0.0 || !mass_kg.is_finite() || !diameter_m.is_finite() {
        return 0.0;
    }

    let frontal_area = std::f64::consts::FRAC_PI_4 * diameter_m * diameter_m;
    let raw_length = mass_kg / (BULLET_LENGTH_RHO_EFF_KG_M3 * frontal_area);

    // Clamp the length-to-diameter ratio (calibers) into a physically sensible band.
    let length_calibers = (raw_length / diameter_m).clamp(MIN_LENGTH_CALIBERS, MAX_LENGTH_CALIBERS);
    length_calibers * diameter_m
}

/// Synthesize a default barrel twist (inches per turn) for a bullet whose twist the shooter did not
/// supply, by inverting the Miller stability formula for the twist that yields
/// [`DEFAULT_TWIST_SG_TARGET`] (MBA-1135).
///
/// This replaces the mass-blind fixed `1:12"` fallback: heavier / longer bullets correctly demand a
/// faster twist. The bullet length is taken from [`estimate_bullet_length_m`] (so the whole default
/// is caliber + weight aware), the velocity term is `(v_fps / 2800)^(1/3)` at the supplied muzzle
/// velocity, and the density correction is `1.0` (sea-level ICAO standard).
///
/// Inverting `Sg = 30 m_gr * vel_corr / (t_cal^2 * d_in^3 * L_cal * (1 + L_cal^2))`:
///
/// `t_cal = sqrt( 30 m_gr * vel_corr / (Sg_target * d_in^3 * L_cal * (1 + L_cal^2)) )`,
/// `twist_in = t_cal * d_in`.
///
/// `Sg_target` was selected by comparing 1.5 vs 2.0 against common factory twists
/// (velocities: .308/175 gr @ 2600 fps, .223/55 gr @ 3240 fps, .224/77 gr @ 2750 fps,
/// 6.5mm/140 gr @ 2700 fps):
///
/// | bullet          | factory | Sg=1.5 | Sg=2.0 |
/// |-----------------|---------|--------|--------|
/// | .308 / 175 gr   | 1:10-12"| 1:12.9"| 1:11.2"|
/// | .223 / 55 gr    | ~1:12"  | 1:11.7"| 1:10.1"|
/// | .224 / 77 gr    | ~1:8"   | 1:8.4" | 1:7.2" |
/// | 6.5mm / 140 gr  | ~1:8"   | 1:8.9" | 1:7.7" |
///
/// 2.0 wins: Sg=1.5 gives 1:12.9" for the 175 gr match bullet (no real .308 barrel is that slow and
/// it under-stabilizes 175 gr), whereas 2.0's 1:11.2" / 1:8" are exactly the standard twists for the
/// twist-sensitive high-BC bullets, and slight over-stabilization is the safe direction for a default.
///
/// Guards degenerate inputs (`diameter_m`, `mass_kg`, or `muzzle_velocity_mps` <= 0) by returning the
/// historical `12.0` (1:12") fallback. Always returns a positive inches-per-turn value.
pub fn default_twist_inches(diameter_m: f64, mass_kg: f64, muzzle_velocity_mps: f64) -> f64 {
    const FALLBACK_TWIST_IN: f64 = 12.0;

    if diameter_m <= 0.0
        || mass_kg <= 0.0
        || muzzle_velocity_mps <= 0.0
        || !diameter_m.is_finite()
        || !mass_kg.is_finite()
        || !muzzle_velocity_mps.is_finite()
    {
        return FALLBACK_TWIST_IN;
    }

    let length_m = estimate_bullet_length_m(diameter_m, mass_kg);
    if length_m <= 0.0 {
        return FALLBACK_TWIST_IN;
    }

    let d_in = diameter_m / METERS_PER_INCH;
    let m_gr = mass_kg * GRAINS_PER_KG;
    let l_cal = length_m / diameter_m;
    let v_fps = muzzle_velocity_mps * MPS_TO_FPS;
    let velocity_correction = (v_fps / MILLER_VEL_REF_FPS).powf(1.0 / 3.0);

    let geom = d_in.powi(3) * l_cal * (1.0 + l_cal * l_cal);
    let denom = DEFAULT_TWIST_SG_TARGET * geom;
    if denom <= 0.0 {
        return FALLBACK_TWIST_IN;
    }

    let t_cal_sq = 30.0 * m_gr * velocity_correction / denom;
    if !t_cal_sq.is_finite() || t_cal_sq <= 0.0 {
        return FALLBACK_TWIST_IN;
    }

    let twist_in = t_cal_sq.sqrt() * d_in;
    if twist_in.is_finite() && twist_in > 0.0 {
        twist_in
    } else {
        FALLBACK_TWIST_IN
    }
}

/// Preserve a supplied twist (already converted to inches per turn), or synthesize the
/// caliber/weight/velocity-aware default used by frontend argument parsers.
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) fn resolve_twist_inches(
    explicit_twist_inches: Option<f64>,
    diameter_m: f64,
    mass_kg: f64,
    muzzle_velocity_mps: f64,
) -> f64 {
    explicit_twist_inches
        .unwrap_or_else(|| default_twist_inches(diameter_m, mass_kg, muzzle_velocity_mps))
}

/// Calculate spin drift in meters using Litz approximation.
///
/// # Arguments
/// * `time_s` - Time of flight in seconds
/// * `stability` - Stability coefficient
/// * `twist_rate` - Twist rate in inches (calibers per turn)
/// * `is_twist_right` - True for right-hand twist, false for left-hand
///
/// # Returns
/// * Spin drift in meters
pub fn compute_spin_drift(
    time_s: f64,
    stability: f64,
    twist_rate: f64,
    is_twist_right: bool,
) -> f64 {
    compute_spin_drift_with_decay(time_s, stability, twist_rate, is_twist_right, None)
}

/// Calculate spin drift with optional spin decay modeling
pub fn compute_spin_drift_with_decay(
    time_s: f64,
    stability: f64,
    twist_rate: f64,
    is_twist_right: bool,
    decay_factor: Option<f64>, // Optional spin decay factor (0-1)
) -> f64 {
    if stability == 0.0 || time_s <= 0.0 || twist_rate == 0.0 {
        return 0.0;
    }

    let sign = if is_twist_right { 1.0 } else { -1.0 };

    // Litz empirical spin drift: inches = 1.25 * (SG + 1.2) * TOF^1.83.
    // Keep this summary/API path consistent with cli_api::apply_spin_drift.
    let scaling_factor = 1.25;
    let base_drift = sign * scaling_factor * (stability + 1.2) * time_s.powf(1.83);

    // Apply spin decay if provided
    let effective_drift = if let Some(decay) = decay_factor {
        base_drift * decay.max(0.0).min(1.0)
    } else {
        base_drift
    };

    // Convert inches to meters
    effective_drift * 0.0254
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_inputs() -> BallisticInputs {
        BallisticInputs {
            muzzle_velocity: 823.0, // 2700 fps in m/s
            bc_value: 0.5,
            bullet_mass: 0.0109,      // 168 grains in kg
            bullet_diameter: 0.00782, // 0.308 inches in meters
            bullet_length: 0.033,     // in meters (1.3 inches)
            twist_rate: 10.0,
            ..Default::default()
        }
    }

    #[test]
    fn test_compute_stability_coefficient() {
        let inputs = create_test_inputs();
        let atmo_params = (0.0, 15.0, 1013.25, 1.0); // Standard conditions

        let stability = compute_stability_coefficient(&inputs, atmo_params);

        // Debug output
        println!("Computed stability: {}", stability);

        // Should be a reasonable stability value
        assert!(stability > 0.0);
        assert!(stability < 10.0); // Sanity check

        // Test with typical values should give stability around 1.5-2.5
        assert!(stability > 1.0);
        assert!(stability < 3.0);
    }

    #[test]
    fn test_compute_stability_coefficient_zero_cases() {
        let mut inputs = create_test_inputs();
        let atmo_params = (0.0, 15.0, 1013.25, 1.0);

        // Test with zero twist rate
        inputs.twist_rate = 0.0;
        assert_eq!(compute_stability_coefficient(&inputs, atmo_params), 0.0);

        // Test with zero bullet length
        inputs = create_test_inputs();
        inputs.bullet_length = 0.0;
        assert_eq!(compute_stability_coefficient(&inputs, atmo_params), 0.0);

        // Test with zero bullet diameter
        inputs = create_test_inputs();
        inputs.bullet_diameter = 0.0;
        assert_eq!(compute_stability_coefficient(&inputs, atmo_params), 0.0);
    }

    #[test]
    fn test_compute_stability_coefficient_atmospheric_effects() {
        let inputs = create_test_inputs();

        // Standard conditions
        let standard_atmo = (0.0, 15.0, 1013.25, 1.0);
        let standard_stability = compute_stability_coefficient(&inputs, standard_atmo);

        // High altitude (lower pressure, lower temperature)
        let high_alt_atmo = (3000.0, 5.0, 900.0, 1.0);
        let high_alt_stability = compute_stability_coefficient(&inputs, high_alt_atmo);

        // High altitude should have higher stability due to lower air density
        assert!(high_alt_stability > standard_stability);

        // Hot conditions (higher temperature)
        let hot_atmo = (0.0, 35.0, 1013.25, 1.0);
        let hot_stability = compute_stability_coefficient(&inputs, hot_atmo);

        // Hot conditions should have higher stability due to lower air density
        assert!(hot_stability > standard_stability);
    }

    #[test]
    fn test_compute_spin_drift() {
        let time_s = 1.5;
        let stability = 2.0;
        let twist_rate = 10.0;

        // Test right-hand twist
        let drift_right = compute_spin_drift(time_s, stability, twist_rate, true);
        assert!(drift_right > 0.0); // Should drift to the right (positive)

        // Test left-hand twist
        let drift_left = compute_spin_drift(time_s, stability, twist_rate, false);
        assert!(drift_left < 0.0); // Should drift to the left (negative)
        assert!((drift_left + drift_right).abs() < 1e-10); // Should be equal magnitude

        // Litz spin drift remains a small correction for this flight time.
        assert!(drift_right.abs() < 0.25); // Less than 25cm for 1.5s flight
    }

    #[test]
    fn test_compute_spin_drift_zero_cases() {
        // Test with zero stability
        assert_eq!(compute_spin_drift(1.5, 0.0, 10.0, true), 0.0);

        // Test with zero time
        assert_eq!(compute_spin_drift(0.0, 2.0, 10.0, true), 0.0);

        // Test with negative time
        assert_eq!(compute_spin_drift(-1.0, 2.0, 10.0, true), 0.0);

        // Test with zero twist rate
        assert_eq!(compute_spin_drift(1.5, 2.0, 0.0, true), 0.0);
    }

    // --- MBA-1135: mass-based length + Miller-inverse twist defaults ---

    const GR_TO_KG: f64 = 0.00006479891;
    const IN_TO_M: f64 = 0.0254;

    fn len_in(diameter_in: f64, mass_gr: f64) -> f64 {
        estimate_bullet_length_m(diameter_in * IN_TO_M, mass_gr * GR_TO_KG) / IN_TO_M
    }

    fn twist_in(diameter_in: f64, mass_gr: f64, v_fps: f64) -> f64 {
        default_twist_inches(diameter_in * IN_TO_M, mass_gr * GR_TO_KG, v_fps * 0.3048)
    }

    #[test]
    fn test_estimate_bullet_length_reference_bullets() {
        // (diameter_in, mass_gr, expected_in, tolerance_in)
        let cases = [
            (0.308, 175.0, 1.24, 0.06), // .308 175 gr match  -> L/d ~4.0
            (0.224, 77.0, 0.99, 0.06),  // .224 77 gr         -> long for caliber
            (0.338, 300.0, 1.74, 0.06), // .338 300 gr
            (0.224, 55.0, 0.72, 0.05),  // .224 55 gr varmint -> L/d ~3.2
            (0.510, 750.0, 1.90, 0.08), // .510 750 gr
        ];
        for (d, m, expected, tol) in cases {
            let got = len_in(d, m);
            assert!(
                (got - expected).abs() < tol,
                ".{d}/{m}gr length: expected ~{expected}\", got {got:.4}\" (L/d {:.2})",
                got / d
            );
        }
    }

    #[test]
    fn test_estimate_bullet_length_degenerate_inputs() {
        assert_eq!(estimate_bullet_length_m(0.00782, 0.0), 0.0);
        assert_eq!(estimate_bullet_length_m(0.00782, -1.0), 0.0);
        assert_eq!(estimate_bullet_length_m(0.0, 0.011), 0.0);
        assert_eq!(estimate_bullet_length_m(-0.1, 0.011), 0.0);
        assert_eq!(estimate_bullet_length_m(f64::NAN, 0.011), 0.0);
    }

    #[test]
    fn test_estimate_bullet_length_clamps_ld_ratio() {
        // Absurdly heavy for a tiny caliber -> clamp at 6.5 calibers.
        let d = 0.172 * IN_TO_M;
        let long = estimate_bullet_length_m(d, 0.05); // way too heavy
        assert!((long / d - 6.5).abs() < 1e-9, "expected clamp at 6.5 cal, got {}", long / d);
        // Absurdly light -> clamp at 2.5 calibers.
        let short = estimate_bullet_length_m(0.510 * IN_TO_M, 0.001);
        assert!(
            (short / (0.510 * IN_TO_M) - 2.5).abs() < 1e-9,
            "expected clamp at 2.5 cal, got {}",
            short / (0.510 * IN_TO_M)
        );
    }

    #[test]
    fn test_default_twist_reference_bullets() {
        // With DEFAULT_TWIST_SG_TARGET = 2.0, these should land near common factory twists.
        // (diameter_in, mass_gr, v_fps, expected_twist_in, tolerance_in)
        let cases = [
            (0.308, 175.0, 2600.0, 11.2, 1.2), // .308 175 gr -> ~1:10-12"
            (0.224, 55.0, 3240.0, 10.1, 1.5),  // .223 55 gr  -> ~1:10-12"
            (0.224, 77.0, 2750.0, 7.2, 1.2),   // .224 77 gr  -> ~1:7-8"
            (0.264, 140.0, 2700.0, 7.7, 1.2),  // 6.5mm 140 gr -> ~1:8"
        ];
        for (d, m, v, expected, tol) in cases {
            let got = twist_in(d, m, v);
            assert!(got > 0.0, ".{d}/{m}gr twist must be positive, got {got}");
            assert!(
                (got - expected).abs() < tol,
                ".{d}/{m}gr twist: expected ~1:{expected}\", got 1:{got:.2}\"",
            );
        }
    }

    #[test]
    fn test_default_twist_yields_target_sg() {
        // By construction, at sea-level standard the synthesized twist should reproduce
        // DEFAULT_TWIST_SG_TARGET when the same estimated length is used for the Sg check.
        let d_m = 0.308 * IN_TO_M;
        let m_kg = 175.0 * GR_TO_KG;
        let v_mps = 2600.0 * 0.3048;
        let twist = default_twist_inches(d_m, m_kg, v_mps);
        let inputs = BallisticInputs {
            muzzle_velocity: v_mps,
            bullet_mass: m_kg,
            bullet_diameter: d_m,
            bullet_length: estimate_bullet_length_m(d_m, m_kg),
            twist_rate: twist,
            ..Default::default()
        };
        let sg = compute_stability_coefficient(&inputs, (0.0, 15.0, 1013.25, 1.0));
        assert!(
            (sg - DEFAULT_TWIST_SG_TARGET).abs() < 0.02,
            "expected Sg ~{DEFAULT_TWIST_SG_TARGET}, got {sg}"
        );
    }

    #[test]
    fn test_default_twist_degenerate_inputs_fall_back() {
        assert_eq!(default_twist_inches(0.0, 0.011, 800.0), 12.0);
        assert_eq!(default_twist_inches(0.00782, 0.0, 800.0), 12.0);
        assert_eq!(default_twist_inches(0.00782, 0.011, 0.0), 12.0);
        assert_eq!(default_twist_inches(f64::NAN, 0.011, 800.0), 12.0);
    }

    #[test]
    fn test_resolve_twist_preserves_explicit_or_uses_miller_default() {
        let diameter_m = 0.224 * IN_TO_M;
        let mass_kg = 77.0 * GR_TO_KG;
        let velocity_mps = 2750.0 * 0.3048;
        let expected_default = default_twist_inches(diameter_m, mass_kg, velocity_mps);

        assert_eq!(
            resolve_twist_inches(Some(9.5), diameter_m, mass_kg, velocity_mps).to_bits(),
            9.5_f64.to_bits()
        );
        assert_eq!(
            resolve_twist_inches(None, diameter_m, mass_kg, velocity_mps).to_bits(),
            expected_default.to_bits()
        );
        assert_ne!(expected_default.to_bits(), 12.0_f64.to_bits());

        let metric_default = resolve_twist_inches(
            None,
            5.6896 * 0.001,
            4.98951607 * 0.001,
            838.2,
        );
        assert!((metric_default - expected_default).abs() < 1e-12);
    }

    #[test]
    fn test_heavier_bullet_needs_faster_twist() {
        // Same caliber + velocity: the heavier (longer) bullet must get a faster (smaller) twist.
        let light = twist_in(0.224, 55.0, 2900.0);
        let heavy = twist_in(0.224, 77.0, 2900.0);
        assert!(heavy < light, "77gr twist {heavy} should be faster than 55gr {light}");
    }

    #[test]
    fn test_print_reference_estimates() {
        // Diagnostic dump for the MBA-1135 report: `cargo test print_reference_estimates -- --nocapture`.
        println!("\n=== MBA-1135 estimate_bullet_length_m ===");
        for (d, m) in [
            (0.308, 175.0),
            (0.224, 77.0),
            (0.338, 300.0),
            (0.224, 55.0),
            (0.510, 750.0),
            (0.264, 140.0),
        ] {
            let l = len_in(d, m);
            println!(".{d}/{m}gr -> {l:.4}\" (L/d {:.3})", l / d);
        }
        println!("\n=== MBA-1135 default_twist_inches (SG=1.5 vs 2.0) ===");
        for (d, m, v) in [
            (0.308, 175.0, 2600.0),
            (0.224, 55.0, 3240.0),
            (0.224, 77.0, 2750.0),
            (0.264, 140.0, 2700.0),
        ] {
            // Reconstruct both targets by scaling: t_cal ~ 1/sqrt(Sg) so twist(1.5)=twist(2.0)*sqrt(2.0/1.5).
            let t20 = twist_in(d, m, v);
            let t15 = t20 * (DEFAULT_TWIST_SG_TARGET / 1.5_f64).sqrt();
            println!(".{d}/{m}gr @ {v}fps -> SG1.5: 1:{t15:.2}\"  SG2.0: 1:{t20:.2}\"");
        }
        println!();
    }

    #[test]
    fn test_compute_spin_drift_scaling() {
        let stability = 2.0;
        let twist_rate = 10.0;

        // Test time scaling
        let drift_1s = compute_spin_drift(1.0, stability, twist_rate, true);
        let drift_2s = compute_spin_drift(2.0, stability, twist_rate, true);

        // Drift should increase with time (non-linearly due to 1.83 exponent)
        assert!(drift_2s > drift_1s);

        // Test stability scaling
        let drift_low_stability = compute_spin_drift(1.5, 1.0, twist_rate, true);
        let drift_high_stability = compute_spin_drift(1.5, 3.0, twist_rate, true);

        // Higher stability should produce more drift
        assert!(drift_high_stability > drift_low_stability);
    }
}
