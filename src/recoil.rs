//! SAAMI free-recoil calculator (MBA-1372).
//!
//! Source: SAAMI "Gun Recoil - Technical" (Rev. 7/9/2018), the Sporting Arms and
//! Ammunition Manufacturers' Institute's freely downloadable momentum-balance formula:
//! <https://saami.org/wp-content/uploads/2025/03/Gun-Recoil-Formulae-2018-07-9.pdf>
//!
//! ## The physics
//!
//! Recoil is conservation of momentum: the momentum of the free-recoiling firearm is
//! equal and opposite to the momentum of everything that leaves the muzzle -- the bullet
//! plus the propellant gas. SAAMI's document equates the propellant gas *mass* to the
//! powder charge weight ("because the propellant gases are extremely difficult to
//! weigh"), giving:
//!
//! ```text
//! firearm_mass * V_recoil = bullet_mass * V_muzzle + charge_mass * V_gas
//! ```
//!
//! and free recoil energy is then simply `0.5 * firearm_mass * V_recoil^2`.
//!
//! ## Gas-velocity convention (a deliberate modelling choice, not an implementation detail)
//!
//! The propellant gas leaves the muzzle *faster* than the bullet, but its true velocity
//! is hard to measure directly. Two conventions are in circulation:
//!
//! 1. **SAAMI's own type-keyed multiplier** (page 2 of the source PDF): `V_gas = f *
//!    V_muzzle`, with `f` sourced from 1929-era British ballistic testing and keyed to
//!    firearm class: high-powered rifle `f = 1.75`, pistol/revolver `f = 1.50`,
//!    average-length shotgun `f = 1.50`, long-barrel shotgun `f = 1.25`.
//! 2. **A fixed constant** (commonly ~4700 fps / ~1433 m/s for smokeless powder),
//!    popularized by several reloading references, independent of muzzle velocity.
//!
//! **This module defaults to the SAAMI type-keyed factor** (`GasVelocityModel::Saami`)
//! rather than the fixed constant, because: it is the cited, auditable industry-standard
//! source with a worked example we reproduce in the tests below; it scales with muzzle
//! velocity, so it stays physically sane across drastically different loads (a fixed
//! constant over-states gas momentum for a slow subsonic load and under-states it for a
//! magnum); and it is keyed to firearm class the way the industry actually publishes it.
//! The fixed-velocity model is still exposed (`GasVelocityModel::Fixed`) for callers who
//! want to reproduce that older convention or plug in a chronographed gas velocity --
//! "expose both," per the brief, rather than silently picking one and hiding the other.
//!
//! ## Units
//!
//! [`free_recoil`] takes and returns pure SI (kg, m/s, J, N*s) -- the same internal
//! convention `cli_api`/the rest of this crate uses; CLI/WASM front ends convert to/from
//! display units (imperial: grains for bullet/charge, pounds for firearm weight, fps;
//! metric: grams, kilograms, m/s).
//!
//! Note this deliberately differs from SAAMI's own metric appendix (page 4 of the source
//! PDF), which divides a *weight-in-kilograms* by `g = 9.8` to produce energy in legacy
//! "kilogram-meters" (kgf*m) rather than joules. This crate already treats every other
//! mass input (bullet grams/grains, in `cli_api`) as a true SI mass rather than a
//! gravity-dependent "weight," so `free_recoil` performs genuine `F = m*a` physics in
//! kg/m/s/J throughout. No `g` factor is needed for the recoil *velocity* at all (it
//! cancels in the momentum ratio -- see the SAAMI PDF's own algebra, page 1); `g` only
//! ever entered SAAMI's imperial derivation to convert *weight* (pounds) to *mass*
//! (slugs) for the energy step, which plain kilograms sidesteps entirely. The tests below
//! confirm this reproduces the SAAMI worked example (page 3) to within its own rounding.

/// Firearm type, selecting SAAMI's empirical propellant-gas-velocity multiplier `f`
/// (`V_gas = f * V_muzzle`). See the module docs for the source and rationale.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FirearmType {
    /// High-powered rifle: `f = 1.75` (SAAMI's own headline worked case).
    Rifle,
    /// Pistol or revolver: `f = 1.50`.
    Pistol,
    /// Shotgun, average barrel length: `f = 1.50`.
    ShotgunAverage,
    /// Shotgun, long barrel: `f = 1.25`.
    ShotgunLong,
}

impl FirearmType {
    /// SAAMI's `f` multiplier: `V_gas = f * V_muzzle`.
    pub fn saami_gas_factor(self) -> f64 {
        match self {
            FirearmType::Rifle => 1.75,
            FirearmType::Pistol => 1.50,
            FirearmType::ShotgunAverage => 1.50,
            FirearmType::ShotgunLong => 1.25,
        }
    }
}

/// How the propellant gas velocity (SAAMI's `V_PG`) is resolved from the muzzle velocity.
/// See the module docs for why the SAAMI type-keyed factor is the default.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GasVelocityModel {
    /// `V_gas = firearm_type.saami_gas_factor() * V_muzzle` -- the default/recommended path.
    Saami(FirearmType),
    /// `V_gas = factor * V_muzzle`: an explicit override of the SAAMI per-type factor.
    Factor(f64),
    /// `V_gas` is an absolute constant (m/s), independent of muzzle velocity -- e.g. the
    /// popular non-SAAMI rule of thumb of roughly 4700 fps (~1433 m/s) for smokeless
    /// powder, or a chronographed/measured gas velocity.
    Fixed(f64),
}

impl GasVelocityModel {
    /// Resolve the propellant gas velocity (m/s) for a given muzzle velocity (m/s).
    pub fn resolve_mps(self, muzzle_velocity_mps: f64) -> f64 {
        match self {
            GasVelocityModel::Saami(t) => t.saami_gas_factor() * muzzle_velocity_mps,
            GasVelocityModel::Factor(f) => f * muzzle_velocity_mps,
            GasVelocityModel::Fixed(v) => v,
        }
    }
}

/// SI inputs (kg, m/s) to [`free_recoil`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FreeRecoilInputs {
    pub bullet_mass_kg: f64,
    pub charge_mass_kg: f64,
    pub muzzle_velocity_mps: f64,
    /// Firearm weight/mass in kg, including all attachments (scope, suppressor, etc.),
    /// per SAAMI's own definition of `M`.
    pub firearm_mass_kg: f64,
    /// Propellant gas velocity (m/s) -- see [`GasVelocityModel::resolve_mps`].
    pub gas_velocity_mps: f64,
}

/// Free recoil result, SI units throughout.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FreeRecoilResult {
    pub recoil_velocity_mps: f64,
    pub recoil_energy_j: f64,
    /// Recoil impulse == the firearm's momentum magnitude == the ejecta's momentum
    /// magnitude (Newton's third law), in kg*m/s (== N*s).
    pub impulse_ns: f64,
}

/// SAAMI free-recoil momentum balance (see module docs):
/// `firearm_mass * V_recoil = bullet_mass * V_muzzle + charge_mass * V_gas`,
/// `FRE = 0.5 * firearm_mass * V_recoil^2`, `impulse = firearm_mass * V_recoil`.
pub fn free_recoil(inputs: FreeRecoilInputs) -> Result<FreeRecoilResult, String> {
    let FreeRecoilInputs {
        bullet_mass_kg,
        charge_mass_kg,
        muzzle_velocity_mps,
        firearm_mass_kg,
        gas_velocity_mps,
    } = inputs;

    if !firearm_mass_kg.is_finite() || firearm_mass_kg <= 0.0 {
        return Err("firearm weight must be a positive, finite number".to_string());
    }
    if !bullet_mass_kg.is_finite() || bullet_mass_kg <= 0.0 {
        return Err("bullet weight must be a positive, finite number".to_string());
    }
    if !charge_mass_kg.is_finite() || charge_mass_kg < 0.0 {
        return Err("charge weight must be a non-negative, finite number".to_string());
    }
    if !muzzle_velocity_mps.is_finite() || muzzle_velocity_mps <= 0.0 {
        return Err("muzzle velocity must be a positive, finite number".to_string());
    }
    if !gas_velocity_mps.is_finite() || gas_velocity_mps < 0.0 {
        return Err("gas velocity must be a non-negative, finite number".to_string());
    }

    let recoil_velocity_mps = (bullet_mass_kg * muzzle_velocity_mps
        + charge_mass_kg * gas_velocity_mps)
        / firearm_mass_kg;
    let recoil_energy_j = 0.5 * firearm_mass_kg * recoil_velocity_mps * recoil_velocity_mps;
    let impulse_ns = firearm_mass_kg * recoil_velocity_mps;

    Ok(FreeRecoilResult {
        recoil_velocity_mps,
        recoil_energy_j,
        impulse_ns,
    })
}

/// Avoirdupois pound in kilograms, exact by definition (1 lb = 0.45359237 kg). Used to
/// convert firearm weight (pounds, imperial units) to the SI mass [`free_recoil`] needs --
/// the only pound-denominated input in this crate (bullet/charge weight use grains, via
/// `constants::GRAINS_TO_KG`).
pub const POUNDS_TO_KG: f64 = 0.45359237;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::GRAINS_TO_KG;

    const FPS_TO_MPS: f64 = 0.3048;
    const J_TO_FTLB: f64 = 0.737562;

    /// SAAMI's own worked example (source PDF, page 3): a 7 lb average-length shotgun
    /// firing a 12ga 2-3/4"-3-1/2" No. 4 shot load.
    ///   WF  = 7 lb (firearm weight)
    ///   WE  = 1.25 oz shot (546.9 gr) + 43 gr wads = 589.9 gr (ejecta weight)
    ///   WPC = 33.4 gr (propellant charge weight)
    ///   VE  = 1275 fps (ejecta/shot velocity)
    ///   f   = 1.50 (average-length shotgun)
    /// SAAMI's own (rounded) arithmetic: V = (589.9*1275 + 33.4*1275*1.5) / (7000*7)
    ///   = (752122 + 63877.5) / 49000 = 816001/49000 = 16.653 fps (rounded to 16.65 fps
    ///   before squaring in the source document); FRE = (7/64.34)*16.65^2 = 30.22 ft-lb
    ///   ("or about 30 ft-lb due to the uncertainty of the exact shot charge weight and
    ///   velocity" -- SAAMI's own words). Converting through exact SI (kg/m/s, no
    ///   intermediate rounding) instead of SAAMI's g=32.17-and-round-16.65 shortcut gives
    ///   16.653 fps / 30.17 ft-lb -- the same result within that acknowledged rounding
    ///   noise, which is exactly the cross-check this test performs.
    #[test]
    fn saami_shotgun_worked_example() {
        let wf_kg = 7.0 * POUNDS_TO_KG;
        let we_kg = 589.9 * GRAINS_TO_KG;
        let wpc_kg = 33.4 * GRAINS_TO_KG;
        let ve_mps = 1275.0 * FPS_TO_MPS;
        let gas_model = GasVelocityModel::Saami(FirearmType::ShotgunAverage);
        let vpg_mps = gas_model.resolve_mps(ve_mps);

        let result = free_recoil(FreeRecoilInputs {
            bullet_mass_kg: we_kg,
            charge_mass_kg: wpc_kg,
            muzzle_velocity_mps: ve_mps,
            firearm_mass_kg: wf_kg,
            gas_velocity_mps: vpg_mps,
        })
        .expect("valid inputs");

        let recoil_velocity_fps = result.recoil_velocity_mps / FPS_TO_MPS;
        let recoil_energy_ftlb = result.recoil_energy_j * J_TO_FTLB;

        assert!(
            (recoil_velocity_fps - 16.653).abs() < 0.01,
            "expected ~16.653 fps, got {recoil_velocity_fps}"
        );
        // SAAMI's own document lands on 30.22 ft-lb by rounding V to 16.65 fps before
        // squaring; carrying full precision through gives ~30.17. Both are "about 30
        // ft-lb" per the source's own caveat, so allow that documented rounding slack.
        assert!(
            (recoil_energy_ftlb - 30.17).abs() < 0.1,
            "expected ~30.17-30.22 ft-lb, got {recoil_energy_ftlb}"
        );
    }

    /// SAAMI's own per-type gas-velocity factors (source PDF, page 2).
    #[test]
    fn saami_gas_factors_by_firearm_type() {
        assert_eq!(FirearmType::Rifle.saami_gas_factor(), 1.75);
        assert_eq!(FirearmType::Pistol.saami_gas_factor(), 1.50);
        assert_eq!(FirearmType::ShotgunAverage.saami_gas_factor(), 1.50);
        assert_eq!(FirearmType::ShotgunLong.saami_gas_factor(), 1.25);
    }

    #[test]
    fn gas_velocity_model_variants_resolve() {
        let muzzle = 800.0;
        assert_eq!(
            GasVelocityModel::Saami(FirearmType::Rifle).resolve_mps(muzzle),
            1.75 * muzzle
        );
        assert_eq!(GasVelocityModel::Factor(2.0).resolve_mps(muzzle), 2.0 * muzzle);
        // Fixed is independent of the muzzle velocity passed in.
        assert_eq!(GasVelocityModel::Fixed(1433.0).resolve_mps(muzzle), 1433.0);
        assert_eq!(GasVelocityModel::Fixed(1433.0).resolve_mps(1.0), 1433.0);
    }

    /// Zero charge weight degenerates to a bullet-only momentum balance
    /// (no propellant-gas contribution at all).
    #[test]
    fn zero_charge_weight_is_bullet_only_momentum() {
        let result = free_recoil(FreeRecoilInputs {
            bullet_mass_kg: 0.01,
            charge_mass_kg: 0.0,
            muzzle_velocity_mps: 800.0,
            firearm_mass_kg: 4.0,
            gas_velocity_mps: 1500.0, // must be ignored: charge mass is zero
        })
        .expect("valid inputs");
        let expected_v = 0.01 * 800.0 / 4.0;
        assert!((result.recoil_velocity_mps - expected_v).abs() < 1e-12);
    }

    /// Doubling the firearm mass alone must halve recoil energy scaled by velocity^2:
    /// V halves, so E = 0.5*M*V^2 -> (2M)*(V/2)^2 = 0.5*M*V^2 -- i.e. energy also halves.
    #[test]
    fn heavier_firearm_reduces_recoil_energy() {
        let gas_velocity_mps = GasVelocityModel::Saami(FirearmType::Rifle).resolve_mps(823.0);
        let base = FreeRecoilInputs {
            bullet_mass_kg: 0.011,
            charge_mass_kg: 0.0028,
            muzzle_velocity_mps: 823.0,
            firearm_mass_kg: 3.0,
            gas_velocity_mps,
        };
        let light = free_recoil(base).expect("valid inputs");
        let heavy = free_recoil(FreeRecoilInputs {
            firearm_mass_kg: 6.0,
            ..base
        })
        .expect("valid inputs");
        assert!(heavy.recoil_velocity_mps < light.recoil_velocity_mps);
        assert!(heavy.recoil_energy_j < light.recoil_energy_j);
        assert!(
            (heavy.recoil_energy_j - light.recoil_energy_j / 2.0).abs() < 1e-9,
            "doubling firearm mass should exactly halve recoil energy"
        );
    }

    #[test]
    fn impulse_matches_firearm_momentum() {
        let firearm_mass_kg = 3.86;
        let result = free_recoil(FreeRecoilInputs {
            bullet_mass_kg: 0.0109,
            charge_mass_kg: 0.0028,
            muzzle_velocity_mps: 823.0,
            firearm_mass_kg,
            gas_velocity_mps: GasVelocityModel::Saami(FirearmType::Rifle).resolve_mps(823.0),
        })
        .expect("valid inputs");
        let expected_impulse = firearm_mass_kg * result.recoil_velocity_mps;
        assert!(
            (result.impulse_ns - expected_impulse).abs() < 1e-9,
            "impulse must equal firearm_mass * recoil_velocity"
        );
    }

    #[test]
    fn rejects_non_positive_firearm_weight() {
        let base = FreeRecoilInputs {
            bullet_mass_kg: 0.01,
            charge_mass_kg: 0.002,
            muzzle_velocity_mps: 800.0,
            firearm_mass_kg: 0.0,
            gas_velocity_mps: 1200.0,
        };
        assert!(free_recoil(base).is_err());
        assert!(free_recoil(FreeRecoilInputs {
            firearm_mass_kg: -1.0,
            ..base
        })
        .is_err());
    }

    #[test]
    fn rejects_non_positive_bullet_weight() {
        assert!(free_recoil(FreeRecoilInputs {
            bullet_mass_kg: 0.0,
            charge_mass_kg: 0.002,
            muzzle_velocity_mps: 800.0,
            firearm_mass_kg: 4.0,
            gas_velocity_mps: 1200.0,
        })
        .is_err());
    }

    #[test]
    fn rejects_non_positive_muzzle_velocity() {
        assert!(free_recoil(FreeRecoilInputs {
            bullet_mass_kg: 0.01,
            charge_mass_kg: 0.002,
            muzzle_velocity_mps: 0.0,
            firearm_mass_kg: 4.0,
            gas_velocity_mps: 1200.0,
        })
        .is_err());
    }

    #[test]
    fn rejects_non_finite_inputs() {
        assert!(free_recoil(FreeRecoilInputs {
            bullet_mass_kg: f64::NAN,
            charge_mass_kg: 0.002,
            muzzle_velocity_mps: 800.0,
            firearm_mass_kg: 4.0,
            gas_velocity_mps: 1200.0,
        })
        .is_err());
        assert!(free_recoil(FreeRecoilInputs {
            bullet_mass_kg: 0.01,
            charge_mass_kg: 0.002,
            muzzle_velocity_mps: f64::INFINITY,
            firearm_mass_kg: 4.0,
            gas_velocity_mps: 1200.0,
        })
        .is_err());
    }

    #[test]
    fn accepts_zero_gas_velocity() {
        // A degenerate but not invalid input (e.g. an airgun-style GasVelocityModel::Fixed(0.0)).
        assert!(free_recoil(FreeRecoilInputs {
            bullet_mass_kg: 0.01,
            charge_mass_kg: 0.0,
            muzzle_velocity_mps: 800.0,
            firearm_mass_kg: 4.0,
            gas_velocity_mps: 0.0,
        })
        .is_ok());
    }
}
