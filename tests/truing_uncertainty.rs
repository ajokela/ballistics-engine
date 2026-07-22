use std::sync::OnceLock;

use ballistics_engine::truing::{DragModelArg, DropUnit, TruingModelInputsV1};
use ballistics_engine::truing_plan::{plan_truing_experiment_v1, TruingExperimentPlanRequestV1};
#[cfg(feature = "validation")]
use ballistics_engine::truing_uncertainty::TruingMapConvergenceCriterionV1;
use ballistics_engine::truing_uncertainty::{
    run_uncertainty_truing_v1, NormalPriorV1, TruingApproximationFailureCodeV1,
    TruingApproximationV1, TruingPredictionRequestV1, TruingPriorsV1,
    TruingUncertaintyWarningCodeV1, UncertaintyTruingErrorV1, UncertaintyTruingRequestV1,
    WeightedTruingObservationV1,
};

const RANGES_YD: [f64; 5] = [200.0, 300.0, 600.0, 900.0, 1_000.0];

fn model() -> TruingModelInputsV1 {
    TruingModelInputsV1 {
        muzzle_velocity_fps: 2_700.0,
        ballistic_coefficient: 0.475,
        drag_model: DragModelArg::G1,
        mass_gr: 168.0,
        diameter_in: 0.308,
        zero_distance_yd: 100.0,
        sight_height_in: 2.0,
        temperature_f: 59.0,
        pressure_inhg: 29.92,
        humidity_pct: 50.0,
        altitude_ft: 0.0,
    }
}

/// Use the MBA-1346 planner's nominal station predictions to construct exact,
/// self-consistent observations.  The planner and uncertainty fitter share the
/// production forward model and finite-difference Jacobian; caching avoids
/// repeating these relatively expensive trajectories in every test.
fn nominal_drops() -> &'static [(f64, f64)] {
    static DROPS: OnceLock<Vec<(f64, f64)>> = OnceLock::new();
    DROPS
        .get_or_init(|| {
            let plan = plan_truing_experiment_v1(&TruingExperimentPlanRequestV1 {
                model: model(),
                candidate_ranges_yd: RANGES_YD.to_vec(),
                observation_count: RANGES_YD.len(),
                minimum_separation_yd: 0.0,
                measurement_sigma_1sd: 0.02,
                drop_unit: DropUnit::Mil,
            })
            .expect("nominal station predictions");
            assert_eq!(plan.selected_stations.len(), RANGES_YD.len());
            plan.selected_stations
                .into_iter()
                .map(|station| (station.range_yd, station.predicted_drop))
                .collect()
        })
        .as_slice()
}

fn nominal_drop(range_yd: f64) -> f64 {
    nominal_drops()
        .iter()
        .find(|(range, _)| (*range - range_yd).abs() < 1.0e-9)
        .map(|(_, drop)| *drop)
        .expect("range in cached nominal predictions")
}

fn observation(range_yd: f64, sigma: f64) -> WeightedTruingObservationV1 {
    WeightedTruingObservationV1 {
        range_yd,
        drop: nominal_drop(range_yd),
        sigma,
    }
}

fn request(ranges_yd: &[f64], sigma: f64, priors: TruingPriorsV1) -> UncertaintyTruingRequestV1 {
    UncertaintyTruingRequestV1 {
        model: model(),
        drop_unit: DropUnit::Mil,
        observations: ranges_yd
            .iter()
            .map(|range| observation(*range, sigma))
            .collect(),
        priors,
        predictions: vec![TruingPredictionRequestV1 {
            range_yd: 1_100.0,
            future_observation_sigma: Some(0.03),
        }],
    }
}

fn available(
    report: &ballistics_engine::truing_uncertainty::UncertaintyTruingReportV1,
) -> &ballistics_engine::truing_uncertainty::TruingGaussianApproximationV1 {
    match &report.approximation {
        TruingApproximationV1::Available(approximation) => approximation,
        TruingApproximationV1::Unavailable(failure) => {
            panic!("expected Gaussian approximation, got {failure:?}")
        }
    }
}

#[test]
fn covariance_is_psd_and_absolute_sigma_scales_intervals_and_predictive_bands() {
    let mut narrow_request = request(&[300.0, 600.0, 900.0], 0.02, TruingPriorsV1::default());
    narrow_request.predictions.insert(
        0,
        TruingPredictionRequestV1 {
            range_yd: 300.0,
            future_observation_sigma: None,
        },
    );
    let narrow = run_uncertainty_truing_v1(&narrow_request).expect("narrow-sigma fit");
    let wide = run_uncertainty_truing_v1(&request(
        &[300.0, 600.0, 900.0],
        0.04,
        TruingPriorsV1::default(),
    ))
    .expect("wide-sigma fit");

    assert!(narrow.converged);
    assert!((narrow.map_muzzle_velocity_fps - model().muzzle_velocity_fps).abs() < 1.0e-8);
    assert!((narrow.map_ballistic_coefficient - model().ballistic_coefficient).abs() < 1.0e-10);

    let narrow_gaussian = available(&narrow);
    let wide_gaussian = available(&wide);
    let covariance = narrow_gaussian.covariance;
    assert!(covariance.mv_variance_fps2 > 0.0);
    assert!(covariance.bc_variance > 0.0);
    assert!(
        covariance.mv_variance_fps2 * covariance.bc_variance
            - covariance.mv_bc_covariance_fps.powi(2)
            >= -1.0e-12
    );
    assert!(narrow_gaussian.mv_bc_correlation.abs() <= 1.0);
    assert!(narrow_gaussian.muzzle_velocity_interval_95.lower <= model().muzzle_velocity_fps);
    assert!(narrow_gaussian.muzzle_velocity_interval_95.upper >= model().muzzle_velocity_fps);
    assert!(
        narrow_gaussian.ballistic_coefficient_interval_95.lower <= model().ballistic_coefficient
    );
    assert!(
        narrow_gaussian.ballistic_coefficient_interval_95.upper >= model().ballistic_coefficient
    );

    // With absolute known observation sigmas and no priors, doubling every
    // sigma doubles parameter standard deviations (no residual-RMS rescaling).
    let mv_scale = wide_gaussian.muzzle_velocity_interval_95.standard_deviation
        / narrow_gaussian
            .muzzle_velocity_interval_95
            .standard_deviation;
    let bc_scale = wide_gaussian
        .ballistic_coefficient_interval_95
        .standard_deviation
        / narrow_gaussian
            .ballistic_coefficient_interval_95
            .standard_deviation;
    assert!((mv_scale - 2.0).abs() < 1.0e-6, "MV sigma scale={mv_scale}");
    assert!((bc_scale - 2.0).abs() < 1.0e-6, "BC sigma scale={bc_scale}");

    let in_domain = narrow.predictive_bands[0]
        .latent_interval_95
        .expect("in-domain latent band");
    let prediction = narrow.predictive_bands[1];
    let latent = prediction.latent_interval_95.expect("latent band");
    let future = prediction
        .future_observation_interval_95
        .expect("future-observation band");
    assert!(future.standard_deviation > latent.standard_deviation);
    assert!(
        (future.standard_deviation.powi(2) - latent.standard_deviation.powi(2) - 0.03_f64.powi(2))
            .abs()
            < 1.0e-12
    );
    assert!(
        latent.standard_deviation > in_domain.standard_deviation,
        "extrapolated latent band should widen for this design"
    );
    assert!(narrow.warnings.iter().any(|warning| {
        warning.code == TruingUncertaintyWarningCodeV1::PredictionOutsideObservedDomain
    }));

    let json = serde_json::to_value(&narrow).expect("report serializes to JSON");
    assert_eq!(json["schema_version"], 1);
    assert_eq!(json["drop_unit"], "mil");
    assert_eq!(json["approximation"]["status"], "available");
}

#[test]
fn adding_a_consistent_long_range_observation_contracts_both_marginals() {
    let base = run_uncertainty_truing_v1(&request(
        &[300.0, 600.0, 900.0],
        0.03,
        TruingPriorsV1::default(),
    ))
    .expect("base fit");
    let augmented = run_uncertainty_truing_v1(&request(
        &[300.0, 600.0, 900.0, 1_000.0],
        0.03,
        TruingPriorsV1::default(),
    ))
    .expect("augmented fit");
    let base = available(&base);
    let augmented = available(&augmented);

    assert!(
        augmented.muzzle_velocity_interval_95.standard_deviation
            < base.muzzle_velocity_interval_95.standard_deviation
    );
    assert!(
        augmented
            .ballistic_coefficient_interval_95
            .standard_deviation
            < base.ballistic_coefficient_interval_95.standard_deviation
    );
}

#[test]
fn short_range_replicates_leave_bc_broad_and_explicitly_prior_dominated() {
    let bc_prior = NormalPriorV1 {
        mean: model().ballistic_coefficient,
        sigma: 0.20,
    };
    let report = run_uncertainty_truing_v1(&request(
        &[200.0, 200.0, 200.0],
        0.02,
        TruingPriorsV1 {
            muzzle_velocity_fps: None,
            ballistic_coefficient: Some(bc_prior),
        },
    ))
    .expect("prior-regularized short-range fit");
    let gaussian = available(&report);

    assert!(
        gaussian
            .ballistic_coefficient_interval_95
            .standard_deviation
            >= 0.99 * bc_prior.sigma
    );
    assert!(report
        .warnings
        .iter()
        .any(|warning| { warning.code == TruingUncertaintyWarningCodeV1::BcPriorDominated }));
    assert!(report
        .warnings
        .iter()
        .any(|warning| { warning.code == TruingUncertaintyWarningCodeV1::IllConditionedData }));
}

#[test]
fn rank_deficient_data_preserve_map_but_withhold_gaussian_approximation() {
    let report = run_uncertainty_truing_v1(&request(
        &[200.0, 200.0, 200.0],
        0.02,
        TruingPriorsV1::default(),
    ))
    .expect("rank-deficient fit still returns MAP");

    assert!((report.map_muzzle_velocity_fps - model().muzzle_velocity_fps).abs() < 1.0e-8);
    assert!((report.map_ballistic_coefficient - model().ballistic_coefficient).abs() < 1.0e-10);
    match report.approximation {
        TruingApproximationV1::Unavailable(failure) => assert_eq!(
            failure.code,
            TruingApproximationFailureCodeV1::RankDeficientInformation
        ),
        TruingApproximationV1::Available(_) => panic!("rank-deficient information must not invert"),
    }
    assert!(report.predictive_bands[0].latent_interval_95.is_none());
}

#[test]
fn observation_specific_sigmas_control_conflicting_measurement_weight() {
    let nominal_600 = nominal_drop(600.0);
    let precise_drop = nominal_600 + 0.12;
    let imprecise_drop = nominal_600 - 0.12;
    let mut weighted = request(
        &[300.0, 600.0, 600.0, 900.0],
        0.03,
        TruingPriorsV1::default(),
    );
    weighted.observations[1].drop = precise_drop;
    weighted.observations[1].sigma = 0.01;
    weighted.observations[2].drop = imprecise_drop;
    weighted.observations[2].sigma = 0.20;

    let report = run_uncertainty_truing_v1(&weighted).expect("heteroscedastic fit");
    let predicted_600 = report.observations[1].predicted_drop;
    assert!(
        (predicted_600 - precise_drop).abs() < (predicted_600 - imprecise_drop).abs(),
        "prediction {predicted_600} should favor precise {precise_drop} over imprecise {imprecise_drop}"
    );
}

#[test]
fn half_sigma_equals_four_identical_full_sigma_observations() {
    const FULL_SIGMA: f64 = 0.04;
    let shifted_drop = nominal_drop(600.0) + 0.08;

    let mut half_sigma = request(
        &[300.0, 600.0, 900.0],
        FULL_SIGMA,
        TruingPriorsV1::default(),
    );
    half_sigma.observations[1].drop = shifted_drop;
    half_sigma.observations[1].sigma = FULL_SIGMA / 2.0;

    let mut four_replicates = request(
        &[300.0, 600.0, 600.0, 600.0, 600.0, 900.0],
        FULL_SIGMA,
        TruingPriorsV1::default(),
    );
    for observation in &mut four_replicates.observations[1..5] {
        observation.drop = shifted_drop;
    }

    let half_sigma = run_uncertainty_truing_v1(&half_sigma).expect("half-sigma fit");
    let four_replicates = run_uncertainty_truing_v1(&four_replicates).expect("four-replicate fit");
    assert!(
        (half_sigma.map_muzzle_velocity_fps - four_replicates.map_muzzle_velocity_fps).abs()
            < 1.0e-6
    );
    assert!(
        (half_sigma.map_ballistic_coefficient - four_replicates.map_ballistic_coefficient).abs()
            < 1.0e-10
    );

    let half_sigma = available(&half_sigma);
    let four_replicates = available(&four_replicates);
    let covariance_pairs = [
        (
            half_sigma.covariance.mv_variance_fps2,
            four_replicates.covariance.mv_variance_fps2,
        ),
        (
            half_sigma.covariance.mv_bc_covariance_fps,
            four_replicates.covariance.mv_bc_covariance_fps,
        ),
        (
            half_sigma.covariance.bc_variance,
            four_replicates.covariance.bc_variance,
        ),
    ];
    for (half_sigma_value, replicate_value) in covariance_pairs {
        let scale = half_sigma_value
            .abs()
            .max(replicate_value.abs())
            .max(1.0e-20);
        assert!(
            (half_sigma_value - replicate_value).abs() <= 1.0e-10 * scale,
            "inverse-information/covariance mismatch: half-sigma={half_sigma_value:.16e}, replicates={replicate_value:.16e}"
        );
    }
    let information_terms =
        |approximation: &ballistics_engine::truing_uncertainty::TruingGaussianApproximationV1| {
            let covariance = approximation.covariance;
            let determinant = covariance.mv_variance_fps2 * covariance.bc_variance
                - covariance.mv_bc_covariance_fps.powi(2);
            [
                covariance.bc_variance / determinant,
                -covariance.mv_bc_covariance_fps / determinant,
                covariance.mv_variance_fps2 / determinant,
            ]
        };
    for (half_sigma_value, replicate_value) in information_terms(half_sigma)
        .into_iter()
        .zip(information_terms(four_replicates))
    {
        let scale = half_sigma_value
            .abs()
            .max(replicate_value.abs())
            .max(1.0e-20);
        assert!(
            (half_sigma_value - replicate_value).abs() <= 1.0e-9 * scale,
            "posterior-information mismatch: half-sigma={half_sigma_value:.16e}, replicates={replicate_value:.16e}"
        );
    }
    assert!(
        (half_sigma.scaled_information_condition_number
            - four_replicates.scaled_information_condition_number)
            .abs()
            <= 1.0e-10 * half_sigma.scaled_information_condition_number,
        "posterior information condition differs"
    );
}

#[test]
fn invalid_observation_sigma_is_rejected_before_the_forward_model_runs() {
    let request = UncertaintyTruingRequestV1 {
        model: model(),
        drop_unit: DropUnit::Mil,
        observations: vec![
            WeightedTruingObservationV1 {
                range_yd: 300.0,
                drop: 1.0,
                sigma: 0.0,
            },
            WeightedTruingObservationV1 {
                range_yd: 600.0,
                drop: 4.0,
                sigma: 0.02,
            },
        ],
        priors: TruingPriorsV1::default(),
        predictions: Vec::new(),
    };
    let error = run_uncertainty_truing_v1(&request).expect_err("zero sigma must fail");
    assert!(matches!(error, UncertaintyTruingErrorV1::InvalidInput(_)));
    assert!(error.to_string().contains("sigma"));
}

#[test]
fn clustered_short_range_no_prior_never_claims_narrow_bc_uncertainty() {
    let mut starting_model = model();
    starting_model.muzzle_velocity_fps = 3_000.0;
    let request = UncertaintyTruingRequestV1 {
        model: starting_model,
        drop_unit: DropUnit::Mil,
        observations: [
            (200.0, 0.489_804_514),
            (250.0, 0.848_459_369),
            (300.0, 1.247_736_416),
        ]
        .into_iter()
        .map(|(range_yd, drop)| WeightedTruingObservationV1 {
            range_yd,
            drop,
            sigma: 0.03,
        })
        .collect(),
        priors: TruingPriorsV1::default(),
        predictions: Vec::new(),
    };

    let report = run_uncertainty_truing_v1(&request).expect("short-range fit");
    match &report.approximation {
        TruingApproximationV1::Available(gaussian) => assert!(
            gaussian.ballistic_coefficient_interval_95.upper
                - gaussian.ballistic_coefficient_interval_95.lower
                > 0.10,
            "short-range 95% BC interval should remain broad: {gaussian:#?}"
        ),
        TruingApproximationV1::Unavailable(failure) => {
            assert!(matches!(
                failure.code,
                TruingApproximationFailureCodeV1::OptimizerDidNotConverge
                    | TruingApproximationFailureCodeV1::RankDeficientInformation
            ));
            assert!(report.warnings.iter().any(|warning| {
                warning.code == TruingUncertaintyWarningCodeV1::GaussianApproximationUnavailable
            }));
        }
    }
}

/// MBA-1353: synthetic, well-conditioned Monte-Carlo interval calibration.
///
/// Truth is the shared `model()` fixture (MV 2700 fps, BC 0.475). Each trial
/// perturbs the noise-free `nominal_drop()` predictions at the suite's
/// well-conditioned `RANGES_YD` stations (200/300/600/900/1000 yd -- a mix of
/// near and far ranges that separates MV from BC, unlike the short-range-only
/// clusters exercised elsewhere in this file) with iid Gaussian noise and
/// checks whether the resulting 95% intervals contain the true MV/BC. This
/// closes MBA-1353's last acceptance criterion: synthetic well-conditioned
/// cases contain known MV/BC at the advertised interval rate within
/// simulation tolerance.
///
/// Requests are built with `request()`, i.e. `TruingPriorsV1::default()` --
/// no priors -- so the intervals under test are entirely data-dominated.
///
/// The noise generator is a hand-rolled 64-bit LCG feeding Box-Muller, not the
/// `rand` crate: coverage depends on a bit-for-bit reproducible noise sequence
/// across every platform and Rust toolchain this crate is built with, and
/// `rand`'s algorithm/output sequence is not a stability guarantee across its
/// own crate versions (see the identical pattern used in the validation-gated
/// `synthetic_gaussian_coverage_tracks_the_advertised_95_percent_level` below).
///
/// Runtime: ~150s single-threaded for the 40 fits at this seed; the harness
/// runs it in parallel with the rest of this binary, so suite wall-time grows
/// far less. This is the always-on smoke tier; the validation-gated sibling
/// is the high-N rigorous tier.
#[test]
fn intervals_cover_the_true_parameters_at_the_advertised_rate() {
    fn uniform_open(state: &mut u64) -> f64 {
        *state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (((*state >> 11) as f64) + 0.5) / ((1_u64 << 53) as f64)
    }

    fn normal_pair(state: &mut u64) -> (f64, f64) {
        let radius = (-2.0 * uniform_open(state).ln()).sqrt();
        let angle = 2.0 * std::f64::consts::PI * uniform_open(state);
        (radius * angle.cos(), radius * angle.sin())
    }

    // sd = 0.05 mil: large enough that noise dominates solver/finite-difference
    // error, small enough that the 5-station, near+far design stays
    // well-conditioned (matches `nominal_drops()`'s planner sigma of 0.02 in
    // spirit, but scaled up to make Monte-Carlo noise the dominant effect).
    const SIGMA_MIL: f64 = 0.05;
    const TRIALS: usize = 40;
    let mut seed: u64 = 0x1353_2026_0722_c0de_u64;

    let mut available_count = 0usize;
    let mut mv_hits = 0usize;
    let mut bc_hits = 0usize;
    let mut mv_half_width_sum = 0.0_f64;
    let mut bc_half_width_sum = 0.0_f64;

    for trial in 0..TRIALS {
        let mut noise = Vec::with_capacity(RANGES_YD.len());
        while noise.len() < RANGES_YD.len() {
            let (a, b) = normal_pair(&mut seed);
            noise.push(a * SIGMA_MIL);
            if noise.len() < RANGES_YD.len() {
                noise.push(b * SIGMA_MIL);
            }
        }

        let mut trial_request = request(&RANGES_YD, SIGMA_MIL, TruingPriorsV1::default());
        for (observation, delta) in trial_request.observations.iter_mut().zip(&noise) {
            observation.drop += *delta;
        }

        let report = run_uncertainty_truing_v1(&trial_request)
            .unwrap_or_else(|error| panic!("coverage trial {trial} failed: {error}"));
        let gaussian = match &report.approximation {
            TruingApproximationV1::Available(gaussian) => {
                available_count += 1;
                gaussian
            }
            TruingApproximationV1::Unavailable(failure) => panic!(
                "trial {trial} was expected to be well-conditioned but had no \
                 Gaussian approximation: {failure:?}"
            ),
        };

        let mv_interval = gaussian.muzzle_velocity_interval_95;
        let bc_interval = gaussian.ballistic_coefficient_interval_95;
        if mv_interval.lower <= model().muzzle_velocity_fps
            && mv_interval.upper >= model().muzzle_velocity_fps
        {
            mv_hits += 1;
        }
        if bc_interval.lower <= model().ballistic_coefficient
            && bc_interval.upper >= model().ballistic_coefficient
        {
            bc_hits += 1;
        }
        mv_half_width_sum += (mv_interval.upper - mv_interval.lower) / 2.0;
        bc_half_width_sum += (bc_interval.upper - bc_interval.lower) / 2.0;
    }

    let mv_mean_half_width = mv_half_width_sum / TRIALS as f64;
    let bc_mean_half_width = bc_half_width_sum / TRIALS as f64;
    eprintln!(
        "MBA-1353 coverage: available {available_count}/{TRIALS}, MV hits {mv_hits}/{TRIALS} \
         (mean 95% half-width {mv_mean_half_width:.3} fps), BC hits {bc_hits}/{TRIALS} \
         (mean 95% half-width {bc_mean_half_width:.5})"
    );

    // The 5-station near+far design is well-conditioned by construction: every
    // trial must produce a Gaussian approximation.
    assert_eq!(
        available_count, TRIALS,
        "expected every well-conditioned trial to yield a Gaussian approximation"
    );

    // Binomial tail bound: for n=40 independent trials at the advertised true
    // coverage p=0.95, P(hits <= 32) = 7.115e-4 (exact binomial CDF). That is
    // small enough to treat >= 33 as a stable, non-flaky floor: material
    // undercoverage (e.g. p=0.85, where P(hits<=32) ~ 0.30) would fail
    // reliably, while nominal p=0.95 sampling noise essentially never does.
    assert!(mv_hits >= 33, "MV coverage {mv_hits}/{TRIALS} is too low");
    assert!(bc_hits >= 33, "BC coverage {bc_hits}/{TRIALS} is too low");

    // Guard against a degenerate "perfect coverage via absurdly wide
    // intervals" pass. Bounds are pinned at ~2x the values measured for this
    // exact committed seed (mean MV half-width 34.84 fps, mean BC half-width
    // 0.02097): wide enough that fit-noise across toolchains cannot trip them,
    // tight enough that intervals inflated ~2x would fail.
    assert!(
        mv_mean_half_width < 70.0,
        "MV 95% interval implausibly wide on average: {mv_mean_half_width:.3} fps"
    );
    assert!(
        bc_mean_half_width < 0.045,
        "BC 95% interval implausibly wide on average: {bc_mean_half_width:.5}"
    );
}

/// End-to-end coverage is intentionally in the opt-in validation suite: every
/// replicate runs the real nonlinear trajectory and MAP solver.
#[cfg(feature = "validation")]
#[test]
fn synthetic_gaussian_coverage_tracks_the_advertised_95_percent_level() {
    fn uniform_open(seed: &mut u64) -> f64 {
        *seed = seed
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (((*seed >> 11) as f64) + 0.5) / ((1_u64 << 53) as f64)
    }

    fn normal_pair(seed: &mut u64) -> (f64, f64) {
        let radius = (-2.0 * uniform_open(seed).ln()).sqrt();
        let angle = 2.0 * std::f64::consts::PI * uniform_open(seed);
        (radius * angle.cos(), radius * angle.sin())
    }

    const TRIALS: usize = 32;
    const SIGMA_MIL: f64 = 0.03;
    let ranges = [300.0, 600.0, 900.0, 1_000.0];
    let mut seed = 0x1353_2026_0718_C0DE_u64;
    let mut mv_hits = 0;
    let mut bc_hits = 0;
    let mut gradient_convergences = 0;
    let mut objective_mesh_convergences = 0;

    for trial in 0..TRIALS {
        let mut noise = Vec::with_capacity(ranges.len());
        while noise.len() < ranges.len() {
            let (a, b) = normal_pair(&mut seed);
            noise.push(a * SIGMA_MIL);
            if noise.len() < ranges.len() {
                noise.push(b * SIGMA_MIL);
            }
        }
        let mut starting_model = model();
        starting_model.muzzle_velocity_fps = 3_000.0;
        starting_model.ballistic_coefficient = 0.45;
        let request = UncertaintyTruingRequestV1 {
            model: starting_model,
            drop_unit: DropUnit::Mil,
            observations: ranges
                .iter()
                .enumerate()
                .map(|(index, range_yd)| WeightedTruingObservationV1 {
                    range_yd: *range_yd,
                    drop: nominal_drop(*range_yd) + noise[index],
                    sigma: SIGMA_MIL,
                })
                .collect(),
            priors: TruingPriorsV1::default(),
            predictions: Vec::new(),
        };
        let report = run_uncertainty_truing_v1(&request)
            .unwrap_or_else(|error| panic!("coverage trial {trial} failed: {error}"));
        let approximation = match &report.approximation {
            TruingApproximationV1::Available(value) => *value,
            TruingApproximationV1::Unavailable(failure) => {
                panic!(
                    "coverage trial {trial} had no approximation: {failure:?}; report={report:#?}"
                )
            }
        };
        match report.diagnostics.map_convergence_criterion {
            Some(TruingMapConvergenceCriterionV1::ScaledGradient) => gradient_convergences += 1,
            Some(TruingMapConvergenceCriterionV1::ObjectiveMesh) => {
                objective_mesh_convergences += 1
            }
            None => panic!("coverage trial {trial} did not converge: {report:#?}"),
        }
        if approximation.muzzle_velocity_interval_95.lower <= model().muzzle_velocity_fps
            && approximation.muzzle_velocity_interval_95.upper >= model().muzzle_velocity_fps
        {
            mv_hits += 1;
        }
        if approximation.ballistic_coefficient_interval_95.lower <= model().ballistic_coefficient
            && approximation.ballistic_coefficient_interval_95.upper
                >= model().ballistic_coefficient
        {
            bc_hits += 1;
        }
    }

    // For n=32 and true p=.95, P(X <= 26) is about 0.0046. Requiring 27 hits
    // catches material undercoverage while retaining a documented allowance
    // for deterministic binomial sampling variation.
    eprintln!(
        "deterministic coverage from 3000 fps / 0.45 BC: MV {mv_hits}/{TRIALS}, BC {bc_hits}/{TRIALS}; convergence: gradient {gradient_convergences}, objective mesh {objective_mesh_convergences}"
    );
    assert!(mv_hits >= 27, "MV coverage {mv_hits}/{TRIALS} is too low");
    assert!(bc_hits >= 27, "BC coverage {bc_hits}/{TRIALS} is too low");
}
