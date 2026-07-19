// MBA-1346: observation-range design for identifiable MV/BC truing.

use ballistics_engine::truing::{
    run_multi_observation_truing_core, DragModelArg, DropUnit, TruingModelInputsV1,
    TruingObservation,
};
use ballistics_engine::truing_plan::{
    plan_truing_experiment_v1, TruingCandidateRejectionReasonV1, TruingExperimentPlanRequestV1,
    TruingPlanErrorCodeV1, TruingPlanModeV1,
};

fn model() -> TruingModelInputsV1 {
    TruingModelInputsV1 {
        muzzle_velocity_fps: 2700.0,
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

#[test]
fn physically_unreachable_candidate_is_excluded_with_a_reason() {
    let plan = plan_truing_experiment_v1(&request(
        vec![300.0, 600.0, 900.0, 1_000.0, 20_000.0],
        3,
        100.0,
        0.03,
    ))
    .expect("reachable candidates still form a plan");

    assert!(plan.rejected_candidates.iter().any(|candidate| {
        candidate.range_yd == 20_000.0
            && candidate.reason == TruingCandidateRejectionReasonV1::Unreachable
            && !candidate.detail.is_empty()
    }));
    assert!(plan
        .selected_stations
        .iter()
        .all(|station| station.range_yd != 20_000.0));
}

fn request(
    candidates: Vec<f64>,
    count: usize,
    separation: f64,
    sigma: f64,
) -> TruingExperimentPlanRequestV1 {
    TruingExperimentPlanRequestV1 {
        model: model(),
        candidate_ranges_yd: candidates,
        observation_count: count,
        minimum_separation_yd: separation,
        measurement_sigma_1sd: sigma,
        drop_unit: DropUnit::Mil,
    }
}

#[test]
fn long_range_design_is_informative_and_recovers_synthetic_parameters() {
    let supplied = vec![
        200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 300.0, -1.0,
    ];
    let plan = plan_truing_experiment_v1(&request(supplied.clone(), 3, 100.0, 0.03)).expect("plan");

    assert_eq!(plan.mode, TruingPlanModeV1::JointMvBc);
    assert_eq!(plan.selected_stations.len(), 3);
    assert_eq!(plan.measurement_sigma_1sd, 0.03);
    assert_eq!(plan.measurement_drop_unit, "mil");
    assert!(plan.information.condition_number.is_some());
    assert!(plan.information.minimum_singular_value > 0.0);
    assert!(plan.information.weak_axis_fractional_sigma_1sd.is_some());
    assert!(plan
        .selected_stations
        .iter()
        .all(|station| supplied.contains(&station.range_yd)));
    for pair in plan.selected_stations.windows(2) {
        assert!(pair[1].range_yd - pair[0].range_yd >= 100.0 - 1.0e-9);
    }
    assert!(plan
        .selected_stations
        .iter()
        .all(|station| station.leave_one_out_information_gain_nats >= 0.0));
    assert!(plan
        .rejected_candidates
        .iter()
        .any(|candidate| { candidate.reason == TruingCandidateRejectionReasonV1::DuplicateRange }));
    assert!(plan
        .rejected_candidates
        .iter()
        .any(|candidate| { candidate.reason == TruingCandidateRejectionReasonV1::InvalidRange }));

    // The planner's nominal predictions are self-consistent synthetic
    // observations.  Starting from a deliberately low BC must recover the
    // profile that was used to design the experiment.
    let observations: Vec<TruingObservation> = plan
        .selected_stations
        .iter()
        .map(|station| TruingObservation {
            range_yd: station.range_yd,
            drop: station.predicted_drop,
        })
        .collect();
    let fitted = run_multi_observation_truing_core(
        &observations,
        DropUnit::Mil,
        0.45,
        DragModelArg::G1,
        168.0,
        0.308,
        100.0,
        2.0,
        59.0,
        29.92,
        50.0,
        0.0,
        &None,
    )
    .expect("synthetic fit");
    assert!(fitted.bc_fitted, "selected design should support joint fit");
    assert!((fitted.fitted_mv_fps - 2700.0).abs() / 2700.0 < 5.0e-4);
    assert!((fitted.fitted_bc - 0.475).abs() / 0.475 < 5.0e-3);
}

#[test]
fn clustered_short_range_facility_is_explicitly_mv_only() {
    let short = plan_truing_experiment_v1(&request(vec![200.0, 250.0, 300.0], 3, 25.0, 0.03))
        .expect("short-range plan");
    assert_eq!(short.mode, TruingPlanModeV1::MvOnly);
    assert!(short
        .warnings
        .iter()
        .any(|warning| warning.contains("MV-only")));

    let long = plan_truing_experiment_v1(&request(vec![300.0, 600.0, 900.0], 3, 100.0, 0.03))
        .expect("long-range plan");
    assert_eq!(long.mode, TruingPlanModeV1::JointMvBc);
    let short_condition = short.information.condition_number.unwrap_or(f64::INFINITY);
    let long_condition = long.information.condition_number.unwrap_or(f64::INFINITY);
    assert!(long_condition < short_condition);
    assert!(
        long.information.minimum_singular_value > short.information.minimum_singular_value,
        "long range should improve the weak information axis"
    );
}

#[test]
fn declared_resolution_scales_weak_axis_uncertainty() {
    let candidates = vec![300.0, 600.0, 900.0];
    let fine = plan_truing_experiment_v1(&request(candidates.clone(), 3, 100.0, 0.03))
        .expect("fine-resolution plan");
    let coarse = plan_truing_experiment_v1(&request(candidates, 3, 100.0, 0.06))
        .expect("coarse-resolution plan");

    let fine_sigma = fine
        .information
        .weak_axis_fractional_sigma_1sd
        .expect("full-rank fine plan");
    let coarse_sigma = coarse
        .information
        .weak_axis_fractional_sigma_1sd
        .expect("full-rank coarse plan");
    assert!((coarse_sigma / fine_sigma - 2.0).abs() < 1.0e-9);
    assert!(
        (coarse.information.minimum_singular_value / fine.information.minimum_singular_value - 0.5)
            .abs()
            < 1.0e-9
    );
    assert!(fine
        .warnings
        .iter()
        .any(|warning| warning.contains("0.5x sigma") && warning.contains("2x sigma")));
}

#[test]
fn planner_detects_when_resolution_changes_the_selected_stations() {
    let candidates = (2..=20).map(|step| step as f64 * 50.0).collect();
    let plan = plan_truing_experiment_v1(&request(candidates, 2, 0.0, 1.0))
        .expect("resolution-sensitive plan");

    assert!(plan.warnings.iter().any(|warning| {
        warning.contains("recommendation is sensitive to measurement resolution")
            && warning.contains("2.0x the declared sigma")
    }));
}

#[test]
fn impossible_minimum_separation_is_structured_error() {
    let error = plan_truing_experiment_v1(&request(vec![100.0, 150.0, 200.0], 3, 100.0, 0.03))
        .expect_err("three stations cannot fit");
    assert_eq!(error.code, TruingPlanErrorCodeV1::NoFeasibleDesign);
}
