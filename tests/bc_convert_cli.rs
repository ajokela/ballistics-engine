//! Native CLI contract for MBA-1375's G1/G7 BC-family converter.

use std::process::{Command, Output};

use ballistics_engine::bc_conversion::{
    analyze_bc_segments, convert_bc_at_mach, convert_bc_at_velocity, format_bc_conversion_report,
    BcConversionFormat, BcConversionReportV1,
};
use ballistics_engine::{BCSegmentData, DragModel};

const BIN: &str = env!("CARGO_BIN_EXE_ballistics");

fn run(args: &[&str]) -> Output {
    Command::new(BIN)
        .args(args)
        .output()
        .expect("run ballistics")
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

#[test]
fn help_discloses_modes_units_and_phase_one_models() {
    let output = run(&["bc-convert", "--help"]);
    assert!(output.status.success(), "{}", stderr(&output));
    let help = String::from_utf8_lossy(&output.stdout);
    for expected in [
        "--source-model <SOURCE_MODEL>",
        "--target-model <TARGET_MODEL>",
        "--bc <BC>",
        "--mach <MACH>",
        "--velocity <VELOCITY>",
        "--bc-segment <VMIN:VMAX:BC>",
        "--speed-of-sound <SPEED_OF_SOUND>",
        "1116.437 fps / 340.29 m/s",
        "Possible values:",
        "g1",
        "g7",
    ] {
        assert!(
            help.contains(expected),
            "help missing {expected:?}:\n{help}"
        );
    }
    assert!(
        !help.to_ascii_lowercase().contains("pdf"),
        "bc-convert help advertises an unsupported PDF output:\n{help}"
    );
}

#[test]
fn phase_one_rejects_an_unshipped_conversion_family_at_parse_time() {
    let output = run(&[
        "bc-convert",
        "--source-model",
        "g2",
        "--target-model",
        "g7",
        "--bc",
        "0.5",
        "--mach",
        "2",
    ]);
    assert!(!output.status.success());
    let err = stderr(&output);
    assert!(err.contains("invalid value 'g2'"), "{err}");
    assert!(err.contains("possible values: g1, g7"), "{err}");
}

#[test]
fn scalar_mode_requires_exactly_one_position() {
    let missing = run(&[
        "bc-convert",
        "--source-model",
        "g1",
        "--target-model",
        "g7",
        "--bc",
        "0.5",
    ]);
    assert!(!missing.status.success());
    assert!(
        stderr(&missing).contains("scalar conversion requires exactly one of --mach or --velocity"),
        "{}",
        stderr(&missing)
    );

    let both = run(&[
        "bc-convert",
        "--source-model",
        "g1",
        "--target-model",
        "g7",
        "--bc",
        "0.5",
        "--mach",
        "2",
        "--velocity",
        "2200",
    ]);
    assert!(!both.status.success());
    let err = stderr(&both);
    assert!(err.contains("cannot be used with"), "{err}");
    assert!(
        err.contains("--mach") || err.contains("--velocity"),
        "{err}"
    );

    let repeated = run(&[
        "bc-convert",
        "--source-model",
        "g1",
        "--target-model",
        "g7",
        "--bc",
        "0.4",
        "--bc",
        "0.5",
        "--mach",
        "2",
    ]);
    assert!(!repeated.status.success());
    assert!(
        stderr(&repeated).contains("cannot be used multiple times"),
        "{}",
        stderr(&repeated)
    );
}

#[test]
fn scalar_and_banded_arguments_conflict_explicitly() {
    let scalar_plus_band = run(&[
        "bc-convert",
        "--source-model",
        "g1",
        "--target-model",
        "g7",
        "--bc",
        "0.5",
        "--mach",
        "2",
        "--bc-segment",
        "1800:2600:0.5",
    ]);
    assert!(!scalar_plus_band.status.success());
    assert!(stderr(&scalar_plus_band).contains("cannot be used with"));

    let band_plus_position = run(&[
        "bc-convert",
        "--source-model",
        "g1",
        "--target-model",
        "g7",
        "--bc-segment",
        "1800:2600:0.5",
        "--velocity",
        "2200",
    ]);
    assert!(!band_plus_position.status.success());
    assert!(stderr(&band_plus_position).contains("cannot be used with"));
}

#[test]
fn speed_of_sound_is_positive_and_does_not_apply_to_explicit_mach() {
    let zero = run(&[
        "bc-convert",
        "--source-model",
        "g1",
        "--target-model",
        "g7",
        "--bc",
        "0.5",
        "--velocity",
        "2200",
        "--speed-of-sound",
        "0",
    ]);
    assert!(!zero.status.success());
    assert!(stderr(&zero).contains("not in range"), "{}", stderr(&zero));

    let with_mach = run(&[
        "bc-convert",
        "--source-model",
        "g1",
        "--target-model",
        "g7",
        "--bc",
        "0.5",
        "--mach",
        "2",
        "--speed-of-sound",
        "1116.437",
    ]);
    assert!(!with_mach.status.success());
    assert!(stderr(&with_mach).contains("cannot be used with"));
}

#[test]
fn scalar_outputs_are_the_shared_formatter_verbatim() {
    let result = convert_bc_at_mach(0.505, DragModel::G1, DragModel::G7, 2.33).unwrap();
    let report = BcConversionReportV1::Scalar { result };
    for (flag, format) in [
        ("table", BcConversionFormat::Table),
        ("csv", BcConversionFormat::Csv),
        ("JSON", BcConversionFormat::Json),
    ] {
        let output = run(&[
            "bc-convert",
            "--source-model",
            "G1",
            "--target-model",
            "G7",
            "--bc",
            "0.505",
            "--mach",
            "2.33",
            "-o",
            flag,
        ]);
        assert!(output.status.success(), "{flag}: {}", stderr(&output));
        assert_eq!(
            String::from_utf8_lossy(&output.stdout),
            format_bc_conversion_report(&report, format).unwrap(),
            "native scalar {flag} output diverged from the shared formatter"
        );
    }
}

#[test]
fn velocity_mode_uses_standard_sound_speed_and_metric_conversion_exactly_once() {
    const METRIC_TO_FPS: f64 = 3.280_839_895;
    let standard_sos_fps =
        ballistics_engine::constants::SPEED_OF_SOUND_MPS / ballistics_engine::constants::FPS_TO_MPS;

    let imperial = run(&[
        "bc-convert",
        "--source-model",
        "g1",
        "--target-model",
        "g7",
        "--bc",
        "0.505",
        "--velocity",
        "2600",
        "-o",
        "json",
    ]);
    assert!(imperial.status.success(), "{}", stderr(&imperial));
    let imperial_expected = BcConversionReportV1::Scalar {
        result: convert_bc_at_velocity(
            0.505,
            DragModel::G1,
            DragModel::G7,
            2600.0,
            standard_sos_fps,
        )
        .unwrap(),
    };
    assert_eq!(
        String::from_utf8_lossy(&imperial.stdout),
        format_bc_conversion_report(&imperial_expected, BcConversionFormat::Json).unwrap()
    );

    let metric = run(&[
        "--units",
        "metric",
        "bc-convert",
        "--source-model",
        "g7",
        "--target-model",
        "g1",
        "--bc",
        "0.26",
        "--velocity",
        "800",
        "--speed-of-sound",
        "343",
        "-o",
        "json",
    ]);
    assert!(metric.status.success(), "{}", stderr(&metric));
    let metric_expected = BcConversionReportV1::Scalar {
        result: convert_bc_at_velocity(
            0.26,
            DragModel::G7,
            DragModel::G1,
            800.0 * METRIC_TO_FPS,
            343.0 * METRIC_TO_FPS,
        )
        .unwrap(),
    };
    assert_eq!(
        String::from_utf8_lossy(&metric.stdout),
        format_bc_conversion_report(&metric_expected, BcConversionFormat::Json).unwrap(),
        "metric velocity or speed of sound was converted more than once"
    );
}

#[test]
fn banded_outputs_include_target_ladder_and_fixed_g1_g7_recommendation() {
    let source_segments = vec![
        BCSegmentData {
            velocity_min: 2400.0,
            velocity_max: 3200.0,
            bc_value: 0.505,
        },
        BCSegmentData {
            velocity_min: 1800.0,
            velocity_max: 2400.0,
            bc_value: 0.480,
        },
        BCSegmentData {
            velocity_min: 1200.0,
            velocity_max: 1800.0,
            bc_value: 0.430,
        },
    ];
    let standard_sos_fps =
        ballistics_engine::constants::SPEED_OF_SOUND_MPS / ballistics_engine::constants::FPS_TO_MPS;
    let result = analyze_bc_segments(
        &source_segments,
        DragModel::G1,
        DragModel::G7,
        &[DragModel::G1, DragModel::G7],
        standard_sos_fps,
    )
    .unwrap();
    assert_eq!(result.conversion.segments.len(), 3);
    assert_eq!(result.recommendation.fits.len(), 2);
    let report = BcConversionReportV1::Banded { result };

    for (flag, format) in [
        ("table", BcConversionFormat::Table),
        ("csv", BcConversionFormat::Csv),
        ("json", BcConversionFormat::Json),
    ] {
        let output = run(&[
            "bc-convert",
            "--source-model",
            "g1",
            "--target-model",
            "g7",
            "--bc-segment",
            "2400:3200:0.505",
            "--bc-segment",
            "1800:2400:0.480",
            "--bc-segment",
            "1200:1800:0.430",
            "-o",
            flag,
        ]);
        assert!(output.status.success(), "{flag}: {}", stderr(&output));
        assert_eq!(
            String::from_utf8_lossy(&output.stdout),
            format_bc_conversion_report(&report, format).unwrap(),
            "native banded {flag} output diverged from the shared formatter"
        );
    }
}

#[test]
fn malformed_bands_pdf_and_out_of_domain_velocity_fail_loudly() {
    let malformed = run(&[
        "bc-convert",
        "--source-model",
        "g1",
        "--target-model",
        "g7",
        "--bc-segment",
        "2400:1800:0.5",
    ]);
    assert!(!malformed.status.success());
    assert!(
        stderr(&malformed).contains("VMIN must be < VMAX"),
        "{}",
        stderr(&malformed)
    );

    let pdf = run(&[
        "bc-convert",
        "--source-model",
        "g1",
        "--target-model",
        "g7",
        "--bc",
        "0.5",
        "--mach",
        "2",
        "-o",
        "pdf",
    ]);
    assert!(!pdf.status.success());
    assert!(
        stderr(&pdf).contains("invalid value 'pdf'"),
        "{}",
        stderr(&pdf)
    );

    let outside = run(&[
        "bc-convert",
        "--source-model",
        "g1",
        "--target-model",
        "g7",
        "--bc",
        "0.5",
        "--velocity",
        "6000",
    ]);
    assert!(!outside.status.success());
    assert!(
        stderr(&outside).contains("outside the G1 reference table domain"),
        "{}",
        stderr(&outside)
    );
}
