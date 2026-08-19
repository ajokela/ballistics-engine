// Unit tests for WASM bindings
#[cfg(test)]
#[cfg(target_arch = "wasm32")]
mod tests {
    use crate::wasm::*;
    use wasm_bindgen_test::*;

    wasm_bindgen_test_configure!(run_in_browser);

    #[wasm_bindgen_test]
    fn test_wasm_ballistics_creation() {
        let wasm = WasmBallistics::new();
        assert!(wasm.run_command("help").is_ok());
    }

    #[wasm_bindgen_test]
    fn test_help_command() {
        let wasm = WasmBallistics::new();
        let result = wasm.run_command("help").unwrap();
        assert!(result.contains("Ballistics Engine"));
        assert!(result.contains("trajectory"));
        assert!(result.contains("zero"));
        assert!(result.contains("monte-carlo"));
    }

    #[wasm_bindgen_test]
    fn bc_convert_help_discloses_both_modes_and_phase_one_models() {
        let wasm = WasmBallistics::new();
        let help = wasm.run_command("help").unwrap();
        let index = help
            .split("Global Options:")
            .next()
            .expect("help has a command index");
        assert!(
            index.contains("\n  bc-convert "),
            "bc-convert is missing from the command index: {index}"
        );
        for expected in [
            "BC Convert Command:",
            "--source-model <MODEL>",
            "--target-model <MODEL>",
            "--bc <BC>",
            "--mach <MACH>",
            "--velocity <VEL>",
            "--bc-segment <VMIN:VMAX:BC>",
            "--speed-of-sound <VEL>",
            "g1|g7",
            "table/json/csv",
        ] {
            assert!(help.contains(expected), "help missing {expected:?}: {help}");
        }
        let command_help = wasm.run_command("bc-convert --help").unwrap();
        assert!(command_help.contains("BC Convert Command:"));
        assert!(command_help.contains("--source-model <MODEL>"));
    }

    #[wasm_bindgen_test]
    fn bc_convert_scalar_json_is_the_shared_formatter_verbatim() {
        let actual = WasmBallistics::new()
            .run_command(
                "bc-convert --source-model g1 --target-model g7 --bc 0.475 \
                 --velocity 2700 -o JSON",
            )
            .unwrap();
        let result = crate::bc_conversion::convert_bc_at_velocity(
            0.475,
            crate::DragModel::G1,
            crate::DragModel::G7,
            2700.0,
            crate::constants::SPEED_OF_SOUND_MPS / crate::constants::FPS_TO_MPS,
        )
        .unwrap();
        let expected = crate::bc_conversion::format_bc_conversion_report(
            &crate::bc_conversion::BcConversionReportV1::Scalar { result },
            crate::bc_conversion::BcConversionFormat::Json,
        )
        .unwrap();
        assert_eq!(actual, expected);
        let _: serde_json::Value =
            serde_json::from_str(&actual).expect("bc-convert JSON must be a pure document");
    }

    #[wasm_bindgen_test]
    fn bc_convert_metric_velocity_and_speed_of_sound_reach_the_canonical_fps_core() {
        let actual = WasmBallistics::new()
            .run_command(
                "--units metric bc-convert --source-model g7 --target-model g1 --bc 0.260 \
                 --velocity 800 --speed-of-sound 343 -o json",
            )
            .unwrap();
        let result = crate::bc_conversion::convert_bc_at_velocity(
            0.260,
            crate::DragModel::G7,
            crate::DragModel::G1,
            800.0 * 3.280_839_895,
            343.0 * 3.280_839_895,
        )
        .unwrap();
        let expected = crate::bc_conversion::format_bc_conversion_report(
            &crate::bc_conversion::BcConversionReportV1::Scalar { result },
            crate::bc_conversion::BcConversionFormat::Json,
        )
        .unwrap();
        assert_eq!(actual, expected);
    }

    #[wasm_bindgen_test]
    fn bc_convert_banded_json_includes_the_shared_g1_g7_recommendation() {
        let actual = WasmBallistics::new()
            .run_command(
                "bc-convert --source-model g1 --target-model g7 \
                 --bc-segment 2400:3200:0.505 --bc-segment 1800:2400:0.480 \
                 --bc-segment 1200:1800:0.430 -o json",
            )
            .unwrap();
        let segments = vec![
            crate::BCSegmentData {
                velocity_min: 2400.0,
                velocity_max: 3200.0,
                bc_value: 0.505,
            },
            crate::BCSegmentData {
                velocity_min: 1800.0,
                velocity_max: 2400.0,
                bc_value: 0.480,
            },
            crate::BCSegmentData {
                velocity_min: 1200.0,
                velocity_max: 1800.0,
                bc_value: 0.430,
            },
        ];
        let result = crate::bc_conversion::analyze_bc_segments(
            &segments,
            crate::DragModel::G1,
            crate::DragModel::G7,
            &[crate::DragModel::G1, crate::DragModel::G7],
            crate::constants::SPEED_OF_SOUND_MPS / crate::constants::FPS_TO_MPS,
        )
        .unwrap();
        let expected = crate::bc_conversion::format_bc_conversion_report(
            &crate::bc_conversion::BcConversionReportV1::Banded { result },
            crate::bc_conversion::BcConversionFormat::Json,
        )
        .unwrap();
        assert_eq!(actual, expected);
        assert!(
            actual.contains("\"recommendation\""),
            "banded JSON omitted the family recommendation: {actual}"
        );
        let _: serde_json::Value =
            serde_json::from_str(&actual).expect("banded bc-convert JSON must be pure JSON");
    }

    #[wasm_bindgen_test]
    fn bc_convert_rejects_invalid_modes_models_and_values() {
        let wasm = WasmBallistics::new();
        for (command, expected) in [
            (
                "bc-convert --source-model g2 --target-model g7 --bc 0.5 --mach 2",
                "expected g1 or g7",
            ),
            (
                "bc-convert --source-model g1 --target-model g7 --bc 0.5",
                "exactly one of --mach or --velocity",
            ),
            (
                "bc-convert --source-model g1 --target-model g7 --bc 0.5 --mach 2 \
                 --velocity 2200",
                "exactly one of --mach or --velocity",
            ),
            (
                "bc-convert --source-model g1 --target-model g7 --bc 0.5 --mach 2 \
                 --bc-segment 1800:2600:0.5",
                "--bc cannot be used with --bc-segment",
            ),
            (
                "bc-convert --source-model g1 --target-model g7 \
                 --bc-segment 1800:2600:0.5 --velocity 2200",
                "--bc-segment cannot be used with --mach or --velocity",
            ),
            (
                "bc-convert --source-model g1 --target-model g7 --bc 0.5 --mach 2 \
                 --speed-of-sound 1116.437",
                "--speed-of-sound cannot be used with --mach",
            ),
            (
                "bc-convert --source-model g1 --target-model g7 --bc 0.4 --bc 0.5 --mach 2",
                "--bc cannot be used multiple times",
            ),
            (
                "bc-convert --source-model g1 --target-model g7 \
                 --bc-segment 2600:1800:0.5",
                "VMIN must be finite and less than VMAX",
            ),
            (
                "bc-convert --source-model g1 --target-model g7 --bc 0.5 --mach 2 -o pdf",
                "has no PDF form",
            ),
        ] {
            let error = wasm.run_command(command).unwrap_err();
            let message = error.as_string().unwrap_or_default();
            assert!(
                message.contains(expected),
                "`{command}`: expected {expected:?}, got {message:?}"
            );
        }
    }

    #[wasm_bindgen_test]
    fn test_basic_trajectory_command() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308")
            .unwrap();
        assert!(result.contains("Trajectory Calculation Results"));
        assert!(result.contains("Range"));
        assert!(result.contains("Drop"));
        assert!(result.contains("Velocity"));
    }

    #[wasm_bindgen_test]
    fn test_trajectory_with_auto_zero() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --auto-zero 200")
            .unwrap();
        assert!(result.contains("Rifle zeroed at 200 yards"));
        assert!(result.contains("MOA"));
        assert!(result.contains("mrad"));
    }

    #[wasm_bindgen_test]
    fn test_metric_units() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command("--units metric trajectory -v 823 -b 0.475 -m 10.9 -d 7.82")
            .unwrap();
        assert!(result.contains("m/s"));
        assert!(result.contains("meters"));
        assert!(!result.contains("yards"));
        assert!(!result.contains("fps"));
    }

    #[wasm_bindgen_test]
    fn test_json_output_format() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 -o json")
            .unwrap();
        assert!(result.contains("\"trajectory\""));
        assert!(result.contains("\"summary\""));
        assert!(result.contains("range_yards"));
        assert!(result.contains("drop_inches"));
        // MBA-1315: self-describing units/axes metadata, additive third top-level key
        // alongside the pre-existing "trajectory"/"summary" keys.
        assert!(result.contains("\"legend\""));
        assert!(result.contains("\"axes\""));
        assert!(result.contains("\"drift\""));
    }

    /// MBA-1315 axis-doc-vs-behavior (WASM parity with the native CLI legacy JSON test of the
    /// same name): `legend.axes.drift` must describe the SAME sign convention the formatter
    /// actually produces, not an assumed one. Wind FROM the left (`--wind-direction 270`)
    /// must drift `drift_inches` positive, and FROM the right (`90`) must drift it negative,
    /// exactly as the legend states -- identical in sign and source (`position.z`) to the
    /// native CLI legacy JSON's `x`.
    #[wasm_bindgen_test]
    fn json_legend_drift_axis_matches_observed_crosswind_sign() {
        let wasm = WasmBallistics::new();
        let out = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 -o json")
            .unwrap();
        let doc: serde_json::Value = serde_json::from_str(&out).unwrap();
        let drift_text = doc["legend"]["axes"]["drift"]
            .as_str()
            .expect("legend.axes.drift string");
        assert!(
            drift_text.contains("shooter's") && drift_text.contains("right"),
            "legend.axes.drift wording changed; update this test's sign assumptions: {drift_text}"
        );
        assert!(drift_text.contains("positive"));

        fn last_drift_inches(json: &str) -> f64 {
            let doc: serde_json::Value = serde_json::from_str(json).unwrap();
            doc["trajectory"]
                .as_array()
                .expect("trajectory array")
                .last()
                .expect("at least one trajectory point")["drift_inches"]
                .as_f64()
                .expect("drift_inches is a number")
        }

        let baseline = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --full -o json")
            .unwrap();
        let from_left = wasm
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --wind-speed 10 \
                 --wind-direction 270 --full -o json",
            )
            .unwrap();
        let from_right = wasm
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --wind-speed 10 \
                 --wind-direction 90 --full -o json",
            )
            .unwrap();

        assert_eq!(
            last_drift_inches(&baseline),
            0.0,
            "a no-wind run must have zero lateral drift"
        );
        assert!(
            last_drift_inches(&from_left) > 0.0,
            "wind FROM the left (--wind-direction 270) must drift right (positive) per legend.axes.drift"
        );
        assert!(
            last_drift_inches(&from_right) < 0.0,
            "wind FROM the right (--wind-direction 90) must drift left (negative) per legend.axes.drift"
        );
    }

    #[wasm_bindgen_test]
    fn test_csv_output_format() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 -o csv")
            .unwrap();
        assert!(result.contains("Range(yards),Drop(inches)"));
        assert!(result.contains(","));
        assert!(!result.contains("Trajectory Calculation Results")); // CSV shouldn't have the header
    }

    #[wasm_bindgen_test]
    fn test_zero_command() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command("zero -v 2700 -b 0.475 -m 168 -d 0.308 --target-distance 300")
            .unwrap();
        assert!(result.contains("Zero Calculation Results"));
        assert!(result.contains("Target Distance"));
        assert!(result.contains("MOA Adjustment"));
        assert!(result.contains("Mrad Adjustment"));
    }

    /// MBA-1427: the per-point effective drag coefficient, in the build its requester actually
    /// runs. MBA-1423 shipped it native-only; `with_drag_coefficient` had 6 references in
    /// main.rs and 0 in wasm.rs.
    #[wasm_bindgen_test]
    fn test_with_drag_coefficient_adds_the_key_to_json_points() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command(
                "trajectory -v 2700 -b 0.243 -m 175 -d 0.308 --drag-model g7 --max-range 300 \
                 -o json --with-drag-coefficient",
            )
            .unwrap();
        assert!(
            result.contains("\"drag_coefficient\""),
            "the flag parsed but the key never reached the JSON: {result}"
        );
    }

    /// Off by default: without the flag the JSON must not gain the key, so an existing consumer
    /// sees an unchanged document.
    #[wasm_bindgen_test]
    fn test_drag_coefficient_is_absent_without_the_flag() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command(
                "trajectory -v 2700 -b 0.243 -m 175 -d 0.308 --drag-model g7 --max-range 300 \
                 -o json",
            )
            .unwrap();
        assert!(
            !result.contains("drag_coefficient"),
            "drag_coefficient leaked into unflagged JSON: {result}"
        );
    }

    /// MBA-1426 item 2: the reference drag curves, from the terminal. The output IS the shared
    /// formatter's string — asserted by equality, not by spot-checking fragments, so this
    /// surface cannot drift from native the way the recoil CSV header did.
    #[wasm_bindgen_test]
    fn test_drag_curve_output_is_the_shared_formatter_verbatim() {
        let wasm = WasmBallistics::new();
        for (args, format) in [
            ("drag-curve", crate::drag::ReferenceDragCurveFormat::Table),
            ("drag-curve -o csv", crate::drag::ReferenceDragCurveFormat::Csv),
            ("drag-curve -o json", crate::drag::ReferenceDragCurveFormat::Json),
            (
                "drag-curve --drag-model gs -o csv",
                crate::drag::ReferenceDragCurveFormat::Csv,
            ),
        ] {
            let model = if args.contains("gs") {
                crate::DragModel::GS
            } else {
                crate::DragModel::G7
            };
            assert_eq!(
                wasm.run_command(args).unwrap(),
                crate::drag::format_reference_drag_curve(&model, format),
                "`{args}` diverged from the shared formatter"
            );
        }
    }

    /// A typo'd model must be a hard error, not a silent fall-through to the default — the
    /// silently-ignored-argument class the last two releases spent effort closing.
    #[wasm_bindgen_test]
    fn test_drag_curve_rejects_an_unknown_model_and_argument() {
        let wasm = WasmBallistics::new();
        assert!(wasm.run_command("drag-curve --drag-model g9").is_err());
        assert!(wasm.run_command("drag-curve --bogus-flag").is_err());
        // A junk format must get the invalid-format message, not be told PDF has no form —
        // blaming a format the user never mentioned (review finding).
        assert!(wasm.run_command("drag-curve -o yaml").is_err());
        assert!(wasm.run_command("drag-curve -o pdf").is_err());
    }

    /// Both additions must be discoverable, or they repeat --from-angle's fate (MBA-1418).
    ///
    /// The drag-curve assertion is anchored to the Commands INDEX row, not a bare substring:
    /// the review found the first draft's `contains("drag-curve")` was vacuous, satisfied by a
    /// pre-existing ".drg vendor drag-curve text" phrase in the loadDragTable Host API prose —
    /// it passed on a build with no drag-curve help at all, and the index row was in fact
    /// missing.
    #[wasm_bindgen_test]
    fn test_help_lists_the_new_surfaces() {
        let wasm = WasmBallistics::new();
        let help = wasm.run_command("help").unwrap();
        assert!(help.contains("--with-drag-coefficient"));
        let index = help
            .split("Global Options:")
            .next()
            .expect("help has a Commands index before Global Options");
        assert!(
            index.contains("\n  drag-curve "),
            "drag-curve is missing from the Commands index a user actually scans: {index}"
        );
        assert!(help.contains("Drag Curve Command:"));
    }

    /// MBA-1402 parity, reported by an external consumer against the 0.29.0 WASM build: native
    /// emits the solved zero angle on its JSON surface, the browser build emitted it only in the
    /// text banner. Someone parsing the WASM JSON therefore could not read back the angle the
    /// same run had just printed.
    #[wasm_bindgen_test]
    fn test_auto_zero_angle_is_present_in_json() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command(
                "trajectory -v 2650 -b 0.19 -m 77 -d 0.224 --drag-model g7 --sight-height 2.48 \
                 --auto-zero 100 --max-range 300 --ignore-ground-impact -o json",
            )
            .unwrap();
        assert!(
            result.contains("zero_angle_degrees"),
            "auto-zero JSON is missing zero_angle_degrees: {result}"
        );
    }

    /// And absent — not null — when auto-zero did not run, so a consumer that never uses it sees
    /// an unchanged document.
    #[wasm_bindgen_test]
    fn test_zero_angle_is_absent_without_auto_zero() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command(
                "trajectory -v 2650 -b 0.19 -m 77 -d 0.224 --drag-model g7 --max-range 300 -o json",
            )
            .unwrap();
        assert!(
            !result.contains("zero_angle_degrees"),
            "zero_angle_degrees must be absent on a run with no auto-zero: {result}"
        );
    }

    /// MBA-1426 items 3+4: the terminal's zero previously REJECTED every atmosphere flag
    /// (both its solve branches ran against AtmosphericConditions::default()), while native
    /// zero has carried them since MBA-1397/MBA-1422. These pin that the flags are accepted
    /// AND reach the physics.
    #[wasm_bindgen_test]
    fn test_zero_atmosphere_flags_reach_the_solve() {
        let wasm = WasmBallistics::new();
        let base = "zero -v 2700 -b 0.243 -m 175 -d 0.308 --target-distance 300";
        let sea = wasm.run_command(base).unwrap();
        let high = wasm
            .run_command(&format!("{base} --altitude 5000"))
            .unwrap();
        assert_ne!(sea, high, "--altitude parsed but never reached the zero solve");

        let qnh = wasm
            .run_command(&format!("{base} --pressure 24.90 --altitude 5000 --pressure-type qnh"))
            .unwrap();
        let abs = wasm
            .run_command(&format!(
                "{base} --pressure 24.90 --altitude 5000 --pressure-type absolute"
            ))
            .unwrap();
        assert_ne!(qnh, abs, "--pressure-type parsed but inert on zero");
    }

    /// Density altitude supersedes pressure on the terminal's zero, with the notice spliced
    /// into the output text (this surface has no stderr).
    #[wasm_bindgen_test]
    fn test_zero_density_altitude_supersedes_and_says_so() {
        let wasm = WasmBallistics::new();
        let out = wasm
            .run_command(
                "zero -v 2700 -b 0.243 -m 175 -d 0.308 --target-distance 300 \
                 --density-altitude 5000 --pressure 29.92",
            )
            .unwrap();
        assert!(
            out.contains("--density-altitude supersedes"),
            "supersede notice missing: {out}"
        );
    }

    /// MBA-1426 item 5: the from-angle output names which root the single zero-range value
    /// refers to, matching native's primary_crossing disclosure.
    #[wasm_bindgen_test]
    fn test_from_angle_names_the_primary_crossing() {
        let wasm = WasmBallistics::new();
        let out = wasm
            .run_command("zero -v 2700 -b 0.243 -m 175 -d 0.308 --from-angle 0.14")
            .unwrap();
        assert!(
            out.contains("Primary Crossing:"),
            "primary crossing line missing: {out}"
        );
    }

    /// MBA-1426 item 4: the zero-day density altitude on trajectory --auto-zero.
    #[wasm_bindgen_test]
    fn test_trajectory_zero_day_density_altitude_changes_the_zero() {
        let wasm = WasmBallistics::new();
        let base = "trajectory -v 2700 -b 0.243 -m 175 -d 0.308 --drag-model g7 \
                    --max-range 500 --auto-zero 300";
        let plain = wasm.run_command(base).unwrap();
        let zday = wasm
            .run_command(&format!("{base} --zero-density-altitude 5000"))
            .unwrap();
        assert_ne!(plain, zday, "--zero-density-altitude parsed but inert");
    }

    /// MBA-1426 item 3: --pressure-type on the calculators that take --pressure.
    #[wasm_bindgen_test]
    fn test_calculators_honor_pressure_type() {
        let wasm = WasmBallistics::new();
        for (qnh, abs) in [
            (
                "true-velocity --measured-drop 30 --range 500 -b 0.243 -m 175 -d 0.308 \
                 --pressure 24.90 --altitude 5000 --pressure-type qnh --offline",
                "true-velocity --measured-drop 30 --range 500 -b 0.243 -m 175 -d 0.308 \
                 --pressure 24.90 --altitude 5000 --pressure-type absolute --offline",
            ),
            (
                "lead --target-speed 10 -v 2700 -b 0.243 -m 175 -d 0.308 \
                 --pressure 24.90 --altitude 5000 --pressure-type qnh",
                "lead --target-speed 10 -v 2700 -b 0.243 -m 175 -d 0.308 \
                 --pressure 24.90 --altitude 5000 --pressure-type absolute",
            ),
        ] {
            let a = wasm.run_command(qnh);
            let b = wasm.run_command(abs);
            match (a, b) {
                (Ok(x), Ok(y)) => assert_ne!(x, y, "pressure-type inert on: {qnh}"),
                (a, b) => panic!("calculator run failed: {a:?} / {b:?}"),
            }
        }
    }

    /// MBA-1418: `zero --from-angle` works in the browser terminal but had no tests here and
    /// was absent from the help text, so it was effectively undiscoverable. The inverse reports
    /// BOTH line-of-sight crossings a bore angle produces, and names which one `zero_range` is.
    #[wasm_bindgen_test]
    fn test_zero_from_angle_reports_both_crossings() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command("zero -v 2700 -b 0.243 -m 175 -d 0.308 --from-angle 0.14")
            .unwrap();
        assert!(result.contains("Zero"), "unexpected output: {result}");
        // A bore angle generally crosses the line of sight twice; both must be reported rather
        // than the tool silently picking one.
        assert!(
            result.contains("Near Zero") || result.contains("near"),
            "the near crossing is missing: {result}"
        );
        assert!(
            result.contains("Far Zero") || result.contains("far"),
            "the far crossing is missing: {result}"
        );
    }

    /// The two modes are mutually exclusive on native, and the terminal must agree rather than
    /// silently honouring one of them.
    #[wasm_bindgen_test]
    fn test_zero_from_angle_conflicts_with_target_distance() {
        let wasm = WasmBallistics::new();
        let result = wasm.run_command(
            "zero -v 2700 -b 0.243 -m 175 -d 0.308 --from-angle 0.14 --target-distance 300",
        );
        assert!(
            result.is_err(),
            "--from-angle with --target-distance must be refused, got: {result:?}"
        );
    }

    /// `--from-angle` must appear in the terminal help, or a browser user has no way to find it.
    #[wasm_bindgen_test]
    fn test_help_lists_from_angle() {
        let wasm = WasmBallistics::new();
        let help = wasm.run_command("help").unwrap();
        assert!(
            help.contains("--from-angle"),
            "the zero command's --from-angle is missing from the terminal help"
        );
    }

    #[wasm_bindgen_test]
    fn test_monte_carlo_command() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command("monte-carlo -v 2700 -b 0.475 -m 168 -d 0.308 -n 100")
            .unwrap();
        assert!(result.contains("Monte Carlo Simulation Results"));
        assert!(result.contains("Simulations Run: 100"));
        assert!(result.contains("Range Statistics"));
        assert!(result.contains("Impact Velocity Statistics"));
    }

    #[wasm_bindgen_test]
    fn test_estimate_bc_command() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command("estimate-bc -v 2700 -m 168 -d 0.308 --data \"100,2.5;200,10.2;300,23.5\"")
            .unwrap();
        assert!(result.contains("BC Estimation Results"));
        assert!(result.contains("Estimated BC"));
        // Default --drag-model is "both", so a G1 and a G7 row are printed, each fit to the
        // 3-point drop series.
        assert!(result.contains("G1"));
        assert!(result.contains("G7"));
        assert!(result.contains("drop (3 pts)"));
    }

    #[wasm_bindgen_test]
    fn test_estimate_bc_velocity_and_drag_model() {
        let wasm = WasmBallistics::new();
        // G7-only, fit against a velocity-retention series.
        let result = wasm
            .run_command(
                "estimate-bc -v 2650 -m 77 -d 0.224 \
                 --velocity-data \"200,2270;400,1930;600,1610\" --drag-model g7",
            )
            .unwrap();
        assert!(result.contains("BC Estimation Results"));
        assert!(result.contains("G7"));
        assert!(result.contains("velocity (3 pts)"));
        // g7-only must NOT print a G1 row.
        assert!(!result.contains("G1"));
    }

    #[wasm_bindgen_test]
    fn test_advanced_physics_flags() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 \
             --enable-magnus --enable-coriolis --enable-spin-drift \
             --twist-rate 10 --latitude 45",
            )
            .unwrap();
        assert!(result.contains("Trajectory Calculation Results"));
        // The trajectory should complete without errors when physics flags are enabled
    }

    #[wasm_bindgen_test]
    fn omitted_twist_matches_miller_default_and_explicit_override_wins() {
        let wasm = WasmBallistics::new();
        let final_drift = |output: &str, key: &str| {
            let json: serde_json::Value = serde_json::from_str(output).unwrap();
            json["trajectory"].as_array().unwrap().last().unwrap()[key]
                .as_f64()
                .unwrap()
        };
        let command = "trajectory -v 2750 -b 0.219 -m 77 -d 0.224 --drag-model G7 \
                       --enable-spin-drift --max-range 500 -o json";
        let default_twist = crate::stability::default_twist_inches(
            0.224 * 0.0254,
            77.0 * crate::constants::GRAINS_TO_KG,
            2750.0 * 0.3048,
        );

        let omitted = wasm.run_command(command).unwrap();
        let explicit_default = wasm
            .run_command(&format!("{command} --twist-rate {default_twist:.17}"))
            .unwrap();
        let explicit_twelve = wasm
            .run_command(&format!("{command} --twist-rate 12"))
            .unwrap();

        let omitted_drift = final_drift(&omitted, "drift_inches");
        let default_drift = final_drift(&explicit_default, "drift_inches");
        let twelve_drift = final_drift(&explicit_twelve, "drift_inches");
        assert!((omitted_drift - default_drift).abs() < 1e-9);
        assert!((omitted_drift - twelve_drift).abs() > 0.1);

        let metric_command = "trajectory --units metric -v 838.2 -b 0.219 -m 4.98951607 \
                              -d 5.6896 --drag-model G7 --enable-spin-drift \
                              --max-range 457.2 -o json";
        let metric_default_twist =
            crate::stability::default_twist_inches(5.6896 * 0.001, 4.98951607 * 0.001, 838.2);
        let omitted_metric = wasm.run_command(metric_command).unwrap();
        let explicit_default_metric = wasm
            .run_command(&format!(
                "{metric_command} --twist-rate {:.17}",
                metric_default_twist * 25.4
            ))
            .unwrap();
        assert!(
            (final_drift(&omitted_metric, "drift_cm")
                - final_drift(&explicit_default_metric, "drift_cm"))
            .abs()
                < 1e-9
        );
    }

    #[wasm_bindgen_test]
    fn test_environmental_conditions() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 \
             --temperature 32 --pressure 25.0 --humidity 80 --altitude 5000 \
             --wind-speed 10 --wind-direction 90",
            )
            .unwrap();
        assert!(result.contains("Trajectory Calculation Results"));
    }

    #[wasm_bindgen_test]
    fn test_full_trajectory_output() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --full --max-range 500")
            .unwrap();
        // With --full flag, we should see more data points
        let lines: Vec<&str> = result.lines().collect();
        let data_lines = lines.iter().filter(|l| l.contains("yd")).count();
        assert!(data_lines > 5); // Should have many data points with --full
    }

    #[wasm_bindgen_test]
    fn test_powder_temperature_sensitivity() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 \
             --use-powder-sensitivity --powder-temp 90 --powder-temp-sensitivity 1.5",
            )
            .unwrap();
        assert!(result.contains("Trajectory Calculation Results"));
    }

    #[wasm_bindgen_test]
    fn explicit_zero_velocity_ignores_linear_powder_adjustment() {
        let wasm = WasmBallistics::new();
        // Table format: the "Rifle zeroed at ..." banner is a table-only human header
        // (MBA-1294(a) keeps JSON/CSV pure), and this test reads the banner's adjustment line.
        let command = "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 \
                       --auto-zero 300 --zero-velocity 2400 --zero-temperature 20 \
                       --max-range 100";
        let zero_line = |output: &str| output.lines().next().unwrap().to_string();

        let without_linear = zero_line(&wasm.run_command(command).unwrap());
        let with_linear = zero_line(
            &wasm
                .run_command(&format!(
                    "{command} --use-powder-sensitivity --powder-temp-sensitivity 5"
                ))
                .unwrap(),
        );

        assert!(without_linear.starts_with("Rifle zeroed at 300 yards"));
        assert_eq!(with_linear, without_linear);
    }

    #[wasm_bindgen_test]
    fn auto_zero_inherits_explicit_shot_day_powder_temperature() {
        let wasm = WasmBallistics::new();
        // Table format: reads the auto-zero banner's adjustment (banner is table-only
        // since MBA-1294(a); JSON/CSV stay pure machine output).
        let command = "trajectory -v 2650 -b 0.19 -m 77 -d 0.224 --drag-model g7 \
                       --temperature 32 --powder-temp 68 \
                       --powder-temp-curve 32:2650,77:2720 \
                       --auto-zero 100 --max-range 100";
        let adjustment = |output: &str| {
            output
                .lines()
                .next()
                .and_then(|line| line.split_once("(Adjustment: "))
                .map(|(_, value)| value.to_string())
                .expect("successful auto-zero adjustment")
        };

        let inherited = adjustment(&wasm.run_command(command).unwrap());
        let explicit = adjustment(
            &wasm
                .run_command(&format!("{command} --zero-powder-temp 68"))
                .unwrap(),
        );

        assert_eq!(inherited, explicit);

        let zero_air = adjustment(
            &wasm
                .run_command(&format!("{command} --zero-temperature 77"))
                .unwrap(),
        );
        let explicit_zero_air = adjustment(
            &wasm
                .run_command(&format!(
                    "{command} --zero-temperature 77 --zero-powder-temp 77"
                ))
                .unwrap(),
        );

        assert_eq!(zero_air, explicit_zero_air);
    }

    #[wasm_bindgen_test]
    fn metric_default_powder_temperature_represents_70_fahrenheit() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command(
                "trajectory --units metric -v 800 -b 0.3 -m 10 -d 7.62 \
                 --temperature 30 --use-powder-sensitivity \
                 --powder-temp-sensitivity 1 --max-range 1 -o json",
            )
            .unwrap();
        let json: serde_json::Value = serde_json::from_str(&result).unwrap();
        let actual_velocity = json["trajectory"][0]["velocity_mps"].as_f64().unwrap();
        let reference_temp_c = (70.0 - 32.0) * 5.0 / 9.0;
        let expected_velocity = 800.0 + (30.0 - reference_temp_c);

        assert!((actual_velocity - expected_velocity).abs() < 1e-9);
    }

    #[wasm_bindgen_test]
    fn test_shooting_angle() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --shooting-angle 15")
            .unwrap();
        assert!(result.contains("Trajectory Calculation Results"));
    }

    #[wasm_bindgen_test]
    fn test_drag_models() {
        let wasm = WasmBallistics::new();

        // Test G1 model
        let g1_result = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --drag-model G1")
            .unwrap();
        assert!(g1_result.contains("Trajectory Calculation Results"));

        // Test G7 model
        let g7_result = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --drag-model G7")
            .unwrap();
        assert!(g7_result.contains("Trajectory Calculation Results"));
    }

    #[wasm_bindgen_test]
    fn test_invalid_command() {
        let wasm = WasmBallistics::new();
        let result = wasm.run_command("invalid-command").unwrap();
        assert!(result.contains("Error: Unknown command"));
        assert!(result.contains("help"));
    }

    #[wasm_bindgen_test]
    fn test_invalid_parameters() {
        let wasm = WasmBallistics::new();
        let result = wasm.run_command("trajectory -v not-a-number -b 0.475 -m 168 -d 0.308");
        assert!(result.is_err());
    }

    #[wasm_bindgen_test]
    fn test_metric_monte_carlo() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command("--units metric monte-carlo -v 823 -b 0.475 -m 10.9 -d 7.82 -n 50")
            .unwrap();
        assert!(result.contains("m/s"));
        assert!(result.contains("meters"));
        assert!(!result.contains("yards"));
        assert!(!result.contains("fps"));
    }

    #[wasm_bindgen_test]
    fn test_metric_zero_calculation() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command(
                "--units metric zero -v 823 -b 0.475 -m 10.9 -d 7.82 --target-distance 200",
            )
            .unwrap();
        assert!(result.contains("200 meters"));
        assert!(!result.contains("yards"));
    }

    #[wasm_bindgen_test]
    fn test_bc_segments() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --use-bc-segments")
            .unwrap();
        assert!(result.contains("Trajectory Calculation Results"));
    }

    #[wasm_bindgen_test]
    fn test_trajectory_sampling() {
        let wasm = WasmBallistics::new();
        let result = wasm.run_command(
            "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --sample-trajectory --sample-interval 25"
        ).unwrap();
        assert!(result.contains("Trajectory Calculation Results"));
    }

    #[wasm_bindgen_test]
    fn test_integration_methods() {
        let wasm = WasmBallistics::new();

        // Test RK4 (default)
        let rk4_result = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308")
            .unwrap();
        assert!(rk4_result.contains("Trajectory Calculation Results"));

        // Test Euler
        let euler_result = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --use-euler")
            .unwrap();
        assert!(euler_result.contains("Trajectory Calculation Results"));
    }

    #[wasm_bindgen_test]
    fn test_estimate_bc_no_data() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command("estimate-bc -v 2700 -m 168 -d 0.308")
            .unwrap();
        assert!(result.contains("Error: No data provided"));
        assert!(result.contains("Example"));
    }

    #[wasm_bindgen_test]
    fn test_zero_with_sight_height() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command(
                "zero -v 2700 -b 0.475 -m 168 -d 0.308 --target-distance 200 --sight-height 2.5",
            )
            .unwrap();
        assert!(result.contains("Sight Height: 2.5 inches"));
    }

    #[wasm_bindgen_test]
    fn test_monte_carlo_with_variations() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command(
                "monte-carlo -v 2700 -b 0.475 -m 168 -d 0.308 -n 50 \
             --velocity-std 15 --angle-std 0.2 --bc-std 0.02 \
             --wind-speed-std 3 --wind-dir-std 10",
            )
            .unwrap();
        assert!(result.contains("Monte Carlo Simulation Results"));
        assert!(result.contains("Std Dev"));
    }

    #[wasm_bindgen_test]
    fn test_command_parsing() {
        let wasm = WasmBallistics::new();

        // Test that ballistics prefix is handled
        let result1 = wasm
            .run_command("ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308")
            .unwrap();
        assert!(result1.contains("Trajectory Calculation Results"));

        // Test without prefix
        let result2 = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308")
            .unwrap();
        assert!(result2.contains("Trajectory Calculation Results"));

        // Test with ./ prefix
        let result3 = wasm
            .run_command("./ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308")
            .unwrap();
        assert!(result3.contains("Trajectory Calculation Results"));
    }

    #[wasm_bindgen_test]
    fn test_empty_command() {
        let wasm = WasmBallistics::new();
        let result = wasm.run_command("").unwrap();
        assert!(result.contains("Ballistics Engine"));
    }

    #[wasm_bindgen_test]
    fn test_help_variants() {
        let wasm = WasmBallistics::new();

        let help1 = wasm.run_command("help").unwrap();
        let help2 = wasm.run_command("--help").unwrap();
        let help3 = wasm.run_command("-h").unwrap();

        assert!(help1.contains("Ballistics Engine"));
        assert!(help2.contains("Ballistics Engine"));
        assert!(help3.contains("Ballistics Engine"));
    }

    #[wasm_bindgen_test]
    fn ignore_ground_impact_reaches_full_range_where_default_truncates() {
        // Reads the number out of the table's trailing "Max Range: <N> yards" line —
        // more robust than the table's own "Bullet struck ground" heuristic, which only
        // fires when the last recorded point lands within 0.01m of the ground and can
        // miss a truncation that stopped a hair above that (as it does for these inputs).
        fn parse_max_range_yards(output: &str) -> f64 {
            output
                .lines()
                .find_map(|l| l.strip_prefix("Max Range: "))
                .and_then(|rest| rest.split_whitespace().next())
                .and_then(|n| n.parse::<f64>().ok())
                .expect("Max Range line present in trajectory output")
        }

        let wasm = WasmBallistics::new();
        // Flat-fire (no auto-zero) drop at 1000 yd far exceeds the default 60in bore
        // height, so the default run truncates at ground impact well short of --max-range.
        let default_result = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 1000")
            .unwrap();
        let default_max_range = parse_max_range_yards(&default_result);
        assert!(
            default_max_range < 900.0,
            "expected the default run to truncate well short of 1000 yd, got {default_max_range}"
        );

        // --ignore-ground-impact (bero-feedback 0.23.0) disables that early termination,
        // so the same run reaches the full requested range instead.
        let ignore_result = wasm
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 1000 \
                 --ignore-ground-impact",
            )
            .unwrap();
        let ignore_max_range = parse_max_range_yards(&ignore_result);
        assert!(
            ignore_max_range >= 999.0,
            "expected --ignore-ground-impact to reach ~1000 yd, got {ignore_max_range}"
        );
    }

    #[wasm_bindgen_test]
    fn muzzle_height_warning_appears_only_above_the_1000m_threshold() {
        let wasm = WasmBallistics::new();
        // ~25 km "muzzle height" (defeating ground truncation the wrong way) must warn.
        let warned = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --muzzle-height 1000000")
            .unwrap();
        assert!(warned.starts_with("WARNING"));
        assert!(warned.contains("--altitude"));
        assert!(warned.contains("--ignore-ground-impact"));

        // An ordinary elevated stand (2400 in ~= 61 m) must not trigger the warning.
        let unwarned = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --muzzle-height 2400")
            .unwrap();
        assert!(!unwarned.contains("WARNING"));
    }

    #[wasm_bindgen_test]
    fn lead_output_json_parses_with_sane_intercept_range() {
        let wasm = WasmBallistics::new();
        let output = wasm
            .run_command("lead --target-speed 10 --target-angle 0 --range 500 -o json")
            .unwrap();
        let json: serde_json::Value = serde_json::from_str(&output).unwrap();
        let range = json["range"].as_f64().unwrap();
        let intercept_range = json["intercept_range"].as_f64().unwrap();
        // angle 0 = directly away (outbound): the target recedes during the time of
        // flight, so the corrected intercept range must be at or beyond the initial range.
        assert!(intercept_range >= range);
        assert_eq!(json["distance_unit"], "yd");
        assert_eq!(json["adjustment_unit"], "mil");
    }

    #[wasm_bindgen_test]
    fn test_all_physics_flags_combined() {
        let wasm = WasmBallistics::new();

        // Test all physics flags together
        let result = wasm
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 \
             --enable-magnus --enable-coriolis --enable-spin-drift \
             --enable-wind-shear --enable-pitch-damping --enable-precession \
             --sample-trajectory --use-bc-segments --use-powder-sensitivity \
             --twist-rate 10 --latitude 45",
            )
            .unwrap();
        assert!(result.contains("Trajectory Calculation Results"));
    }

    /// MBA-1297 (field report, Bero): at 90 degrees of cant the vertical and
    /// lateral misses at the zero distance MUST be equal — the same rotation
    /// moves both the bore offset and the zero tilt. The WASM formatters used
    /// points[0].y + sight_height as the LOS, double-counting the canted
    /// bore's vertical rise (reported drop = true miss + one sight height).
    #[wasm_bindgen_test]
    fn canted_90_vertical_and_lateral_misses_are_symmetric() {
        let wasm = WasmBallistics::new();
        let out = wasm
            .run_command(
                "trajectory -b 0.19 -m 77 -d 0.224 --drag-model g7 --max-range 200 \
                 --sight-height 2.48 --ignore-ground-impact -v 2650 --auto-zero 100 \
                 --cant 90 --full -o csv",
            )
            .unwrap();
        // CSV rows: distance,drop_in,drift_in,... — find the row nearest 100 yd
        // and the first (muzzle) row.
        let mut muzzle_drop: Option<f64> = None;
        let mut drop100: Option<f64> = None;
        let mut drift100: Option<f64> = None;
        for line in out.lines() {
            let cols: Vec<&str> = line.split(',').collect();
            if cols.len() < 3 {
                continue;
            }
            if let (Ok(dist), Ok(drop), Ok(drift)) = (
                cols[0].trim().parse::<f64>(),
                cols[1].trim().parse::<f64>(),
                cols[2].trim().parse::<f64>(),
            ) {
                if dist.abs() < 1.0 && muzzle_drop.is_none() {
                    muzzle_drop = Some(drop);
                }
                if (dist - 100.0).abs() < 6.0 && drop100.is_none() {
                    drop100 = Some(drop);
                    drift100 = Some(drift);
                }
            }
        }
        let muzzle_drop = muzzle_drop.expect("muzzle row");
        let (d, w) = (drop100.expect("100yd drop"), drift100.expect("100yd drift"));
        // At 90 degrees the bore is BESIDE the scope: t=0 vertical offset from the
        // LOS is sh*cos(90) = 0, not the full sight height.
        assert!(
            muzzle_drop.abs() < 0.3,
            "muzzle drop must be ~0 at 90 deg cant (bore beside scope), got {muzzle_drop}"
        );
        // Symmetry: |vertical miss| == |lateral miss| within 5%.
        let ratio = d.abs() / w.abs();
        assert!(
            (0.95..=1.05).contains(&ratio),
            "90-deg cant symmetry broken: drop {d} vs drift {w} (ratio {ratio})"
        );
    }

    // -----------------------------------------------------------------------------
    // MBA-1325: mover ring (`trajectory --target-speed`) WASM parity.
    // -----------------------------------------------------------------------------

    #[wasm_bindgen_test]
    fn target_speed_adds_ring_column_to_trajectory_table() {
        let wasm = WasmBallistics::new();
        let without = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --full --max-range 500")
            .unwrap();
        assert!(!without.contains("Ring"), "Ring column must not appear without --target-speed");

        let with_ring = wasm
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --full --max-range 500 \
                 --target-speed 5",
            )
            .unwrap();
        assert!(with_ring.contains("Ring"), "table should carry a Ring column:\n{with_ring}");
        assert!(with_ring.contains("mil"), "ring values should be labeled mil");
    }

    #[wasm_bindgen_test]
    fn target_speed_adds_mover_ring_m_to_trajectory_json() {
        let wasm = WasmBallistics::new();
        let without = wasm
            .run_command(
                "trajectory --units metric -v 823 -b 0.475 -m 10.9 -d 7.82 --full \
                 --max-range 400 -o json",
            )
            .unwrap();
        assert!(
            !without.contains("mover_ring"),
            "mover_ring fields must not appear without --target-speed"
        );

        let result = wasm
            .run_command(
                "trajectory --units metric -v 823 -b 0.475 -m 10.9 -d 7.82 --full \
                 --max-range 400 --target-speed 3 -o json",
            )
            .unwrap();
        let json: serde_json::Value = serde_json::from_str(&result).unwrap();
        let points = json["trajectory"].as_array().expect("trajectory array");
        assert!(points.len() > 5);

        // mover_ring_m is EXACTLY target_speed_mps * time at every point (metric run,
        // so --target-speed 3 is already 3 m/s, matching the native CLI convention).
        let mut checked = 0;
        for p in points {
            let time = p["time_seconds"].as_f64().expect("time_seconds");
            let ring_m = p["mover_ring_m"].as_f64().expect("mover_ring_m present");
            let expected = 3.0 * time;
            assert!(
                (ring_m - expected).abs() < 1e-9,
                "ring_m {ring_m} != target_speed_mps*time {expected} at t={time}"
            );
            checked += 1;
        }
        assert!(checked > 5);

        // Muzzle point: mover_ring_m is 0 (present), mover_ring_mil is omitted (0/0
        // is degenerate) — matches the native CLI contract exactly.
        let muzzle = &points[0];
        assert_eq!(muzzle["range_meters"].as_f64().unwrap(), 0.0);
        assert_eq!(muzzle["mover_ring_m"].as_f64().unwrap(), 0.0);
        assert!(muzzle.get("mover_ring_mil").is_none());
    }

    // -----------------------------------------------------------------------------
    // MBA-1325: `lead` powder-temperature parity with `trajectory`.
    // -----------------------------------------------------------------------------

    #[wasm_bindgen_test]
    fn lead_powder_curve_matches_manually_resolved_velocity() {
        let wasm = WasmBallistics::new();
        // WASM's `lead` has no --temperature flag (it always solves at the default
        // atmosphere — a pre-existing gap unrelated to MBA-1325, which only adds the
        // POWDER flags), so drive the curve lookup directly via --powder-temp rather
        // than the ambient air temperature. 55F is the midpoint of the 40/70 curve
        // points, so this is a genuine correction, not a no-op.
        let resolved_json = wasm
            .run_command(
                "trajectory -v 2700 -b 0.19 -m 77 -d 0.224 --drag-model g7 --max-range 5 \
                 --sight-height 2.48 --powder-temp-curve 40:2620,70:2700,100:2760 \
                 --powder-temp 55 --full -o json",
            )
            .unwrap();
        let resolved_json: serde_json::Value = serde_json::from_str(&resolved_json).unwrap();
        let resolved_velocity = resolved_json["trajectory"][0]["velocity_fps"].as_f64().unwrap();
        assert!(
            (resolved_velocity - 2660.0).abs() < 0.5,
            "sanity: 55F should interpolate to ~2660 fps, got {resolved_velocity}"
        );

        let with_curve = wasm
            .run_command(
                "lead -v 2700 -b 0.19 -m 77 -d 0.224 --drag-model g7 --sight-height 2.48 \
                 --powder-temp-curve 40:2620,70:2700,100:2760 --powder-temp 55 \
                 --target-speed 5 --range 300 -o json",
            )
            .unwrap();
        let manual = wasm
            .run_command(&format!(
                "lead -v {resolved_velocity} -b 0.19 -m 77 -d 0.224 --drag-model g7 \
                 --sight-height 2.48 --target-speed 5 --range 300 -o json"
            ))
            .unwrap();

        let with_curve: serde_json::Value = serde_json::from_str(&with_curve).unwrap();
        let manual: serde_json::Value = serde_json::from_str(&manual).unwrap();
        // `resolved_velocity` was already round-tripped once through WASM's imprecise
        // fps<->m/s constants (m/s -> fps via *3.28084, not the exact reciprocal of
        // the *0.3048 used on the way in — see lead_without_powder_flags_uses_
        // velocity_verbatim above); feeding it back in via `-v` round-trips it AGAIN,
        // so the two lead solutions agree to within float noise (~1e-7), not bit-for-
        // bit. Compare numerically with a tolerance instead of assert_eq! on the raw
        // JSON values (iterations is an integer count, compared exactly).
        for field in ["lead_mil", "lead_moa", "lead", "tof_s", "intercept_range"] {
            let a = with_curve[field].as_f64().unwrap_or_else(|| panic!("{field} missing"));
            let b = manual[field].as_f64().unwrap_or_else(|| panic!("{field} missing"));
            assert!(
                (a - b).abs() < 1e-5,
                "field '{field}' differs between curve-corrected and manually-resolved lead runs: {a} vs {b}"
            );
        }
        assert_eq!(
            with_curve["iterations"], manual["iterations"],
            "iterations should match exactly (integer count)"
        );
    }

    // -----------------------------------------------------------------------------
    // MBA-1325 env-flags addendum: `lead` gains --temperature/--pressure/--humidity/
    // --altitude/--wind-speed/--wind-direction (native parity).
    // -----------------------------------------------------------------------------

    /// With NO environmental flags, `lead` output must be byte-identical to the build
    /// from before the flags existed. These literals were captured VERBATIM from the
    /// pre-env-flags build (commit c805229) via a temporary probe test run under
    /// `wasm-pack test --safari --headless` — a true golden-compare, not a
    /// same-code-both-sides tautology. JSON strings have no trailing newline
    /// (serde_json::to_string_pretty); the table format string ends with exactly one.
    #[wasm_bindgen_test]
    fn lead_no_env_flags_output_is_byte_identical_to_pre_env_build() {
        const GOLDEN_IMPERIAL_JSON: &str = r#"{
  "adjustment_unit": "mil",
  "distance_unit": "yd",
  "intercept_range": 399.9987936,
  "iterations": 0,
  "lead": 1.267051911310233,
  "lead_mil": 3.1676393318758076,
  "lead_moa": 10.890344022989026,
  "range": 400.0,
  "target_angle_deg": 90.0,
  "target_speed": 5.0,
  "target_speed_unit": "mph",
  "tof_s": 0.5183409815796777
}"#;
        const GOLDEN_METRIC_JSON: &str = r#"{
  "adjustment_unit": "mil",
  "distance_unit": "m",
  "intercept_range": 350.0,
  "iterations": 0,
  "lead": 1.4774105669002737,
  "lead_mil": 4.221173048286496,
  "lead_moa": 14.512392940008974,
  "range": 350.0,
  "target_angle_deg": 90.0,
  "target_speed": 3.0,
  "target_speed_unit": "m/s",
  "tof_s": 0.49247018896675787
}"#;
        const GOLDEN_IMPERIAL_TABLE: &str = "Moving-Target Lead\n\
===================\n\
Target: 5.0 mph at 90\u{b0} (0=away, 90=left-to-right, 180=toward, 270=right-to-left;\n\
positive lead = hold in direction of travel)\n\
\n\
Range: 400 yd\n\
Time of Flight: 0.518 s\n\
Lead: 1.27 yd (3.17 MIL / 10.89 MOA)\n\
Intercept Range: 400.0 yd\n\
Iterations: 0\n";

        let wasm = WasmBallistics::new();
        let g1 = wasm
            .run_command(
                "lead -v 2700 -b 0.475 -m 168 -d 0.308 --target-speed 5 --range 400 -o json",
            )
            .unwrap();
        assert_eq!(g1, GOLDEN_IMPERIAL_JSON, "imperial JSON drifted from pre-env-flags build");
        let g2 = wasm
            .run_command(
                "--units metric lead -v 823 -b 0.475 -m 10.9 -d 7.82 --target-speed 3 \
                 --range 350 -o json",
            )
            .unwrap();
        assert_eq!(g2, GOLDEN_METRIC_JSON, "metric JSON drifted from pre-env-flags build");
        let g3 = wasm
            .run_command("lead -v 2700 -b 0.475 -m 168 -d 0.308 --target-speed 5 --range 400")
            .unwrap();
        assert_eq!(g3, GOLDEN_IMPERIAL_TABLE, "imperial table drifted from pre-env-flags build");
    }

    /// Hot air (100F) vs the standard default (59F) must change the lead solution —
    /// lower density, less drag, shorter TOF, smaller hold — proving --temperature is
    /// actually plumbed into the solve, not just parsed. (Fails on the pre-addendum
    /// build with "Unknown flag: --temperature".)
    #[wasm_bindgen_test]
    fn lead_temperature_is_live_hot_air_shortens_lead() {
        let wasm = WasmBallistics::new();
        let base = "lead -v 2700 -b 0.475 -m 168 -d 0.308 --target-speed 5 --range 400 -o json";
        let std_out = wasm.run_command(base).unwrap();
        let hot_out = wasm
            .run_command(&format!("{base} --temperature 100"))
            .unwrap();
        let std_json: serde_json::Value = serde_json::from_str(&std_out).unwrap();
        let hot_json: serde_json::Value = serde_json::from_str(&hot_out).unwrap();
        let std_tof = std_json["tof_s"].as_f64().unwrap();
        let hot_tof = hot_json["tof_s"].as_f64().unwrap();
        let std_mil = std_json["lead_mil"].as_f64().unwrap();
        let hot_mil = hot_json["lead_mil"].as_f64().unwrap();
        assert!(
            hot_tof < std_tof,
            "hot air must shorten TOF: hot {hot_tof} vs standard {std_tof}"
        );
        assert!(
            hot_mil < std_mil && (std_mil - hot_mil) / std_mil > 1e-4,
            "hot air must shrink the hold by more than float noise: hot {hot_mil} vs standard {std_mil}"
        );
    }

    /// A 20 mph headwind (--wind-direction 0 = wind-FROM dead ahead) adds drag along
    /// the flight path, lengthening TOF and growing the hold vs calm air — proving the
    /// wind flags reach the WindConditions the solve actually uses.
    #[wasm_bindgen_test]
    fn lead_headwind_is_live_and_lengthens_tof() {
        let wasm = WasmBallistics::new();
        let base = "lead -v 2700 -b 0.475 -m 168 -d 0.308 --target-speed 5 --range 400 -o json";
        let calm_out = wasm.run_command(base).unwrap();
        let wind_out = wasm
            .run_command(&format!("{base} --wind-speed 20 --wind-direction 0"))
            .unwrap();
        let calm_json: serde_json::Value = serde_json::from_str(&calm_out).unwrap();
        let wind_json: serde_json::Value = serde_json::from_str(&wind_out).unwrap();
        let calm_tof = calm_json["tof_s"].as_f64().unwrap();
        let wind_tof = wind_json["tof_s"].as_f64().unwrap();
        assert!(
            wind_tof > calm_tof && (wind_tof - calm_tof) / calm_tof > 1e-4,
            "a 20 mph headwind must lengthen TOF beyond float noise: wind {wind_tof} vs calm {calm_tof}"
        );
    }

    #[wasm_bindgen_test]
    fn lead_without_powder_flags_uses_velocity_verbatim() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command(
                "trajectory -v 2700 -b 0.19 -m 77 -d 0.224 --drag-model g7 --max-range 5 \
                 --sight-height 2.48 --temperature 100 --full -o json",
            )
            .unwrap();
        let json: serde_json::Value = serde_json::from_str(&result).unwrap();
        let velocity = json["trajectory"][0]["velocity_fps"].as_f64().unwrap();
        // Tolerance, not exact equality: the WASM imperial round-trip (m/s -> fps via
        // a fixed *3.28084, not the exact reciprocal of the fps -> m/s *0.3048 used on
        // the way in) carries an inherent ~3e-8 relative bias (~0.0001 fps at 2700),
        // pre-existing and unrelated to MBA-1325. 1e-3 is >100x looser than that noise
        // floor while still being >1000x tighter than any real powder correction
        // (which is many fps), so this still catches a genuine verbatim-velocity
        // regression.
        assert!(
            (velocity - 2700.0).abs() < 1e-3,
            "no powder flags should keep -v 2700 verbatim regardless of --temperature, got {velocity}"
        );
    }

    // -----------------------------------------------------------------------------
    // MBA-1325 review fix pass: M1 (target-speed bound), M2 (curve parser parity),
    // M3 (ring table --adjustment-unit).
    // -----------------------------------------------------------------------------

    /// M1: --target-speed enforces the same [0, 300] bound the native flags carry via
    /// f64_range(0.0, 300.0), on BOTH WASM commands — rejected with a clear error, not
    /// silently clamped or let through.
    #[wasm_bindgen_test]
    fn target_speed_out_of_range_is_rejected_on_both_commands() {
        let wasm = WasmBallistics::new();
        let err = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --target-speed 301")
            .unwrap_err();
        assert!(
            err.as_string().unwrap_or_default().contains("--target-speed must be between 0 and 300"),
            "trajectory should reject 301: {err:?}"
        );
        let err = wasm
            .run_command(
                "lead -v 2700 -b 0.475 -m 168 -d 0.308 --target-speed 1000000000 --range 400",
            )
            .unwrap_err();
        assert!(
            err.as_string().unwrap_or_default().contains("--target-speed must be between 0 and 300"),
            "lead should reject 1e9 mph: {err:?}"
        );
        // Boundary values stay legal (300 is lead's documented cap).
        assert!(wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --target-speed 300 --max-range 200")
            .is_ok());
    }

    /// M2: duplicate temperatures in --powder-temp-curve are rejected on both WASM
    /// commands with the native parser's message (they now share one parser, so this
    /// also guards against the two copies drifting again).
    #[wasm_bindgen_test]
    fn powder_temp_curve_duplicate_temperatures_rejected() {
        let wasm = WasmBallistics::new();
        for cmd in [
            "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 5 \
             --powder-temp-curve 40:2620,40:2700,100:2760",
            "lead -v 2700 -b 0.475 -m 168 -d 0.308 --target-speed 5 --range 300 \
             --powder-temp-curve 40:2620,40:2700,100:2760",
        ] {
            let err = wasm.run_command(cmd).unwrap_err();
            assert!(
                err.as_string()
                    .unwrap_or_default()
                    .contains("--powder-temp-curve has duplicate temperatures"),
                "duplicate curve temps must be rejected for `{cmd}`: {err:?}"
            );
        }
    }

    /// M3: the trajectory table's Ring cell honors --adjustment-unit, and the MOA
    /// figure is exactly mil x 3.438 (the crate's MBA-724 printed-table dial constant,
    /// shared with the native table). --target-speed 300 makes ring_mil large enough
    /// (~190) that a wrong constant (exact-angle 3437.7467/1000) would shift the value
    /// by > 0.012 — past the 0.01 print quantum, so the 2dp string match below pins
    /// the constant, not just the ballpark.
    #[wasm_bindgen_test]
    fn ring_table_honors_moa_adjustment_unit() {
        let wasm = WasmBallistics::new();
        let base = "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 400 \
                    --target-speed 300 --full";

        // Full-precision mil from JSON (adjustment-unit-invariant by contract).
        let json_out = wasm.run_command(&format!("{base} -o json")).unwrap();
        let json: serde_json::Value = serde_json::from_str(&json_out).unwrap();
        let mil = json["trajectory"]
            .as_array()
            .expect("points")
            .last()
            .expect("last point")["mover_ring_mil"]
            .as_f64()
            .expect("mover_ring_mil");
        assert!(mil > 160.0, "sensitivity precondition: ring_mil > 160, got {mil}");

        // Last ring cell of a table, as the printed numeric string before " <unit>".
        let last_ring_cell = |table: &str, unit: &str| -> String {
            let suffix = format!(" {unit}");
            table
                .lines()
                .filter_map(|line| {
                    let cell = line.rsplit('|').next()?.trim();
                    cell.strip_suffix(suffix.as_str()).map(str::to_string)
                })
                .last()
                .unwrap_or_else(|| panic!("no ring cell ending in '{suffix}' found:\n{table}"))
        };

        let moa_table = wasm
            .run_command(&format!("{base} --adjustment-unit moa"))
            .unwrap();
        assert!(moa_table.contains(" moa"), "moa run must label ring cells moa");
        assert!(!moa_table.contains(" mil"), "moa run must not label ring cells mil");
        assert_eq!(
            last_ring_cell(&moa_table, "moa"),
            format!("{:.2}", mil * 3.438),
            "Ring moa cell must be exactly mil x 3.438, 2dp-rounded"
        );

        let mil_table = wasm.run_command(base).unwrap();
        assert_eq!(
            last_ring_cell(&mil_table, "mil"),
            format!("{:.2}", mil),
            "default (mil) Ring cell must match JSON mover_ring_mil, 2dp-rounded"
        );
    }

    // -----------------------------------------------------------------------------
    // MBA-1355: WASM terminal parity for smoa/iphy/clicks turret adjustment units
    // (native CLI got these in main.rs; this brings `trajectory`'s Ring column and
    // `lead`'s adjustment display up to the same surface). The real flag is
    // `--adjustment-unit` (NOT `--units`, which selects imperial/metric).
    // -----------------------------------------------------------------------------

    /// Ring column accepts smoa/iphy (real constant-factor conversion, same
    /// smoa_per_mil() ratio as moa/mil) and clicks (whole-integer rounding via a
    /// resolved elevation click graduation) — mirroring native's Ring(smoa)/Ring(iphy)/
    /// Ring(clicks) headers exactly (CLI_USAGE.md's documented "(mil)/(moa) convention
    /// ... e.g. the mover Ring column reads Ring(clicks)"). A bare --adjustment-unit
    /// clicks (no graduation from any source — WASM has no --profile) must fail fast,
    /// naming --elevation-click-value.
    #[wasm_bindgen_test]
    fn test_units_smoa_and_clicks_accepted() {
        let base = "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 300 --target-speed 3";

        let smoa = WasmBallistics::new()
            .run_command(&format!("{base} --adjustment-unit smoa"))
            .unwrap();
        assert!(smoa.contains("Ring(smoa)"), "{smoa}");

        let iphy = WasmBallistics::new()
            .run_command(&format!("{base} --adjustment-unit iphy"))
            .unwrap();
        assert!(iphy.contains("Ring(iphy)"), "{iphy}");
        // smoa and iphy are numerically identical (both 3600 factor) — only the label differs.
        let smoa_json = WasmBallistics::new()
            .run_command(&format!("{base} --adjustment-unit smoa -o json"))
            .unwrap();
        let iphy_json = WasmBallistics::new()
            .run_command(&format!("{base} --adjustment-unit iphy -o json"))
            .unwrap();
        assert_eq!(smoa_json, iphy_json, "smoa/iphy JSON is unit-in-name-only (mil-based fields)");

        let clicks = WasmBallistics::new()
            .run_command(&format!("{base} --adjustment-unit clicks --elevation-click-value 0.25moa"))
            .unwrap();
        assert!(clicks.contains("Ring(clicks)"), "{clicks}");
        // Clicks cells are whole integers with no per-cell unit suffix (native's
        // RingUnit::Clicks arm: `format!("{:>8}", clicks_for(...))`, no trailing label) —
        // unlike the Factor arms, which append " mil"/" moa"/" smoa"/" iphy" per cell.
        let last_ring_cell = clicks
            .lines()
            .filter(|l| l.contains('|') && !l.starts_with("Range") && !l.starts_with("---"))
            .last()
            .map(|l| l.rsplit('|').next().unwrap().trim().to_string())
            .unwrap_or_default();
        assert!(
            last_ring_cell.parse::<i64>().is_ok(),
            "clicks Ring cell must be a bare whole integer, got {last_ring_cell:?} in:\n{clicks}"
        );

        let missing = WasmBallistics::new()
            .run_command(&format!("{base} --adjustment-unit clicks"))
            .unwrap_err();
        let msg = missing.as_string().unwrap_or_default();
        assert!(
            msg.contains("--adjustment-unit clicks requires a turret elevation graduation"),
            "{msg}"
        );
        assert!(msg.contains("--elevation-click-value"), "{msg}");
    }

    /// A malformed --elevation-click-value (missing unit suffix, non-positive
    /// magnitude, etc.) surfaces parse_click_value's own error message, same parser
    /// the native CLI uses (src/adjustment.rs, shared by both binaries).
    #[wasm_bindgen_test]
    fn clicks_rejects_malformed_elevation_click_value() {
        let err = WasmBallistics::new()
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 300 --target-speed 3 \
                 --adjustment-unit clicks --elevation-click-value notaclick",
            )
            .unwrap_err();
        let msg = err.as_string().unwrap_or_default();
        assert!(msg.contains("needs a unit suffix"), "{msg}");
    }

    /// --windage-click-value is validated even though the Ring column never displays
    /// it (it always reuses the elevation graduation, same "accepted but inert"
    /// contract native's come-ups documents for its own windage flag) — a typo'd
    /// windage value must still error, not be silently ignored.
    #[wasm_bindgen_test]
    fn clicks_validates_windage_click_value_even_though_unused_by_ring() {
        let err = WasmBallistics::new()
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 300 --target-speed 3 \
                 --adjustment-unit clicks --elevation-click-value 0.25moa \
                 --windage-click-value bogus",
            )
            .unwrap_err();
        let msg = err.as_string().unwrap_or_default();
        assert!(msg.contains("needs a unit suffix"), "{msg}");

        // A well-formed windage value alongside elevation still resolves normally.
        let ok = WasmBallistics::new()
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 300 --target-speed 3 \
                 --adjustment-unit clicks --elevation-click-value 0.25moa \
                 --windage-click-value 0.1mil",
            )
            .unwrap();
        assert!(ok.contains("Ring(clicks)"), "{ok}");
    }

    /// An unrecognized --adjustment-unit value lists all five accepted choices, not
    /// just the pre-MBA-1355 mil/moa pair.
    #[wasm_bindgen_test]
    fn invalid_adjustment_unit_lists_all_five_choices() {
        let err = WasmBallistics::new()
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 300 --target-speed 3 \
                 --adjustment-unit bogus",
            )
            .unwrap_err();
        let msg = err.as_string().unwrap_or_default();
        assert!(msg.contains("mil, moa, smoa, iphy, or clicks"), "{msg}");
    }

    /// CSV/JSON ring fields stay mil-only regardless of --adjustment-unit (pre-existing
    /// contract, unaffected by MBA-1355): ring_mil (csv) / mover_ring_mil (json) are
    /// always in mil even when the table's Ring column is requesting clicks.
    #[wasm_bindgen_test]
    fn ring_csv_and_json_fields_stay_mil_only_under_clicks() {
        let base = "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 300 --target-speed 3 \
                    --adjustment-unit clicks --elevation-click-value 0.25moa";

        let csv = WasmBallistics::new().run_command(&format!("{base} -o csv")).unwrap();
        assert!(csv.contains("Ring(mil)"), "{csv}");

        let json = WasmBallistics::new().run_command(&format!("{base} -o json")).unwrap();
        let json: serde_json::Value = serde_json::from_str(&json).unwrap();
        let points = json["trajectory"].as_array().expect("trajectory array");
        assert!(points.iter().any(|p| p.get("mover_ring_mil").is_some()));
    }

    /// `lead` gains real smoa/iphy display (mirrors native handle_lead's
    /// Smoa|Iphy => sol.lead_mil * smoa_per_mil() arm). MBA-1410: clicks now resolves too,
    /// via --windage-click-value (lead's single column is windage-like) -- see
    /// `lead_accepts_clicks_with_windage_click_value` below.
    #[wasm_bindgen_test]
    fn lead_accepts_smoa_and_iphy_with_real_values() {
        let smoa_json = WasmBallistics::new()
            .run_command(
                "lead -v 2700 -b 0.475 -m 168 -d 0.308 --target-speed 5 --range 300 \
                 --adjustment-unit smoa -o json",
            )
            .unwrap();
        let json: serde_json::Value = serde_json::from_str(&smoa_json).unwrap();
        let lead_mil = json["lead_mil"].as_f64().unwrap();
        let lead_smoa = json["lead_smoa"].as_f64().unwrap();
        assert!(
            (lead_smoa - lead_mil * 3.6).abs() < 1e-6,
            "lead_smoa must be lead_mil * 3.6 exactly: {lead_smoa} vs {}",
            lead_mil * 3.6
        );
        assert_eq!(json["adjustment_unit"], "smoa");

        let smoa_table = WasmBallistics::new()
            .run_command(
                "lead -v 2700 -b 0.475 -m 168 -d 0.308 --target-speed 5 --range 300 \
                 --adjustment-unit smoa",
            )
            .unwrap();
        assert!(smoa_table.contains("SMOA"), "{smoa_table}");

        let iphy_table = WasmBallistics::new()
            .run_command(
                "lead -v 2700 -b 0.475 -m 168 -d 0.308 --target-speed 5 --range 300 \
                 --adjustment-unit iphy",
            )
            .unwrap();
        assert!(iphy_table.contains("IPHY"), "{iphy_table}");
    }

    /// MBA-1410: `lead --adjustment-unit clicks` now resolves for real via
    /// --windage-click-value (mirrors native handle_lead's Clicks arm, which resolves
    /// through the WINDAGE graduation since lead's single column is windage-like). A bare
    /// `--adjustment-unit clicks` with no graduation (WASM has no --profile, so there is
    /// no fallback source) still fails fast, naming --windage-click-value.
    #[wasm_bindgen_test]
    fn lead_accepts_clicks_with_windage_click_value() {
        let json = WasmBallistics::new()
            .run_command(
                "lead -v 2700 -b 0.475 -m 168 -d 0.308 --target-speed 5 --range 300 \
                 --adjustment-unit clicks --windage-click-value 0.25moa -o json",
            )
            .unwrap();
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(v["adjustment_unit"], "clicks");
        assert!(v["lead_clicks"].as_i64().is_some(), "{json}");

        let table = WasmBallistics::new()
            .run_command(
                "lead -v 2700 -b 0.475 -m 168 -d 0.308 --target-speed 5 --range 300 \
                 --adjustment-unit clicks --windage-click-value 0.25moa",
            )
            .unwrap();
        assert!(table.contains("clicks"), "{table}");

        let missing = WasmBallistics::new()
            .run_command(
                "lead -v 2700 -b 0.475 -m 168 -d 0.308 --target-speed 5 --range 300 \
                 --adjustment-unit clicks",
            )
            .unwrap_err();
        let msg = missing.as_string().unwrap_or_default();
        assert!(msg.contains("--windage-click-value"), "{msg}");
    }

    // -----------------------------------------------------------------------------
    // MBA-1328: custom Mach:Cd drag table, mirroring the BC5D bytes-loader pattern
    // (load_bc5d_table / has_bc5d_table above). loadDragTable(bytes) parses the SAME
    // CSV format native --drag-table accepts (DragTable::from_csv_str) and, once
    // loaded, is applied automatically to every trajectory/zero/lead/monte-carlo run
    // — no --use-* gate flag needed (unlike BC5D, which needs --use-bc-segments).
    // -----------------------------------------------------------------------------

    /// A deliberately non-G1/non-G7 curve: flat Cd 0.5 across the whole Mach range.
    /// 6 points, comfortably above try_new's 2-point minimum, spanning Mach 0-3.
    const FLAT_CD_CSV: &str = "mach,cd\n0.0,0.5\n0.5,0.5\n1.0,0.5\n1.5,0.5\n2.0,0.5\n3.0,0.5\n";

    #[wasm_bindgen_test]
    fn drag_table_not_loaded_by_default() {
        let wasm = WasmBallistics::new();
        assert!(!wasm.has_drag_table());
    }

    #[wasm_bindgen_test]
    fn load_drag_table_reports_point_count_and_mach_range() {
        let wasm = WasmBallistics::new();
        let summary = wasm.load_drag_table(FLAT_CD_CSV.as_bytes()).unwrap();
        assert!(wasm.has_drag_table());
        assert!(summary.contains("6 points"), "got: {summary}");
        assert!(summary.contains("0.000-3.000"), "got: {summary}");
    }

    #[wasm_bindgen_test]
    fn load_drag_table_replaces_previous_table() {
        let wasm = WasmBallistics::new();
        wasm.load_drag_table(FLAT_CD_CSV.as_bytes()).unwrap();
        assert!(wasm.has_drag_table());
        // A second, different (2-point) table must fully replace the first, not merge
        // with or append to it.
        let summary = wasm
            .load_drag_table(b"mach,cd\n0.5,0.30\n2.5,0.25\n")
            .unwrap();
        assert!(wasm.has_drag_table());
        assert!(summary.contains("2 points"), "got: {summary}");
    }

    /// Live-change: a table loaded via loadDragTable() must actually change the
    /// physics of a trajectory run, not just be stored inertly. Compares the flat,
    /// deliberately non-G7-shaped Cd curve above against the unloaded (G1-model) run
    /// of otherwise-identical args.
    #[wasm_bindgen_test]
    fn drag_table_live_change_alters_trajectory_output() {
        let command = "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 100 -o json";

        let unloaded = WasmBallistics::new();
        let baseline = unloaded.run_command(command).unwrap();

        let loaded = WasmBallistics::new();
        loaded.load_drag_table(FLAT_CD_CSV.as_bytes()).unwrap();
        let with_table = loaded.run_command(command).unwrap();

        assert_ne!(
            baseline, with_table,
            "loading a drag table must change trajectory output for identical args"
        );
    }

    /// Live-change also applies to `zero`, `lead`, and `monte-carlo` (deliberately
    /// broader than the BC5D auto-apply, which only wires into `trajectory` — see
    /// handle_trajectory_command's BC5D block vs the unconditional drag_table checks
    /// added to all four handlers).
    #[wasm_bindgen_test]
    fn drag_table_live_change_alters_zero_lead_and_monte_carlo() {
        let zero_cmd = "zero -v 2700 -b 0.475 -m 168 -d 0.308 --target-distance 300";
        let lead_cmd =
            "lead -v 2700 -b 0.475 -m 168 -d 0.308 --target-speed 10 --range 300 -o json";
        let mc_cmd = "monte-carlo -v 2700 -b 0.475 -m 168 -d 0.308 -n 50";

        for command in [zero_cmd, lead_cmd, mc_cmd] {
            let unloaded = WasmBallistics::new();
            let baseline = unloaded.run_command(command).unwrap();

            let loaded = WasmBallistics::new();
            loaded.load_drag_table(FLAT_CD_CSV.as_bytes()).unwrap();
            let with_table = loaded.run_command(command).unwrap();

            assert_ne!(
                baseline, with_table,
                "loading a drag table must change `{command}` output"
            );
        }
    }

    /// Clearing reverts the physics: after load → clear, a run must be byte-identical to a
    /// never-loaded instance's run of the same args — not merely "different from the CDM
    /// run". This is exactly the property Bero needs for a G7-vs-CDM comparison on ONE
    /// instance (load → solve CDM → clear → solve G7), the inverse of the live-change test
    /// above.
    #[wasm_bindgen_test]
    fn clear_drag_table_reverts_trajectory_output_to_g7_baseline() {
        let command = "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 100 -o json";

        let baseline = WasmBallistics::new().run_command(command).unwrap();

        let engine = WasmBallistics::new();
        engine.load_drag_table(FLAT_CD_CSV.as_bytes()).unwrap();
        let with_table = engine.run_command(command).unwrap();
        assert_ne!(
            baseline, with_table,
            "sanity: the loaded CDM run must differ from the G7 baseline"
        );

        assert!(
            engine.clear_drag_table(),
            "clear must report it unloaded a table"
        );
        assert!(!engine.has_drag_table());
        let after_clear = engine.run_command(command).unwrap();
        assert_eq!(
            baseline, after_clear,
            "after clearDragTable the solve must revert byte-for-byte to the G7-BC baseline"
        );
    }

    /// Return value + idempotence: clear reports whether it actually unloaded a table, and
    /// a second clear (or a clear on a fresh instance) is a harmless `false` no-op.
    #[wasm_bindgen_test]
    fn clear_drag_table_reports_prior_state_and_is_idempotent() {
        let engine = WasmBallistics::new();
        assert!(
            !engine.clear_drag_table(),
            "no table loaded → clear returns false"
        );

        engine.load_drag_table(FLAT_CD_CSV.as_bytes()).unwrap();
        assert!(engine.has_drag_table());
        assert!(
            engine.clear_drag_table(),
            "a loaded table → clear returns true"
        );
        assert!(
            !engine.has_drag_table(),
            "after clear, no table remains loaded"
        );
        assert!(
            !engine.clear_drag_table(),
            "a second clear is an idempotent false no-op"
        );
    }

    /// Golden-unloaded: with NO table loaded, trajectory output for a fixed command
    /// must remain byte-identical to a literal captured from the pre-MBA-1328 build
    /// (harvested via a temporary probe test run under `wasm-pack test --node`,
    /// following the MBA-1325 golden-compare methodology — see
    /// `lead_no_env_flags_output_is_byte_identical_to_pre_env_build` above). Guards
    /// against the new drag_table plumbing silently perturbing the no-table path.
    #[wasm_bindgen_test]
    fn drag_table_unloaded_output_is_byte_identical_to_pre_change_build() {
        const GOLDEN: &str = "Trajectory Calculation Results\n\
==============================\n\
\n\
Range | Drop | Drift | Velocity | Energy | Time\n\
------|------|-------|----------|--------|------\n\
0 yd   | 2.0 in | 0.0 in  | 2700 fps   | 2719 ft-lb | 0.000 s\n\
55 yd  | 2.7 in | 0.0 in  | 2595 fps   | 2512 ft-lb | 0.062 s\n\
100 yd | 4.5 in | 0.0 in  | 2510 fps   | 2350 ft-lb | 0.115 s\n\
\n\
Max Range: 100 yards\n\
Max Height: 60.0 inches\n\
Time of Flight: 0.12 seconds\n\
Impact Velocity: 2510 fps\n";

        let wasm = WasmBallistics::new();
        assert!(!wasm.has_drag_table());
        let out = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 100")
            .unwrap();
        assert_eq!(
            out, GOLDEN,
            "unloaded trajectory output drifted from the pre-MBA-1328 build"
        );
    }

    #[wasm_bindgen_test]
    fn load_drag_table_rejects_malformed_csv_with_parser_message() {
        let wasm = WasmBallistics::new();
        let err = wasm
            .load_drag_table(b"0.5,0.23\n1.0,notanumber\n")
            .unwrap_err();
        let msg = err.as_string().unwrap_or_default();
        assert!(
            msg.contains("line 2"),
            "expected the DragTable::from_csv_str parser's line-numbered message, got: {msg}"
        );
        assert!(
            !wasm.has_drag_table(),
            "a failed load must not leave a table installed"
        );
    }

    #[wasm_bindgen_test]
    fn load_drag_table_rejects_invalid_utf8_cleanly() {
        let wasm = WasmBallistics::new();
        // 0xFF is never a valid UTF-8 lead byte.
        let err = wasm
            .load_drag_table(&[0xFF, 0xFE, 0x00, 0x00])
            .unwrap_err();
        let msg = err.as_string().unwrap_or_default();
        assert!(
            msg.contains("UTF-8"),
            "expected a clean UTF-8 decode error, got: {msg}"
        );
        assert!(!wasm.has_drag_table());
    }

    /// Deliverable #4: the WASM surface never hard-requires `-b`/`--bc` (every handler
    /// defaults it, e.g. `let mut bc = 0.475;` in handle_trajectory_command) — unlike native
    /// `zero`/`monte-carlo`, which clap-require it even though it's ignored once a table is
    /// active (CLI_USAGE.md). Confirms the solve path actually succeeds end-to-end with a
    /// loaded table and the default bc, omitting `-b` entirely, across all four commands.
    #[wasm_bindgen_test]
    fn drag_table_solve_succeeds_without_explicit_bc_flag() {
        let wasm = WasmBallistics::new();
        wasm.load_drag_table(FLAT_CD_CSV.as_bytes()).unwrap();

        let trajectory = wasm
            .run_command("trajectory -v 2700 -m 168 -d 0.308 --max-range 100")
            .unwrap();
        assert!(trajectory.contains("Trajectory Calculation Results"));

        let zero = wasm
            .run_command("zero -v 2700 -m 168 -d 0.308 --target-distance 300")
            .unwrap();
        assert!(zero.contains("Zero Calculation Results"));

        let lead = wasm
            .run_command("lead -v 2700 -m 168 -d 0.308 --target-speed 10 --range 300")
            .unwrap();
        assert!(lead.contains("Moving-Target Lead"));

        let monte_carlo = wasm
            .run_command("monte-carlo -v 2700 -m 168 -d 0.308 -n 50")
            .unwrap();
        assert!(monte_carlo.contains("Monte Carlo Simulation Results"));
    }

    // -----------------------------------------------------------------------------
    // MBA-1356: `--cd-scale <FACTOR>` — whole-curve drag scale for a loaded custom drag table.
    // WASM parity with the native CLI's `--cd-scale` on trajectory/zero/monte-carlo: same
    // pairing requirement with a loaded table (Err(JsValue) here instead of exit(1)), same
    // out-of-range nudge text, prepended table-only. 0.28.1 sweep: unlike native (whose
    // `trajectory` has a `--bc-adjustment` flag to suggest), the WASM terminal has no
    // `--bc-adjustment` flag on ANY surface, so its pairing text never suggests one.
    // -----------------------------------------------------------------------------

    /// Without a loaded drag table, `--cd-scale` must be rejected before any solve. Unlike
    /// native's `trajectory` pairing error, this never suggests `--bc-adjustment` -- the
    /// WASM terminal has no such flag on any surface (0.28.1 sweep).
    #[wasm_bindgen_test]
    fn cd_scale_without_a_loaded_drag_table_errors_on_every_command() {
        for command in [
            "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 100 --cd-scale 1.1",
            "zero -v 2700 -b 0.475 -m 168 -d 0.308 --target-distance 300 --cd-scale 1.1",
            "monte-carlo -v 2700 -b 0.475 -m 168 -d 0.308 -n 50 --cd-scale 1.1",
        ] {
            let err = WasmBallistics::new().run_command(command).unwrap_err();
            let msg = err.as_string().unwrap_or_default();
            assert_eq!(
                msg,
                "--cd-scale requires --drag-table",
                "command: {command}"
            );
        }
    }

    /// The same pairing requirement applies to the `--wez` sweep path (a separate compute
    /// core, `crate::wez::compute_wez`, from the base monte-carlo path above).
    #[wasm_bindgen_test]
    fn cd_scale_without_a_loaded_drag_table_errors_on_wez() {
        let err = WasmBallistics::new()
            .run_command(
                "monte-carlo -v 2700 -b 0.475 -m 168 -d 0.308 --wez --target-size 18x30 \
                 -n 20 --wez-start 200 --wez-end 300 --wez-step 100 --cd-scale 1.1",
            )
            .unwrap_err();
        let msg = err.as_string().unwrap_or_default();
        assert_eq!(msg, "--cd-scale requires --drag-table");
    }

    /// With a table loaded, `--cd-scale` is accepted and actually shifts trajectory output
    /// (the plumb-through smoke): a higher scale means more drag, so a `--cd-scale 1.1` run
    /// must differ from the neutral `1.0` run of otherwise-identical args.
    #[wasm_bindgen_test]
    fn cd_scale_accepted_with_a_loaded_table_and_shifts_trajectory_output() {
        let wasm = WasmBallistics::new();
        wasm.load_drag_table(FLAT_CD_CSV.as_bytes()).unwrap();

        let neutral = wasm
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 100 --cd-scale 1.0 -o json",
            )
            .unwrap();
        let scaled = wasm
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 100 --cd-scale 1.1 -o json",
            )
            .unwrap();

        assert_ne!(
            neutral, scaled,
            "--cd-scale 1.1 must change trajectory output relative to the neutral 1.0"
        );

        // Omitting --cd-scale entirely must be byte-identical to the explicit neutral 1.0 (the
        // engine default), on the SAME loaded table.
        let omitted = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 100 -o json")
            .unwrap();
        assert_eq!(
            omitted, neutral,
            "omitting --cd-scale must be byte-identical to an explicit --cd-scale 1.0"
        );
    }

    /// Same shift, on `zero` and `monte-carlo` (the other two commands `--cd-scale` is wired
    /// into), and on the WEZ sweep path — a separate compute core from the other three.
    #[wasm_bindgen_test]
    fn cd_scale_shifts_zero_and_monte_carlo_and_wez_output() {
        let zero_cmd_neutral =
            "zero -v 2700 -b 0.475 -m 168 -d 0.308 --target-distance 300 --cd-scale 1.0";
        let zero_cmd_scaled =
            "zero -v 2700 -b 0.475 -m 168 -d 0.308 --target-distance 300 --cd-scale 1.1";
        let mc_cmd_neutral = "monte-carlo -v 2700 -b 0.475 -m 168 -d 0.308 -n 50 --cd-scale 1.0";
        let mc_cmd_scaled = "monte-carlo -v 2700 -b 0.475 -m 168 -d 0.308 -n 50 --cd-scale 1.1";
        let wez_cmd_neutral = "monte-carlo -v 2700 -b 0.475 -m 168 -d 0.308 --wez \
             --target-size 18x30 -n 20 --wez-start 200 --wez-end 300 --wez-step 100 \
             --cd-scale 1.0";
        let wez_cmd_scaled = "monte-carlo -v 2700 -b 0.475 -m 168 -d 0.308 --wez \
             --target-size 18x30 -n 20 --wez-start 200 --wez-end 300 --wez-step 100 \
             --cd-scale 1.1";

        for (neutral_cmd, scaled_cmd) in [
            (zero_cmd_neutral, zero_cmd_scaled),
            (mc_cmd_neutral, mc_cmd_scaled),
            (wez_cmd_neutral, wez_cmd_scaled),
        ] {
            let wasm = WasmBallistics::new();
            wasm.load_drag_table(FLAT_CD_CSV.as_bytes()).unwrap();
            let neutral = wasm.run_command(neutral_cmd).unwrap();
            let scaled = wasm.run_command(scaled_cmd).unwrap();
            assert_ne!(
                neutral, scaled,
                "`{scaled_cmd}` must differ from `{neutral_cmd}`"
            );
        }
    }

    /// A `--cd-scale` far outside the typical truing range (outside `[0.5, 2.0]`) must still
    /// solve successfully (the engine's own gate is only finite && > 0) but warn once, table-only
    /// — the warning text must not appear in JSON output.
    #[wasm_bindgen_test]
    fn cd_scale_out_of_range_warns_table_only() {
        let wasm = WasmBallistics::new();
        wasm.load_drag_table(FLAT_CD_CSV.as_bytes()).unwrap();

        let table_output = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 100 --cd-scale 3.0")
            .unwrap();
        assert!(
            table_output.contains("--cd-scale 3 is far outside the typical truing range (0.90-1.10)"),
            "got: {table_output}"
        );

        let json_output = wasm
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 100 --cd-scale 3.0 -o json",
            )
            .unwrap();
        assert!(
            !json_output.contains("far outside the typical truing range"),
            "the range warning must never contaminate JSON output, got: {json_output}"
        );
    }

    /// MBA-1366 review, Important #1: `--density-altitude` derives an ABSOLUTE station
    /// pressure, so the declared `--pressure-type` no longer describes it. The shot-day
    /// atmosphere was reset correctly but the mode was not, leaving the zero-day solve to
    /// inherit a stale `Qnh` and reduce an already-absolute pressure a second time (~7%
    /// pressure error at 2000 ft DA — shot day right, zero day silently wrong). Native
    /// always reset the mode; the browser terminal did not mirror it.
    #[wasm_bindgen_test]
    fn density_altitude_supersedes_the_pressure_mode_on_the_zero_day_too() {
        let wasm = WasmBallistics::new();
        const BASE: &str = "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 500 \
                            --density-altitude 2000 --auto-zero 100 -o json";

        // Declaring a QNH pressure alongside --density-altitude must change nothing at all:
        // density altitude supersedes both the value and the mode, on both solves.
        let with_qnh = wasm
            .run_command(&format!("{BASE} --pressure 1030 --pressure-type qnh"))
            .unwrap();
        let da_only = wasm.run_command(BASE).unwrap();
        assert_eq!(
            with_qnh, da_only,
            "--density-altitude must fully supersede --pressure/--pressure-type; a residual \
             mode would re-reduce the derived absolute pressure on the zero-day solve"
        );
    }

    /// MBA-1414: on the browser terminal `--adjustment-unit` reaches only the mover Ring
    /// column, which exists only for a nonzero `--target-speed`. A tester set it, got a
    /// byte-identical table and no signal at all — so the inert case now says so, table-only.
    #[wasm_bindgen_test]
    fn inert_adjustment_unit_warns_table_only() {
        let wasm = WasmBallistics::new();
        const BASE: &str = "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 300";

        // Supplied with no mover ring: warn, and leave the table below it unchanged.
        let inert = wasm
            .run_command(&format!("{BASE} --adjustment-unit moa"))
            .unwrap();
        assert!(
            inert.contains("--adjustment-unit only affects the mover Ring column"),
            "got: {inert}"
        );

        // A mover ring consumes the unit: no warning, and the column carries it.
        let with_ring = wasm
            .run_command(&format!("{BASE} --adjustment-unit moa --target-speed 10"))
            .unwrap();
        assert!(
            !with_ring.contains("--adjustment-unit only affects"),
            "a rendered Ring column means the flag did something: {with_ring}"
        );
        assert!(with_ring.contains("Ring(moa)"), "got: {with_ring}");

        // Never supplied: silence (the default unit is not a user choice we can detect).
        let untouched = wasm.run_command(BASE).unwrap();
        assert!(
            !untouched.contains("--adjustment-unit only affects"),
            "an unsupplied flag must not warn: {untouched}"
        );

        // Machine formats stay pure.
        for fmt in ["json", "csv"] {
            let machine = wasm
                .run_command(&format!("{BASE} --adjustment-unit moa -o {fmt}"))
                .unwrap();
            assert!(
                !machine.contains("--adjustment-unit only affects"),
                "the warning must never contaminate {fmt} output, got: {machine}"
            );
        }
    }

    // -----------------------------------------------------------------------------
    // MBA-1409: `loadDragTable` accepts `.drg` vendor drag-curve text as a fallback when the
    // bytes don't parse as CSV. WASM has no filesystem and thus no file extension to key off
    // (unlike native `--drag-table`, which dispatches on the `.drg`/other-extension split), so
    // the fallback signal is `looks_like_drg(text)` on the CSV-parse-failure text itself.
    // -----------------------------------------------------------------------------

    /// A synthetic .drg deck (invented values; no vendor data) with a leading name/header line
    /// and tab-separated (mach, cd) rows -- structurally like the real vendor format. Its points
    /// are deliberately identical to `FLAT_CD_CSV` above so a solve can be compared byte-for-byte
    /// against the plain-CSV path.
    const DRG_TEXT_SAME_POINTS_AS_FLAT_CD_CSV: &str =
        "SYNTH TEST DECK, invented values (MBA-1409 WASM test)\n\
         0.00\t0.5\n\
         0.50\t0.5\n\
         1.00\t0.5\n\
         1.50\t0.5\n\
         2.00\t0.5\n\
         3.00\t0.5\n";

    #[wasm_bindgen_test]
    fn test_load_drag_table_accepts_drg_text() {
        let wasm = WasmBallistics::new();
        let summary = wasm
            .load_drag_table(DRG_TEXT_SAME_POINTS_AS_FLAT_CD_CSV.as_bytes())
            .expect(".drg-shaped text should be accepted as a fallback");
        assert!(wasm.has_drag_table());
        assert!(summary.contains("6 points"), "got: {summary}");
        assert!(summary.contains("0.000-3.000"), "got: {summary}");

        // The .drg fallback must drive a trajectory run exactly as the CSV path does for the
        // same underlying (mach, cd) points -- not merely "be accepted".
        let command = "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 100 -o json";
        let via_drg = wasm.run_command(command).unwrap();

        let via_csv = WasmBallistics::new();
        via_csv.load_drag_table(FLAT_CD_CSV.as_bytes()).unwrap();
        let via_csv_out = via_csv.run_command(command).unwrap();

        assert_eq!(
            via_drg, via_csv_out,
            ".drg-loaded and CSV-loaded identical-point tables must solve identically"
        );
    }

    #[wasm_bindgen_test]
    fn test_load_drag_table_junk_error_names_both_formats() {
        let wasm = WasmBallistics::new();
        let err = wasm
            .load_drag_table(b"this is not a drag table at all\nit is just some prose\n")
            .unwrap_err();
        let msg = err.as_string().unwrap_or_default();
        assert!(
            msg.to_lowercase().contains("csv"),
            "combined error should name the csv format, got: {msg}"
        );
        assert!(
            msg.to_lowercase().contains("drg"),
            "combined error should name the .drg format, got: {msg}"
        );
        assert!(!wasm.has_drag_table());
    }

    /// MBA-1294(a): the auto-zero "Rifle zeroed at ..." banner is a table-only human header;
    /// JSON output must stay pure machine output (a text prefix makes it unparseable).
    #[wasm_bindgen_test]
    fn test_json_output_pure_with_auto_zero() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --auto-zero 200 -o json")
            .unwrap();
        assert!(
            !result.contains("Rifle zeroed at"),
            "auto-zero banner leaked into JSON output: {result}"
        );
        // The prefix, if present, would break this parse.
        let _: serde_json::Value =
            serde_json::from_str(&result).expect("auto-zero JSON output must be pure JSON");
    }

    /// MBA-1294(a): same purity guarantee for CSV output.
    #[wasm_bindgen_test]
    fn test_csv_output_pure_with_auto_zero() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --auto-zero 200 -o csv")
            .unwrap();
        assert!(
            !result.contains("Rifle zeroed at"),
            "auto-zero banner leaked into CSV output: {result}"
        );
    }

    /// MBA-1294(c): --print-bc-segments appends the BC ladder in TABLE view only; it must
    /// never contaminate a JSON payload.
    #[wasm_bindgen_test]
    fn test_print_bc_segments_table_only() {
        let wasm = WasmBallistics::new();
        let table = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --print-bc-segments")
            .unwrap();
        // A plain scalar BC yields the "none active" note; either way the report block is present.
        assert!(
            table.contains("BC segments") || table.contains("BC Segments"),
            "bc-segments report missing from table view: {table}"
        );
        let json = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --print-bc-segments -o json")
            .unwrap();
        assert!(
            !json.contains("BC Segments (active)") && !json.contains("No BC segments active"),
            "bc-segments report leaked into JSON output: {json}"
        );
        let _: serde_json::Value =
            serde_json::from_str(&json).expect("--print-bc-segments must not break JSON");
    }

    /// MBA-1294(c): the `lead` command rejects an unknown -o value instead of silently
    /// defaulting to a format; a valid value is accepted.
    #[wasm_bindgen_test]
    fn test_lead_invalid_output_rejected() {
        let wasm = WasmBallistics::new();
        assert!(
            wasm.run_command("lead -v 2700 -m 168 -d 0.308 --target-speed 10 --range 300 -o json")
                .is_ok(),
            "lead -o json should be accepted"
        );
        let err = wasm
            .run_command("lead -v 2700 -m 168 -d 0.308 --target-speed 10 --range 300 -o bogus")
            .expect_err("lead must reject an invalid -o value");
        assert!(
            format!("{err:?}").contains("Invalid --output"),
            "unexpected lead -o error: {err:?}"
        );
    }

    /// MBA-1339: the WASM CLI accepts --bore-height as an alias for --muzzle-height on the
    /// same inches/mm units (unifying with the native CLI), so the two must produce identical
    /// trajectories.
    #[wasm_bindgen_test]
    fn test_bore_height_aliases_muzzle_height() {
        let wasm = WasmBallistics::new();
        let via_muzzle = wasm
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --muzzle-height 40 --max-range 300",
            )
            .unwrap();
        let via_bore = wasm
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --bore-height 40 --max-range 300",
            )
            .unwrap();
        assert_eq!(
            via_muzzle, via_bore,
            "--bore-height must alias --muzzle-height in the WASM CLI"
        );
    }

    // --- MBA-737: `powder` command (parity with the native CLI handler) ---

    #[wasm_bindgen_test]
    fn test_auto_zero_converges_on_inclined_shot() {
        // Field report (PRS use): --auto-zero with a nonzero --shooting-angle failed
        // ("not bracketed") because the zero solve inherited the shot-day incline.
        // The zero is torn on a level range: same adjustment at any shooting angle.
        let cmd = |ang: u32| {
            WasmBallistics::new()
                .run_command(&format!(
                    "trajectory -v 2650 -b 0.19 -m 77 -d 0.224 --drag-model g7 \
                     --sight-height 2.48 --auto-zero 100 --shooting-angle {ang} --max-range 300"
                ))
                .unwrap()
        };
        let flat = cmd(0);
        let incline = cmd(5);
        assert!(!incline.contains("Error calculating zero"), "{}", incline);
        let zero_line =
            |s: &str| s.lines().find(|l| l.contains("Rifle zeroed")).map(str::to_string);
        assert_eq!(
            zero_line(&flat),
            zero_line(&incline),
            "same rifle zero regardless of shot incline"
        );
    }

    #[wasm_bindgen_test]
    fn test_auto_zero_ignores_shot_coriolis() {
        // Fix-half of MBA-1384: the zero is torn on a known range without Coriolis
        // (native CLI behavior — its zero literal defaults enable_coriolis=false).
        // Before the fix, auto-zero inherited the shot's Coriolis conditions and the
        // browser terminal disagreed with the CLI on the rifle zero.
        let cmd = |coriolis: &str| {
            WasmBallistics::new()
                .run_command(&format!(
                    "trajectory -v 2650 -b 0.19 -m 77 -d 0.224 --drag-model g7 \
                     --sight-height 2.48 --auto-zero 100 --max-range 1000 {coriolis}"
                ))
                .unwrap()
        };
        let plain = cmd("");
        let coriolis = cmd("--enable-coriolis --latitude 45 --shot-direction 90");
        let zero_line =
            |s: &str| s.lines().find(|l| l.contains("Rifle zeroed")).map(str::to_string);
        assert_eq!(
            zero_line(&plain),
            zero_line(&coriolis),
            "rifle zero must not change with shot-day Coriolis conditions"
        );
    }

    #[wasm_bindgen_test]
    fn test_trajectory_plot_appends_charts() {
        // MBA-1337 p3: --plot appends both charts after the table (table-only).
        // MBA-1394 adds two more panels (velocity, energy) after drop/drift, mirroring
        // the native CLI's four-panel --plot output.
        let out = WasmBallistics::new()
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 300 --plot")
            .unwrap();
        assert!(out.contains("Drop vs Range:"), "{}", out);
        assert!(out.contains("Lateral Drift vs Range:"), "{}", out);
        assert!(out.contains("Velocity vs Range:"), "{}", out);
        assert!(out.contains("Energy vs Range:"), "{}", out);
        assert!(out.contains("velocity (fps)"), "{}", out);
        assert!(out.contains("energy (ft-lb)"), "{}", out);
        let ascii = WasmBallistics::new()
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 300 --plot ascii",
            )
            .unwrap();
        assert!(ascii.contains("Drop vs Range:"), "{}", ascii);
        assert!(ascii.contains("Velocity vs Range:"), "{}", ascii);
        assert!(ascii.contains("Energy vs Range:"), "{}", ascii);
        // JSON stays pure machine output even with --plot.
        let json = WasmBallistics::new()
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 300 --plot -o json",
            )
            .unwrap();
        assert!(!json.contains("Drop vs Range:"), "{}", json);
        assert!(!json.contains("Velocity vs Range:"), "{}", json);
        assert!(!json.contains("Energy vs Range:"), "{}", json);
    }

    #[wasm_bindgen_test]
    fn test_trajectory_plot_velocity_energy_labels_follow_units_metric() {
        // MBA-1394: under --units metric, the new panels must switch labels to m/s
        // and J, matching the metric labels the rest of this command uses.
        let out = WasmBallistics::new()
            .run_command(
                "trajectory --units metric -v 823 -b 0.475 -m 168 -d 0.308 --max-range 300 --plot",
            )
            .unwrap();
        assert!(out.contains("Velocity vs Range:"), "{}", out);
        assert!(out.contains("Energy vs Range:"), "{}", out);
        assert!(out.contains("velocity (m/s)"), "{}", out);
        assert!(out.contains("energy (J)"), "{}", out);
        assert!(!out.contains("velocity (fps)"), "{}", out);
        assert!(!out.contains("energy (ft-lb)"), "{}", out);
    }

    #[wasm_bindgen_test]
    fn test_powder_linear_model_mba1296_pin() {
        // The MBA-1296 repro: 1.0 fps/degF, 40F day, 70F reference -> 2770, never 8400.
        let out = WasmBallistics::new()
            .run_command("powder -v 2800 --temperature 40")
            .unwrap();
        assert!(out.contains("2770.0 fps"), "{}", out);
        assert!(out.contains("(-30.0)"), "{}", out);
    }

    #[wasm_bindgen_test]
    fn test_powder_curve_interpolates_at_powder_temp() {
        let out = WasmBallistics::new()
            .run_command("powder --powder-temp-curve 40:2620,70:2700,100:2760 --powder-temp 55")
            .unwrap();
        assert!(out.contains("2660.0 fps"), "{}", out);
        assert!(out.contains("measured curve, 3 points"), "{}", out);
    }

    #[wasm_bindgen_test]
    fn test_powder_sweep_csv_includes_end_row() {
        let out = WasmBallistics::new()
            .run_command("powder -v 2700 --sweep 20:110:30 -o csv")
            .unwrap();
        assert!(out.contains("2650.0,-50.0"), "{}", out);
        // The END row must not be dropped to float drift in the sweep expansion.
        assert!(out.contains("2740.0,40.0"), "{}", out);
    }

    #[wasm_bindgen_test]
    fn test_powder_json_resolved_object() {
        let out = WasmBallistics::new()
            .run_command("powder -v 2800 --temperature 40 -o json")
            .unwrap();
        assert!(out.contains("\"velocity\": 2770.0"), "{}", out);
        assert!(out.contains("\"shift\": -30.0"), "{}", out);
    }

    #[wasm_bindgen_test]
    fn test_powder_metric_resolves_in_celsius() {
        // 0.54864 m/s/degC x (5 - 21) degC = -8.778 m/s -> 814.2
        let out = WasmBallistics::new()
            .run_command("--units metric powder -v 823 --temperature 5 --powder-temp 21")
            .unwrap();
        assert!(out.contains("814.2 m/s"), "{}", out);
    }

    #[wasm_bindgen_test]
    fn test_powder_missing_inputs_error() {
        assert!(WasmBallistics::new()
            .run_command("powder --temperature 40")
            .is_err());
        assert!(WasmBallistics::new().run_command("powder -v 2800").is_err());
    }

    #[wasm_bindgen_test]
    fn test_powder_energy_with_mass() {
        // 168 gr at 2700 fps -> ~2719 ft-lb
        let out = WasmBallistics::new()
            .run_command("powder -v 2700 --temperature 70 -m 168")
            .unwrap();
        assert!(out.contains("2719 ft\u{b7}lb"), "{}", out);
    }

    // -----------------------------------------------------------------------------
    // MBA-1372: `recoil` (SAAMI free-recoil momentum balance) and `power-factor`
    // (USPSA/IDPA/SASS rulebook pass/fail). Numeric expectations here are the SAME
    // numbers asserted for the native CLI (see the doctests/tests in main.rs's
    // handle_recoil/handle_power_factor coverage and CLI_USAGE.md's worked examples) --
    // this is the WASM-parity cross-check: both front ends call into the identical
    // `ballistics_engine::recoil`/`ballistics_engine::power_factor` core.
    // -----------------------------------------------------------------------------

    #[wasm_bindgen_test]
    fn test_recoil_saami_rifle_default_matches_native_worked_example() {
        // .308 Win: 168 gr bullet, 43 gr charge, 2700 fps, 8.5 lb rifle, default
        // --firearm-type rifle (SAAMI f=1.75). Same numbers as the native CLI's
        // `recoil -b 168 -c 43 -v 2700 -f 8.5` smoke test.
        let out = WasmBallistics::new()
            .run_command("recoil -b 168 -c 43 -v 2700 -f 8.5")
            .unwrap();
        assert!(out.contains("saami-rifle"), "{}", out);
        assert!(out.contains("4725.0 fps"), "{}", out); // gas velocity = 1.75 * 2700
        assert!(out.contains("11.04 fps"), "{}", out); // recoil velocity
        assert!(out.contains("16.09 ft-lb"), "{}", out); // recoil energy
        assert!(out.contains("2.916 lb-s"), "{}", out); // recoil impulse
    }

    #[wasm_bindgen_test]
    fn test_recoil_metric_round_trip_matches_native() {
        // Same load as above, expressed in metric display units.
        let out = WasmBallistics::new()
            .run_command(
                "--units metric recoil -b 10.89 -c 2.79 -v 823.0 -f 3.86",
            )
            .unwrap();
        assert!(out.contains("1440.2 m/s") || out.contains("1440.3 m/s"), "{}", out);
        assert!(out.contains("3.36 m/s"), "{}", out);
        assert!(out.contains("21.8"), "{}", out); // ~21.82-21.83 J
        assert!(out.contains("12.98"), "{}", out); // ~12.97-12.98 N-s
    }

    #[wasm_bindgen_test]
    fn test_recoil_firearm_type_selects_saami_factor() {
        let pistol = WasmBallistics::new()
            .run_command("recoil -b 230 -c 6 -v 830 -f 2.5 --firearm-type pistol")
            .unwrap();
        assert!(pistol.contains("saami-pistol"), "{}", pistol);
        assert!(pistol.contains("1245.0 fps"), "{}", pistol); // 1.50 * 830

        let shotgun_long = WasmBallistics::new()
            .run_command("recoil -b 490 -c 30 -v 1300 -f 7.5 --firearm-type shotgun-long")
            .unwrap();
        assert!(shotgun_long.contains("saami-shotgun-long"), "{}", shotgun_long);
        assert!(shotgun_long.contains("1625.0 fps"), "{}", shotgun_long); // 1.25 * 1300
    }

    #[wasm_bindgen_test]
    fn test_recoil_gas_velocity_overrides_take_precedence() {
        // --gas-velocity-factor overrides --firearm-type.
        let factor = WasmBallistics::new()
            .run_command("recoil -b 168 -c 43 -v 2700 -f 8.5 --gas-velocity-factor 2.0")
            .unwrap();
        assert!(factor.contains("factor(2.000)"), "{}", factor);
        assert!(factor.contains("5400.0 fps"), "{}", factor); // 2.0 * 2700

        // --gas-velocity overrides both --gas-velocity-factor and --firearm-type, and
        // is independent of muzzle velocity.
        let fixed = WasmBallistics::new()
            .run_command(
                "recoil -b 168 -c 43 -v 2700 -f 8.5 --gas-velocity-factor 2.0 --gas-velocity 4700",
            )
            .unwrap();
        assert!(fixed.contains("fixed"), "{}", fixed);
        assert!(fixed.contains("4700.0 fps"), "{}", fixed);
    }

    #[wasm_bindgen_test]
    fn test_recoil_zero_charge_weight_is_bullet_only() {
        let out = WasmBallistics::new()
            .run_command("recoil -b 168 -c 0 -v 2700 -f 8.5")
            .unwrap();
        assert!(out.contains("Charge weight:       0.0 gr"), "{}", out);
    }

    #[wasm_bindgen_test]
    fn test_recoil_json_fields() {
        let out = WasmBallistics::new()
            .run_command("recoil -b 168 -c 43 -v 2700 -f 8.5 -o json")
            .unwrap();
        assert!(out.contains("\"command\": \"recoil\""), "{}", out);
        assert!(out.contains("\"gas_velocity_model\": \"saami-rifle\""), "{}", out);
    }

    #[wasm_bindgen_test]
    fn test_recoil_missing_required_flags_error() {
        assert!(WasmBallistics::new()
            .run_command("recoil -c 43 -v 2700 -f 8.5")
            .is_err());
        assert!(WasmBallistics::new()
            .run_command("recoil -b 168 -v 2700 -f 8.5")
            .is_err());
        assert!(WasmBallistics::new()
            .run_command("recoil -b 168 -c 43 -f 8.5")
            .is_err());
        assert!(WasmBallistics::new()
            .run_command("recoil -b 168 -c 43 -v 2700")
            .is_err());
    }

    #[wasm_bindgen_test]
    fn test_recoil_invalid_firearm_type_rejected() {
        assert!(WasmBallistics::new()
            .run_command("recoil -b 168 -c 43 -v 2700 -f 8.5 --firearm-type revolver")
            .is_err());
    }

    #[wasm_bindgen_test]
    fn test_recoil_nonpositive_firearm_weight_rejected() {
        // f64_range enforces > 0 (min 0.1), same as the native clap definition.
        assert!(WasmBallistics::new()
            .run_command("recoil -b 168 -c 43 -v 2700 -f 0")
            .is_err());
    }

    #[wasm_bindgen_test]
    fn test_power_factor_9mm_matches_native_worked_example() {
        // 147 gr @ 900 fps -> raw 132.30, scored 132: passes USPSA Minor/IDPA
        // SSP-ESP-CO, fails USPSA Major/IDPA PCC/Enhanced Revolver/CDP.
        let out = WasmBallistics::new()
            .run_command("power-factor -w 147 -v 900")
            .unwrap();
        assert!(out.contains("132.30"), "{}", out);
        assert!(out.contains("Power factor (scored): 132"), "{}", out);
        assert!(out.contains("USPSA     Minor                            125    PASS"), "{}", out);
        assert!(out.contains("USPSA     Major                            165    FAIL"), "{}", out);
        assert!(out.contains("IDPA      PCC                              135    FAIL"), "{}", out);
    }

    #[wasm_bindgen_test]
    fn test_power_factor_every_boundary_both_sides() {
        // USPSA Minor floor (125): weight 1000 makes velocity numerically equal PF.
        let at = WasmBallistics::new()
            .run_command("power-factor -w 1000 -v 125 --organization uspsa")
            .unwrap();
        assert!(at.contains("Minor                            125    PASS"), "{}", at);
        let below = WasmBallistics::new()
            .run_command("power-factor -w 1000 -v 124 --organization uspsa")
            .unwrap();
        assert!(below.contains("Minor                            125    FAIL"), "{}", below);

        // SASS floor (60) and minimum velocity (400 fps) boundary, isolated from each
        // other with a 200gr bullet (PF stays >= 60 across 399-400 fps).
        let sass_min_v_pass = WasmBallistics::new()
            .run_command("power-factor -w 200 -v 400 --organization sass")
            .unwrap();
        assert!(sass_min_v_pass.contains("PASS"), "{}", sass_min_v_pass);
        let sass_min_v_fail = WasmBallistics::new()
            .run_command("power-factor -w 200 -v 399 --organization sass")
            .unwrap();
        assert!(sass_min_v_fail.contains("FAIL"), "{}", sass_min_v_fail);
    }

    #[wasm_bindgen_test]
    fn test_power_factor_metric_converts_to_grains_fps_internally() {
        // 9.526 g == 147.0 gr (9.526 / 0.06479891); 900 fps == 274.32 m/s.
        let out = WasmBallistics::new()
            .run_command("--units metric power-factor -w 9.526 -v 274.32")
            .unwrap();
        assert!(out.contains("Power factor (raw):  132."), "{}", out);
    }

    #[wasm_bindgen_test]
    fn test_power_factor_organization_filter() {
        let out = WasmBallistics::new()
            .run_command("power-factor -w 147 -v 900 --organization sass")
            .unwrap();
        assert!(out.contains("SASS"), "{}", out);
        assert!(!out.contains("USPSA"), "{}", out);
        assert!(!out.contains("IDPA"), "{}", out);
    }

    #[wasm_bindgen_test]
    fn test_power_factor_unknown_organization_errors() {
        assert!(WasmBallistics::new()
            .run_command("power-factor -w 147 -v 900 --organization xyz")
            .is_err());
    }

    #[wasm_bindgen_test]
    fn test_power_factor_json_fields() {
        let out = WasmBallistics::new()
            .run_command("power-factor -w 147 -v 900 -o json")
            .unwrap();
        assert!(out.contains("\"power_factor_scored\": 132.0"), "{}", out);
        assert!(out.contains("\"organization\": \"SASS\""), "{}", out);
    }

    #[wasm_bindgen_test]
    fn test_power_factor_csv_has_no_embedded_commas_in_class_field() {
        // Every `class` value must be CSV-safe (no literal comma) since this handler
        // does no field quoting/escaping.
        let out = WasmBallistics::new()
            .run_command("power-factor -w 147 -v 900 -o csv")
            .unwrap();
        for line in out.lines().skip(1) {
            assert_eq!(
                line.split(',').count(),
                8,
                "row does not have exactly 8 CSV fields: {line}"
            );
        }
    }

    // -----------------------------------------------------------------------------
    // MBA-1343 Phase C: `true-velocity` (single- and multi-observation truing) and
    // `monte-carlo --wez`, both sharing the native compute cores
    // (ballistics_engine::truing / ballistics_engine::wez) with renderers that
    // replicate the native printers byte-for-byte.
    // -----------------------------------------------------------------------------

    /// Single-observation smoke: recovers an effective velocity via the shared
    /// binary-search core and prints the native table's headline strings. The
    /// pinned 2355.5 fps matches the native binary's `--offline` output for these
    /// exact args (deterministic solve, native parity pin).
    #[wasm_bindgen_test]
    fn true_velocity_single_observation_recovers_a_velocity() {
        let out = WasmBallistics::new()
            .run_command(
                "true-velocity --range 300 --measured-drop 1.8 --bc 0.475 -m 168 -d 0.308 \
                 --zero-distance 100 --chrono-velocity 2700",
            )
            .unwrap();
        assert!(out.contains("VELOCITY TRUING RESULTS"), "{}", out);
        assert!(out.contains("Effective Muzzle Velocity:"), "{}", out);
        assert!(out.contains("2355.5"), "{}", out);
        assert!(out.contains("Adjustment from Chrono:"), "{}", out);
    }

    /// Multi-observation native parity pin: the 300/600/900 mil set must run the
    /// joint MV+BC calibration and land on the native binary's fitted MV (2514.5
    /// fps) — the compute core is shared and deterministic, so any drift here is a
    /// real cross-surface divergence.
    #[wasm_bindgen_test]
    fn true_velocity_multi_observation_matches_native_joint_fit() {
        let out = WasmBallistics::new()
            .run_command(
                "true-velocity --range 300 --measured-drop 1.8 --observed 600:5.1 \
                 --observed 900:10.9 --bc 0.475 -m 168 -d 0.308 --zero-distance 100 \
                 --chrono-velocity 2700",
            )
            .unwrap();
        assert!(out.contains("Joint MV+BC fit"), "{}", out);
        assert!(out.contains("2514.5"), "{}", out);
    }

    /// MBA-1405 Task 2: the MV-calibration window table line is byte-parity with
    /// the native CLI's `display_multi_truing_result` (same 300/600/900 fixture,
    /// same pinned 656.7-729.7 yd window — see
    /// `tests/truing_multi_obs.rs::mv_calibration_window_line_appears_for_a_transonic_fixture`
    /// for the native-side pin this mirrors) — plus the "run plan-truing"
    /// cross-reference line, always present in table mode.
    #[wasm_bindgen_test]
    fn true_velocity_multi_observation_mv_window_line_matches_native() {
        let out = WasmBallistics::new()
            .run_command(
                "true-velocity --range 300 --measured-drop 1.8 --observed 600:5.1 \
                 --observed 900:10.9 --bc 0.475 -m 168 -d 0.308 --zero-distance 100 \
                 --chrono-velocity 2700",
            )
            .unwrap();
        assert!(
            out.contains("MV-calibration window: 656.7-729.7 yd (90-100% of the Mach 1.2 distance)"),
            "{out}"
        );
        assert!(
            out.contains("for optimal observation ranges run: ballistics plan-truing"),
            "{out}"
        );
    }

    /// The "no MV window" note, table only — mirrors the native CLI's fully
    /// supersonic branch (same fixture/pin as
    /// `tests/truing_multi_obs.rs::no_mv_window_note_for_a_fully_supersonic_fixture`).
    #[wasm_bindgen_test]
    fn true_velocity_multi_observation_no_window_note_matches_native() {
        let out = WasmBallistics::new()
            .run_command(
                "true-velocity --range 300 --measured-drop 0.3 --observed 600:0.65 \
                 --bc 2.0 -m 750 -d 0.510 --zero-distance 300 --chrono-velocity 3300",
            )
            .unwrap();
        assert!(
            out.contains(
                "note: no MV window: trajectory is supersonic through 3109.4 yd; MV is identifiable at any range"
            ),
            "{out}"
        );
        assert!(!out.contains("MV-calibration window:"), "{out}");
    }

    /// MBA-1405 Task 2 fix pass: the OTHER "no MV window" cause — a load that
    /// launches BELOW Mach 1.2 and so never crosses downward at all (as opposed to
    /// `true_velocity_multi_observation_no_window_note_matches_native` above, which
    /// is still supersonic through the whole window envelope). These two causes
    /// need opposite report text; this pins the never-supersonic branch byte-parity
    /// with native (same fixture/pin as
    /// `tests/truing_multi_obs.rs::no_mv_window_note_for_a_never_supersonic_fixture`).
    #[wasm_bindgen_test]
    fn true_velocity_multi_observation_never_supersonic_note_matches_native() {
        let out = WasmBallistics::new()
            .run_command(
                "true-velocity --range 100.0 --measured-drop 1.834371 --observed \
                 400.0:17.880184 --bc 0.3 -m 220 -d 0.308 --zero-distance 25 \
                 --chrono-velocity 1050",
            )
            .unwrap();
        assert!(
            out.contains(
                "note: no MV window: trajectory never reaches Mach 1.2; calibrate muzzle \
                 velocity with a chronograph, then collect DSF points"
            ),
            "{out}"
        );
        assert!(!out.contains("MV-calibration window:"), "{out}");
        assert!(!out.contains("supersonic through"), "{out}");
    }

    /// JSON purity for the never-supersonic branch too — still additive-only null
    /// fields, no note text of either flavor (matches the native CLI exactly).
    #[wasm_bindgen_test]
    fn true_velocity_multi_observation_never_supersonic_json_window_fields_null() {
        let out = WasmBallistics::new()
            .run_command(
                "true-velocity --range 100.0 --measured-drop 1.834371 --observed \
                 400.0:17.880184 --bc 0.3 -m 220 -d 0.308 --zero-distance 25 \
                 --chrono-velocity 1050 -o json",
            )
            .unwrap();
        assert!(!out.to_lowercase().contains("calibration window"), "{out}");
        assert!(!out.to_lowercase().contains("mach 1.2"), "{out}");
        let v: serde_json::Value = serde_json::from_str(&out).expect("json");
        assert!(v["mv_window_start_m"].is_null(), "{v}");
        assert!(v["mv_window_end_m"].is_null(), "{v}");
    }

    /// JSON additive-only purity: the window fields are present and numeric for the
    /// transonic fixture, and no note text ever leaks into JSON (matches the native
    /// CLI's JSON purity rule exactly).
    #[wasm_bindgen_test]
    fn true_velocity_multi_observation_json_has_additive_window_fields_only() {
        let out = WasmBallistics::new()
            .run_command(
                "true-velocity --range 300 --measured-drop 1.8 --observed 600:5.1 \
                 --observed 900:10.9 --bc 0.475 -m 168 -d 0.308 --zero-distance 100 \
                 --chrono-velocity 2700 -o json",
            )
            .unwrap();
        assert!(!out.to_lowercase().contains("calibration window"), "{out}");
        let v: serde_json::Value = serde_json::from_str(&out).expect("json");
        let lo = v["mv_window_start_m"].as_f64().expect("mv_window_start_m present");
        let hi = v["mv_window_end_m"].as_f64().expect("mv_window_end_m present");
        assert!((lo - 656.7 * 0.9144).abs() < 0.5, "lo={lo}");
        assert!((hi - 729.7 * 0.9144).abs() < 0.5, "hi={hi}");
    }

    /// JSON window fields are `null` (not absent, not a note string) when there is
    /// no window — mirrors the native CLI exactly.
    #[wasm_bindgen_test]
    fn true_velocity_multi_observation_json_window_fields_null_when_absent() {
        let out = WasmBallistics::new()
            .run_command(
                "true-velocity --range 300 --measured-drop 0.3 --observed 600:0.65 \
                 --bc 2.0 -m 750 -d 0.510 --zero-distance 300 --chrono-velocity 3300 -o json",
            )
            .unwrap();
        let v: serde_json::Value = serde_json::from_str(&out).expect("json");
        assert!(v["mv_window_start_m"].is_null(), "{v}");
        assert!(v["mv_window_end_m"].is_null(), "{v}");
    }

    /// WEZ summary parity: the seeded sweep (fixed per-step seed in the shared core)
    /// is deterministic, so the sims/step banner and the 200 yd P(hit) row must
    /// match the native binary's output for identical args verbatim.
    #[wasm_bindgen_test]
    fn wez_summary_reports_sims_per_step_and_p_hit_rows() {
        let out = WasmBallistics::new()
            .run_command(
                "monte-carlo -v 2700 -b 0.475 -m 168 -d 0.308 --wez --target-size 18x30 \
                 -n 300 --wez-start 200 --wez-end 500 --wez-step 100",
            )
            .unwrap();
        assert!(out.contains("WEZ sweep: 300 sims/step"), "{}", out);
        // A full P(hit) table row, byte-identical to the native summary printer.
        assert!(
            out.contains("│      200.0 │    49.0% │ other         │      0.0% │      0.0% │    100.0% │"),
            "{}",
            out
        );
    }

    /// `--wez` without `--target-size` must fail with the native dispatch's message.
    #[wasm_bindgen_test]
    fn wez_without_target_size_errors_like_native() {
        let err = WasmBallistics::new()
            .run_command("monte-carlo -v 2700 -b 0.475 -m 168 -d 0.308 --wez -n 10")
            .unwrap_err();
        assert!(
            err.as_string()
                .unwrap_or_default()
                .contains("--target-size is required with --wez (e.g. --target-size 18x30)"),
            "{err:?}"
        );
    }

    /// `--drag-model` must actually reach the WEZ sweep (MBA-1343 review: it was
    /// parsed and then dropped, so a G7 user silently got a G1 sweep). G1 and G7
    /// referenced to the same BC value give materially different trajectories, so
    /// the two summaries must differ; the G1 run must equal the default-model run.
    #[wasm_bindgen_test]
    fn wez_drag_model_reaches_the_sweep() {
        let base = "monte-carlo -v 2700 -b 0.475 -m 168 -d 0.308 --wez --target-size 18x30 \
                    -n 100 --wez-start 200 --wez-end 400 --wez-step 100";
        let default_out = WasmBallistics::new().run_command(base).unwrap();
        let g1_out = WasmBallistics::new()
            .run_command(&format!("{base} --drag-model g1"))
            .unwrap();
        let g7_out = WasmBallistics::new()
            .run_command(&format!("{base} --drag-model g7"))
            .unwrap();
        assert_eq!(default_out, g1_out, "explicit g1 must equal the default");
        assert_ne!(
            g1_out, g7_out,
            "a G7 sweep must not silently degrade to G1:\n{g7_out}"
        );
    }

    /// Bare (non-flag) tokens and dangling value-taking flags are hard errors on
    /// the truing/monte-carlo surfaces (MBA-1343 review: both were silently
    /// ignored, corrupting fits — `--observed 600:4.8 700:5.9` dropped the second
    /// point; a dangling `--observed` fell back to single-observation truing).
    #[wasm_bindgen_test]
    fn true_velocity_rejects_bare_tokens_and_dangling_flags() {
        let wasm = WasmBallistics::new();
        let err = wasm
            .run_command(
                "true-velocity --range 300 --measured-drop 1.8 --observed 600:4.8 700:5.9 \
                 --bc 0.475 -m 168 -d 0.308",
            )
            .unwrap_err();
        assert!(
            err.as_string()
                .unwrap_or_default()
                .contains("unexpected argument '700:5.9'"),
            "{err:?}"
        );
        let err = wasm
            .run_command("true-velocity --range 300 --measured-drop 1.8 --bc 0.475 --observed")
            .unwrap_err();
        assert!(
            err.as_string()
                .unwrap_or_default()
                .contains("a value is required for '--observed'"),
            "{err:?}"
        );
        let err = wasm
            .run_command("monte-carlo -v 2700 -b 0.475 -m 168 -d 0.308 300")
            .unwrap_err();
        assert!(
            err.as_string()
                .unwrap_or_default()
                .contains("unexpected argument '300'"),
            "{err:?}"
        );
        let err = wasm
            .run_command("monte-carlo -v 2700 -b 0.475 -m 168 -d 0.308 --velocity-std")
            .unwrap_err();
        assert!(
            err.as_string()
                .unwrap_or_default()
                .contains("a value is required for '--velocity-std'"),
            "{err:?}"
        );
    }

    /// Native clap's f64_range bounds for the true-velocity bullet parameters
    /// (-b 0.001..=2, -m 0.1..=2000, -d 0.01..=60) now also apply on this surface.
    #[wasm_bindgen_test]
    fn true_velocity_range_validates_bullet_parameters() {
        let wasm = WasmBallistics::new();
        for (args, needle) in [
            ("-b 5 -m 168 -d 0.308", "for '--bc'"),
            ("-b 0.475 -m 0.05 -d 0.308", "for '--mass'"),
            ("-b 0.475 -m 168 -d 99", "for '--diameter'"),
        ] {
            let err = wasm
                .run_command(&format!(
                    "true-velocity --range 300 --measured-drop 1.8 {args}"
                ))
                .unwrap_err();
            assert!(
                err.as_string().unwrap_or_default().contains(needle),
                "{args}: {err:?}"
            );
        }
    }

    /// MBA-1386: every `DragModel` family — including G2/G5/GI/GS and the new RA4 —
    /// now has a real reference table (src/drag.rs); none of them fall back to G1
    /// anymore, so the fix-half "using the G1 curve" note this test used to check
    /// for is gone (`g1_fallback_note`/`is_g1_fallback` were retired). This
    /// replaces that check with proof g5 actually solves with its own table: no
    /// fallback text anywhere in the output, and its result differs from a g1 run
    /// under identical inputs (a silent G1 alias would make them identical).
    #[wasm_bindgen_test]
    fn test_g5_solves_with_its_own_real_table() {
        let g5 = WasmBallistics::new()
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --drag-model g5 --max-range 300")
            .unwrap();
        assert!(!g5.contains("using the G1 curve"), "{}", g5);
        let g1 = WasmBallistics::new()
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --drag-model g1 --max-range 300")
            .unwrap();
        assert!(!g1.contains("using the G1 curve"), "{}", g1);
        assert_ne!(g5, g1, "g5 must apply its own drag table, not silently alias g1");
    }

    /// MBA-1386: RA4 (McCoy's British RA 1929 reference drag function) is the new
    /// ninth `DragModel` family. The WASM terminal's `--drag-model` parsing already
    /// goes through the shared `DragModel::from_str`, so `ra4` must be accepted and
    /// solve like any other family, with no fallback text in the output.
    #[wasm_bindgen_test]
    fn test_ra4_drag_model_accepted_and_solves() {
        let out = WasmBallistics::new()
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --drag-model ra4 --max-range 300")
            .unwrap();
        assert!(!out.contains("using the G1 curve"), "{}", out);
        assert!(out.contains("Range"), "expected a trajectory table, got: {}", out);
    }

    // MBA-1386: `test_g1_fallback_note_absent_from_structured_outputs` (formerly here)
    // guarded the table-only gating of the retired `g1_fallback_note` machinery across
    // trajectory/lead/monte-carlo --wez's JSON/CSV/table branches. That machinery no
    // longer exists (every family has a real table, so there is nothing left to gate),
    // so the test is deleted rather than adapted — there is no fallback note in any
    // output mode to assert absent.

    /// Serialize a minimal single-cell BC5D table into the v2 `.bin` byte layout so the
    /// coercion-warning tests below can exercise `loadBc5dTable` without an external file.
    /// Mirrors `bc_table_5d.rs`'s own (private-to-its-module) `serialize_test_table` test
    /// helper — duplicated here rather than shared since that one lives inside a `mod tests`
    /// this file cannot reach.
    fn synthetic_bc5d_bytes(correction: f32) -> Vec<u8> {
        synthetic_bc5d_bytes_for_caliber(correction, 0.308)
    }

    /// Same fixture with a chosen HEADER caliber, for the caliber-identity guard test.
    fn synthetic_bc5d_bytes_for_caliber(correction: f32, caliber: f32) -> Vec<u8> {
        let weight_bins = [168.0f32];
        let bc_bins = [0.4f32];
        let muzzle_vel_bins = [2500.0f32];
        let current_vel_bins = [2000.0f32];
        let data = [correction];

        let mut out = Vec::new();
        out.extend_from_slice(b"BC5D");
        out.extend_from_slice(&2u32.to_le_bytes()); // version
        out.extend_from_slice(&caliber.to_le_bytes());
        out.extend_from_slice(&0u32.to_le_bytes()); // flags
        out.extend_from_slice(&0u32.to_le_bytes()); // padding
        out.extend_from_slice(&(weight_bins.len() as u32).to_le_bytes());
        out.extend_from_slice(&(bc_bins.len() as u32).to_le_bytes());
        out.extend_from_slice(&(muzzle_vel_bins.len() as u32).to_le_bytes());
        out.extend_from_slice(&(current_vel_bins.len() as u32).to_le_bytes());
        out.extend_from_slice(&1u32.to_le_bytes()); // num_drag_types
        out.extend_from_slice(&0u64.to_le_bytes()); // timestamp

        let mut checksum_data = Vec::new();
        for v in weight_bins
            .iter()
            .chain(&bc_bins)
            .chain(&muzzle_vel_bins)
            .chain(&current_vel_bins)
            .chain(&data)
        {
            checksum_data.extend_from_slice(&v.to_le_bytes());
        }
        out.extend_from_slice(&crate::bc_table_5d::crc32_ieee(&checksum_data).to_le_bytes());

        let mut api = [0u8; 16];
        api[..4].copy_from_slice(b"test");
        out.extend_from_slice(&api);
        out.extend_from_slice(&[0u8; 12]); // reserved

        for v in weight_bins
            .iter()
            .chain(&bc_bins)
            .chain(&muzzle_vel_bins)
            .chain(&current_vel_bins)
            .chain(&data)
        {
            out.extend_from_slice(&v.to_le_bytes());
        }
        out
    }

    /// `loadBc5dTable` takes raw BYTES, so nothing before the solve ties the loaded table
    /// to the shot: a `.224` table applied to a `.308` bullet still clamped to its edge
    /// bins and silently biased every row. It must be refused, naming both calibers — and
    /// the same bytes must still be accepted for the caliber they actually describe.
    #[wasm_bindgen_test]
    fn bc5d_table_for_another_caliber_is_refused() {
        let wasm = WasmBallistics::new();
        wasm.load_bc5d_table(&synthetic_bc5d_bytes_for_caliber(0.9, 0.224))
            .unwrap();
        let error = wasm
            .run_command(
                "trajectory -v 2500 -b 0.4 -m 168 -d 0.308 --drag-model g1 --use-bc-segments \
                 --max-range 300",
            )
            .expect_err("a .224 table must not correct a .308 shot");
        let message = error.as_string().unwrap_or_default();
        assert!(
            message.contains("table is for 0.224, shot is 0.308"),
            "{message}"
        );

        // The same table is fine for a .224 shot.
        assert!(wasm
            .run_command(
                "trajectory -v 2500 -b 0.4 -m 55 -d 0.224 --drag-model g1 --use-bc-segments \
                 --max-range 300",
            )
            .is_ok());
    }

    /// 0.28.1 sweep: the BC5D G1/G7-coercion warning (MBA-1386, bcdd213) was only carried
    /// by the `Ok` arm of `match solver.solve()` in `handle_trajectory_command`; a solve
    /// that succeeds must show it in the table (and never in JSON/CSV).
    #[wasm_bindgen_test]
    fn bc5d_coercion_warning_survives_a_successful_solve() {
        let wasm = WasmBallistics::new();
        wasm.load_bc5d_table(&synthetic_bc5d_bytes(0.9)).unwrap();
        let table = wasm
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --drag-model g5 --use-bc-segments \
                 --max-range 300",
            )
            .unwrap();
        assert!(
            table.contains("BC5D correction tables support G1/G7 only"),
            "{table}"
        );
        assert!(table.contains("treating drag model 'G5' as G1"), "{table}");

        let json = wasm
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --drag-model g5 --use-bc-segments \
                 --max-range 300 -o json",
            )
            .unwrap();
        assert!(!json.contains("BC5D correction tables"), "{json}");
    }

    /// The same warning must survive the `Err` arm too: a solve that fails (here, an
    /// invalid `--bc` rejected by `validate_for_solve`) still computed the coercion warning
    /// up front, but the `Err` arm was dropping it on the floor instead of carrying it like
    /// the `Ok` arm does. Table-only, same as every warning in this formatter.
    #[wasm_bindgen_test]
    fn bc5d_coercion_warning_survives_a_failed_solve() {
        let wasm = WasmBallistics::new();
        wasm.load_bc5d_table(&synthetic_bc5d_bytes(0.9)).unwrap();
        let table = wasm
            .run_command(
                "trajectory -v 2700 -b -1 -m 168 -d 0.308 --drag-model g5 --use-bc-segments \
                 --max-range 300",
            )
            .unwrap();
        assert!(
            table.contains("BC5D correction tables support G1/G7 only"),
            "{table}"
        );
        assert!(table.contains("Error:"), "{table}");

        let json = wasm
            .run_command(
                "trajectory -v 2700 -b -1 -m 168 -d 0.308 --drag-model g5 --use-bc-segments \
                 --max-range 300 -o json",
            )
            .unwrap();
        assert!(!json.contains("BC5D correction tables"), "{json}");
        assert!(json.contains("Error:"), "{json}");
    }

    // -----------------------------------------------------------------------------
    // MBA-1411: `--dsf-point <MACH:DSF>` on `trajectory` — per-call parity with the
    // native CLI's saved-profile-carried DSF (drop-scale-factor) table (MBA-1357). WASM
    // has no profile storage, so the table is a repeatable per-call argument instead
    // (mirrors the MBA-1343 per-call pattern). Parsing/format errors are this command's
    // own; bounds/cap validation (0 < mach < 1.2, 0.5 < dsf < 2.0, <= 6 points) is NOT
    // reimplemented — `DsfTable::from_points`'s error text IS the user-facing error,
    // identical to native.
    // -----------------------------------------------------------------------------

    /// A malformed `MACH:DSF` token (wrong field count) is rejected with a message naming
    /// the bad token, before any solve.
    #[wasm_bindgen_test]
    fn dsf_point_bad_token_errors() {
        let err = WasmBallistics::new()
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --dsf-point abc")
            .unwrap_err();
        let msg = err.as_string().unwrap_or_default();
        assert!(
            msg.contains("--dsf-point expects MACH:DSF") && msg.contains("'abc'"),
            "{msg}"
        );
    }

    /// A non-numeric MACH or DSF field is rejected naming which field and the offending
    /// token, before any solve.
    #[wasm_bindgen_test]
    fn dsf_point_non_numeric_fields_error() {
        let wasm = WasmBallistics::new();
        let err = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --dsf-point x:1.0")
            .unwrap_err();
        let msg = err.as_string().unwrap_or_default();
        assert!(
            msg.contains("--dsf-point: invalid MACH in 'x:1.0'"),
            "{msg}"
        );

        let err = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --dsf-point 0.5:x")
            .unwrap_err();
        let msg = err.as_string().unwrap_or_default();
        assert!(
            msg.contains("--dsf-point: invalid DSF in '0.5:x'"),
            "{msg}"
        );
    }

    /// A Mach at/above the DSF ceiling (1.2 — MV-truing territory, not DSF's) surfaces
    /// `DsfTable::from_points`'s own error text verbatim, not a WASM-specific rewording.
    #[wasm_bindgen_test]
    fn dsf_point_out_of_range_mach_uses_engine_error_text() {
        let err = WasmBallistics::new()
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --dsf-point 1.5:1.1")
            .unwrap_err();
        let msg = err.as_string().unwrap_or_default();
        assert!(
            msg.contains("DSF point Mach 1.5 is out of range") && msg.contains("0 < mach < 1.2"),
            "{msg}"
        );
    }

    /// A DSF ratio outside `(0.5, 2.0)` likewise surfaces the engine's own error text.
    #[wasm_bindgen_test]
    fn dsf_point_out_of_range_dsf_uses_engine_error_text() {
        let err = WasmBallistics::new()
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --dsf-point 0.5:3.0")
            .unwrap_err();
        let msg = err.as_string().unwrap_or_default();
        assert!(
            msg.contains("DSF value 3 is out of range") && msg.contains("0.5 < dsf < 2"),
            "{msg}"
        );
    }

    /// More than 6 distinct `--dsf-point` occurrences hits `DsfTable::from_points`'s cap,
    /// naming the 6-point limit — the WASM parser does not reimplement this check.
    #[wasm_bindgen_test]
    fn dsf_point_more_than_six_uses_engine_cap_error() {
        let err = WasmBallistics::new()
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 \
                 --dsf-point 0.1:1.0 --dsf-point 0.2:1.0 --dsf-point 0.3:1.0 \
                 --dsf-point 0.4:1.0 --dsf-point 0.5:1.0 --dsf-point 0.6:1.0 \
                 --dsf-point 0.7:1.0",
            )
            .unwrap_err();
        let msg = err.as_string().unwrap_or_default();
        assert!(
            msg.contains("DSF table supports at most 6 points") && msg.contains("got 7"),
            "{msg}"
        );
    }

    /// Drop-only invariant through the WASM command path (mirrors
    /// `tests/dsf_workflow.rs`'s end-to-end beat): an active DSF table changes drop and
    /// nothing else — velocity/energy/time/range/drift stay byte-identical, and the point
    /// objects' key sets are unchanged (no new field, matching the native purity rule).
    /// `--max-range 1500 --ignore-ground-impact` gets a .308/168gr/2700fps/0.475 G1 shot
    /// well down into the transonic/subsonic band (crosses Mach 1.2 around 860 yd, Mach
    /// 1.0 around 1090 yd — verified against the native CLI's own JSON at the same inputs)
    /// so the DSF point at Mach 0.9 has real points to correct.
    #[wasm_bindgen_test]
    fn dsf_point_applies_drop_only_correction() {
        let wasm = WasmBallistics::new();
        let base_cmd = "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 1500 \
             --ignore-ground-impact -o json";
        let dsf_cmd = "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 1500 \
             --ignore-ground-impact --dsf-point 0.9:1.3 -o json";

        let baseline = wasm.run_command(base_cmd).unwrap();
        let corrected = wasm.run_command(dsf_cmd).unwrap();
        assert_ne!(
            baseline, corrected,
            "an active DSF table must change trajectory JSON output"
        );

        let base_json: serde_json::Value = serde_json::from_str(&baseline).unwrap();
        let corr_json: serde_json::Value = serde_json::from_str(&corrected).unwrap();

        // Drop-only: the summary block (max_range/max_height/time_of_flight/impact
        // velocity) never reaches DSF at all.
        assert_eq!(
            base_json["summary"], corr_json["summary"],
            "summary must stay identical"
        );
        assert_eq!(
            base_json["legend"], corr_json["legend"],
            "legend must stay identical"
        );

        let base_points = base_json["trajectory"].as_array().unwrap();
        let corr_points = corr_json["trajectory"].as_array().unwrap();
        assert!(!base_points.is_empty(), "expected a non-empty trajectory");
        assert_eq!(base_points.len(), corr_points.len());
        let mut any_drop_changed = false;
        for (bp, cp) in base_points.iter().zip(corr_points.iter()) {
            assert_eq!(
                bp["range_yards"], cp["range_yards"],
                "range must stay identical"
            );
            assert_eq!(
                bp["drift_inches"], cp["drift_inches"],
                "drift must stay identical"
            );
            assert_eq!(
                bp["velocity_fps"], cp["velocity_fps"],
                "velocity must stay identical"
            );
            assert_eq!(
                bp["energy_ftlb"], cp["energy_ftlb"],
                "energy must stay identical"
            );
            assert_eq!(
                bp["time_seconds"], cp["time_seconds"],
                "time must stay identical"
            );
            // No new field: DSF is a pure value correction, not an additive one.
            let bp_keys: std::collections::BTreeSet<_> = bp.as_object().unwrap().keys().collect();
            let cp_keys: std::collections::BTreeSet<_> = cp.as_object().unwrap().keys().collect();
            assert_eq!(bp_keys, cp_keys, "point key set must stay identical");

            let bd = bp["drop_inches"].as_f64().unwrap();
            let cd = cp["drop_inches"].as_f64().unwrap();
            if (bd - cd).abs() > 1e-9 {
                any_drop_changed = true;
            }
        }
        assert!(
            any_drop_changed,
            "at least one point's drop must change under an active DSF table"
        );
    }

    /// The "DSF table active (...)" note appears in table output only, with the same
    /// wording native main.rs prints for a saved profile's table — and never leaks into
    /// JSON, which also keeps an unchanged key set (additive-field purity, same rule the
    /// cd_scale/bc5d warnings above already follow).
    #[wasm_bindgen_test]
    fn dsf_point_note_is_table_output_only() {
        let wasm = WasmBallistics::new();
        let table = wasm
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 1500 \
                 --ignore-ground-impact --dsf-point 0.9:1.3",
            )
            .unwrap();
        assert!(
            table.contains("DSF table active (1 points, Mach 0.90-0.90)"),
            "{table}"
        );

        let json_without = wasm
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 1500 \
                 --ignore-ground-impact -o json",
            )
            .unwrap();
        let json_with = wasm
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 1500 \
                 --ignore-ground-impact --dsf-point 0.9:1.3 -o json",
            )
            .unwrap();
        assert!(
            !json_with.contains("DSF table active"),
            "note text must never leak into JSON: {json_with}"
        );

        // Unchanged top-level key set between the two JSON payloads (only values differ).
        let without_json: serde_json::Value = serde_json::from_str(&json_without).unwrap();
        let with_json: serde_json::Value = serde_json::from_str(&json_with).unwrap();
        let without_keys: std::collections::BTreeSet<_> =
            without_json.as_object().unwrap().keys().collect();
        let with_keys: std::collections::BTreeSet<_> =
            with_json.as_object().unwrap().keys().collect();
        assert_eq!(without_keys, with_keys, "top-level JSON key set must stay identical");
    }

    // -----------------------------------------------------------------------------
    // MBA-1411 (carried from the MBA-1356 review, "WASM lead untrued-curve gap"):
    // `lead` already applies a loaded custom drag table unconditionally, but never wired
    // up `--cd-scale` alongside it, unlike trajectory/zero/monte-carlo — a table trued
    // via `--cd-scale` on those commands would silently revert to its untrued curve on
    // `lead`. This section pins the fix: same pairing requirement, same shift, same
    // table-only warning as the existing cd_scale tests above.
    // -----------------------------------------------------------------------------

    /// Without a loaded drag table, `lead --cd-scale` is rejected before any solve —
    /// same pairing requirement, same error text as trajectory/zero/monte-carlo.
    #[wasm_bindgen_test]
    fn lead_cd_scale_without_a_loaded_drag_table_errors() {
        let err = WasmBallistics::new()
            .run_command(
                "lead -v 2700 -m 168 -d 0.308 --target-speed 10 --range 300 --cd-scale 1.1",
            )
            .unwrap_err();
        let msg = err.as_string().unwrap_or_default();
        assert_eq!(msg, "--cd-scale requires --drag-table");
    }

    /// With a table loaded, `--cd-scale` is accepted and actually shifts `lead`'s output
    /// (more drag changes time of flight, which changes the lead) — the plumb-through
    /// smoke, mirroring `cd_scale_accepted_with_a_loaded_table_and_shifts_trajectory_output`.
    #[wasm_bindgen_test]
    fn lead_cd_scale_accepted_with_a_loaded_table_and_shifts_output() {
        let wasm = WasmBallistics::new();
        wasm.load_drag_table(FLAT_CD_CSV.as_bytes()).unwrap();

        let neutral = wasm
            .run_command(
                "lead -v 2700 -m 168 -d 0.308 --target-speed 10 --range 300 --cd-scale 1.0 \
                 -o json",
            )
            .unwrap();
        let scaled = wasm
            .run_command(
                "lead -v 2700 -m 168 -d 0.308 --target-speed 10 --range 300 --cd-scale 1.1 \
                 -o json",
            )
            .unwrap();
        assert_ne!(
            neutral, scaled,
            "--cd-scale 1.1 must change lead output relative to the neutral 1.0"
        );

        let omitted = wasm
            .run_command("lead -v 2700 -m 168 -d 0.308 --target-speed 10 --range 300 -o json")
            .unwrap();
        assert_eq!(
            omitted, neutral,
            "omitting --cd-scale must be byte-identical to an explicit --cd-scale 1.0"
        );
    }

    /// A `--cd-scale` far outside the typical truing range warns once, table-only — the
    /// warning must never contaminate `lead`'s JSON output.
    #[wasm_bindgen_test]
    fn lead_cd_scale_out_of_range_warns_table_only() {
        let wasm = WasmBallistics::new();
        wasm.load_drag_table(FLAT_CD_CSV.as_bytes()).unwrap();

        let table_output = wasm
            .run_command(
                "lead -v 2700 -m 168 -d 0.308 --target-speed 10 --range 300 --cd-scale 3.0",
            )
            .unwrap();
        assert!(
            table_output.contains("--cd-scale 3 is far outside the typical truing range (0.90-1.10)"),
            "{table_output}"
        );

        let json_output = wasm
            .run_command(
                "lead -v 2700 -m 168 -d 0.308 --target-speed 10 --range 300 --cd-scale 3.0 \
                 -o json",
            )
            .unwrap();
        assert!(
            !json_output.contains("far outside the typical truing range"),
            "the range warning must never contaminate JSON output: {json_output}"
        );
    }

    // ---- MBA-1397: --pressure-type <absolute|qnh> / --zero-pressure-type ------------------

    const QNH_TRAJECTORY_BASE: &str = "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 \
         --units metric --max-range 300 --ignore-ground-impact --altitude 1500 \
         --pressure 1030.0 -o json";

    #[wasm_bindgen_test]
    fn trajectory_omitted_pressure_type_is_identical_to_explicit_absolute() {
        let wasm = WasmBallistics::new();
        let omitted = wasm.run_command(QNH_TRAJECTORY_BASE).unwrap();
        let explicit = wasm
            .run_command(&format!("{QNH_TRAJECTORY_BASE} --pressure-type absolute"))
            .unwrap();
        assert_eq!(
            omitted, explicit,
            "omitted --pressure-type must be byte-identical to --pressure-type absolute"
        );
    }

    #[wasm_bindgen_test]
    fn trajectory_qnh_mode_yields_a_higher_impact_velocity_than_absolute() {
        let wasm = WasmBallistics::new();
        let absolute_out = wasm.run_command(QNH_TRAJECTORY_BASE).unwrap();
        let qnh_out = wasm
            .run_command(&format!("{QNH_TRAJECTORY_BASE} --pressure-type qnh"))
            .unwrap();
        let absolute: serde_json::Value = serde_json::from_str(&absolute_out).unwrap();
        let qnh: serde_json::Value = serde_json::from_str(&qnh_out).unwrap();
        let v_absolute = absolute["impact_velocity"].as_f64().unwrap();
        let v_qnh = qnh["impact_velocity"].as_f64().unwrap();
        assert!(
            v_qnh > v_absolute + 10.0,
            "QNH-reduced (lower) pressure must retain velocity better than treating the same \
             reading as absolute station pressure: absolute={v_absolute} qnh={v_qnh}"
        );
    }

    #[wasm_bindgen_test]
    fn zero_omitted_pressure_type_is_identical_to_explicit_absolute() {
        let wasm = WasmBallistics::new();
        let base = "zero -v 2700 -b 0.475 -m 168 -d 0.308 --units metric \
             --target-distance 300 --altitude 1500 --pressure 1030.0 -o json";
        let omitted = wasm.run_command(base).unwrap();
        let explicit = wasm
            .run_command(&format!("{base} --pressure-type absolute"))
            .unwrap();
        assert_eq!(
            omitted, explicit,
            "omitted --pressure-type must be byte-identical to --pressure-type absolute on zero"
        );
    }

    #[wasm_bindgen_test]
    fn trajectory_auto_zero_omitted_zero_pressure_type_is_identical_to_explicit_absolute() {
        let wasm = WasmBallistics::new();
        let base = "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --units metric \
             --max-range 300 --ignore-ground-impact --altitude 1500 --pressure 1013.25 \
             --auto-zero 300 --zero-pressure 1030.0 -o json";
        let omitted = wasm.run_command(base).unwrap();
        let explicit = wasm
            .run_command(&format!("{base} --zero-pressure-type absolute"))
            .unwrap();
        assert_eq!(
            omitted, explicit,
            "omitted --zero-pressure-type must be byte-identical to --zero-pressure-type absolute"
        );
    }

    // ---- MBA-1366: --density-altitude -----------------------------------------------------

    const DA_TRAJECTORY_BASE: &str = "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 \
         --units metric --max-range 500 --ignore-ground-impact -o json";

    #[wasm_bindgen_test]
    fn trajectory_without_density_altitude_is_unaffected_by_the_new_flag_existing() {
        // The mere existence of --density-altitude parsing must not change output for a run
        // that never references it -- omitting it twice must be byte-identical.
        let wasm = WasmBallistics::new();
        let a = wasm.run_command(DA_TRAJECTORY_BASE).unwrap();
        let b = wasm.run_command(DA_TRAJECTORY_BASE).unwrap();
        assert_eq!(a, b);
    }

    #[wasm_bindgen_test]
    fn trajectory_density_altitude_yields_lower_impact_velocity_than_sea_level() {
        let wasm = WasmBallistics::new();
        let sea_level = wasm.run_command(DA_TRAJECTORY_BASE).unwrap();
        let with_da = wasm
            .run_command(&format!("{DA_TRAJECTORY_BASE} --density-altitude 2000"))
            .unwrap();
        let sea_level: serde_json::Value = serde_json::from_str(&sea_level).unwrap();
        let with_da: serde_json::Value = serde_json::from_str(&with_da).unwrap();
        let v_sea = sea_level["impact_velocity"].as_f64().unwrap();
        let v_da = with_da["impact_velocity"].as_f64().unwrap();
        // Thinner air at altitude -> less drag -> higher retained velocity.
        assert!(
            v_da > v_sea,
            "density altitude must thin the air relative to sea level: sea={v_sea} da={v_da}"
        );
    }

    #[wasm_bindgen_test]
    fn trajectory_density_altitude_supersedes_altitude_and_pressure() {
        let wasm = WasmBallistics::new();
        let da_only = wasm
            .run_command(&format!("{DA_TRAJECTORY_BASE} --density-altitude 2000"))
            .unwrap();
        let da_with_overrides = wasm
            .run_command(&format!(
                "{DA_TRAJECTORY_BASE} --density-altitude 2000 --altitude 50 --pressure 1013.25"
            ))
            .unwrap();
        assert_eq!(
            da_only, da_with_overrides,
            "--density-altitude must supersede --altitude/--pressure entirely"
        );
    }

    #[wasm_bindgen_test]
    fn trajectory_density_altitude_supersedes_pressure_type_and_notes_it_table_only() {
        let wasm = WasmBallistics::new();
        let table_base = "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --units metric \
             --max-range 500 --ignore-ground-impact --density-altitude 2000 --pressure 1030.0 \
             --pressure-type qnh";
        let out = wasm.run_command(table_base).unwrap();
        assert!(
            out.contains("--density-altitude supersedes --pressure/--pressure-type"),
            "table output must carry the supersede note (WASM has no visible stderr): {out}"
        );

        // JSON output must stay pure machine output -- no note text.
        let json_out = wasm
            .run_command(&format!("{table_base} -o json"))
            .unwrap();
        assert!(
            !json_out.contains("supersedes"),
            "JSON output must never be contaminated with the human-readable note: {json_out}"
        );
    }

    #[wasm_bindgen_test]
    fn trajectory_explicit_temperature_overrides_isa_default_under_density_altitude() {
        let wasm = WasmBallistics::new();
        let default_temp = wasm
            .run_command(&format!("{DA_TRAJECTORY_BASE} --density-altitude 3000"))
            .unwrap();
        let explicit_temp = wasm
            .run_command(&format!(
                "{DA_TRAJECTORY_BASE} --density-altitude 3000 --temperature 35"
            ))
            .unwrap();
        assert_ne!(
            default_temp, explicit_temp,
            "an explicit --temperature must change the resolved atmosphere under \
             --density-altitude"
        );
    }

    // -----------------------------------------------------------------------------
    // MBA-1377: `--chrono-distance` (correct a downrange chronograph reading back to
    // muzzle velocity). Native parity with `tests/chrono_distance_correction.rs`: the
    // hand-computed solver math itself is covered natively
    // (`truing::chrono_correction_tests`, shared by both surfaces); this section
    // proves the WASM terminal's independent hand-rolled arg parser wires the flag
    // through identically -- same additive no-op guarantee, same unit conversion,
    // same error surfacing.
    // -----------------------------------------------------------------------------

    const CHRONO_DISTANCE_BASE: &str = "true-velocity --range 300 --measured-drop 1.8 -b 0.475 \
         --drag-model g7 -m 168 -d 0.308 --chrono-velocity 2680";

    #[wasm_bindgen_test]
    fn chrono_distance_absent_is_byte_identical_to_explicit_zero() {
        let wasm = WasmBallistics::new();
        let absent = wasm.run_command(CHRONO_DISTANCE_BASE).unwrap();
        let zero = wasm
            .run_command(&format!("{CHRONO_DISTANCE_BASE} --chrono-distance 0"))
            .unwrap();
        assert_eq!(
            absent, zero,
            "omitting --chrono-distance must be byte-identical to passing an explicit 0"
        );
    }

    #[wasm_bindgen_test]
    fn chrono_distance_nonzero_changes_the_displayed_adjustment() {
        let wasm = WasmBallistics::new();
        let baseline = wasm.run_command(CHRONO_DISTANCE_BASE).unwrap();
        let corrected = wasm
            .run_command(&format!("{CHRONO_DISTANCE_BASE} --chrono-distance 15"))
            .unwrap();
        assert_ne!(
            baseline, corrected,
            "a nonzero --chrono-distance must change the displayed comparison"
        );
    }

    #[wasm_bindgen_test]
    fn chrono_distance_requires_chrono_velocity() {
        let wasm = WasmBallistics::new();
        let err = wasm
            .run_command(
                "true-velocity --range 300 --measured-drop 1.8 -b 0.475 -m 168 -d 0.308 \
                 --chrono-distance 15",
            )
            .unwrap_err();
        assert!(
            err.as_string()
                .unwrap_or_default()
                .contains("--chrono-distance requires --chrono-velocity"),
            "{err:?}"
        );
    }

    #[wasm_bindgen_test]
    fn chrono_distance_out_of_band_is_rejected() {
        let wasm = WasmBallistics::new();
        for bad in ["-5", "0.05", "1000"] {
            let err = wasm
                .run_command(&format!("{CHRONO_DISTANCE_BASE} --chrono-distance {bad}"))
                .unwrap_err();
            assert!(
                err.as_string().unwrap_or_default().contains("out of range"),
                "distance {bad}: {err:?}"
            );
        }
    }
}
