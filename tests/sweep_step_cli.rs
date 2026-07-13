// MBA-1288: sweep subcommands must reject non-positive --step at parse time.
// Before the fix, `lead` and `wind-card` spun forever on --step 0 (and negative
// steps walked the sweep backward forever while solving a trajectory per
// iteration); `come-ups`/`range-table` only escaped via accidental downstream
// guards with inconsistent errors. All four now share a parse-time floor.
use std::process::Command;

const BIN: &str = env!("CARGO_BIN_EXE_ballistics");

/// Common ballistic args each subcommand needs to get past required-arg checks.
fn base_args(subcommand: &'static str) -> Vec<&'static str> {
    let mut v = vec![subcommand, "-v", "2700", "-m", "168", "-d", "0.308", "-b", "0.5"];
    match subcommand {
        "come-ups" | "wind-card" | "range-table" => v.extend(["--zero-distance", "100"]),
        "lead" => v.extend(["--target-speed", "3"]),
        _ => {}
    }
    v
}

fn run_with_step(subcommand: &'static str, step: &str) -> std::process::Output {
    Command::new(BIN)
        .args(base_args(subcommand))
        .arg(format!("--step={step}"))
        .output()
        .expect("spawn ballistics")
}

#[test]
fn zero_step_is_rejected_at_parse_time_by_all_sweep_subcommands() {
    for sub in ["lead", "wind-card", "come-ups", "range-table"] {
        let out = run_with_step(sub, "0");
        assert!(
            !out.status.success(),
            "{sub} --step 0 must fail (was an infinite loop before MBA-1288)"
        );
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(
            stderr.contains("not in range"),
            "{sub} --step 0 must be a clap range rejection, got: {stderr}"
        );
    }
}

#[test]
fn negative_step_is_rejected_at_parse_time_by_all_sweep_subcommands() {
    for sub in ["lead", "wind-card", "come-ups", "range-table"] {
        let out = run_with_step(sub, "-50");
        assert!(
            !out.status.success(),
            "{sub} --step=-50 must fail (walked the sweep backward forever before MBA-1288)"
        );
    }
}

#[test]
fn valid_step_still_accepted() {
    // Positive control: the floor must not over-reject a normal sweep.
    let out = run_with_step("lead", "200");
    assert!(
        out.status.success(),
        "lead --step=200 must still work: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(stdout.contains("Lead"), "expected a lead table, got: {stdout}");
}
