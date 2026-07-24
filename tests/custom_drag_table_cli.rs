// MBA-1285: CLI-level coverage for --drag-table on the `trajectory` subcommand. Drives the real
// built binary (no assert_cmd in this repo's dev-dependencies) via the cargo-provided
// CARGO_BIN_EXE_<name> path: a well-formed deck must run to completion, and a malformed deck
// must fail cleanly (non-zero exit, no panic) rather than silently falling back to the G-model.
//
// MBA-1356 adds `--cd-scale` coverage at the bottom of this file: the pairing requirement with
// --drag-table, and the plumb-through smoke (a higher scale must measurably increase drag).
use std::io::Write;
use std::process::Command;

fn write_temp_csv(name: &str, contents: &[u8]) -> std::path::PathBuf {
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let path =
        std::env::temp_dir().join(format!("mba1285_{name}_{}_{nonce}.csv", std::process::id()));
    std::fs::File::create(&path)
        .unwrap()
        .write_all(contents)
        .unwrap();
    path
}

#[test]
fn cli_drag_table_runs_on_good_deck() {
    let bin = env!("CARGO_BIN_EXE_ballistics");
    let good = write_temp_csv("good", b"mach,cd\n0.5,0.23\n1.0,0.40\n2.0,0.30\n3.0,0.26\n");

    let out = Command::new(bin)
        .args([
            "trajectory",
            "--velocity",
            "850",
            "--bc",
            "0.5",
            "--mass",
            "168",
            "--diameter",
            "0.308",
            "--max-range",
            "400",
            "--drag-table",
        ])
        .arg(&good)
        .output()
        .unwrap();

    std::fs::remove_file(&good).ok();

    assert!(
        out.status.success(),
        "good deck should run: stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(
        stdout.contains("TRAJECTORY") || stdout.contains("Range"),
        "expected trajectory output, got:\n{stdout}"
    );
}

#[test]
fn cli_drag_table_rejects_malformed_deck() {
    let bin = env!("CARGO_BIN_EXE_ballistics");
    let bad = write_temp_csv("bad", b"0.5,0.23\n1.0,nope\n");

    let out = Command::new(bin)
        .args([
            "trajectory",
            "--velocity",
            "850",
            "--bc",
            "0.5",
            "--mass",
            "168",
            "--diameter",
            "0.308",
            "--max-range",
            "400",
            "--drag-table",
        ])
        .arg(&bad)
        .output()
        .unwrap();

    std::fs::remove_file(&bad).ok();

    assert!(
        !out.status.success(),
        "malformed deck should fail, not silently succeed"
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("Error"),
        "expected a clear error message, got stderr:\n{stderr}\nstdout:\n{}",
        String::from_utf8_lossy(&out.stdout)
    );
}

// ---- MBA-1356: --cd-scale -------------------------------------------------------------------

/// `--cd-scale` without `--drag-table` must be rejected BEFORE any solve, with the exact
/// pairing-error text named in the MBA-1356 spec (naming --bc-adjustment as the G1/G7
/// alternative), not a generic parse failure or a silent no-op.
#[test]
fn cli_cd_scale_requires_drag_table_errors() {
    let bin = env!("CARGO_BIN_EXE_ballistics");

    let out = Command::new(bin)
        .args([
            "trajectory",
            "--velocity",
            "850",
            "--bc",
            "0.5",
            "--mass",
            "168",
            "--diameter",
            "0.308",
            "--max-range",
            "400",
            "--cd-scale",
            "1.1",
        ])
        .output()
        .unwrap();

    assert!(
        !out.status.success(),
        "--cd-scale without --drag-table must fail, not silently run"
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains(
            "error: --cd-scale requires --drag-table (for G1/G7 use --bc-adjustment instead)"
        ),
        "expected the exact pairing-error text, got stderr:\n{stderr}"
    );
}

/// Same pairing requirement on `monte-carlo` and `zero` — `--cd-scale` is not trajectory-only.
/// Neither subcommand has a `--bc-adjustment` flag (only `trajectory` does), so unlike the
/// `trajectory` pairing error above, theirs must NOT suggest a flag the surface doesn't
/// accept (0.28.1 sweep).
#[test]
fn cli_cd_scale_requires_drag_table_errors_on_monte_carlo_and_zero() {
    let bin = env!("CARGO_BIN_EXE_ballistics");

    let mc_out = Command::new(bin)
        .args([
            "monte-carlo",
            "--velocity",
            "850",
            "--bc",
            "0.5",
            "--mass",
            "168",
            "--diameter",
            "0.308",
            "--num-sims",
            "10",
            "--cd-scale",
            "1.1",
        ])
        .output()
        .unwrap();
    assert!(!mc_out.status.success(), "monte-carlo: must fail");
    let mc_stderr = String::from_utf8_lossy(&mc_out.stderr);
    assert!(
        mc_stderr.contains("--cd-scale requires --drag-table"),
        "monte-carlo stderr:\n{mc_stderr}"
    );
    assert!(
        !mc_stderr.contains("--bc-adjustment"),
        "monte-carlo has no --bc-adjustment flag to suggest: {mc_stderr}"
    );

    let zero_out = Command::new(bin)
        .args([
            "zero",
            "--velocity",
            "850",
            "--bc",
            "0.5",
            "--mass",
            "168",
            "--diameter",
            "0.308",
            "--target-distance",
            "300",
            "--cd-scale",
            "1.1",
        ])
        .output()
        .unwrap();
    assert!(!zero_out.status.success(), "zero: must fail");
    let zero_stderr = String::from_utf8_lossy(&zero_out.stderr);
    assert!(
        zero_stderr.contains("--cd-scale requires --drag-table"),
        "zero stderr:\n{zero_stderr}"
    );
    assert!(
        !zero_stderr.contains("--bc-adjustment"),
        "zero has no --bc-adjustment flag to suggest: {zero_stderr}"
    );
}

/// Run `trajectory --drag-table <deck> --cd-scale <scale> -o json` and return the parsed
/// `impact_velocity` field (fps, imperial default).
fn trajectory_impact_velocity_with_scale(deck: &std::path::Path, cd_scale: &str) -> f64 {
    let bin = env!("CARGO_BIN_EXE_ballistics");
    let out = Command::new(bin)
        .args([
            "trajectory",
            "--velocity",
            "2700",
            "--bc",
            "0.5",
            "--mass",
            "168",
            "--diameter",
            "0.308",
            "--max-range",
            "500",
            "--drag-table",
        ])
        .arg(deck)
        .args(["--cd-scale", cd_scale, "-o", "json"])
        .output()
        .unwrap();

    assert!(
        out.status.success(),
        "trajectory --cd-scale {cd_scale} should run: stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    let json: serde_json::Value =
        serde_json::from_str(&stdout).unwrap_or_else(|e| panic!("invalid JSON ({e}):\n{stdout}"));
    json["impact_velocity"]
        .as_f64()
        .unwrap_or_else(|| panic!("no impact_velocity in JSON:\n{stdout}"))
}

/// Plumb-through smoke: `--cd-scale 1.1` must actually reach `BallisticInputs.cd_scale` and
/// increase drag over the whole flight, measurably lowering impact velocity relative to the
/// explicit neutral `1.0` on the identical deck/shot.
#[test]
fn cli_cd_scale_plumbs_through_and_increases_drag() {
    let deck = write_temp_csv(
        "cd_scale_plumb",
        b"mach,cd\n0.5,0.220\n0.8,0.230\n1.0,0.520\n1.2,0.480\n1.5,0.400\n2.0,0.330\n2.5,0.300\n",
    );

    let baseline = trajectory_impact_velocity_with_scale(&deck, "1.0");
    let scaled = trajectory_impact_velocity_with_scale(&deck, "1.1");

    std::fs::remove_file(&deck).ok();

    assert!(
        scaled < baseline,
        "--cd-scale 1.1 must increase drag -> lower impact velocity: baseline={baseline} scaled={scaled}"
    );
}

/// `--cd-scale 3.0` (far outside the documented typical truing range 0.90-1.10) must warn once
/// on stderr but still solve successfully — the engine's own gate is only finite && > 0.
#[test]
fn cli_cd_scale_out_of_range_warns_once_but_still_solves() {
    let deck = write_temp_csv(
        "cd_scale_warn",
        b"mach,cd\n0.5,0.23\n1.0,0.40\n2.0,0.30\n3.0,0.26\n",
    );

    let out = Command::new(env!("CARGO_BIN_EXE_ballistics"))
        .args([
            "trajectory",
            "--velocity",
            "850",
            "--bc",
            "0.5",
            "--mass",
            "168",
            "--diameter",
            "0.308",
            "--max-range",
            "400",
            "--drag-table",
        ])
        .arg(&deck)
        .args(["--cd-scale", "3.0"])
        .output()
        .unwrap();

    std::fs::remove_file(&deck).ok();

    assert!(
        out.status.success(),
        "an out-of-range (but finite, positive) cd_scale must still solve: stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    let occurrences = stderr.matches("far outside the typical truing range").count();
    assert_eq!(
        occurrences, 1,
        "expected exactly one out-of-range warning, got {occurrences} in stderr:\n{stderr}"
    );
    assert!(stderr.contains("--cd-scale 3"), "{stderr}");
}
