// MBA-1285: CLI-level coverage for --drag-table on the `trajectory` subcommand. Drives the real
// built binary (no assert_cmd in this repo's dev-dependencies) via the cargo-provided
// CARGO_BIN_EXE_<name> path: a well-formed deck must run to completion, and a malformed deck
// must fail cleanly (non-zero exit, no panic) rather than silently falling back to the G-model.
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
