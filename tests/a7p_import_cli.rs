//! End-to-end test of `ballistics profile import` against a synthetic .a7p.
#![cfg(feature = "profile-import")]

use std::process::Command;

fn cli() -> Command {
    Command::new(env!("CARGO_BIN_EXE_ballistics"))
}

// Spec-derived encoder (same shape as the unit-test helpers).
fn enc_varint(mut v: u64, out: &mut Vec<u8>) {
    loop {
        let byte = (v & 0x7f) as u8;
        v >>= 7;
        if v == 0 {
            out.push(byte);
            return;
        }
        out.push(byte | 0x80);
    }
}
fn enc_i32(number: u32, value: i64, out: &mut Vec<u8>) {
    enc_varint(u64::from(number) << 3, out);
    enc_varint(value as u64, out);
}
fn enc_str(number: u32, s: &str, out: &mut Vec<u8>) {
    enc_varint((u64::from(number) << 3) | 2, out);
    enc_varint(s.len() as u64, out);
    out.extend_from_slice(s.as_bytes());
}
fn enc_bytes(number: u32, payload: &[u8], out: &mut Vec<u8>) {
    enc_varint((u64::from(number) << 3) | 2, out);
    enc_varint(payload.len() as u64, out);
    out.extend_from_slice(payload);
}

fn synthetic_a7p() -> Vec<u8> {
    let mut p = Vec::new();
    enc_str(1, "CLI 6.5CM", &mut p);
    enc_str(3, "140 ELD-M", &mut p);
    enc_i32(9, 50, &mut p); // 50 mm sight height
    enc_i32(10, 800, &mut p); // 8.00 in twist
    enc_i32(11, 8230, &mut p); // 823.0 m/s
    enc_i32(15, 15, &mut p);
    enc_i32(16, 10132, &mut p); // 1013.2 hPa
    enc_i32(17, 50, &mut p);
    enc_i32(20, 264, &mut p); // 0.264 in
    enc_i32(21, 1400, &mut p); // 140 gr
    enc_i32(22, 1360, &mut p); // 1.36 in
    enc_i32(24, 1, &mut p); // G7
    let mut packed = Vec::new();
    enc_varint(10_000, &mut packed);
    enc_bytes(26, &packed, &mut p);
    let mut row = Vec::new();
    enc_i32(1, 3260, &mut row); // G7 BC 0.326
    enc_i32(2, 8230, &mut row);
    enc_bytes(27, &row, &mut p);
    let mut payload = Vec::new();
    enc_bytes(1, &p, &mut payload);
    ballistics_engine::profile_import::wrap_payload(&payload)
}

#[test]
fn import_dry_run_saves_nothing_and_prints_report() {
    let home = tempfile_dir();
    let a7p_path = home.join("test.a7p");
    std::fs::write(&a7p_path, synthetic_a7p()).unwrap();

    let output = cli()
        .env("HOME", &home)
        .args(["profile", "import"])
        .arg(&a7p_path)
        .arg("--dry-run")
        .output()
        .unwrap();
    assert!(output.status.success(), "{output:?}");
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("Imported fields:"), "{stdout}");
    assert!(stdout.contains("823.0 m/s"), "{stdout}");
    assert!(stdout.contains("Dry run"), "{stdout}");
    assert!(
        !home.join(".ballistics").join("profiles").join("CLI 6.5CM.json").exists(),
        "dry run must not write"
    );
}

#[test]
fn import_then_show_roundtrip() {
    let home = tempfile_dir();
    let a7p_path = home.join("test.a7p");
    std::fs::write(&a7p_path, synthetic_a7p()).unwrap();

    let output = cli()
        .env("HOME", &home)
        .args(["profile", "import"])
        .arg(&a7p_path)
        .args(["--name", "imported-cm"])
        .output()
        .unwrap();
    assert!(output.status.success(), "{output:?}");
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("Saved profile"), "{stdout}");

    let show = cli()
        .env("HOME", &home)
        .args(["profile", "show", "imported-cm"])
        .output()
        .unwrap();
    assert!(show.status.success(), "{show:?}");
    let shown = String::from_utf8_lossy(&show.stdout);
    assert!(shown.contains("imported-cm"), "{shown}");
    assert!(shown.contains("G7"), "{shown}");
}

#[test]
fn corrupted_checksum_warns_by_default_and_fails_with_strict() {
    let home = tempfile_dir();
    let mut bytes = synthetic_a7p();
    bytes[0] = if bytes[0] == b'0' { b'1' } else { b'0' };
    let a7p_path = home.join("bad.a7p");
    std::fs::write(&a7p_path, &bytes).unwrap();

    let lenient = cli()
        .env("HOME", &home)
        .args(["profile", "import"])
        .arg(&a7p_path)
        .arg("--dry-run")
        .output()
        .unwrap();
    assert!(lenient.status.success(), "{lenient:?}");
    let stdout = String::from_utf8_lossy(&lenient.stdout);
    assert!(stdout.contains("checksum mismatch"), "{stdout}");

    let strict = cli()
        .env("HOME", &home)
        .args(["profile", "import"])
        .arg(&a7p_path)
        .args(["--dry-run", "--strict"])
        .output()
        .unwrap();
    assert!(!strict.status.success(), "strict must fail on checksum mismatch");
}

#[test]
fn garbage_file_is_a_clean_error() {
    let home = tempfile_dir();
    let path = home.join("junk.a7p");
    std::fs::write(&path, b"not a profile").unwrap();
    let output = cli()
        .env("HOME", &home)
        .args(["profile", "import"])
        .arg(&path)
        .output()
        .unwrap();
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("a7p"), "{stderr}");
}

#[test]
fn import_with_slash_name_is_sanitized_with_note() {
    let home = tempfile_dir();
    let a7p_path = home.join("test.a7p");
    std::fs::write(&a7p_path, synthetic_a7p()).unwrap();

    let output = cli()
        .env("HOME", &home)
        .args(["profile", "import"])
        .arg(&a7p_path)
        .args(["--name", "match/338"])
        .output()
        .unwrap();
    assert!(output.status.success(), "{output:?}");
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        stdout.contains("note: profile name sanitized to 'match_338' for the profile store"),
        "{stdout}"
    );

    let show = cli()
        .env("HOME", &home)
        .args(["profile", "show", "match_338"])
        .output()
        .unwrap();
    assert!(show.status.success(), "{show:?}");
}

/// Unique temp dir per test without external crates.
fn tempfile_dir() -> std::path::PathBuf {
    use std::sync::atomic::{AtomicU32, Ordering};
    static N: AtomicU32 = AtomicU32::new(0);
    let dir = std::env::temp_dir().join(format!(
        "a7p-import-test-{}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}
