// MBA-953: the BCCR loader previously read but never verified the stored CRC32, so corrupt-but-
// in-range files loaded silently. The checksum is a standard CRC-32 (IEEE/zlib) over the DATA
// section only (the f32 cells) — confirmed against the in-tree reference data/bc_corrections_
// transonic.bin, whose stored checksum is 0xb0a9edc9. These tests guard that the reference passes
// and that a single corrupted data byte is rejected.
use ballistics_engine::bc_table::{BcCorrectionTable, BcTableError};
use std::path::PathBuf;

fn ref_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("data/bc_corrections_transonic.bin")
}

#[test]
fn bccr_reference_passes_crc() {
    let table =
        BcCorrectionTable::load(ref_path()).expect("reference BCCR should load and pass its CRC");
    assert!(table.total_cells() > 0, "reference table should have cells");
}

#[test]
fn bccr_corrupt_data_byte_rejected() {
    let mut bytes = std::fs::read(ref_path()).expect("read reference file");
    // Flip a byte deep in the DATA section (data starts after the 60-byte header + 4 bin arrays
    // = offset 428; the file is ~1.5 MB, so len-100 is well inside the data).
    let i = bytes.len() - 100;
    bytes[i] ^= 0xFF;

    let tmp = std::env::temp_dir().join("bccr_corrupt_953.bin");
    std::fs::write(&tmp, &bytes).expect("write corrupt copy");
    let result = BcCorrectionTable::load(&tmp);
    let _ = std::fs::remove_file(&tmp);

    assert!(
        matches!(result, Err(BcTableError::ChecksumMismatch)),
        "a corrupt data byte must be rejected with ChecksumMismatch (load now verifies the CRC)"
    );
}
