// MBA-953: the BCCR loader previously read but never verified the stored CRC32, so corrupt-but-
// in-range files loaded silently. The checksum is a standard CRC-32 (IEEE/zlib) over the DATA
// section only (the f32 cells). Build the smallest valid table in the test so release-source
// archives do not depend on the ignored 1.5 MB benchmark/reference table.
use ballistics_engine::bc_table::{BcCorrectionTable, BcTableError};
use std::path::PathBuf;

fn push_f32(output: &mut Vec<u8>, value: f32) {
    output.extend_from_slice(&value.to_le_bytes());
}

fn checksummed_fixture() -> Vec<u8> {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(b"BCCR");
    bytes.extend_from_slice(&1_u32.to_le_bytes());
    bytes.extend_from_slice(&0_u32.to_le_bytes()); // flags
    for dimension in [1_u32; 5] {
        bytes.extend_from_slice(&dimension.to_le_bytes());
    }
    bytes.extend_from_slice(&0_u64.to_le_bytes()); // timestamp

    // IEEE/zlib CRC-32 of the little-endian bytes for the sole 0.8f32 data cell:
    // cd cc 4c 3f -> d77cabd7. Keep this independent constant so the test does not reproduce
    // the loader's checksum implementation to generate its own expected value.
    bytes.extend_from_slice(&0xd77c_abd7_u32.to_le_bytes());
    bytes.extend_from_slice(&[0_u8; 16]);

    push_f32(&mut bytes, 0.4); // BC bin
    push_f32(&mut bytes, 168.0); // mass bin (grains)
    push_f32(&mut bytes, 1.0); // length bin (inches)
    push_f32(&mut bytes, 2_800.0); // velocity bin (fps)
    push_f32(&mut bytes, 0.8); // correction data cell covered by the checksum
    assert_eq!(bytes.len(), 80);
    bytes
}

fn temp_path(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!("bccr_{label}_{}.bin", std::process::id()))
}

fn write_fixture(label: &str, bytes: &[u8]) -> PathBuf {
    let path = temp_path(label);
    std::fs::write(&path, bytes).expect("write BCCR fixture");
    path
}

#[test]
fn bccr_checksummed_fixture_passes_crc() {
    let path = write_fixture("valid_953", &checksummed_fixture());
    let table = BcCorrectionTable::load(&path).expect("valid BCCR fixture should pass its CRC");
    let _ = std::fs::remove_file(path);
    assert!(table.total_cells() > 0, "fixture table should have cells");
}

#[test]
fn bccr_corrupt_data_byte_rejected() {
    let mut bytes = checksummed_fixture();
    // The final four bytes are the sole data cell; corrupt one without changing dimensions.
    let data_byte = bytes.len() - 1;
    bytes[data_byte] ^= 0xFF;

    let path = write_fixture("corrupt_953", &bytes);
    let result = BcCorrectionTable::load(&path);
    let _ = std::fs::remove_file(path);

    assert!(
        matches!(result, Err(BcTableError::ChecksumMismatch)),
        "a corrupt data byte must be rejected with ChecksumMismatch (load now verifies the CRC)"
    );
}
