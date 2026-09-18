//! Minimal proto3 wire-format WRITER — just enough to encode `.a7p` payloads.
//!
//! Written from the protobuf wire-format specification, like its reading
//! counterpart in `profile_import::wire`; no upstream code or schema files are
//! vendored.
//!
//! Deliberately NOT a mirror image of the reader that shares its helpers: the
//! reader is `pub(crate)` inside the `profile_import` subtree and stays there.
//! Keeping the writer independent means the round-trip test in `a7p.rs` is a
//! genuine cross-check of two separately-derived implementations rather than a
//! function proving it can undo itself — the same reason the reader's own tests
//! carry spec-derived encoders instead of reusing a shared one.

/// Base-128 varint, little-endian groups of seven bits, high bit = continue.
pub(crate) fn write_varint(mut value: u64, out: &mut Vec<u8>) {
    loop {
        let byte = (value & 0x7f) as u8;
        value >>= 7;
        if value == 0 {
            out.push(byte);
            return;
        }
        out.push(byte | 0x80);
    }
}

fn write_key(number: u32, wire_type: u8, out: &mut Vec<u8>) {
    write_varint((u64::from(number) << 3) | u64::from(wire_type), out);
}

/// proto3 `int32`: the value is sign-extended to 64 bits before encoding, so a
/// negative int32 occupies the full ten bytes. This is the encoder half of
/// `profile_import::wire::varint_to_i32`'s truncating cast — do not "optimize"
/// it to `value as u32`, which would produce a five-byte varint that canonical
/// decoders read back as a large positive number.
pub(crate) fn write_i32_field(number: u32, value: i32, out: &mut Vec<u8>) {
    write_key(number, 0, out);
    write_varint(i64::from(value) as u64, out);
}

pub(crate) fn write_bytes_field(number: u32, payload: &[u8], out: &mut Vec<u8>) {
    write_key(number, 2, out);
    write_varint(payload.len() as u64, out);
    out.extend_from_slice(payload);
}

pub(crate) fn write_string_field(number: u32, value: &str, out: &mut Vec<u8>) {
    write_bytes_field(number, value.as_bytes(), out);
}

/// `repeated int32` in the PACKED encoding (proto3's default for repeated
/// scalars). The reader accepts both packed and unpacked; we emit packed
/// because that is what a canonical proto3 serializer produces and therefore
/// what other tools in the ecosystem are most likely to have been tested
/// against.
pub(crate) fn write_packed_i32_field(number: u32, values: &[i32], out: &mut Vec<u8>) {
    let mut payload = Vec::with_capacity(values.len() * 2);
    for &value in values {
        write_varint(i64::from(value) as u64, &mut payload);
    }
    write_bytes_field(number, &payload, out);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn varint_matches_the_spec_examples() {
        let cases: [(u64, &[u8]); 4] = [
            (0, &[0x00]),
            (1, &[0x01]),
            (127, &[0x7f]),
            (300, &[0xac, 0x02]),
        ];
        for (value, expected) in cases {
            let mut out = Vec::new();
            write_varint(value, &mut out);
            assert_eq!(out, expected, "varint {value}");
        }
    }

    #[test]
    fn negative_int32_is_sign_extended_to_ten_bytes() {
        // The trap this pins: `value as u32` would encode -5 in five bytes and
        // decode as 4294967291 in any canonical reader.
        let mut out = Vec::new();
        write_i32_field(12, -5, &mut out);
        assert_eq!(out.len(), 1 + 10, "key byte + ten varint bytes");
    }

    #[test]
    fn packed_repeated_is_one_length_delimited_field() {
        let mut out = Vec::new();
        write_packed_i32_field(26, &[10_000, 20_000], &mut out);
        // Spelled out byte for byte, because the two-byte KEY is the part that is
        // easy to get wrong: field 26 with wire type 2 is 210, which does not fit
        // in one varint byte.
        assert_eq!(
            out,
            vec![
                0xd2, 0x01, // key: (26 << 3) | 2 = 210
                0x05, // payload length
                0x90, 0x4e, // 10000
                0xa0, 0x9c, 0x01, // 20000
            ]
        );
    }
}
