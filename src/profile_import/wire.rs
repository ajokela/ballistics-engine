//! Minimal proto3 wire-format reader — just enough to decode `.a7p` payloads.
//!
//! Implemented from the protobuf wire-format specification for
//! interoperability; no upstream code or schema files are vendored. Unknown
//! fields are represented, not dropped, so the importer can report what it
//! did not understand.

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum WireValue<'a> {
    /// Wire type 0 — raw varint value (int32/int64/enum/bool).
    Varint(u64),
    /// Wire type 1 — fixed 64-bit payload.
    Fixed64([u8; 8]),
    /// Wire type 2 — length-delimited: strings, bytes, sub-messages, packed repeats.
    Bytes(&'a [u8]),
    /// Wire type 5 — fixed 32-bit payload.
    Fixed32([u8; 4]),
}

#[derive(Debug)]
pub(crate) struct WireField<'a> {
    pub number: u32,
    pub value: WireValue<'a>,
}

#[derive(Debug, PartialEq)]
pub(crate) enum WireError {
    /// Input ended inside a field.
    Truncated,
    /// A varint ran past 64 bits.
    VarintOverflow,
    /// Deprecated group wire types (3/4) or reserved values (6/7).
    UnsupportedWireType(u8),
}

impl std::fmt::Display for WireError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            WireError::Truncated => write!(f, "input truncated mid-field"),
            WireError::VarintOverflow => write!(f, "varint exceeds 64 bits"),
            WireError::UnsupportedWireType(t) => write!(f, "unsupported wire type {t}"),
        }
    }
}

impl std::error::Error for WireError {}

pub(crate) fn read_varint(buf: &[u8], pos: &mut usize) -> Result<u64, WireError> {
    let mut value = 0u64;
    let mut shift = 0u32;
    loop {
        let byte = *buf.get(*pos).ok_or(WireError::Truncated)?;
        *pos += 1;
        if shift >= 64 {
            return Err(WireError::VarintOverflow);
        }
        value |= u64::from(byte & 0x7f) << shift;
        if byte & 0x80 == 0 {
            return Ok(value);
        }
        shift += 7;
    }
}

/// Split one message level into its fields. Callers descend into
/// `WireValue::Bytes` for the sub-messages they know by field number.
pub(crate) fn parse_message(buf: &[u8]) -> Result<Vec<WireField<'_>>, WireError> {
    let mut pos = 0usize;
    let mut fields = Vec::new();
    while pos < buf.len() {
        let key = read_varint(buf, &mut pos)?;
        let number = (key >> 3) as u32;
        let wire_type = (key & 7) as u8;
        let value = match wire_type {
            0 => WireValue::Varint(read_varint(buf, &mut pos)?),
            1 => {
                let end = pos.checked_add(8).ok_or(WireError::Truncated)?;
                if end > buf.len() {
                    return Err(WireError::Truncated);
                }
                let mut b = [0u8; 8];
                b.copy_from_slice(&buf[pos..end]);
                pos = end;
                WireValue::Fixed64(b)
            }
            2 => {
                let len = read_varint(buf, &mut pos)? as usize;
                let end = pos.checked_add(len).ok_or(WireError::Truncated)?;
                if end > buf.len() {
                    return Err(WireError::Truncated);
                }
                let v = &buf[pos..end];
                pos = end;
                WireValue::Bytes(v)
            }
            5 => {
                let end = pos.checked_add(4).ok_or(WireError::Truncated)?;
                if end > buf.len() {
                    return Err(WireError::Truncated);
                }
                let mut b = [0u8; 4];
                b.copy_from_slice(&buf[pos..end]);
                pos = end;
                WireValue::Fixed32(b)
            }
            other => return Err(WireError::UnsupportedWireType(other)),
        };
        fields.push(WireField { number, value });
    }
    Ok(fields)
}

/// proto3 `int32` semantics: the wire value is the 64-bit two's complement;
/// truncating to i32 recovers the original (matches canonical decoders).
pub(crate) fn varint_to_i32(raw: u64) -> i32 {
    raw as i32
}

/// Collect a `repeated int32`, accepting both packed (length-delimited) and
/// unpacked (repeated varint) encodings — proto3 serializers may emit either.
pub(crate) fn collect_repeated_i32(
    fields: &[WireField<'_>],
    number: u32,
) -> Result<Vec<i32>, WireError> {
    let mut out = Vec::new();
    for field in fields.iter().filter(|f| f.number == number) {
        match &field.value {
            WireValue::Varint(v) => out.push(varint_to_i32(*v)),
            WireValue::Bytes(b) => {
                let mut pos = 0usize;
                while pos < b.len() {
                    let v = read_varint(b, &mut pos)?;
                    out.push(varint_to_i32(v));
                }
            }
            _ => {}
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Test-only encoder helpers. Deliberately written against the wire SPEC,
    /// independently of the decoder, so the two cross-check each other.
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
    fn enc_field_varint(number: u32, value: u64, out: &mut Vec<u8>) {
        enc_varint(u64::from(number) << 3, out);
        enc_varint(value, out);
    }
    fn enc_field_bytes(number: u32, payload: &[u8], out: &mut Vec<u8>) {
        enc_varint((u64::from(number) << 3) | 2, out);
        enc_varint(payload.len() as u64, out);
        out.extend_from_slice(payload);
    }

    #[test]
    fn varint_roundtrip_boundaries() {
        for v in [0u64, 1, 127, 128, 300, u64::from(u32::MAX), u64::MAX] {
            let mut buf = Vec::new();
            enc_varint(v, &mut buf);
            let mut pos = 0;
            assert_eq!(read_varint(&buf, &mut pos).unwrap(), v);
            assert_eq!(pos, buf.len());
        }
    }

    #[test]
    fn truncated_varint_is_an_error() {
        let mut pos = 0;
        assert_eq!(read_varint(&[0x80], &mut pos), Err(WireError::Truncated));
    }

    #[test]
    fn parses_mixed_fields_and_preserves_order() {
        let mut buf = Vec::new();
        enc_field_varint(11, 7920, &mut buf); // c_muzzle_velocity-like
        enc_field_bytes(1, b"hello", &mut buf); // string
        let fields = parse_message(&buf).unwrap();
        assert_eq!(fields.len(), 2);
        assert_eq!(fields[0].number, 11);
        assert_eq!(fields[0].value, WireValue::Varint(7920));
        assert_eq!(fields[1].number, 1);
        assert_eq!(fields[1].value, WireValue::Bytes(b"hello"));
    }

    #[test]
    fn negative_int32_arrives_as_ten_byte_varint() {
        // proto3 int32 -5 is encoded as the 64-bit two's complement varint.
        let mut buf = Vec::new();
        enc_field_varint(12, (-5i64) as u64, &mut buf);
        let fields = parse_message(&buf).unwrap();
        let WireValue::Varint(raw) = fields[0].value else {
            panic!("expected varint")
        };
        assert_eq!(varint_to_i32(raw), -5);
    }

    #[test]
    fn repeated_i32_accepts_packed_and_unpacked() {
        // packed: one length-delimited blob of varints
        let mut packed_payload = Vec::new();
        for v in [10_000u64, 20_000, 25_000] {
            enc_varint(v, &mut packed_payload);
        }
        let mut buf = Vec::new();
        enc_field_bytes(26, &packed_payload, &mut buf);
        // unpacked: one more bare varint occurrence of the same field
        enc_field_varint(26, 30_000, &mut buf);
        let fields = parse_message(&buf).unwrap();
        assert_eq!(
            collect_repeated_i32(&fields, 26).unwrap(),
            vec![10_000, 20_000, 25_000, 30_000]
        );
    }

    #[test]
    fn fixed32_and_fixed64_are_consumed_not_confused() {
        let mut buf = Vec::new();
        enc_varint((3u64 << 3) | 5, &mut buf); // field 3, fixed32
        buf.extend_from_slice(&[1, 2, 3, 4]);
        enc_varint((4u64 << 3) | 1, &mut buf); // field 4, fixed64
        buf.extend_from_slice(&[5, 6, 7, 8, 9, 10, 11, 12]);
        let fields = parse_message(&buf).unwrap();
        assert_eq!(fields[0].value, WireValue::Fixed32([1, 2, 3, 4]));
        assert_eq!(fields[1].value, WireValue::Fixed64([5, 6, 7, 8, 9, 10, 11, 12]));
    }

    #[test]
    fn truncated_length_delimited_is_an_error() {
        let mut buf = Vec::new();
        enc_varint((1u64 << 3) | 2, &mut buf);
        enc_varint(100, &mut buf); // claims 100 bytes
        buf.extend_from_slice(b"short");
        assert_eq!(parse_message(&buf).unwrap_err(), WireError::Truncated);
    }

    #[test]
    fn group_wire_types_are_rejected() {
        let mut buf = Vec::new();
        enc_varint((1u64 << 3) | 3, &mut buf); // deprecated start-group
        assert_eq!(
            parse_message(&buf).unwrap_err(),
            WireError::UnsupportedWireType(3)
        );
    }
}
