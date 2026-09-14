//! C ABI for the JSON command bridge — the entire mobile-facing surface.
//!
//! Three symbols. Ownership rules (also documented in the shipped header):
//!
//! - The caller owns the input buffer; the engine never retains it past the call.
//! - Every pointer returned by `ballistics_bridge_call` / `ballistics_bridge_call_n`
//!   is heap-allocated by the engine and must be released exactly once with
//!   `ballistics_bridge_free`. Freeing NULL is a no-op.
//! - The calls NEVER return NULL: every failure mode — invalid UTF-8, bad JSON,
//!   unknown command, command failure, internal panic — is an in-band
//!   `{"ok":false,...}` envelope.
//! - Calls are independent and thread-safe; there is no shared mutable state.
//!
//! Android callers want `ballistics_bridge_call_n` over a `ByteArray`, not
//! `ballistics_bridge_call` over a `jstring`: JNI's string conversions speak modified UTF-8,
//! which `call_with_bytes` below refuses. `docs/ANDROID_JNI_BRIDGE.md` has the worked example
//! and the reasoning; `tests::a_request_in_jni_modified_utf8_is_refused` pins the behaviour it
//! argues from.

use std::ffi::{c_char, CStr, CString};

use super::{bridge_call, BridgeErrorCode, BRIDGE_API_VERSION};

fn respond(json: String) -> *mut c_char {
    // A JSON string can legally contain no interior NULs (serde_json escapes
    // control characters), so this only fails on memory exhaustion-class bugs;
    // fall back to a static-shaped minimal envelope built from a clean literal.
    CString::new(json)
        .unwrap_or_else(|_| {
            CString::new(format!(
                r#"{{"ok":false,"api_version":{BRIDGE_API_VERSION},"error":{{"code":"internal_error","message":"response contained an interior NUL"}}}}"#
            ))
            .expect("fallback envelope contains no NUL")
        })
        .into_raw()
}

fn error_envelope(code: BridgeErrorCode, message: &str) -> String {
    // Reuse the bridge's own serializer by round-tripping through it would drag
    // private helpers into the ABI layer; a literal is simpler and stable.
    let code = serde_json::to_string(&code).unwrap_or_else(|_| "\"internal_error\"".into());
    format!(
        r#"{{"ok":false,"api_version":{BRIDGE_API_VERSION},"error":{{"code":{code},"message":"{message}"}}}}"#
    )
}

/// Process one bridge request (NUL-terminated UTF-8 JSON envelope).
///
/// # Safety
/// `request_json` must be NULL or a valid NUL-terminated C string. The returned
/// pointer must be freed with [`ballistics_bridge_free`] exactly once.
#[no_mangle]
pub unsafe extern "C" fn ballistics_bridge_call(request_json: *const c_char) -> *mut c_char {
    if request_json.is_null() {
        return respond(error_envelope(
            BridgeErrorCode::InvalidJson,
            "request pointer is NULL",
        ));
    }
    let bytes = CStr::from_ptr(request_json).to_bytes();
    call_with_bytes(bytes)
}

/// Length-explicit variant for callers whose buffers are not NUL-terminated.
///
/// # Safety
/// `request` must be NULL or point to at least `len` readable bytes. The returned
/// pointer must be freed with [`ballistics_bridge_free`] exactly once.
#[no_mangle]
pub unsafe extern "C" fn ballistics_bridge_call_n(
    request: *const u8,
    len: usize,
) -> *mut c_char {
    if request.is_null() {
        return respond(error_envelope(
            BridgeErrorCode::InvalidJson,
            "request pointer is NULL",
        ));
    }
    let bytes = std::slice::from_raw_parts(request, len);
    call_with_bytes(bytes)
}

fn call_with_bytes(bytes: &[u8]) -> *mut c_char {
    match std::str::from_utf8(bytes) {
        Ok(text) => respond(bridge_call(text)),
        Err(_) => respond(error_envelope(
            BridgeErrorCode::InvalidJson,
            "request is not valid UTF-8",
        )),
    }
}

/// Release a response produced by the bridge calls. NULL is a no-op.
///
/// # Safety
/// `response` must be NULL or a pointer previously returned by
/// [`ballistics_bridge_call`] / [`ballistics_bridge_call_n`] that has not already
/// been freed.
#[no_mangle]
pub unsafe extern "C" fn ballistics_bridge_free(response: *mut c_char) {
    if !response.is_null() {
        drop(CString::from_raw(response));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::ffi::CString;

    fn roundtrip(request: &str) -> serde_json::Value {
        let c_request = CString::new(request).unwrap();
        let raw = unsafe { ballistics_bridge_call(c_request.as_ptr()) };
        assert!(!raw.is_null(), "bridge must never return NULL");
        let text = unsafe { CStr::from_ptr(raw) }.to_str().unwrap().to_owned();
        unsafe { ballistics_bridge_free(raw) };
        serde_json::from_str(&text).expect("bridge output is JSON")
    }

    #[test]
    fn c_abi_smoke_meta_version() {
        let out = roundtrip(r#"{"api_version":1,"command":"meta.version"}"#);
        assert_eq!(out["ok"], true);
        assert_eq!(out["result"]["engine_version"], env!("CARGO_PKG_VERSION"));
    }

    #[test]
    fn null_request_is_an_envelope() {
        let raw = unsafe { ballistics_bridge_call(std::ptr::null()) };
        let text = unsafe { CStr::from_ptr(raw) }.to_str().unwrap().to_owned();
        unsafe { ballistics_bridge_free(raw) };
        let out: serde_json::Value = serde_json::from_str(&text).unwrap();
        assert_eq!(out["ok"], false);
        assert_eq!(out["error"]["code"], "invalid_json");
    }

    #[test]
    fn length_variant_handles_non_terminated_buffers() {
        let payload = br#"{"api_version":1,"command":"meta.version"}"#;
        let raw = unsafe { ballistics_bridge_call_n(payload.as_ptr(), payload.len()) };
        let text = unsafe { CStr::from_ptr(raw) }.to_str().unwrap().to_owned();
        unsafe { ballistics_bridge_free(raw) };
        let out: serde_json::Value = serde_json::from_str(&text).unwrap();
        assert_eq!(out["ok"], true);
    }

    #[test]
    fn invalid_utf8_is_an_envelope() {
        let payload = [0xFFu8, 0xFE, 0x00];
        let raw = unsafe { ballistics_bridge_call_n(payload.as_ptr(), payload.len()) };
        let text = unsafe { CStr::from_ptr(raw) }.to_str().unwrap().to_owned();
        unsafe { ballistics_bridge_free(raw) };
        let out: serde_json::Value = serde_json::from_str(&text).unwrap();
        assert_eq!(out["error"]["code"], "invalid_json");
    }

    #[test]
    fn free_null_is_a_no_op() {
        unsafe { ballistics_bridge_free(std::ptr::null_mut()) };
    }

    /// Encode `text` the way JNI's `GetStringUTFChars` and `GetStringUTFRegion` do: MODIFIED
    /// UTF-8, which is UTF-8 applied to UTF-16 code units rather than to scalar values. A
    /// supplementary character therefore comes out as its two surrogates, each in the 3-byte
    /// form (six bytes, where UTF-8 uses four), and U+0000 comes out as `C0 80`.
    fn modified_utf8(text: &str) -> Vec<u8> {
        let mut out = Vec::new();
        for unit in text.encode_utf16() {
            match unit {
                0 => out.extend_from_slice(&[0xC0, 0x80]),
                0x0001..=0x007F => out.push(unit as u8),
                0x0080..=0x07FF => {
                    out.push(0xC0 | (unit >> 6) as u8);
                    out.push(0x80 | (unit & 0x3F) as u8);
                }
                _ => {
                    out.push(0xE0 | (unit >> 12) as u8);
                    out.push(0x80 | ((unit >> 6) & 0x3F) as u8);
                    out.push(0x80 | (unit & 0x3F) as u8);
                }
            }
        }
        out
    }

    /// MBA-1543: the bridge refuses JNI's modified UTF-8, which is what makes
    /// `ballistics_bridge_call_n` over a `ByteArray` the Android entry point rather than a
    /// stylistic preference. `docs/ANDROID_JNI_BRIDGE.md` argues from exactly this behaviour,
    /// so it is driven here rather than asserted there.
    ///
    /// The discriminator is the error CODE, not the fact of an error: the command in this
    /// request does not exist either way, so getting as far as `unknown_command` proves the
    /// bytes were decoded, and `invalid_json` proves they were not.
    #[test]
    fn a_request_in_jni_modified_utf8_is_refused() {
        // U+1F3AF is supplementary, which is where the two encodings diverge.
        let request = r#"{"api_version":1,"command":"no.such.command.🎯"}"#;
        let standard = request.as_bytes();
        let modified = modified_utf8(request);
        assert_ne!(
            modified.as_slice(),
            standard,
            "the fixture must actually exercise the divergence"
        );

        let decoded = |bytes: &[u8]| {
            let raw = unsafe { ballistics_bridge_call_n(bytes.as_ptr(), bytes.len()) };
            let text = unsafe { CStr::from_ptr(raw) }.to_str().unwrap().to_owned();
            unsafe { ballistics_bridge_free(raw) };
            serde_json::from_str::<serde_json::Value>(&text).unwrap()
        };

        let out = decoded(standard);
        assert_eq!(out["error"]["code"], "unknown_command", "{out}");

        let out = decoded(&modified);
        assert_eq!(out["error"]["code"], "invalid_json", "{out}");
        assert_eq!(out["error"]["message"], "request is not valid UTF-8", "{out}");
    }

    /// The scope of that refusal, so the docs do not overstate it: the two encodings agree on
    /// everything from U+0001 to U+FFFF, so ordinary accented text goes through JNI's own
    /// conversion unharmed and the defect hides until someone types outside the BMP.
    #[test]
    fn modified_utf8_and_utf8_agree_below_the_supplementary_planes() {
        for text in ["ascii", "München", "Ω", "中文", "\u{FFFF}"] {
            assert_eq!(
                modified_utf8(text),
                text.as_bytes(),
                "{text} must encode identically in both"
            );
        }
        for text in ["🎯", "\u{10000}", "\u{10FFFF}"] {
            assert_ne!(modified_utf8(text), text.as_bytes(), "{text} must diverge");
        }
        // And the other member of the divergence: an interior NUL, which modified UTF-8 exists
        // to keep out of the byte stream in the first place.
        assert_eq!(modified_utf8("\0"), vec![0xC0, 0x80]);
    }
}
