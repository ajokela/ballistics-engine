//! Process transport for the `ballistics solve-json` command.
//!
//! This module deliberately contains only stdin/stdout framing. Request semantics remain in the
//! transport-free library service, [`ballistics_engine::solve_v1()`].

use ballistics_engine::solve_json::{
    decode_solve_request_v1, SolveErrorCodeV1, SolveErrorEnvelopeV1, SolveErrorV1, SolveSuccessV1,
};
use std::io::{self, Read, Write};
use std::panic::{catch_unwind, AssertUnwindSafe};

const MAX_INPUT_BYTES: usize = 1024 * 1024;
const EXIT_SUCCESS: i32 = 0;
const EXIT_IO_OR_INTERNAL: i32 = 1;
const EXIT_REQUEST: i32 = 2;
const EXIT_RESOURCE_OR_SOLVE: i32 = 3;

type ProtocolResult = Result<SolveSuccessV1, SolveErrorEnvelopeV1>;

/// Run one solve-json exchange on the process standard streams.
pub(crate) fn run_stdio() -> i32 {
    let stdin = io::stdin();
    let stdout = io::stdout();
    run_transport(stdin.lock(), stdout.lock(), process_request)
}

fn process_request(input: &[u8]) -> ProtocolResult {
    let input = std::str::from_utf8(input).map_err(|error| invalid_utf8_error(input, error))?;
    let request = decode_solve_request_v1(input)?;
    ballistics_engine::solve_v1(request)
}

fn run_transport<R, W, F>(mut reader: R, mut writer: W, processor: F) -> i32
where
    R: Read,
    W: Write,
    F: FnOnce(&[u8]) -> ProtocolResult,
{
    // Serialize completely before the first write so a serialization failure can never append an
    // error envelope after a partial success document. AssertUnwindSafe is appropriate here: all
    // captured state is process-local and discarded after this single exchange.
    let guarded = catch_unwind(AssertUnwindSafe(|| {
        let result = read_request(&mut reader).and_then(|bytes| processor(&bytes));
        let exit_code = result
            .as_ref()
            .map_or_else(error_exit_code, |_| EXIT_SUCCESS);
        match encode_result(&result) {
            Ok(payload) => (payload, exit_code),
            Err(_) => (internal_error_payload(), EXIT_IO_OR_INTERNAL),
        }
    }));

    let (mut payload, exit_code) = match guarded {
        Ok(response) => response,
        Err(_) => (internal_error_payload(), EXIT_IO_OR_INTERNAL),
    };
    payload.push(b'\n');

    if writer
        .write_all(&payload)
        .and_then(|()| writer.flush())
        .is_err()
    {
        // If stdout itself is unavailable, delivering a protocol envelope is physically
        // impossible. The process still returns the documented I/O status.
        return EXIT_IO_OR_INTERNAL;
    }

    exit_code
}

fn read_request(reader: &mut impl Read) -> Result<Vec<u8>, SolveErrorEnvelopeV1> {
    let mut bytes = Vec::new();
    reader
        .take((MAX_INPUT_BYTES + 1) as u64)
        .read_to_end(&mut bytes)
        .map_err(|error| {
            protocol_error(
                SolveErrorCodeV1::IoError,
                format!("failed to read solve-json request from stdin: {error}"),
            )
        })?;

    if bytes.len() > MAX_INPUT_BYTES {
        return Err(protocol_error(
            SolveErrorCodeV1::ResourceLimit,
            format!("solve-json input exceeds the {MAX_INPUT_BYTES}-byte limit"),
        ));
    }

    Ok(bytes)
}

fn encode_result(result: &ProtocolResult) -> Result<Vec<u8>, serde_json::Error> {
    match result {
        Ok(success) => serde_json::to_vec(success),
        Err(error) => serde_json::to_vec(error),
    }
}

fn protocol_error(code: SolveErrorCodeV1, message: impl Into<String>) -> SolveErrorEnvelopeV1 {
    SolveErrorEnvelopeV1::new(SolveErrorV1::new(code, message))
}

fn invalid_utf8_error(input: &[u8], error: std::str::Utf8Error) -> SolveErrorEnvelopeV1 {
    let valid_prefix = &input[..error.valid_up_to()];
    let line = valid_prefix.iter().filter(|&&byte| byte == b'\n').count() + 1;
    let current_line = valid_prefix
        .rsplit(|&byte| byte == b'\n')
        .next()
        .unwrap_or_default();
    // `valid_up_to` guarantees that this prefix is valid UTF-8.
    let column = std::str::from_utf8(current_line)
        .map(|line| line.chars().count() + 1)
        .unwrap_or(1);
    let error = SolveErrorV1::new(
        SolveErrorCodeV1::InvalidJson,
        "solve-json input is not valid UTF-8 JSON",
    )
    .at_location(line, column)
    .unwrap_or_else(|_| {
        SolveErrorV1::new(
            SolveErrorCodeV1::InternalError,
            "failed to locate invalid UTF-8 input",
        )
    });
    SolveErrorEnvelopeV1::new(error)
}

fn error_exit_code(error: &SolveErrorEnvelopeV1) -> i32 {
    match error.error.code {
        SolveErrorCodeV1::InvalidJson
        | SolveErrorCodeV1::UnsupportedSchemaVersion
        | SolveErrorCodeV1::UnknownField
        | SolveErrorCodeV1::MissingField
        | SolveErrorCodeV1::InvalidValue
        | SolveErrorCodeV1::ConflictingFields => EXIT_REQUEST,
        SolveErrorCodeV1::ResourceLimit | SolveErrorCodeV1::SolveFailed => EXIT_RESOURCE_OR_SOLVE,
        SolveErrorCodeV1::IoError | SolveErrorCodeV1::InternalError => EXIT_IO_OR_INTERNAL,
    }
}

fn internal_error_payload() -> Vec<u8> {
    let envelope = protocol_error(
        SolveErrorCodeV1::InternalError,
        "solve-json command failed unexpectedly",
    );
    serde_json::to_vec(&envelope).unwrap_or_else(|_| {
        br#"{"schema_version":1,"status":"error","error":{"code":"internal_error","message":"solve-json command failed unexpectedly","path":null,"line":null,"column":null}}"#
            .to_vec()
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    struct FailingReader;

    impl Read for FailingReader {
        fn read(&mut self, _buffer: &mut [u8]) -> io::Result<usize> {
            Err(io::Error::other("injected read failure"))
        }
    }

    fn parse_output(output: &[u8]) -> serde_json::Value {
        assert!(output.ends_with(b"\n"), "envelope needs one terminator");
        assert_eq!(
            output.iter().filter(|&&byte| byte == b'\n').count(),
            1,
            "transport must write exactly one compact envelope"
        );
        serde_json::from_slice(output).expect("transport output must be one JSON envelope")
    }

    #[test]
    fn stdin_failure_is_an_io_envelope() {
        let mut output = Vec::new();
        let code = run_transport(FailingReader, &mut output, process_request);
        assert_eq!(code, EXIT_IO_OR_INTERNAL);
        let wire = parse_output(&output);
        assert_eq!(wire["status"], "error");
        assert_eq!(wire["error"]["code"], "io_error");
        assert!(wire["error"]["message"]
            .as_str()
            .is_some_and(|message| message.contains("injected read failure")));
    }

    #[test]
    fn processor_panic_is_a_generic_internal_error() {
        let mut output = Vec::new();
        let code = run_transport(&b"{}"[..], &mut output, |_| {
            panic!("private panic payload and backtrace marker")
        });
        assert_eq!(code, EXIT_IO_OR_INTERNAL);
        let wire = parse_output(&output);
        assert_eq!(wire["error"]["code"], "internal_error");
        let serialized = String::from_utf8(output).expect("UTF-8 envelope");
        assert!(!serialized.contains("private panic payload"));
        assert!(!serialized.contains("backtrace"));
    }

    #[test]
    fn every_protocol_error_code_has_a_documented_exit_class() {
        for (code, expected) in [
            (SolveErrorCodeV1::InvalidJson, EXIT_REQUEST),
            (SolveErrorCodeV1::UnsupportedSchemaVersion, EXIT_REQUEST),
            (SolveErrorCodeV1::UnknownField, EXIT_REQUEST),
            (SolveErrorCodeV1::MissingField, EXIT_REQUEST),
            (SolveErrorCodeV1::InvalidValue, EXIT_REQUEST),
            (SolveErrorCodeV1::ConflictingFields, EXIT_REQUEST),
            (SolveErrorCodeV1::ResourceLimit, EXIT_RESOURCE_OR_SOLVE),
            (SolveErrorCodeV1::SolveFailed, EXIT_RESOURCE_OR_SOLVE),
            (SolveErrorCodeV1::IoError, EXIT_IO_OR_INTERNAL),
            (SolveErrorCodeV1::InternalError, EXIT_IO_OR_INTERNAL),
        ] {
            assert_eq!(error_exit_code(&protocol_error(code, "test")), expected);
        }
    }
}
