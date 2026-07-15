//! Import of third-party ballistic profile files.
//!
//! Currently supports the ArcherBC2 `.a7p` format (MD5-hex envelope + proto3
//! payload). The wire format is implemented here from the protobuf wire
//! specification for interoperability; no upstream schema files or code are
//! vendored (the a7p project is LGPL-3.0, this crate is MIT OR Apache-2.0).

mod a7p;
mod md5;
mod wire;

pub use a7p::{
    parse_a7p, wrap_payload, A7pBcType, A7pDocument, A7pError, A7pProfile, EnvelopeStatus,
    UnknownField,
};
