//! Import of third-party ballistic profile files.
//!
//! Currently supports the ArcherBC2 `.a7p` format (MD5-hex envelope + proto3
//! payload). The wire format is implemented here from the protobuf wire
//! specification for interoperability; no upstream schema files or code are
//! vendored (the a7p project is LGPL-3.0, this crate is MIT OR Apache-2.0).

// Task 1 (MBA-1323) lands this module ahead of its consumer: the a7p parser
// added in Task 2 is what calls md5_hex/parse_message/etc. Until then these
// pub(crate) items have no caller, so silence dead_code rather than let the
// lint gate flag intentionally-not-yet-wired scaffolding.
#[allow(dead_code)]
mod md5;
#[allow(dead_code)]
mod wire;
