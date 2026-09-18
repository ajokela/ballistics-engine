//! Export of saved profiles into third-party ballistic profile files.
//!
//! Currently the ArcherBC2 `.a7p` format — the outbound half of
//! [`crate::profile_import`], written from the same protobuf wire specification
//! with no upstream schema files or code vendored (the a7p project is LGPL-3.0,
//! this crate is MIT OR Apache-2.0).
//!
//! Export is LOSSY in one direction the importer never had to worry about: a
//! saved profile holds a good deal the format has no slot for. Every export
//! therefore returns the bytes AND [`NotCarried`] — the full list of
//! `ProfileData` fields that did not make it, and whether this particular
//! profile had anything in them. See `a7p.rs`'s module doc for why the list is
//! unconditional.

mod a7p;
mod wire;

pub use a7p::{export_a7p, A7pExport, A7pExportError, NotCarried, CARRIED_FIELDS};
