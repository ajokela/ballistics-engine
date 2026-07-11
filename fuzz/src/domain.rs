//! Shared byte->input generators and invariant assertions for the fuzz harnesses.

use arbitrary::{Result, Unstructured};

/// Map raw bytes to a finite f64 in `[lo, hi]` (inclusive-ish, deterministic).
pub fn ranged(u: &mut Unstructured, lo: f64, hi: f64) -> Result<f64> {
    let bits: u64 = u.arbitrary()?;
    let frac = (bits as f64) / (u64::MAX as f64); // [0, 1]
    Ok(lo + frac * (hi - lo))
}
