# Multi-BC Support via `bridge_call()` — Required Engine Changes

## Overview

`bridge_call()` accepts and returns JSON strings through the **solve-json v1 schema**
(`src/solve_json.rs` + `src/solve_v1.rs`). This document describes every change required
inside the `ballistics-engine` library to make velocity-segmented ballistic coefficients
(multi-BC) reachable through that JSON path — without modifying `bridge_call()` itself.

---

## Root Cause — Why Multi-BC Is Currently Unreachable

`solve_v1.rs` lines 650–652 always hard-code:

```rust
use_bc_segments: false,
bc_segments: None,
bc_segments_data: None,
```

The `apply_bc5d_correction()` function (line 787–870) is the **only code path** that ever
sets `use_bc_segments = true` and populates `bc_segments_data`. That function reads from
`request.corrections.bc5d_table_path` — a local filesystem path to a proprietary binary
`.bin` file. There is no JSON field anywhere in the v1 schema that lets a caller supply
BC segments directly as data values.

---

## Required Changes — Engine Library (`ballistics-engine`)

### 1. `src/solve_json.rs` — New struct `BcSegmentV1`

Add a new DTO representing one velocity band. Place it near `DragTableV1` / `DragPointV1`:

```rust
/// One velocity-segmented ballistic coefficient for multi-BC flight (JSON bridge path).
///
/// Segments must be supplied in **descending** order of `velocity_from_mps`
/// (highest velocity first), matching `BCSegmentData` conventions internally.
///
/// When `projectile.bc_segments` is present, the scalar `ballistic_coefficient`
/// is used only as the gap-coverage fallback; it is otherwise not flown.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BcSegmentV1 {
    /// Lower bound of this velocity band, m/s (inclusive). Must be finite and >= 0.
    pub velocity_from_mps: f64,
    /// Upper bound, m/s (exclusive). Must be finite and > velocity_from_mps.
    pub velocity_to_mps: f64,
    /// Ballistic coefficient for this band (lb/in²). Must be finite and > 0.
    pub ballistic_coefficient: f64,
}
```

---

### 2. `src/solve_json.rs` — Extend `ProjectileV1`

Add the optional `bc_segments` field to the request DTO:

```rust
pub struct ProjectileV1 {
    pub mass_kg: f64,
    pub diameter_m: f64,
    pub length_m: Option<f64>,
    pub drag_model: DragModelV1,
    pub ballistic_coefficient: f64,
    pub drag_table: Option<DragTableV1>,           // existing (MBA-1597)
    /// Velocity-segmented BCs. Mutually exclusive with `drag_table`
    /// and `corrections.bc5d_table_path`.
    pub bc_segments: Option<Vec<BcSegmentV1>>,     // NEW
}
```

---

### 3. `src/solve_json.rs` — Extend `ResolvedProjectileV1`

Echo the field on the resolved DTO so the perturbation kernel can re-solve correctly:

```rust
pub struct ResolvedProjectileV1 {
    pub mass_kg: f64,
    pub diameter_m: f64,
    pub length_m: Option<f64>,
    pub drag_model: DragModelV1,
    pub ballistic_coefficient: f64,
    pub drag_table: Option<DragTableV1>,           // existing
    /// Verbatim echo of `bc_segments` when supplied, so a resolved-request
    /// re-solve (perturbation, round-trip) re-applies the same ladder.
    pub bc_segments: Option<Vec<BcSegmentV1>>,     // NEW
}
```

---

### 4. `src/solve_json.rs` — Shape validation

#### 4a. Extend `validate_projectile()`

Add `"bc_segments"` to the allowed-member list and call the new shape checker:

```rust
validate_members(
    object, path,
    &["mass_kg","diameter_m","length_m","drag_model",
      "ballistic_coefficient","drag_table",
      "bc_segments"],       // NEW entry
    &["mass_kg","diameter_m","drag_model","ballistic_coefficient"],
)?;
// ...existing validates unchanged...
if let Some(bc_segments) = object.get("bc_segments") {
    validate_bc_segments_shape(bc_segments)?;     // NEW call
}
```

#### 4b. New function `validate_bc_segments_shape()`

```rust
fn validate_bc_segments_shape(value: &Value) -> Result<(), SolveErrorEnvelopeV1> {
    let path = "$.projectile.bc_segments";
    let Some(segments) = value.as_array() else {
        return Err(protocol_error(SolveErrorCodeV1::InvalidValue,
            "bc_segments must be an array", path));
    };
    for (index, seg) in segments.iter().enumerate() {
        let seg_path = format!("{path}[{index}]");
        let seg = require_object(seg, &seg_path)?;
        validate_members(seg, &seg_path,
            &["velocity_from_mps","velocity_to_mps","ballistic_coefficient"],
            &["velocity_from_mps","velocity_to_mps","ballistic_coefficient"])?;
        validate_required_numbers(seg, &seg_path,
            &["velocity_from_mps","velocity_to_mps","ballistic_coefficient"])?;
    }
    Ok(())
}
```

---

### 5. `src/solve_json.rs` — Value-range validation in `validate_request_ranges()`

Add a block that checks segment values and enforces mutual-exclusion rules:

```rust
if let Some(segments) = &request.projectile.bc_segments {
    if segments.is_empty() {
        return Err(invalid_value("$.projectile.bc_segments",
            "must contain at least one segment"));
    }
    for (i, seg) in segments.iter().enumerate() {
        let base = format!("$.projectile.bc_segments[{i}]");
        if !seg.velocity_from_mps.is_finite() || seg.velocity_from_mps < 0.0 {
            return Err(invalid_value(&format!("{base}.velocity_from_mps"),
                "must be finite and >= 0"));
        }
        if !seg.velocity_to_mps.is_finite()
            || seg.velocity_to_mps <= seg.velocity_from_mps
        {
            return Err(invalid_value(&format!("{base}.velocity_to_mps"),
                "must be finite and > velocity_from_mps"));
        }
        if !seg.ballistic_coefficient.is_finite() || seg.ballistic_coefficient <= 0.0 {
            return Err(invalid_value(&format!("{base}.ballistic_coefficient"),
                "must be finite and > 0"));
        }
    }
    // Mutually exclusive with a custom drag_table
    if request.projectile.drag_table.is_some() {
        return Err(conflicting_fields(
            "$.projectile.bc_segments", "$.projectile.drag_table",
            "bc_segments cannot be combined with a custom drag_table"));
    }
    // Mutually exclusive with a BC5D correction table
    if request.corrections
        .as_ref()
        .and_then(|c| c.bc5d_table_path.as_deref())
        .is_some()
    {
        return Err(conflicting_fields(
            "$.projectile.bc_segments", "$.corrections.bc5d_table_path",
            "bc_segments cannot be combined with a BC5D correction table"));
    }
}
```

---

### 6. `src/solve_v1.rs` — Echo in `resolve_projectile()`

The function that builds `ResolvedProjectileV1` must carry the segments through:

```rust
let resolved_projectile = ResolvedProjectileV1 {
    mass_kg: projectile.mass_kg,
    diameter_m: projectile.diameter_m,
    length_m: projectile.length_m,
    drag_model: projectile.drag_model,
    ballistic_coefficient: projectile.ballistic_coefficient,
    drag_table: projectile.drag_table.clone(),
    bc_segments: projectile.bc_segments.clone(),   // NEW
};
```

---

### 7. `src/solve_v1.rs` — New function `apply_bc_segments()`

Add directly after `apply_bc5d_correction()`. Converts m/s → fps (the unit
`BCSegmentData` uses internally) and writes into `BallisticInputs`:

```rust
/// Apply `projectile.bc_segments` from the resolved request to the prepared engine
/// inputs, mirroring the segment-injection that `apply_bc5d_correction` performs for
/// BC5D tables. The caller supplies SI (m/s); `BCSegmentData` stores fps.
///
/// When the field is absent this is a no-op — byte-identical to every prior request.
fn apply_bc_segments(
    resolved_request: &ResolvedSolveRequestV1,
    inputs: &mut BallisticInputs,
) -> Result<(), SolveErrorEnvelopeV1> {
    let Some(segments) = resolved_request.projectile.bc_segments.as_deref() else {
        return Ok(());
    };
    const MPS_TO_FPS: f64 = 1.0 / 0.3048;
    let engine_segments: Vec<BCSegmentData> = segments
        .iter()
        .map(|s| BCSegmentData {
            velocity_min: s.velocity_from_mps * MPS_TO_FPS,
            velocity_max: s.velocity_to_mps   * MPS_TO_FPS,
            bc_value:     s.ballistic_coefficient,
        })
        .collect();
    inputs.use_bc_segments = true;
    inputs.bc_segments_data = Some(engine_segments);
    Ok(())
}
```

---

### 8. `src/solve_v1.rs` — Call `apply_bc_segments()` from `prepare_request()`

Add immediately after the existing `apply_bc5d_correction()` call (line 716):

```rust
apply_bc5d_correction(request, &resolved_request, &mut inputs, &mut warnings)?;
apply_bc_segments(&resolved_request, &mut inputs)?;   // NEW
```

Because `prepare_request()` feeds **both** the zero-angle search (Pass 0) and the
trajectory solve through the same `BallisticInputs`, the BC segment ladder is
automatically active for all passes — the pinned muzzle angle from the zero solve
uses the same drag schedule as the flight, which is required for correct drops.

---

## Summary Table — Exact Files and Sections

| File | Location | Change |
|---|---|---|
| `src/solve_json.rs` | New struct `BcSegmentV1` | Define the per-segment wire DTO |
| `src/solve_json.rs` | `ProjectileV1` | Add `pub bc_segments: Option<Vec<BcSegmentV1>>` |
| `src/solve_json.rs` | `ResolvedProjectileV1` | Add `pub bc_segments: Option<Vec<BcSegmentV1>>` |
| `src/solve_json.rs` | `validate_projectile()` | Allow `"bc_segments"` key; call `validate_bc_segments_shape()` |
| `src/solve_json.rs` | New fn `validate_bc_segments_shape()` | Shape-check the segments array |
| `src/solve_json.rs` | `validate_request_ranges()` | Range-check values; conflict-check vs `drag_table` and `bc5d_table_path` |
| `src/solve_v1.rs` | `resolve_projectile()` | Echo `bc_segments` onto `ResolvedProjectileV1` |
| `src/solve_v1.rs` | `prepare_request()` | Call `apply_bc_segments()` after `apply_bc5d_correction()` |
| `src/solve_v1.rs` | New fn `apply_bc_segments()` | Convert m/s → fps; write `use_bc_segments` + `bc_segments_data` |

---

## What Does NOT Need to Change

| Item | Reason |
|---|---|
| `BCSegmentData` struct in `lib.rs` | Already correct and public |
| `TrajectorySolver` internals | Already handles `bc_segments_data` correctly when set |
| `bridge_call()` itself | Calls `solve_v1()` → `prepare_request()`; the chain works once `prepare_request` populates the fields |
| Zero-angle search (Pass 0) | Goes through the same `prepare_request()` → same `BallisticInputs`; segments active automatically |
| `CorrectionsV1` / `SolveRequestV1` root | No root-level additions needed; everything lives inside `projectile` |
| `apply_bc5d_correction()` | Unchanged; BC5D and direct `bc_segments` are additive-exclusive alternatives |

---

## Resulting JSON Wire Format

Single-BC request (unchanged, backward-compatible — absent `bc_segments` is a no-op):

```json
{
  "schema_version": 1,
  "projectile": {
    "mass_kg": 0.000648,
    "diameter_m": 0.00455,
    "drag_model": "G1",
    "ballistic_coefficient": 0.040
  }
}
```

Multi-BC request with a two-segment ladder:

```json
{
  "schema_version": 1,
  "projectile": {
    "mass_kg": 0.000648,
    "diameter_m": 0.00455,
    "drag_model": "G1",
    "ballistic_coefficient": 0.040,
    "bc_segments": [
      { "velocity_from_mps": 274.3, "velocity_to_mps": 999.0, "ballistic_coefficient": 0.042 },
      { "velocity_from_mps": 0.0,   "velocity_to_mps": 274.3, "ballistic_coefficient": 0.038 }
    ]
  }
}
```

When `bc_segments` is absent the solve is byte-identical to every prior request.
When present, `ballistic_coefficient` becomes the gap-coverage fallback only,
mirroring exactly what `apply_bc5d_correction()` does for BC5D tables.

---

## Key Design Decisions

1. **Field location: `projectile.bc_segments`** — logically belongs to the projectile,
   additive-optional, no existing request breaks. The root-level member list in
   `validate_request_shape()` does not need updating since the field lives one level down
   inside `projectile`.

2. **Units: m/s on the wire, fps in the engine** — SI on the wire is consistent with
   every other velocity field in solve-json v1 (`muzzle_velocity_mps`, wind speeds, etc.).
   The `MPS_TO_FPS` conversion is applied once in `apply_bc_segments()`, not scattered
   across callers.

3. **`ballistic_coefficient` as gap fallback** — identical to the BC5D table behaviour.
   The scalar BC in `inputs.bc_value` is left untouched; the engine uses it automatically
   for any velocity gap not covered by a segment.

4. **Mutual exclusion with `drag_table` and `bc5d_table_path`** — all three features
   drive the same `bc_segments_data` / `custom_drag_table` slots in `BallisticInputs`;
   combining any two would silently apply whichever ran last. Hard errors at validation
   time prevent that.

5. **Zero-search correctness guaranteed by design** — `prepare_request()` is called once
   and its `BallisticInputs` is reused for both the elevation search and the trajectory.
   No separate injection step is needed for the zero pass. This is the structural advantage
   over a `run_pass_direct()` approach that would require explicit per-pass injection.

