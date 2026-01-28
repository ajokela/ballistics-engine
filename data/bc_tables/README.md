# BC5D Correction Tables

## Overview

This directory contains BC5D (5-Dimensional BC Correction) tables that provide ML-derived, velocity-dependent BC corrections for specific calibers. These tables capture how ballistic coefficient changes throughout the flight envelope.

## Available Tables

| File | Caliber | Size |
|------|---------|------|
| bc5d_224.bin | .224 (5.56mm) | 1.0 MB |
| bc5d_243.bin | .243 (6mm) | 1.1 MB |
| bc5d_264.bin | .264 (6.5mm) | 1.2 MB |
| bc5d_277.bin | .277 (6.8mm) | 1.0 MB |
| bc5d_284.bin | .284 (7mm) | 1.1 MB |
| bc5d_308.bin | .308 (7.62mm) | 1.5 MB |
| bc5d_338.bin | .338 (8.6mm) | 1.1 MB |

## Binary Format (v2)

### Header (80 bytes)

| Offset | Size | Type | Field |
|--------|------|------|-------|
| 0 | 4 | bytes | Magic: `BC5D` |
| 4 | 4 | uint32 | Version (2) |
| 8 | 4 | float32 | Caliber (inches) |
| 12 | 4 | uint32 | Flags |
| 16 | 4 | uint32 | Padding |
| 20 | 4 | uint32 | dim_weight |
| 24 | 4 | uint32 | dim_bc |
| 28 | 4 | uint32 | dim_muzzle_vel |
| 32 | 4 | uint32 | dim_current_vel |
| 36 | 4 | uint32 | dim_drag_types |
| 40 | 8 | uint64 | Timestamp (Unix) |
| 48 | 4 | uint32 | CRC32 (data section) |
| 52 | 16 | string | API version (null-padded) |
| 68 | 12 | bytes | Reserved |

### Bin Definitions

Following the header:
- Weight bins: `dim_weight × 4 bytes` (float32)
- BC bins: `dim_bc × 4 bytes` (float32)
- Muzzle velocity bins: `dim_muzzle_vel × 4 bytes` (float32)
- Current velocity bins: `dim_current_vel × 4 bytes` (float32)

### Data Section

Correction factors as float32 values:
- Layout: `[drag_type][weight][bc][muzzle_vel][current_vel]`
- Total cells: `drag_types × weight × bc × muzzle_vel × current_vel`
- Values represent correction ratios (typically 0.5 to 1.5)

## Usage

### CLI

```bash
# Use tables from this directory
./ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --bc-table-dir ./data/bc_tables/

# Auto-download tables (caches locally)
./ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --bc-table-auto
```

### Library

```rust
use ballistics_engine::bc_table_5d::Bc5dTableManager;

let mut manager = Bc5dTableManager::new("./data/bc_tables/");
let table = manager.get_table(0.308)?;

// Get effective BC at current velocity
let effective_bc = table.get_effective_bc(
    168.0,   // weight (grains)
    0.475,   // base BC
    2700.0,  // muzzle velocity (fps)
    2000.0,  // current velocity (fps)
    "G1",    // drag model
);
```

## Online Distribution

These tables are also available for download at:
- URL: `https://ballistics.tools/downloads/bc5d/`
- Manifest: `https://ballistics.tools/downloads/bc5d/manifest.json`

The CLI can auto-download tables using `--bc-table-auto`.

## Generation

Tables are generated from ML models trained on doppler-derived ballistics data. The generation process:

1. Sample the 5D parameter space (weight, BC, muzzle velocity, current velocity, drag model)
2. Query the ML model for each sample point
3. Store correction factors as a binary lookup table
4. Apply CRC32 validation for data integrity

## License

Tables are free to download and use. The ML models used to generate them are proprietary.
