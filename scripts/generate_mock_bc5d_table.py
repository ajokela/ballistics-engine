#!/usr/bin/env python3
"""
Generate a mock BC5D table for testing the Rust loader.

This creates a small table with known values to verify the binary format
is correctly parsed by the Rust code.
"""

import struct
import zlib
from pathlib import Path

# Configuration
MAGIC = b"BC5D"
VERSION = 2
OUTPUT_FILE = Path("data/bc_tables/bc5d_308.bin")

# Small dimensions for testing
WEIGHT_BINS = [150.0, 175.0, 200.0]
BC_BINS = [0.4, 0.5, 0.6]
MUZZLE_VEL_BINS = [2700.0, 3000.0]
CURRENT_VEL_BINS = [1000.0, 1500.0, 2000.0, 2500.0]
DRAG_TYPES = 2  # G1, G7


def main():
    print(f"Generating mock BC5D table...")
    print(f"  Weight bins: {WEIGHT_BINS}")
    print(f"  BC bins: {BC_BINS}")
    print(f"  Muzzle vel bins: {MUZZLE_VEL_BINS}")
    print(f"  Current vel bins: {CURRENT_VEL_BINS}")
    print(f"  Drag types: {DRAG_TYPES}")

    # Calculate total cells
    total_cells = (
        DRAG_TYPES
        * len(WEIGHT_BINS)
        * len(BC_BINS)
        * len(MUZZLE_VEL_BINS)
        * len(CURRENT_VEL_BINS)
    )
    print(f"  Total cells: {total_cells}")

    # Generate correction values with a pattern for verification
    # Correction = 1.0 + 0.01 * (velocity_idx) - 0.05 * (bc_idx > 1)
    # This creates predictable values we can verify in tests
    corrections = []
    for drag_idx in range(DRAG_TYPES):
        for weight_idx in range(len(WEIGHT_BINS)):
            for bc_idx in range(len(BC_BINS)):
                for muzzle_idx in range(len(MUZZLE_VEL_BINS)):
                    for vel_idx in range(len(CURRENT_VEL_BINS)):
                        # Simple pattern: varies with velocity index
                        # At low velocity (idx 0) = 0.95 (transonic degradation)
                        # At high velocity (idx 3) = 1.0 (no correction)
                        correction = 0.95 + 0.0125 * vel_idx
                        if drag_idx == 1:  # G7 has slightly less degradation
                            correction += 0.02
                        corrections.append(correction)

    # Prepare data buffer (bins + corrections)
    data_buffer = bytearray()

    # Write bins
    for w in WEIGHT_BINS:
        data_buffer.extend(struct.pack("<f", w))
    for b in BC_BINS:
        data_buffer.extend(struct.pack("<f", b))
    for m in MUZZLE_VEL_BINS:
        data_buffer.extend(struct.pack("<f", m))
    for c in CURRENT_VEL_BINS:
        data_buffer.extend(struct.pack("<f", c))

    # Write corrections
    for corr in corrections:
        data_buffer.extend(struct.pack("<f", corr))

    # Calculate checksum
    checksum = zlib.crc32(data_buffer) & 0xFFFFFFFF

    # Create header
    api_version = b"mock_v1.0\x00\x00\x00\x00\x00\x00\x00"[:16]
    reserved = b"\x00" * 12
    timestamp = 1706400000  # Fixed timestamp for reproducibility

    header = struct.pack(
        "<4sIfIIIIIIIQI16s12s",
        MAGIC,
        VERSION,
        0.308,  # caliber
        0,  # flags
        0,  # padding
        len(WEIGHT_BINS),
        len(BC_BINS),
        len(MUZZLE_VEL_BINS),
        len(CURRENT_VEL_BINS),
        DRAG_TYPES,
        timestamp,
        checksum,
        api_version,
        reserved,
    )

    # Ensure output directory exists
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)

    # Write file
    with open(OUTPUT_FILE, "wb") as f:
        f.write(header)
        f.write(data_buffer)

    print(f"\nWritten: {OUTPUT_FILE} ({OUTPUT_FILE.stat().st_size} bytes)")
    print(f"Checksum: {checksum:08x}")

    # Verify by reading back
    print("\nVerification:")
    with open(OUTPUT_FILE, "rb") as f:
        magic = f.read(4)
        print(f"  Magic: {magic}")
        assert magic == MAGIC, f"Magic mismatch: {magic}"

        version = struct.unpack("<I", f.read(4))[0]
        print(f"  Version: {version}")
        assert version == VERSION, f"Version mismatch: {version}"

        caliber = struct.unpack("<f", f.read(4))[0]
        print(f"  Caliber: {caliber}")
        assert abs(caliber - 0.308) < 0.001, f"Caliber mismatch: {caliber}"

    print("\nMock table generated successfully!")
    print(f"\nTest with: cargo run -- trajectory -v 3000 -b 0.5 -m 175 -d 0.308 --bc-table-dir data/bc_tables --max-range 1000")


if __name__ == "__main__":
    main()
