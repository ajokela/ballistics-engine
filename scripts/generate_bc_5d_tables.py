#!/usr/bin/env python3
"""
5D BC Correction Table Generator

Queries the production API to capture ML model predictions and stores them
as compact binary lookup tables for offline use.

Dimensions:
  1. Weight (caliber-specific ranges)
  2. Base BC (0.05 - 1.2)
  3. Muzzle Velocity (2000 - 4000 fps)
  4. Current Velocity (500 - 4000 fps)
  5. Drag Model (G1, G7)

Usage:
  python generate_bc_5d_tables.py --caliber 0.308 --output data/bc_tables/
  python generate_bc_5d_tables.py --all --output data/bc_tables/
"""

import argparse
import json
import os
import struct
import sys
import time
import zlib
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import requests

# Configuration
API_BASE_URL = "https://api.ballistics.7.62x51mm.sh"
BC_SEGMENTS_ENDPOINT = "/v1/bc_segments/estimate"
RATE_LIMIT_DELAY = 0.67  # 1.5 req/sec = 0.67 seconds between requests (for production API)
LOCAL_RATE_LIMIT_DELAY = 0.01  # 100 req/sec for local server

# Binary format constants
MAGIC = b"BC5D"
VERSION = 2
HEADER_SIZE = 80

# Caliber configurations with weight ranges and bins
CALIBER_CONFIGS = {
    0.224: {"weight_range": (40, 90), "weight_bins": 14, "file": "bc5d_224.bin"},
    0.243: {"weight_range": (55, 115), "weight_bins": 15, "file": "bc5d_243.bin"},
    0.264: {"weight_range": (95, 160), "weight_bins": 16, "file": "bc5d_264.bin"},
    0.277: {"weight_range": (115, 170), "weight_bins": 14, "file": "bc5d_277.bin"},
    0.284: {"weight_range": (140, 200), "weight_bins": 15, "file": "bc5d_284.bin"},
    0.308: {"weight_range": (110, 230), "weight_bins": 20, "file": "bc5d_308.bin"},
    0.338: {"weight_range": (200, 300), "weight_bins": 15, "file": "bc5d_338.bin"},
}

# Fixed dimension bins
BC_BINS = 24  # 0.05 to 1.2
MUZZLE_VEL_BINS = 11  # 2000 to 4000 fps
CURRENT_VEL_BINS = 36  # 500 to 4000 fps (dense in transonic)
DRAG_TYPES = 2  # G1, G7


@dataclass
class BC5DHeader:
    """Header structure for BC5D binary format."""

    magic: bytes
    version: int
    caliber: float
    flags: int
    dim_weight: int
    dim_bc: int
    dim_muzzle_vel: int
    dim_current_vel: int
    dim_drag_types: int
    timestamp: int
    checksum: int
    api_version: str
    reserved: bytes

    def to_bytes(self) -> bytes:
        """Serialize header to bytes."""
        api_version_bytes = self.api_version.encode("utf-8")[:16].ljust(16, b"\x00")
        return struct.pack(
            "<4sIfIIIIIIIQI16s12s",
            self.magic,
            self.version,
            self.caliber,
            self.flags,
            0,  # padding
            self.dim_weight,
            self.dim_bc,
            self.dim_muzzle_vel,
            self.dim_current_vel,
            self.dim_drag_types,
            self.timestamp,
            self.checksum,
            api_version_bytes,
            self.reserved,
        )

    @classmethod
    def from_bytes(cls, data: bytes) -> "BC5DHeader":
        """Deserialize header from bytes."""
        unpacked = struct.unpack("<4sIfIIIIIIIQI16s12s", data[:HEADER_SIZE])
        return cls(
            magic=unpacked[0],
            version=unpacked[1],
            caliber=unpacked[2],
            flags=unpacked[3],
            dim_weight=unpacked[5],
            dim_bc=unpacked[6],
            dim_muzzle_vel=unpacked[7],
            dim_current_vel=unpacked[8],
            dim_drag_types=unpacked[9],
            timestamp=unpacked[10],
            checksum=unpacked[11],
            api_version=unpacked[12].rstrip(b"\x00").decode("utf-8"),
            reserved=unpacked[13],
        )


def generate_weight_bins(config: dict) -> list[float]:
    """Generate weight bins for a caliber."""
    min_w, max_w = config["weight_range"]
    n_bins = config["weight_bins"]
    step = (max_w - min_w) / (n_bins - 1)
    return [min_w + i * step for i in range(n_bins)]


def generate_bc_bins() -> list[float]:
    """Generate BC bins (0.05 to 1.2)."""
    return [0.05 + i * (1.15 / (BC_BINS - 1)) for i in range(BC_BINS)]


def generate_muzzle_vel_bins() -> list[float]:
    """Generate muzzle velocity bins (2000 to 4000 fps)."""
    return [2000 + i * (2000 / (MUZZLE_VEL_BINS - 1)) for i in range(MUZZLE_VEL_BINS)]


def generate_current_vel_bins() -> list[float]:
    """
    Generate current velocity bins (500 to 4000 fps).
    Dense sampling in transonic region (900-1400 fps).
    """
    bins = []

    # Subsonic region: 500-900 fps (sparse)
    for v in range(500, 900, 100):
        bins.append(float(v))

    # Transonic region: 900-1400 fps (dense - every 25 fps)
    for v in range(900, 1400, 25):
        bins.append(float(v))

    # Supersonic region: 1400-4000 fps (moderate)
    for v in range(1400, 4050, 100):
        bins.append(float(v))

    return sorted(set(bins))[:CURRENT_VEL_BINS]


def query_api_bc_segments(
    caliber: float,
    weight: float,
    bc_value: float,
    bc_type: str,
    muzzle_velocity: float,
    session: requests.Session,
    api_url: str = API_BASE_URL,
) -> Optional[list[dict]]:
    """
    Query the API to get BC segments for a given configuration.

    Returns list of segments: [{"velocity_min": f, "velocity_max": f, "bc_value": f}, ...]
    """
    payload = {
        "manufacturer": "Generic",  # Required field - use placeholder
        "model": f"{weight:.0f}gr Match",  # Required field - descriptive placeholder
        "caliber": caliber,
        "weight": weight,
        "bc_value": bc_value,
        "bc_type": bc_type,
        "muzzle_velocity": muzzle_velocity,
        "bullet_type": "match",  # Use match for most accurate predictions
    }

    try:
        resp = session.post(
            f"{api_url}{BC_SEGMENTS_ENDPOINT}",
            json=payload,
            timeout=10,
        )
        if resp.status_code == 200:
            data = resp.json()
            return data.get("segments", [])
        else:
            print(f"API error {resp.status_code}: {resp.text[:100]}", file=sys.stderr)
            return None
    except requests.RequestException as e:
        print(f"Request failed: {e}", file=sys.stderr)
        return None


def get_bc_at_velocity(segments: list[dict], velocity: float) -> float:
    """
    Get BC at a specific velocity from segments.
    Uses linear interpolation between segment boundaries.
    """
    if not segments:
        return 1.0  # No correction

    # Sort segments by velocity_min
    sorted_segs = sorted(segments, key=lambda s: s["velocity_min"])

    # Find containing segment
    for seg in sorted_segs:
        if seg["velocity_min"] <= velocity <= seg["velocity_max"]:
            return seg["bc_value"]

    # Below all segments - use lowest
    if velocity < sorted_segs[0]["velocity_min"]:
        return sorted_segs[0]["bc_value"]

    # Above all segments - use highest
    if velocity > sorted_segs[-1]["velocity_max"]:
        return sorted_segs[-1]["bc_value"]

    # Shouldn't reach here, but fallback
    return sorted_segs[-1]["bc_value"]


def generate_caliber_table(
    caliber: float,
    config: dict,
    output_dir: Path,
    verbose: bool = False,
    api_url: str = API_BASE_URL,
) -> bool:
    """
    Generate the 5D BC correction table for a single caliber.

    Returns True on success, False on failure.
    """
    print(f"\n{'=' * 60}")
    print(f"Generating table for caliber {caliber}")
    print(f"{'=' * 60}")

    # Generate bin arrays
    weight_bins = generate_weight_bins(config)
    bc_bins = generate_bc_bins()
    muzzle_vel_bins = generate_muzzle_vel_bins()
    current_vel_bins = generate_current_vel_bins()
    drag_types = ["G1", "G7"]

    print(f"  Weight bins: {len(weight_bins)} ({weight_bins[0]:.0f} - {weight_bins[-1]:.0f} gr)")
    print(f"  BC bins: {len(bc_bins)} ({bc_bins[0]:.3f} - {bc_bins[-1]:.3f})")
    print(f"  Muzzle vel bins: {len(muzzle_vel_bins)} ({muzzle_vel_bins[0]:.0f} - {muzzle_vel_bins[-1]:.0f} fps)")
    print(f"  Current vel bins: {len(current_vel_bins)} ({current_vel_bins[0]:.0f} - {current_vel_bins[-1]:.0f} fps)")
    print(f"  Drag types: {drag_types}")

    # Calculate total API queries needed
    total_queries = len(drag_types) * len(weight_bins) * len(bc_bins) * len(muzzle_vel_bins)
    print(f"\n  Total API queries: {total_queries:,}")
    print(f"  Estimated time: {total_queries * RATE_LIMIT_DELAY / 60:.1f} minutes")

    # Initialize correction array
    # Shape: [drag_type][weight][bc][muzzle_vel][current_vel]
    corrections = []
    for _ in range(len(drag_types)):
        drag_arr = []
        for _ in range(len(weight_bins)):
            weight_arr = []
            for _ in range(len(bc_bins)):
                bc_arr = []
                for _ in range(len(muzzle_vel_bins)):
                    bc_arr.append([1.0] * len(current_vel_bins))
                weight_arr.append(bc_arr)
            drag_arr.append(weight_arr)
        corrections.append(drag_arr)

    # Query API and populate corrections
    session = requests.Session()
    query_count = 0
    start_time = time.time()

    for drag_idx, drag_type in enumerate(drag_types):
        for weight_idx, weight in enumerate(weight_bins):
            for bc_idx, base_bc in enumerate(bc_bins):
                for muzzle_idx, muzzle_vel in enumerate(muzzle_vel_bins):
                    query_count += 1

                    if verbose or query_count % 100 == 0:
                        elapsed = time.time() - start_time
                        rate = query_count / elapsed if elapsed > 0 else 0
                        eta = (total_queries - query_count) / rate if rate > 0 else 0
                        print(
                            f"\r  Progress: {query_count}/{total_queries} "
                            f"({100 * query_count / total_queries:.1f}%) "
                            f"Rate: {rate:.1f}/s ETA: {eta / 60:.1f}m",
                            end="",
                            flush=True,
                        )

                    # Query API
                    segments = query_api_bc_segments(
                        caliber=caliber,
                        weight=weight,
                        bc_value=base_bc,
                        bc_type=drag_type,
                        muzzle_velocity=muzzle_vel,
                        session=session,
                        api_url=api_url,
                    )

                    if segments:
                        # Extract BC at each current velocity and compute correction
                        for vel_idx, current_vel in enumerate(current_vel_bins):
                            bc_at_vel = get_bc_at_velocity(segments, current_vel)
                            correction = bc_at_vel / base_bc if base_bc > 0 else 1.0
                            # Clip to reasonable range
                            correction = max(0.5, min(1.5, correction))
                            corrections[drag_idx][weight_idx][bc_idx][muzzle_idx][vel_idx] = correction

                    # Rate limiting (skip for local server)
                    if "localhost" in api_url or "127.0.0.1" in api_url:
                        time.sleep(LOCAL_RATE_LIMIT_DELAY)
                    else:
                        time.sleep(RATE_LIMIT_DELAY)

    print(f"\n  Completed {query_count} queries in {time.time() - start_time:.1f}s")

    # Write binary file
    output_path = output_dir / config["file"]
    write_binary_table(
        output_path=output_path,
        caliber=caliber,
        weight_bins=weight_bins,
        bc_bins=bc_bins,
        muzzle_vel_bins=muzzle_vel_bins,
        current_vel_bins=current_vel_bins,
        corrections=corrections,
    )

    print(f"  Written: {output_path} ({output_path.stat().st_size:,} bytes)")
    return True


def write_binary_table(
    output_path: Path,
    caliber: float,
    weight_bins: list[float],
    bc_bins: list[float],
    muzzle_vel_bins: list[float],
    current_vel_bins: list[float],
    corrections: list,
) -> None:
    """Write the 5D correction table to binary format."""
    import time as time_module

    # Prepare data buffer
    data_buffer = bytearray()

    # Write bins
    for w in weight_bins:
        data_buffer.extend(struct.pack("<f", w))
    for b in bc_bins:
        data_buffer.extend(struct.pack("<f", b))
    for m in muzzle_vel_bins:
        data_buffer.extend(struct.pack("<f", m))
    for c in current_vel_bins:
        data_buffer.extend(struct.pack("<f", c))

    # Write corrections as flat array
    # Order: [drag_type][weight][bc][muzzle_vel][current_vel]
    for drag_arr in corrections:
        for weight_arr in drag_arr:
            for bc_arr in weight_arr:
                for muzzle_arr in bc_arr:
                    for correction in muzzle_arr:
                        data_buffer.extend(struct.pack("<f", correction))

    # Calculate checksum
    checksum = zlib.crc32(data_buffer) & 0xFFFFFFFF

    # Create header
    header = BC5DHeader(
        magic=MAGIC,
        version=VERSION,
        caliber=caliber,
        flags=0,
        dim_weight=len(weight_bins),
        dim_bc=len(bc_bins),
        dim_muzzle_vel=len(muzzle_vel_bins),
        dim_current_vel=len(current_vel_bins),
        dim_drag_types=DRAG_TYPES,
        timestamp=int(time_module.time()),
        checksum=checksum,
        api_version="v0.31.x",
        reserved=b"\x00" * 12,
    )

    # Write file
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "wb") as f:
        f.write(header.to_bytes())
        f.write(data_buffer)


def verify_table(table_path: Path) -> bool:
    """Verify a generated table file."""
    print(f"\nVerifying {table_path}...")

    with open(table_path, "rb") as f:
        header_bytes = f.read(HEADER_SIZE)
        header = BC5DHeader.from_bytes(header_bytes)

        if header.magic != MAGIC:
            print(f"  ERROR: Invalid magic bytes: {header.magic}")
            return False

        print(f"  Magic: {header.magic}")
        print(f"  Version: {header.version}")
        print(f"  Caliber: {header.caliber}")
        print(f"  Dimensions: {header.dim_weight} x {header.dim_bc} x {header.dim_muzzle_vel} x {header.dim_current_vel} x {header.dim_drag_types}")
        print(f"  API Version: {header.api_version}")

        # Read and verify checksum
        data = f.read()
        calculated_checksum = zlib.crc32(data) & 0xFFFFFFFF
        if calculated_checksum != header.checksum:
            print(f"  ERROR: Checksum mismatch! Expected {header.checksum}, got {calculated_checksum}")
            return False
        print(f"  Checksum: OK ({header.checksum:08x})")

    return True


def main():
    parser = argparse.ArgumentParser(description="Generate 5D BC correction tables from API")
    parser.add_argument("--caliber", type=float, help="Single caliber to generate (e.g., 0.308)")
    parser.add_argument("--all", action="store_true", help="Generate all calibers")
    parser.add_argument("--output", type=str, default="data/bc_tables", help="Output directory")
    parser.add_argument("--verify", type=str, help="Verify an existing table file")
    parser.add_argument("--verbose", action="store_true", help="Verbose output")
    parser.add_argument("--api-url", type=str, help=f"API base URL (default: {API_BASE_URL})")

    args = parser.parse_args()

    # Update API_BASE_URL if custom URL provided
    api_url = args.api_url if args.api_url else API_BASE_URL

    output_dir = Path(args.output)

    if args.verify:
        verify_table(Path(args.verify))
        return

    if args.all:
        calibers = list(CALIBER_CONFIGS.keys())
    elif args.caliber:
        if args.caliber not in CALIBER_CONFIGS:
            print(f"Unknown caliber {args.caliber}. Available: {list(CALIBER_CONFIGS.keys())}")
            sys.exit(1)
        calibers = [args.caliber]
    else:
        parser.print_help()
        sys.exit(1)

    print(f"BC5D Table Generator")
    print(f"API: {api_url}")
    print(f"Output: {output_dir}")
    print(f"Calibers: {calibers}")

    for caliber in calibers:
        config = CALIBER_CONFIGS[caliber]
        success = generate_caliber_table(
            caliber=caliber,
            config=config,
            output_dir=output_dir,
            verbose=args.verbose,
            api_url=api_url,
        )
        if success:
            verify_table(output_dir / config["file"])


if __name__ == "__main__":
    main()
