#!/usr/bin/env python3
"""
Benchmark BC5D offline tables versus online API predictions.

Compares the BC correction factors from the precomputed tables
against fresh API queries to validate accuracy.
"""

import struct
import sys
import requests
from pathlib import Path

API_URL = "http://127.0.0.1:5001"
BC_SEGMENTS_ENDPOINT = "/v1/bc_segments/estimate"
TABLE_DIR = Path(__file__).parent.parent / "data/bc_tables"

# Test configurations: (caliber, weight, bc, muzzle_vel, drag_type)
TEST_CASES = [
    # .224 (5.56mm)
    (0.224, 55, 0.243, 3200, "G1"),
    (0.224, 77, 0.372, 2750, "G1"),
    # .243 (6mm)
    (0.243, 95, 0.430, 3000, "G1"),
    (0.243, 105, 0.530, 2900, "G7"),
    # .264 (6.5mm)
    (0.264, 140, 0.625, 2900, "G1"),
    (0.264, 147, 0.697, 2700, "G7"),
    # .277 (270)
    (0.277, 130, 0.450, 3050, "G1"),
    (0.277, 150, 0.574, 2850, "G1"),
    # .284 (7mm)
    (0.284, 162, 0.625, 2900, "G1"),
    (0.284, 175, 0.672, 2800, "G7"),
    # .308 (30 cal)
    (0.308, 168, 0.462, 2700, "G1"),
    (0.308, 175, 0.505, 2800, "G1"),
    (0.308, 190, 0.533, 2650, "G7"),
    # .338
    (0.338, 250, 0.670, 2800, "G1"),
    (0.338, 300, 0.818, 2650, "G7"),
]

# Velocity test points (fps)
VELOCITY_POINTS = [3000, 2500, 2000, 1500, 1200, 1000, 800]


def load_bc5d_table(caliber: float) -> dict:
    """Load a BC5D table and return bin/data info."""
    caliber_int = int(caliber * 1000)
    table_path = TABLE_DIR / f"bc5d_{caliber_int}.bin"

    if not table_path.exists():
        return None

    with open(table_path, "rb") as f:
        # Read header
        magic = f.read(4)
        if magic != b"BC5D":
            return None

        version = struct.unpack("<I", f.read(4))[0]
        cal = struct.unpack("<f", f.read(4))[0]
        flags = struct.unpack("<I", f.read(4))[0]
        padding = struct.unpack("<I", f.read(4))[0]

        dim_weight = struct.unpack("<I", f.read(4))[0]
        dim_bc = struct.unpack("<I", f.read(4))[0]
        dim_muzzle = struct.unpack("<I", f.read(4))[0]
        dim_current = struct.unpack("<I", f.read(4))[0]
        dim_drag = struct.unpack("<I", f.read(4))[0]

        timestamp = struct.unpack("<Q", f.read(8))[0]
        checksum = struct.unpack("<I", f.read(4))[0]
        api_version = f.read(16).rstrip(b"\x00").decode()
        reserved = f.read(12)

        # Read bins
        weight_bins = [struct.unpack("<f", f.read(4))[0] for _ in range(dim_weight)]
        bc_bins = [struct.unpack("<f", f.read(4))[0] for _ in range(dim_bc)]
        muzzle_bins = [struct.unpack("<f", f.read(4))[0] for _ in range(dim_muzzle)]
        current_bins = [struct.unpack("<f", f.read(4))[0] for _ in range(dim_current)]

        # Read data
        total = dim_drag * dim_weight * dim_bc * dim_muzzle * dim_current
        data = [struct.unpack("<f", f.read(4))[0] for _ in range(total)]

    return {
        "weight_bins": weight_bins,
        "bc_bins": bc_bins,
        "muzzle_bins": muzzle_bins,
        "current_bins": current_bins,
        "dim_drag": dim_drag,
        "data": data,
    }


def interp_idx(value: float, bins: list) -> tuple:
    """Find interpolation index and weight."""
    if value <= bins[0]:
        return 0, 0.0
    if value >= bins[-1]:
        return len(bins) - 2, 1.0

    for i in range(len(bins) - 1):
        if bins[i] <= value <= bins[i + 1]:
            w = (value - bins[i]) / (bins[i + 1] - bins[i])
            return i, w
    return len(bins) - 2, 1.0


def lookup_bc5d(table: dict, weight: float, bc: float, muzzle_vel: float, current_vel: float, drag_type: str) -> float:
    """Look up correction factor from BC5D table with 4D interpolation."""
    drag_idx = 1 if drag_type.upper() == "G7" else 0

    wi, ww = interp_idx(weight, table["weight_bins"])
    bi, bw = interp_idx(bc, table["bc_bins"])
    mi, mw = interp_idx(muzzle_vel, table["muzzle_bins"])
    ci, cw = interp_idx(current_vel, table["current_bins"])

    n_weight = len(table["weight_bins"])
    n_bc = len(table["bc_bins"])
    n_muzzle = len(table["muzzle_bins"])
    n_current = len(table["current_bins"])

    def flat_idx(d, w, b, m, c):
        return (d * n_weight * n_bc * n_muzzle * n_current +
                w * n_bc * n_muzzle * n_current +
                b * n_muzzle * n_current +
                m * n_current +
                c)

    result = 0.0
    for dw in [0, 1]:
        for db in [0, 1]:
            for dm in [0, 1]:
                for dc in [0, 1]:
                    weight_factor = (1 - ww) if dw == 0 else ww
                    bc_factor = (1 - bw) if db == 0 else bw
                    muzzle_factor = (1 - mw) if dm == 0 else mw
                    current_factor = (1 - cw) if dc == 0 else cw

                    w_idx = min(wi + dw, n_weight - 1)
                    b_idx = min(bi + db, n_bc - 1)
                    m_idx = min(mi + dm, n_muzzle - 1)
                    c_idx = min(ci + dc, n_current - 1)

                    idx = flat_idx(drag_idx, w_idx, b_idx, m_idx, c_idx)
                    result += weight_factor * bc_factor * muzzle_factor * current_factor * table["data"][idx]

    return max(0.5, min(1.5, result))


def query_api_bc(caliber: float, weight: float, bc: float, muzzle_vel: float, drag_type: str) -> list:
    """Query API for BC segments."""
    payload = {
        "manufacturer": "Generic",
        "model": f"{weight:.0f}gr Match",
        "caliber": caliber,
        "weight": weight,
        "bc_value": bc,
        "bc_type": drag_type,
        "muzzle_velocity": muzzle_vel,
        "bullet_type": "match",
    }

    try:
        resp = requests.post(f"{API_URL}{BC_SEGMENTS_ENDPOINT}", json=payload, timeout=10)
        if resp.status_code == 200:
            return resp.json().get("segments", [])
    except Exception as e:
        print(f"API error: {e}")
    return []


def get_api_bc_at_velocity(segments: list, velocity: float) -> float:
    """Get BC at velocity from API segments."""
    if not segments:
        return None

    for seg in sorted(segments, key=lambda s: s["velocity_min"]):
        if seg["velocity_min"] <= velocity <= seg["velocity_max"]:
            return seg["bc_value"]

    # Fallback to nearest
    if velocity < segments[0]["velocity_min"]:
        return min(segments, key=lambda s: s["velocity_min"])["bc_value"]
    return max(segments, key=lambda s: s["velocity_max"])["bc_value"]


def main():
    print("=" * 80)
    print("BC5D Offline vs Online API Benchmark")
    print("=" * 80)
    print()

    total_comparisons = 0
    total_error = 0.0
    max_error = 0.0
    max_error_case = None
    errors_above_1pct = 0
    errors_above_2pct = 0
    errors_above_5pct = 0

    tables = {}

    for caliber, weight, bc, muzzle_vel, drag_type in TEST_CASES:
        # Load table if not cached
        if caliber not in tables:
            tables[caliber] = load_bc5d_table(caliber)
            if tables[caliber] is None:
                print(f"Warning: No table for caliber {caliber}")
                continue

        table = tables[caliber]
        if table is None:
            continue

        # Query API
        segments = query_api_bc(caliber, weight, bc, muzzle_vel, drag_type)
        if not segments:
            print(f"Warning: No API response for {caliber}/{weight}gr/{bc}/{drag_type}")
            continue

        print(f"\n{caliber} cal, {weight}gr, BC={bc} {drag_type}, MV={muzzle_vel}fps:")
        print(f"  {'Velocity':>8} {'API BC':>10} {'Table BC':>10} {'Diff':>8} {'Error%':>8}")
        print(f"  {'-'*8} {'-'*10} {'-'*10} {'-'*8} {'-'*8}")

        for vel in VELOCITY_POINTS:
            if vel > muzzle_vel:
                continue

            # Get API BC at this velocity
            api_bc = get_api_bc_at_velocity(segments, vel)
            if api_bc is None:
                continue

            # Get table BC at this velocity
            correction = lookup_bc5d(table, weight, bc, muzzle_vel, vel, drag_type)
            table_bc = bc * correction

            # Calculate error
            diff = table_bc - api_bc
            error_pct = abs(diff / api_bc) * 100 if api_bc > 0 else 0

            total_comparisons += 1
            total_error += error_pct

            if error_pct > max_error:
                max_error = error_pct
                max_error_case = (caliber, weight, bc, drag_type, muzzle_vel, vel, api_bc, table_bc)

            if error_pct > 1:
                errors_above_1pct += 1
            if error_pct > 2:
                errors_above_2pct += 1
            if error_pct > 5:
                errors_above_5pct += 1

            status = "OK" if error_pct < 2 else "!!" if error_pct < 5 else "XXX"
            print(f"  {vel:>8} {api_bc:>10.5f} {table_bc:>10.5f} {diff:>+8.5f} {error_pct:>7.2f}% {status}")

    # Summary
    print()
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Total comparisons: {total_comparisons}")
    print(f"Average error: {total_error / total_comparisons:.3f}%" if total_comparisons > 0 else "N/A")
    print(f"Max error: {max_error:.3f}%")
    if max_error_case:
        cal, wt, bc, dt, mv, vel, api, tbl = max_error_case
        print(f"  at: {cal} cal, {wt}gr, BC={bc} {dt}, MV={mv}, vel={vel}")
        print(f"      API={api:.5f}, Table={tbl:.5f}")
    print()
    print(f"Errors > 1%: {errors_above_1pct} ({100*errors_above_1pct/total_comparisons:.1f}%)" if total_comparisons > 0 else "N/A")
    print(f"Errors > 2%: {errors_above_2pct} ({100*errors_above_2pct/total_comparisons:.1f}%)" if total_comparisons > 0 else "N/A")
    print(f"Errors > 5%: {errors_above_5pct} ({100*errors_above_5pct/total_comparisons:.1f}%)" if total_comparisons > 0 else "N/A")

    # Pass/fail
    print()
    if total_comparisons > 0 and max_error < 5 and errors_above_2pct / total_comparisons < 0.1:
        print("RESULT: PASS - Tables match API within acceptable tolerance")
        return 0
    else:
        print("RESULT: FAIL - Tables have significant deviation from API")
        return 1


if __name__ == "__main__":
    sys.exit(main())
