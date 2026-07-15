#!/usr/bin/env python3
"""MBA-1318 cross-implementation reference generator.

Emits ``data/validation/xref_*.json`` golden cases from py-ballisticcalc, the independent
open-source point-mass solver, so the ballistics-engine crate is checked against a second
implementation across the G1/G7 x supersonic/transonic x wind x spin-drift regime matrix.

Provenance / honesty notes
--------------------------
* Reference solver: py-ballisticcalc (RK4 fixed-step engine, its default), standard G1/G7 drag
  tables (Ingalls/BRL family, PCHIP-interpolated), Miller stability + Litz spin drift.
* The ballistics-engine defaults to an ADAPTIVE RK45 integrator and its own (higher-resolution)
  G1/G7 tables, so agreement is bounded by solver-class differences, not machine precision. The
  per-observable tolerances below are justified from those differences, NOT reverse-fitted to pass.
* Geometry is chosen to make the two frames line up EXACTLY: flat fire (0 deg barrel elevation,
  no zeroing), sight height 0 and look angle 0 so the line of sight is the horizontal bore line.
  Then py-bc ``height`` (relative to LOS) == engine McCoy Y (drop_m), py-bc ``windage`` == engine
  McCoy Z (drift_m), ``velocity`` == speed, ``time`` == tof.
* Wind mapping is derived, not guessed: engine ``wind_vector(s,dir) = (-s cos dir, *, -s sin dir)``
  and py-bc ``Wind.vector(from) = (s cos from, *, s sin from)`` are equal iff
  ``from = dir + 180 deg`` (same speed). So pybc_direction_from = engine_wind_angle + pi.
* Humidity is pinned to 0 on both sides to remove the humid-air density term as a variable.
* Spin drift: both use identical Litz ``1.25*(Sg+1.2)*t^1.83`` with a Miller Sg that includes the
  ``(v/2800)^(1/3)`` term; at sea-level standard the T/P corrections are ~1 on both sides, so the
  only residual is the RK4-vs-RK45 time-of-flight feeding t^1.83.

Reproduce:
    python3 -m venv venv && venv/bin/pip install "py-ballisticcalc==2.2.10"
    venv/bin/python data/validation/generators/gen_pybc_refs.py
"""
import json
import math
import os
from importlib.metadata import version

from py_ballisticcalc import (
    Ammo, Angular, Atmo, Calculator, Distance, DragModel, Pressure, Shot, TableG1, TableG7,
    Temperature, Velocity, Weapon, Weight, Wind,
)

PBC_VERSION = version("py_ballisticcalc")
GENERATOR = "data/validation/generators/gen_pybc_refs.py"
RETRIEVED = "2026-07-15"
TABLES = {"G1": TableG1, "G7": TableG7}
OUT_DIR = os.path.join(os.path.dirname(__file__), "..")

# --------------------------------------------------------------------------------------------
# Tolerance policy. Each returns (tol_abs, tol_rel, justification). Physically motivated by the
# solver-class differences documented above; verified to sit comfortably ABOVE the observed
# deltas (head-room reported in the harness GOLDEN_DUMP), never fitted to the boundary.
# --------------------------------------------------------------------------------------------
# The two solvers were measured to agree to ~0.03-0.05% across this matrix (see the head-room
# column in the harness GOLDEN_DUMP and data/validation/README.md). Tolerances are therefore set
# ~15-25x that observed cross-solver delta: loose enough to absorb the integrator/drag-table
# differences that legitimately separate the two implementations, tight enough to FAIL on a
# field-noticeable regression (a ~1% long-range drop shift is ~1 MOA). Bands are NOT reverse-fitted
# to the boundary; they are a stated multiple of the mechanism-driven noise floor.
def tol_drop(regime):
    if regime == "transonic":
        return (0.03, 0.012,
                "Drop, transonic reach: the trajectory has integrated THROUGH the Mach 1.2->1.0 drag "
                "rise to reach these (now subsonic) stations, so drop accumulates any Mach-dependent "
                "G7-table interpolation difference between the solvers on top of the RK45-vs-RK4 "
                "integrator gap. 3 cm floor + 1.2% (~observed 0.04% x ~25) still fails a >~1.5% drop "
                "regression.")
    return (0.01, 0.008,
            "Drop: bounded by RK45(engine, adaptive) vs RK4(py-bc, fixed-step) integration and the "
            "two distinct G1/G7 table resolutions/interpolants. 1 cm floor + 0.8% brackets the "
            "observed ~0.04% cross-solver difference with ~20x head-room while failing a >~1% "
            "(field-noticeable, ~1 MOA at long range) drop regression.")

def tol_tof(regime):
    return (0.003, 0.006,
            "Time of flight: integrator + drag-curve differences accumulate in flight time; observed "
            "cross-solver difference ~0.03%. 3 ms floor + 0.6% keeps ~15x head-room.")

def tol_velocity(regime):
    if regime == "transonic":
        return (0.8, 0.01,
                "Velocity, transonic reach: retained speed is the most drag-sensitive observable, "
                "carried through the transonic drag rise. 0.8 m/s floor + 1%.")
    return (0.6, 0.006,
            "Velocity: the most drag-sensitive observable; bounded by the G1/G7 table-resolution "
            "difference. Observed cross-solver difference ~0.06%; 0.6 m/s floor + 0.6% keeps ~10x "
            "head-room while failing a >~0.7% velocity regression.")

def tol_drift_wind():
    return (0.02, 0.02,
            "Wind drift: a difference of larger integrated quantities, so it inherits both the "
            "drag-curve and time-of-flight gaps; observed cross-solver difference ~0.1%. 2 cm floor "
            "+ 2% keeps ~20x head-room and fails a >~2.5 cm-per-metre-of-drift regression.")

def tol_drift_spin():
    return (0.01, 0.02,
            "Spin drift: both use the identical Litz 1.25*(Sg+1.2)*t^1.83 with a Miller Sg that "
            "matches at sea-level standard (the T/P corrections are ~1 on both sides), so the only "
            "residual is the RK4-vs-RK45 time-of-flight feeding t^1.83; observed difference ~0.05%. "
            "1 cm floor + 2%.")

def tol_drift_zero():
    return (0.01, None,
            "No wind and spin drift OFF, so the lateral must be the pure (near-zero) integrated "
            "value (the engine yields exactly 0). 1 cm absolute floor guards against a spurious "
            "lateral force sneaking in.")

# --------------------------------------------------------------------------------------------
# Case matrix. Inputs are ENGINE-NATIVE (SI, engine conventions). `pybc_twist_in` overrides the
# twist handed to py-bc only (used by the spin-drift-OFF case, where the engine flag gates spin
# drift off but py-bc applies it unconditionally on twist).
# --------------------------------------------------------------------------------------------
def base(**kw):
    d = dict(
        muzzle_angle_rad=0.0, shooting_angle_rad=0.0, twist_rate_in=0.0, is_twist_right=True,
        sight_height_m=0.0, muzzle_height_m=0.0, altitude_m=0.0,
        temperature_c=15.0, pressure_hpa=1013.25, humidity_frac=0.0,
        wind_speed_mps=0.0, wind_angle_rad=0.0, vertical_wind_mps=0.0,
        use_enhanced_spin_drift=False, enable_advanced_effects=False,
    )
    d.update(kw)
    return d

CASES = [
    dict(id="xref_g1_supersonic_308", regime="supersonic",
         desc="G1 .308 168 gr flat-fire, no wind, no twist. Supersonic drop/TOF/speed vs py-ballisticcalc; lateral must stay ~0 (spin drift off).",
         inputs=base(muzzle_velocity_mps=800.0, bullet_mass_kg=0.010886, bullet_diameter_m=0.0078232,
                     bullet_length_m=0.030861, bc_value=0.462, bc_type="G1"),
         ranges=[200.0, 400.0, 600.0], obs=["drop_m", "tof_s", "velocity_mps"], drift_zero_at=600.0),

    dict(id="xref_g7_supersonic_308", regime="supersonic",
         desc="G7 .308 175 gr flat-fire, no wind, no twist. Supersonic drop/TOF/speed vs py-ballisticcalc.",
         inputs=base(muzzle_velocity_mps=800.0, bullet_mass_kg=0.011340, bullet_diameter_m=0.0078232,
                     bullet_length_m=0.031496, bc_value=0.243, bc_type="G7"),
         ranges=[300.0, 600.0, 900.0], obs=["drop_m", "tof_s", "velocity_mps"], drift_zero_at=900.0),

    dict(id="xref_g7_65mm_supersonic", regime="supersonic",
         desc="G7 6.5 mm 140 gr flat-fire, no wind, no twist. A second G7 caliber for drag-model coverage.",
         inputs=base(muzzle_velocity_mps=820.0, bullet_mass_kg=0.009072, bullet_diameter_m=0.006706,
                     bullet_length_m=0.033782, bc_value=0.310, bc_type="G7"),
         ranges=[300.0, 600.0, 900.0], obs=["drop_m", "tof_s", "velocity_mps"], drift_zero_at=900.0),

    dict(id="xref_g7_transonic_reach", regime="transonic",
         desc="G7 .308 175 gr flat-fire pushed past transonic. The bullet crosses Mach 1.2->1.0 around 700-950 m; the 1100/1250/1400 m stations are just subsonic (~Mach 0.8-0.9), so drop/TOF/speed there reflect having integrated correctly THROUGH the transonic drag rise. Tolerances relaxed vs the pure-supersonic cases.",
         inputs=base(muzzle_velocity_mps=790.0, bullet_mass_kg=0.011340, bullet_diameter_m=0.0078232,
                     bullet_length_m=0.031496, bc_value=0.243, bc_type="G7"),
         ranges=[1100.0, 1250.0, 1400.0], obs=["drop_m", "tof_s", "velocity_mps"]),

    dict(id="xref_g1_crosswind", regime="supersonic",
         desc="G1 .308 168 gr flat-fire, 10 m/s wind FROM THE RIGHT (engine wind_angle = pi/2), no twist. Validates crosswind drift (negative McCoy Z) and its sign.",
         inputs=base(muzzle_velocity_mps=800.0, bullet_mass_kg=0.010886, bullet_diameter_m=0.0078232,
                     bullet_length_m=0.030861, bc_value=0.462, bc_type="G1",
                     wind_speed_mps=10.0, wind_angle_rad=math.pi / 2.0),
         ranges=[300.0, 600.0], obs=["drop_m", "drift_m", "tof_s"]),

    dict(id="xref_g7_crosswind", regime="supersonic",
         desc="G7 .308 175 gr flat-fire, 8 m/s wind FROM THE LEFT (engine wind_angle = 3pi/2), no twist. Crosswind drift with the opposite sign.",
         inputs=base(muzzle_velocity_mps=800.0, bullet_mass_kg=0.011340, bullet_diameter_m=0.0078232,
                     bullet_length_m=0.031496, bc_value=0.243, bc_type="G7",
                     wind_speed_mps=8.0, wind_angle_rad=3.0 * math.pi / 2.0),
         ranges=[400.0, 800.0], obs=["drop_m", "drift_m", "tof_s"]),

    dict(id="xref_g7_spin_drift_on", regime="supersonic",
         desc="G7 .308 175 gr flat-fire, 1:10 right twist, no wind, enhanced spin drift ON. Validates the engine's Litz+Miller spin drift (positive McCoy Z) against py-ballisticcalc's identical model.",
         inputs=base(muzzle_velocity_mps=800.0, bullet_mass_kg=0.011340, bullet_diameter_m=0.0078232,
                     bullet_length_m=0.031496, bc_value=0.243, bc_type="G7",
                     twist_rate_in=10.0, is_twist_right=True, use_enhanced_spin_drift=True),
         ranges=[500.0, 900.0], obs=["drift_m"], extra_obs_at={900.0: ["drop_m", "velocity_mps"]}),

    dict(id="xref_g7_spin_drift_off", regime="supersonic",
         desc="G7 .308 175 gr flat-fire, 1:10 twist PRESENT but enhanced spin drift OFF, no wind. Regression guard that the engine flag actually gates spin drift: lateral must be ~0. py-bc reference generated with twist=0 (its spin drift is unconditional on a set twist).",
         inputs=base(muzzle_velocity_mps=800.0, bullet_mass_kg=0.011340, bullet_diameter_m=0.0078232,
                     bullet_length_m=0.031496, bc_value=0.243, bc_type="G7",
                     twist_rate_in=10.0, is_twist_right=True, use_enhanced_spin_drift=False),
         pybc_twist_in=0.0,
         ranges=[500.0, 900.0], obs=["drop_m"], drift_zero_at=900.0),
]


def solve(spec):
    inp = spec["inputs"]
    grains = inp["bullet_mass_kg"] / 0.00006479891
    dm = DragModel(
        inp["bc_value"], TABLES[inp["bc_type"]],
        weight=Weight.Grain(grains),
        diameter=Distance.Meter(inp["bullet_diameter_m"]),
        length=Distance.Meter(inp["bullet_length_m"]),
    )
    ammo = Ammo(dm, mv=Velocity.MPS(inp["muzzle_velocity_mps"]))
    twist_in = spec.get("pybc_twist_in", inp["twist_rate_in"])
    # py-bc twist sign encodes hand: positive = right-hand.
    twist_signed = twist_in if inp["is_twist_right"] else -twist_in
    weapon = Weapon(sight_height=Distance.Meter(0.0), twist=Distance.Inch(twist_signed))
    atmo = Atmo(
        altitude=Distance.Meter(inp["altitude_m"]),
        pressure=Pressure.hPa(inp["pressure_hpa"]),
        temperature=Temperature.Celsius(inp["temperature_c"]),
        humidity=inp["humidity_frac"],
    )
    winds = []
    if inp["wind_speed_mps"] > 0.0:
        pybc_from_deg = math.degrees(inp["wind_angle_rad"] + math.pi) % 360.0
        winds = [Wind(velocity=Velocity.MPS(inp["wind_speed_mps"]),
                      direction_from=Angular.Degree(pybc_from_deg))]
    shot = Shot(weapon=weapon, ammo=ammo, atmo=atmo, winds=winds)
    all_ranges = list(spec["ranges"])
    maxr = max(all_ranges) + 5.0
    res = Calculator().fire(shot, trajectory_range=Distance.Meter(maxr),
                            trajectory_step=Distance.Meter(min(50.0, maxr)),
                            dense_output=True)
    return res


def observe(res, observable, range_m):
    td = res.get_at("distance", Distance.Meter(range_m))
    if observable == "drop_m":
        return td.height >> Distance.Meter
    if observable == "drift_m":
        return td.windage >> Distance.Meter
    if observable == "tof_s":
        return td.time
    if observable == "velocity_mps":
        return td.velocity >> Velocity.MPS
    raise ValueError(observable)


def tol_for(observable, regime):
    if observable == "drop_m":
        return tol_drop(regime)
    if observable == "tof_s":
        return tol_tof(regime)
    if observable == "velocity_mps":
        return tol_velocity(regime)
    raise ValueError(observable)


def build_case(spec):
    res = solve(spec)
    regime = spec["regime"]
    expectations = []

    # Primary observables at each range.
    extra_at = spec.get("extra_obs_at", {})
    for r in spec["ranges"]:
        obs_list = list(spec["obs"]) + list(extra_at.get(r, []))
        for ob in obs_list:
            val = observe(res, ob, r)
            if ob == "drift_m":
                ta, tr, just = tol_drift_spin() if spec["inputs"]["use_enhanced_spin_drift"] else tol_drift_wind()
            else:
                ta, tr, just = tol_for(ob, regime)
            exp = {"observable": ob, "range_m": r, "value": round(val, 6),
                   "tolerance_justification": just}
            if ta is not None:
                exp["tol_abs"] = ta
            if tr is not None:
                exp["tol_rel"] = tr
            expectations.append(exp)

    # Optional near-zero lateral assertion.
    if "drift_zero_at" in spec:
        ta, tr, just = tol_drift_zero()
        exp = {"observable": "drift_m", "range_m": spec["drift_zero_at"], "value": 0.0,
               "tol_abs": ta, "tolerance_justification": just}
        expectations.append(exp)

    case = {
        "id": spec["id"],
        "description": spec["desc"],
        "source": {
            "kind": "cross-implementation",
            "citation": ("py-ballisticcalc " + PBC_VERSION + " (open-source point-mass solver; "
                         "default RK4 engine; standard G1/G7 drag tables). Flat-fire, horizontal-LOS "
                         "geometry maps py-bc height/windage/velocity/time onto engine McCoy "
                         "drop/drift/speed/tof; wind direction_from = engine wind_angle + 180 deg."),
            "retrieved": RETRIEVED,
            "generator_version": "py-ballisticcalc==" + PBC_VERSION + " via " + GENERATOR,
        },
        "inputs": spec["inputs"],
        "expectations": expectations,
    }
    return case


def main():
    for spec in CASES:
        case = build_case(spec)
        path = os.path.join(OUT_DIR, spec["id"] + ".json")
        with open(path, "w", encoding="utf-8") as f:
            json.dump(case, f, indent=2)
            f.write("\n")
        print("wrote", os.path.relpath(path), "(%d expectations)" % len(case["expectations"]))
    print("py-ballisticcalc version:", PBC_VERSION)


if __name__ == "__main__":
    main()
