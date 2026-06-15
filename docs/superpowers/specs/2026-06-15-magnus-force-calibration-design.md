# Magnus Force Calibration — Design

**Date:** 2026-06-15
**Status:** Approved (design), pending implementation plan

## Problem

The Magnus side force in the cli_api `TrajectorySolver` (and the parallel
`derivatives.rs` solver) is unphysically large — ~156″ of lateral drift at
1000 yd for a .308, roughly 10× the entire gyroscopic spin drift. The current
model applies a **steady** side-force coefficient with **no yaw-of-repose
term** and an empirical `MAGNUS_CALIBRATION_FACTOR = 1.8` tuned to "4–6″ at
200 yd" — i.e. it was hand-fit to spin-drift magnitudes, conflating the small
aerodynamic Magnus force with the dominant gyroscopic drift.

Current (cli_api `calculate_acceleration`):
```
F = MAGNUS_CALIBRATION_FACTOR · ½ρV² · S · (spin_param · C_la)   // no yaw term
```

## Research basis

- Lateral drift of a spin-stabilized bullet is dominated by **gyroscopic spin
  drift** (lift acting on the yaw of repose). The Litz formula
  `SD = 1.25·(Sg+1.2)·t^1.83` — already used by `--enable-spin-drift` — bundles
  this, including the Magnus contribution.
- The Magnus **force** effect on the trajectory is "usually insignificant"; its
  real role is the Magnus **moment** (stability). McCoy Eq. 10.83 gives the
  yaw of repose `βR = P·G/M` (~0.1°, small).
- Conclusion: a *correct* Magnus force is a small secondary effect —
  sub-inch to ~1″ at 1000 yd, far below the ~9″ spin drift.

Sources: ballistics.guide lesson 4; McCoy *Modern Exterior Ballistics*
(Eq. 10.83); Military Wiki / UAF physics on Magnus force significance.

## Design

Replace the steady-coefficient model with the proper yaw-dependent McCoy
Magnus force in cli_api `calculate_acceleration`:

```
F_Magnus = ½·ρ·V² · S · C_Npα · (p·d / V) · sin(α_R)
```
- `S = π d² / 4` (reference area, d in meters)
- `(p·d / V)` = spin parameter (p = spin rate rad/s)
- `α_R` = yaw of repose
- `C_Npα` = Magnus force coefficient

**Decisions (approved):**

1. **Yaw of repose — reuse existing approximation (A-i).** Use the engine's
   `calculate_yaw_of_repose` (spin_drift.rs). It needs the gyroscopic stability
   factor `Sg`; compute `Sg` via the same Miller formula already used in
   `apply_spin_drift` (factor that out into a shared helper to avoid
   duplication). Magnus is tiny regardless, so the extra rigor of a full
   `Ix`/`C_Mα` closed form is not worth the added infrastructure.

2. **C_Npα — McCoy-typical, validated.** Reuse `calculate_magnus_moment_coefficient`
   (already small, Mach-dependent: 0.015–0.030) as the force-coefficient proxy,
   and **drop the 1.8 fudge factor**. The new `sin(α_R)` term (~0.001–0.002)
   is what brings the magnitude down to the realistic band. Validate the
   resulting 1000-yd Magnus drift is sub-inch-to-~1″ and ≪ the Litz spin drift;
   tune `C_Npα` only if validation falls outside that band.

3. **Direction unchanged.** Keep the McCoy-frame horizontal-perpendicular
   direction (`v_unit × up`, twist-signed).

4. **Gating unchanged.** Stays behind `enable_magnus` (Phase 1 decoupling).

## Scope

- **cli_api.rs** `calculate_acceleration` Magnus block — rewrite per above;
  add a shared `miller_stability(...)` helper used by both this and
  `apply_spin_drift`.
- **derivatives.rs** Magnus block (secondary solver) — apply the same proper
  formula for consistency (it already has `calculate_yaw_of_repose` available
  and shares the convention).
- Remove `MAGNUS_CALIBRATION_FACTOR = 1.8` from both.

## Verification

- New test: Magnus-only lateral drift for .308/168gr at 1000 yd is in
  (0, ~2″] and strictly less than the `--enable-spin-drift` drift at the same
  range.
- Existing CLI/FFI golden masters: Magnus-enabled runs will **change** (this is
  an intended behavior change, not output-preserving). Re-baseline the
  Magnus-enabled golden cases after validating the new magnitude; all
  non-Magnus cases must stay identical.
- `cargo test`, `cargo clippy` clean (no new findings).

## Notes / non-goals

- Small double-count remains between `--enable-magnus` (now a small proper
  force) and `--enable-spin-drift` (Litz, which already bundles Magnus). It is
  negligible at the new magnitude; documented, not "fixed".
- Not implementing a full 6-DOF yaw-of-repose (Ix/C_Mα) model.
