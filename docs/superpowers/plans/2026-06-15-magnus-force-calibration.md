# Magnus Force Calibration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the unphysical steady-coefficient Magnus model (`1.8` fudge factor, ~156″ drift at 1000 yd) with the proper yaw-of-repose-dependent McCoy Magnus force, giving a realistic sub-inch secondary effect in both solvers.

**Architecture:** The Magnus force becomes `F = ½·ρ·V²·S·C_Npα·(pd/2V)·sin(α_R)`, where `α_R` (yaw of repose) comes from the engine's existing `calculate_yaw_of_repose`, fed by the gyroscopic stability `Sg` from a new shared `miller_stability` helper. Both the cli_api solver and the `derivatives.rs` solver are updated; the `1.8` factor is removed from both.

**Tech Stack:** Rust, nalgebra (Vector3), existing `spin_drift.rs` / `derivatives.rs` physics.

---

### Task 1: Shared `miller_stability` helper

Extract the Miller gyroscopic stability calculation (currently inline in
cli_api's `apply_spin_drift`) into a reusable free function so the Magnus code
in both solvers can compute `Sg`.

**Files:**
- Modify: `src/spin_drift.rs` (add free function near the top, after imports)
- Test: `src/spin_drift.rs` (in the existing `#[cfg(test)] mod tests`)

- [ ] **Step 1: Write the failing test**

Add to the `tests` module in `src/spin_drift.rs`:

```rust
    #[test]
    fn test_miller_stability_308_168gr() {
        // .308, 168 gr, 1:12 twist, ~1.215 in length -> Sg ~ 1.1-1.3 (matches CLI "SG: 1.17")
        let sg = miller_stability(0.308, 168.0, 12.0, 1.215);
        assert!(sg > 1.0 && sg < 1.4, "expected Sg ~1.17, got {}", sg);
    }

    #[test]
    fn test_miller_stability_invalid_inputs_zero() {
        assert_eq!(miller_stability(0.0, 168.0, 12.0, 1.2), 0.0);
        assert_eq!(miller_stability(0.308, 0.0, 12.0, 1.2), 0.0);
        assert_eq!(miller_stability(0.308, 168.0, 0.0, 1.2), 0.0);
    }
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test --lib miller_stability 2>&1 | tail -20`
Expected: FAIL — `cannot find function miller_stability in this scope`.

- [ ] **Step 3: Write minimal implementation**

Add to `src/spin_drift.rs` (after the `use` lines, before `calculate_dynamic_stability`):

```rust
/// Base Miller gyroscopic stability factor (no velocity/density correction).
/// All inputs imperial: caliber/length in inches, mass in grains, twist in
/// inches-per-turn. Returns 0.0 for non-positive inputs.
///   Sg = 30 m / (t^2 d^3 l (1 + l^2)),  t,l in calibers, d in inches, m in grains
pub(crate) fn miller_stability(
    caliber_in: f64,
    weight_gr: f64,
    twist_in: f64,
    length_in: f64,
) -> f64 {
    if caliber_in <= 0.0 || weight_gr <= 0.0 || twist_in <= 0.0 || length_in <= 0.0 {
        return 0.0;
    }
    let twist_cal = twist_in / caliber_in;
    let l_cal = length_in / caliber_in;
    let denom = twist_cal * twist_cal * caliber_in.powi(3) * l_cal * (1.0 + l_cal * l_cal);
    if denom == 0.0 {
        return 0.0;
    }
    30.0 * weight_gr / denom
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test --lib miller_stability 2>&1 | tail -20`
Expected: PASS (2 tests).

- [ ] **Step 5: Commit**

```bash
git add src/spin_drift.rs
git commit -m "Add shared miller_stability helper"
```

---

### Task 2: Refactor cli_api `apply_spin_drift` to use the helper (no behavior change)

**Files:**
- Modify: `src/cli_api.rs` — `apply_spin_drift` (around lines 408–438)

- [ ] **Step 1: Replace the inline Miller math with the helper**

In `apply_spin_drift`, replace this:

```rust
        let twist_cal = twist_in / d_in;
        // Real length-to-diameter ratio when available, else 4.5 cal (typical match bullet).
        let l_cal = if self.inputs.bullet_length > 0.0 {
            self.inputs.bullet_length / self.inputs.bullet_diameter
        } else {
            4.5
        };
        // Miller stability factor (Sg).
        let sg = 30.0 * m_gr
            / (twist_cal * twist_cal * d_in.powi(3) * l_cal * (1.0 + l_cal * l_cal));
```

with:

```rust
        // Real length when available, else 4.5 cal (typical match bullet).
        let length_in = if self.inputs.bullet_length > 0.0 {
            self.inputs.bullet_length / 0.0254
        } else {
            4.5 * d_in
        };
        let sg = crate::spin_drift::miller_stability(d_in, m_gr, twist_in, length_in);
```

(`d_in`, `m_gr`, `twist_in` are already defined above this block.)

- [ ] **Step 2: Verify the CLI golden master is unchanged**

Establish a fresh baseline from the current committed branch, then confirm
this refactor changes nothing:

```bash
cargo build --release
# capture baseline BEFORE this change is built in — if not already saved:
# (the verification script from earlier sessions lives at /tmp/mccoy_golden/verify.sh)
/tmp/mccoy_golden/verify.sh 2>&1 | head -3
```
Expected: `GOLDEN: IDENTICAL ✓` (spin-drift math is numerically identical; the
old code used `length/diameter` in meters = `length_in/d_in`, the helper uses
`length_in` and `caliber_in` — same ratio, same `Sg`).

If `/tmp/mccoy_golden` no longer exists, rebuild the baseline from `git stash`:
`git stash && cargo build --release && <run the trajectory commands to a file> && git stash pop`, then diff. See the spec's Verification section for the command set (table/json/csv × imperial/metric × effects × wind × zero).

- [ ] **Step 3: Run the full test suite**

Run: `cargo test 2>&1 | grep -E "test result:|FAILED"`
Expected: all pass.

- [ ] **Step 4: Commit**

```bash
git add src/cli_api.rs
git commit -m "Refactor apply_spin_drift to use miller_stability helper"
```

---

### Task 3: Proper Magnus force in the cli_api solver

**Files:**
- Modify: `src/cli_api.rs` — Magnus block in `calculate_acceleration` (around lines 1301–1338)
- Test: `tests/magnus_calibration.rs` (new integration test)

- [ ] **Step 1: Write the failing test**

Create `tests/magnus_calibration.rs`:

```rust
// Magnus force must be a small secondary effect: lateral drift far below the
// gyroscopic spin drift, and below ~2 inches at 1000 yd.
use ballistics_engine::{
    AtmosphericConditions, BallisticInputs, TrajectorySolver, WindConditions,
};

fn lateral_drift_m(enable_magnus: bool, enable_spin: bool) -> f64 {
    let mut inputs = BallisticInputs::default();
    inputs.muzzle_velocity = 823.0; // m/s (~2700 fps)
    inputs.bullet_mass = 168.0 * 0.00006479891; // kg
    inputs.bullet_diameter = 0.308 * 0.0254; // m
    inputs.bullet_length = 1.215 * 0.0254; // m
    inputs.bc_value = 0.475;
    inputs.twist_rate = 12.0;
    inputs.is_twist_right = true;
    inputs.enable_magnus = enable_magnus;
    inputs.use_enhanced_spin_drift = enable_spin;
    // McCoy frame: Z is lateral
    let mut solver =
        TrajectorySolver::new(inputs, WindConditions::default(), AtmosphericConditions::default());
    solver.set_max_range(914.4); // 1000 yd
    let r = solver.solve().unwrap();
    r.points.last().unwrap().position.z
}

#[test]
fn magnus_is_small_secondary_effect() {
    let magnus = lateral_drift_m(true, false).abs();
    let spin = lateral_drift_m(false, true).abs();
    assert!(magnus > 0.0, "Magnus should produce some drift, got {magnus}");
    assert!(
        magnus < 0.05,
        "Magnus drift should be < ~2in (0.05m) at 1000yd, got {magnus} m"
    );
    assert!(
        magnus < spin,
        "Magnus drift ({magnus} m) should be far below spin drift ({spin} m)"
    );
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test --test magnus_calibration 2>&1 | tail -20`
Expected: FAIL — the current `1.8`-factor Magnus produces ~3–4 m at 1000 yd, so `magnus < 0.05` fails.

- [ ] **Step 3: Rewrite the Magnus block**

In `src/cli_api.rs` `calculate_acceleration`, replace the body of the
`if self.inputs.enable_magnus && ...` block (the part from `let (_, spin_rad_s)`
down to `let magnus_force = MAGNUS_CALIBRATION_FACTOR ... * c_l;`) with:

```rust
            let (_, spin_rad_s) =
                crate::spin_drift::calculate_spin_rate(velocity_magnitude, self.inputs.twist_rate);
            let temp_k = self.atmosphere.temperature + 273.15;
            let speed_of_sound = (1.4 * 287.05 * temp_k).sqrt();
            let mach = velocity_magnitude / speed_of_sound;

            // Imperial conversions for the stability / yaw-of-repose helpers.
            let d_in = self.inputs.bullet_diameter / 0.0254;
            let m_gr = self.inputs.bullet_mass / 0.00006479891;
            let l_in = if self.inputs.bullet_length > 0.0 {
                self.inputs.bullet_length / 0.0254
            } else {
                4.5 * d_in
            };
            let sg = crate::spin_drift::miller_stability(d_in, m_gr, self.inputs.twist_rate, l_in);

            // Yaw of repose (radians); zero for unstable bullets (Sg <= 1).
            let (yaw_rad, _) = crate::spin_drift::calculate_yaw_of_repose(
                sg,
                velocity_magnitude,
                spin_rad_s,
                0.0, // crosswind handled elsewhere
                0.0, // pitch rate not tracked
                air_density,
                d_in,
                l_in,
                m_gr,
                mach,
                "match",
                false,
            );

            // Proper McCoy Magnus FORCE: F = q S C_Npa (pd/2V) sin(alpha_R).
            let diameter_m = self.inputs.bullet_diameter; // already meters
            let spin_param = spin_rad_s * diameter_m / (2.0 * velocity_magnitude);
            let c_np = crate::derivatives::calculate_magnus_moment_coefficient(mach);
            let area = std::f64::consts::PI * (diameter_m / 2.0).powi(2);
            let magnus_force = 0.5
                * air_density
                * velocity_magnitude.powi(2)
                * area
                * c_np
                * spin_param
                * yaw_rad.sin();
```

Leave the direction/sign/apply code below it unchanged (the `velocity_unit`,
`dir`, `if dir_norm > 1e-12 && magnus_force.abs() > 1e-12 { ... }` block).

- [ ] **Step 4: Run the Magnus test to verify it passes**

Run: `cargo test --test magnus_calibration 2>&1 | tail -20`
Expected: PASS — Magnus drift is sub-inch and well under spin drift.

- [ ] **Step 5: Run full suite**

Run: `cargo test 2>&1 | grep -E "test result:|FAILED"`
Expected: all pass.

- [ ] **Step 6: Commit**

```bash
git add src/cli_api.rs tests/magnus_calibration.rs
git commit -m "Use proper yaw-dependent McCoy Magnus force in cli_api solver"
```

---

### Task 4: Proper Magnus force in the `derivatives.rs` solver

Apply the same model to the secondary solver's Magnus block so both agree.
This solver's `inputs` use imperial fields (`caliber_inches`, `weight_grains`)
and `bullet_diameter` is in inches there.

**Files:**
- Modify: `src/derivatives.rs` — Magnus block (around lines 284–315)

- [ ] **Step 1: Rewrite the force magnitude**

In `src/derivatives.rs`, inside the Magnus `if` block, replace:

```rust
            let c_la = calculate_magnus_moment_coefficient(mach);
```
...through...
```rust
            const MAGNUS_CALIBRATION_FACTOR: f64 = 1.8; // Calibrated to produce 4-6 inches drift at 200 yards
            let magnus_force_magnitude =
                MAGNUS_CALIBRATION_FACTOR * 0.5 * air_density * speed_air.powi(2) * area * c_l;
```

with (keeping the existing `diameter_m`, `spin_param`, `area`, `spin_rate_rad_s`
locals — recompute `c_l` removed):

```rust
            let c_np = calculate_magnus_moment_coefficient(mach);

            // Yaw of repose for the proper Magnus force. Imperial fields.
            let d_in = inputs.caliber_inches;
            let m_gr = inputs.weight_grains;
            let l_in = if inputs.bullet_length > 0.0 {
                inputs.bullet_length
            } else {
                4.5 * d_in.max(1e-9)
            };
            let sg = crate::spin_drift::miller_stability(d_in, m_gr, inputs.twist_rate, l_in);
            let (yaw_rad, _) = crate::spin_drift::calculate_yaw_of_repose(
                sg,
                speed_air,
                spin_rate_rad_s,
                0.0,
                0.0,
                air_density,
                d_in,
                l_in,
                m_gr,
                mach,
                "match",
                false,
            );

            // Proper McCoy Magnus FORCE: F = q S C_Npa (pd/2V) sin(alpha_R).
            let magnus_force_magnitude =
                0.5 * air_density * speed_air.powi(2) * area * c_np * spin_param * yaw_rad.sin();
```

Note: in this file `inputs.bullet_length` is in inches (see the existing
construction in `trajectory_integration.rs`), so no unit conversion is needed.
Add `use crate::spin_drift::{calculate_yaw_of_repose, miller_stability};` to the
imports at the top of `derivatives.rs` if not already importable via the
`crate::spin_drift::` path used inline above (inline path needs no `use`).

- [ ] **Step 2: Build**

Run: `cargo build 2>&1 | grep -E "^error" ; echo done`
Expected: `done` with no error lines. Fix any unused-variable warnings for
`c_l`/`spin_param` by removing now-dead locals if the compiler flags them.

- [ ] **Step 3: Run full suite**

Run: `cargo test 2>&1 | grep -E "test result:|FAILED"`
Expected: all pass. The `derivatives.rs` coriolis/magnus tests assert magnitude
(`> 1e-3`) or norms, which still hold for a small Magnus.

- [ ] **Step 4: Commit**

```bash
git add src/derivatives.rs
git commit -m "Use proper yaw-dependent McCoy Magnus force in derivatives solver"
```

---

### Task 5: Re-baseline Magnus golden cases, final verification

The Magnus magnitude intentionally changed, so Magnus-enabled golden outputs
differ now; all non-Magnus output must remain identical.

**Files:** none (verification + version bump decision is separate)

- [ ] **Step 1: Confirm non-Magnus output is still identical**

Run the golden comparison; expect differences ONLY in the `--enable-magnus`
cases:

```bash
cargo build --release
/tmp/mccoy_golden/verify.sh 2>&1 | sed -n '1,40p'
```
Expected: differences appear, and every differing block's command header
contains `--enable-magnus`. Spin-drift-only, Coriolis-only, wind, zero, and
plain runs must be byte-identical.

- [ ] **Step 2: Spot-check the new Magnus magnitude via CLI**

```bash
./target/release/ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 \
  --twist-rate 12 --max-range 1000 -o csv --full --enable-magnus 2>/dev/null | tail -1
```
Expected: the `x_yd` (drift) column is now sub-inch-scale (well under 1 yd),
not the previous ~4 m / multi-yard value.

- [ ] **Step 3: Clippy**

Run: `cargo clippy 2>&1 | grep -A2 "^error" | grep -oE "src/[a-z_]+\.rs" | sort | uniq -c`
Expected: only `src/ffi.rs` (the 10 pre-existing errors); no new files.

- [ ] **Step 4: Update docs**

In `README.md` / `CLI_USAGE.md`, if `--enable-magnus` is documented with a
magnitude or example, update it to reflect the realistic small effect. If no
magnitude is documented, no change needed.

- [ ] **Step 5: Final commit**

```bash
git add -A
git commit -m "Re-baseline Magnus golden cases; note realistic Magnus magnitude"
```

---

## Self-Review

- **Spec coverage:** proper McCoy force (Tasks 3,4) ✓; reuse yaw-of-repose A-i (Tasks 3,4) ✓; C_Npα from `calculate_magnus_moment_coefficient`, drop 1.8 (Tasks 3,4) ✓; shared Miller helper (Task 1) ✓; both solvers (Tasks 3,4) ✓; verification Magnus≪spin (Task 3 test) + golden re-baseline (Task 5) ✓; direction/gating unchanged ✓.
- **Placeholders:** none — all steps have concrete code/commands.
- **Type consistency:** `miller_stability(caliber_in, weight_gr, twist_in, length_in)` used identically in Tasks 1–4; `calculate_yaw_of_repose` arg order matches its signature in `spin_drift.rs`; `calculate_magnus_moment_coefficient(mach)` is `pub(crate)`, reachable via `crate::derivatives::` (cli_api) and bare (derivatives).
