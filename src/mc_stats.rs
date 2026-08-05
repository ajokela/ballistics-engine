//! Streaming moments and Bernoulli intervals for Monte Carlo hit statistics.
//!
//! [`Welford`] accumulates a running mean and variance in constant memory, one trial at a
//! time, so an adaptive Monte Carlo driver never needs to retain the full trial history just
//! to report a standard deviation. [`wilson_interval`] is the classic fixed-`n` confidence
//! interval for a Bernoulli hit/miss proportion, evaluated at one of three pinned confidence
//! levels ([`ConfidenceLevel`]). The anytime-valid beta-binomial-mixture confidence sequence
//! that lets an adaptive driver peek at partial results without inflating its error rate lives
//! in this same module too (Task 3 of Plan C, MBA-1352) -- it is added separately.
//!
//! No randomness lives here: this module consumes trial outcomes a caller already produced
//! (the existing Monte Carlo trial loop in `src/cli_api.rs`) and is pure `std` math only --
//! no `rand`, no `fs`, no `clap` -- so it compiles for `wasm32-unknown-unknown` unconditionally.

/// Online (streaming) mean and variance via Welford's algorithm.
///
/// Naive "sum of squares minus square of sum" variance loses precision catastrophically once
/// the values share a large common offset (see this module's tests for a worked example);
/// Welford's algorithm instead updates the mean and the sum of squared deviations from the
/// *running* mean one sample at a time, and never needs to retain the samples themselves --
/// exactly the shape an adaptive Monte Carlo driver needs: an unbounded trial count held in
/// `O(1)` memory.
#[derive(Debug, Clone, Default)]
pub struct Welford {
    n: u64,
    mean: f64,
    m2: f64,
}

impl Welford {
    /// A fresh accumulator with no observations.
    pub fn new() -> Self {
        Self::default()
    }

    /// Folds one more observation into the running mean and variance.
    pub fn push(&mut self, x: f64) {
        self.n += 1;
        let delta = x - self.mean;
        self.mean += delta / self.n as f64;
        self.m2 += delta * (x - self.mean);
    }

    /// Number of observations folded in so far.
    pub fn count(&self) -> u64 {
        self.n
    }

    /// Running mean. `0.0` when no observations have been pushed yet.
    pub fn mean(&self) -> f64 {
        self.mean
    }

    /// Unbiased (Bessel-corrected) sample variance, `m2 / (n - 1)`.
    ///
    /// `0.0` when fewer than two observations have been pushed -- a single point has no
    /// defined spread.
    pub fn sample_variance(&self) -> f64 {
        if self.n < 2 {
            0.0
        } else {
            self.m2 / (self.n - 1) as f64
        }
    }

    /// Population variance, `m2 / n`. `0.0` when no observations have been pushed yet.
    pub fn population_variance(&self) -> f64 {
        if self.n == 0 {
            0.0
        } else {
            self.m2 / self.n as f64
        }
    }

    /// Sample standard deviation, `sqrt(sample_variance())`.
    pub fn sample_std(&self) -> f64 {
        self.sample_variance().sqrt()
    }
}

/// A confidence level pinned to one of three fixed two-sided critical values.
///
/// A probit implementation (inverse normal CDF) is deliberately out of scope: three levels
/// cover every caller in this train, so the module states the constants directly rather than
/// deriving them numerically (Plan C spec 9.2).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConfidenceLevel {
    /// 90% two-sided confidence.
    P90,
    /// 95% two-sided confidence.
    P95,
    /// 99% two-sided confidence.
    P99,
}

impl ConfidenceLevel {
    /// The two-sided critical value `z` such that `P(-z <= Z <= z) == 1 - alpha` for a
    /// standard normal `Z`.
    pub fn z(&self) -> f64 {
        match self {
            Self::P90 => 1.644_853_626_951_472_2,
            Self::P95 => 1.959_963_984_540_054,
            Self::P99 => 2.575_829_303_548_900_4,
        }
    }

    /// The significance level `alpha = 1 - confidence`.
    pub fn alpha(&self) -> f64 {
        match self {
            Self::P90 => 0.10,
            Self::P95 => 0.05,
            Self::P99 => 0.01,
        }
    }

    /// The confidence level as a whole-number percentage (`90`, `95`, or `99`).
    pub fn as_percent(&self) -> u32 {
        match self {
            Self::P90 => 90,
            Self::P95 => 95,
            Self::P99 => 99,
        }
    }
}

/// Wilson score interval for a Bernoulli proportion at fixed `n`.
///
/// `n == 0` returns `(0.0, 1.0)` -- total ignorance, not an error: with zero trials nothing is
/// known about the proportion beyond it lying in `[0, 1]`, so the widest possible interval is
/// the honest answer rather than a divide-by-zero or a panic.
///
/// `successes` above `trials` is not a domain error either: it saturates at `trials` (`p =
/// 1.0`), so the call reports exactly the interval a fully-successful run of that size would.
/// This is deliberate, not merely tolerated: for `p` only just above `1.0` the radicand below
/// does not reliably go negative (`z^2/(4n^2)` can outweigh a small negative `p(1-p)/n`), so an
/// unclamped implementation would silently return a plausible-looking wrong interval rather
/// than a loud `NaN` -- saturating `p` at the input removes that failure mode entirely instead
/// of documenting around it.
///
/// Unlike the naive Wald interval (`p +- z * sqrt(p(1-p)/n)`), the Wilson interval stays
/// inside `[0, 1]` and has correct coverage even at small `n` or `p` near `0` or `1` --
/// exactly the regime a hit-probability Monte Carlo run starts in before enough trials have
/// accumulated. Computed directly in "center +- spread" form (Wilson 1927; worked example in
/// Newcombe 1998, cross-checked by this module's tests):
///
/// ```text
/// center = (p + z^2 / (2n)) / (1 + z^2 / n)
/// spread = (z / (1 + z^2 / n)) * sqrt(p(1-p)/n + z^2/(4n^2))
/// ```
///
/// where `p = min(successes / trials, 1.0)` and `z` is [`ConfidenceLevel::z`]. The result is
/// additionally clamped to `[0, 1]`: floating-point rounding in `center +- spread` can overshoot
/// by up to a few ULP right at the `p == 0` / `p == 1` edges, and a probability bound outside
/// `[0, 1]` is a worse answer than one that is merely maximally uninformative.
pub fn wilson_interval(successes: u64, trials: u64, level: ConfidenceLevel) -> (f64, f64) {
    if trials == 0 {
        return (0.0, 1.0);
    }
    let n = trials as f64;
    let p = (successes as f64 / n).min(1.0);
    let z = level.z();
    let z2 = z * z;
    let denom = 1.0 + z2 / n;
    let center = (p + z2 / (2.0 * n)) / denom;
    let spread = (z / denom) * (p * (1.0 - p) / n + z2 / (4.0 * n * n)).sqrt();
    ((center - spread).max(0.0), (center + spread).min(1.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn welford_matches_two_pass_moments_on_a_hostile_fixture() {
        // Catastrophic-cancellation fixture: large offset, small spread — the case
        // naive sum-of-squares gets wrong and Welford exists to protect.
        let xs: Vec<f64> = (0..1000).map(|i| 1.0e9 + (i % 7) as f64 * 0.25).collect();
        let mut w = Welford::new();
        for &x in &xs { w.push(x); }
        let n = xs.len() as f64;
        let mean = xs.iter().sum::<f64>() / n;
        let m2 = xs.iter().map(|x| (x - mean).powi(2)).sum::<f64>();
        assert_eq!(w.count(), 1000);
        assert!((w.mean() - mean).abs() < 1e-6);
        assert!((w.population_variance() - m2 / n).abs() < 1e-6);
        assert!((w.sample_variance() - m2 / (n - 1.0)).abs() < 1e-6);
    }

    #[test]
    fn welford_empty_and_single_are_defined() {
        let mut w = Welford::new();
        assert_eq!(w.count(), 0);
        assert_eq!(w.mean(), 0.0);
        assert_eq!(w.sample_variance(), 0.0);
        assert_eq!(w.population_variance(), 0.0);
        w.push(42.0);
        assert_eq!(w.mean(), 42.0);
        assert_eq!(w.sample_variance(), 0.0); // n < 2
        assert_eq!(w.population_variance(), 0.0);
    }

    #[test]
    fn confidence_level_accessors_are_pinned_exactly() {
        // z(): bit-for-bit against the three literals pinned in the plan's Global Constraints
        // (docs/superpowers/plans/2026-08-04-decision-support-plan-c.md). The 95% value is the
        // same number as `truing_uncertainty::NORMAL_95_TWO_SIDED_Z`
        // (src/truing_uncertainty.rs:36), but that `const` is module-private (not `pub` or
        // `pub(crate)`), so it is not reachable from here; the digits are restated rather than
        // imported. Without this test, z() was previously anchored only indirectly (P95 via
        // the Newcombe reference, P99 only by an inequality, P90 by nothing at all) -- review I3.
        assert_eq!(ConfidenceLevel::P90.z().to_bits(), 1.644_853_626_951_472_2_f64.to_bits());
        assert_eq!(ConfidenceLevel::P95.z().to_bits(), 1.959_963_984_540_054_f64.to_bits());
        assert_eq!(ConfidenceLevel::P99.z().to_bits(), 2.575_829_303_548_900_4_f64.to_bits());

        assert_eq!(ConfidenceLevel::P90.alpha(), 0.10);
        assert_eq!(ConfidenceLevel::P95.alpha(), 0.05);
        assert_eq!(ConfidenceLevel::P99.alpha(), 0.01);

        assert_eq!(ConfidenceLevel::P90.as_percent(), 90);
        assert_eq!(ConfidenceLevel::P95.as_percent(), 95);
        assert_eq!(ConfidenceLevel::P99.as_percent(), 99);
    }

    #[test]
    fn wilson_matches_the_newcombe_canonical_example() {
        // Newcombe (1998), worked example: 81/263 at 95% → (0.2553, 0.3662) to 4 dp.
        let (lo, hi) = wilson_interval(81, 263, ConfidenceLevel::P95);
        assert!((lo - 0.2553).abs() < 5e-4, "lo = {lo}");
        assert!((hi - 0.3662).abs() < 5e-4, "hi = {hi}");
    }

    #[test]
    fn wilson_two_algebraic_forms_agree() {
        // center ± spread form (production) vs the quadratic-root form (independent oracle,
        // structurally different): solving (p_hat - p)^2 = z^2 * p(1-p)/n directly for p gives
        // the quadratic p^2*(n + z^2) - p*(2*n*p_hat + z^2) + n*p_hat^2 == 0, whose roots are
        //   p = [ (2*n*p_hat + z^2) +- z * sqrt(z^2 + 4*n*p_hat*(1 - p_hat)) ] / (2*(n + z^2))
        // A mis-transcription of either form fails this cross-check (review I2: the previous
        // version of this test restated the production expression term-for-term and could not).
        for &(k, n) in &[(0_u64, 10_u64), (10, 10), (1, 10), (8, 10), (500, 1000), (3, 7)] {
            for level in [ConfidenceLevel::P90, ConfidenceLevel::P95, ConfidenceLevel::P99] {
                let (lo, hi) = wilson_interval(k, n, level);
                let z = level.z();
                let (kf, nf) = (k as f64, n as f64);
                let p_hat = kf / nf;
                let z2 = z * z;
                let a = nf + z2;
                let b = 2.0 * nf * p_hat + z2;
                let root_term = z * (z2 + 4.0 * nf * p_hat * (1.0 - p_hat)).sqrt();
                let lo_root = (b - root_term) / (2.0 * a);
                let hi_root = (b + root_term) / (2.0 * a);
                assert!((lo - lo_root).abs() < 1e-12, "lo mismatch at k={k} n={n}: {lo} vs {lo_root}");
                assert!((hi - hi_root).abs() < 1e-12, "hi mismatch at k={k} n={n}: {hi} vs {hi_root}");
            }
        }
    }

    #[test]
    // The brief's bound checks are intentionally spelled out as explicit `>=`/`<` comparisons
    // (matching the two-sided prose "lower bound must be 0" / "upper bound must be 1") rather
    // than `Range::contains`, which clippy's default lint set does not recognize as equivalent
    // here.
    #[allow(clippy::manual_range_contains)]
    fn wilson_edges_and_ordering_properties() {
        assert_eq!(wilson_interval(0, 0, ConfidenceLevel::P95), (0.0, 1.0));
        // k=0 / k=n bounds must land in [0, 1] EXACTLY (unclamped center +- spread can overshoot
        // by a few ULP at these edges; `wilson_interval` clamps for it) and stay tight against
        // 0 / 1 (exact in real arithmetic for every n, at every level). Swept rather than probed
        // at a single n: n = 20 at P95 happens to round to exactly 0.0 / 1.0 even unclamped,
        // which gave false assurance that the property held in general (review I1).
        for n in 1_u64..=200 {
            for level in [ConfidenceLevel::P90, ConfidenceLevel::P95, ConfidenceLevel::P99] {
                let (lo0, _) = wilson_interval(0, n, level);
                assert!(lo0 >= 0.0 && lo0 < 1e-9, "k=0 n={n} {level:?}: lo = {lo0:e}, want [0, 1e-9)");
                let (_, hi_n) = wilson_interval(n, n, level);
                assert!(hi_n <= 1.0 && hi_n > 1.0 - 1e-9, "k=n n={n} {level:?}: hi = {hi_n:e}, want (1-1e-9, 1]");
            }
        }
        // Wider at higher confidence, narrower at larger n.
        let w = |k, n, l| { let (a, b) = wilson_interval(k, n, l); b - a };
        assert!(w(40, 100, ConfidenceLevel::P99) > w(40, 100, ConfidenceLevel::P95));
        assert!(w(40, 100, ConfidenceLevel::P95) > w(400, 1000, ConfidenceLevel::P95));
        // Interval always contains the point estimate.
        let (lo, hi) = wilson_interval(3, 7, ConfidenceLevel::P90);
        let p = 3.0 / 7.0;
        assert!(lo < p && p < hi);
    }

    #[test]
    fn wilson_interval_saturates_when_successes_exceeds_trials() {
        // successes > trials is not a panic and not a silent NaN: it saturates at trials
        // (p = 1.0), reporting exactly the k=n interval, for every level (review I4).
        for level in [ConfidenceLevel::P90, ConfidenceLevel::P95, ConfidenceLevel::P99] {
            assert_eq!(wilson_interval(37, 20, level), wilson_interval(20, 20, level));
            assert_eq!(wilson_interval(u64::MAX, 20, level), wilson_interval(20, 20, level));
        }
    }
}
