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
/// where `p = successes / trials` and `z` is [`ConfidenceLevel::z`].
pub fn wilson_interval(successes: u64, trials: u64, level: ConfidenceLevel) -> (f64, f64) {
    if trials == 0 {
        return (0.0, 1.0);
    }
    let n = trials as f64;
    let p = successes as f64 / n;
    let z = level.z();
    let z2 = z * z;
    let denom = 1.0 + z2 / n;
    let center = (p + z2 / (2.0 * n)) / denom;
    let spread = (z / denom) * (p * (1.0 - p) / n + z2 / (4.0 * n * n)).sqrt();
    (center - spread, center + spread)
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
    fn wilson_matches_the_newcombe_canonical_example() {
        // Newcombe (1998), worked example: 81/263 at 95% → (0.2553, 0.3662) to 4 dp.
        let (lo, hi) = wilson_interval(81, 263, ConfidenceLevel::P95);
        assert!((lo - 0.2553).abs() < 5e-4, "lo = {lo}");
        assert!((hi - 0.3662).abs() < 5e-4, "hi = {hi}");
    }

    #[test]
    fn wilson_two_algebraic_forms_agree() {
        // center ± spread form vs the quadratic-root form — an implementation
        // transposition in either fails this cross-check.
        for &(k, n) in &[(0_u64, 10_u64), (10, 10), (1, 10), (8, 10), (500, 1000), (3, 7)] {
            for level in [ConfidenceLevel::P90, ConfidenceLevel::P95, ConfidenceLevel::P99] {
                let (lo, hi) = wilson_interval(k, n, level);
                let z = level.z();
                let (kf, nf) = (k as f64, n as f64);
                let p = kf / nf;
                let z2 = z * z;
                let denom = 1.0 + z2 / nf;
                let center = (p + z2 / (2.0 * nf)) / denom;
                let spread = (z / denom) * ((p * (1.0 - p) / nf) + z2 / (4.0 * nf * nf)).sqrt();
                assert!((lo - (center - spread)).abs() < 1e-12);
                assert!((hi - (center + spread)).abs() < 1e-12);
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
        let (lo0, _) = wilson_interval(0, 20, ConfidenceLevel::P95);
        let (_, hi_n) = wilson_interval(20, 20, ConfidenceLevel::P95);
        assert!(lo0 >= 0.0 && lo0 < 1e-9, "k=0 lower bound must be 0, got {lo0}");
        assert!(hi_n <= 1.0 && hi_n > 1.0 - 1e-9, "k=n upper bound must be 1");
        // Wider at higher confidence, narrower at larger n.
        let w = |k, n, l| { let (a, b) = wilson_interval(k, n, l); b - a };
        assert!(w(40, 100, ConfidenceLevel::P99) > w(40, 100, ConfidenceLevel::P95));
        assert!(w(40, 100, ConfidenceLevel::P95) > w(400, 1000, ConfidenceLevel::P95));
        // Interval always contains the point estimate.
        let (lo, hi) = wilson_interval(3, 7, ConfidenceLevel::P90);
        let p = 3.0 / 7.0;
        assert!(lo < p && p < hi);
    }
}
