//! Streaming moments and Bernoulli intervals for Monte Carlo hit statistics.
//!
//! [`Welford`] accumulates a running mean and variance in constant memory, one trial at a
//! time, so an adaptive Monte Carlo driver never needs to retain the full trial history just
//! to report a standard deviation. [`wilson_interval`] is the classic fixed-`n` confidence
//! interval for a Bernoulli hit/miss proportion, evaluated at one of three pinned confidence
//! levels ([`ConfidenceLevel`]). [`BernoulliConfidenceSequence`] is the anytime-valid
//! beta-binomial-mixture counterpart to [`wilson_interval`]: it lets an adaptive driver peek at
//! partial results after every trial, and stop on a data-dependent rule, without inflating its
//! error rate (Task 3 of Plan C, MBA-1352).
//!
//! No randomness lives in the production path: this module consumes trial outcomes a caller
//! already produced (the existing Monte Carlo trial loop in `src/cli_api.rs`) and is pure `std`
//! math only -- no `rand`, no `fs`, no `clap` -- so it compiles for `wasm32-unknown-unknown`
//! unconditionally. (One `#[cfg(test)]` test, the empirical coverage check, does draw from a
//! seeded `rand` generator; nothing outside `#[cfg(test)]` does.)

use crate::special::ln_beta;

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

/// Bracket floor for the endpoint search: roots are hunted on `(EPS, p_hat]` and
/// `[p_hat, 1 - EPS)` rather than on the open unit interval, because `ln M_n` diverges at both
/// ends and a bracket endpoint has to be a finite number. `1e-15` is ~9 ULP away from `0.0` on
/// the low side and ~9 ULP from `1.0` on the high side, so the excluded slivers are narrower
/// than any proportion this crate could report meaningfully; a root inside a sliver is reported
/// as the saturated `0.0` / `1.0`, which widens the interval and therefore cannot break
/// coverage.
const BRACKET_EPS: f64 = 1e-15;

/// Bisection steps taken per endpoint by [`BernoulliConfidenceSequence::bounds`]. Fixed, not
/// tolerance-driven -- see that method's docs for why.
const BISECTION_ITERS: u32 = 200;

/// `ln M_n(p)`, the log mixture martingale, given the precomputed prior term
/// `ln_b = ln B(S+1, n-S+1)` (constant across an endpoint search, so it is hoisted out).
///
/// The two guards apply the `0 * ln(0) == 0` convention, which makes this a total function at
/// `p == 0` and `p == 1` (where `S == 0` / `S == n` respectively kill the divergent term) rather
/// than a `0.0 * -INFINITY == NaN`. The bisection itself never evaluates outside
/// `[BRACKET_EPS, 1 - BRACKET_EPS]`, so the guards are a correctness margin, not a hot path.
///
/// `(1.0 - p).ln()` is used rather than `(-p).ln_1p()`: for `p >= 0.5` the subtraction is exact
/// (Sterbenz), and for `p < 0.5` the operand is at least `0.5`, so its half-ULP rounding is a
/// relative `1e-16` that contributes an absolute `~1e-16` to a logarithm -- far below the
/// `ln_beta` error budget quantified in [`BernoulliConfidenceSequence::bounds`].
fn ln_mixture(ln_b: f64, successes: f64, failures: f64, p: f64) -> f64 {
    let mut ln_m = ln_b;
    if successes > 0.0 {
        ln_m -= successes * p.ln();
    }
    if failures > 0.0 {
        ln_m -= failures * (1.0 - p).ln();
    }
    ln_m
}

/// Bisects for the endpoint where `ln M_n(p) == threshold`, on a bracket whose `outside` end is
/// known to exceed the threshold and whose `inside` end is known not to.
///
/// Direction-agnostic on purpose: the lower endpoint passes `outside = BRACKET_EPS,
/// inside = p_hat` and the upper passes `outside = 1 - BRACKET_EPS, inside = p_hat`, so both
/// endpoints run the identical loop and cannot drift apart under maintenance. Midpoints stay
/// strictly inside the initial bracket, so `p` never reaches `0.0` or `1.0`.
fn bisect_endpoint(
    ln_b: f64,
    successes: f64,
    failures: f64,
    threshold: f64,
    mut outside: f64,
    mut inside: f64,
) -> f64 {
    for _ in 0..BISECTION_ITERS {
        let mid = 0.5 * (outside + inside);
        if ln_mixture(ln_b, successes, failures, mid) > threshold {
            outside = mid;
        } else {
            inside = mid;
        }
    }
    0.5 * (outside + inside)
}

/// Anytime-valid confidence sequence for a Bernoulli proportion: Robbins' beta-binomial mixture
/// with a uniform `Beta(1, 1)` prior.
///
/// # Why not just re-run [`wilson_interval`]
///
/// A fixed-`n` interval is only valid if `n` was fixed *before* the data. An adaptive driver
/// that watches the interval after every trial and stops when it looks tight enough is doing
/// optional stopping, and a fixed-`n` interval checked that way has no coverage guarantee at
/// all: the error rate grows with the number of peeks. A confidence sequence is instead a whole
/// family of intervals, one per `n`, that are *simultaneously* valid:
///
/// ```text
/// P( there exists n >= 1 with p_true outside CS_n )  <=  alpha
/// ```
///
/// so the caller may look as often as it likes, stop on any rule it likes -- including a
/// data-dependent one -- and the interval it stops on still covers at the advertised rate. The
/// price is width: at every `n` this interval is strictly wider than the Wilson interval on the
/// same counts (pinned by this module's tests). Paying it is the point.
///
/// # The mixture martingale
///
/// Robbins (1970), "Statistical methods related to the law of the iterated logarithm"; the
/// modern treatment is Howard, Ramdas, McAuliffe & Sekhon (2021), "Time-uniform, nonasymptotic,
/// nonparametric confidence sequences". Against a candidate value `p`, with `S` successes in
/// `n` trials, mix the likelihood ratio over a uniform `Beta(1, 1)` prior on the alternative:
///
/// ```text
/// M_n(p) = integral_0^1 q^S (1-q)^(n-S) dq / [ p^S (1-p)^(n-S) ]
///        = B(S+1, n-S+1) / [ p^S (1-p)^(n-S) ]
/// ```
///
/// The prior's own normaliser `1 / B(1, 1)` is `1` (`ln B(1,1) == 0`), so it drops out. In logs,
/// via [`crate::special::ln_beta`]:
///
/// ```text
/// ln M_n(p) = ln B(S+1, n-S+1) - S*ln(p) - (n-S)*ln(1-p)
/// ```
///
/// When `p` is the true success probability, `M_n` is a nonnegative martingale with `M_0 = 1`,
/// so Ville's inequality bounds it uniformly over all time: `P(exists n : M_n >= 1/alpha) <=
/// alpha`. Inverting that test -- keeping every `p` the data has not yet ruled out -- gives the
/// confidence set
///
/// ```text
/// CS_n = { p in (0, 1) : ln M_n(p) <= ln(1/alpha) }
/// ```
///
/// whose miss probability, over the whole infinite sequence at once, is at most `alpha`.
///
/// # Asymptotics
///
/// Substituting Stirling into `ln B(S+1, n-S+1) = -ln[(n+1) * C(n, S)]` gives
/// `ln M_n(p) ~ n*KL(p_hat || p) - ln(n+1) + 0.5*ln(2*pi*n*p_hat*(1-p_hat))`, so the half-width
/// shrinks like `sqrt((ln(1/alpha) + 0.5*ln n) * p_hat*(1-p_hat) * 2 / n)` -- the familiar
/// `sqrt(log(n)/n)` rate of a confidence sequence, a `sqrt(log n)` factor wider than the
/// `sqrt(1/n)` of a fixed-`n` interval. That extra factor *is* the optional-stopping licence.
#[derive(Debug, Clone)]
pub struct BernoulliConfidenceSequence {
    successes: u64,
    trials: u64,
    alpha: f64,
}

impl BernoulliConfidenceSequence {
    /// A fresh sequence with no observations, at the given confidence level.
    pub fn new(level: ConfidenceLevel) -> Self {
        Self {
            successes: 0,
            trials: 0,
            alpha: level.alpha(),
        }
    }

    /// Folds in one more trial.
    pub fn update(&mut self, hit: bool) {
        self.trials = self.trials.saturating_add(1);
        if hit {
            self.successes = self.successes.saturating_add(1);
        }
    }

    /// Folds in `trials` more trials of which `hits` succeeded.
    ///
    /// `hits` above `trials` saturates at `trials` rather than panicking or asserting, matching
    /// [`wilson_interval`]'s posture for the same malformed input so the two functions cannot
    /// disagree about what a bad count means. That choice is load-bearing here rather than
    /// merely tidy: `S > n` would make the second argument of `ln B(S+1, n-S+1)` non-positive,
    /// [`crate::special::ln_gamma`] returns `NaN` out of domain, and `NaN > threshold` is
    /// `false` -- so `bisect_endpoint` would take its `else` arm 200 times running and return
    /// a perfectly plausible-looking number pinned to the bracket end instead of signalling
    /// anything. Saturating at the input removes that silent-wrong-answer path entirely. A
    /// `debug_assert!` was considered and rejected for the reason recorded on `wilson_interval`:
    /// it would panic under plain `cargo test`, leaving the documented behaviour unpinnable by
    /// an ordinary test.
    ///
    /// Both counters saturate at `u64::MAX` too, which preserves the `successes <= trials`
    /// invariant (`min(S+h, MAX) <= min(n+t, MAX)` whenever `S <= n` and `h <= t`).
    pub fn update_batch(&mut self, hits: u64, trials: u64) {
        self.trials = self.trials.saturating_add(trials);
        self.successes = self.successes.saturating_add(hits.min(trials));
    }

    /// Successes observed so far.
    pub fn successes(&self) -> u64 {
        self.successes
    }

    /// Trials observed so far.
    pub fn trials(&self) -> u64 {
        self.trials
    }

    /// The current interval: every `p` not yet ruled out at level `alpha`.
    ///
    /// `(0.0, 1.0)` when no trials have been folded in -- total ignorance, the same answer
    /// [`wilson_interval`] gives at `n == 0`, and not an error.
    ///
    /// # How the endpoints are found
    ///
    /// `d^2/dp^2 ln M_n(p) = S/p^2 + (n-S)/(1-p)^2 > 0`, so `ln M_n` is strictly convex on
    /// `(0, 1)` with its unique minimum at the MLE `p_hat = S/n`. The confidence set is
    /// therefore a genuine interval, bounded by the two solutions of
    /// `ln M_n(p) = ln(1/alpha)`, one on each side of `p_hat`.
    ///
    /// `p_hat` is always strictly inside the set, so both brackets are always valid: the
    /// minimum value is `ln M_n(p_hat) = -ln[(n+1) * C(n,S) * p_hat^S * (1-p_hat)^(n-S)]`, and
    /// the method-of-types bound (Cover & Thomas, *Elements of Information Theory*, Thm 11.1.4:
    /// a type's probability under its own maximum-likelihood distribution is at least
    /// `1/(n+1)^(|alphabet|-1)`, here `1/(n+1)`) makes that bracketed product at least `1`,
    /// hence `ln M_n(p_hat) <= 0 < ln(1/alpha)` for every `n >= 1` and every level.
    ///
    /// Each endpoint is then bisected with exactly `BISECTION_ITERS` (200) steps -- a fixed count,
    /// not a tolerance test. Every step halves the bracket, so after 200 the bracket is
    /// `2^-200 ~ 6e-61` of a starting width below `1.0`: the loop has been sitting on two
    /// adjacent doubles since roughly step 55, and the remaining steps are no-ops that keep
    /// `bounds()` a pure function of `(S, n, alpha)` with no data-dependent trip count. Since
    /// the search resolves `p` to the last bit, accuracy is set by `ln M_n` itself, i.e. by
    /// [`crate::special::ln_beta`]: its measured `~7e-16` relative error at `n = 1e6` is an
    /// absolute `~5e-10` on an `ln B` of magnitude `~7e5`, and dividing by the slope of
    /// `ln M_n` at the endpoint (`~9e3` there) puts the endpoint error near `1e-13` -- ten
    /// orders of magnitude below the `~2e-3` half-width at that `n`.
    ///
    /// # Saturation
    ///
    /// Each side saturates for either of two independent reasons.
    ///
    /// **Rule 1, the `p_hat` rule.** The lower root can never sit above `p_hat`, so
    /// `p_hat <= BRACKET_EPS` (which includes `S == 0`, where `p_hat` is exactly `0`) means the
    /// root is below anything the bracket can resolve: the bound is reported as exactly `0.0`.
    /// The upper side mirrors it at `p_hat >= 1 - BRACKET_EPS` (including `S == n`, reported as
    /// exactly `1.0`). This rule is what handles every small-`n` case, on both sides.
    ///
    /// **Rule 2, the endpoint-value guard.** Independently, `ln M_n` evaluated *at* the bracket
    /// end can already sit at or below the threshold, meaning no root exists inside the bracket
    /// at all; the bound then saturates rather than bisecting for something that is not there.
    /// Despite the "tiny `n`" intuition this guard is **not** load-bearing at small `n` -- at
    /// `S = 1, n = 1` and 95% confidence, `ln M_n(eps) = 33.85` against a threshold of `2.996`,
    /// 30 nats away from firing, and the low bound there is a genuine root at `0.025`. Every
    /// small-`n` case that does satisfy this guard is already caught by rule 1. It first becomes
    /// non-redundant in the opposite regime -- large `n` with an extreme `p_hat`, where the true
    /// root falls below `BRACKET_EPS` while `p_hat` itself does not. Measured first firings:
    /// `S = 1` at `n ~ 1e7` (giving `(0.0, 2.22e-6)`), and `S = 2` at `n ~ 1e10`. Do not read it
    /// as dead code because a small-`n` probe never reaches it.
    ///
    /// Both saturations widen the interval, so neither can cost coverage.
    pub fn bounds(&self) -> (f64, f64) {
        if self.trials == 0 {
            return (0.0, 1.0);
        }
        let n = self.trials as f64;
        // The invariant is maintained by both update paths; re-imposed here so a future
        // constructor cannot leak an S > n state into ln_gamma's NaN domain.
        let successes = self.successes.min(self.trials) as f64;
        let failures = n - successes;
        let p_hat = successes / n;
        let threshold = -self.alpha.ln(); // ln(1/alpha)
        let ln_b = ln_beta(successes + 1.0, failures + 1.0);

        let lower = if p_hat <= BRACKET_EPS
            || ln_mixture(ln_b, successes, failures, BRACKET_EPS) <= threshold
        {
            0.0
        } else {
            bisect_endpoint(ln_b, successes, failures, threshold, BRACKET_EPS, p_hat)
        };

        let upper = if p_hat >= 1.0 - BRACKET_EPS
            || ln_mixture(ln_b, successes, failures, 1.0 - BRACKET_EPS) <= threshold
        {
            1.0
        } else {
            bisect_endpoint(
                ln_b,
                successes,
                failures,
                threshold,
                1.0 - BRACKET_EPS,
                p_hat,
            )
        };

        (lower, upper)
    }

    /// Half the width of the current interval; `0.5` before any trials.
    pub fn half_width(&self) -> f64 {
        let (lo, hi) = self.bounds();
        (hi - lo) / 2.0
    }
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

    #[test]
    fn cs_starts_ignorant_and_shrinks_monotonically_in_n() {
        let mut cs = BernoulliConfidenceSequence::new(ConfidenceLevel::P95);
        assert_eq!(cs.bounds(), (0.0, 1.0));
        let mut prev = 1.0_f64;
        // Alternate hits/misses so p̂ stays near 0.5 while n grows.
        for i in 0..2000 {
            cs.update(i % 2 == 0);
            if i % 100 == 99 {
                let hw = cs.half_width();
                assert!(hw <= prev + 1e-12, "half-width grew at n={}: {hw} > {prev}", i + 1);
                prev = hw;
            }
        }
        let (lo, hi) = cs.bounds();
        assert!(lo > 0.0 && hi < 1.0 && lo < 0.5 && 0.5 < hi);
    }

    #[test]
    fn cs_edges_saturate_exactly() {
        let mut cs = BernoulliConfidenceSequence::new(ConfidenceLevel::P95);
        cs.update_batch(0, 50);
        let (lo, hi) = cs.bounds();
        assert_eq!(lo, 0.0);
        assert!(hi < 0.25 && hi > 0.0);
        let mut cs2 = BernoulliConfidenceSequence::new(ConfidenceLevel::P95);
        cs2.update_batch(50, 50);
        let (lo2, hi2) = cs2.bounds();
        assert_eq!(hi2, 1.0);
        assert!(lo2 > 0.75 && lo2 < 1.0);
    }

    #[test]
    fn cs_is_wider_than_wilson_at_the_same_n() {
        // The honesty property that motivates anytime-validity: the price of
        // optional stopping is a wider interval than the fixed-n Wilson at every n.
        for &(k, n) in &[(10_u64, 40_u64), (81, 263), (500, 1000)] {
            let mut cs = BernoulliConfidenceSequence::new(ConfidenceLevel::P95);
            cs.update_batch(k, n);
            let (clo, chi) = cs.bounds();
            let (wlo, whi) = wilson_interval(k, n, ConfidenceLevel::P95);
            assert!(chi - clo > whi - wlo, "CS not wider than Wilson at k={k}, n={n}");
        }
    }

    #[test]
    fn cs_bounds_are_deterministic_and_batch_order_invariant() {
        let mut a = BernoulliConfidenceSequence::new(ConfidenceLevel::P99);
        let mut b = BernoulliConfidenceSequence::new(ConfidenceLevel::P99);
        a.update_batch(30, 100);
        for i in 0..100 { b.update(i < 30); }
        assert_eq!(a.bounds(), b.bounds()); // state is (S, n) only — exact equality
    }

    /// The structural preconditions `bounds()` relies on, swept rather than argued: the
    /// endpoints bracket the MLE (so the two bisections really did straddle the root, and
    /// `ln M_n(p_hat) <= 0 < ln(1/alpha)` really does hold at every `(S, n)` — the
    /// method-of-types claim in the `bounds()` docs), the interval stays inside `[0, 1]`,
    /// and lowering alpha only ever widens it. A sign error or a swapped bracket in the
    /// endpoint search fails this even where the brief's fixtures happen to look sane.
    #[test]
    fn cs_brackets_the_mle_and_nests_by_confidence_level() {
        for n in 1_u64..=60 {
            for s in 0..=n {
                let p_hat = s as f64 / n as f64;
                let mut widths = Vec::new();
                for level in [ConfidenceLevel::P90, ConfidenceLevel::P95, ConfidenceLevel::P99] {
                    let mut cs = BernoulliConfidenceSequence::new(level);
                    cs.update_batch(s, n);
                    let (lo, hi) = cs.bounds();
                    assert!(lo.is_finite() && hi.is_finite(), "non-finite bound at s={s} n={n}");
                    assert!((0.0..=1.0).contains(&lo), "lo={lo} out of [0,1] at s={s} n={n}");
                    assert!((0.0..=1.0).contains(&hi), "hi={hi} out of [0,1] at s={s} n={n}");
                    assert!(lo <= p_hat, "lo={lo} above p_hat={p_hat} at s={s} n={n} {level:?}");
                    assert!(hi >= p_hat, "hi={hi} below p_hat={p_hat} at s={s} n={n} {level:?}");
                    assert!((hi - lo - 2.0 * cs.half_width()).abs() < 1e-15);
                    widths.push(hi - lo);
                }
                assert!(widths[0] <= widths[1], "P90 wider than P95 at s={s} n={n}");
                assert!(widths[1] <= widths[2], "P95 wider than P99 at s={s} n={n}");
            }
        }
        // Malformed counts saturate instead of reaching ln_gamma's NaN domain, matching
        // `wilson_interval`'s posture (see `update_batch`'s docs).
        let mut bad = BernoulliConfidenceSequence::new(ConfidenceLevel::P95);
        bad.update_batch(u64::MAX, 20);
        assert_eq!(bad.successes(), 20);
        assert_eq!(bad.trials(), 20);
        let (lo, hi) = bad.bounds();
        assert!(lo > 0.0 && !lo.is_nan(), "lo={lo}");
        assert_eq!(hi, 1.0);
    }

    /// Exact binomial CDF `P(X <= k)` for `X ~ Bin(n, p)`, used to set the coverage floor
    /// below. Forward pmf recursion — `pmf(0) = (1-p)^n`, `pmf(i) = pmf(i-1) * ((n-i+1)/i) *
    /// (p/(1-p))` — so no factorial ever has to be represented, and no special function from
    /// this crate is involved (the floor stays an independent oracle).
    fn binomial_cdf(k: usize, n: usize, p: f64) -> f64 {
        let ratio = p / (1.0 - p);
        let mut term = (1.0 - p).powi(n as i32);
        let mut total = term;
        for i in 1..=k.min(n) {
            term *= ratio * (n - i + 1) as f64 / i as f64;
            total += term;
        }
        total.min(1.0)
    }

    /// Empirical coverage under OPTIONAL STOPPING — the spec's acceptance test,
    /// following the exact-binomial-tail methodology of tests/truing_uncertainty.rs
    /// (fixed committed seed; a floor that nominal coverage essentially never
    /// violates but material undercoverage reliably does; an anti-width guard so
    /// trivially-wide intervals cannot pass).
    #[test]
    fn cs_coverage_survives_optional_stopping() {
        // `RngExt`, not the brief's `Rng`: this crate pins rand 0.10, which moved the
        // `random::<T>()` extension method off `Rng` onto `RngExt`.
        use rand::{rngs::StdRng, RngExt, SeedableRng};
        const TRIALS: usize = 200;
        const P_TRUE: f64 = 0.30;
        let mut rng = StdRng::seed_from_u64(0x1352_C0FF_EE00_0001);
        let mut covered = 0_usize;
        let mut final_half_widths = 0.0_f64;
        for _ in 0..TRIALS {
            let mut cs = BernoulliConfidenceSequence::new(ConfidenceLevel::P95);
            // Optional stopping: stop the moment the half-width crosses 0.08,
            // or at 5000 draws. This is exactly the usage pattern that breaks
            // a repeated-Wilson check.
            for _ in 0..5000 {
                cs.update(rng.random::<f64>() < P_TRUE);
                if cs.trials() >= 50 && cs.half_width() <= 0.08 { break; }
            }
            let (lo, hi) = cs.bounds();
            if lo <= P_TRUE && P_TRUE <= hi { covered += 1; }
            final_half_widths += cs.half_width();
        }
        let mean_hw = final_half_widths / TRIALS as f64;
        eprintln!("MBA-1352 CS coverage: {covered}/{TRIALS}, mean final half-width {mean_hw:.5}");

        // Floor by exact binomial tail at the *advertised* level, p = 0.95 — the worst case a
        // correct anytime-valid construction has to survive (its true coverage is >= 0.95, and
        // measurably higher than that here, so this is conservative twice over). Recomputed by
        // `binomial_cdf` above rather than quoted:
        //   P(X <= 179 | n = 200, p = 0.95) = 1.1599e-3   <- this floor's false-failure rate
        // the same order as truing_uncertainty.rs's ">= 33 of 40" precedent (7.115e-4).
        //
        // NOTE, deliberate divergence from the plan: the plan sketched "P(covered <= 183) ~
        // 8.4e-4" and a floor of 184. The exact tail at 183 is 2.3799e-2 — 28x the quoted
        // figure, i.e. a 1-in-42 false-failure rate at nominal, well outside the cited
        // precedent. 8.4e-4 actually lands between k=178 (4.8111e-4) and k=179 (1.1599e-3), so
        // it is the tail belonging to a floor of 179-180, not 184. The computed value governs.
        //
        // What this floor does and does not catch, measured by mutating `bounds()` at this
        // seed rather than assumed: dropping the `ln B(S+1, n-S+1)` prior term collapses
        // coverage to 26/200 and fails loudly. But substituting a repeated Wilson interval
        // scores 185/200 and swapping alpha for the confidence in the threshold scores
        // 192/200 — both pass, at this floor *and* at the plan's 184. That is a property of
        // the stopping rule, not of the floor: stopping at a fixed half-width of 0.08 fixes
        // the final interval width, so every construction is judged at whatever `n` it needs
        // to reach that width, which equalizes coverage. Those two failure modes are pinned
        // elsewhere, and each has exactly ONE robust detector — they are not redundantly
        // covered, so neither leg is safe to prune:
        //   * repeated Wilson  -> `cs_is_wider_than_wilson_at_the_same_n`, definitionally: the
        //     two widths become bit-identical, so its strict `>` cannot hold. (The nesting
        //     sweep does also fail on this mutant today, but only through a 1-ULP artifact of
        //     Wilson's k=n upper bound landing on 0.9999999999999999 — an accident, not a
        //     property. Do not rely on it.)
        //   * alpha/confidence swap -> `cs_brackets_the_mle_and_nests_by_confidence_level`,
        //     structurally: the threshold ordering inverts, so P90 comes out wider than P95 in
        //     all 1890 swept cells. The width test catches this one at only 1 of its 3
        //     fixtures — (10, 40), by 1.2% — because at (81, 263) and (500, 1000) the
        //     alpha-swapped interval is still 1.19x and 1.31x wider than Wilson. Dropping the
        //     (10, 40) fixture would silently remove that leg.
        // Read this test as "the sequence does not undercover under optional stopping", not as
        // the sole guard on the construction.
        const COVERAGE_FLOOR: usize = 180;
        let false_failure_rate = binomial_cdf(COVERAGE_FLOOR - 1, TRIALS, 0.95);
        assert!(
            (false_failure_rate - 1.1599e-3).abs() < 1e-6,
            "exact tail P(X <= {}) = {false_failure_rate:.6e}, expected ~1.1599e-3",
            COVERAGE_FLOOR - 1
        );
        assert!(
            covered >= COVERAGE_FLOOR,
            "coverage {covered}/{TRIALS} below the exact-tail floor {COVERAGE_FLOOR}"
        );
        // Anti-width guard: mean final half-width must show real convergence.
        assert!(mean_hw < 0.12, "mean final half-width {mean_hw} — intervals not converging");
    }
}
