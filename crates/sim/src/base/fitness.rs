use std::fmt;
use std::ops::{Add, AddAssign, Mul, MulAssign, Deref};
use std::iter::Sum;
// no type-state marker; simplify types
use serde::{Serialize, Deserialize};

/// A fitness value representing non-negative fitness. Values are not normalized
/// by default and may exceed 1.0. Use `.ln()` for log-space operations.
///
/// Construction/assignment from NaN is forbidden and will cause a panic to
/// prevent NaNs from entering the numeric system; this keeps invariants simple
/// and prevents hidden propagation of NaN through calculations.
///
/// # Examples
///
/// Basic construction and operations:
///
/// ```rust
/// # use centrevo_sim::base::fitness::FitnessValue;
/// let f1 = FitnessValue::new(0.2);
/// let f2 = FitnessValue::new(0.3);
/// let sum = f1 + f2; // addition
/// let sum_val: f64 = sum.into();
/// assert!((sum_val - 0.5).abs() < 1e-12);
/// ```
///
/// Manual normalization by total weight
///
/// ```rust
/// # use centrevo_sim::base::fitness::FitnessValue;
/// let w1 = FitnessValue::new(2.0);
/// let w2 = FitnessValue::new(3.0);
/// let total = [w1, w2].into_iter().sum::<FitnessValue>();
/// let norm_w1: f64 = (*w1 as f64) / (*total as f64);
/// assert!((norm_w1 - 0.4).abs() < 1e-12);
/// ```
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct FitnessValue(f64);

impl FitnessValue {
    /// Lethal fitness value constant.
    /// This represents no fitness (0.0).
    pub const LETHAL_FITNESS: Self = Self(0.0);

    /// Neutral fitness value constant.
    /// This represents neutral fitness (1.0).
    pub const NEUTRAL_FITNESS: Self = Self(1.0);

    /// Creates a new FitnessValue representing a non-negative weight.
    ///
    /// Values are in [0.0, +inf) and may act as weights for selection or scoring.
    /// 
    /// # Arguments
    /// * `value` - The fitness value (must be finite and non-negative).
    /// 
    /// # Panics
    /// * Panics if `value` is NaN or infinite.
    /// * Panics if `value` is negative.
    pub fn new(value: f64) -> Self {
        // Panic on NaN or infinite values to maintain the invariant that FitnessValue contains only finite numeric values.
        assert!(value.is_finite(), "FitnessValue must be finite");
        // Panic on negative values to maintain the invariant that FitnessValue is in [0.0, +inf)
        assert!(value >= 0.0, "FitnessValue cannot be negative");
        Self(value)
    }

    /// Converts to the log-scale fitness value.
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use centrevo_sim::base::fitness::FitnessValue;
    /// let f = FitnessValue::new(0.25);
    /// let log = f.ln();
    /// let back: f64 = log.exp().into();
    /// assert!((back - 0.25).abs() < 1e-12);
    /// ```
    pub fn ln(self) -> LogFitnessValue {
        LogFitnessValue(self.0.ln())
    }

    /// Returns true if this fitness value is lethal (0.0).
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use centrevo_sim::base::fitness::FitnessValue;
    /// let lethal = FitnessValue::LETHAL_FITNESS;
    /// assert!(lethal.is_lethal());
    /// let neutral = FitnessValue::NEUTRAL_FITNESS;
    /// assert!(!neutral.is_lethal());
    /// ```
    pub fn is_lethal(self) -> bool {
        self.0 == 0.0
    }
}

impl Deref for FitnessValue {
    type Target = f64;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl AsRef<f64> for FitnessValue {
    fn as_ref(&self) -> &f64 {
        &self.0
    }
}

impl From<FitnessValue> for f64 {
    fn from(fitness: FitnessValue) -> Self {
        fitness.0
    }
}

impl From<FitnessValue> for LogFitnessValue {
    fn from(fitness: FitnessValue) -> Self {
        fitness.ln()
    }
}

impl PartialEq for FitnessValue {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl PartialOrd for FitnessValue {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.0.partial_cmp(&other.0)
    }
}

impl fmt::Display for FitnessValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl Default for FitnessValue {
    fn default() -> Self {
        Self::NEUTRAL_FITNESS // 1.0
    }
}

impl Add for FitnessValue {
    type Output = FitnessValue;

    /// Adds two linear-scale fitness values (weight addition).
    /// Result can exceed 1.0, remains in [0.0, +inf).
    fn add(self, rhs: Self) -> Self::Output {
        FitnessValue(self.0 + rhs.0)
    }
}

impl AddAssign for FitnessValue {
    /// Adds and assigns two linear-scale fitness values.
    fn add_assign(&mut self, rhs: Self) {
        *self = Self(self.0 + rhs.0);
    }
}

impl Sum for FitnessValue {
    /// Sums an iterator of linear-scale FitnessValues.
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let total: f64 = iter.map(|f| f.0).sum();
        Self(total)
    }
}

impl<'a> Sum<&'a FitnessValue> for FitnessValue {
    /// Sums an iterator of references to linear-scale FitnessValues.
    fn sum<I: Iterator<Item = &'a FitnessValue>>(iter: I) -> Self {
        let total: f64 = iter.map(|f| f.0).sum();
        Self(total)
    }
}

/// Multiplication for `FitnessValue` uses log-space for numerical stability.
///
/// Multiplication is implemented via `ln(a) + ln(b)` to prevent underflow
/// when multiplying very small values. The implementation also contains
/// fast-paths for `0.0` and `1.0` to avoid unnecessary conversions.
///
/// # Examples
///
/// ```rust
/// # use centrevo_sim::base::fitness::FitnessValue;
/// let f1 = FitnessValue::new(0.5);
/// let f2 = FitnessValue::new(0.5);
/// let res = f1 * f2;
/// assert!((*res - 0.25).abs() < 1e-12);
/// ```
impl Mul for FitnessValue {
    type Output = Self;

    /// Multiplies two linear-scale fitness values using log-space addition for numerical stability.
    ///
    /// Converts to log scale, adds, then converts back: a × b = exp(ln(a) + ln(b))
    /// This helps prevent underflow when multiplying very small fitness values.
    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        // Fast path for common cases
        match (self.0, rhs.0) {
            (1.0, _) => return rhs,
            (_, 1.0) => return self,
            (0.0, _) | (_, 0.0) => return Self(0.0),
            _ => {}  // fall through to log-space multiplication
        }
        let log_self = LogFitnessValue::from(self);
        let log_rhs = LogFitnessValue::from(rhs);
        let product = log_self.add(log_rhs).exp();
        Self(product.0)
    }
}

impl MulAssign for FitnessValue {
    /// Multiplies and assigns using log-space multiplication for numerical stability.
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        // Use log-space multiplication to avoid underflow, consistent with `Mul`.
        match (self.0, rhs.0) {
            (1.0, _) => *self = rhs,
            (_, 1.0) => {},
            (0.0, _) | (_, 0.0) => *self = FitnessValue(0.0),
            _ => {
                let log_self = LogFitnessValue(self.0.ln());
                let log_rhs = LogFitnessValue(rhs.0.ln());
                *self = log_self.add(log_rhs).exp();
            }
        }
    }
}

/// A log-scale fitness value (natural logarithm of fitness).
///
/// Note: Construction/assignment from NaN is forbidden and will cause a panic.
///
/// # Examples
///
/// ```rust
/// # use centrevo_sim::base::fitness::LogFitnessValue;
/// let log_val = LogFitnessValue::new(-std::f64::consts::LN_2); // ln(0.5)
/// let lin = log_val.exp();
/// assert!((*lin - 0.5).abs() < 1e-12);
/// ```
/// This is useful for numerical stability when working with very small fitness values,
/// as it avoids underflow issues. The value represents ln(fitness).
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct LogFitnessValue(f64);

impl LogFitnessValue {
    /// Creates a new LogFitnessValue from a log-scale value.
    /// Since fitness is typically in [0.0, 1.0], log-scale values are typically in [-∞, 0.0].
    /// 
    /// # Arguments
    /// * `log_value` - The natural logarithm of the fitness value.
    /// 
    /// # Panics
    /// * Panics if `log_value` is NaN.
    /// * Panics if `log_value` is positive infinity or non-finite (only negative infinity is allowed).
    pub fn new(log_value: f64) -> Self {
        // Logs must not be NaN. We allow negative infinity to represent zero
        // fitness, but disallow positive infinity or other non-finite values.
        assert!(!log_value.is_nan(), "LogFitnessValue cannot be NaN");
        assert!(log_value.is_finite() || log_value == f64::NEG_INFINITY,
                "LogFitnessValue must be finite or NEG_INFINITY");
        Self(log_value)
    }
}

impl LogFitnessValue {
    /// Lethal log fitness value constant.
    /// This represents no fitness (log(0.0) = -∞).
    pub const LETHAL_FITNESS: Self = Self(f64::NEG_INFINITY);

    /// Neutral log fitness value constant.
    /// This represents neutral fitness (log(1.0) = 0.0).
    pub const NEUTRAL_FITNESS: Self = Self(0.0);

    /// Exponentiates to get the linear-scale fitness.
    pub fn exp(self) -> FitnessValue {
        FitnessValue(self.0.exp())
    }

    /// Returns true if this represents lethal fitness (log = -∞).
    ///
    /// # Examples
    /// ```rust
    /// # use centrevo_sim::base::fitness::LogFitnessValue;
    /// let lethal_log = LogFitnessValue::LETHAL_FITNESS;
    /// assert!(lethal_log.is_lethal());
    /// let neutral_log = LogFitnessValue::NEUTRAL_FITNESS;
    /// assert!(!neutral_log.is_lethal());
    /// ```
    pub fn is_lethal(self) -> bool {
        self.0.is_infinite() && self.0.is_sign_negative()
    }
}

impl Deref for LogFitnessValue {
    type Target = f64;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl AsRef<f64> for LogFitnessValue {
    fn as_ref(&self) -> &f64 {
        &self.0
    }
}

impl From<LogFitnessValue> for f64 {
    fn from(log_fitness: LogFitnessValue) -> Self {
        log_fitness.0
    }
}

impl From<LogFitnessValue> for FitnessValue {
    fn from(log_fitness: LogFitnessValue) -> Self {
        log_fitness.exp()
    }
}

impl fmt::Display for LogFitnessValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl From<f64> for LogFitnessValue {
    fn from(value: f64) -> Self {
        Self::new(value)
    }
}


impl Default for LogFitnessValue {
    fn default() -> Self {
        Self::NEUTRAL_FITNESS // ln(1.0) = 0.0
    }
}

impl Add for LogFitnessValue {
    type Output = Self;

    /// Adds two log fitness values (equivalent to multiplying linear fitness values).
    ///
    /// In log space: ln(a × b) = ln(a) + ln(b)
    /// Correctly handles -∞: if either fitness is zero, the result is zero.
    fn add(self, other: Self) -> Self::Output {
        Self(self.0 + other.0)
    }
}

impl AddAssign for LogFitnessValue {
    /// Adds and assigns two log fitness values.
    fn add_assign(&mut self, rhs: Self) {
        *self = Self(self.0 + rhs.0);
    }
}

impl Sum for LogFitnessValue {
    /// Sums an iterator of log FitnessValues.
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let total_log: f64 = iter.map(|f| f.0).sum();
        Self(total_log)
    }
}

impl<'a> Sum<&'a LogFitnessValue> for LogFitnessValue {
    /// Sums an iterator of references to log FitnessValues.
    fn sum<I: Iterator<Item = &'a LogFitnessValue>>(iter: I) -> Self {
        let total_log: f64 = iter.map(|f| f.0).sum();
        Self(total_log)
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    /// Helper for comparing f64 values with relative tolerance.
    fn approx_eq(a: f64, b: f64, eps: f64) -> bool {
        if a.is_nan() || b.is_nan() {
            return false;
        }
        let diff = (a - b).abs();
        if a == 0.0 || b == 0.0 {
            return diff < eps;
        }
        diff / a.abs().max(b.abs()) < eps
    }

    // #[test]
    // fn test_new_clamps_negative_to_zero() {
    //     let f = FitnessValue::new(-1.0);
    //     assert!(approx_eq(f.get(), 0.0, 1e-12));
    // }

    #[test]
    fn test_new_preserves_zero() {
        let f = FitnessValue::new(0.0);
        assert!(approx_eq(*f, 0.0, 1e-12));
    }

    #[test]
    fn test_new_preserves_midrange() {
        let f = FitnessValue::new(0.5);
        assert!(approx_eq(*f, 0.5, 1e-12));
    }

    #[test]
    fn test_new_preserves_one() {
        let f = FitnessValue::new(1.0);
        assert!(approx_eq(*f, 1.0, 1e-12));
    }

    // #[test]
    // fn test_new_clamps_above_one_to_one() {
    //     let f = FitnessValue::new(2.0);
    //     assert!(approx_eq(f.get(), 1.0, 1e-12));
    // }

    // #[test]
    // fn test_from_f64_clamps_and_preserves_values() {
    //     let f_from_neg: FitnessValue = (-1.0).into();
    //     assert!(approx_eq(f_from_neg.get(), 0.0, 1e-12));

    //     let f_from_pos: FitnessValue = 0.75.into();
    //     assert!(approx_eq(f_from_pos.get(), 0.75, 1e-12));

    //     let f_from_big: FitnessValue = 10.0.into();
    //     assert!(approx_eq(f_from_big.get(), 1.0, 1e-12));
    // }

    #[test]
    #[should_panic(expected = "FitnessValue must be finite")]
    fn test_new_nan_panics() {
        let nan_val = f64::NAN;
        let _nan_f = FitnessValue::new(nan_val);
    }

    #[test]
    #[should_panic(expected = "FitnessValue must be finite")]
    fn test_new_infinite_panics() {
        let _ = FitnessValue::new(f64::INFINITY);
    }

    #[test]
    fn test_ln_and_exp_roundtrip_zero() {
        let input = 0.0;
        let f= FitnessValue::new(input);
        let log = f.ln();
        let log_val = *log;
        let back = log.exp();
        assert!(log_val.is_infinite() && log_val.is_sign_negative());
        assert!(*back == 0.0);
    }

    #[test]
    fn test_ln_and_exp_roundtrip_one() {
        let input = 1.0;
        let f = FitnessValue::new(input);
        let log = f.ln();
        let back = log.exp();
        assert!(approx_eq(input, *back, 1e-12));
    }

    #[test]
    fn test_ln_and_exp_roundtrip_mid() {
        let input = 0.5;
        let f = FitnessValue::new(input);
        let log = f.ln();
        let back = log.exp();
        assert!(approx_eq(input, *back, 1e-12));
    }

    #[test]
    fn test_ln_and_exp_roundtrip_very_small() {
        let input = 1e-308;
        let f = FitnessValue::new(input);
        let log = f.ln();
        let back = log.exp();
        assert!(approx_eq(input, *back, 1e-12));
    }

    #[test]
    #[should_panic(expected = "FitnessValue must be finite")]
    fn test_ln_and_exp_roundtrip_nan() {
        // Creating a FitnessValue from NaN should panic; therefore this test checks
        // that such construction triggers a panic as part of the new NaN handling.
        let _nan_f= FitnessValue::new(f64::NAN);
    }

    #[test]
    fn test_log_exp_of_ln_half() {
        let log_val = LogFitnessValue::new(-std::f64::consts::LN_2); // ln(0.5)
        let lin = log_val.exp();
        assert!(approx_eq(*lin, 0.5, 1e-12));
    }

    #[test]
    fn test_log_from_f64_preserves_value() {
        let log_from_f: LogFitnessValue = (-1.234).into();
        assert!(approx_eq(*log_from_f, -1.234, 1e-12));
    }

    #[test]
    fn test_conversion_roundtrip_for_quarter() {
        let f = FitnessValue::new(0.25);
        let log = LogFitnessValue::from(f);
        let back: FitnessValue = log.into();
        assert!(approx_eq(0.25, *back, 1e-12));
    }

    #[test]
    fn test_defaults_are_correct() {
        let default_f: FitnessValue = Default::default();
        assert!(approx_eq(*default_f, 1.0, 1e-12));

        let default_f: FitnessValue = Default::default();
        assert!(approx_eq(*default_f, 1.0, 1e-12));

        let default_log: LogFitnessValue = Default::default();
        assert!(approx_eq(*default_log, 0.0, 1e-12));

        let default_log: LogFitnessValue = Default::default();
        assert!(approx_eq(*default_log, 0.0, 1e-12));
    }

    #[test]
    fn test_display_parsable_for_default() {
        let default_f: FitnessValue = Default::default();
        let disp = default_f.to_string();
        let parsed: f64 = disp.parse().expect("display should be parsable as f64");
        assert!(approx_eq(parsed, 1.0, 1e-12));
    }

    // #[test]
    // fn test_log_new_clamps_positive_to_zero() {
    //     let log_val = LogFitnessValue::new(5.0);
    //     assert!(approx_eq(log_val.get(), 0.0, 1e-12));
    // }

    #[test]
    fn test_log_new_preserves_negative() {
        let log_val = LogFitnessValue::new(-2.5);
        assert!(approx_eq(*log_val, -2.5, 1e-12));
    }

    #[test]
    fn test_log_new_preserves_neg_infinity() {
        let log = LogFitnessValue::new(f64::NEG_INFINITY);
        let log_val = *log;
        assert!(log_val.is_infinite() && log_val.is_sign_negative());
    }

    #[test]
    #[should_panic(expected = "LogFitnessValue must be finite or NEG_INFINITY")]
    fn test_log_new_positive_infinite_panics() {
        let _ = LogFitnessValue::new(f64::INFINITY);
    }

    #[test]
    fn test_is_zero_fitness_detects_neg_infinity() {
        let zero_log = LogFitnessValue::new(f64::NEG_INFINITY);
        assert!(zero_log.is_lethal());

        let nonzero_log = LogFitnessValue::new(-1.0);
        assert!(!nonzero_log.is_lethal());
    }

    #[test]
    fn test_add_combines_log_fitness() {
        let log1 = LogFitnessValue::new(-std::f64::consts::LN_2); // ln(0.5)
        let log2 = LogFitnessValue::new(-std::f64::consts::LN_2); // ln(0.5)
        let combined = log1.add(log2);
        // ln(0.5) + ln(0.5) = ln(0.25)
        let result = combined.exp();
        assert!(approx_eq(*result, 0.25, 1e-12));
    }

    #[test]
    fn test_add_with_zero_fitness_gives_zero() {
        let zero_log = LogFitnessValue::new(f64::NEG_INFINITY);
        let normal_log = LogFitnessValue::new(-1.0);

        let result1 = zero_log.add(normal_log);
        assert!(result1.is_lethal());

        let result2 = normal_log.add(zero_log);
        assert!(result2.is_lethal());
    }

    #[test]
    fn test_mul_fitness_values() {
        let f1 = FitnessValue::new(0.5);
        let f2 = FitnessValue::new(0.5);
        let result = f1 * f2;
        assert!(approx_eq(*result, 0.25, 1e-12));
    }

    #[test]
    fn test_mul_with_one_preserves_value() {
        let f1 = FitnessValue::new(0.7);
        let f2 = FitnessValue::new(1.0);
        let result = f1 * f2;
        assert!(approx_eq(*result, 0.7, 1e-12));
    }

    #[test]
    fn test_mul_with_zero_gives_zero() {
        let f1 = FitnessValue::new(0.8);
        let f2 = FitnessValue::new(0.0);
        let result = f1 * f2;
        assert!(approx_eq(*result, 0.0, 1e-12));
    }

    #[test]
    fn test_mul_chain_three_values() {
        let f1 = FitnessValue::new(0.5);
        let f2 = FitnessValue::new(0.5);
        let f3 = FitnessValue::new(0.5);
        let result = f1 * f2 * f3;
        assert!(approx_eq(*result, 0.125, 1e-12));
    }

    #[test]
    fn test_mul_very_small_values() {
        let f1 = FitnessValue::new(1e-100);
        let f2 = FitnessValue::new(1e-100);
        let result = f1 * f2;
        // Direct multiplication would underflow, but log-space handles it
        assert!(*result > 0.0);
        assert!(approx_eq(*result, 1e-200, 1e-12));
    }

    #[test]
    fn test_add_fitness_values() {
        let f1 = FitnessValue::new(0.2);
        let f2 = FitnessValue::new(0.3);
        let result = f1 + f2;
        assert!(approx_eq(*result, 0.5, 1e-12));
    }

    // #[test]
    // fn test_add_fitness_values_clamp() {
    //     let f1 = FitnessValue::new(0.7);
    //     let f2 = FitnessValue::new(0.5);
    //     let result = f1 + f2;
    //     assert!(approx_eq(result.get(), 1.0, 1e-12));
    // }

    #[test]
    fn test_add_log_values() {
        let log1 = LogFitnessValue::new(-std::f64::consts::LN_2); // ln(0.5)
        let log2 = LogFitnessValue::new(-std::f64::consts::LN_2); // ln(0.5)
        let combined = log1 + log2;
        let result = combined.exp();
        assert!(approx_eq(*result, 0.25, 1e-12));
    }

    #[test]
    fn test_add_log_with_zero_gives_zero() {
        let zero_log = LogFitnessValue::new(f64::NEG_INFINITY);
        let normal_log = LogFitnessValue::new(-1.0);
        let result1 = zero_log + normal_log;
        assert!(result1.is_lethal());
        let result2 = normal_log + zero_log;
        assert!(result2.is_lethal());
    }

    #[test]
    fn test_add_assign_basic_sum() {
        let mut f = FitnessValue::new(0.2);
        f += FitnessValue::new(0.3);
        assert!(approx_eq(*f, 0.5, 1e-12));
    }

    // #[test]
    // fn test_add_assign_clamp() {
    //     let mut f = FitnessValue::new(0.7);
    //     f += FitnessValue::new(0.5);
    //     assert!(approx_eq(f.get(), 1.0, 1e-12));
    // }

    #[test]
    fn test_add_assign_with_zero() {
        let mut f = FitnessValue::new(0.0);
        f += FitnessValue::new(0.5);
        assert!(approx_eq(*f, 0.5, 1e-12));

        let mut g = FitnessValue::new(0.5);
        g += FitnessValue::new(0.0);
        assert!(approx_eq(*g, 0.5, 1e-12));
    }

    #[test]
    #[should_panic(expected = "FitnessValue must be finite")]
    fn test_add_assign_with_nan() {
        // Constructing a FitnessValue with NaN should panic; add-assign tries to
        // create such a value (by converting from NaN) so it should panic at the
        // point of construction.
        let _ = FitnessValue::new(f64::NAN);
    }

    #[test]
    fn test_add_assign_chained() {
        let mut f = FitnessValue::new(0.2);
        f += FitnessValue::new(0.3);
        f += FitnessValue::new(0.4);
        assert!(approx_eq(*f, 0.9, 1e-12));
    }

    // MulAssign tests
    #[test]
    fn test_mul_assign_basic() {
        let mut f = FitnessValue::new(0.5);
        f *= FitnessValue::new(0.5);
        assert!(approx_eq(*f, 0.25, 1e-12));
    }

    #[test]
    fn test_mul_assign_one_preserves_value() {
        let mut f = FitnessValue::new(0.7);
        f *= FitnessValue::new(1.0);
        assert!(approx_eq(*f, 0.7, 1e-12));

        let mut g = FitnessValue::new(1.0);
        g *= FitnessValue::new(0.7);
        assert!(approx_eq(*g, 0.7, 1e-12));
    }

    #[test]
    fn test_mul_assign_zero_gives_zero() {
        let mut f = FitnessValue::new(0.8);
        f *= FitnessValue::new(0.0);
        assert!(approx_eq(*f, 0.0, 1e-12));

        let mut g = FitnessValue::new(0.0);
        g *= FitnessValue::new(0.8);
        assert!(approx_eq(*g, 0.0, 1e-12));
    }

    #[test]
    fn test_mul_assign_chained() {
        let mut f = FitnessValue::new(0.5);
        f *= FitnessValue::new(0.5);
        f *= FitnessValue::new(0.5);
        assert!(approx_eq(*f, 0.125, 1e-12));
    }

    #[test]
    fn test_mul_assign_very_small_values() {
        let mut f = FitnessValue::new(1e-100);
        f *= FitnessValue::new(1e-100);
        assert!(*f > 0.0);
        assert!(approx_eq(*f, 1e-200, 1e-12));
    }
}
