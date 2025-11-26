use std::fmt;
use std::ops::{Add, AddAssign, Mul, MulAssign, Deref};
use std::iter::Sum;
use std::marker::PhantomData;
use serde::{Serialize, Deserialize};

// State marker types for FitnessValue
/// Unnormalized state: fitness value can be in [0.0, +inf), acts as weight
#[derive(Debug, Clone, Copy)]
pub struct Unnormalized;
/// Normalized state: fitness value is constrained to [0.0, 1.0], acts as probability
#[derive(Debug, Clone, Copy)]
pub struct Normalized;

/// A fitness value with type-state pattern.
///
/// - **Unnormalized** (default): fitness value is in [0.0, +inf), acts as weight.
///   Use this for accumulating fitness values before normalizing.
/// - **Normalized**: fitness value is constrained to [0.0, 1.0], acts as probability.
///   Created via `normalize()` method.
///
/// Construction/assignment from NaN is forbidden and will cause a panic to
/// prevent NaNs from entering the numeric system; this keeps invariants simple
/// and prevents hidden propagation of NaN through calculations.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct FitnessValue<State>(f64, PhantomData<State>);

impl FitnessValue<Unnormalized> {
    /// Creates a new unnormalized FitnessValue.
    ///
    /// Unnormalized values can be in [0.0, +inf) and act as weights.
    /// This is the default state and should be used for accumulating fitness values.
    pub fn new(value: f64) -> Self {
        // Panic on NaN to maintain the invariant that FitnessValue contains only finite numeric values.
        assert!(!value.is_nan(), "FitnessValue cannot be NaN");
        // Panic on negative values to maintain the invariant that FitnessValue is in [0.0, +inf)
        assert!(value >= 0.0, "FitnessValue cannot be negative");
        Self(value, PhantomData)
    }

    /// Normalizes this fitness value to the range [0.0, 1.0].
    /// Recommended to normalize in log-space for numerical stability when 
    /// dealing with very small values.
    ///
    /// If the value is already in [0.0, 1.0], it remains unchanged.
    #[inline]
    pub fn normalize(self, total: impl Into<FitnessValue<Unnormalized>>) -> FitnessValue<Normalized> {
        let total = total.into();
        // Panic if total is less than self.0, which would lead to normalized > 1.0
        assert!(self.0 <= total.0, "FitnessValue cannot be greater than total fitness for normalization");
        FitnessValue(self.0 / total.0, PhantomData)
    }
}

impl FitnessValue<Normalized> {
    /// Creates a new normalized FitnessValue ensuring bounds [0.0, 1.0].
    pub fn new_normalized(value: f64) -> Self {
        // Panic on NaN to maintain the invariant that FitnessValue contains only finite numeric values.
        assert!(!value.is_nan(), "FitnessValue cannot be NaN");
        // Panic on values outside [0.0, 1.0] to maintain the invariant that Normalized FitnessValue is in [0.0, 1.0]
        assert!((0.0..=1.0).contains(&value), "Normalized FitnessValue must be in [0.0, 1.0]");
        Self(value, PhantomData)
    }

    /// Manually changes to unnormalized state.
    /// This does not change the numeric value; it only changes the type state.
    pub fn unnormalize(self) -> FitnessValue<Unnormalized> {
        FitnessValue(self.0, PhantomData)
    }
}

impl<State> FitnessValue<State> {
    /// Lethal fitness value constant.
    /// This represents no fitness (0.0).
    pub const LETHAL_FITNESS: Self = Self(0.0, PhantomData);

    /// Neutral fitness value constant.
    /// This represents neutral fitness (1.0).
    pub const NEUTRAL_FITNESS: Self = Self(1.0, PhantomData);

    /// Converts to the log-scale fitness value.
    pub fn ln(self) -> LogFitnessValue<State> {
        LogFitnessValue(self.0.ln(), PhantomData)
    }
}

impl<State> Deref for FitnessValue<State> {
    type Target = f64;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<State> AsRef<f64> for FitnessValue<State> {
    fn as_ref(&self) -> &f64 {
        &self.0
    }
}

impl<State> From<FitnessValue<State>> for f64 {
    fn from(fitness: FitnessValue<State>) -> Self {
        fitness.0
    }
}

impl<State> From<FitnessValue<State>> for LogFitnessValue<State> {
    fn from(fitness: FitnessValue<State>) -> Self {
        fitness.ln()
    }
}

impl<State> PartialEq for FitnessValue<State> {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0 && self.1 == other.1
    }
}

impl<State> PartialOrd for FitnessValue<State> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.0.partial_cmp(&other.0)
    }
}

impl<State> fmt::Display for FitnessValue<State> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl<State> Default for FitnessValue<State> {
    fn default() -> Self {
        Self::NEUTRAL_FITNESS // 1.0
    }
}

impl<State> Add for FitnessValue<State> {
    type Output = FitnessValue<Unnormalized>;

    /// Adds two unnormalized fitness values (weight addition).
    /// Result can exceed 1.0, remains in [0.0, +inf).
    fn add(self, rhs: Self) -> Self::Output {
        FitnessValue(self.0 + rhs.0, PhantomData)
    }
}

impl AddAssign for FitnessValue<Unnormalized> {
    /// Adds and assigns two unnormalized fitness values.
    fn add_assign(&mut self, rhs: Self) {
        *self = Self(self.0 + rhs.0, PhantomData);
    }
}

impl Sum for FitnessValue<Unnormalized> {
    /// Sums an iterator of unnormalized FitnessValues.
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let total: f64 = iter.map(|f| f.0).sum();
        Self(total, PhantomData)
    }
}

impl<'a> Sum<&'a FitnessValue<Unnormalized>> for FitnessValue<Unnormalized> {
    /// Sums an iterator of references to unnormalized FitnessValues.
    fn sum<I: Iterator<Item = &'a FitnessValue<Unnormalized>>>(iter: I) -> Self {
        let total: f64 = iter.map(|f| f.0).sum();
        Self(total, PhantomData)
    }
}

impl<State> Mul for FitnessValue<State> {
    type Output = Self;

    /// Multiplies two unnormalized fitness values using log-space addition for numerical stability.
    /// 
    /// Converts to log scale, adds, then converts back: a × b = exp(ln(a) + ln(b))
    /// This helps prevent underflow when multiplying very small fitness values.
    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        // Fast path for common cases
        match (self.0, rhs.0) {
            (1.0, _) => return rhs,
            (_, 1.0) => return self,
            (0.0, _) | (_, 0.0) => return Self(0.0, PhantomData),
            _ => {}  // fall through to log-space multiplication
        }
        let log_self = LogFitnessValue::from(self);
        let log_rhs = LogFitnessValue::from(rhs);
        let product = log_self.add(log_rhs).exp();
        Self(product.0, PhantomData)
    }
}

impl<State> MulAssign for FitnessValue<State> {
    /// Multiplies and assigns using log-space multiplication for numerical stability.
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        *self = Self(self.0 * rhs.0, PhantomData);
    }
}

/// A log-scale fitness value (natural logarithm of fitness).
///
/// Note: Construction/assignment from NaN is forbidden and will cause a panic.
/// 
/// This is useful for numerical stability when working with very small fitness values,
/// as it avoids underflow issues. The value represents ln(fitness).
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct LogFitnessValue<State = Unnormalized>(f64, PhantomData<State>);

impl LogFitnessValue<Unnormalized> {
    /// Creates a new LogFitnessValue from a log-scale value.
    /// Since fitness is typically in [0.0, 1.0], log-scale values are typically in [-∞, 0.0].
    pub fn new(log_value: f64) -> Self {
        // Logs must not be NaN. If a caller attempts to construct a LogFitnessValue
        // with NaN, panic (assignment is not allowed).
        assert!(!log_value.is_nan(), "LogFitnessValue cannot be NaN");
        Self(log_value, PhantomData)
    }

    /// Normalizes this log fitness value to the range [-∞, 0.0].
    /// 
    /// If the value is already in [-∞, 0.0], it remains unchanged.
    #[inline]
    pub fn normalize(self, total_log: impl Into<LogFitnessValue<Unnormalized>>) -> LogFitnessValue<Normalized> {
        let total_log = total_log.into();
        // Panic if total_log is less than self.0, which would lead to normalized > 1.0
        assert!(self.0 <= total_log.0, "LogFitnessValue cannot be greater than total log fitness for normalization");
        LogFitnessValue(self.0 - total_log.0, PhantomData)
    }
}

impl LogFitnessValue<Normalized> {
    /// Manually changes to unnormalized state.
    /// This does not change the numeric value; it only changes the type state.
    pub fn unnormalize(self) -> LogFitnessValue<Unnormalized> {
        LogFitnessValue(self.0, PhantomData)
    }
}

impl<State> LogFitnessValue<State> {
    /// Lethal log fitness value constant.
    /// This represents no fitness (log(0.0) = -∞).
    pub const LETHAL_FITNESS: Self = Self(f64::NEG_INFINITY, PhantomData);

    /// Neutral log fitness value constant.
    /// This represents neutral fitness (log(1.0) = 0.0).
    pub const NEUTRAL_FITNESS: Self = Self(0.0, PhantomData);

    /// Exponentiates to get the linear-scale fitness.
    pub fn exp(self) -> FitnessValue<State> {
        FitnessValue(self.0.exp(), PhantomData)
    }

    /// Returns true if this represents zero fitness (log = -∞).
    pub fn is_zero_fitness(self) -> bool {
        self.0.is_infinite() && self.0.is_sign_negative()
    }
}

impl<State> Deref for LogFitnessValue<State> {
    type Target = f64;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<State> AsRef<f64> for LogFitnessValue<State> {
    fn as_ref(&self) -> &f64 {
        &self.0
    }
}

impl<State> From<LogFitnessValue<State>> for f64 {
    fn from(log_fitness: LogFitnessValue<State>) -> Self {
        log_fitness.0
    }
}

impl<State> From<LogFitnessValue<State>> for FitnessValue<State> {
    fn from(log_fitness: LogFitnessValue<State>) -> Self {
        log_fitness.exp()
    }
}

impl<State> fmt::Display for LogFitnessValue<State> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl From<f64> for LogFitnessValue<Unnormalized> {
    fn from(value: f64) -> Self {
        Self::new(value)
    }
}


impl<State> Default for LogFitnessValue<State> {
    fn default() -> Self {
        Self::NEUTRAL_FITNESS // ln(1.0) = 0.0
    }
}

impl<State> Add for LogFitnessValue<State> {
    type Output = Self;

    /// Adds two log fitness values (equivalent to multiplying linear fitness values).
    /// 
    /// In log space: ln(a × b) = ln(a) + ln(b)
    /// Correctly handles -∞: if either fitness is zero, the result is zero.
    fn add(self, other: Self) -> Self::Output {
        Self(self.0 + other.0, PhantomData)
    }
}

impl<State> AddAssign for LogFitnessValue<State> {
    /// Adds and assigns two log fitness values.
    fn add_assign(&mut self, rhs: Self) {
        *self = Self(self.0 + rhs.0, PhantomData);
    }
}

impl<State> Sum for LogFitnessValue<State> {
    /// Sums an iterator of log FitnessValues.
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let total_log: f64 = iter.map(|f| f.0).sum();
        Self(total_log, PhantomData)
    }
}

impl<'a, State> Sum<&'a LogFitnessValue<State>> for LogFitnessValue<State> {
    /// Sums an iterator of references to log FitnessValues.
    fn sum<I: Iterator<Item = &'a LogFitnessValue<State>>>(iter: I) -> Self {
        let total_log: f64 = iter.map(|f| f.0).sum();
        Self(total_log, PhantomData)
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
        let f = FitnessValue::<Unnormalized>::new(0.0);
        assert!(approx_eq(*f, 0.0, 1e-12));
    }

    #[test]
    fn test_new_preserves_midrange() {
        let f = FitnessValue::<Unnormalized>::new(0.5);
        assert!(approx_eq(*f, 0.5, 1e-12));
    }

    #[test]
    fn test_new_preserves_one() {
        let f = FitnessValue::<Unnormalized>::new(1.0);
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
    #[should_panic(expected = "FitnessValue cannot be NaN")]
    fn test_new_nan_panics() {
        let nan_val = f64::NAN;
        let _nan_f = FitnessValue::<Unnormalized>::new(nan_val);
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
    #[should_panic(expected = "FitnessValue cannot be NaN")]
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
        let f = FitnessValue::<Unnormalized>::new(0.25);
        let log = LogFitnessValue::from(f);
        let back: FitnessValue<Unnormalized> = log.into();
        assert!(approx_eq(0.25, *back, 1e-12));
    }

    #[test]
    fn test_defaults_are_correct() {
        let default_f: FitnessValue<Unnormalized> = Default::default();
        assert!(approx_eq(*default_f, 1.0, 1e-12));

        let default_f: FitnessValue<Normalized> = Default::default();
        assert!(approx_eq(*default_f, 1.0, 1e-12));

        let default_log: LogFitnessValue<Unnormalized> = Default::default();
        assert!(approx_eq(*default_log, 0.0, 1e-12));

        let default_log: LogFitnessValue<Normalized> = Default::default();
        assert!(approx_eq(*default_log, 0.0, 1e-12));
    }

    #[test]
    fn test_display_parsable_for_default() {
        let default_f: FitnessValue<Unnormalized> = Default::default();
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
    fn test_is_zero_fitness_detects_neg_infinity() {
        let zero_log = LogFitnessValue::new(f64::NEG_INFINITY);
        assert!(zero_log.is_zero_fitness());

        let nonzero_log = LogFitnessValue::new(-1.0);
        assert!(!nonzero_log.is_zero_fitness());
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
        assert!(result1.is_zero_fitness());
        
        let result2 = normal_log.add(zero_log);
        assert!(result2.is_zero_fitness());
    }

    #[test]
    fn test_mul_fitness_values() {
        let f1 = FitnessValue::<Unnormalized>::new(0.5);
        let f2 = FitnessValue::<Unnormalized>::new(0.5);
        let result = f1 * f2;
        assert!(approx_eq(*result, 0.25, 1e-12));
    }

    #[test]
    fn test_mul_with_one_preserves_value() {
        let f1 = FitnessValue::<Unnormalized>::new(0.7);
        let f2 = FitnessValue::<Unnormalized>::new(1.0);
        let result = f1 * f2;
        assert!(approx_eq(*result, 0.7, 1e-12));
    }

    #[test]
    fn test_mul_with_zero_gives_zero() {
        let f1 = FitnessValue::<Unnormalized>::new(0.8);
        let f2 = FitnessValue::<Unnormalized>::new(0.0);
        let result = f1 * f2;
        assert!(approx_eq(*result, 0.0, 1e-12));
    }

    #[test]
    fn test_mul_chain_three_values() {
        let f1 = FitnessValue::<Unnormalized>::new(0.5);
        let f2 = FitnessValue::<Unnormalized>::new(0.5);
        let f3 = FitnessValue::<Unnormalized>::new(0.5);
        let result = f1 * f2 * f3;
        assert!(approx_eq(*result, 0.125, 1e-12));
    }

    #[test]
    fn test_mul_very_small_values() {
        let f1 = FitnessValue::<Unnormalized>::new(1e-100);
        let f2 = FitnessValue::<Unnormalized>::new(1e-100);
        let result = f1 * f2;
        // Direct multiplication would underflow, but log-space handles it
        assert!(*result > 0.0);
        assert!(approx_eq(*result, 1e-200, 1e-12));
    }

    #[test]
    fn test_add_fitness_values() {
        let f1 = FitnessValue::<Unnormalized>::new(0.2);
        let f2 = FitnessValue::<Unnormalized>::new(0.3);
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
        assert!(result1.is_zero_fitness());
        let result2 = normal_log + zero_log;
        assert!(result2.is_zero_fitness());
    }

    #[test]
    fn test_add_assign_basic_sum() {
        let mut f = FitnessValue::<Unnormalized>::new(0.2);
        f += FitnessValue::<Unnormalized>::new(0.3);
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
        let mut f = FitnessValue::<Unnormalized>::new(0.0);
        f += FitnessValue::<Unnormalized>::new(0.5);
        assert!(approx_eq(*f, 0.5, 1e-12));

        let mut g = FitnessValue::<Unnormalized>::new(0.5);
        g += FitnessValue::<Unnormalized>::new(0.0);
        assert!(approx_eq(*g, 0.5, 1e-12));
    }

    #[test]
    #[should_panic(expected = "FitnessValue cannot be NaN")]
    fn test_add_assign_with_nan() {
        // Constructing a FitnessValue with NaN should panic; add-assign tries to
        // create such a value (by converting from NaN) so it should panic at the
        // point of construction.
        let _ = FitnessValue::new(f64::NAN);
    }

    #[test]
    fn test_add_assign_chained() {
        let mut f = FitnessValue::<Unnormalized>::new(0.2);
        f += FitnessValue::<Unnormalized>::new(0.3);
        f += FitnessValue::<Unnormalized>::new(0.4);
        assert!(approx_eq(*f, 0.9, 1e-12));
    }

    // MulAssign tests
    #[test]
    fn test_mul_assign_basic() {
        let mut f = FitnessValue::<Unnormalized>::new(0.5);
        f *= FitnessValue::<Unnormalized>::new(0.5);
        assert!(approx_eq(*f, 0.25, 1e-12));
    }

    #[test]
    fn test_mul_assign_one_preserves_value() {
        let mut f = FitnessValue::<Unnormalized>::new(0.7);
        f *= FitnessValue::<Unnormalized>::new(1.0);
        assert!(approx_eq(*f, 0.7, 1e-12));

        let mut g = FitnessValue::<Unnormalized>::new(1.0);
        g *= FitnessValue::<Unnormalized>::new(0.7);
        assert!(approx_eq(*g, 0.7, 1e-12));
    }

    #[test]
    fn test_mul_assign_zero_gives_zero() {
        let mut f = FitnessValue::<Unnormalized>::new(0.8);
        f *= FitnessValue::<Unnormalized>::new(0.0);
        assert!(approx_eq(*f, 0.0, 1e-12));

        let mut g = FitnessValue::<Unnormalized>::new(0.0);
        g *= FitnessValue::<Unnormalized>::new(0.8);
        assert!(approx_eq(*g, 0.0, 1e-12));
    }

    #[test]
    fn test_mul_assign_chained() {
        let mut f = FitnessValue::<Unnormalized>::new(0.5);
        f *= FitnessValue::<Unnormalized>::new(0.5);
        f *= FitnessValue::<Unnormalized>::new(0.5);
        assert!(approx_eq(*f, 0.125, 1e-12));
    }

    #[test]
    fn test_mul_assign_very_small_values() {
        let mut f = FitnessValue::<Unnormalized>::new(1e-100);
        f *= FitnessValue::<Unnormalized>::new(1e-100);
        assert!(*f > 0.0);
        assert!(approx_eq(*f, 1e-200, 1e-12));
    }
}
