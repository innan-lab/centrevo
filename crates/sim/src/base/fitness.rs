use std::fmt;
use std::ops::{Add, Mul};

use serde::{Serialize, Deserialize};

/// A fitness value constrained to the range [0.0, 1.0].
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct FitnessValue(f64);

impl FitnessValue {
    /// Creates a new FitnessValue, clamping the input to [0.0, 1.0].
    pub fn new(value: f64) -> Self {
        Self(value.clamp(0.0, 1.0))
    }

    /// Returns the inner f64 value.
    pub fn get(self) -> f64 {
        self.0
    }

    /// Returns the inner f64 value as a reference.
    pub fn as_f64(&self) -> f64 {
        self.0
    }

    /// Converts to the log-scale fitness value.
    pub fn ln(self) -> LogFitnessValue {
        LogFitnessValue::from(self)
    }
}

impl From<FitnessValue> for f64 {
    fn from(fitness: FitnessValue) -> Self {
        fitness.0
    }
}

impl From<f64> for FitnessValue {
    fn from(value: f64) -> Self {
        Self::new(value)
    }
}

impl Default for FitnessValue {
    fn default() -> Self {
        Self(1.0)
    }
}

impl fmt::Display for FitnessValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl Add for FitnessValue {
    type Output = Self;

    /// Adds two fitness values, clamping the result to [0.0, 1.0].
    fn add(self, rhs: Self) -> Self::Output {
        FitnessValue::new(self.0 + rhs.0)
    }
}

impl Mul for FitnessValue {
    type Output = Self;

    /// Multiplies two fitness values using log-space addition for numerical stability.
    /// 
    /// Converts to log scale, adds, then converts back: a × b = exp(ln(a) + ln(b))
    fn mul(self, rhs: Self) -> Self::Output {
        let log_self = LogFitnessValue::from(self);
        let log_rhs = LogFitnessValue::from(rhs);
        log_self.add(log_rhs).exp()
    }
}

/// A log-scale fitness value (natural logarithm of fitness).
/// 
/// This is useful for numerical stability when working with very small fitness values,
/// as it avoids underflow issues. The value represents ln(fitness).
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct LogFitnessValue(f64);

impl LogFitnessValue {
    /// Creates a new LogFitnessValue from a log-scale value, clamping to [-∞, 0.0].
    /// 
    /// Since fitness is in [0.0, 1.0], the log-scale value must be in [-∞, 0.0].
    /// Values above 0.0 are clamped to 0.0.
    pub fn new(log_value: f64) -> Self {
        Self(log_value.min(0.0))
    }

    /// Returns the inner log-scale f64 value.
    pub fn get(self) -> f64 {
        self.0
    }

    /// Returns the inner log-scale f64 value as a reference.
    pub fn as_f64(&self) -> f64 {
        self.0
    }

    /// Converts to the linear-scale fitness value.
    pub fn exp(self) -> FitnessValue {
        FitnessValue::new(self.0.exp())
    }

    /// Returns true if this represents zero fitness (log = -∞).
    pub fn is_zero_fitness(self) -> bool {
        self.0.is_infinite() && self.0.is_sign_negative()
    }

    /// Adds two log fitness values (equivalent to multiplying linear fitness values).
    /// 
    /// In log space: ln(a × b) = ln(a) + ln(b)
    /// Correctly handles -∞: if either fitness is zero, the result is zero.
    pub fn add(self, other: Self) -> Self {
        Self::new(self.0 + other.0)
    }
}

impl From<LogFitnessValue> for f64 {
    fn from(log_fitness: LogFitnessValue) -> Self {
        log_fitness.0
    }
}

impl From<f64> for LogFitnessValue {
    fn from(value: f64) -> Self {
        Self::new(value)
    }
}

impl From<FitnessValue> for LogFitnessValue {
    fn from(fitness: FitnessValue) -> Self {
        Self(fitness.get().ln())
    }
}

impl From<LogFitnessValue> for FitnessValue {
    fn from(log_fitness: LogFitnessValue) -> Self {
        log_fitness.exp()
    }
}

impl Default for LogFitnessValue {
    fn default() -> Self {
        Self(0.0) // ln(1.0) = 0.0
    }
}

impl fmt::Display for LogFitnessValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl Add for LogFitnessValue {
    type Output = Self;

    /// Adds two log fitness values (equivalent to multiplying linear fitness values).
    /// 
    /// In log space: ln(a × b) = ln(a) + ln(b)
    /// Correctly handles -∞: if either fitness is zero, the result is zero.
    fn add(self, other: Self) -> Self::Output {
        Self::new(self.0 + other.0)
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

    #[test]
    fn test_new_clamps_negative_to_zero() {
        let f = FitnessValue::new(-1.0);
        assert!(approx_eq(f.get(), 0.0, 1e-12));
    }

    #[test]
    fn test_new_preserves_zero() {
        let f = FitnessValue::new(0.0);
        assert!(approx_eq(f.get(), 0.0, 1e-12));
    }

    #[test]
    fn test_new_preserves_midrange() {
        let f = FitnessValue::new(0.5);
        assert!(approx_eq(f.get(), 0.5, 1e-12));
    }

    #[test]
    fn test_new_preserves_one() {
        let f = FitnessValue::new(1.0);
        assert!(approx_eq(f.get(), 1.0, 1e-12));
    }

    #[test]
    fn test_new_clamps_above_one_to_one() {
        let f = FitnessValue::new(2.0);
        assert!(approx_eq(f.get(), 1.0, 1e-12));
    }

    #[test]
    fn test_from_f64_clamps_and_preserves_values() {
        let f_from_neg: FitnessValue = (-1.0).into();
        assert!(approx_eq(f_from_neg.get(), 0.0, 1e-12));

        let f_from_pos: FitnessValue = 0.75.into();
        assert!(approx_eq(f_from_pos.get(), 0.75, 1e-12));

        let f_from_big: FitnessValue = 10.0.into();
        assert!(approx_eq(f_from_big.get(), 1.0, 1e-12));
    }

    #[test]
    fn test_new_nan_passes_through() {
        let nan_val = f64::NAN;
        let nan_f = FitnessValue::new(nan_val);
        assert!(nan_f.get().is_nan());
    }

    #[test]
    fn test_ln_and_exp_roundtrip_zero() {
        let input = 0.0;
        let f = FitnessValue::from(input);
        let log = f.ln();
        let back = log.exp();
        assert!(log.get().is_infinite() && log.get().is_sign_negative());
        assert!(back.get() == 0.0);
    }

    #[test]
    fn test_ln_and_exp_roundtrip_one() {
        let input = 1.0;
        let f = FitnessValue::from(input);
        let log = f.ln();
        let back = log.exp();
        assert!(approx_eq(f.get(), back.get(), 1e-12));
    }

    #[test]
    fn test_ln_and_exp_roundtrip_mid() {
        let input = 0.5;
        let f = FitnessValue::from(input);
        let log = f.ln();
        let back = log.exp();
        assert!(approx_eq(f.get(), back.get(), 1e-12));
    }

    #[test]
    fn test_ln_and_exp_roundtrip_very_small() {
        let input = 1e-308;
        let f = FitnessValue::from(input);
        let log = f.ln();
        let back = log.exp();
        assert!(approx_eq(f.get(), back.get(), 1e-12));
    }

    #[test]
    fn test_ln_and_exp_roundtrip_nan() {
        let nan_f = FitnessValue::from(f64::NAN);
        let log_nan = nan_f.ln();
        assert!(log_nan.get().is_nan());
        let back = log_nan.exp();
        assert!(back.get().is_nan());
    }

    #[test]
    fn test_log_exp_of_ln_half() {
        let log_val = LogFitnessValue::new(-std::f64::consts::LN_2); // ln(0.5)
        let lin = log_val.exp();
        assert!(approx_eq(lin.get(), 0.5, 1e-12));
    }

    #[test]
    fn test_log_from_f64_preserves_value() {
        let log_from_f: LogFitnessValue = (-1.234).into();
        assert!(approx_eq(log_from_f.get(), -1.234, 1e-12));
    }

    #[test]
    fn test_conversion_roundtrip_for_quarter() {
        let f = FitnessValue::new(0.25);
        let log = LogFitnessValue::from(f);
        let back: FitnessValue = log.into();
        assert!(approx_eq(f.get(), back.get(), 1e-12));
    }

    #[test]
    fn test_defaults_are_correct() {
        let default_f = FitnessValue::default();
        assert!(approx_eq(default_f.get(), 1.0, 1e-12));

        let default_log = LogFitnessValue::default();
        assert!(approx_eq(default_log.get(), 0.0, 1e-12));
    }

    #[test]
    fn test_display_parsable_for_default() {
        let default_f = FitnessValue::default();
        let disp = default_f.to_string();
        let parsed: f64 = disp.parse().expect("display should be parsable as f64");
        assert!(approx_eq(parsed, 1.0, 1e-12));
    }

    #[test]
    fn test_log_new_clamps_positive_to_zero() {
        let log_val = LogFitnessValue::new(5.0);
        assert!(approx_eq(log_val.get(), 0.0, 1e-12));
    }

    #[test]
    fn test_log_new_preserves_negative() {
        let log_val = LogFitnessValue::new(-2.5);
        assert!(approx_eq(log_val.get(), -2.5, 1e-12));
    }

    #[test]
    fn test_log_new_preserves_neg_infinity() {
        let log_val = LogFitnessValue::new(f64::NEG_INFINITY);
        assert!(log_val.get().is_infinite() && log_val.get().is_sign_negative());
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
        let expected = FitnessValue::new(0.25);
        let result = combined.exp();
        assert!(approx_eq(result.get(), expected.get(), 1e-12));
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
        let f1 = FitnessValue::new(0.5);
        let f2 = FitnessValue::new(0.5);
        let result = f1 * f2;
        assert!(approx_eq(result.get(), 0.25, 1e-12));
    }

    #[test]
    fn test_mul_with_one_preserves_value() {
        let f1 = FitnessValue::new(0.7);
        let f2 = FitnessValue::new(1.0);
        let result = f1 * f2;
        assert!(approx_eq(result.get(), 0.7, 1e-12));
    }

    #[test]
    fn test_mul_with_zero_gives_zero() {
        let f1 = FitnessValue::new(0.8);
        let f2 = FitnessValue::new(0.0);
        let result = f1 * f2;
        assert!(approx_eq(result.get(), 0.0, 1e-12));
    }

    #[test]
    fn test_mul_chain_three_values() {
        let f1 = FitnessValue::new(0.5);
        let f2 = FitnessValue::new(0.5);
        let f3 = FitnessValue::new(0.5);
        let result = f1 * f2 * f3;
        assert!(approx_eq(result.get(), 0.125, 1e-12));
    }

    #[test]
    fn test_mul_very_small_values() {
        let f1 = FitnessValue::new(1e-100);
        let f2 = FitnessValue::new(1e-100);
        let result = f1 * f2;
        // Direct multiplication would underflow, but log-space handles it
        assert!(result.get() > 0.0);
        assert!(approx_eq(result.get(), 1e-200, 1e-12));
    }

    #[test]
    fn test_add_fitness_values() {
        let f1 = FitnessValue::new(0.2);
        let f2 = FitnessValue::new(0.3);
        let result = f1 + f2;
        assert!(approx_eq(result.get(), 0.5, 1e-12));
    }

    #[test]
    fn test_add_fitness_values_clamp() {
        let f1 = FitnessValue::new(0.7);
        let f2 = FitnessValue::new(0.5);
        let result = f1 + f2;
        assert!(approx_eq(result.get(), 1.0, 1e-12));
    }

    #[test]
    fn test_add_log_values() {
        let log1 = LogFitnessValue::new(-std::f64::consts::LN_2); // ln(0.5)
        let log2 = LogFitnessValue::new(-std::f64::consts::LN_2); // ln(0.5)
        let combined = log1 + log2;
        let expected = FitnessValue::new(0.25);
        let result = combined.exp();
        assert!(approx_eq(result.get(), expected.get(), 1e-12));
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
}
