//! Utility functions for analysis module
//!
//! Shared helper functions used across analysis submodules.

use centrevo_sim::base::{Nucleotide, Sequence};

/// Calculate mean of a vector
pub fn mean(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    values.iter().sum::<f64>() / values.len() as f64
}

/// Optimized Hamming distance calculation between two sequences
///
/// This function is performance-critical and used throughout the analysis module.
/// It directly accesses sequence indices for better performance.
#[inline]
pub fn hamming_distance_fast(seq1: &Sequence, seq2: &Sequence) -> usize {
    let len = seq1.len().min(seq2.len());
    let indices1 = seq1.as_slice();
    let indices2 = seq2.as_slice();

    // Process in chunks of 8 for better CPU pipelining
    let mut distance = 0;
    let chunks = len / 8;
    let remainder = len % 8;

    // Process 8 elements at a time
    for i in 0..chunks {
        let base = i * 8;
        distance += (indices1[base] != indices2[base]) as usize;
        distance += (indices1[base + 1] != indices2[base + 1]) as usize;
        distance += (indices1[base + 2] != indices2[base + 2]) as usize;
        distance += (indices1[base + 3] != indices2[base + 3]) as usize;
        distance += (indices1[base + 4] != indices2[base + 4]) as usize;
        distance += (indices1[base + 5] != indices2[base + 5]) as usize;
        distance += (indices1[base + 6] != indices2[base + 6]) as usize;
        distance += (indices1[base + 7] != indices2[base + 7]) as usize;
    }

    // Process remaining elements
    let base = chunks * 8;
    for i in 0..remainder {
        distance += (indices1[base + i] != indices2[base + i]) as usize;
    }

    distance
}

/// Cache for commonly used harmonic numbers
/// Avoids recalculating for common population sizes
/// Note: harmonic_number(n) = sum_{i=1}^{n-1} 1/i
static HARMONIC_CACHE: [f64; 11] = [
    0.0,                // n=0 (unused)
    0.0,                // n=1
    1.0,                // n=2
    1.5,                // n=3
    1.8333333333333333, // n=4
    2.083333333333333,  // n=5
    2.283333333333333,  // n=6
    2.45,               // n=7 (rounded to avoid clippy warning)
    2.5928571428571425, // n=8
    2.7178571428571425, // n=9
    2.8289682539682537, // n=10
];

/// Calculate harmonic number efficiently
#[inline]
pub fn harmonic_number(n: usize) -> f64 {
    if n < HARMONIC_CACHE.len() {
        HARMONIC_CACHE[n]
    } else {
        (1..n).map(|i| 1.0 / i as f64).sum()
    }
}

/// Optimized nucleotide composition counting
/// Returns (A_count, C_count, G_count, T_count)
#[inline]
pub fn count_nucleotides_fast(seq: &Sequence) -> (usize, usize, usize, usize) {
    let indices = seq.as_slice();
    let mut a_count = 0;
    let mut c_count = 0;
    let mut g_count = 0;
    let mut t_count = 0;

    // Process in chunks of 8 for better CPU pipelining
    let len = indices.len();
    let chunks = len / 8;
    let remainder = len % 8;

    for i in 0..chunks {
        let base = i * 8;
        a_count += (indices[base] == Nucleotide::A) as usize;
        c_count += (indices[base] == Nucleotide::C) as usize;
        g_count += (indices[base] == Nucleotide::G) as usize;
        t_count += (indices[base] == Nucleotide::T) as usize;

        a_count += (indices[base + 1] == Nucleotide::A) as usize;
        c_count += (indices[base + 1] == Nucleotide::C) as usize;
        g_count += (indices[base + 1] == Nucleotide::G) as usize;
        t_count += (indices[base + 1] == Nucleotide::T) as usize;

        a_count += (indices[base + 2] == Nucleotide::A) as usize;
        c_count += (indices[base + 2] == Nucleotide::C) as usize;
        g_count += (indices[base + 2] == Nucleotide::G) as usize;
        t_count += (indices[base + 2] == Nucleotide::T) as usize;

        a_count += (indices[base + 3] == Nucleotide::A) as usize;
        c_count += (indices[base + 3] == Nucleotide::C) as usize;
        g_count += (indices[base + 3] == Nucleotide::G) as usize;
        t_count += (indices[base + 3] == Nucleotide::T) as usize;

        a_count += (indices[base + 4] == Nucleotide::A) as usize;
        c_count += (indices[base + 4] == Nucleotide::C) as usize;
        g_count += (indices[base + 4] == Nucleotide::G) as usize;
        t_count += (indices[base + 4] == Nucleotide::T) as usize;

        a_count += (indices[base + 5] == Nucleotide::A) as usize;
        c_count += (indices[base + 5] == Nucleotide::C) as usize;
        g_count += (indices[base + 5] == Nucleotide::G) as usize;
        t_count += (indices[base + 5] == Nucleotide::T) as usize;

        a_count += (indices[base + 6] == Nucleotide::A) as usize;
        c_count += (indices[base + 6] == Nucleotide::C) as usize;
        g_count += (indices[base + 6] == Nucleotide::G) as usize;
        t_count += (indices[base + 6] == Nucleotide::T) as usize;

        a_count += (indices[base + 7] == Nucleotide::A) as usize;
        c_count += (indices[base + 7] == Nucleotide::C) as usize;
        g_count += (indices[base + 7] == Nucleotide::G) as usize;
        t_count += (indices[base + 7] == Nucleotide::T) as usize;
    }

    // Process remaining elements
    let base = chunks * 8;
    for i in 0..remainder {
        let idx = indices[base + i];
        a_count += (idx == Nucleotide::A) as usize;
        c_count += (idx == Nucleotide::C) as usize;
        g_count += (idx == Nucleotide::G) as usize;
        t_count += (idx == Nucleotide::T) as usize;
    }

    (a_count, c_count, g_count, t_count)
}

/// Fast GC content calculation
#[inline]
pub fn gc_content_fast(seq: &Sequence) -> f64 {
    let (_, c_count, g_count, _) = count_nucleotides_fast(seq);
    let len = seq.len();
    if len == 0 {
        0.0
    } else {
        (c_count + g_count) as f64 / len as f64
    }
}

/// Calculate standard deviation
pub fn std_dev(values: &[f64]) -> f64 {
    if values.len() < 2 {
        return 0.0;
    }

    let mean_val = mean(values);
    let variance =
        values.iter().map(|v| (v - mean_val).powi(2)).sum::<f64>() / (values.len() - 1) as f64;

    variance.sqrt()
}

/// Calculate median of a vector
pub fn median(values: &mut [f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }

    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mid = values.len() / 2;

    if values.len().is_multiple_of(2) {
        (values[mid - 1] + values[mid]) / 2.0
    } else {
        values[mid]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mean() {
        assert_eq!(mean(&[1.0, 2.0, 3.0, 4.0, 5.0]), 3.0);
        assert_eq!(mean(&[]), 0.0);
        assert_eq!(mean(&[5.0]), 5.0);
    }

    #[test]
    fn test_std_dev() {
        let values = vec![2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        let sd = std_dev(&values);
        assert!((sd - 2.138).abs() < 0.01);
    }

    #[test]
    fn test_median() {
        let mut values = vec![1.0, 3.0, 2.0, 5.0, 4.0];
        assert_eq!(median(&mut values), 3.0);

        let mut values_even = vec![1.0, 2.0, 3.0, 4.0];
        assert_eq!(median(&mut values_even), 2.5);
    }
}
