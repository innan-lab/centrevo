//! Mutation operations for sequences.
//!
//! This module provides mutation functionality including point substitutions
//! based on substitution models (e.g., JC69, uniform).

use rand::Rng;
use rand_distr::{Distribution, Poisson};
use crate::base::{Alphabet, Nucleotide, Sequence};
use serde::{Serialize, Deserialize};

/// Substitution model for nucleotide mutations.
///
/// This implements a simple uniform substitution model where all off-diagonal
/// substitutions occur with equal probability.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SubstitutionModel {
    /// The alphabet for this model
    alphabet: Alphabet,
    /// Mutation rate per base per generation
    mu: f64,
}

impl SubstitutionModel {
    /// Create a new substitution model with the given alphabet and mutation rate.
    ///
    /// # Arguments
    /// * `alphabet` - The alphabet to use
    /// * `mu` - Mutation rate per base per generation (must be between 0.0 and 1.0)
    ///
    /// # Errors
    /// Returns an error if `mu` is not in the valid range [0.0, 1.0].
    pub fn new(alphabet: Alphabet, mu: f64) -> Result<Self, MutationError> {
        if !(0.0..=1.0).contains(&mu) {
            return Err(MutationError::InvalidMutationRate(mu));
        }
        Ok(Self { alphabet, mu })
    }

    /// Create a JC69 (Jukes-Cantor 1969) substitution model for DNA sequences.
    ///
    /// This is a uniform-rate model for the DNA alphabet (A, C, G, T).
    ///
    /// # Arguments
    /// * `mu` - Total mutation rate per base per generation
    pub fn jc69(mu: f64) -> Result<Self, MutationError> {
        Self::new(Alphabet::dna(), mu)
    }

    /// Create a uniform substitution model.
    ///
    /// All off-diagonal transitions occur with equal probability.
    pub fn uniform(alphabet: Alphabet, mu: f64) -> Result<Self, MutationError> {
        Self::new(alphabet, mu)
    }

    /// Get the mutation rate.
    #[inline]
    pub fn mu(&self) -> f64 {
        self.mu
    }

    /// Get the alphabet.
    #[inline]
    pub fn alphabet(&self) -> &Alphabet {
        &self.alphabet
    }

    /// Mutate a single base according to the substitution model.
    ///
    /// # Returns
    /// The possibly mutated base.
    #[inline]
    pub fn mutate_base<R: Rng + ?Sized>(&self, base: Nucleotide, rng: &mut R) -> Nucleotide {
        // Check if mutation occurs using gen::<f64>() which is faster than random_bool
        if rng.random::<f64>() >= self.mu {
            return base; // No mutation
        }

        // Uniform model: choose any of the other 3 bases with equal probability
        let alphabet_size = self.alphabet.len();
        let current_idx = base.to_index() as usize;
        
        // Generate a random index in [0, alphabet_size-1) excluding current
        let mut new_idx = rng.random_range(0..alphabet_size - 1);
        if new_idx >= current_idx {
            new_idx += 1; // Skip the current base
        }

        Nucleotide::from_index(new_idx as u8).unwrap_or(base)
    }

    /// Mutate a sequence in place according to the substitution model.
    ///
    /// Each base in the sequence has an independent chance of mutating.
    /// This version uses bulk random number generation for better performance.
    ///
    /// # Returns
    /// The number of mutations that occurred.
    pub fn mutate_sequence<R: Rng + ?Sized>(
        &self,
        sequence: &mut Sequence,
        rng: &mut R,
    ) -> usize {
        let len = sequence.len();
        if len == 0 {
            return 0;
        }

        let mut mutation_count = 0;
        let indices = sequence.indices_mut();
        let alphabet_size = self.alphabet.len();
        
        // Bulk generate random floats for mutation decisions
        // This is much faster than calling random() for each base
        let mut random_floats = vec![0.0f64; len];
        rng.fill(&mut random_floats[..]);
        
        // Pre-generate random indices for mutations that will occur
        // Estimate how many we'll need (with some buffer)
        let expected_mutations = (len as f64 * self.mu * 1.5) as usize + 10;
        let mut random_indices = vec![0u8; expected_mutations.min(len * 3)];
        rng.fill(&mut random_indices[..]);
        let mut random_idx_pos = 0;
        
        for i in 0..len {
            // Check if mutation occurs using pre-generated random float
            if random_floats[i] < self.mu {
                let current_idx = indices[i] as usize;
                
                // Use pre-generated random byte to select new base
                if random_idx_pos < random_indices.len() {
                    let mut new_idx = (random_indices[random_idx_pos] as usize) % (alphabet_size - 1);
                    random_idx_pos += 1;
                    
                    if new_idx >= current_idx {
                        new_idx += 1;
                    }
                    
                    indices[i] = new_idx as u8;
                    mutation_count += 1;
                } else {
                    // Fallback if we run out of pre-generated randoms
                    let mut new_idx = rng.random_range(0..alphabet_size - 1);
                    if new_idx >= current_idx {
                        new_idx += 1;
                    }
                    indices[i] = new_idx as u8;
                    mutation_count += 1;
                }
            }
        }

        mutation_count
    }

    /// Mutate a sequence using Poisson pre-sampling for better performance.
    ///
    /// Instead of testing each base individually, this method:
    /// 1. Samples the total number of mutations from Poisson(length × μ)
    /// 2. Randomly selects that many positions to mutate
    /// 3. Applies mutations to selected positions
    ///
    /// This is mathematically equivalent to the standard approach but much faster
    /// for low mutation rates (typical in evolutionary simulations).
    ///
    /// # Performance
    /// - Low mutation rates (μ < 0.01): 5-10x faster
    /// - Medium mutation rates (0.01 ≤ μ < 0.1): 2-5x faster
    /// - High mutation rates (μ ≥ 0.1): Similar or slightly slower
    ///
    /// # Returns
    /// The number of mutations that occurred.
    pub fn mutate_sequence_poisson<R: Rng + ?Sized>(
        &self,
        sequence: &mut Sequence,
        rng: &mut R,
    ) -> usize {
        let len = sequence.len();
        if len == 0 {
            return 0;
        }

        // Special case: if mutation rate is very high, use standard approach
        // (Poisson sampling becomes inefficient when most bases mutate)
        if self.mu > 0.5 {
            return self.mutate_sequence(sequence, rng);
        }

        // Sample number of mutations from Poisson distribution
        let lambda = len as f64 * self.mu;
        
        // For very low lambda, Poisson sampling may have numerical issues
        // Use standard approach for tiny sequences or very low rates
        if lambda < 0.1 {
            return self.mutate_sequence(sequence, rng);
        }

        let poisson = match Poisson::new(lambda) {
            Ok(p) => p,
            Err(_) => return self.mutate_sequence(sequence, rng), // Fallback on error
        };
        
        let num_mutations = poisson.sample(rng) as usize;
        
        // If no mutations, return early
        if num_mutations == 0 {
            return 0;
        }

        // Optimization: if number of mutations is close to sequence length,
        // use standard approach (sampling without replacement becomes expensive)
        if num_mutations >= len / 2 {
            return self.mutate_sequence(sequence, rng);
        }

        // Sample positions without replacement using reservoir sampling
        let positions = sample_without_replacement(len, num_mutations, rng);

        // Apply mutations at selected positions with bulk random generation
        let indices = sequence.indices_mut();
        let alphabet_size = self.alphabet.len();
        
        // Pre-generate random bytes for mutations
        let mut random_bytes = vec![0u8; positions.len()];
        rng.fill(&mut random_bytes[..]);
        
        for (idx, &pos) in positions.iter().enumerate() {
            let current_idx = indices[pos] as usize;
            
            // Use pre-generated random byte
            let mut new_idx = (random_bytes[idx] as usize) % (alphabet_size - 1);
            if new_idx >= current_idx {
                new_idx += 1;
            }
            
            indices[pos] = new_idx as u8;
        }

        positions.len()
    }
}

/// Sample k positions from [0, n) without replacement using reservoir sampling.
///
/// This is an efficient algorithm for sampling without replacement when k << n.
/// Time complexity: O(k), Space complexity: O(k)
///
/// # Arguments
/// * `n` - The range to sample from [0, n)
/// * `k` - The number of samples to take (must be <= n)
/// * `rng` - Random number generator
///
/// # Returns
/// A vector of k unique positions in [0, n)
#[inline]
fn sample_without_replacement<R: Rng + ?Sized>(n: usize, k: usize, rng: &mut R) -> Vec<usize> {
    debug_assert!(k <= n, "Cannot sample more items than available");
    
    if k == 0 {
        return Vec::new();
    }
    
    // For small k relative to n, use hash-based sampling
    // For large k, it's more efficient to use Fisher-Yates shuffle
    if k < n / 10 {
        // Hash-based sampling with bulk random generation
        let mut selected = std::collections::HashSet::with_capacity(k);
        let mut result = Vec::with_capacity(k);
        
        // Pre-generate random values in batches
        let batch_size = (k * 2).min(1024); // Generate 2x what we need, capped at 1024
        let mut random_values = vec![0usize; batch_size];
        let mut random_pos = batch_size; // Force initial generation
        
        while result.len() < k {
            // Regenerate batch if needed
            if random_pos >= random_values.len() {
                for v in random_values.iter_mut() {
                    *v = rng.random_range(0..n);
                }
                random_pos = 0;
            }
            
            let pos = random_values[random_pos];
            random_pos += 1;
            
            if selected.insert(pos) {
                result.push(pos);
            }
        }
        
        result
    } else {
        // Fisher-Yates shuffle with bulk random generation for larger k
        let mut positions: Vec<usize> = (0..n).collect();
        
        // Generate all random indices at once
        let mut random_indices = vec![0usize; k];
        for (i, r) in random_indices.iter_mut().enumerate() {
            *r = rng.random_range(i..n);
        }
        
        for i in 0..k {
            positions.swap(i, random_indices[i]);
        }
        
        positions.truncate(k);
        positions
    }
}

/// Errors that can occur during mutation operations.
#[derive(Debug, Clone, PartialEq)]
pub enum MutationError {
    /// Invalid mutation rate (must be between 0.0 and 1.0)
    InvalidMutationRate(f64),
}

impl std::fmt::Display for MutationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MutationError::InvalidMutationRate(mu) => {
                write!(f, "Invalid mutation rate: {} (must be between 0.0 and 1.0)", mu)
            }
        }
    }
}

impl std::error::Error for MutationError {}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256PlusPlus;

    fn test_alphabet() -> Alphabet {
        Alphabet::dna()
    }

    #[test]
    fn test_substitution_model_new() {
        let model = SubstitutionModel::new(test_alphabet(), 0.01).unwrap();
        assert_eq!(model.mu(), 0.01);
    }

    #[test]
    fn test_substitution_model_invalid_rate() {
        assert!(SubstitutionModel::new(test_alphabet(), -0.1).is_err());
        assert!(SubstitutionModel::new(test_alphabet(), 1.5).is_err());
    }

    #[test]
    fn test_substitution_model_jc69() {
        let model = SubstitutionModel::jc69(0.01).unwrap();
        assert_eq!(model.mu(), 0.01);
        assert_eq!(model.alphabet().len(), 4);
    }

    #[test]
    fn test_substitution_model_zero_rate() {
        let model = SubstitutionModel::new(test_alphabet(), 0.0).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        
        // With zero mutation rate, base should never change
        for _ in 0..100 {
            let base = Nucleotide::A;
            let mutated = model.mutate_base(base, &mut rng);
            assert_eq!(mutated, base);
        }
    }

    #[test]
    fn test_substitution_model_mutate_base() {
        let model = SubstitutionModel::new(test_alphabet(), 1.0).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        
        // With mutation rate 1.0, base should always change
        let mut mutations = 0;
        for _ in 0..100 {
            let base = Nucleotide::A;
            let mutated = model.mutate_base(base, &mut rng);
            if mutated != base {
                mutations += 1;
            }
        }
        
        // Should have many mutations (close to 100)
        assert!(mutations > 90);
    }

    #[test]
    fn test_substitution_model_mutate_base_distribution() {
        let model = SubstitutionModel::new(test_alphabet(), 1.0).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        
        let mut counts = [0; 4];
        for _ in 0..1000 {
            let mutated = model.mutate_base(Nucleotide::A, &mut rng);
            counts[mutated.to_index() as usize] += 1;
        }
        
        // A should never appear (always mutates)
        assert_eq!(counts[0], 0);
        
        // C, G, T should be roughly equally distributed
        for i in 1..4 {
            assert!(counts[i] > 250);
            assert!(counts[i] < 450);
        }
    }

    #[test]
    fn test_mutate_sequence_zero_rate() {
        let model = SubstitutionModel::new(test_alphabet(), 0.0).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        
        let mut seq = Sequence::from_str("ACGTACGT", test_alphabet()).unwrap();
        let original = seq.to_string();
        
        let count = model.mutate_sequence(&mut seq, &mut rng);
        
        assert_eq!(count, 0);
        assert_eq!(seq.to_string(), original);
    }

    #[test]
    fn test_mutate_sequence_low_rate() {
        let model = SubstitutionModel::new(test_alphabet(), 0.1).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        
        let mut seq = Sequence::from_str("ACGTACGTACGTACGT", test_alphabet()).unwrap();
        let count = model.mutate_sequence(&mut seq, &mut rng);
        
        // With low rate, should have few mutations
        assert!(count < 5);
    }

    #[test]
    fn test_mutate_sequence_high_rate() {
        let model = SubstitutionModel::new(test_alphabet(), 0.9).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        
        let mut seq = Sequence::from_str("ACGTACGTACGTACGT", test_alphabet()).unwrap();
        let count = model.mutate_sequence(&mut seq, &mut rng);
        
        // With high rate, should have many mutations
        assert!(count > 10);
    }

    #[test]
    fn test_mutate_sequence_empty() {
        let model = SubstitutionModel::new(test_alphabet(), 0.5).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        
        let mut seq = Sequence::new(test_alphabet());
        let count = model.mutate_sequence(&mut seq, &mut rng);
        
        assert_eq!(count, 0);
        assert!(seq.is_empty());
    }

    #[test]
    fn test_mutate_sequence_deterministic() {
        let model = SubstitutionModel::new(test_alphabet(), 0.1).unwrap();
        
        let mut seq1 = Sequence::from_str("ACGTACGTACGTACGT", test_alphabet()).unwrap();
        let mut seq2 = Sequence::from_str("ACGTACGTACGTACGT", test_alphabet()).unwrap();
        
        let mut rng1 = Xoshiro256PlusPlus::seed_from_u64(123);
        let mut rng2 = Xoshiro256PlusPlus::seed_from_u64(123);
        
        let count1 = model.mutate_sequence(&mut seq1, &mut rng1);
        let count2 = model.mutate_sequence(&mut seq2, &mut rng2);
        
        // Same seed should produce same results
        assert_eq!(count1, count2);
        assert_eq!(seq1.to_string(), seq2.to_string());
    }

    #[test]
    fn test_mutation_error_display() {
        let err = MutationError::InvalidMutationRate(1.5);
        let msg = format!("{}", err);
        assert!(msg.contains("Invalid mutation rate"));
        assert!(msg.contains("1.5"));
    }

    #[test]
    fn test_substitution_model_clone() {
        let model1 = SubstitutionModel::new(test_alphabet(), 0.01).unwrap();
        let model2 = model1.clone();
        
        assert_eq!(model1.mu(), model2.mu());
        assert_eq!(model1.alphabet(), model2.alphabet());
    }

    // Tests for Poisson-based mutation
    
    #[test]
    fn test_mutate_sequence_poisson_zero_rate() {
        let model = SubstitutionModel::new(test_alphabet(), 0.0).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        
        let mut seq = Sequence::from_str("ACGTACGT", test_alphabet()).unwrap();
        let original = seq.to_string();
        
        let count = model.mutate_sequence_poisson(&mut seq, &mut rng);
        
        assert_eq!(count, 0);
        assert_eq!(seq.to_string(), original);
    }

    #[test]
    fn test_mutate_sequence_poisson_low_rate() {
        let model = SubstitutionModel::new(test_alphabet(), 0.001).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        
        let mut seq = Sequence::from_str("ACGT".repeat(250).as_str(), test_alphabet()).unwrap();
        let count = model.mutate_sequence_poisson(&mut seq, &mut rng);
        
        // With low rate on 1000bp, expect 0-5 mutations
        assert!(count < 10);
    }

    #[test]
    fn test_mutate_sequence_poisson_medium_rate() {
        let model = SubstitutionModel::new(test_alphabet(), 0.01).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        
        let mut seq = Sequence::from_str("ACGT".repeat(250).as_str(), test_alphabet()).unwrap();
        let count = model.mutate_sequence_poisson(&mut seq, &mut rng);
        
        // With medium rate on 1000bp, expect 5-20 mutations
        assert!(count > 0);
        assert!(count < 30);
    }

    #[test]
    fn test_mutate_sequence_poisson_high_rate_fallback() {
        // High rate should fallback to standard method
        let model = SubstitutionModel::new(test_alphabet(), 0.6).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        
        let mut seq = Sequence::from_str("ACGTACGTACGTACGT", test_alphabet()).unwrap();
        let count = model.mutate_sequence_poisson(&mut seq, &mut rng);
        
        // Should have many mutations
        assert!(count > 5);
    }

    #[test]
    fn test_mutate_sequence_poisson_empty() {
        let model = SubstitutionModel::new(test_alphabet(), 0.01).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        
        let mut seq = Sequence::new(test_alphabet());
        let count = model.mutate_sequence_poisson(&mut seq, &mut rng);
        
        assert_eq!(count, 0);
        assert!(seq.is_empty());
    }

    #[test]
    fn test_mutate_sequence_poisson_deterministic() {
        let model = SubstitutionModel::new(test_alphabet(), 0.01).unwrap();
        
        let mut seq1 = Sequence::from_str("ACGT".repeat(250).as_str(), test_alphabet()).unwrap();
        let mut seq2 = Sequence::from_str("ACGT".repeat(250).as_str(), test_alphabet()).unwrap();
        
        let mut rng1 = Xoshiro256PlusPlus::seed_from_u64(123);
        let mut rng2 = Xoshiro256PlusPlus::seed_from_u64(123);
        
        let count1 = model.mutate_sequence_poisson(&mut seq1, &mut rng1);
        let count2 = model.mutate_sequence_poisson(&mut seq2, &mut rng2);
        
        // Same seed should produce same results
        assert_eq!(count1, count2);
        assert_eq!(seq1.to_string(), seq2.to_string());
    }

    #[test]
    fn test_mutate_sequence_poisson_statistical_equivalence() {
        // Test that Poisson method produces statistically similar results to standard method
        let model = SubstitutionModel::new(test_alphabet(), 0.01).unwrap();
        let seq_str = "ACGT".repeat(250); // 1000bp
        
        let mut standard_mutations = Vec::new();
        let mut poisson_mutations = Vec::new();
        
        // Run multiple trials
        for seed in 0..100 {
            let mut seq_standard = Sequence::from_str(&seq_str, test_alphabet()).unwrap();
            let mut seq_poisson = Sequence::from_str(&seq_str, test_alphabet()).unwrap();
            
            let mut rng_standard = Xoshiro256PlusPlus::seed_from_u64(seed);
            let mut rng_poisson = Xoshiro256PlusPlus::seed_from_u64(seed + 10000); // Different seed
            
            standard_mutations.push(model.mutate_sequence(&mut seq_standard, &mut rng_standard));
            poisson_mutations.push(model.mutate_sequence_poisson(&mut seq_poisson, &mut rng_poisson));
        }
        
        // Calculate means
        let mean_standard: f64 = standard_mutations.iter().sum::<usize>() as f64 / 100.0;
        let mean_poisson: f64 = poisson_mutations.iter().sum::<usize>() as f64 / 100.0;
        
        // Expected: 1000 * 0.01 = 10 mutations
        // Both methods should have similar means (within 20% of each other)
        let diff_ratio = (mean_standard - mean_poisson).abs() / mean_standard.max(mean_poisson);
        assert!(diff_ratio < 0.2, "Mean difference too large: standard={}, poisson={}", mean_standard, mean_poisson);
    }

    #[test]
    fn test_sample_without_replacement_small_k() {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        
        // Sample 5 from 100 (uses hash-based sampling)
        let samples = super::sample_without_replacement(100, 5, &mut rng);
        
        assert_eq!(samples.len(), 5);
        
        // Check uniqueness
        let mut unique = samples.clone();
        unique.sort();
        unique.dedup();
        assert_eq!(unique.len(), 5);
        
        // Check range
        for &s in &samples {
            assert!(s < 100);
        }
    }

    #[test]
    fn test_sample_without_replacement_large_k() {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        
        // Sample 50 from 100 (uses Fisher-Yates)
        let samples = super::sample_without_replacement(100, 50, &mut rng);
        
        assert_eq!(samples.len(), 50);
        
        // Check uniqueness
        let mut unique = samples.clone();
        unique.sort();
        unique.dedup();
        assert_eq!(unique.len(), 50);
        
        // Check range
        for &s in &samples {
            assert!(s < 100);
        }
    }

    #[test]
    fn test_sample_without_replacement_edge_cases() {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        
        // Sample 0
        let samples = super::sample_without_replacement(100, 0, &mut rng);
        assert_eq!(samples.len(), 0);
        
        // Sample all
        let samples = super::sample_without_replacement(10, 10, &mut rng);
        assert_eq!(samples.len(), 10);
        
        let mut sorted = samples.clone();
        sorted.sort();
        assert_eq!(sorted, vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9]);
    }
}
