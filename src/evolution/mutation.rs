//! Mutation operations for sequences.
//!
//! This module provides mutation functionality including point substitutions
//! based on substitution models (e.g., JC69, uniform).

use rand::Rng;
use crate::base::{Alphabet, Nucleotide, Sequence};

/// Substitution model for nucleotide mutations.
///
/// This implements a simple uniform substitution model where all off-diagonal
/// substitutions occur with equal probability.
#[derive(Debug, Clone)]
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
    pub fn mutate_base<R: Rng + ?Sized>(&self, base: Nucleotide, rng: &mut R) -> Nucleotide {
        // Check if mutation occurs
        if !rng.random_bool(self.mu) {
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
    ///
    /// # Returns
    /// The number of mutations that occurred.
    pub fn mutate_sequence<R: Rng + ?Sized>(
        &self,
        sequence: &mut Sequence,
        rng: &mut R,
    ) -> usize {
        let mut mutation_count = 0;

        for i in 0..sequence.len() {
            if let Some(base) = sequence.get(i) {
                let new_base = self.mutate_base(base, rng);
                if new_base != base {
                    sequence.set(i, new_base).expect("valid index");
                    mutation_count += 1;
                }
            }
        }

        mutation_count
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
    use rand::rngs::StdRng;

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
        let mut rng = StdRng::seed_from_u64(42);
        
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
        let mut rng = StdRng::seed_from_u64(42);
        
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
        let mut rng = StdRng::seed_from_u64(42);
        
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
        let mut rng = StdRng::seed_from_u64(42);
        
        let mut seq = Sequence::from_str("ACGTACGT", test_alphabet()).unwrap();
        let original = seq.to_string();
        
        let count = model.mutate_sequence(&mut seq, &mut rng);
        
        assert_eq!(count, 0);
        assert_eq!(seq.to_string(), original);
    }

    #[test]
    fn test_mutate_sequence_low_rate() {
        let model = SubstitutionModel::new(test_alphabet(), 0.1).unwrap();
        let mut rng = StdRng::seed_from_u64(42);
        
        let mut seq = Sequence::from_str("ACGTACGTACGTACGT", test_alphabet()).unwrap();
        let count = model.mutate_sequence(&mut seq, &mut rng);
        
        // With low rate, should have few mutations
        assert!(count < 5);
    }

    #[test]
    fn test_mutate_sequence_high_rate() {
        let model = SubstitutionModel::new(test_alphabet(), 0.9).unwrap();
        let mut rng = StdRng::seed_from_u64(42);
        
        let mut seq = Sequence::from_str("ACGTACGTACGTACGT", test_alphabet()).unwrap();
        let count = model.mutate_sequence(&mut seq, &mut rng);
        
        // With high rate, should have many mutations
        assert!(count > 10);
    }

    #[test]
    fn test_mutate_sequence_empty() {
        let model = SubstitutionModel::new(test_alphabet(), 0.5).unwrap();
        let mut rng = StdRng::seed_from_u64(42);
        
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
        
        let mut rng1 = StdRng::seed_from_u64(123);
        let mut rng2 = StdRng::seed_from_u64(123);
        
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
}
