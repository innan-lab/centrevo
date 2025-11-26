//! Recombination operations for sequences.
//!
//! This module provides recombination functionality including crossover
//! and gene conversion events.

use crate::base::Sequence;
pub use crate::errors::RecombinationError;
use rand::Rng;
use rand_distr::Geometric;
use serde::{Deserialize, Serialize};

/// Type of recombination event that can occur on a sequence.
///
/// Variants describe whether no recombination occurred, a crossover at a
/// specific position, or a gene conversion tract (start..end).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum RecombinationType {
    /// No recombination occurs
    None,
    /// Crossover (recombination) at a specific position
    Crossover { position: usize },
    /// Gene conversion from start to end position
    GeneConversion { start: usize, end: usize },
}

/// Parameters controlling recombination behavior.
///
/// These parameters determine the per-base break probability and how breaks
/// are resolved (crossover vs gene conversion) as well as the extension
/// behavior of gene conversion tracts. Use `RecombinationParams::new` to
/// validate values before using sampling/operation helpers.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RecombinationParams {
    /// Probability of a double-strand break per base per generation
    break_prob: f64,
    /// Probability that a break is repaired via crossover (vs gene conversion)
    crossover_prob: f64,
    /// Probability that gene conversion extends by one more base
    gc_extension_prob: f64,
}

impl RecombinationParams {
    /// Create new recombination parameters.
    ///
    /// # Arguments
    /// * `break_prob` - Probability of break per base per generation [0.0, 1.0]
    /// * `crossover_prob` - Probability break is crossover vs gene conversion [0.0, 1.0]
    /// * `gc_extension_prob` - Probability gene conversion extends one more base [0.0, 1.0]
    ///
    /// # Errors
    /// Returns an error if any probability is outside [0.0, 1.0].
    pub fn new(
        break_prob: f64,
        crossover_prob: f64,
        gc_extension_prob: f64,
    ) -> Result<Self, RecombinationError> {
        if !(0.0..=1.0).contains(&break_prob) {
            return Err(RecombinationError::InvalidProbability(
                "break_prob",
                break_prob,
            ));
        }
        if !(0.0..=1.0).contains(&crossover_prob) {
            return Err(RecombinationError::InvalidProbability(
                "crossover_prob",
                crossover_prob,
            ));
        }
        if !(0.0..=1.0).contains(&gc_extension_prob) {
            return Err(RecombinationError::InvalidProbability(
                "gc_extension_prob",
                gc_extension_prob,
            ));
        }

        Ok(Self {
            break_prob,
            crossover_prob,
            gc_extension_prob,
        })
    }

    /// Get the break probability.
    #[inline]
    pub fn break_prob(&self) -> f64 {
        self.break_prob
    }

    /// Get the crossover probability.
    #[inline]
    pub fn crossover_prob(&self) -> f64 {
        self.crossover_prob
    }

    /// Get the gene conversion extension probability.
    #[inline]
    pub fn gc_extension_prob(&self) -> f64 {
        self.gc_extension_prob
    }

    /// Sample recombination events for a sequence of given length.
    ///
    /// Returns a list of events that occurred, sorted by position.
    /// This uses a geometric distribution to correctly simulate the per-base
    /// break probability, ensuring that the probability of at least one break
    /// scales correctly with sequence length (unlike the previous implementation).
    pub fn sample_events<R: Rng + ?Sized>(
        &self,
        length: usize,
        rng: &mut R,
    ) -> Vec<RecombinationType> {
        let mut events = Vec::new();
        if length == 0 || self.break_prob <= 0.0 {
            return events;
        }

        // If break_prob is 1.0, we have a break at every position.
        // This is an edge case, but we should handle it.
        if self.break_prob >= 1.0 {
            for position in 0..length {
                self.add_event_at(position, length, rng, &mut events);
            }
            return events;
        }

        // Use geometric distribution to skip bases between breaks.
        // Geometric(p) gives the number of failures (non-breaks) before the first success (break).
        // So if we are at `current_pos`, the next break is at `current_pos + k` where k ~ Geometric(p).
        // Note: rand_distr::Geometric is defined on k = 0, 1, 2... (failures before success).
        // So if k=0, the break is at current_pos.
        let geo = Geometric::new(self.break_prob).unwrap();

        let mut current_pos = 0;
        loop {
            let skip: u64 = rng.sample(geo);
            current_pos += skip as usize;

            if current_pos >= length {
                break;
            }

            self.add_event_at(current_pos, length, rng, &mut events);

            // Move past the current break for the next iteration
            current_pos += 1;
        }

        events
    }

    // Helper to add an event at a specific position
    fn add_event_at<R: Rng + ?Sized>(
        &self,
        position: usize,
        length: usize,
        rng: &mut R,
        events: &mut Vec<RecombinationType>,
    ) {
        // Determine if crossover or gene conversion
        if rng.random::<f64>() < self.crossover_prob {
            events.push(RecombinationType::Crossover { position });
        } else {
            // Gene conversion: sample tract length
            // Tract length starts at 1 and extends with probability gc_extension_prob
            // This is effectively a geometric distribution on the extension (k failures before success of stopping?)
            // Or just a loop as before.
            let mut tract_length = 1;
            while (position + tract_length) < length && rng.random::<f64>() < self.gc_extension_prob
            {
                tract_length += 1;
            }

            let end = position + tract_length;
            events.push(RecombinationType::GeneConversion {
                start: position,
                end,
            });
        }
    }

    /// Perform crossover between two sequences at the given position.
    ///
    /// Creates two offspring by swapping sequence content after the crossover point.
    ///
    /// # Returns
    /// A tuple of two new sequences (offspring1, offspring2).
    ///
    /// # Errors
    /// Returns an error if the sequences have different lengths.
    pub fn crossover(
        &self,
        seq1: &Sequence,
        seq2: &Sequence,
        position: usize,
    ) -> Result<(Sequence, Sequence), RecombinationError> {
        if seq1.len() != seq2.len() {
            return Err(RecombinationError::LengthMismatch {
                len1: seq1.len(),
                len2: seq2.len(),
            });
        }

        if position >= seq1.len() {
            return Err(RecombinationError::InvalidPosition {
                position,
                length: seq1.len(),
            });
        }

        // Create offspring by combining parts
        let mut offspring1_data = Vec::with_capacity(seq1.len());
        let mut offspring2_data = Vec::with_capacity(seq2.len());

        // Copy first part
        offspring1_data.extend_from_slice(&seq1.as_slice()[..position]);
        offspring2_data.extend_from_slice(&seq2.as_slice()[..position]);

        // Copy second part (swapped)
        offspring1_data.extend_from_slice(&seq2.as_slice()[position..]);
        offspring2_data.extend_from_slice(&seq1.as_slice()[position..]);

        Ok((
            Sequence::from_nucleotides(offspring1_data),
            Sequence::from_nucleotides(offspring2_data),
        ))
    }

    /// Perform gene conversion by copying a tract from donor to recipient.
    ///
    /// Copies the sequence content from `donor` to `recipient` in the range [start, end).
    ///
    /// # Returns
    /// A new sequence with the converted tract.
    ///
    /// # Errors
    /// Returns an error if sequences have different lengths or indices are invalid.
    pub fn gene_conversion(
        &self,
        recipient: &Sequence,
        donor: &Sequence,
        start: usize,
        end: usize,
    ) -> Result<Sequence, RecombinationError> {
        if recipient.len() != donor.len() {
            return Err(RecombinationError::LengthMismatch {
                len1: recipient.len(),
                len2: donor.len(),
            });
        }

        if start >= recipient.len() {
            return Err(RecombinationError::InvalidPosition {
                position: start,
                length: recipient.len(),
            });
        }

        if end > recipient.len() {
            return Err(RecombinationError::InvalidPosition {
                position: end,
                length: recipient.len(),
            });
        }

        if start >= end {
            return Err(RecombinationError::InvalidRange { start, end });
        }

        // Create new sequence with converted tract
        let mut new_data = recipient.as_slice().to_vec();
        new_data[start..end].copy_from_slice(&donor.as_slice()[start..end]);

        Ok(Sequence::from_nucleotides(new_data))
    }
}

// Removed RecombinationError definition, imported from crate::errors

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256PlusPlus;
    use std::str::FromStr;

    #[test]
    fn test_recombination_params_new() {
        let params = RecombinationParams::new(0.01, 0.5, 0.1).unwrap();
        assert_eq!(params.break_prob(), 0.01);
        assert_eq!(params.crossover_prob(), 0.5);
        assert_eq!(params.gc_extension_prob(), 0.1);
    }

    #[test]
    fn test_recombination_params_invalid() {
        assert!(RecombinationParams::new(-0.1, 0.5, 0.1).is_err());
        assert!(RecombinationParams::new(0.01, 1.5, 0.1).is_err());
        assert!(RecombinationParams::new(0.01, 0.5, -0.1).is_err());
    }

    #[test]
    fn test_sample_event_zero_break_prob() {
        let params = RecombinationParams::new(0.0, 0.5, 0.1).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

        for _ in 0..100 {
            let events = params.sample_events(100, &mut rng);
            assert!(events.is_empty());
        }
    }

    #[test]
    fn test_sample_event_empty_sequence() {
        let params = RecombinationParams::new(0.5, 0.5, 0.1).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

        let events = params.sample_events(0, &mut rng);
        assert!(events.is_empty());
    }

    #[test]
    fn test_sample_event_types() {
        let params = RecombinationParams::new(1.0, 1.0, 0.0).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

        // With break_prob=1.0 and crossover_prob=1.0, should always get crossover
        // For length 10, we expect 10 events (one at each position)
        let events = params.sample_events(10, &mut rng);
        assert_eq!(events.len(), 10);
        for event in events {
            match event {
                RecombinationType::Crossover { .. } => {}
                _ => panic!("Expected crossover event"),
            }
        }
    }

    #[test]
    fn test_sample_event_gene_conversion() {
        let params = RecombinationParams::new(1.0, 0.0, 0.0).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

        // With break_prob=1.0, crossover_prob=0.0, should get gene conversion
        let events = params.sample_events(10, &mut rng);
        assert!(!events.is_empty());

        for event in events {
            match event {
                RecombinationType::GeneConversion { start, end } => {
                    assert!(start < end);
                    assert!(end <= 10);
                }
                _ => panic!("Expected gene conversion event"),
            }
        }
    }

    #[test]
    fn test_crossover_basic() {
        let params = RecombinationParams::new(0.01, 0.5, 0.1).unwrap();
        let seq1 = Sequence::from_str("AAAA").unwrap();
        let seq2 = Sequence::from_str("TTTT").unwrap();

        let (offspring1, offspring2) = params.crossover(&seq1, &seq2, 2).unwrap();

        assert_eq!(offspring1.to_string(), "AATT");
        assert_eq!(offspring2.to_string(), "TTAA");
    }

    #[test]
    fn test_crossover_at_start() {
        let params = RecombinationParams::new(0.01, 0.5, 0.1).unwrap();
        let seq1 = Sequence::from_str("AAAA").unwrap();
        let seq2 = Sequence::from_str("TTTT").unwrap();

        let (offspring1, offspring2) = params.crossover(&seq1, &seq2, 0).unwrap();

        assert_eq!(offspring1.to_string(), "TTTT");
        assert_eq!(offspring2.to_string(), "AAAA");
    }

    #[test]
    fn test_crossover_length_mismatch() {
        let params = RecombinationParams::new(0.01, 0.5, 0.1).unwrap();
        let seq1 = Sequence::from_str("AAA").unwrap();
        let seq2 = Sequence::from_str("TTTT").unwrap();

        let result = params.crossover(&seq1, &seq2, 1);
        assert!(result.is_err());
    }

    #[test]
    fn test_crossover_invalid_position() {
        let params = RecombinationParams::new(0.01, 0.5, 0.1).unwrap();
        let seq1 = Sequence::from_str("AAAA").unwrap();
        let seq2 = Sequence::from_str("TTTT").unwrap();

        let result = params.crossover(&seq1, &seq2, 10);
        assert!(result.is_err());
    }

    #[test]
    fn test_gene_conversion_basic() {
        let params = RecombinationParams::new(0.01, 0.5, 0.1).unwrap();
        let recipient = Sequence::from_str("AAAA").unwrap();
        let donor = Sequence::from_str("TTTT").unwrap();

        let result = params.gene_conversion(&recipient, &donor, 1, 3).unwrap();

        assert_eq!(result.to_string(), "ATTA");
    }

    #[test]
    fn test_gene_conversion_full() {
        let params = RecombinationParams::new(0.01, 0.5, 0.1).unwrap();
        let recipient = Sequence::from_str("AAAA").unwrap();
        let donor = Sequence::from_str("TTTT").unwrap();

        let result = params.gene_conversion(&recipient, &donor, 0, 4).unwrap();

        assert_eq!(result.to_string(), "TTTT");
    }

    #[test]
    fn test_gene_conversion_single_base() {
        let params = RecombinationParams::new(0.01, 0.5, 0.1).unwrap();
        let recipient = Sequence::from_str("AAAA").unwrap();
        let donor = Sequence::from_str("TTTT").unwrap();

        let result = params.gene_conversion(&recipient, &donor, 2, 3).unwrap();

        assert_eq!(result.to_string(), "AATA");
    }

    #[test]
    fn test_gene_conversion_length_mismatch() {
        let params = RecombinationParams::new(0.01, 0.5, 0.1).unwrap();
        let recipient = Sequence::from_str("AAA").unwrap();
        let donor = Sequence::from_str("TTTT").unwrap();

        let result = params.gene_conversion(&recipient, &donor, 1, 2);
        assert!(result.is_err());
    }

    #[test]
    fn test_gene_conversion_invalid_range() {
        let params = RecombinationParams::new(0.01, 0.5, 0.1).unwrap();
        let recipient = Sequence::from_str("AAAA").unwrap();
        let donor = Sequence::from_str("TTTT").unwrap();

        // start >= end
        let result = params.gene_conversion(&recipient, &donor, 2, 2);
        assert!(result.is_err());

        let result = params.gene_conversion(&recipient, &donor, 3, 2);
        assert!(result.is_err());
    }

    #[test]
    fn test_gene_conversion_invalid_position() {
        let params = RecombinationParams::new(0.01, 0.5, 0.1).unwrap();
        let recipient = Sequence::from_str("AAAA").unwrap();
        let donor = Sequence::from_str("TTTT").unwrap();

        let result = params.gene_conversion(&recipient, &donor, 10, 11);
        assert!(result.is_err());

        let result = params.gene_conversion(&recipient, &donor, 0, 10);
        assert!(result.is_err());
    }

    #[test]
    fn test_recombination_type_equality() {
        assert_eq!(RecombinationType::None, RecombinationType::None);
        assert_eq!(
            RecombinationType::Crossover { position: 10 },
            RecombinationType::Crossover { position: 10 }
        );
        assert_ne!(
            RecombinationType::Crossover { position: 10 },
            RecombinationType::Crossover { position: 20 }
        );
    }

    #[test]
    fn test_error_display() {
        let err = RecombinationError::InvalidProbability("test", 1.5);
        let msg = format!("{err}");
        assert!(msg.contains("Invalid probability"));

        let err = RecombinationError::LengthMismatch { len1: 10, len2: 20 };
        let msg = format!("{err}");
        assert!(msg.contains("length mismatch"));
    }

    #[test]
    fn test_params_clone() {
        let params1 = RecombinationParams::new(0.01, 0.5, 0.1).unwrap();
        let params2 = params1.clone();

        assert_eq!(params1.break_prob(), params2.break_prob());
        assert_eq!(params1.crossover_prob(), params2.crossover_prob());
        assert_eq!(params1.gc_extension_prob(), params2.gc_extension_prob());
    }

    #[test]
    fn test_sample_events_multiple() {
        // High break prob, long sequence -> multiple events
        let params = RecombinationParams::new(0.1, 0.5, 0.1).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

        let events = params.sample_events(1000, &mut rng);
        // Expected events ~ 100
        assert!(events.len() > 50);
        assert!(events.len() < 150);

        // Verify sorted positions
        let mut last_pos = 0;
        for event in events {
            let pos = match event {
                RecombinationType::Crossover { position } => position,
                RecombinationType::GeneConversion { start, .. } => start,
                RecombinationType::None => panic!("Should not have None in events list"),
            };
            assert!(pos >= last_pos);
            last_pos = pos;
        }
    }

    #[test]
    fn test_sample_event_probability_scaling() {
        // break_prob is 0.01 per base.
        // For length 100, expected breaks = 1.
        // Probability of at least one break = 1 - (0.99)^100 ≈ 0.63

        let params = RecombinationParams::new(0.01, 0.5, 0.1).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

        let mut event_count = 0;
        let trials = 10000;
        for _ in 0..trials {
            let events = params.sample_events(100, &mut rng);
            if !events.is_empty() {
                event_count += 1;
            }
        }

        let frequency = event_count as f64 / trials as f64;
        println!("Frequency of events for L=100, p=0.01: {}", frequency);

        // Now we expect correct scaling
        assert!(
            frequency > 0.60 && frequency < 0.66,
            "Frequency {} should be around 0.63",
            frequency
        );
    }
}
