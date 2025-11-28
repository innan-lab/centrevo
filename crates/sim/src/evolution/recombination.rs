//! Recombination operations for sequences.
//!
//! This module provides recombination functionality including crossover
//! and gene conversion events.

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
    /// Crossover (recombination) at specific positions in both chromosomes
    Crossover { pos1: usize, pos2: usize },
    /// Gene conversion from start to end position
    // `donor_start` is the position in the first chromosome passed to `sample_events`
    // (source/donor). `recipient_start` is the position in the second chromosome passed
    // to `sample_events` (target/recipient).
    GeneConversion {
        donor_start: usize,
        recipient_start: usize,
        length: usize,
    },
}

/// Parameters controlling recombination behavior.
///
/// These parameters determine the per-base break probability and how breaks
/// are resolved (crossover vs gene conversion) as well as the extension
/// behavior of gene conversion tracts. Use `RecombinationModel::builder()` to
/// construct a new instance with validation.
///
/// # Examples
///
/// ```
/// use centrevo_sim::evolution::RecombinationModel;
///
/// // Create a model with custom parameters
/// let model = RecombinationModel::builder()
///     .break_prob(0.01)
///     .crossover_prob(0.8)
///     .gc_extension_prob(0.5)
///     .build()
///     .unwrap();
///
/// assert_eq!(model.break_prob(), 0.01);
/// assert_eq!(model.crossover_prob(), 0.8);
/// ```
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RecombinationModel {
    /// Probability of a double-strand break per base per generation
    break_prob: f64,
    /// Probability that a break is repaired via crossover (vs gene conversion)
    crossover_prob: f64,
    /// Probability that gene conversion extends by one more base
    gc_extension_prob: f64,
    /// Strength of homology preference [0.0, infinity).
    /// 0.0 = Random selection within window.
    /// High = Deterministic best match.
    homology_strength: f64,
    /// Window size (in RUs) to search for homologous sites around the syntenic position.
    search_window: usize,
    /// K-mer size for similarity calculation.
    kmer_size: usize,
}

impl Default for RecombinationModel {
    fn default() -> Self {
        Self::builder().build().unwrap()
    }
}

/// Builder for `RecombinationModel`.
///
/// Allows constructing a `RecombinationModel` with validation.
///
/// # Examples
///
/// ```
/// use centrevo_sim::evolution::RecombinationModel;
///
/// let model = RecombinationModel::builder()
///     .break_prob(0.05)
///     .homology_strength(2.0)
///     .build()
///     .expect("Valid parameters");
/// ```
#[derive(Debug, Clone)]
pub struct RecombinationModelBuilder {
    break_prob: f64,
    crossover_prob: f64,
    gc_extension_prob: f64,
    homology_strength: f64,
    search_window: usize,
    kmer_size: usize,
}

impl Default for RecombinationModelBuilder {
    fn default() -> Self {
        Self {
            break_prob: 0.0,
            crossover_prob: 0.5,
            gc_extension_prob: 0.5,
            homology_strength: 0.0,
            search_window: 100,
            kmer_size: 7,
        }
    }
}

impl RecombinationModelBuilder {
    /// Create a new builder with default values.
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the break probability per base per generation [0.0, 1.0].
    pub fn break_prob(mut self, prob: f64) -> Self {
        self.break_prob = prob;
        self
    }

    /// Set the crossover probability (vs gene conversion) [0.0, 1.0].
    pub fn crossover_prob(mut self, prob: f64) -> Self {
        self.crossover_prob = prob;
        self
    }

    /// Set the gene conversion extension probability [0.0, 1.0].
    pub fn gc_extension_prob(mut self, prob: f64) -> Self {
        self.gc_extension_prob = prob;
        self
    }

    /// Set the homology strength (>= 0.0).
    pub fn homology_strength(mut self, strength: f64) -> Self {
        self.homology_strength = strength;
        self
    }

    /// Set the search window size in RUs.
    pub fn search_window(mut self, window: usize) -> Self {
        self.search_window = window;
        self
    }

    /// Set the k-mer size for similarity calculation (> 0).
    pub fn kmer_size(mut self, size: usize) -> Self {
        self.kmer_size = size;
        self
    }

    /// Build the `RecombinationModel` with validation.
    ///
    /// # Errors
    /// Returns an error if any probability is outside [0.0, 1.0] or other params invalid.
    pub fn build(self) -> Result<RecombinationModel, RecombinationError> {
        if !(0.0..=1.0).contains(&self.break_prob) {
            return Err(RecombinationError::InvalidProbability(
                "break_prob",
                self.break_prob,
            ));
        }
        if !(0.0..=1.0).contains(&self.crossover_prob) {
            return Err(RecombinationError::InvalidProbability(
                "crossover_prob",
                self.crossover_prob,
            ));
        }
        if !(0.0..=1.0).contains(&self.gc_extension_prob) {
            return Err(RecombinationError::InvalidProbability(
                "gc_extension_prob",
                self.gc_extension_prob,
            ));
        }
        if self.homology_strength < 0.0 {
            return Err(RecombinationError::InvalidHomologyStrength(
                self.homology_strength,
            ));
        }
        if self.kmer_size == 0 {
            return Err(RecombinationError::InvalidKmerSize(0));
        }

        Ok(RecombinationModel {
            break_prob: self.break_prob,
            crossover_prob: self.crossover_prob,
            gc_extension_prob: self.gc_extension_prob,
            homology_strength: self.homology_strength,
            search_window: self.search_window,
            kmer_size: self.kmer_size,
        })
    }
}

impl RecombinationModel {
    /// Create a new builder for `RecombinationModel`.
    pub fn builder() -> RecombinationModelBuilder {
        RecombinationModelBuilder::default()
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

    /// Get homology strength.
    #[inline]
    pub fn homology_strength(&self) -> f64 {
        self.homology_strength
    }

    /// Get search window size.
    #[inline]
    pub fn search_window(&self) -> usize {
        self.search_window
    }

    /// Get kmer size.
    #[inline]
    pub fn kmer_size(&self) -> usize {
        self.kmer_size
    }

    /// Sample recombination events for a sequence of given length.
    ///
    /// Returns a list of events that occurred, sorted by position.
    /// This uses a geometric distribution to correctly simulate the per-base
    /// break probability, ensuring that the probability of at least one break
    /// scales correctly with sequence length.
    ///
    /// # Arguments
    /// * `seq1` - The first chromosome (recipient/break source)
    /// * `seq2` - The second chromosome (donor/repair template)
    /// * `rng` - Random number generator
    ///
    /// # Examples
    ///
    /// ```
    /// use centrevo_sim::evolution::RecombinationModel;
    /// use centrevo_sim::genome::{Chromosome, RepeatMap};
    /// use centrevo_sim::base::Sequence;
    /// use rand::SeedableRng;
    /// use rand_xoshiro::Xoshiro256PlusPlus;
    ///
    /// let model = RecombinationModel::builder()
    ///     .break_prob(0.1) // High probability for example
    ///     .build()
    ///     .unwrap();
    ///
    /// let seq: Sequence = "AAAAA".parse().unwrap();
    /// let map = RepeatMap::uniform(1, 1, 5);
    /// let chr1 = Chromosome::new("chr1", seq.clone(), map.clone());
    /// let chr2 = Chromosome::new("chr2", seq, map);
    ///
    /// let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
    /// let events = model.sample_events(&chr1, &chr2, &mut rng);
    /// ```
    pub fn sample_events<R: Rng + ?Sized>(
        &self,
        seq1: &crate::genome::Chromosome,
        seq2: &crate::genome::Chromosome,
        rng: &mut R,
    ) -> Vec<RecombinationType> {
        let mut events = Vec::new();
        let length = seq1.len();
        if length == 0 || self.break_prob <= 0.0 {
            return events;
        }

        // If break_prob is 1.0, we have a break at every position.
        if self.break_prob >= 1.0 {
            for position in 0..length {
                self.add_event_at(position, seq1, seq2, rng, &mut events);
            }
            return events;
        }

        let geo = Geometric::new(self.break_prob).unwrap();

        let mut current_pos = 0;
        loop {
            let skip: u64 = rng.sample(geo);
            current_pos += skip as usize;

            if current_pos >= length {
                break;
            }

            self.add_event_at(current_pos, seq1, seq2, rng, &mut events);

            // Move past the current break for the next iteration
            current_pos += 1;
        }

        events
    }

    // Helper to add an event at a specific position
    fn add_event_at<R: Rng + ?Sized>(
        &self,
        pos1: usize,
        seq1: &crate::genome::Chromosome,
        seq2: &crate::genome::Chromosome,
        rng: &mut R,
        events: &mut Vec<RecombinationType>,
    ) {
        // Find homologous position in seq2
        // NOTE: `pos1` refers to a position in `seq1` (source/donor) and `pos2`
        // refers to the mapped position in `seq2` (target/recipient). The
        // `RecombinationType::GeneConversion` event stores `donor_start` and `recipient_start`
        // corresponding to these positions respectively.
        let pos2 = self.find_homologous_site(seq1, pos1, seq2, rng);

        // Determine if crossover or gene conversion
        if rng.random::<f64>() < self.crossover_prob {
            events.push(RecombinationType::Crossover { pos1, pos2 });
        } else {
            // Gene conversion: sample tract length
            let mut tract_length = 1;
            // Check boundaries for both chromosomes
            while (pos1 + tract_length) < seq1.len()
                && (pos2 + tract_length) < seq2.len()
                && rng.random::<f64>() < self.gc_extension_prob
            {
                tract_length += 1;
            }

            // Store gene conversion event with donor_start=pos1 and recipient_start=pos2
            events.push(RecombinationType::GeneConversion {
                donor_start: pos1,
                recipient_start: pos2,
                length: tract_length,
            });
        }
    }

    /// Find a homologous site in the target chromosome.
    fn find_homologous_site<R: Rng + ?Sized>(
        &self,
        source: &crate::genome::Chromosome,
        source_pos: usize,
        target: &crate::genome::Chromosome,
        rng: &mut R,
    ) -> usize {
        // 1. Identify Source RU
        let source_ru_idx = match source.map().find_ru_index(source_pos) {
            Some(idx) => idx,
            None => return (source_pos * target.len()) / source.len(), // Fallback: proportional mapping
        };

        // 2. Calculate Syntenic Center in Target
        // Map source RU index to target RU index proportionally
        let source_num_rus = source.map().num_rus();
        let target_num_rus = target.map().num_rus();

        if source_num_rus == 0 || target_num_rus == 0 {
            return 0;
        }

        let syntenic_ru_idx = (source_ru_idx * target_num_rus) / source_num_rus;

        // 3. Define Window
        let window = self.search_window;
        let start_idx = syntenic_ru_idx.saturating_sub(window);
        let end_idx = (syntenic_ru_idx + window).min(target_num_rus - 1);

        // 4. Score Candidates
        let mut candidates = Vec::with_capacity(end_idx - start_idx + 1);
        let mut total_weight = 0.0;

        for idx in start_idx..=end_idx {
            let similarity =
                source.calculate_similarity(source_ru_idx, target, idx, self.kmer_size);

            // Weight = Similarity ^ Strength
            let weight = if self.homology_strength == 0.0 {
                1.0 // Uniform
            } else {
                similarity.powf(self.homology_strength)
            };

            if weight > 0.0 {
                candidates.push((idx, weight));
                total_weight += weight;
            }
        }

        // 5. Sample Target RU
        let target_ru_idx = if candidates.is_empty() {
            syntenic_ru_idx // Fallback
        } else {
            let mut choice = rng.random::<f64>() * total_weight;
            let mut selected = candidates.last().unwrap().0;
            for (idx, weight) in candidates {
                choice -= weight;
                if choice <= 0.0 {
                    selected = idx;
                    break;
                }
            }
            selected
        };

        // 6. Map offset within Source RU to Target RU
        let (s_start, s_end) = source.map().get_ru_interval(source_ru_idx).unwrap();
        let (t_start, t_end) = target.map().get_ru_interval(target_ru_idx).unwrap();

        let offset = source_pos - s_start;
        let s_len = s_end - s_start;
        let t_len = t_end - t_start;

        if s_len == 0 {
            return t_start;
        }

        // Proportional offset mapping
        let t_offset = (offset * t_len) / s_len;
        t_start + t_offset
    }

    /// Perform crossover between two sequences at the given position.
    ///
    /// Creates two offspring by swapping sequence content after the crossover point.
    ///
    /// # Returns
    /// A tuple of two new sequences (offspring1, offspring2).
    ///
    /// # Errors
    /// Returns an error if the provided positions are invalid for the given
    /// sequences. Crossover does not require the sequences to be the same total length.
    pub fn crossover(
        &self,
        seq1: &crate::genome::Chromosome,
        seq2: &crate::genome::Chromosome,
        pos1: usize,
        pos2: usize,
    ) -> Result<(crate::genome::Chromosome, crate::genome::Chromosome), RecombinationError> {
        seq1.crossover(seq2, pos1, pos2).map_err(|_| {
            RecombinationError::InvalidPosition {
                position: pos1,
                length: seq1.len(),
            }
        }) // Simplified error mapping
    }

    /// Perform gene conversion by copying a tract from donor to recipient.
    ///
    /// Copies the sequence content from `donor` to `recipient` in the range [start, end).
    ///
    /// # Returns
    /// A new sequence with the converted tract.
    ///
    /// # Errors
    /// Returns an error if the provided indices are invalid for the given
    /// chromosomes. Gene conversion does not require the chromosomes to be
    /// the same total length.
    pub fn gene_conversion(
        &self,
        recipient: &crate::genome::Chromosome,
        donor: &crate::genome::Chromosome,
        recipient_start: usize,
        donor_start: usize,
        length: usize,
    ) -> Result<crate::genome::Chromosome, RecombinationError> {
        recipient
            .gene_conversion(donor, recipient_start, donor_start, length)
            .map_err(|_| RecombinationError::InvalidRange {
                start: recipient_start,
                end: recipient_start + length,
            }) // Simplified error mapping
    }
}

// Removed RecombinationError definition, imported from crate::errors

#[cfg(test)]
mod tests {
    use super::*;
    use crate::base::Sequence;
    use crate::genome::{Chromosome, RepeatMap};
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256PlusPlus;
    use std::str::FromStr;

    fn make_test_chr(seq_str: &str) -> Chromosome {
        let seq = Sequence::from_str(seq_str).unwrap();
        // Uniform map 1bp RUs for simplicity
        let map = RepeatMap::uniform(1, 1, seq_str.len());
        Chromosome::new("chr1", seq, map)
    }

    #[test]
    fn test_recombination_model_builder() {
        let model = RecombinationModel::builder()
            .break_prob(0.01)
            .crossover_prob(0.5)
            .gc_extension_prob(0.1)
            .homology_strength(1.0)
            .search_window(5)
            .kmer_size(7)
            .build()
            .unwrap();

        assert_eq!(model.break_prob(), 0.01);
        assert_eq!(model.crossover_prob(), 0.5);
        assert_eq!(model.gc_extension_prob(), 0.1);
        assert_eq!(model.homology_strength(), 1.0);
        assert_eq!(model.search_window(), 5);
        assert_eq!(model.kmer_size(), 7);
    }

    #[test]
    fn test_recombination_model_invalid() {
        assert!(
            RecombinationModel::builder()
                .break_prob(-0.1)
                .build()
                .is_err()
        );
        assert!(
            RecombinationModel::builder()
                .crossover_prob(1.5)
                .build()
                .is_err()
        );
        assert!(
            RecombinationModel::builder()
                .gc_extension_prob(-0.1)
                .build()
                .is_err()
        );
        assert!(
            RecombinationModel::builder()
                .homology_strength(-1.0)
                .build()
                .is_err()
        );
        assert!(RecombinationModel::builder().kmer_size(0).build().is_err());
    }

    #[test]
    fn test_sample_event_zero_break_prob() {
        let model = RecombinationModel::builder()
            .break_prob(0.0)
            .crossover_prob(0.5)
            .gc_extension_prob(0.1)
            .build()
            .unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        let chr1 = make_test_chr("A".repeat(100).as_str());
        let chr2 = make_test_chr("T".repeat(100).as_str());

        for _ in 0..100 {
            let events = model.sample_events(&chr1, &chr2, &mut rng);
            assert!(events.is_empty());
        }
    }

    #[test]
    fn test_sample_event_empty_sequence() {
        let model = RecombinationModel::builder()
            .break_prob(0.5)
            .build()
            .unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        let chr1 = make_test_chr("");
        let chr2 = make_test_chr("");

        let events = model.sample_events(&chr1, &chr2, &mut rng);
        assert!(events.is_empty());
    }

    #[test]
    fn test_sample_event_types() {
        let model = RecombinationModel::builder()
            .break_prob(1.0)
            .crossover_prob(1.0)
            .build()
            .unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        let chr1 = make_test_chr("AAAAAAAAAA"); // Len 10
        let chr2 = make_test_chr("TTTTTTTTTT");

        // With break_prob=1.0 and crossover_prob=1.0, should always get crossover
        // For length 10, we expect 10 events (one at each position)
        let events = model.sample_events(&chr1, &chr2, &mut rng);
        assert_eq!(events.len(), 10);
        for event in events {
            match event {
                RecombinationType::Crossover { pos1, pos2 } => {
                    // With homology_strength=0, pos2 is random or proportional.
                    // Here lengths are same, so likely pos1==pos2 or close.
                    assert!(pos1 < 10);
                    assert!(pos2 < 10);
                }
                _ => panic!("Expected crossover event"),
            }
        }
    }

    #[test]
    fn test_sample_event_gene_conversion() {
        let model = RecombinationModel::builder()
            .break_prob(1.0)
            .crossover_prob(0.0)
            .build()
            .unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        let chr1 = make_test_chr("AAAAAAAAAA");
        let chr2 = make_test_chr("TTTTTTTTTT");

        // With break_prob=1.0, crossover_prob=0.0, should get gene conversion
        let events = model.sample_events(&chr1, &chr2, &mut rng);
        assert!(!events.is_empty());

        for event in events {
            match event {
                RecombinationType::GeneConversion {
                    donor_start,
                    recipient_start,
                    length,
                } => {
                    assert!(donor_start + length <= 10);
                    assert!(recipient_start + length <= 10);
                }
                _ => panic!("Expected gene conversion event"),
            }
        }
    }

    #[test]
    fn test_crossover_basic() {
        let model = RecombinationModel::builder()
            .break_prob(0.01)
            .crossover_prob(0.5)
            .build()
            .unwrap();
        let chr1 = make_test_chr("AAAA");
        let chr2 = make_test_chr("TTTT");

        let (offspring1, offspring2) = model.crossover(&chr1, &chr2, 2, 2).unwrap();

        assert_eq!(offspring1.to_string(), "AATT");
        assert_eq!(offspring2.to_string(), "TTAA");
    }

    #[test]
    fn test_crossover_at_start() {
        let model = RecombinationModel::builder()
            .break_prob(0.01)
            .crossover_prob(0.5)
            .build()
            .unwrap();
        let chr1 = make_test_chr("AAAA");
        let chr2 = make_test_chr("TTTT");

        let (offspring1, offspring2) = model.crossover(&chr1, &chr2, 0, 0).unwrap();

        assert_eq!(offspring1.to_string(), "TTTT");
        assert_eq!(offspring2.to_string(), "AAAA");
    }

    #[test]
    fn test_crossover_length_mismatch() {
        let model = RecombinationModel::builder()
            .break_prob(0.01)
            .crossover_prob(1.0)
            .build()
            .unwrap();
        let seq1 = make_test_chr("AAAAA"); // 5
        let seq2 = make_test_chr("TTTTTTTTTT"); // 10

        // Crossover at 2 in seq1 and 2 in seq2
        // seq1: AA|AAA, seq2: TT|TTTTTTTT
        // new1: AA|TTTTTTTT (10)
        // new2: TT|AAA (5)
        let result = model.crossover(&seq1, &seq2, 2, 2);
        assert!(result.is_ok());
        let (new1, new2) = result.unwrap();
        assert_eq!(new1.len(), 10);
        assert_eq!(new2.len(), 5);
        assert_eq!(new1.sequence().to_string(), "AATTTTTTTT");
        assert_eq!(new2.sequence().to_string(), "TTAAA");
    }

    #[test]
    fn test_crossover_invalid_position() {
        let model = RecombinationModel::builder()
            .break_prob(0.01)
            .crossover_prob(0.5)
            .build()
            .unwrap();
        let chr1 = make_test_chr("AAAA");
        let chr2 = make_test_chr("TTTT");

        let result = model.crossover(&chr1, &chr2, 10, 10);
        assert!(result.is_err());
    }

    #[test]
    fn test_gene_conversion_basic() {
        let model = RecombinationModel::builder()
            .break_prob(0.01)
            .crossover_prob(0.0)
            .build()
            .unwrap();
        let recipient = make_test_chr("AAAAA");
        let donor = make_test_chr("TTTTT");

        // Replace index 1 (len 2) with donor's index 1
        // Recipient: A|AA|AA -> A|TT|AA
        let result = model.gene_conversion(&recipient, &donor, 1, 1, 2);
        assert!(result.is_ok());
        let new_seq = result.unwrap();
        assert_eq!(new_seq.sequence().to_string(), "ATTAA");
    }

    #[test]
    fn test_gene_conversion_full() {
        let model = RecombinationModel::builder()
            .break_prob(0.01)
            .crossover_prob(0.0)
            .build()
            .unwrap();
        let recipient = make_test_chr("AAAAA");
        let donor = make_test_chr("TTTTT");

        let result = model.gene_conversion(&recipient, &donor, 0, 0, 5);
        assert!(result.is_ok());
        assert_eq!(result.unwrap().sequence().to_string(), "TTTTT");
    }

    #[test]
    fn test_gene_conversion_single_base() {
        let model = RecombinationModel::builder()
            .break_prob(0.01)
            .crossover_prob(0.0)
            .build()
            .unwrap();
        let recipient = make_test_chr("AAAAA");
        let donor = make_test_chr("TTTTT");

        let result = model.gene_conversion(&recipient, &donor, 2, 2, 1);
        assert!(result.is_ok());
        assert_eq!(result.unwrap().sequence().to_string(), "AATAA");
    }

    #[test]
    fn test_gene_conversion_length_mismatch() {
        let model = RecombinationModel::builder()
            .break_prob(0.01)
            .crossover_prob(0.0)
            .build()
            .unwrap();
        let recipient = make_test_chr("AAAAA"); // 5
        let donor = make_test_chr("TTTTTTTTTT"); // 10

        // Replace index 1 (len 2) in recipient with index 5 (len 2) from donor
        // Recipient: A|AA|AA -> A|TT|AA
        // Donor: TTTTT|TT|TTT
        let result = model.gene_conversion(&recipient, &donor, 1, 5, 2);
        assert!(result.is_ok());
        assert_eq!(result.unwrap().sequence().to_string(), "ATTAA");
    }

    #[test]
    fn test_gene_conversion_invalid_range() {
        let model = RecombinationModel::builder()
            .break_prob(0.01)
            .crossover_prob(0.0)
            .build()
            .unwrap();
        let recipient = make_test_chr("AAAAA"); // len 5
        let donor = make_test_chr("TTTTT"); // len 5

        // start + length > recipient.len()
        // 4 + 2 = 6 > 5
        let result = model.gene_conversion(&recipient, &donor, 4, 0, 2);
        assert!(result.is_err());

        // start + length > donor.len()
        let recipient_long = make_test_chr("AAAAA");
        let donor_short = make_test_chr("TTT");
        let result = model.gene_conversion(&recipient_long, &donor_short, 1, 1, 3);
        assert!(result.is_err());
    }

    #[test]
    fn test_gene_conversion_invalid_position() {
        let model = RecombinationModel::builder()
            .break_prob(0.01)
            .crossover_prob(0.5)
            .build()
            .unwrap();
        let recipient = make_test_chr("AAAA");
        let donor = make_test_chr("TTTT");

        let result = model.gene_conversion(&recipient, &donor, 10, 10, 1);
        assert!(result.is_err());

        let result = model.gene_conversion(&recipient, &donor, 0, 10, 1);
        assert!(result.is_err());
    }

    #[test]
    fn test_recombination_type_equality() {
        assert_eq!(RecombinationType::None, RecombinationType::None);
        assert_eq!(
            RecombinationType::Crossover { pos1: 10, pos2: 10 },
            RecombinationType::Crossover { pos1: 10, pos2: 10 }
        );
        assert_ne!(
            RecombinationType::Crossover { pos1: 10, pos2: 10 },
            RecombinationType::Crossover { pos1: 10, pos2: 20 }
        );
        assert_eq!(
            RecombinationType::GeneConversion {
                donor_start: 1,
                recipient_start: 1,
                length: 5
            },
            RecombinationType::GeneConversion {
                donor_start: 1,
                recipient_start: 1,
                length: 5
            }
        );
        assert_ne!(
            RecombinationType::GeneConversion {
                donor_start: 1,
                recipient_start: 1,
                length: 5
            },
            RecombinationType::GeneConversion {
                donor_start: 1,
                recipient_start: 2,
                length: 5
            }
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
        let model1 = RecombinationModel::builder()
            .break_prob(0.01)
            .crossover_prob(0.5)
            .gc_extension_prob(0.1)
            .homology_strength(1.0)
            .search_window(5)
            .kmer_size(7)
            .build()
            .unwrap();
        let model2 = model1.clone();

        assert_eq!(model1.break_prob(), model2.break_prob());
        assert_eq!(model1.crossover_prob(), model2.crossover_prob());
        assert_eq!(model1.gc_extension_prob(), model2.gc_extension_prob());
    }

    #[test]
    fn test_sample_events_multiple() {
        // High break prob, long sequence -> multiple events
        let model = RecombinationModel::builder()
            .break_prob(0.1)
            .crossover_prob(0.5)
            .build()
            .unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        let chr1 = make_test_chr("A".repeat(1000).as_str());
        let chr2 = make_test_chr("T".repeat(1000).as_str());

        let events = model.sample_events(&chr1, &chr2, &mut rng);
        // Expected events ~ 100
        assert!(events.len() > 50);
        assert!(events.len() < 150);

        // Verify sorted positions
        let mut last_pos = 0;
        for event in events {
            let pos = match event {
                RecombinationType::Crossover { pos1, .. } => pos1,
                RecombinationType::GeneConversion { donor_start, .. } => donor_start,
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

        let model = RecombinationModel::builder()
            .break_prob(0.01)
            .crossover_prob(0.5)
            .build()
            .unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        let chr1 = make_test_chr("A".repeat(100).as_str());
        let chr2 = make_test_chr("T".repeat(100).as_str());

        let mut event_count = 0;
        let trials = 10000;
        for _ in 0..trials {
            let events = model.sample_events(&chr1, &chr2, &mut rng);
            if !events.is_empty() {
                event_count += 1;
            }
        }

        let frequency = event_count as f64 / trials as f64;
        println!("Frequency of events for L=100, p=0.01: {frequency}");

        // Now we expect correct scaling
        assert!(
            frequency > 0.60 && frequency < 0.66,
            "Frequency {frequency} should be around 0.63"
        );
    }
}
