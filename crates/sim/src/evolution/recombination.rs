//! Recombination operations for sequences.
//!
//! This module simulates sexual reproduction through recombination, the process by
//! which genetic material from two parents is shuffled to create offspring with new
//! combinations of alleles. This is a fundamental mechanism driving genetic diversity
//! and the efficacy of natural selection.
//!
//! ## Recombination Types
//!
//! ### Crossover (Reciprocal Recombination)
//! DNA strands from two chromosomes break at the same (or different) positions and
//! rejoin in crossed patterns. This results in both chromosomes having new combinations
//! of alleles. Crossovers are balanced - both parents lose sequence at the break point
//! but gain sequence from the other parent.
//!
//! ### Gene Conversion (Non-Reciprocal Recombination)
//! A segment from the donor chromosome is copied to the recipient chromosome, replacing
//! the local alleles. Unlike crossover, this is non-reciprocal: the donor is unchanged.
//! Gene conversion is thought to be the mechanism for "homogenization" of tandem repeats.
//!
//! ## How Recombination Works in This Model
//!
//! 1. **Break Sites**: Meiosis creates double-strand breaks at random positions (Poisson process)
//! 2. **Homology Search**: The break site searches for matching sequence in the homologous chromosome
//! 3. **Repair Choice**: The break is repaired as either crossover or gene conversion
//! 4. **Tract Extension**: Gene conversion tracts extend stochastically along the chromosome

pub use crate::errors::RecombinationError;
use rand::Rng;
use rand_distr::Geometric;
use serde::{Deserialize, Serialize};

/// Type of recombination event that can occur on a sequence.
///
/// During meiosis, when two homologous chromosomes pair up, recombination events
/// can occur to shuffle genetic material. This enum describes what happened at each
/// recombination point:
///
/// - **None**: No recombination (used for neutral regions or modeling)
/// - **Crossover**: Reciprocal exchange at two positions (one per chromosome)
/// - **Gene Conversion**: Unidirectional copying from donor to recipient chromosome
///
/// Biologically, these represent different mechanistic outcomes of meiotic DSBs
/// (double-strand breaks) during prophase I of meiosis.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum RecombinationType {
    /// No recombination occurs
    None,
    /// Crossover (recombination) at specific positions in both chromosomes
    Crossover { pos1: usize, pos2: usize },
    /// Gene conversion from start to end position
    ///
    /// Biologically: A "copy-paste" event where sequence from the **donor** is copied
    /// to the **recipient**, overwriting its original sequence. This homogenizes sequences,
    /// reducing differences between repeats.
    ///
    /// - `donor_start`: Start position in the donor chromosome
    /// - `recipient_start`: Start position in the recipient chromosome
    /// - `length`: Length of the tract copied
    GeneConversion {
        donor_start: usize,
        recipient_start: usize,
        length: usize,
    },
}

/// Parameters controlling recombination behavior.
///
/// These parameters model the key biological aspects of meiotic recombination:
/// - How frequently breaks occur (break_prob)
/// - Whether breaks resolve as crossovers or gene conversions (crossover_prob)
/// - How long gene conversion tracts extend (gc_extension_prob)
/// - How sensitive the cell is to homology when choosing repair templates (homology_strength)
/// - Where the cell searches for homologous sequences (search_window, kmer_size)
///
/// The parameters reflect mechanistic details observed in real meiosis:
/// 1. **Break Frequency**: Determined by meiotic machinery; ~1-3 breaks per chromosome arm
/// 2. **Break Resolution**: DSBs can be repaired as crossovers (~1 per arm, constrained) or
///    gene conversions (many per arm, homology-dependent)
/// 3. **Homology Search**: Meiotic complexes search locally for homologous sequence
/// 4. **GC Tract Length**: Real gene conversion tracts average ~1-4 kb in fungi, ~0.5-2 kb in mammals
///
/// Use `RecombinationModel::builder()` to construct a new instance with validation.
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
    /// This simulates meiosis for a chromosome: we generate break sites stochastically
    /// using a Bernoulli process (which approximates a Poisson process for rare events).
    ///
    /// The process mimics real meiosis:
    /// 1. Meiotic DSBs occur at random intervals along the chromosome. We sample the
    ///    distance between breaks using a Geometric distribution, which describes the
    ///    inter-arrival time of a Bernoulli process.
    /// 2. RAD51/DMC1 proteins coat the single-stranded DNA and search for homology
    /// 3. Homology-directed strand invasion initiates repair
    /// 4. Breaks are resolved as crossovers (~50% via BLM helicase dissolution) or
    ///    gene conversions (~50% via nuclease resolution)
    ///
    /// # Arguments
    /// * `seq1` - The first chromosome (recipient for gene conversion, one strand of crossover)
    /// * `seq2` - The second chromosome (donor for gene conversion, homolog for pairing)
    /// * `rng` - Random number generator
    ///
    /// # Returns
    /// A list of RecombinationType events, sorted by position in seq1. The events
    /// represent what actually happened to the sequence during meiosis.
    ///
    /// # Examples
    ///
    /// ```
    /// use centrevo_sim::evolution::RecombinationModel;
    /// use centrevo_sim::genome::{Chromosome, RepeatMap};
    /// use centrevo_sim::base::Sequence;
    /// use centrevo_sim::base::Nucleotide;
    /// use rand::SeedableRng;
    /// use rand_xoshiro::Xoshiro256PlusPlus;
    ///
    /// use centrevo_sim::base::GenomeArena;
    ///
    /// let mut arena = GenomeArena::new();
    /// let model = RecombinationModel::builder()
    ///     .break_prob(0.1) // High probability for example
    ///     .build()
    ///     .unwrap();
    ///
    /// let nucs: Vec<Nucleotide> = "AAAAA".chars()
    ///     .map(|c| Nucleotide::from_ascii(c as u8).unwrap())
    ///     .collect();
    /// let map = RepeatMap::uniform(1, 1, 5);
    /// // Allocate sequence in arena
    /// let data = arena.alloc(&nucs);
    /// let chr1 = Chromosome::new("chr1", data, map.clone());
    /// let chr2 = Chromosome::new("chr2", data, map);
    ///
    /// let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
    /// let events = model.sample_events(&chr1, &chr2, &arena, &mut rng);
    /// ```
    pub fn sample_events<R: Rng + ?Sized>(
        &self,
        seq1: &crate::genome::Chromosome,
        seq2: &crate::genome::Chromosome,
        arena: &crate::base::GenomeArena,
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
                self.add_event_at(position, seq1, seq2, arena, rng, &mut events);
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

            self.add_event_at(current_pos, seq1, seq2, arena, rng, &mut events);

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
        arena: &crate::base::GenomeArena,
        rng: &mut R,
        events: &mut Vec<RecombinationType>,
    ) {
        // Find homologous position in seq2
        // NOTE: `pos1` refers to a position in `seq1` (source/donor) and `pos2`
        // refers to the mapped position in `seq2` (target/recipient). The
        // `RecombinationType::GeneConversion` event stores `donor_start` and `recipient_start`
        // corresponding to these positions respectively.
        let pos2 = self.find_homologous_site(seq1, pos1, seq2, arena, rng);

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
    ///
    /// This models the strand invasion and homology search process during meiotic
    /// recombination. After a DSB, the 5'-3' exonuclease degrades one strand,
    /// creating a 3' single-stranded tail. RAD51 and DMC1 proteins then:
    ///
    /// 1. Load onto the 3' tail
    /// 2. Probe the partner chromosome for homologous sequence
    /// 3. Search a local window around the syntenic position (expected position based
    ///    on chromosome alignment)
    /// 4. Score candidates by similarity (k-mer matching)
    /// 5. Weight candidates by homology strength (biological: stronger homology = preference)
    /// 6. Sample a site proportional to these weights
    /// 7. Invade and prime DNA synthesis
    ///
    /// The algorithm implements weighted random sampling: higher-similarity sequences
    /// are more likely to be chosen, but with randomness reflecting biological variance
    /// in the meiotic machinery.
    fn find_homologous_site<R: Rng + ?Sized>(
        &self,
        source: &crate::genome::Chromosome,
        source_pos: usize,
        target: &crate::genome::Chromosome,
        arena: &crate::base::GenomeArena,
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

        // OPTIMIZATION: Fast path for Random Homology (strength == 0.0)
        // If strength is 0, weights are all 1.0 (Uniform).
        // We can just sample a random index in [start_idx, end_idx]
        if self.homology_strength == 0.0 {
            if start_idx > end_idx {
                return syntenic_ru_idx;
            }
            let target_ru_idx = rng.random_range(start_idx..=end_idx);

            // Map offset within Source RU to Target RU
            return self.map_offset_to_target(
                source,
                source_pos,
                source_ru_idx,
                target,
                target_ru_idx,
            );
        }

        // 4. Score Candidates
        let mut candidates = Vec::with_capacity(end_idx - start_idx + 1);
        let mut total_weight = 0.0;

        for idx in start_idx..=end_idx {
            let similarity =
                source.calculate_similarity(source_ru_idx, target, idx, self.kmer_size, arena);

            // Weight = Similarity ^ Strength
            // We already handled strength == 0.0 above
            let weight = similarity.powf(self.homology_strength);

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
        self.map_offset_to_target(source, source_pos, source_ru_idx, target, target_ru_idx)
    }

    /// Helper to map the offset from a source RU to a target RU.
    fn map_offset_to_target(
        &self,
        source: &crate::genome::Chromosome,
        source_pos: usize,
        source_ru_idx: usize,
        target: &crate::genome::Chromosome,
        target_ru_idx: usize,
    ) -> usize {
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
    /// Crossover (recombination) creates two recombinant chromosomes by reciprocal
    /// exchange of sequence. If seq1 is broken at pos1 and seq2 is broken at pos2,
    /// the results are:
    /// - Offspring 1: [seq1 before pos1] + [seq2 from pos2]
    /// - Offspring 2: [seq2 before pos2] + [seq1 from pos1]
    ///
    /// Biologically: both offspring inherit new combinations of alleles around the
    /// crossover point, increasing genetic diversity. For two parents of equal fitness,
    /// crossover produces balanced offspring (no allele loss).
    ///
    /// # Returns
    /// A tuple of two new sequences (offspring1, offspring2).
    ///
    /// Returns an error if the provided positions are invalid for the given
    /// sequences. Crossover does not require the sequences to be the same total length.
    ///
    /// # Examples
    ///
    /// ```
    /// use centrevo_sim::evolution::RecombinationModel;
    /// use centrevo_sim::genome::{Chromosome, RepeatMap};
    /// use centrevo_sim::base::{Nucleotide, GenomeArena};
    ///
    /// let mut arena = GenomeArena::new();
    /// let model = RecombinationModel::builder().build().unwrap();
    ///
    /// // Create two chromosomes: "AAAA" and "TTTT"
    /// let seq_a: Vec<Nucleotide> = vec![Nucleotide::A; 4];
    /// let seq_t: Vec<Nucleotide> = vec![Nucleotide::T; 4];
    /// let map = RepeatMap::uniform(1, 4, 1);
    ///
    /// let data_a = arena.alloc(&seq_a);
    /// let chr_a = Chromosome::new("A", data_a, map.clone());
    ///
    /// let data_t = arena.alloc(&seq_t);
    /// let chr_t = Chromosome::new("T", data_t, map);
    ///
    /// // Crossover at index 2 (midpoint) implies:
    /// // Offspring 1 takes first 2 of A ("AA") and rest of T ("TT") -> "AATT"
    /// // Offspring 2 takes first 2 of T ("TT") and rest of A ("AA") -> "TTAA"
    /// let (off1, off2) = model.crossover(&chr_a, &chr_t, 2, 2, &mut arena).unwrap();
    ///
    /// assert_eq!(off1.to_formatted_string("", "", &arena), "AATT");
    /// assert_eq!(off2.to_formatted_string("", "", &arena), "TTAA");
    /// ```
    pub fn crossover(
        &self,
        seq1: &crate::genome::Chromosome,
        seq2: &crate::genome::Chromosome,
        pos1: usize,
        pos2: usize,
        arena: &mut crate::base::GenomeArena,
    ) -> Result<(crate::genome::Chromosome, crate::genome::Chromosome), RecombinationError> {
        seq1.crossover(seq2, pos1, pos2, arena)
            .map_err(|_| RecombinationError::InvalidPosition {
                position: pos1,
                length: seq1.len(),
            }) // Simplified error mapping
    }

    /// Perform gene conversion by copying a tract from donor to recipient.
    ///
    /// Gene conversion replaces a segment of the recipient chromosome with the
    /// corresponding segment from the donor. Unlike crossover (reciprocal), only
    /// the recipient changes; the donor is unmodified.
    ///
    /// Biologically: gene conversion arises from asymmetric resolution of recombination
    /// intermediates. The tract length reflects processivity of meiotic exonucleases:
    /// longer tracts = more processive exonucleases; shorter tracts = pausing/dissociation.
    /// Gene conversion can homogenize sequences within
    /// tandem repeats, reducing sequence divergence between copies (important for
    /// maintaining repeat unit similarity).
    ///
    /// # Returns
    /// A new sequence with the converted tract replaced by donor sequence.
    ///
    /// # Errors
    /// Returns an error if the provided indices are invalid for the given
    /// chromosomes. Gene conversion does not require the chromosomes to be
    /// the same total length.
    ///
    /// # Examples
    ///
    /// ```
    /// use centrevo_sim::evolution::RecombinationModel;
    /// use centrevo_sim::genome::{Chromosome, RepeatMap};
    /// use centrevo_sim::base::{Nucleotide, GenomeArena};
    ///
    /// let mut arena = GenomeArena::new();
    /// let model = RecombinationModel::builder().build().unwrap();
    ///
    /// // Recipient: "AAAAA"
    /// let seq_r: Vec<Nucleotide> = vec![Nucleotide::A; 5];
    /// // Donor: "TTTTT"
    /// let seq_d: Vec<Nucleotide> = vec![Nucleotide::T; 5];
    /// let map = RepeatMap::uniform(1, 5, 1);
    ///
    /// let data_r = arena.alloc(&seq_r);
    /// let chr_r = Chromosome::new("Rec", data_r, map.clone());
    ///
    /// let data_d = arena.alloc(&seq_d);
    /// let chr_d = Chromosome::new("Don", data_d, map);
    ///
    /// // Copy 2 bases from Donor starting at index 1 to Recipient starting at index 1
    /// // Recipient becomes: A (orig) + TT (copied) + AA (orig) = "ATTAA"
    /// let res = model.gene_conversion(&chr_r, &chr_d, 1, 1, 2, &mut arena).unwrap();
    ///
    /// assert_eq!(res.to_formatted_string("", "", &arena), "ATTAA");
    /// ```
    pub fn gene_conversion(
        &self,
        recipient: &crate::genome::Chromosome,
        donor: &crate::genome::Chromosome,
        recipient_start: usize,
        donor_start: usize,
        length: usize,
        arena: &mut crate::base::GenomeArena,
    ) -> Result<crate::genome::Chromosome, RecombinationError> {
        recipient
            .gene_conversion(donor, recipient_start, donor_start, length, arena)
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
    use crate::base::{GenomeArena, Nucleotide}; // Added GenomeArena
    use crate::genome::{Chromosome, RepeatMap};
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256PlusPlus;
    // use crate::base::Sequence; // Removed unused

    fn make_test_chr(seq_str: &str, arena: &mut GenomeArena) -> Chromosome {
        let nucs: Vec<Nucleotide> = seq_str
            .chars()
            .map(|c| Nucleotide::from_ascii(c as u8).unwrap())
            .collect();
        let data = arena.alloc(&nucs);
        let map = RepeatMap::uniform(1, 1, seq_str.len());
        Chromosome::new("chr1", data, map)
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
        let mut arena = GenomeArena::new();
        let model = RecombinationModel::builder()
            .break_prob(0.0)
            .crossover_prob(0.5)
            .gc_extension_prob(0.1)
            .build()
            .unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        let chr1 = make_test_chr(&"A".repeat(100), &mut arena); // .as_str() removed
        let chr2 = make_test_chr(&"T".repeat(100), &mut arena);

        for _ in 0..100 {
            let events = model.sample_events(&chr1, &chr2, &arena, &mut rng);
            assert!(events.is_empty());
        }
    }

    #[test]
    fn test_sample_event_empty_sequence() {
        let mut arena = GenomeArena::new();
        let model = RecombinationModel::builder()
            .break_prob(0.5)
            .build()
            .unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        let chr1 = make_test_chr("", &mut arena);
        let chr2 = make_test_chr("", &mut arena);

        let events = model.sample_events(&chr1, &chr2, &arena, &mut rng);
        assert!(events.is_empty());
    }

    #[test]
    fn test_sample_event_types() {
        let mut arena = GenomeArena::new();
        let model = RecombinationModel::builder()
            .break_prob(1.0)
            .crossover_prob(1.0)
            .build()
            .unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        let chr1 = make_test_chr("AAAAAAAAAA", &mut arena); // Len 10
        let chr2 = make_test_chr("TTTTTTTTTT", &mut arena);

        // With break_prob=1.0 and crossover_prob=1.0, should always get crossover
        // For length 10, we expect 10 events (one at each position)
        let events = model.sample_events(&chr1, &chr2, &arena, &mut rng);
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
        let mut arena = GenomeArena::new();
        let model = RecombinationModel::builder()
            .break_prob(1.0)
            .crossover_prob(0.0)
            .build()
            .unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        let chr1 = make_test_chr("AAAAAAAAAA", &mut arena);
        let chr2 = make_test_chr("TTTTTTTTTT", &mut arena);

        // With break_prob=1.0, crossover_prob=0.0, should get gene conversion
        let events = model.sample_events(&chr1, &chr2, &arena, &mut rng);
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
        let mut arena = GenomeArena::new();
        let model = RecombinationModel::builder()
            .break_prob(0.01)
            .crossover_prob(0.5)
            .build()
            .unwrap();
        let chr1 = make_test_chr("AAAA", &mut arena);
        let chr2 = make_test_chr("TTTT", &mut arena);

        let (offspring1, offspring2) = model.crossover(&chr1, &chr2, 2, 2, &mut arena).unwrap();

        assert_eq!(
            offspring1
                .to_formatted_string(" ", " ", &arena)
                .replace(" ", ""),
            "AATT"
        );
        assert_eq!(
            offspring2
                .to_formatted_string(" ", " ", &arena)
                .replace(" ", ""),
            "TTAA"
        );
    }

    #[test]
    fn test_crossover_at_start() {
        let mut arena = GenomeArena::new();
        let model = RecombinationModel::builder()
            .break_prob(0.01)
            .crossover_prob(0.5)
            .build()
            .unwrap();
        let chr1 = make_test_chr("AAAA", &mut arena);
        let chr2 = make_test_chr("TTTT", &mut arena);

        let (offspring1, offspring2) = model.crossover(&chr1, &chr2, 0, 0, &mut arena).unwrap();

        assert_eq!(
            offspring1
                .to_formatted_string(" ", " ", &arena)
                .replace(" ", ""),
            "TTTT"
        );
        assert_eq!(
            offspring2
                .to_formatted_string(" ", " ", &arena)
                .replace(" ", ""),
            "AAAA"
        );
    }

    #[test]
    fn test_crossover_length_mismatch() {
        let mut arena = GenomeArena::new();
        let model = RecombinationModel::builder()
            .break_prob(0.01)
            .crossover_prob(1.0)
            .build()
            .unwrap();
        let seq1 = make_test_chr("AAAAA", &mut arena); // 5
        let seq2 = make_test_chr("TTTTTTTTTT", &mut arena); // 10

        // Crossover at 2 in seq1 and 2 in seq2
        // seq1: AA|AAA, seq2: TT|TTTTTTTT
        // new1: AA|TTTTTTTT (10)
        // new2: TT|AAA (5)
        let result = model.crossover(&seq1, &seq2, 2, 2, &mut arena);
        assert!(result.is_ok());
        let (new1, new2) = result.unwrap();
        assert_eq!(new1.len(), 10);
        assert_eq!(new2.len(), 5);
        assert_eq!(
            new1.to_formatted_string(" ", " ", &arena).replace(" ", ""),
            "AATTTTTTTT"
        );
        assert_eq!(
            new2.to_formatted_string(" ", " ", &arena).replace(" ", ""),
            "TTAAA"
        );
    }

    #[test]
    fn test_crossover_invalid_position() {
        let mut arena = GenomeArena::new();
        let model = RecombinationModel::builder()
            .break_prob(0.01)
            .crossover_prob(0.5)
            .build()
            .unwrap();
        let chr1 = make_test_chr("AAAA", &mut arena);
        let chr2 = make_test_chr("TTTT", &mut arena);

        let result = model.crossover(&chr1, &chr2, 10, 10, &mut arena);
        assert!(result.is_err());
    }

    #[test]
    fn test_gene_conversion_basic() {
        let mut arena = GenomeArena::new();
        let model = RecombinationModel::builder()
            .break_prob(0.01)
            .crossover_prob(0.0)
            .build()
            .unwrap();
        let recipient = make_test_chr("AAAAA", &mut arena);
        let donor = make_test_chr("TTTTT", &mut arena);

        // Replace index 1 (len 2) with donor's index 1
        // Recipient: A|AA|AA -> A|TT|AA
        let result = model.gene_conversion(&recipient, &donor, 1, 1, 2, &mut arena);
        assert!(result.is_ok());
        let new_seq = result.unwrap();
        assert_eq!(
            new_seq
                .to_formatted_string(" ", " ", &arena)
                .replace(" ", ""),
            "ATTAA"
        );
    }

    #[test]
    fn test_gene_conversion_full() {
        let mut arena = GenomeArena::new();
        let model = RecombinationModel::builder()
            .break_prob(0.01)
            .crossover_prob(0.0)
            .build()
            .unwrap();
        let recipient = make_test_chr("AAAAA", &mut arena);
        let donor = make_test_chr("TTTTT", &mut arena);

        let result = model.gene_conversion(&recipient, &donor, 0, 0, 5, &mut arena);
        assert!(result.is_ok());
        assert_eq!(
            result
                .unwrap()
                .to_formatted_string(" ", " ", &arena)
                .replace(" ", ""),
            "TTTTT"
        );
    }

    #[test]
    fn test_gene_conversion_single_base() {
        let mut arena = GenomeArena::new();
        let model = RecombinationModel::builder()
            .break_prob(0.01)
            .crossover_prob(0.0)
            .build()
            .unwrap();
        let recipient = make_test_chr("AAAAA", &mut arena);
        let donor = make_test_chr("TTTTT", &mut arena);

        let result = model.gene_conversion(&recipient, &donor, 2, 2, 1, &mut arena);
        assert!(result.is_ok());
        assert_eq!(
            result
                .unwrap()
                .to_formatted_string(" ", " ", &arena)
                .replace(" ", ""),
            "AATAA"
        );
    }

    #[test]
    fn test_gene_conversion_length_mismatch() {
        let mut arena = GenomeArena::new();
        let model = RecombinationModel::builder()
            .break_prob(0.01)
            .crossover_prob(0.0)
            .build()
            .unwrap();
        let recipient = make_test_chr("AAAAA", &mut arena); // 5
        let donor = make_test_chr("TTTTTTTTTT", &mut arena); // 10

        // Replace index 1 (len 2) in recipient with index 5 (len 2) from donor
        // Recipient: A|AA|AA -> A|TT|AA
        // Donor: TTTTT|TT|TTT
        let result = model.gene_conversion(&recipient, &donor, 1, 5, 2, &mut arena);
        assert!(result.is_ok());
        assert_eq!(
            result
                .unwrap()
                .to_formatted_string(" ", " ", &arena)
                .replace(" ", ""),
            "ATTAA"
        );
    }

    #[test]
    fn test_gene_conversion_invalid_range() {
        let mut arena = GenomeArena::new();
        let model = RecombinationModel::builder()
            .break_prob(0.01)
            .crossover_prob(0.0)
            .build()
            .unwrap();
        let recipient = make_test_chr("AAAAA", &mut arena); // len 5
        let donor = make_test_chr("TTTTT", &mut arena); // len 5

        // start + length > recipient.len()
        // 4 + 2 = 6 > 5
        let result = model.gene_conversion(&recipient, &donor, 4, 0, 2, &mut arena);
        assert!(result.is_err());

        // start + length > donor.len()
        let recipient_long = make_test_chr("AAAAA", &mut arena);
        let donor_short = make_test_chr("TTT", &mut arena);
        let result = model.gene_conversion(&recipient_long, &donor_short, 1, 1, 3, &mut arena);
        assert!(result.is_err());
    }

    #[test]
    fn test_gene_conversion_invalid_position() {
        let mut arena = GenomeArena::new();
        let model = RecombinationModel::builder()
            .break_prob(0.01)
            .crossover_prob(0.5)
            .build()
            .unwrap();
        let recipient = make_test_chr("AAAA", &mut arena);
        let donor = make_test_chr("TTTT", &mut arena);

        let result = model.gene_conversion(&recipient, &donor, 10, 10, 1, &mut arena);
        assert!(result.is_err());

        let result = model.gene_conversion(&recipient, &donor, 0, 10, 1, &mut arena);
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
        let mut arena = GenomeArena::new();
        let chr1 = make_test_chr(&"A".repeat(1000), &mut arena);
        let chr2 = make_test_chr(&"T".repeat(1000), &mut arena);

        let events = model.sample_events(&chr1, &chr2, &arena, &mut rng);
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
        let mut arena = GenomeArena::new();
        let chr1 = make_test_chr(&"A".repeat(100), &mut arena);
        let chr2 = make_test_chr(&"T".repeat(100), &mut arena);

        let mut event_count = 0;
        let trials = 10000;
        for _ in 0..trials {
            let events = model.sample_events(&chr1, &chr2, &arena, &mut rng);
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
