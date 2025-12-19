//! Selection and fitness functions for individuals.
//!
//! This module provides fitness scoring capabilities for simulating natural selection.
//! Fitness determines reproductive success: individuals with higher fitness contribute
//! more offspring to the next generation, driving allele frequencies toward beneficial
//! variants and maintaining harmful variation at mutation-selection balance.
//!
//! ## Fitness Functions
//!
//! ### Haplotype Fitness (Single Chromosome)
//! Fitness of a single chromosome, typically based on sequence properties:
//! - **GC Content**: Many organisms have GC-content-dependent fitness, reflecting
//!   GC-bias mutations or selection for stability (GC bonds are stronger)
//! - **Length**: Longer sequences might be more robust (more genes) or costly (more to replicate)
//!
//! ### Individual Fitness (Diploid Organism)
//! Fitness of a whole organism with two chromosome copies (alleles). Can be:
//! - **Multiplicative**: `fitness = haplotype1_fitness × haplotype2_fitness` (default, models additive effects)
//! - **Interactive**: `fitness = f(haplotype1, haplotype2)` (models dominance, overdominance, etc.)

use crate::base::{FitnessValue, Nucleotide};
use crate::errors::FitnessError;
use crate::genome::{Chromosome, Individual};
use serde::{Deserialize, Serialize};

/// Trait for scoring fitness of a single haplotype (chromosome).
///
/// Implementors should provide a `haplotype_fitness` method that computes a
/// non-negative fitness value for the haplotype. The default implementation
/// returns neutral fitness of 1.0.
///
/// Biologically: haplotype fitness represents the relative survival and
/// reproductive success of individuals carrying a particular allele or chromosome.
/// Fitness is typically measured relative to a reference (often the wild-type),
/// so values > 1 are beneficial and < 1 are deleterious.
pub trait HaplotypeFitness {
    /// Calculate fitness score for a single haplotype.
    ///
    /// Returns a non-negative fitness value. Higher values indicate better fitness.
    /// The default implementation returns neutral fitness of 1.0.
    fn haplotype_fitness(&self, _chromosome: &Chromosome) -> FitnessValue {
        FitnessValue::NEUTRAL_FITNESS
    }
}

/// Trait for scoring fitness of a diploid individual (two haplotypes).
///
/// Implementors should provide an `individual_fitness` method that computes a
/// non-negative fitness value for the whole individual. The default
/// implementation composes the per-haplotype score by multiplying the two
/// haplotype fitness values. Override `individual_fitness` when the fitness
/// depends on interactions between haplotypes (for example sequence similarity
/// or length similarity).
///
/// Biologically: diploid fitness reflects how the two alleles interact:
/// - **Multiplicative (default)**: Assumes additive/independent effects (common for weak selection)
/// - **Dominance**: One allele masks the effect of another (classic Mendelian pattern)
/// - **Overdominance**: Heterozygotes more fit than either homozygote (heterozygote advantage)
/// - **Underdominance**: Heterozygotes less fit than homozygotes (reproductive incompatibility)
///
/// Example (override default):
///
/// ```rust
/// # use centrevo_sim::evolution::IndividualFitness;
/// # use centrevo_sim::genome::Individual;
/// use centrevo_sim::base::FitnessValue;
///
/// struct MyFitness;
/// impl centrevo_sim::evolution::HaplotypeFitness for MyFitness {}
/// impl IndividualFitness for MyFitness {
///     fn individual_fitness(&self, _ind: &Individual) -> FitnessValue {
///         FitnessValue::new(0.42)
///     }
/// }
/// ```
pub trait IndividualFitness: HaplotypeFitness {
    /// Calculate fitness score for a diploid individual.
    ///
    /// The default implementation multiplies the fitness scores of both haplotypes.
    fn individual_fitness(&self, individual: &Individual) -> FitnessValue {
        let fit1 = self.haplotype_fitness(individual.haplotype1().get(0).unwrap());
        let fit2 = self.haplotype_fitness(individual.haplotype2().get(0).unwrap());
        fit1 * fit2
    }
}

/// GC content-based fitness function.
///
/// Fitness follows a Beta distribution centered on the optimal GC content.
/// This models real biological phenomena:
///
/// 1. **GC-Bias Mutations**: Most organisms have biased mutation patterns favoring GC
///    or AT, shifting genome-wide GC content over time
/// 2. **GC-Dependent Fitness**: GC content affects:
///    - Codon usage bias (translation efficiency)
///    - mRNA stability (GC bonds are stronger)
///    - Thermodynamic stability of DNA secondary structures
///    - Horizontal gene transfer barriers (foreign DNA can have unusual GC)
/// 3. **Selection-Mutation Balance**: Optimal GC reflects balance between mutation bias
///    and stabilizing selection
///
/// The Beta distribution provides a biologically realistic fitness landscape:
/// - Smooth optimum: selection is weak near the optimum
/// - Symmetric: allows modeling any optimum from 0 to 1
/// - Peaked at the optimum: strong selection against deviations
/// - Concentration parameter controls the sharpness (higher = stronger selection)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GCContentFitness {
    /// Optimal GC content (between 0.0 and 1.0)
    pub optimum: f64,
    /// Concentration parameter controlling the sharpness of the fitness curve
    pub concentration: f64,
}

impl GCContentFitness {
    /// Create a new GC content fitness function.
    ///
    /// # Arguments
    /// * `optimum` - Optimal GC content (must be between 0.0 and 1.0, inclusive)
    ///   - 0.0 = AT-rich sequences are optimal (e.g., parasites, thermophiles)
    ///   - 0.5 = GC-neutral (no preference)
    ///   - 1.0 = GC-rich sequences are optimal (e.g., some bacteria, mammalian CpG islands)
    /// * `concentration` - Sharpness of fitness curve (must be > 0.0)
    ///   - 0.5 = very weak selection (broad fitness landscape)
    ///   - 1.0 = moderate selection (typical for real organisms)
    ///   - 5.0+ = very strong selection (narrow tolerance for deviations)
    ///
    /// # Errors
    /// Returns an error if parameters are outside valid ranges.
    pub fn new(optimum: f64, concentration: f64) -> Result<Self, FitnessError> {
        if !(0.0..=1.0).contains(&optimum) {
            return Err(FitnessError::InvalidParameter(
                "optimum must be between 0.0 and 1.0 (inclusive)".into(),
            ));
        }
        if concentration <= 0.0 {
            return Err(FitnessError::InvalidParameter(
                "concentration must be greater than 0.0".into(),
            ));
        }
        Ok(Self {
            optimum,
            concentration,
        })
    }

    /// Convert optimum and concentration to alpha/beta parameters for Beta distribution.
    fn to_alpha_beta(&self) -> (f64, f64) {
        let alpha = 1.0 + self.concentration * self.optimum;
        let beta = 1.0 + self.concentration * (1.0 - self.optimum);
        (alpha, beta)
    }
}

impl HaplotypeFitness for GCContentFitness {
    fn haplotype_fitness(&self, chromosome: &Chromosome) -> FitnessValue {
        let gc_content = chromosome.gc_content();

        // Handle edge cases with small epsilon to avoid log(0)
        let gc = gc_content.clamp(1e-10, 1.0 - 1e-10);

        let (alpha, beta) = self.to_alpha_beta();

        // Beta PDF returns values in [0, 1], but the area under the curve is not 1; use raw PDF values as weights
        FitnessValue::new(gc.powf(alpha - 1.0) * (1.0 - gc).powf(beta - 1.0))
    }
}

impl IndividualFitness for GCContentFitness {}

/// Length-based fitness function.
///
/// Fitness follows a log-normal distribution centered around the optimal length.
/// This models real biological trade-offs:
///
/// 1. **Genome Size Constraints**: Organisms face pressure to maintain compact genomes
///    (replication cost) versus having sufficient genes (functional complexity)
/// 2. **Indel Mutations**: Insertions and deletions push sequences away from optimum
/// 3. **Length Polymorphisms**: Real populations show length variation (e.g., STRs)
///    with fitness effects
///
/// The log-normal distribution is biologically motivated:
/// - Asymmetric: sequences much smaller than optimum are more deleterious
///   than sequences much larger (truncated genomes lose genes)
/// - Accounts for both mutation pressure and selection
///
/// Example: If optimum=1000bp and std_dev=0.2, sequences 500-2000bp have reasonable
/// fitness, but 100bp (95% reduction) or 10000bp (10-fold increase) are severely deleterious.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LengthFitness {
    /// Optimal sequence length in bases
    pub optimum: usize,
    /// Standard deviation controlling the width of the fitness curve
    pub std_dev: f64,
}

impl LengthFitness {
    /// Create a new length-based fitness function.
    ///
    /// # Arguments
    /// * `optimum` - Optimal sequence length in bases (must be > 0)
    ///   - Typical values: 1000-10000 for viruses, 100kb-1Mb for bacteria, 1-100 Mb for eukaryotes
    /// * `std_dev` - Standard deviation on log scale (must be > 0.0)
    ///   - Measures tolerance for length deviations: log scale reflects exponential effects
    ///   - 0.1 = very tight length constraint (±10% has significant fitness cost)
    ///   - 0.5 = moderate tolerance (sequences 0.6-1.6× optimum are reasonably fit)
    ///   - 1.0 = wide tolerance (even 0.4-2.5× optimum viable)
    ///
    /// # Errors
    /// Returns an error if parameters are outside valid ranges.
    pub fn new(optimum: usize, std_dev: f64) -> Result<Self, FitnessError> {
        if optimum == 0 {
            return Err(FitnessError::InvalidParameter(
                "optimum must be greater than 0".into(),
            ));
        }
        if std_dev <= 0.0 {
            return Err(FitnessError::InvalidParameter(
                "std_dev must be greater than 0.0".into(),
            ));
        }
        Ok(Self { optimum, std_dev })
    }
}

impl HaplotypeFitness for LengthFitness {
    fn haplotype_fitness(&self, chromosome: &Chromosome) -> FitnessValue {
        if chromosome.is_empty() {
            return FitnessValue::LETHAL_FITNESS;
        }

        let log_length = (chromosome.len() as f64).ln();
        let log_optimum = (self.optimum as f64).ln();
        let deviation = log_length - log_optimum;

        // Log-normal fitness
        let lognorm_fit = -(deviation * deviation) / (2.0 * self.std_dev * self.std_dev);
        FitnessValue::new(lognorm_fit.exp())
    }
}

impl IndividualFitness for LengthFitness {}

/// Sequence similarity fitness function.
///
/// Fitness is higher when two haplotypes have more similar sequences. This models
/// biological scenarios where sequence divergence reduces fitness:
///
/// 1. **Reproductive Incompatibility**: Divergent alleles may not interact properly
/// 2. **Gene Dosage Balance**: Some genes require stoichiometric expression of protein products
/// 3. **Intragenic Complementation**: Within-gene diversity can interfere with splicing or folding
/// 4. **Heterozygote Disadvantage**: Underdominance where heterozygotes are less fit than homozygotes
///
/// The fitness function is based on Hamming distance (counting differences):
/// - Identical sequences: fitness = 1.0 (perfect compatibility)
/// - Completely different sequences: fitness = 0.0 (lethal incompatibility)
/// - Partial difference: fitness between 0 and 1 (proportional to similarity)
/// - Shape parameter controls how sharply fitness declines with divergence
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SequenceSimilarityFitness {
    /// Shape parameter controlling decline rate (higher = steeper decline)
    pub shape: f64,
}

impl SequenceSimilarityFitness {
    /// Create a new sequence similarity fitness function.
    ///
    /// # Arguments
    /// * `shape` - Shape parameter (must be > 0.0)
    ///   - Controls how sharply fitness declines with sequence divergence
    ///   - 0.5 = gradual decline (weak selection against incompatibility)
    ///   - 1.0 = linear decline with distance (moderate selection)
    ///   - 2.0+ = steep decline (strong selection for compatibility)
    ///
    /// # Errors
    /// Returns an error if shape is not positive.
    pub fn new(shape: f64) -> Result<Self, FitnessError> {
        if shape <= 0.0 {
            return Err(FitnessError::InvalidParameter(
                "shape must be greater than 0.0".into(),
            ));
        }
        Ok(Self { shape })
    }

    /// Calculate hamming distance between two sequences.
    fn hamming_distance(seq1: &[Nucleotide], seq2: &[Nucleotide]) -> usize {
        seq1.iter().zip(seq2.iter()).filter(|(a, b)| a != b).count()
    }
}

impl HaplotypeFitness for SequenceSimilarityFitness {
    fn haplotype_fitness(&self, _chromosome: &Chromosome) -> FitnessValue {
        panic!(
            "SequenceSimilarityFitness requires two haplotypes. Use individual_fitness instead."
        );
    }
}

impl IndividualFitness for SequenceSimilarityFitness {
    fn individual_fitness(&self, individual: &Individual) -> FitnessValue {
        let chr1 = individual.haplotype1().get(0).unwrap();
        let chr2 = individual.haplotype2().get(0).unwrap();

        if chr1.is_empty() || chr2.is_empty() {
            return FitnessValue::LETHAL_FITNESS;
        }

        let min_len = chr1.len().min(chr2.len());
        let max_len = chr1.len().max(chr2.len());

        // Compare up to the shorter length (compare actual sequence content)
        // Count differences using Hamming distance
        let distance = Self::hamming_distance(
            &chr1.sequence().as_slice()[..min_len],
            &chr2.sequence().as_slice()[..min_len],
        );

        // Add penalty for length difference
        let length_diff = max_len - min_len;
        let total_diff = distance + length_diff;

        // Normalize by maximum possible length
        let similarity = (max_len - total_diff) as f64 / max_len as f64;

        // Apply shape parameter
        FitnessValue::new(similarity.powf(self.shape))
    }
}

/// Length similarity fitness function.
///
/// Fitness is higher when two haplotypes have more similar lengths.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LengthSimilarityFitness {
    /// Shape parameter controlling decline rate (higher = steeper decline)
    pub shape: f64,
}

impl LengthSimilarityFitness {
    /// Create a new length similarity fitness function.
    ///
    /// # Arguments
    /// * `shape` - Shape parameter (must be > 0.0)
    ///
    /// # Errors
    /// Returns an error if shape is not positive.
    pub fn new(shape: f64) -> Result<Self, FitnessError> {
        if shape <= 0.0 {
            return Err(FitnessError::InvalidParameter(
                "shape must be greater than 0.0".into(),
            ));
        }
        Ok(Self { shape })
    }

    /// Calculate the percent length difference between two sequences.
    /// Identical lengths return 0, completely different lengths return 1.
    fn length_ratio(seq1: &[Nucleotide], seq2: &[Nucleotide]) -> f64 {
        let len1 = seq1.len();
        let len2 = seq2.len();
        let max_len = len1.max(len2);
        if max_len == 0 {
            return 0.0;
        }
        let length_diff = usize::abs_diff(len1, len2);
        length_diff as f64 / max_len as f64
    }
}

impl HaplotypeFitness for LengthSimilarityFitness {
    fn haplotype_fitness(&self, _chromosome: &Chromosome) -> FitnessValue {
        panic!("LengthSimilarityFitness requires two haplotypes. Use individual_fitness instead.");
    }
}

impl IndividualFitness for LengthSimilarityFitness {
    fn individual_fitness(&self, individual: &Individual) -> FitnessValue {
        let chr1 = individual.haplotype1().get(0).unwrap();
        let chr2 = individual.haplotype2().get(0).unwrap();

        if chr1.is_empty() || chr2.is_empty() {
            return FitnessValue::LETHAL_FITNESS;
        }

        // Compare the lengths of the sequences as a ratio
        let similarity =
            1.0 - Self::length_ratio(chr1.sequence().as_slice(), chr2.sequence().as_slice());

        // Apply shape parameter
        FitnessValue::new(similarity.powf(self.shape))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genome::{Haplotype, RepeatMap};
    use std::str::FromStr;

    fn test_chromosome(seq: &str) -> Chromosome {
        let sequence = crate::base::Sequence::from_str(seq).unwrap();
        let total_len = sequence.len();
        let ru_len = if total_len > 0 { total_len } else { 4 };
        let rus_per_hor = 1;
        let num_hors = if total_len > 0 { 1 } else { 0 };
        let map = RepeatMap::uniform(ru_len, rus_per_hor, num_hors);
        Chromosome::new("test", sequence, map)
    }

    // ===== GCContentFitness Tests =====

    #[test]
    fn test_gc_content_fitness_new() {
        let fitness = GCContentFitness::new(0.5, 2.0).unwrap();
        assert_eq!(fitness.optimum, 0.5);
        assert_eq!(fitness.concentration, 2.0);
    }

    #[test]
    fn test_gc_content_fitness_invalid_optimum() {
        assert!(GCContentFitness::new(-0.1, 2.0).is_err());
        assert!(GCContentFitness::new(1.1, 2.0).is_err());
        // 0.0 and 1.0 are now valid
        assert!(GCContentFitness::new(0.0, 2.0).is_ok());
        assert!(GCContentFitness::new(1.0, 2.0).is_ok());
    }

    #[test]
    fn test_gc_content_fitness_invalid_concentration() {
        assert!(GCContentFitness::new(0.5, 0.0).is_err());
        assert!(GCContentFitness::new(0.5, -1.0).is_err());
    }

    #[test]
    fn test_gc_content_fitness_all_gc() {
        let fitness = GCContentFitness::new(0.5, 2.0).unwrap();
        let chr = test_chromosome("GCGCGCGC");
        let score = fitness.haplotype_fitness(&chr);
        assert!(*score > 0.0);
    }

    #[test]
    fn test_gc_content_fitness_all_at() {
        let fitness = GCContentFitness::new(0.5, 2.0).unwrap();
        let chr = test_chromosome("ATATATAT");
        let score = fitness.haplotype_fitness(&chr);
        assert!(*score > 0.0);
    }

    #[test]
    fn test_gc_content_fitness_half() {
        let fitness = GCContentFitness::new(0.5, 2.0).unwrap();
        let chr = test_chromosome("ACGT");
        let score = fitness.haplotype_fitness(&chr);
        assert!(*score > 0.0);
    }

    // ===== LengthFitness Tests =====

    #[test]
    fn test_length_fitness_new() {
        let fitness = LengthFitness::new(100, 0.5).unwrap();
        assert_eq!(fitness.optimum, 100);
        assert_eq!(fitness.std_dev, 0.5);
    }

    #[test]
    fn test_length_fitness_invalid_optimum() {
        assert!(LengthFitness::new(0, 0.5).is_err());
    }

    #[test]
    fn test_length_fitness_invalid_std_dev() {
        assert!(LengthFitness::new(100, 0.0).is_err());
        assert!(LengthFitness::new(100, -1.0).is_err());
    }

    #[test]
    fn test_length_fitness_at_optimum() {
        let fitness = LengthFitness::new(8, 0.5).unwrap();
        let chr = test_chromosome("ACGTACGT");
        let score = fitness.haplotype_fitness(&chr);
        assert!((*score - 1.0).abs() < 0.01); // Should be close to 1.0
    }

    #[test]
    fn test_length_fitness_decreases_with_distance() {
        let fitness = LengthFitness::new(8, 0.5).unwrap();
        let chr_optimal = test_chromosome("ACGTACGT");
        let chr_longer = test_chromosome("ACGTACGTACGTACGT");

        let score_optimal = fitness.haplotype_fitness(&chr_optimal);
        let score_longer = fitness.haplotype_fitness(&chr_longer);

        assert!(*score_optimal > *score_longer);
    }

    // ===== SequenceSimilarityFitness Tests =====

    #[test]
    fn test_sequence_similarity_fitness_new() {
        let fitness = SequenceSimilarityFitness::new(1.0).unwrap();
        assert_eq!(fitness.shape, 1.0);
    }

    #[test]
    fn test_sequence_similarity_fitness_invalid_shape() {
        assert!(SequenceSimilarityFitness::new(0.0).is_err());
        assert!(SequenceSimilarityFitness::new(-1.0).is_err());
    }

    #[test]
    #[should_panic(expected = "requires two haplotypes")]
    fn test_sequence_similarity_panic_on_single_haplotype() {
        let fitness = SequenceSimilarityFitness::new(1.0).unwrap();
        let chr = test_chromosome("ACGT");
        fitness.haplotype_fitness(&chr);
    }

    #[test]
    fn test_sequence_similarity_identical() {
        let fitness = SequenceSimilarityFitness::new(1.0).unwrap();

        let chr = test_chromosome("ACGTACGT");
        let mut hap1 = Haplotype::new();
        let mut hap2 = Haplotype::new();
        hap1.push(chr.clone());
        hap2.push(chr);

        let ind = Individual::new("test", hap1, hap2);
        let score = fitness.individual_fitness(&ind);

        assert_eq!(*score, 1.0); // Identical sequences
    }

    #[test]
    fn test_sequence_similarity_completely_different() {
        let fitness = SequenceSimilarityFitness::new(1.0).unwrap();

        let chr1 = test_chromosome("AAAAAAAA");
        let chr2 = test_chromosome("TTTTTTTT");
        let mut hap1 = Haplotype::new();
        let mut hap2 = Haplotype::new();
        hap1.push(chr1);
        hap2.push(chr2);

        let ind = Individual::new("test", hap1, hap2);
        let score = fitness.individual_fitness(&ind);

        assert_eq!(*score, 0.0); // Completely different
    }

    #[test]
    fn test_sequence_similarity_partial() {
        let fitness = SequenceSimilarityFitness::new(1.0).unwrap();

        let chr1 = test_chromosome("ACGTACGT");
        let chr2 = test_chromosome("ACGTTCGT");
        let mut hap1 = Haplotype::new();
        let mut hap2 = Haplotype::new();
        hap1.push(chr1);
        hap2.push(chr2);

        let ind = Individual::new("test", hap1, hap2);
        let score = fitness.individual_fitness(&ind);

        assert!(*score > 0.0 && *score < 1.0); // Partially similar
    }

    #[test]
    fn test_hamming_distance() {
        let seq1 = vec![Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T];
        let seq2 = vec![Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T];
        assert_eq!(SequenceSimilarityFitness::hamming_distance(&seq1, &seq2), 0);

        let seq3 = vec![Nucleotide::A, Nucleotide::A, Nucleotide::A, Nucleotide::A];
        let seq4 = vec![Nucleotide::C, Nucleotide::C, Nucleotide::C, Nucleotide::C];
        assert_eq!(SequenceSimilarityFitness::hamming_distance(&seq3, &seq4), 4);

        let seq5 = vec![Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T];
        let seq6 = vec![Nucleotide::A, Nucleotide::A, Nucleotide::G, Nucleotide::G];
        assert_eq!(SequenceSimilarityFitness::hamming_distance(&seq5, &seq6), 2);
    }

    #[test]
    fn test_length_similarity_fitness_new() {
        let fitness = LengthSimilarityFitness::new(1.0).unwrap();
        assert_eq!(fitness.shape, 1.0);
    }

    #[test]
    fn test_length_similarity_fitness_invalid_shape() {
        assert!(LengthSimilarityFitness::new(0.0).is_err());
        assert!(LengthSimilarityFitness::new(-1.0).is_err());
    }

    #[test]
    #[should_panic(expected = "requires two haplotypes")]
    fn test_length_similarity_panic_on_single_haplotype() {
        let fitness = LengthSimilarityFitness::new(1.0).unwrap();
        let chr = test_chromosome("ACGT");
        fitness.haplotype_fitness(&chr);
    }

    #[test]
    fn test_length_similarity_identical() {
        let fitness = LengthSimilarityFitness::new(1.0).unwrap();

        let chr = test_chromosome("ACGTACGT");
        let mut hap1 = Haplotype::new();
        let mut hap2 = Haplotype::new();
        hap1.push(chr.clone());
        hap2.push(chr);

        let ind = Individual::new("test", hap1, hap2);
        let score = fitness.individual_fitness(&ind);

        assert_eq!(*score, 1.0); // Identical lengths
    }

    #[test]
    fn test_length_similarity_different_lengths() {
        let fitness = LengthSimilarityFitness::new(1.0).unwrap();

        let chr1 = test_chromosome("ACGT");
        let chr2 = test_chromosome("ACGTACGT");
        let mut hap1 = Haplotype::new();
        let mut hap2 = Haplotype::new();
        hap1.push(chr1);
        hap2.push(chr2);

        let ind = Individual::new("test", hap1, hap2);
        let score = fitness.individual_fitness(&ind);

        // Length diff is 4, max length is 8, so similarity = 1 - 4/8 = 0.5
        assert_eq!(*score, 0.5);
    }

    #[test]
    fn test_length_similarity_completely_different() {
        let fitness = LengthSimilarityFitness::new(1.0).unwrap();

        let chr1 = test_chromosome("A");
        let chr2 = test_chromosome("ACGTACGTACGTACGT");
        let mut hap1 = Haplotype::new();
        let mut hap2 = Haplotype::new();
        hap1.push(chr1);
        hap2.push(chr2);

        let ind = Individual::new("test", hap1, hap2);
        let score = fitness.individual_fitness(&ind);

        // Length diff is 15, max length is 16, so similarity = 1 - 15/16 = 0.0625
        assert!((*score - 0.0625).abs() < 0.001);
    }

    #[test]
    fn test_length_similarity_with_shape() {
        let fitness = LengthSimilarityFitness::new(2.0).unwrap();

        let chr1 = test_chromosome("ACGT");
        let chr2 = test_chromosome("ACGTACGT");
        let mut hap1 = Haplotype::new();
        let mut hap2 = Haplotype::new();
        hap1.push(chr1);
        hap2.push(chr2);

        let ind = Individual::new("test", hap1, hap2);
        let score = fitness.individual_fitness(&ind);

        // Similarity is 0.5, with shape 2.0: 0.5^2.0 = 0.25
        assert_eq!(*score, 0.25);
    }

    #[test]
    fn test_length_ratio() {
        let seq1 = vec![Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T];
        let seq2 = vec![Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T];
        assert_eq!(LengthSimilarityFitness::length_ratio(&seq1, &seq2), 0.0);

        let seq3 = vec![Nucleotide::A, Nucleotide::C];
        let seq4 = vec![Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T];
        assert_eq!(LengthSimilarityFitness::length_ratio(&seq3, &seq4), 0.5);

        let seq5 = vec![Nucleotide::A];
        let seq6 = vec![Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T];
        assert_eq!(LengthSimilarityFitness::length_ratio(&seq5, &seq6), 0.75);

        let empty1: Vec<Nucleotide> = vec![];
        let empty2: Vec<Nucleotide> = vec![];
        assert_eq!(LengthSimilarityFitness::length_ratio(&empty1, &empty2), 0.0);
    }

    #[test]
    fn test_fitness_error_display() {
        let err = FitnessError::InvalidParameter("test message".into());
        let msg = format!("{err}");
        assert!(msg.contains("Invalid fitness parameter"));
        assert!(msg.contains("test message"));
    }

    #[test]
    fn test_gc_content_fitness_clone() {
        let fitness1 = GCContentFitness::new(0.5, 2.0).unwrap();
        let fitness2 = fitness1.clone();

        assert_eq!(fitness1.optimum, fitness2.optimum);
        assert_eq!(fitness1.concentration, fitness2.concentration);
    }
}
