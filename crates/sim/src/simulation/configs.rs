//! Simulation parameters and configuration.
//!
//! This module provides parameter structures for configuring simulations,
//! including mutation rates, recombination rates, fitness parameters, and
//! simulation settings.

use crate::base::Nucleotide;
use crate::errors::FitnessError;
use crate::evolution::{
    GCContentFitness, IndelModel, LengthFitness, LengthSimilarityFitness, RecombinationModel,
    SequenceSimilarityFitness, SubstitutionModel,
};
use centrevo_codec::CodecStrategy;
use serde::{Deserialize, Serialize};
use std::marker::PhantomData;

/// The master configuration struct.
/// Can be deserialized from a file to fully reproduce a simulation setup.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Configuration {
    pub execution: ExecutionConfig,
    pub evolution: EvolutionConfig,
    pub initialization: InitializationConfig,
}

/// High-level simulation parameters.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExecutionConfig {
    /// Number of diploid individuals in population
    pub population_size: usize,
    /// Total number of generations to simulate
    pub total_generations: usize,
    /// Optional RNG seed for reproducibility
    pub seed: Option<u64>,
    /// Encoding strategy for database storage
    #[serde(default = "default_codec_compat")]
    pub codec: CodecStrategy,
}

/// Grouped evolutionary parameters.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EvolutionConfig {
    pub mutation: MutationConfig,
    pub recombination: RecombinationConfig,
    pub fitness: FitnessConfig,
}

/// Configuration for sequence initialization.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum InitializationConfig {
    /// Generate sequences based on a repeat structure.
    Generate {
        /// Repeat structure parameters
        structure: UniformRepeatStructure,
        /// Generation mode (Uniform or Random)
        mode: GenerationMode,
    },
    /// Load sequences from an external source.
    Load {
        /// Source of sequences
        source: SequenceSource,
    },
}

/// Parameters for defining a uniform repeat sequence structure.
///
/// This structure is used for configuration and reproducibility of simulations
/// that start with a uniform repeat structure.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UniformRepeatStructure {
    /// Initial nucleotide to fill sequences
    pub init_base: Nucleotide,
    /// Repeat unit length in bases
    pub ru_length: usize,
    /// Number of repeat units per higher-order repeat
    pub rus_per_hor: usize,
    /// Number of HORs per chromosome
    pub hors_per_chr: usize,
    /// Number of chromosomes per haplotype
    pub chrs_per_hap: usize,
}

impl UniformRepeatStructure {
    /// Create a new uniform repeat structure configuration.
    pub fn new(
        init_base: Nucleotide,
        ru_length: usize,
        rus_per_hor: usize,
        hors_per_chr: usize,
        chrs_per_hap: usize,
    ) -> Self {
        Self {
            init_base,
            ru_length,
            rus_per_hor,
            hors_per_chr,
            chrs_per_hap,
        }
    }

    /// Get total sequence length for a single chromosome.
    pub fn chr_length(&self) -> usize {
        self.ru_length * self.rus_per_hor * self.hors_per_chr
    }
}

/// Mode for sequence generation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum GenerationMode {
    /// Uniform sequences (all repeats identical)
    Uniform,
    /// Random sequences (random bases)
    Random,
}

/// Source for loading sequences.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum SequenceSource {
    /// FASTA file path
    Fasta {
        /// Path to FASTA file
        path: String,
        /// Path to BED file for structure
        bed_path: String,
    },
    /// Simulation database
    Database {
        /// Path to database
        path: String,
        /// Simulation ID
        sim_id: String,
        /// Generation to load (None for last)
        generation: Option<usize>,
    },
    /// Formatted string
    FormattedString {
        /// Sequence string
        sequence: String,
        /// HOR delimiter
        hor_delim: char,
        /// RU delimiter
        ru_delim: char,
    },
}

/// Parameters for mutation processes.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MutationConfig {
    /// Substitution model for point mutations
    #[serde(rename = "substitution")]
    pub substitution: SubstitutionModel,
    /// Indel model for insertions and deletions
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub indel: Option<IndelModel>,
}

impl MutationConfig {
    /// Create new mutation configuration.
    pub fn new(substitution: SubstitutionModel) -> Self {
        Self {
            substitution,
            indel: None,
        }
    }

    /// Create new mutation configuration with indels.
    pub fn with_indels(substitution: SubstitutionModel, indel: IndelModel) -> Self {
        Self {
            substitution,
            indel: Some(indel),
        }
    }

    /// Create with uniform mutation rate.
    pub fn uniform(rate: f64) -> Result<Self, String> {
        let substitution = SubstitutionModel::uniform(rate).map_err(|e| format!("{e}"))?;
        Ok(Self {
            substitution,
            indel: None,
        })
    }
}

/// Parameters for recombination processes.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RecombinationConfig {
    /// Recombination parameters
    pub params: RecombinationModel,
}

impl RecombinationConfig {
    /// Create new recombination configuration.
    pub fn new(params: RecombinationModel) -> Self {
        Self { params }
    }

    /// Create with standard parameters.
    pub fn standard(
        break_prob: f64,
        crossover_prob: f64,
        gc_extension_prob: f64,
    ) -> Result<Self, String> {
        // Default homology parameters: no homology preference (random), no window, k=1 (ignored)
        let params = RecombinationModel::builder()
            .break_prob(break_prob)
            .crossover_prob(crossover_prob)
            .gc_extension_prob(gc_extension_prob)
            .homology_strength(0.0)
            .search_window(0)
            .kmer_size(1)
            .build()
            .map_err(|e| format!("{e}"))?;
        Ok(Self { params })
    }
}

/// Parameters for fitness/selection.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct FitnessConfig {
    /// GC content fitness function
    pub gc_content: Option<GCContentFitness>,
    /// Length-based fitness function
    pub length: Option<LengthFitness>,
    /// Sequence similarity fitness function
    pub seq_similarity: Option<SequenceSimilarityFitness>,
    /// Length similarity fitness function
    pub length_similarity: Option<LengthSimilarityFitness>,
}

impl FitnessConfig {
    /// Create new fitness configuration.
    pub fn new(
        gc_content: Option<GCContentFitness>,
        length: Option<LengthFitness>,
        seq_similarity: Option<SequenceSimilarityFitness>,
        length_similarity: Option<LengthSimilarityFitness>,
    ) -> Self {
        Self {
            gc_content,
            length,
            seq_similarity,
            length_similarity,
        }
    }

    /// Create neutral fitness (no selection).
    pub fn neutral() -> Self {
        Self {
            gc_content: None,
            length: None,
            seq_similarity: None,
            length_similarity: None,
        }
    }

    /// Check if all fitness components are None (neutral).
    pub fn is_neutral(&self) -> bool {
        self.gc_content.is_none()
            && self.length.is_none()
            && self.seq_similarity.is_none()
            && self.length_similarity.is_none()
    }

    /// Compute the total fitness of an individual by multiplying all configured fitness components.
    ///
    /// This method applies all enabled fitness functions (GC content, length, sequence similarity,
    /// length similarity) and returns their product.
    pub fn compute_fitness(
        &self,
        individual: &crate::genome::Individual,
        arena: &crate::base::GenomeArena,
    ) -> crate::base::FitnessValue {
        use crate::evolution::IndividualFitness;

        let mut fitness = crate::base::FitnessValue::default();

        if let Some(gc) = &self.gc_content {
            fitness *= gc.individual_fitness(individual, arena);
        }
        if let Some(len_fit) = &self.length {
            fitness *= len_fit.individual_fitness(individual, arena);
        }
        if let Some(sim) = &self.seq_similarity {
            fitness *= sim.individual_fitness(individual, arena);
        }
        if let Some(len_sim) = &self.length_similarity {
            fitness *= len_sim.individual_fitness(individual, arena);
        }

        fitness
    }

    /// Update cached fitness values for an individual and its haplotypes.
    ///
    /// This computes per-haplotype fitness (for compatible models like GC and Length)
    /// and total individual fitness, caching all values in the provided individual.
    pub fn update_cached_fitness(
        &self,
        individual: &mut crate::genome::Individual,
        arena: &crate::base::GenomeArena,
    ) {
        use crate::evolution::HaplotypeFitness;

        // 1. Compute and cache per-haplotype fitness components
        // We only compute for models that implement HaplotypeFitness (GC and Length)
        let mut h1_fitness = crate::base::FitnessValue::default();
        let mut h2_fitness = crate::base::FitnessValue::default();

        if let Some(gc) = &self.gc_content {
            h1_fitness *= gc.haplotype_fitness(
                individual.haplotype1().get(0).expect("Chromosome missing"),
                arena,
            );
            h2_fitness *= gc.haplotype_fitness(
                individual.haplotype2().get(0).expect("Chromosome missing"),
                arena,
            );
        }
        if let Some(len_fit) = &self.length {
            h1_fitness *= len_fit.haplotype_fitness(
                individual.haplotype1().get(0).expect("Chromosome missing"),
                arena,
            );
            h2_fitness *= len_fit.haplotype_fitness(
                individual.haplotype2().get(0).expect("Chromosome missing"),
                arena,
            );
        }

        individual.haplotype1_mut().set_cached_fitness(h1_fitness);
        individual.haplotype2_mut().set_cached_fitness(h2_fitness);

        // 2. Compute and cache total individual fitness
        let total_fitness = self.compute_fitness(individual, arena);
        individual.set_cached_fitness(total_fitness);
    }
}

/// Type state for an empty builder (no fitness components added yet).
#[derive(Debug, Clone, Copy)]
pub struct BuilderEmpty;

/// Type state for a builder with at least one fitness component.
#[derive(Debug, Clone, Copy)]
pub struct BuilderInitialized;

/// Builder for constructing fitness configurations.
///
/// Uses type states to ensure compile-time safety:
/// - `BuilderEmpty` state: Can call `neutral()` or start with `with_*` methods
/// - `BuilderWithSelection` state: Can chain `with_*` methods and call `build()`
///
/// # Examples
///
/// ```
/// use centrevo_sim::simulation::{FitnessConfigBuilder, BuilderEmpty};
///
/// // Neutral fitness (no selection)
/// let fitness = FitnessConfigBuilder::neutral();
///
/// // Single fitness component
/// let fitness = FitnessConfigBuilder::<BuilderEmpty>::with_gc_content(0.5, 2.0)
///     .unwrap()
///     .build();
///
/// // Multiple fitness components
/// let fitness = FitnessConfigBuilder::<BuilderEmpty>::with_gc_content(0.5, 2.0)
///     .unwrap()
///     .with_length(20000, 0.5)
///     .unwrap()
///     .with_similarity(2.0)
///     .unwrap()
///     .build();
/// ```
#[derive(Debug, Clone)]
pub struct FitnessConfigBuilder<State = BuilderEmpty> {
    gc_content: Option<GCContentFitness>,
    length: Option<LengthFitness>,
    seq_similarity: Option<SequenceSimilarityFitness>,
    length_similarity: Option<LengthSimilarityFitness>,
    _state: PhantomData<State>,
}

impl FitnessConfigBuilder<BuilderEmpty> {
    /// Create neutral fitness configuration (no selection).
    ///
    /// Returns a `FitnessConfig` directly, not a builder.
    pub fn neutral() -> FitnessConfig {
        FitnessConfig::neutral()
    }

    /// Start building with GC content fitness.
    ///
    /// # Arguments
    /// * `optimum` - Optimal GC content (0.0 to 1.0)
    /// * `concentration` - Concentration parameter (> 0.0)
    pub fn with_gc_content(
        optimum: f64,
        concentration: f64,
    ) -> Result<FitnessConfigBuilder<BuilderInitialized>, FitnessError> {
        Ok(FitnessConfigBuilder {
            gc_content: Some(GCContentFitness::new(optimum, concentration)?),
            length: None,
            seq_similarity: None,
            length_similarity: None,
            _state: PhantomData,
        })
    }

    /// Start building with length-based fitness.
    ///
    /// # Arguments
    /// * `optimum` - Optimal sequence length in bases (> 0)
    /// * `std_dev` - Standard deviation (> 0.0)
    pub fn with_length(
        optimum: usize,
        std_dev: f64,
    ) -> Result<FitnessConfigBuilder<BuilderInitialized>, FitnessError> {
        Ok(FitnessConfigBuilder {
            gc_content: None,
            length: Some(LengthFitness::new(optimum, std_dev)?),
            seq_similarity: None,
            length_similarity: None,
            _state: PhantomData,
        })
    }

    /// Start building with sequence similarity fitness.
    ///
    /// # Arguments
    /// * `shape` - Shape parameter (> 0.0)
    pub fn with_similarity(
        shape: f64,
    ) -> Result<FitnessConfigBuilder<BuilderInitialized>, FitnessError> {
        Ok(FitnessConfigBuilder {
            gc_content: None,
            length: None,
            seq_similarity: Some(SequenceSimilarityFitness::new(shape)?),
            length_similarity: None,
            _state: PhantomData,
        })
    }

    /// Start building with length similarity fitness.
    ////
    /// # Arguments
    /// * `shape` - Shape parameter (> 0.0)
    pub fn with_length_similarity(
        shape: f64,
    ) -> Result<FitnessConfigBuilder<BuilderInitialized>, FitnessError> {
        Ok(FitnessConfigBuilder {
            gc_content: None,
            length: None,
            seq_similarity: None,
            length_similarity: Some(LengthSimilarityFitness::new(shape)?),
            _state: PhantomData,
        })
    }
}

impl FitnessConfigBuilder<BuilderInitialized> {
    //TODO: Change name prefix to "and_" to indicate chaining and avoid multiple implementations confusion.
    /// Add GC content fitness to the configuration.
    ///
    /// # Arguments
    /// * `optimum` - Optimal GC content (0.0 to 1.0)
    /// * `concentration` - Concentration parameter (> 0.0)
    ///
    /// # Errors
    /// Returns an error if GC content fitness is already set.
    pub fn with_gc_content(
        mut self,
        optimum: f64,
        concentration: f64,
    ) -> Result<Self, FitnessError> {
        if self.gc_content.is_some() {
            return Err(FitnessError::InvalidParameter(
                "GC content fitness already set".into(),
            ));
        }
        self.gc_content = Some(GCContentFitness::new(optimum, concentration)?);
        Ok(self)
    }

    /// Add length-based fitness to the configuration.
    ///
    /// # Arguments
    /// * `optimum` - Optimal sequence length in bases (> 0)
    /// * `std_dev` - Standard deviation (> 0.0)
    ///
    /// # Errors
    /// Returns an error if length fitness is already set.
    pub fn with_length(mut self, optimum: usize, std_dev: f64) -> Result<Self, FitnessError> {
        if self.length.is_some() {
            return Err(FitnessError::InvalidParameter(
                "Length fitness already set".into(),
            ));
        }
        self.length = Some(LengthFitness::new(optimum, std_dev)?);
        Ok(self)
    }

    /// Add sequence similarity fitness to the configuration.
    ///
    /// # Arguments
    /// * `shape` - Shape parameter (> 0.0)
    ///
    /// # Errors
    /// Returns an error if similarity fitness is already set.
    pub fn with_similarity(mut self, shape: f64) -> Result<Self, FitnessError> {
        if self.seq_similarity.is_some() {
            return Err(FitnessError::InvalidParameter(
                "Similarity fitness already set".into(),
            ));
        }
        self.seq_similarity = Some(SequenceSimilarityFitness::new(shape)?);
        Ok(self)
    }

    /// Add length similarity fitness to the configuration.
    ///
    /// # Arguments
    /// * `shape` - Shape parameter (> 0.0)
    ///
    /// # Errors
    /// Returns an error if length similarity fitness is already set.
    pub fn with_length_similarity(mut self, shape: f64) -> Result<Self, FitnessError> {
        if self.length_similarity.is_some() {
            return Err(FitnessError::InvalidParameter(
                "Length similarity fitness already set".into(),
            ));
        }
        self.length_similarity = Some(LengthSimilarityFitness::new(shape)?);
        Ok(self)
    }

    /// Build the final fitness configuration.
    pub fn build(self) -> FitnessConfig {
        FitnessConfig::new(
            self.gc_content,
            self.length,
            self.seq_similarity,
            self.length_similarity,
        )
    }
}

/// Default codec for backward compatibility (loading old configs).
/// Old configurations used BitPackedRS implicitly.
fn default_codec_compat() -> CodecStrategy {
    CodecStrategy::BitPackedRS
}

impl ExecutionConfig {
    /// Create new simulation configuration.
    pub fn new(population_size: usize, total_generations: usize, seed: Option<u64>) -> Self {
        Self {
            population_size,
            total_generations,
            seed,
            codec: CodecStrategy::default(),
        }
    }

    /// Set the codec strategy.
    pub fn with_codec(mut self, codec: CodecStrategy) -> Self {
        self.codec = codec;
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simulation_config_codec() {
        let config = ExecutionConfig::new(100, 100, None);
        // Default should be UnpackedZRS
        assert!(matches!(config.codec, CodecStrategy::UnpackedZRS));

        let config = config.with_codec(CodecStrategy::UnpackedRS);
        assert!(matches!(config.codec, CodecStrategy::UnpackedRS));
    }

    #[test]
    fn test_repeat_structure_new() {
        let structure = UniformRepeatStructure::new(Nucleotide::A, 171, 12, 100, 1);

        assert_eq!(structure.ru_length, 171);
        assert_eq!(structure.rus_per_hor, 12);
        assert_eq!(structure.hors_per_chr, 100);
        assert_eq!(structure.chrs_per_hap, 1);
    }

    #[test]
    fn test_repeat_structure_chr_length() {
        let structure = UniformRepeatStructure::new(Nucleotide::A, 10, 5, 20, 1);

        // 10 * 5 * 20 = 1000
        assert_eq!(structure.chr_length(), 1000);
    }

    #[test]
    fn test_mutation_config_uniform() {
        let config = MutationConfig::uniform(0.001).unwrap();
        // Just check it was created successfully
        assert!(matches!(config.substitution, SubstitutionModel::Uniform(_)));
        assert!(config.indel.is_none());
    }

    #[test]
    fn test_recombination_config_standard() {
        let config = RecombinationConfig::standard(0.01, 0.7, 0.1).unwrap();
        assert_eq!(config.params.break_prob(), 0.01);
        assert_eq!(config.params.crossover_prob(), 0.7);
    }

    #[test]
    fn test_fitness_config_neutral() {
        let config = FitnessConfig::neutral();
        assert!(config.is_neutral());
        assert!(config.gc_content.is_none());
        assert!(config.length.is_none());
        assert!(config.seq_similarity.is_none());
    }

    #[test]
    fn test_fitness_config_with_components() {
        let gc = GCContentFitness::new(0.5, 2.0).unwrap();
        let config = FitnessConfig::new(Some(gc), None, None, None);

        assert!(!config.is_neutral());
        assert!(config.gc_content.is_some());
        assert!(config.length.is_none());
    }

    #[test]
    fn test_fitness_builder_neutral() {
        let config = FitnessConfigBuilder::<BuilderEmpty>::neutral();
        assert!(config.is_neutral());
    }

    #[test]
    fn test_fitness_builder_single_gc_content() {
        let config = FitnessConfigBuilder::<BuilderEmpty>::with_gc_content(0.5, 2.0)
            .unwrap()
            .build();

        assert!(!config.is_neutral());
        assert!(config.gc_content.is_some());
        assert!(config.length.is_none());
        assert!(config.seq_similarity.is_none());
    }

    #[test]
    fn test_fitness_builder_single_length() {
        let config = FitnessConfigBuilder::<BuilderEmpty>::with_length(20000, 0.5)
            .unwrap()
            .build();

        assert!(!config.is_neutral());
        assert!(config.gc_content.is_none());
        assert!(config.length.is_some());
        assert!(config.seq_similarity.is_none());
    }

    #[test]
    fn test_fitness_builder_multiple_components() {
        let config = FitnessConfigBuilder::<BuilderEmpty>::with_gc_content(0.5, 2.0)
            .unwrap()
            .with_length(20000, 0.5)
            .unwrap()
            .with_similarity(2.0)
            .unwrap()
            .build();

        assert!(!config.is_neutral());
        assert!(config.gc_content.is_some());
        assert!(config.length.is_some());
        assert!(config.seq_similarity.is_some());
    }

    #[test]
    fn test_fitness_builder_duplicate_component_error() {
        let result = FitnessConfigBuilder::<BuilderEmpty>::with_gc_content(0.5, 2.0)
            .unwrap()
            .with_gc_content(0.6, 3.0);

        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("already set"));
    }

    #[test]
    fn test_fitness_builder_invalid_parameters() {
        // Invalid GC content optimum
        assert!(FitnessConfigBuilder::<BuilderEmpty>::with_gc_content(1.5, 2.0).is_err());

        // Invalid concentration
        assert!(FitnessConfigBuilder::<BuilderEmpty>::with_gc_content(0.5, -1.0).is_err());

        // Invalid length optimum
        assert!(FitnessConfigBuilder::<BuilderEmpty>::with_length(0, 0.5).is_err());

        // Invalid std_dev
        assert!(FitnessConfigBuilder::<BuilderEmpty>::with_length(20000, 0.0).is_err());

        // Invalid similarity shape
        assert!(FitnessConfigBuilder::<BuilderEmpty>::with_similarity(0.0).is_err());
    }

    #[test]
    fn test_fitness_builder_order_independence() {
        // Build with different orders
        let config1 = FitnessConfigBuilder::<BuilderEmpty>::with_gc_content(0.5, 2.0)
            .unwrap()
            .with_length(20000, 0.5)
            .unwrap()
            .build();

        let config2 = FitnessConfigBuilder::<BuilderEmpty>::with_length(20000, 0.5)
            .unwrap()
            .with_gc_content(0.5, 2.0)
            .unwrap()
            .build();

        // Both should have the same components
        assert!(config1.gc_content.is_some());
        assert!(config1.length.is_some());
        assert!(config2.gc_content.is_some());
        assert!(config2.length.is_some());
    }

    #[test]
    fn test_fitness_builder_single_length_similarity() {
        let config = FitnessConfigBuilder::<BuilderEmpty>::with_length_similarity(2.0)
            .unwrap()
            .build();

        assert!(!config.is_neutral());
        assert!(config.gc_content.is_none());
        assert!(config.length.is_none());
        assert!(config.seq_similarity.is_none());
        assert!(config.length_similarity.is_some());
    }

    #[test]
    fn test_fitness_builder_with_length_similarity_chained() {
        let config = FitnessConfigBuilder::<BuilderEmpty>::with_gc_content(0.5, 2.0)
            .unwrap()
            .with_length_similarity(3.0)
            .unwrap()
            .build();

        assert!(!config.is_neutral());
        assert!(config.gc_content.is_some());
        assert!(config.length_similarity.is_some());
    }

    #[test]
    fn test_fitness_builder_all_components_including_length_similarity() {
        let config = FitnessConfigBuilder::<BuilderEmpty>::with_gc_content(0.5, 2.0)
            .unwrap()
            .with_length(20000, 0.5)
            .unwrap()
            .with_similarity(2.0)
            .unwrap()
            .with_length_similarity(1.5)
            .unwrap()
            .build();

        assert!(!config.is_neutral());
        assert!(config.gc_content.is_some());
        assert!(config.length.is_some());
        assert!(config.seq_similarity.is_some());
        assert!(config.length_similarity.is_some());
    }

    #[test]
    fn test_fitness_builder_duplicate_length_similarity_error() {
        let result = FitnessConfigBuilder::<BuilderEmpty>::with_length_similarity(2.0)
            .unwrap()
            .with_length_similarity(3.0);

        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("already set"));
    }

    #[test]
    fn test_fitness_builder_invalid_length_similarity_shape() {
        // Invalid shape (zero)
        assert!(FitnessConfigBuilder::<BuilderEmpty>::with_length_similarity(0.0).is_err());

        // Invalid shape (negative)
        assert!(FitnessConfigBuilder::<BuilderEmpty>::with_length_similarity(-1.0).is_err());
    }

    #[test]
    fn test_simulation_config_new() {
        let config = ExecutionConfig::new(100, 1000, Some(42));

        assert_eq!(config.population_size, 100);
        assert_eq!(config.total_generations, 1000);
        assert_eq!(config.seed, Some(42));
    }

    #[test]
    fn test_fitness_config_update_cached_fitness() {
        use crate::base::GenomeArena;
        use crate::base::Nucleotide;
        use crate::genome::{Chromosome, Haplotype, Individual, RepeatMap};

        let gc_fitness = GCContentFitness::new(0.5, 10.0).unwrap();
        let config = FitnessConfig::new(Some(gc_fitness), None, None, None);

        let seq1 = crate::base::Sequence::from_nucleotides(vec![
            Nucleotide::G,
            Nucleotide::C,
            Nucleotide::G,
            Nucleotide::C,
        ]);
        let seq2 = crate::base::Sequence::from_nucleotides(vec![
            Nucleotide::A,
            Nucleotide::T,
            Nucleotide::A,
            Nucleotide::T,
        ]);

        // GC content: seq1=1.0, seq2=0.0. Optimum 0.5.
        // Haplotype fitness will be non-zero (due to clamping in GCContentFitness)

        // Create arena and allocate sequences
        let mut arena = GenomeArena::new();
        let slice1 = arena.alloc(seq1.as_slice());
        let slice2 = arena.alloc(seq2.as_slice());

        let chr1 = Chromosome::new("c1", slice1, RepeatMap::uniform(4, 1, 1));
        let chr2 = Chromosome::new("c2", slice2, RepeatMap::uniform(4, 1, 1));

        let h1 = Haplotype::from_chromosomes(vec![chr1]);
        let h2 = Haplotype::from_chromosomes(vec![chr2]);

        let mut ind = Individual::new("ind1", h1, h2);

        assert!(ind.cached_fitness().is_none());
        assert!(ind.haplotype1().cached_fitness().is_none());
        assert!(ind.haplotype2().cached_fitness().is_none());

        config.update_cached_fitness(&mut ind, &arena);

        assert!(ind.cached_fitness().is_some());
        assert!(ind.haplotype1().cached_fitness().is_some());
        assert!(ind.haplotype2().cached_fitness().is_some());

        // Haplotype 1 (GC=1.0) and Haplotype 2 (GC=0.0) should have same fitness with optimum 0.5 and symmetric Beta
        assert!(
            (*ind.haplotype1().cached_fitness().unwrap()
                - *ind.haplotype2().cached_fitness().unwrap())
            .abs()
                < 1e-10
        );
    }
}
