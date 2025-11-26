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
use serde::{Deserialize, Serialize};
use std::marker::PhantomData;

/// Parameters for defining the initial repeat sequence structure.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RepeatStructure {
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

impl RepeatStructure {
    /// Create a new repeat structure configuration.
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

/// High-level simulation parameters.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationConfig {
    /// Number of diploid individuals in population
    pub population_size: usize,
    /// Total number of generations to simulate
    pub total_generations: usize,
    /// Optional RNG seed for reproducibility
    pub seed: Option<u64>,
}

impl SimulationConfig {
    /// Create new simulation configuration.
    pub fn new(population_size: usize, total_generations: usize, seed: Option<u64>) -> Self {
        Self {
            population_size,
            total_generations,
            seed,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_repeat_structure_new() {
        let structure = RepeatStructure::new(Nucleotide::A, 171, 12, 100, 1);

        assert_eq!(structure.ru_length, 171);
        assert_eq!(structure.rus_per_hor, 12);
        assert_eq!(structure.hors_per_chr, 100);
        assert_eq!(structure.chrs_per_hap, 1);
    }

    #[test]
    fn test_repeat_structure_chr_length() {
        let structure = RepeatStructure::new(Nucleotide::A, 10, 5, 20, 1);

        // 10 * 5 * 20 = 1000
        assert_eq!(structure.chr_length(), 1000);
    }

    #[test]
    fn test_mutation_config_uniform() {
        let config = MutationConfig::uniform(0.001).unwrap();
        // Just check it was created successfully
        assert!(matches!(config.substitution, SubstitutionModel { .. }));
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
        let config = SimulationConfig::new(100, 1000, Some(42));

        assert_eq!(config.population_size, 100);
        assert_eq!(config.total_generations, 1000);
        assert_eq!(config.seed, Some(42));
    }
}
