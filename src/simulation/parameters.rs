//! Simulation parameters and configuration.
//!
//! This module provides parameter structures for configuring simulations,
//! including mutation rates, recombination rates, fitness parameters, and
//! simulation settings.

use crate::evolution::{SubstitutionModel, RecombinationParams, GCContentFitness, LengthFitness, SequenceSimilarityFitness};
use crate::base::{Alphabet, Nucleotide};

/// Parameters for defining the initial repeat sequence structure.
#[derive(Debug, Clone)]
pub struct RepeatStructure {
    /// Alphabet used for sequences
    pub alphabet: Alphabet,
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
        alphabet: Alphabet,
        init_base: Nucleotide,
        ru_length: usize,
        rus_per_hor: usize,
        hors_per_chr: usize,
        chrs_per_hap: usize,
    ) -> Self {
        Self {
            alphabet,
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
#[derive(Debug, Clone)]
pub struct MutationConfig {
    /// Substitution model for point mutations
    pub model: SubstitutionModel,
}

impl MutationConfig {
    /// Create new mutation configuration.
    pub fn new(model: SubstitutionModel) -> Self {
        Self { model }
    }

    /// Create with uniform mutation rate.
    pub fn uniform(alphabet: Alphabet, rate: f64) -> Result<Self, String> {
        let model = SubstitutionModel::uniform(alphabet, rate)
            .map_err(|e| format!("{}", e))?;
        Ok(Self { model })
    }
}

/// Parameters for recombination processes.
#[derive(Debug, Clone)]
pub struct RecombinationConfig {
    /// Recombination parameters
    pub params: RecombinationParams,
}

impl RecombinationConfig {
    /// Create new recombination configuration.
    pub fn new(params: RecombinationParams) -> Self {
        Self { params }
    }

    /// Create with standard parameters.
    pub fn standard(
        break_prob: f64,
        crossover_prob: f64,
        gc_extension_prob: f64,
    ) -> Result<Self, String> {
        let params = RecombinationParams::new(break_prob, crossover_prob, gc_extension_prob)
            .map_err(|e| format!("{}", e))?;
        Ok(Self { params })
    }
}

/// Parameters for fitness/selection.
#[derive(Debug, Clone)]
pub struct FitnessConfig {
    /// GC content fitness function
    pub gc_content: Option<GCContentFitness>,
    /// Length-based fitness function
    pub length: Option<LengthFitness>,
    /// Sequence similarity fitness function
    pub similarity: Option<SequenceSimilarityFitness>,
}

impl FitnessConfig {
    /// Create new fitness configuration.
    pub fn new(
        gc_content: Option<GCContentFitness>,
        length: Option<LengthFitness>,
        similarity: Option<SequenceSimilarityFitness>,
    ) -> Self {
        Self {
            gc_content,
            length,
            similarity,
        }
    }

    /// Create neutral fitness (no selection).
    pub fn neutral() -> Self {
        Self {
            gc_content: None,
            length: None,
            similarity: None,
        }
    }

    /// Check if all fitness components are None (neutral).
    pub fn is_neutral(&self) -> bool {
        self.gc_content.is_none() && self.length.is_none() && self.similarity.is_none()
    }
}

/// High-level simulation parameters.
#[derive(Debug, Clone)]
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
    pub fn new(
        population_size: usize,
        total_generations: usize,
        seed: Option<u64>,
    ) -> Self {
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
        let structure = RepeatStructure::new(
            Alphabet::dna(),
            Nucleotide::A,
            171,
            12,
            100,
            1,
        );
        
        assert_eq!(structure.ru_length, 171);
        assert_eq!(structure.rus_per_hor, 12);
        assert_eq!(structure.hors_per_chr, 100);
        assert_eq!(structure.chrs_per_hap, 1);
    }

    #[test]
    fn test_repeat_structure_chr_length() {
        let structure = RepeatStructure::new(
            Alphabet::dna(),
            Nucleotide::A,
            10,
            5,
            20,
            1,
        );
        
        // 10 * 5 * 20 = 1000
        assert_eq!(structure.chr_length(), 1000);
    }

    #[test]
    fn test_mutation_config_uniform() {
        let config = MutationConfig::uniform(Alphabet::dna(), 0.001).unwrap();
        // Just check it was created successfully
        assert!(matches!(config.model, SubstitutionModel { .. }));
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
        assert!(config.similarity.is_none());
    }

    #[test]
    fn test_fitness_config_with_components() {
        let gc = GCContentFitness::new(0.5, 2.0).unwrap();
        let config = FitnessConfig::new(Some(gc), None, None);
        
        assert!(!config.is_neutral());
        assert!(config.gc_content.is_some());
        assert!(config.length.is_none());
    }

    #[test]
    fn test_simulation_config_new() {
        let config = SimulationConfig::new(100, 1000, Some(42));
        
        assert_eq!(config.population_size, 100);
        assert_eq!(config.total_generations, 1000);
        assert_eq!(config.seed, Some(42));
    }
}
