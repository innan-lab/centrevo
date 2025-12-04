//! Builder pattern for creating simulations.
//!
//! Provides a fluent API for configuring and creating simulations with
//! sensible defaults and comprehensive validation.

use crate::base::Nucleotide;
pub use crate::errors::BuilderError;
use crate::genome::RepeatMap;
use crate::simulation::{
    FitnessConfig, GenerationMode, MutationConfig, RecombinationConfig, SequenceConfig,
    SequenceSource, Simulation, SimulationConfig, UniformRepeatStructure,
};

/// Builder for constructing Simulation instances with a fluent API.
///
/// # Examples
///
/// ```
/// use centrevo_sim::simulation::SimulationBuilder;
/// use centrevo_sim::simulation::UniformRepeatStructure;
/// use centrevo_sim::base::Nucleotide;
///
/// // Simple simulation with defaults
/// let sim = SimulationBuilder::new()
///     .population_size(50)
///     .generations(100)
///     .repeat_structure(171, 12, 10)
///     .build()
///     .unwrap();
///
/// // With mutation and recombination
/// let sim = SimulationBuilder::new()
///     .population_size(50)
///     .generations(100)
///     .repeat_structure(171, 12, 10)
///     .mutation_rate(0.0001)
///     .recombination(0.01, 0.7, 0.1)
///     .seed(42)
///     .build()
///     .unwrap();
/// ```
///
/// ```no_run
/// use centrevo_sim::simulation::SimulationBuilder;
/// use centrevo_sim::simulation::UniformRepeatStructure;
/// use centrevo_sim::base::Nucleotide;
///
/// // From imported sequences (requires sequences.fasta to exist)
/// let structure = UniformRepeatStructure::new(Nucleotide::A, 171, 12, 10, 1);
/// let sim = SimulationBuilder::new()
///     .population_size(50)
///     .generations(100)
///     .init_from_fasta_uniform("sequences.fasta", structure)
///     .mutation_rate(0.0001)
///     .build()
///     .unwrap();
/// ```
#[derive(Clone)]
pub struct SimulationBuilder {
    // Required parameters
    population_size: Option<usize>,
    generations: Option<usize>,

    // Repeat structure (required for random/uniform, ignored for imported)
    ru_length: Option<usize>,
    rus_per_hor: Option<usize>,
    hors_per_chr: Option<usize>,
    chrs_per_hap: usize, // Default: 1

    // Explicit map for init_with_map
    repeat_map: Option<RepeatMap>,

    // Initialization configuration
    config: SequenceConfig,
    // Base for uniform initialization (stored separately as it's part of structure construction)
    init_base: Nucleotide,

    // Evolutionary parameters (with defaults)
    mutation_rate: f64,                           // Default: 0.0 (no mutation)
    recombination_rates: Option<(f64, f64, f64)>, // Default: None (no recombination)
    fitness: FitnessConfig,                       // Default: neutral
    seed: Option<u64>,                            // Default: None (random)
}

impl Default for SimulationBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl SimulationBuilder {
    /// Create a new simulation builder with default values.
    pub fn new() -> Self {
        Self {
            population_size: None,
            generations: None,
            ru_length: None,
            rus_per_hor: None,
            hors_per_chr: None,
            chrs_per_hap: 1,
            repeat_map: None,
            config: SequenceConfig::Generate {
                structure: UniformRepeatStructure::new(Nucleotide::A, 0, 0, 0, 1), // Placeholder
                mode: GenerationMode::Uniform,
            },
            init_base: Nucleotide::A,
            mutation_rate: 0.0,
            recombination_rates: None,
            fitness: FitnessConfig::neutral(),
            seed: None,
        }
    }

    /// Set the population size (required).
    pub fn population_size(mut self, size: usize) -> Self {
        self.population_size = Some(size);
        self
    }

    /// Set the number of generations to run (required).
    pub fn generations(mut self, generations: usize) -> Self {
        self.generations = Some(generations);
        self
    }

    /// Set the repeat structure parameters (required for random/uniform initialization).
    ///
    /// # Arguments
    /// * `ru_length` - Repeat unit length
    /// * `rus_per_hor` - Repeat units per higher-order repeat
    /// * `hors_per_chr` - Higher-order repeats per chromosome
    pub fn repeat_structure(
        mut self,
        ru_length: usize,
        rus_per_hor: usize,
        hors_per_chr: usize,
    ) -> Self {
        self.ru_length = Some(ru_length);
        self.rus_per_hor = Some(rus_per_hor);
        self.hors_per_chr = Some(hors_per_chr);
        self
    }

    /// Set the number of chromosomes per haplotype (default: 1).
    pub fn chromosomes_per_haplotype(mut self, chrs_per_hap: usize) -> Self {
        self.chrs_per_hap = chrs_per_hap;
        self
    }

    /// Initialize with uniform sequences using the specified base.
    pub fn init_uniform(mut self, base: Nucleotide) -> Self {
        self.init_base = base;
        // We set config to Generate Uniform with placeholder, actual structure built in build()
        self.config = SequenceConfig::Generate {
            structure: UniformRepeatStructure::new(base, 0, 0, 0, 1),
            mode: GenerationMode::Uniform,
        };
        self
    }

    /// Initialize with random sequences (each position gets random base from alphabet).
    pub fn init_random(mut self) -> Self {
        self.config = SequenceConfig::Generate {
            structure: UniformRepeatStructure::new(Nucleotide::A, 0, 0, 0, 1), // Placeholder
            mode: GenerationMode::Random,
        };
        self
    }

    /// Initialize from a FASTA file with a BED file for structure.
    pub fn init_from_fasta(mut self, path: impl Into<String>, bed_path: impl Into<String>) -> Self {
        self.config = SequenceConfig::Load {
            source: SequenceSource::Fasta {
                path: path.into(),
                bed_path: Some(bed_path.into()),
            },
            structure: None,
        };
        self
    }

    /// Initialize from a FASTA file with uniform structure.
    pub fn init_from_fasta_uniform(
        mut self,
        path: impl Into<String>,
        structure: UniformRepeatStructure,
    ) -> Self {
        self.config = SequenceConfig::Load {
            source: SequenceSource::Fasta {
                path: path.into(),
                bed_path: None,
            },
            structure: Some(structure),
        };
        self
    }

    /// Initialize from a formatted string.
    pub fn init_from_formatted_string(
        mut self,
        sequence: impl Into<String>,
        hor_delim: char,
        ru_delim: char,
    ) -> Self {
        self.config = SequenceConfig::Load {
            source: SequenceSource::FormattedString {
                sequence: sequence.into(),
                hor_delim,
                ru_delim,
            },
            structure: None,
        };
        self
    }

    /// Initialize with a specific RepeatMap.
    pub fn init_with_map(mut self, map: RepeatMap, _base: Nucleotide) -> Self {
        self.repeat_map = Some(map);
        // TODO: Implement SequenceInput::WithMap or similar
        // For now, we'll error in build() if this is used without proper support
        self
    }

    /// Initialize from JSON (file path or JSON string).
    pub fn init_from_json(mut self, input: impl Into<String>) -> Self {
        self.config = SequenceConfig::Load {
            source: SequenceSource::Json(input.into()),
            structure: None,
        };
        self
    }

    /// Initialize from a checkpoint database (sequences only, new parameters).
    ///
    /// # Arguments
    /// * `db_path` - Path to the checkpoint database
    /// * `sim_id` - Simulation ID to load
    /// * `generation` - Optional generation to load (defaults to last)
    pub fn init_from_checkpoint(
        mut self,
        db_path: impl Into<String>,
        sim_id: impl Into<String>,
        generation: Option<usize>,
    ) -> Self {
        self.config = SequenceConfig::Load {
            source: SequenceSource::Database {
                path: db_path.into(),
                sim_id: sim_id.into(),
                generation,
            },
            structure: None,
        };
        self
    }

    /// Set the mutation rate (default: 0.0).
    pub fn mutation_rate(mut self, rate: f64) -> Self {
        self.mutation_rate = rate;
        self
    }

    /// Set recombination parameters.
    ///
    /// # Arguments
    /// * `break_prob` - Probability of DNA strand break
    /// * `crossover_prob` - Probability of crossover
    /// * `gc_extension_prob` - Probability of gene conversion extension
    pub fn recombination(
        mut self,
        break_prob: f64,
        crossover_prob: f64,
        gc_extension_prob: f64,
    ) -> Self {
        self.recombination_rates = Some((break_prob, crossover_prob, gc_extension_prob));
        self
    }

    /// Set the fitness configuration (default: neutral).
    pub fn fitness(mut self, fitness: FitnessConfig) -> Self {
        self.fitness = fitness;
        self
    }

    /// Set the random seed for reproducibility (default: None = random).
    pub fn seed(mut self, seed: u64) -> Self {
        self.seed = Some(seed);
        self
    }

    /// Build and validate the simulation.
    pub fn build(self) -> Result<Simulation, BuilderError> {
        // Validate required parameters
        let population_size = self
            .population_size
            .ok_or(BuilderError::MissingRequired("population_size"))?;
        let generations = self
            .generations
            .ok_or(BuilderError::MissingRequired("generations"))?;

        // Create simulation config
        let config = SimulationConfig::new(population_size, generations, self.seed);

        // Create mutation config
        let mutation = MutationConfig::uniform(self.mutation_rate)
            .map_err(|e| BuilderError::InvalidParameter(format!("mutation_rate: {e}")))?;

        // Create recombination config
        let recombination = if let Some((break_prob, crossover_prob, gc_extension_prob)) =
            self.recombination_rates
        {
            RecombinationConfig::standard(break_prob, crossover_prob, gc_extension_prob)
                .map_err(|e| BuilderError::InvalidParameter(format!("recombination: {e}")))?
        } else {
            // No recombination (all rates = 0.0)
            RecombinationConfig::standard(0.0, 0.0, 0.0)
                .map_err(|e| BuilderError::InvalidParameter(format!("recombination: {e}")))?
        };

        // Build simulation based on initialization mode
        // Build simulation based on initialization mode
        match &self.config {
            SequenceConfig::Generate { mode, .. } => {
                // For generation (Uniform/Random), we need repeat structure
                let ru_length = self.ru_length.ok_or(BuilderError::MissingRequired(
                    "ru_length (via repeat_structure)",
                ))?;
                let rus_per_hor = self.rus_per_hor.ok_or(BuilderError::MissingRequired(
                    "rus_per_hor (via repeat_structure)",
                ))?;
                let hors_per_chr = self.hors_per_chr.ok_or(BuilderError::MissingRequired(
                    "hors_per_chr (via repeat_structure)",
                ))?;

                let structure = UniformRepeatStructure::new(
                    if *mode == GenerationMode::Uniform {
                        self.init_base
                    } else {
                        Nucleotide::A
                    },
                    ru_length,
                    rus_per_hor,
                    hors_per_chr,
                    self.chrs_per_hap,
                );

                // Update config with actual structure
                let seq_config = SequenceConfig::Generate {
                    structure: structure.clone(),
                    mode: *mode,
                };

                Simulation::from_config(
                    seq_config,
                    Some(structure),
                    mutation,
                    recombination,
                    self.fitness,
                    config,
                )
                .map_err(BuilderError::InvalidParameter)
            }

            SequenceConfig::Load {
                source: _source,
                structure,
            } => {
                // If structure is provided in config (e.g. FastaUniform), use it.
                // Otherwise, check if user provided structure via builder methods (optional for Load).
                // But builder structure params (ru_length etc) are ignored for Load unless we want to enforce validation.
                // Current logic: use structure from config if present.

                // If user tried to use init_with_map (not implemented fully)
                if self.repeat_map.is_some() {
                    return Err(BuilderError::InvalidParameter(
                        "init_with_map not fully implemented yet".into(),
                    ));
                }

                Self::build_from_config_static(
                    self.config.clone(),
                    structure.clone(),
                    mutation,
                    recombination,
                    self.fitness,
                    config,
                )
            }
        }
    }

    /// Helper to build simulation from configuration.
    fn build_from_config_static(
        seq_config: SequenceConfig,
        structure: Option<UniformRepeatStructure>,
        mutation: MutationConfig,
        recombination: RecombinationConfig,
        fitness: FitnessConfig,
        config: SimulationConfig,
    ) -> Result<Simulation, BuilderError> {
        Simulation::from_config(
            seq_config,
            structure,
            mutation,
            recombination,
            fitness,
            config,
        )
        .map_err(|e| BuilderError::SequenceImport(e.to_string()))
    }
}

// Removed BuilderError definition, imported from crate::errors

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_builder_minimal() {
        let sim = SimulationBuilder::new()
            .population_size(10)
            .generations(50)
            .repeat_structure(10, 5, 2)
            .build();

        assert!(sim.is_ok());
        let sim = sim.unwrap();
        assert_eq!(sim.population().size(), 10);
        assert_eq!(sim.generation(), 0);
    }

    #[test]
    fn test_builder_with_mutation() {
        let sim = SimulationBuilder::new()
            .population_size(10)
            .generations(50)
            .repeat_structure(10, 5, 2)
            .mutation_rate(0.001)
            .build();

        assert!(sim.is_ok());
    }

    #[test]
    fn test_builder_with_recombination() {
        let sim = SimulationBuilder::new()
            .population_size(10)
            .generations(50)
            .repeat_structure(10, 5, 2)
            .recombination(0.01, 0.7, 0.1)
            .build();

        assert!(sim.is_ok());
    }

    #[test]
    fn test_builder_with_seed() {
        let sim = SimulationBuilder::new()
            .population_size(10)
            .generations(50)
            .repeat_structure(10, 5, 2)
            .seed(42)
            .build();

        assert!(sim.is_ok());
    }

    #[test]
    fn test_builder_missing_population_size() {
        let sim = SimulationBuilder::new()
            .generations(50)
            .repeat_structure(10, 5, 2)
            .build();

        assert!(sim.is_err());
        match sim.unwrap_err() {
            BuilderError::MissingRequired(param) => {
                assert_eq!(param, "population_size");
            }
            _ => panic!("Expected MissingRequired error"),
        }
    }

    #[test]
    fn test_builder_missing_generations() {
        let sim = SimulationBuilder::new()
            .population_size(10)
            .repeat_structure(10, 5, 2)
            .build();

        assert!(sim.is_err());
        match sim.unwrap_err() {
            BuilderError::MissingRequired(param) => {
                assert_eq!(param, "generations");
            }
            _ => panic!("Expected MissingRequired error"),
        }
    }

    #[test]
    fn test_builder_missing_repeat_structure() {
        let sim = SimulationBuilder::new()
            .population_size(10)
            .generations(50)
            .build();

        assert!(sim.is_err());
        match sim.unwrap_err() {
            BuilderError::MissingRequired(_) => {}
            _ => panic!("Expected MissingRequired error"),
        }
    }

    #[test]
    fn test_builder_invalid_mutation_rate() {
        let sim = SimulationBuilder::new()
            .population_size(10)
            .generations(50)
            .repeat_structure(10, 5, 2)
            .mutation_rate(-0.1) // Invalid: negative
            .build();

        assert!(sim.is_err());
    }

    #[test]
    fn test_builder_uniform_with_different_base() {
        let sim = SimulationBuilder::new()
            .population_size(10)
            .generations(50)
            .repeat_structure(10, 5, 2)
            .init_uniform(Nucleotide::G)
            .build();

        assert!(sim.is_ok());
    }

    #[test]
    fn test_builder_random_init() {
        let sim = SimulationBuilder::new()
            .population_size(10)
            .generations(50)
            .repeat_structure(10, 5, 2)
            .init_random()
            .seed(42)
            .build();

        if let Err(e) = &sim {
            eprintln!("Builder error: {e:?}");
        }
        assert!(sim.is_ok());
    }

    #[test]
    fn test_builder_all_options() {
        let sim = SimulationBuilder::new()
            .population_size(20)
            .generations(100)
            .repeat_structure(171, 12, 10)
            .chromosomes_per_haplotype(2)
            .mutation_rate(0.0001)
            .recombination(0.01, 0.7, 0.1)
            .seed(12345)
            .build();

        assert!(sim.is_ok());
        let sim = sim.unwrap();
        assert_eq!(sim.population().size(), 20);
    }
}
