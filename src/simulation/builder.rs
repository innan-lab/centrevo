//! Builder pattern for creating simulations.
//!
//! Provides a fluent API for configuring and creating simulations with
//! sensible defaults and comprehensive validation.

use crate::base::Nucleotide;
use crate::simulation::initialization::SequenceInput;
use crate::simulation::{
    FitnessConfig, MutationConfig, RecombinationConfig, RepeatStructure, Simulation,
    SimulationConfig,
};

/// Builder for constructing Simulation instances with a fluent API.
///
/// # Examples
///
/// ```
/// use centrevo::simulation::SimulationBuilder;
/// use centrevo::base::Nucleotide;
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
/// # From imported sequences
///
/// ```no_run
/// use centrevo::simulation::SimulationBuilder;
///
/// // From imported sequences (requires sequences.fasta to exist)
/// let sim = SimulationBuilder::new()
///     .population_size(50)
///     .generations(100)
///     .init_from_fasta("sequences.fasta")
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

    // Initialization mode
    init_mode: InitMode,

    // Evolutionary parameters (with defaults)
    mutation_rate: f64,                           // Default: 0.0 (no mutation)
    recombination_rates: Option<(f64, f64, f64)>, // Default: None (no recombination)
    fitness: FitnessConfig,                       // Default: neutral
    seed: Option<u64>,                            // Default: None (random)
}

/// Initialization mode for the simulation.
#[derive(Clone)]
enum InitMode {
    /// Uniform initialization with specified base (default: A)
    Uniform(Nucleotide),
    /// Random initialization (each position gets random base)
    Random,
    /// Initialize from FASTA file
    FromFasta(String),
    /// Initialize from JSON file or string
    FromJson(String),
    /// Initialize from checkpoint (sequences only, not full resume)
    FromCheckpoint {
        db_path: String,
        sim_id: String,
        generation: Option<usize>,
    },
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
            init_mode: InitMode::Uniform(Nucleotide::A),
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
        self.init_mode = InitMode::Uniform(base);
        self
    }

    /// Initialize with random sequences (each position gets random base from alphabet).
    pub fn init_random(mut self) -> Self {
        self.init_mode = InitMode::Random;
        self
    }

    /// Initialize from a FASTA file.
    pub fn init_from_fasta(mut self, path: impl Into<String>) -> Self {
        self.init_mode = InitMode::FromFasta(path.into());
        self
    }

    /// Initialize from JSON (file path or JSON string).
    pub fn init_from_json(mut self, input: impl Into<String>) -> Self {
        self.init_mode = InitMode::FromJson(input.into());
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
        self.init_mode = InitMode::FromCheckpoint {
            db_path: db_path.into(),
            sim_id: sim_id.into(),
            generation,
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
        match self.init_mode {
            InitMode::Uniform(base) => {
                // For uniform, we need repeat structure
                let ru_length = self.ru_length.ok_or(BuilderError::MissingRequired(
                    "ru_length (via repeat_structure)",
                ))?;
                let rus_per_hor = self.rus_per_hor.ok_or(BuilderError::MissingRequired(
                    "rus_per_hor (via repeat_structure)",
                ))?;
                let hors_per_chr = self.hors_per_chr.ok_or(BuilderError::MissingRequired(
                    "hors_per_chr (via repeat_structure)",
                ))?;

                let structure = RepeatStructure::new(
                    base,
                    ru_length,
                    rus_per_hor,
                    hors_per_chr,
                    self.chrs_per_hap,
                );

                Simulation::new(structure, mutation, recombination, self.fitness, config)
                    .map_err(BuilderError::InvalidParameter)
            }

            InitMode::Random => {
                // For random, we need repeat structure
                let ru_length = self.ru_length.ok_or(BuilderError::MissingRequired(
                    "ru_length (via repeat_structure)",
                ))?;
                let rus_per_hor = self.rus_per_hor.ok_or(BuilderError::MissingRequired(
                    "rus_per_hor (via repeat_structure)",
                ))?;
                let hors_per_chr = self.hors_per_chr.ok_or(BuilderError::MissingRequired(
                    "hors_per_chr (via repeat_structure)",
                ))?;

                let structure = RepeatStructure::new(
                    Nucleotide::A, // Doesn't matter for random init
                    ru_length,
                    rus_per_hor,
                    hors_per_chr,
                    self.chrs_per_hap,
                );

                Simulation::new_random(structure, mutation, recombination, self.fitness, config)
                    .map_err(BuilderError::InvalidParameter)
            }

            InitMode::FromFasta(path) => {
                // For imported sequences, infer structure from sequences
                let source = SequenceInput::Fasta(path);
                let chrs_per_hap = self.chrs_per_hap;
                let fitness = self.fitness;
                Self::build_from_source_static(
                    source,
                    chrs_per_hap,
                    mutation,
                    recombination,
                    fitness,
                    config,
                )
            }

            InitMode::FromJson(input) => {
                let source = SequenceInput::Json(input);
                let chrs_per_hap = self.chrs_per_hap;
                let fitness = self.fitness;
                Self::build_from_source_static(
                    source,
                    chrs_per_hap,
                    mutation,
                    recombination,
                    fitness,
                    config,
                )
            }

            InitMode::FromCheckpoint {
                db_path,
                sim_id,
                generation,
            } => {
                let source = SequenceInput::Database {
                    path: db_path,
                    sim_id,
                    generation,
                };
                let chrs_per_hap = self.chrs_per_hap;
                let fitness = self.fitness;
                Self::build_from_source_static(
                    source,
                    chrs_per_hap,
                    mutation,
                    recombination,
                    fitness,
                    config,
                )
            }
        }
    }

    /// Helper to build simulation from imported sequence source.
    fn build_from_source_static(
        source: SequenceInput,
        chrs_per_hap: usize,
        mutation: MutationConfig,
        recombination: RecombinationConfig,
        fitness: FitnessConfig,
        config: SimulationConfig,
    ) -> Result<Simulation, BuilderError> {
        // For imported sequences, we need to infer structure
        // We'll use a dummy structure that will be replaced by actual sequence structure
        let dummy_structure = RepeatStructure::new(
            Nucleotide::A,
            1, // Will be inferred from sequences
            1,
            1,
            chrs_per_hap,
        );

        Simulation::from_sequences(
            source,
            dummy_structure,
            mutation,
            recombination,
            fitness,
            config,
        )
        .map_err(|e| BuilderError::SequenceImport(e.to_string()))
    }
}

/// Errors that can occur during simulation building.
#[derive(Debug)]
pub enum BuilderError {
    /// A required parameter is missing
    MissingRequired(&'static str),
    /// An invalid parameter value was provided
    InvalidParameter(String),
    /// Error importing sequences
    SequenceImport(String),
}

impl std::fmt::Display for BuilderError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::MissingRequired(param) => {
                write!(f, "Missing required parameter: {param}")
            }
            Self::InvalidParameter(msg) => {
                write!(f, "Invalid parameter: {msg}")
            }
            Self::SequenceImport(msg) => {
                write!(f, "Failed to import sequences: {msg}")
            }
        }
    }
}

impl std::error::Error for BuilderError {}

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
