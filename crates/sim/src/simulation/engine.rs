//! Simulation engine for evolutionary processes.
//!
//! This module provides the main simulation loop that orchestrates mutation,
//! recombination, selection, and reproduction across generations.

use crate::errors::RecombinationError;
use crate::evolution::RecombinationType;

use crate::genome::{Chromosome, Individual};
use crate::simulation::{
    Configuration, ExecutionConfig, FitnessConfig, InitializationConfig, MutationConfig,
    Population, RecombinationConfig, SequenceSource, UniformRepeatStructure,
};
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;

/// Main simulation engine.
#[derive(Debug)]
pub struct Simulation {
    /// Current population
    population: Population,
    /// Master configuration
    config: Configuration,
    /// Random number generator (using Xoshiro256++ for better performance)
    rng: Xoshiro256PlusPlus,
}

impl Simulation {
    /// Create a new simulation from configuration.
    ///
    /// This allows initializing a simulation with sequences based on the provided configuration.
    ///
    /// # Arguments
    ///
    /// * `config` - Master simulation configuration
    ///
    /// # Returns
    ///
    /// A `Simulation` instance initialized with the provided configuration.
    pub fn new(config: Configuration) -> Result<Self, String> {
        // Create RNG from seed or thread_rng
        let mut rng = if let Some(seed) = config.execution.seed {
            Xoshiro256PlusPlus::seed_from_u64(seed)
        } else {
            Xoshiro256PlusPlus::from_seed(rand::rng().random())
        };

        // Initialize individuals from config
        let mut individuals = crate::simulation::initialize(
            &config.initialization,
            config.execution.population_size,
            &mut rng,
            &config.execution.codec,
        )
        .map_err(|e| format!("Failed to initialize from sequences: {e}"))?;

        // Compute fitness for initial population
        // This ensures the invariant that all individuals have cached fitness
        Simulation::compute_population_fitness(&mut individuals, &config.evolution.fitness);

        let population = Population::new("pop0", individuals);

        Ok(Self {
            population,
            config,
            rng,
        })
    }

    /// Helper to compute valid initial fitness for a set of individuals.
    fn compute_population_fitness(individuals: &mut [Individual], fitness: &FitnessConfig) {
        individuals.par_iter_mut().for_each(|ind| {
            fitness.update_cached_fitness(ind);
        });
    }
    /// Apply a single recombination `event` to a pair of chromosomes.
    ///
    /// This helper is used by `apply_recombination` and is kept as a separate
    /// private method to make testing and error handling centralized.
    fn apply_event_to_pair(
        params: &crate::evolution::RecombinationModel,
        chr1: &mut Chromosome,
        chr2: &mut Chromosome,
        event: RecombinationType,
    ) -> Result<(), String> {
        match event {
            RecombinationType::None => Ok(()),
            RecombinationType::Crossover { pos1, pos2 } => {
                let (new1, new2) = params
                    .crossover(chr1, chr2, pos1, pos2)
                    .map_err(|e| format!("Crossover failed: {e}"))?;

                *chr1 = new1;
                *chr2 = new2;
                Ok(())
            }
            RecombinationType::GeneConversion {
                donor_start,
                recipient_start,
                length,
            } => {
                match params.gene_conversion(chr2, chr1, recipient_start, donor_start, length) {
                    Ok(new2) => {
                        *chr2 = new2;
                        Ok(())
                    }
                    Err(e) => match e {
                        RecombinationError::InvalidRange { .. }
                        | RecombinationError::InvalidPosition { .. } => {
                            // Try to clamp to remaining lengths
                            if recipient_start >= chr2.len() || donor_start >= chr1.len() {
                                // Nothing to do
                                return Ok(());
                            }
                            let recipient_remaining = chr2.len() - recipient_start;
                            let donor_remaining = chr1.len() - donor_start;
                            let clamped_len = length.min(recipient_remaining).min(donor_remaining);
                            if clamped_len == 0 {
                                return Ok(());
                            }

                            if let Ok(new2) = params.gene_conversion(
                                chr2,
                                chr1,
                                recipient_start,
                                donor_start,
                                clamped_len,
                            ) {
                                *chr2 = new2;
                                Ok(())
                            } else {
                                Ok(())
                            }
                        }
                        other => Err(format!("Gene conversion failed: {other}")),
                    },
                }
            }
        }
    }

    /// Resume a simulation from a checkpoint in the database.
    ///
    /// This loads the complete simulation state from the last recorded checkpoint,
    /// including population state, RNG state, and all configuration parameters.
    /// The simulation can then continue exactly as if it had never stopped.
    ///
    /// # Arguments
    ///
    /// * `db_path` - Path to the SQLite database containing the checkpoint
    /// * `sim_id` - Simulation ID to resume
    ///
    /// # Returns
    ///
    /// A `Simulation` instance ready to continue from the checkpoint, or an error
    /// if the checkpoint is invalid or configuration cannot be loaded.
    pub fn from_checkpoint(
        db_path: impl AsRef<std::path::Path>,
        sim_id: &str,
    ) -> Result<Self, String> {
        use crate::storage::QueryBuilder;

        // Open database for querying
        let query =
            QueryBuilder::new(&db_path).map_err(|e| format!("Failed to open database: {e}"))?;

        // Get the latest checkpoint
        let checkpoint = query
            .get_latest_checkpoint(sim_id)
            .map_err(|e| format!("Failed to load checkpoint: {e}"))?;

        // Verify sim_id matches
        if checkpoint.sim_id != sim_id {
            return Err(format!(
                "Checkpoint sim_id mismatch: expected '{expected}', found '{found}'",
                expected = sim_id,
                found = checkpoint.sim_id
            ));
        }

        // Load complete configuration
        // Load complete configuration
        let mut config = query
            .get_full_config(sim_id)
            .map_err(|e| format!("Failed to load configuration: {e}"))?;

        // Adapt initialization config to load from database
        let db_path_str = db_path.as_ref().to_string_lossy().to_string();
        config.initialization = InitializationConfig::Load {
            source: SequenceSource::Database {
                path: db_path_str,
                sim_id: sim_id.to_string(),
                generation: Some(checkpoint.generation),
            },
        };

        // Create simulation instance (loads population via initialize())
        // Note: We use from_config logic here but we need to advance the population state
        // correctly after initialization.
        // Actually, initialize() for Database source will load the individuals.
        let mut sim = Self::new(config)?;

        // Override population name to match generation
        sim.population = Population::new(
            format!("pop_{gen}", gen = checkpoint.generation),
            sim.population.individuals().to_vec(),
        );

        // Force set the correct generation
        // The default new() starts at 0.
        for _ in 0..checkpoint.generation {
            sim.population.increment_generation();
        }

        // Restore RNG state
        let rng: Xoshiro256PlusPlus = bincode::deserialize(&checkpoint.rng_state)
            .map_err(|e| format!("Failed to restore RNG state: {e}"))?;
        sim.rng = rng;

        // Close query builder
        query
            .close()
            .map_err(|e| format!("Failed to close database: {e}"))?;

        Ok(sim)
    }

    /// Get the current population.
    pub fn population(&self) -> &Population {
        &self.population
    }

    /// Get mutable access to the population.
    pub fn population_mut(&mut self) -> &mut Population {
        &mut self.population
    }

    /// Get the current generation number.
    pub fn generation(&self) -> usize {
        self.population.generation()
    }

    /// Get reference to the full master configuration.
    pub fn configuration(&self) -> &Configuration {
        &self.config
    }

    // Deprecated accessors - kept for compatibility but should be migrated
    /// Get reference to simulation configuration.
    pub fn simulation_config(&self) -> &ExecutionConfig {
        &self.config.execution
    }

    /// Get reference to repeat structure (if available).
    pub fn structure_config(&self) -> Option<&UniformRepeatStructure> {
        match &self.config.initialization {
            InitializationConfig::Generate { structure, .. } => Some(structure),
            InitializationConfig::Load { .. } => None,
        }
    }

    /// Get reference to mutation configuration.
    pub fn mutation_config(&self) -> &MutationConfig {
        &self.config.evolution.mutation
    }

    /// Get reference to recombination configuration.
    pub fn recombination_config(&self) -> &RecombinationConfig {
        &self.config.evolution.recombination
    }

    /// Get reference to fitness configuration.
    pub fn fitness_config(&self) -> &FitnessConfig {
        &self.config.evolution.fitness
    }

    /// Get reference to the full master configuration.
    pub fn config(&self) -> &Configuration {
        &self.config
    }

    /// Get the current RNG state for checkpointing.
    /// Returns the internal state as bytes.
    pub fn rng_state_bytes(&self) -> Vec<u8> {
        bincode::serialize(&self.rng).expect("Failed to serialize RNG state")
    }

    /// Set the RNG state from a checkpoint.
    /// Takes a byte array representing the serialized RNG state.
    pub fn set_rng_from_bytes(&mut self, bytes: &[u8]) -> Result<(), String> {
        self.rng = bincode::deserialize(bytes)
            .map_err(|e| format!("Failed to deserialize RNG state: {e}"))?;
        Ok(())
    }

    /// Produce a single gamete from a parent.
    ///
    /// This process simulates:
    /// 1. Meiosis (Recombination between parent's haplotypes)
    /// 2. Assortment (Random selection of one recombinant haplotype)
    /// 3. Mutation (Applied to the gamete)
    fn produce_gamete(
        &self,
        parent: &Individual,
        rng: &mut Xoshiro256PlusPlus,
    ) -> Result<crate::genome::Haplotype, String> {
        // 1. & 2. Meiosis and Assortment
        // We start by cloning the parent's haplotypes to avoid mutating the population state.
        // We need mutable copies to perform recombination.
        let mut hap1 = parent.haplotype1().clone();
        let mut hap2 = parent.haplotype2().clone();

        // Apply recombination between the two haplotypes
        // Iterate over shared chromosomes
        let num_chromosomes = hap1.len().min(hap2.len());
        for chr_idx in 0..num_chromosomes {
            if let (Some(chr1), Some(chr2)) = (hap1.get_mut(chr_idx), hap2.get_mut(chr_idx)) {
                // Sample recombination events
                let mut events = self
                    .config
                    .evolution
                    .recombination
                    .params
                    .sample_events(chr1, chr2, rng);

                // Process events from right to left (descending position)
                events.reverse();

                // Apply events
                for event in events {
                    Simulation::apply_event_to_pair(
                        &self.config.evolution.recombination.params,
                        chr1,
                        chr2,
                        event,
                    )?;
                }
            }
        }

        // Randomly select one of the recombinant haplotypes to be the gamete
        // Use gen::<f64>() which is faster than random_bool
        let mut gamete = if rng.random::<f64>() < 0.5 {
            hap1
        } else {
            hap2
        };

        // 3. Mutation
        // Apply mutation to the chosen gamete
        for chr in gamete.chromosomes_mut() {
            let seq = chr.sequence_mut();

            // Apply substitutions settings thread-local RNG
            self.config
                .evolution
                .mutation
                .substitution
                .mutate_sequence_sparse(seq, rng);

            // Apply indels if configured
            if let Some(indel_model) = &self.config.evolution.mutation.indel {
                indel_model.apply_indels(seq, rng);
            }
        }

        // Clear cached fitness since the sequence has changed
        gamete.set_cached_fitness(crate::base::FitnessValue::default());

        Ok(gamete)
    }

    /// Advance simulation by one generation.
    ///
    /// The parallel iteration occurs over parent pairs, where each task:
    /// 1. Generates two gametes (one from each parent via meiosis + mutation)
    /// 2. Combines gametes to form a new diploid individual
    /// 3. Computes the individual's fitness based on configured fitness functions
    pub fn step(&mut self) -> Result<(), String> {
        // 1. Select parent pairs using cached fitness from current population
        let parent_pairs = self
            .population
            .select_parents(&mut self.rng, self.config.execution.population_size);

        // Pre-allocate generation string for offspring IDs
        let gen_str = format!("ind_gen{}_", self.generation() + 1);

        // Generate seeds for parallel offspring generation (ensures reproducibility)
        let seeds: Vec<u64> = (0..self.config.execution.population_size)
            .map(|_| self.rng.random())
            .collect();

        // Get immutable references for parallel access
        let population = &self.population;
        let fitness_config = &self.config.evolution.fitness;

        // 2. Generate offspring in parallel: each task processes one parent pair
        let offspring: Result<Vec<Individual>, String> = parent_pairs
            .par_iter()
            .zip(seeds.par_iter())
            .enumerate()
            .map(|(i, ((parent1_idx, parent2_idx), &seed))| {
                // Each parallel task gets its own RNG for reproducibility
                let mut local_rng = Xoshiro256PlusPlus::seed_from_u64(seed);

                let parent1 = population.get(*parent1_idx).unwrap();
                let parent2 = population.get(*parent2_idx).unwrap();

                // Generate two gametes (meiosis + mutation), one from each parent
                let gamete1 = self.produce_gamete(parent1, &mut local_rng)?;
                let gamete2 = self.produce_gamete(parent2, &mut local_rng)?;

                // Combine gametes to form new diploid individual
                let id = format!("{gen_str}{i}");
                let mut new_ind = Individual::new(id, gamete1, gamete2);

                // Compute intrinsic fitness immediately (born with it)
                fitness_config.update_cached_fitness(&mut new_ind);

                Ok(new_ind)
            })
            .collect();

        // 3. Collect offspring (Result -> Vec)
        let offspring = offspring?;

        // 4. Replace population with new generation
        self.population.set_individuals(offspring);

        // 5. Increment generation counter
        self.population.increment_generation();

        Ok(())
    }

    /// Run simulation for the configured number of generations.
    pub fn run(&mut self) -> Result<(), String> {
        for _ in 0..self.config.execution.total_generations {
            self.step()?;
        }
        Ok(())
    }

    /// Run simulation for a specific number of generations.
    pub fn run_for(&mut self, generations: usize) -> Result<(), String> {
        for _ in 0..generations {
            self.step()?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::base::Nucleotide;
    use crate::evolution::SubstitutionModel;
    use crate::simulation::{
        Configuration, EvolutionConfig, ExecutionConfig, FitnessConfig, GenerationMode,
        SimulationBuilder,
    };

    /// Helper function to create a test simulation with standard configuration.
    ///
    /// Creates a simulation with:
    /// - Population size: 10
    /// - Generations: 5
    /// - Repeat structure: 10 ru_length, 5 rus_per_hor, 10 hors_per_chr
    /// - Mutation rate: 0.001
    /// - Recombination: standard(0.01, 0.7, 0.1)
    /// - Seed: 42 (reproducible)
    fn create_test_simulation() -> Simulation {
        SimulationBuilder::new()
            .population_size(10)
            .generations(5)
            .repeat_structure(10, 5, 10)
            .mutation_rate(0.001)
            .recombination(0.01, 0.7, 0.1)
            .seed(42)
            .build()
            .unwrap()
    }

    // ============================================================================
    // Tests for private helper methods
    // ============================================================================

    #[test]
    fn test_apply_event_to_pair_gene_conversion_clamps() {
        // Create a small simulation with trivial configs
        let sim = SimulationBuilder::new()
            .population_size(1)
            .generations(1)
            .repeat_structure(1, 1, 6) // chr len 6
            .mutation_rate(0.0)
            .recombination(0.0, 0.5, 0.5)
            .seed(42)
            .build()
            .unwrap();

        // Chron 1: donor (length 6), Chron 2: recipient (length 4)
        let mut donor = Chromosome::uniform("d", Nucleotide::A, 1, 1, 6);
        let mut recipient = Chromosome::uniform("r", Nucleotide::T, 1, 1, 4);

        // Create a gene conversion event where length exceeds recipient remaining length
        let event = crate::evolution::RecombinationType::GeneConversion {
            donor_start: 0,     // donor start
            recipient_start: 1, // recipient start
            length: 5,          // length > recipient remaining (3)
        };

        // Apply via helper — should clamp to 3 and succeed
        let params = sim.config.evolution.recombination.params.clone();
        Simulation::apply_event_to_pair(&params, &mut donor, &mut recipient, event).unwrap();

        assert_eq!(recipient.to_string(), "TAAA");
    }

    #[test]
    fn test_apply_event_to_pair_gene_conversion_skip_when_out_of_bounds() {
        let sim = SimulationBuilder::new()
            .population_size(1)
            .generations(1)
            .repeat_structure(1, 1, 4)
            .mutation_rate(0.0)
            .recombination(0.0, 0.5, 0.5)
            .seed(42)
            .build()
            .unwrap();

        let mut donor = Chromosome::uniform("d", Nucleotide::A, 1, 1, 4);
        let mut recipient = Chromosome::uniform("r", Nucleotide::T, 1, 1, 4);

        // The recipient start is out of bounds
        let event = crate::evolution::RecombinationType::GeneConversion {
            donor_start: 0,
            recipient_start: 10, // invalid
            length: 2,
        };

        let params = sim.config.evolution.recombination.params.clone();
        Simulation::apply_event_to_pair(&params, &mut donor, &mut recipient, event).unwrap();

        // No changes to recipient - still TTTT
        assert_eq!(recipient.to_string(), "TTTT");
    }

    // ============================================================================
    // Tests for from_config
    // ============================================================================

    #[test]
    fn test_from_config_creates_simulation() {
        let sim = create_test_simulation();

        assert_eq!(sim.population().size(), 10);
        assert_eq!(sim.generation(), 0);
    }

    #[test]
    fn test_from_config_creates_correct_initial_population_structure() {
        let sim = create_test_simulation();

        // Check that all individuals have the correct structure
        for individual in sim.population().individuals() {
            assert_eq!(individual.haplotype1().len(), 1);
            assert_eq!(individual.haplotype2().len(), 1);

            let chr1 = individual.haplotype1().get(0).unwrap();
            let chr2 = individual.haplotype2().get(0).unwrap();

            assert_eq!(chr1.len(), 500); // 10 * 5 * 10
            assert_eq!(chr2.len(), 500);
        }
    }

    #[test]
    fn test_from_config_with_uniform_initialization() {
        let sim = SimulationBuilder::new()
            .population_size(1)
            .generations(1)
            .repeat_structure(10, 5, 10)
            .chromosomes_per_haplotype(2)
            .init_uniform(Nucleotide::C)
            .mutation_rate(0.0)
            .recombination(0.0, 0.0, 0.0)
            .build()
            .unwrap();

        let ind = sim.population().get(0).unwrap();

        assert_eq!(ind.id(), "ind_0");
        assert_eq!(ind.haplotype1().len(), 2);
        assert_eq!(ind.haplotype2().len(), 2);

        // Check that all bases are C
        for chr in ind.haplotype1().chromosomes() {
            for i in 0..chr.len() {
                assert_eq!(chr.sequence().get(i), Some(Nucleotide::C));
            }
        }
    }

    #[test]
    fn test_from_config_golden_path() {
        let structure = UniformRepeatStructure::new(Nucleotide::A, 5, 10, 20, 1);
        let seq_config = InitializationConfig::Generate {
            structure: structure.clone(),
            mode: GenerationMode::Uniform,
        };
        let mutation = MutationConfig::uniform(0.001).unwrap();
        let recombination = RecombinationConfig::standard(0.01, 0.7, 0.1).unwrap();
        let fitness = FitnessConfig::neutral();
        let config = ExecutionConfig::new(100, 10, Some(123));

        let evolution = crate::simulation::EvolutionConfig {
            mutation,
            recombination,
            fitness,
        };

        let configuration = Configuration {
            execution: config,
            evolution,
            initialization: seq_config,
        };

        let sim = Simulation::new(configuration).unwrap();

        assert_eq!(sim.population().size(), 100);
        assert_eq!(sim.generation(), 0);
        assert_eq!(sim.simulation_config().population_size, 100);
        assert_eq!(sim.simulation_config().total_generations, 10);
    }

    #[test]
    fn test_from_config_with_seed_is_reproducible() {
        let structure = UniformRepeatStructure::new(Nucleotide::A, 5, 10, 20, 1);
        let seq_config1 = InitializationConfig::Generate {
            structure: structure.clone(),
            mode: GenerationMode::Uniform,
        };
        let seq_config2 = InitializationConfig::Generate {
            structure: structure.clone(),
            mode: GenerationMode::Uniform,
        };
        let mutation1 = MutationConfig::uniform(0.001).unwrap();
        let mutation2 = MutationConfig::uniform(0.001).unwrap();
        let recombination1 = RecombinationConfig::standard(0.01, 0.7, 0.1).unwrap();
        let recombination2 = RecombinationConfig::standard(0.01, 0.7, 0.1).unwrap();
        let fitness1 = FitnessConfig::neutral();
        let fitness2 = FitnessConfig::neutral();
        let config1 = ExecutionConfig::new(10, 5, Some(999));
        let config2 = ExecutionConfig::new(10, 5, Some(999));

        let sim1 = Simulation::new(Configuration {
            execution: config1,
            evolution: crate::simulation::EvolutionConfig {
                mutation: mutation1,
                recombination: recombination1,
                fitness: fitness1,
            },
            initialization: seq_config1,
        })
        .unwrap();

        let sim2 = Simulation::new(Configuration {
            execution: config2,
            evolution: crate::simulation::EvolutionConfig {
                mutation: mutation2,
                recombination: recombination2,
                fitness: fitness2,
            },
            initialization: seq_config2,
        })
        .unwrap();

        // Check that initial populations are identical
        for i in 0..10 {
            let ind1 = sim1.population().get(i).unwrap();
            let ind2 = sim2.population().get(i).unwrap();

            assert_eq!(
                ind1.haplotype1().get(0).unwrap().sequence().to_string(),
                ind2.haplotype1().get(0).unwrap().sequence().to_string()
            );
        }
    }

    #[test]
    fn test_from_config_different_seeds_produce_different_populations() {
        let structure = UniformRepeatStructure::new(Nucleotide::A, 5, 10, 20, 1);
        let seq_config1 = InitializationConfig::Generate {
            structure: structure.clone(),
            mode: GenerationMode::Random,
        };
        let seq_config2 = InitializationConfig::Generate {
            structure: structure.clone(),
            mode: GenerationMode::Random,
        };
        let mutation1 = MutationConfig::uniform(0.001).unwrap();
        let mutation2 = MutationConfig::uniform(0.001).unwrap();
        let recombination1 = RecombinationConfig::standard(0.01, 0.7, 0.1).unwrap();
        let recombination2 = RecombinationConfig::standard(0.01, 0.7, 0.1).unwrap();
        let fitness1 = FitnessConfig::neutral();
        let fitness2 = FitnessConfig::neutral();
        let config1 = ExecutionConfig::new(10, 5, Some(111));
        let config2 = ExecutionConfig::new(10, 5, Some(222));

        let sim1 = Simulation::new(Configuration {
            execution: config1,
            evolution: crate::simulation::EvolutionConfig {
                mutation: mutation1,
                recombination: recombination1,
                fitness: fitness1,
            },
            initialization: seq_config1,
        })
        .unwrap();

        let sim2 = Simulation::new(Configuration {
            execution: config2,
            evolution: crate::simulation::EvolutionConfig {
                mutation: mutation2,
                recombination: recombination2,
                fitness: fitness2,
            },
            initialization: seq_config2,
        })
        .unwrap();

        // Check that at least one individual is different
        let mut found_difference = false;
        for i in 0..10 {
            let ind1 = sim1.population().get(i).unwrap();
            let ind2 = sim2.population().get(i).unwrap();

            if ind1.haplotype1().get(0).unwrap().sequence().to_string()
                != ind2.haplotype1().get(0).unwrap().sequence().to_string()
            {
                found_difference = true;
                break;
            }
        }

        assert!(
            found_difference,
            "Different seeds should produce different populations"
        );
    }

    #[test]
    fn test_from_config_computes_initial_fitness() {
        use crate::evolution::GCContentFitness;

        let structure = UniformRepeatStructure::new(Nucleotide::A, 5, 10, 20, 1);
        let seq_config = InitializationConfig::Generate {
            structure: structure.clone(),
            mode: GenerationMode::Random,
        };
        let mutation = MutationConfig::uniform(0.0).unwrap();
        let recombination = RecombinationConfig::standard(0.0, 0.0, 0.0).unwrap();

        // Create fitness config with GC content selection
        let gc_fitness = GCContentFitness::new(0.5, 1.0).unwrap();
        let fitness = FitnessConfig::new(Some(gc_fitness), None, None, None);

        let execution = ExecutionConfig::new(10, 5, Some(42));

        let config = Configuration {
            execution,
            evolution: EvolutionConfig {
                mutation,
                recombination,
                fitness,
            },
            initialization: seq_config,
        };

        let sim = Simulation::new(config).unwrap();

        // All individuals should have fitness computed
        for ind in sim.population().individuals() {
            assert!(*ind.cached_fitness().unwrap_or_default() > 0.0);
        }
    }

    // ============================================================================
    // Tests for population accessors
    // ============================================================================

    #[test]
    fn test_population_returns_current_population() {
        let sim = create_test_simulation();
        let pop = sim.population();

        assert_eq!(pop.size(), 10);
        assert_eq!(pop.generation(), 0);
    }

    #[test]
    fn test_population_mut_allows_modification() {
        let mut sim = create_test_simulation();

        // Get mutable reference and verify we can access it
        let pop_mut = sim.population_mut();
        assert_eq!(pop_mut.size(), 10);

        // Verify the changes are reflected
        assert_eq!(sim.population().size(), 10);
    }

    #[test]
    fn test_generation_tracks_correctly() {
        let mut sim = create_test_simulation();

        assert_eq!(sim.generation(), 0);

        sim.step().unwrap();
        assert_eq!(sim.generation(), 1);

        sim.step().unwrap();
        assert_eq!(sim.generation(), 2);

        sim.run_for(3).unwrap();
        assert_eq!(sim.generation(), 5);
    }

    // ============================================================================
    // Tests for configuration getters
    // ============================================================================

    #[test]
    fn test_simulation_config_getter() {
        let sim = create_test_simulation();
        let config = sim.simulation_config();

        assert_eq!(config.population_size, 10);
        assert_eq!(config.total_generations, 5);
        assert_eq!(config.seed, Some(42));
    }

    #[test]
    fn test_structure_config_getter() {
        let sim = create_test_simulation();
        let structure = sim.structure_config();

        assert!(structure.is_some());
        let s = structure.unwrap();
        assert_eq!(s.ru_length, 10);
        assert_eq!(s.rus_per_hor, 5);
        assert_eq!(s.hors_per_chr, 10);
    }

    #[test]
    fn test_mutation_config_getter() {
        let sim = create_test_simulation();
        let mutation = sim.mutation_config();

        assert!(matches!(
            mutation.substitution,
            SubstitutionModel::Uniform(_)
        ));
    }

    #[test]
    fn test_recombination_config_getter() {
        let sim = create_test_simulation();
        let _recombination = sim.recombination_config();

        // Just verify we can access the config
        // RecombinationConfig has a params field which is a RecombinationModel
    }

    #[test]
    fn test_fitness_config_getter() {
        let sim = create_test_simulation();
        let fitness = sim.fitness_config();

        // Test simulation uses neutral fitness
        assert!(fitness.is_neutral());
    }

    // ============================================================================
    // Tests for RNG state serialization
    // ============================================================================

    #[test]
    fn test_rng_state_bytes_can_be_serialized() {
        let sim = create_test_simulation();
        let state = sim.rng_state_bytes();

        // Should return non-empty bytes
        assert!(!state.is_empty());
    }

    #[test]
    fn test_set_rng_from_bytes_golden_path() {
        let mut sim1 = create_test_simulation();
        let mut sim2 = create_test_simulation();

        // Capture sim1's initial RNG state
        let state = sim1.rng_state_bytes();

        // Set sim2's RNG to match sim1
        sim2.set_rng_from_bytes(&state).unwrap();

        // Now both should produce same random sequences when stepped
        sim1.step().unwrap();
        sim2.step().unwrap();

        // Check that populations evolved similarly (same size and generation)
        assert_eq!(sim1.population().size(), sim2.population().size());
        assert_eq!(sim1.generation(), sim2.generation());
    }

    #[test]
    fn test_set_rng_from_bytes_invalid_data_fails() {
        let mut sim = create_test_simulation();

        // Try to set RNG state with invalid bytes
        let invalid_bytes = vec![1, 2, 3, 4, 5];
        let result = sim.set_rng_from_bytes(&invalid_bytes);

        assert!(result.is_err());
        assert!(result.unwrap_err().contains("Failed to deserialize"));
    }

    #[test]
    fn test_set_rng_from_bytes_empty_data_fails() {
        let mut sim = create_test_simulation();

        let empty_bytes = vec![];
        let result = sim.set_rng_from_bytes(&empty_bytes);

        assert!(result.is_err());
    }

    // ============================================================================
    // Tests for step
    // ============================================================================

    #[test]
    fn test_step_advances_one_generation() {
        let mut sim = create_test_simulation();

        assert_eq!(sim.generation(), 0);

        sim.step().unwrap();

        assert_eq!(sim.generation(), 1);
        assert_eq!(sim.population().size(), 10);
    }

    #[test]
    fn test_step_maintains_population_size() {
        let mut sim = create_test_simulation();
        let initial_size = sim.population().size();

        for _ in 0..10 {
            sim.step().unwrap();
            assert_eq!(sim.population().size(), initial_size);
        }
    }

    #[test]
    fn test_step_increments_generation() {
        let mut sim = create_test_simulation();

        for i in 0..5 {
            assert_eq!(sim.generation(), i);
            sim.step().unwrap();
            assert_eq!(sim.generation(), i + 1);
        }
    }

    #[test]
    fn test_step_produces_new_individuals() {
        let mut sim = create_test_simulation();

        let gen0_id = sim.population().get(0).unwrap().id().to_string();

        sim.step().unwrap();

        let gen1_id = sim.population().get(0).unwrap().id().to_string();

        // IDs should be different (gen0 vs gen1)
        assert_ne!(gen0_id, gen1_id);
        assert!(gen1_id.contains("gen1"));
    }

    #[test]
    fn test_step_with_mutation_changes_sequences() {
        let mut sim = SimulationBuilder::new()
            .population_size(10)
            .generations(5)
            .repeat_structure(10, 10, 10)
            .init_uniform(Nucleotide::A)
            .mutation_rate(0.1) // High mutation rate
            .recombination(0.0, 0.0, 0.0)
            .seed(42)
            .build()
            .unwrap();

        let gen0_seq = sim
            .population()
            .get(0)
            .unwrap()
            .haplotype1()
            .get(0)
            .unwrap()
            .sequence()
            .to_string();

        sim.step().unwrap();

        let gen1_seq = sim
            .population()
            .get(0)
            .unwrap()
            .haplotype1()
            .get(0)
            .unwrap()
            .sequence()
            .to_string();

        // With high mutation rate, sequences should differ
        assert_ne!(gen0_seq, gen1_seq);
    }

    #[test]
    fn test_step_without_mutation_with_recombination() {
        let mut sim = SimulationBuilder::new()
            .population_size(10)
            .generations(5)
            .repeat_structure(10, 10, 10)
            .mutation_rate(0.0)
            .recombination(0.5, 0.7, 0.3) // High recombination
            .seed(42)
            .build()
            .unwrap();

        // Should complete without errors even with high recombination
        sim.step().unwrap();

        assert_eq!(sim.generation(), 1);
        assert_eq!(sim.population().size(), 10);
    }

    #[test]
    fn test_step_computes_offspring_fitness() {
        use crate::evolution::GCContentFitness;

        let structure = UniformRepeatStructure::new(Nucleotide::A, 5, 10, 20, 1);
        let seq_config = InitializationConfig::Generate {
            structure: structure.clone(),
            mode: GenerationMode::Random,
        };
        let mutation = MutationConfig::uniform(0.01).unwrap();
        let recombination = RecombinationConfig::standard(0.01, 0.7, 0.1).unwrap();

        let gc_fitness = GCContentFitness::new(0.5, 1.0).unwrap();
        let fitness = FitnessConfig::new(Some(gc_fitness), None, None, None);

        let execution = ExecutionConfig::new(10, 5, Some(42));

        let sim_config = Configuration {
            execution,
            evolution: EvolutionConfig {
                mutation,
                recombination,
                fitness,
            },
            initialization: seq_config,
        };
        let mut sim = Simulation::new(sim_config.clone()).unwrap();

        sim.step().unwrap();

        // All offspring should have fitness computed
        for ind in sim.population().individuals() {
            assert!(*ind.cached_fitness().unwrap_or_default() > 0.0);
        }
    }

    // ============================================================================
    // Tests for run and run_for
    // ============================================================================

    #[test]
    fn test_run_for_advances_specified_generations() {
        let mut sim = create_test_simulation();

        assert_eq!(sim.generation(), 0);

        sim.run_for(3).unwrap();

        assert_eq!(sim.generation(), 3);
        assert_eq!(sim.population().size(), 10);
    }

    #[test]
    fn test_run_advances_configured_generations() {
        let mut sim = create_test_simulation();

        sim.run().unwrap();

        // Should run for configured number of generations (5)
        assert_eq!(sim.generation(), 5);
    }

    #[test]
    fn test_run_for_zero_generations() {
        let mut sim = create_test_simulation();

        sim.run_for(0).unwrap();

        assert_eq!(sim.generation(), 0);
        assert_eq!(sim.population().size(), 10);
    }

    #[test]
    fn test_run_for_multiple_generations() {
        let mut sim = create_test_simulation();

        sim.run_for(10).unwrap();

        assert_eq!(sim.generation(), 10);
        assert_eq!(sim.population().size(), 10);
    }

    #[test]
    fn test_run_for_is_reproducible_with_same_seed() {
        let mut sim1 = SimulationBuilder::new()
            .population_size(5)
            .generations(10)
            .repeat_structure(5, 5, 5)
            .mutation_rate(0.01)
            .recombination(0.01, 0.7, 0.1)
            .seed(777)
            .build()
            .unwrap();

        let mut sim2 = SimulationBuilder::new()
            .population_size(5)
            .generations(10)
            .repeat_structure(5, 5, 5)
            .mutation_rate(0.01)
            .recombination(0.01, 0.7, 0.1)
            .seed(777)
            .build()
            .unwrap();

        sim1.run_for(5).unwrap();
        sim2.run_for(5).unwrap();

        // Both should be at generation 5
        assert_eq!(sim1.generation(), sim2.generation());

        // Check that at least some sequences match
        let seq1 = sim1
            .population()
            .get(0)
            .unwrap()
            .haplotype1()
            .get(0)
            .unwrap()
            .sequence()
            .to_string();

        let seq2 = sim2
            .population()
            .get(0)
            .unwrap()
            .haplotype1()
            .get(0)
            .unwrap()
            .sequence()
            .to_string();

        assert_eq!(seq1, seq2, "Same seed should produce same evolution");
    }

    #[test]
    fn test_run_uses_configured_generations() {
        let mut sim = SimulationBuilder::new()
            .population_size(10)
            .generations(7) // Config specifies 7 generations
            .repeat_structure(5, 5, 5)
            .mutation_rate(0.001)
            .recombination(0.01, 0.7, 0.1)
            .seed(42)
            .build()
            .unwrap();

        sim.run().unwrap();

        assert_eq!(sim.generation(), 7);
    }

    #[test]
    fn test_run_with_zero_configured_generations() {
        let mut sim = SimulationBuilder::new()
            .population_size(10)
            .generations(0) // No generations
            .repeat_structure(5, 5, 5)
            .mutation_rate(0.001)
            .recombination(0.01, 0.7, 0.1)
            .seed(42)
            .build()
            .unwrap();

        sim.run().unwrap();

        assert_eq!(sim.generation(), 0);
    }

    #[test]
    fn test_run_then_run_for_continues_correctly() {
        let mut sim = SimulationBuilder::new()
            .population_size(10)
            .generations(3)
            .repeat_structure(5, 5, 5)
            .mutation_rate(0.001)
            .recombination(0.01, 0.7, 0.1)
            .seed(42)
            .build()
            .unwrap();

        sim.run().unwrap();
        assert_eq!(sim.generation(), 3);

        sim.run_for(2).unwrap();
        assert_eq!(sim.generation(), 5);
    }

    // ============================================================================
    // Tests for from_checkpoint
    // ============================================================================
    // Note: Comprehensive checkpoint tests are in tests/checkpoint_resume.rs
    // These are integration tests that use AsyncRecorder for storage.

    #[test]
    fn test_from_checkpoint_nonexistent_file() {
        let result = Simulation::from_checkpoint("nonexistent.db", "test_sim");

        assert!(result.is_err());
        // The actual error message varies, just check it's an error
    }
}
