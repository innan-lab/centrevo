//! Simulation engine for evolutionary processes.
//!
//! This module provides the main simulation loop that orchestrates mutation,
//! recombination, selection, and reproduction across generations.

use crate::base::Sequence;
use crate::genome::repeat_map::RepeatMap;
use crate::genome::{Chromosome, Haplotype, Individual};
use crate::simulation::{
    FitnessConfig, MutationConfig, Population, RecombinationConfig, RepeatStructure,
    SimulationConfig,
};
use crate::evolution::RecombinationType;
use crate::errors::RecombinationError;
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;
use std::sync::Arc;

/// Main simulation engine.
#[derive(Debug)]
pub struct Simulation {
    /// Current population
    population: Population,
    /// Repeat structure parameters
    #[allow(dead_code)] // Will be used in future iterations
    structure: RepeatStructure,
    /// Mutation configuration
    mutation: MutationConfig,
    /// Recombination configuration
    recombination: RecombinationConfig,
    /// Fitness configuration
    fitness: FitnessConfig,
    /// Simulation configuration
    config: SimulationConfig,
    /// Random number generator (using Xoshiro256++ for better performance)
    rng: Xoshiro256PlusPlus,
}

impl Simulation {
    /// Create a new simulation with uniform initial sequences.
    pub fn new(
        structure: RepeatStructure,
        mutation: MutationConfig,
        recombination: RecombinationConfig,
        fitness: FitnessConfig,
        config: SimulationConfig,
    ) -> Result<Self, String> {
        // Create RNG from seed or thread_rng
        // Using Xoshiro256++ which is 2-3x faster than StdRng
        let rng = if let Some(seed) = config.seed {
            Xoshiro256PlusPlus::seed_from_u64(seed)
        } else {
            Xoshiro256PlusPlus::from_seed(rand::rng().random())
        };

        // Create initial population with uniform sequences
        let individuals = Self::create_initial_population(&structure, config.population_size)?;
        let population = Population::new("pop0", individuals);

        Ok(Self {
            population,
            structure,
            mutation,
            recombination,
            fitness,
            config,
            rng,
        })
    }

    /// Create a new simulation with random initial sequences.
    ///
    /// Each position in each sequence will be initialized with a random base
    /// from the alphabet. This is useful for starting simulations from a
    /// completely random state.
    ///
    /// # Arguments
    ///
    /// * `structure` - Repeat structure parameters
    /// * `mutation` - Mutation configuration
    /// * `recombination` - Recombination configuration
    /// * `fitness` - Fitness configuration
    /// * `config` - Simulation configuration
    ///
    /// # Returns
    ///
    /// A `Simulation` instance initialized with random sequences.
    pub fn new_random(
        structure: RepeatStructure,
        mutation: MutationConfig,
        recombination: RecombinationConfig,
        fitness: FitnessConfig,
        config: SimulationConfig,
    ) -> Result<Self, String> {
        // Create RNG from seed or thread_rng
        let mut rng = if let Some(seed) = config.seed {
            Xoshiro256PlusPlus::seed_from_u64(seed)
        } else {
            Xoshiro256PlusPlus::from_seed(rand::rng().random())
        };

        // Create initial population with random sequences
        let individuals =
            Self::create_random_population(&structure, config.population_size, &mut rng)?;
        let population = Population::new("pop0", individuals);

        Ok(Self {
            population,
            structure,
            mutation,
            recombination,
            fitness,
            config,
            rng,
        })
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
                let (new1, new2) = params.crossover(chr1, chr2, pos1, pos2)
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
                            let clamped_len = length
                                .min(recipient_remaining)
                                .min(donor_remaining);
                            if clamped_len == 0 {
                                return Ok(());
                            }

                            if let Ok(new2) = params.gene_conversion(
                                chr2, chr1, recipient_start, donor_start, clamped_len,
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

    /// Create a new simulation from custom sequences.
    ///
    /// This allows initializing a simulation with sequences from FASTA files,
    /// JSON data, or a previous simulation database. The sequences are validated
    /// to ensure they match the expected parameters.
    ///
    /// # Arguments
    ///
    /// * `source` - Source of initial sequences (FASTA, JSON, or Database)
    /// * `structure` - Repeat structure parameters
    /// * `mutation` - Mutation configuration
    /// * `recombination` - Recombination configuration
    /// * `fitness` - Fitness configuration
    /// * `config` - Simulation configuration
    ///
    /// # Returns
    ///
    /// A `Simulation` instance initialized with the provided sequences, or an error
    /// if the sequences are invalid or don't match the parameters.
    pub fn from_sequences(
        source: crate::simulation::initialization::SequenceInput,
        structure: RepeatStructure,
        mutation: MutationConfig,
        recombination: RecombinationConfig,
        fitness: FitnessConfig,
        config: SimulationConfig,
    ) -> Result<Self, String> {
        // Create RNG from seed or thread_rng
        let rng = if let Some(seed) = config.seed {
            Xoshiro256PlusPlus::seed_from_u64(seed)
        } else {
            Xoshiro256PlusPlus::from_seed(rand::rng().random())
        };

        // Initialize individuals from source
        let individuals = crate::simulation::initialization::initialize_from_source(
            source,
            &structure,
            config.population_size,
        )
        .map_err(|e| format!("Failed to initialize from sequences: {e}"))?;

        let population = Population::new("pop0", individuals);

        Ok(Self {
            population,
            structure,
            mutation,
            recombination,
            fitness,
            config,
            rng,
        })
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
            QueryBuilder::new(db_path).map_err(|e| format!("Failed to open database: {e}"))?;

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
        let snapshot = query
            .get_full_config(sim_id)
            .map_err(|e| format!("Failed to load configuration: {e}"))?;

        // Load population state at checkpoint generation
        let snapshots = query
            .get_generation(sim_id, checkpoint.generation)
            .map_err(|e| format!("Failed to load population: {e}"))?;

        if snapshots.is_empty() {
            return Err(format!(
                "No population data found for generation {generation}",
                generation = checkpoint.generation
            ));
        }

        // Reconstruct individuals from snapshots
        let individuals: Result<Vec<_>, String> = snapshots
            .iter()
            .map(|snap| snap.to_individual(&snapshot.structure))
            .collect();
        let individuals = individuals?;

        // Validate population size matches configuration
        if individuals.len() != snapshot.config.population_size {
            return Err(format!(
                "Population size mismatch: config expects {expected}, checkpoint has {actual}",
                expected = snapshot.config.population_size,
                actual = individuals.len()
            ));
        }

        // Create population with correct generation number
        let mut population = Population::new(format!("pop_{gen}", gen = checkpoint.generation), individuals);

        // Set generation counter to checkpoint generation
        for _ in 0..checkpoint.generation {
            population.increment_generation();
        }

        // Restore RNG state
        let rng = bincode::deserialize(&checkpoint.rng_state)
            .map_err(|e| format!("Failed to restore RNG state: {e}"))?;

        // Close query builder
        query
            .close()
            .map_err(|e| format!("Failed to close database: {e}"))?;

        Ok(Self {
            population,
            structure: snapshot.structure,
            mutation: snapshot.mutation,
            recombination: snapshot.recombination,
            fitness: snapshot.fitness,
            config: snapshot.config,
            rng,
        })
    }

    /// Create initial population with uniform sequences.
    fn create_initial_population(
        structure: &RepeatStructure,
        pop_size: usize,
    ) -> Result<Vec<Individual>, String> {
        let mut individuals = Vec::with_capacity(pop_size);

        for i in 0..pop_size {
            let ind = Self::create_uniform_individual(format!("ind_{i}"), structure)?;
            individuals.push(ind);
        }

        Ok(individuals)
    }

    /// Create initial population with random sequences.
    fn create_random_population(
        structure: &RepeatStructure,
        pop_size: usize,
        rng: &mut impl Rng,
    ) -> Result<Vec<Individual>, String> {
        let mut individuals = Vec::with_capacity(pop_size);

        for i in 0..pop_size {
            let ind = Self::create_random_individual(format!("ind_{i}"), structure, rng)?;
            individuals.push(ind);
        }

        Ok(individuals)
    }

    /// Create a single individual with uniform sequences.
    fn create_uniform_individual(
        id: impl Into<Arc<str>>,
        structure: &RepeatStructure,
    ) -> Result<Individual, String> {
        // Create haplotypes with uniform chromosomes
        let mut hap1 = Haplotype::with_capacity(structure.chrs_per_hap);
        let mut hap2 = Haplotype::with_capacity(structure.chrs_per_hap);

        for chr_idx in 0..structure.chrs_per_hap {
            // Create chromosomes
            let chr1 = Chromosome::uniform(
                format!("chr{chr_idx}"),
                structure.init_base,
                structure.ru_length,
                structure.rus_per_hor,
                structure.hors_per_chr,
            );
            let chr2 = Chromosome::uniform(
                format!("chr{chr_idx}"),
                structure.init_base,
                structure.ru_length,
                structure.rus_per_hor,
                structure.hors_per_chr,
            );

            hap1.push(chr1);
            hap2.push(chr2);
        }

        Ok(Individual::new(id, hap1, hap2))
    }

    /// Create a single individual with random sequences.
    fn create_random_individual(
        id: impl Into<Arc<str>>,
        structure: &RepeatStructure,
        rng: &mut impl Rng,
    ) -> Result<Individual, String> {
        use crate::base::Nucleotide;

        let chr_length = structure.chr_length();
        let alphabet_size = 4; // DNA has 4 bases

        // Create haplotypes with random chromosomes
        let mut hap1 = Haplotype::with_capacity(structure.chrs_per_hap);
        let mut hap2 = Haplotype::with_capacity(structure.chrs_per_hap);

        for chr_idx in 0..structure.chrs_per_hap {
            // Create random sequence for haplotype 1
            let mut seq1 = Sequence::with_capacity(chr_length);
            for _ in 0..chr_length {
                let idx = rng.random_range(0..alphabet_size);
                let base = Nucleotide::from_index(idx as u8)
                    .ok_or_else(|| format!("Invalid nucleotide index: {idx}"))?;
                seq1.push(base);
            }

            // Create random sequence for haplotype 2
            let mut seq2 = Sequence::with_capacity(chr_length);
            for _ in 0..chr_length {
                let idx = rng.random_range(0..alphabet_size);
                let base = Nucleotide::from_index(idx as u8)
                    .ok_or_else(|| format!("Invalid nucleotide index: {idx}"))?;
                seq2.push(base);
            }

            // Create map
            let map = RepeatMap::uniform(
                structure.ru_length,
                structure.rus_per_hor,
                structure.hors_per_chr,
            );

            // Create chromosomes
            let chr1 = Chromosome::new(format!("chr{chr_idx}"), seq1, map.clone());
            let chr2 = Chromosome::new(format!("chr{chr_idx}"), seq2, map);

            hap1.push(chr1);
            hap2.push(chr2);
        }

        Ok(Individual::new(id, hap1, hap2))
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

    /// Get reference to simulation configuration.
    pub fn config(&self) -> &SimulationConfig {
        &self.config
    }

    /// Get reference to repeat structure.
    pub fn structure(&self) -> &RepeatStructure {
        &self.structure
    }

    /// Get reference to mutation configuration.
    pub fn mutation(&self) -> &MutationConfig {
        &self.mutation
    }

    /// Get reference to recombination configuration.
    pub fn recombination(&self) -> &RecombinationConfig {
        &self.recombination
    }

    /// Get reference to fitness configuration.
    pub fn fitness(&self) -> &FitnessConfig {
        &self.fitness
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

    /// Apply mutation to all individuals in the population.
    fn apply_mutation(&mut self) -> Result<(), String> {
        let pop_size = self.population.size();

        // Generate seeds for each individual (using faster RNG)
        let seeds: Vec<u64> = (0..pop_size).map(|_| self.rng.random()).collect();

        // Parallel mutation with independent RNGs per individual
        self.population
            .individuals_mut()
            .par_iter_mut()
            .zip(seeds.par_iter())
            .for_each(|(individual, &seed)| {
                // Use Xoshiro256++ for thread-local RNG (faster than StdRng)
                let mut local_rng = Xoshiro256PlusPlus::seed_from_u64(seed);

                // Mutate first haplotype using Poisson pre-sampling
                for chr in individual.haplotype1_mut().chromosomes_mut() {
                    let seq = chr.sequence_mut();
                    // Apply substitutions
                    self.mutation
                        .substitution
                        .mutate_sequence_poisson(seq, &mut local_rng);

                    // Apply indels if configured
                    if let Some(indel_model) = &self.mutation.indel {
                        indel_model.apply_indels(seq, &mut local_rng);
                    }
                }

                // Mutate second haplotype using Poisson pre-sampling
                for chr in individual.haplotype2_mut().chromosomes_mut() {
                    let seq = chr.sequence_mut();
                    // Apply substitutions
                    self.mutation
                        .substitution
                        .mutate_sequence_poisson(seq, &mut local_rng);

                    // Apply indels if configured
                    if let Some(indel_model) = &self.mutation.indel {
                        indel_model.apply_indels(seq, &mut local_rng);
                    }
                }
            });

        Ok(())
    }

    /// Apply recombination to all individuals in the population.
    ///
    /// Note: Some recombination events may become invalid while being applied
    /// (for example, if previous events change chromosome lengths). For gene
    /// conversion events where the specified range becomes invalid, the engine
    /// attempts to clamp the tract length to the available space and retry the
    /// operation. If clamping cannot resolve the issue, the event is skipped.
    fn apply_recombination(&mut self) -> Result<(), String> {
        let pop_size = self.population.size();

        // Generate seeds for each individual
        let seeds: Vec<u64> = (0..pop_size).map(|_| self.rng.random()).collect();

        // Clone params to avoid immutable borrow of `self` in the parallel closure
        let params = self.recombination.params.clone();

        // Parallel recombination with independent RNGs per individual
        self.population
            .individuals_mut()
            .par_iter_mut()
            .zip(seeds.par_iter())
            .try_for_each(|(individual, &seed)| -> Result<(), String> {
                let mut local_rng = Xoshiro256PlusPlus::seed_from_u64(seed);
                let (hap1, hap2) = individual.haplotypes_mut();

                // Recombine corresponding chromosomes
                for chr_idx in 0..hap1.len().min(hap2.len()) {
                    if let (Some(chr1), Some(chr2)) = (hap1.get_mut(chr_idx), hap2.get_mut(chr_idx))
                    {
                        // Sample recombination events
                        let mut events =
                            self.recombination
                                .params
                                .sample_events(chr1, chr2, &mut local_rng);

                        // Process events from right to left (descending position)
                        // This ensures that changes to the tail (right side) do not invalidate
                        // indices for subsequent events on the left side.
                        events.reverse();

                        // Apply the recombination events sequentially
                        for event in events {
                            Simulation::apply_event_to_pair(&params, chr1, chr2, event)?;
                        }
                    }
                }

                Ok(())
            })?;

        Ok(())
    }

    /// Select parents and generate offspring for the next generation.
    fn generate_offspring(&mut self) -> Result<Vec<Individual>, String> {
        // Compute fitness for current population
        let fitness_values = self.population.compute_fitness(&self.fitness);

        // Select parent pairs
        let pairs = self.population.select_parents(
            &mut self.rng,
            &fitness_values,
            self.config.population_size,
        );

        // Pre-allocate generation string to avoid repeated allocations
        let gen_str = format!("ind_gen{}_", self.generation() + 1);

        // Generate seeds for each offspring
        let seeds: Vec<u64> = (0..self.config.population_size)
            .map(|_| self.rng.random())
            .collect();

        // Get immutable reference to population for parallel access
        let population = &self.population;

        // Generate offspring from parent pairs in parallel
        let offspring: Vec<Individual> = pairs
            .par_iter()
            .zip(seeds.par_iter())
            .enumerate()
            .map(|(i, ((parent1_idx, parent2_idx), &seed))| {
                let mut local_rng = Xoshiro256PlusPlus::seed_from_u64(seed);

                let parent1 = population.get(*parent1_idx).unwrap();
                let parent2 = population.get(*parent2_idx).unwrap();

                // Each parent contributes one gamete (haplotype)
                // Randomly pick one haplotype from each parent
                // Use gen::<f64>() which is faster than random_bool
                let hap1 = if local_rng.random::<f64>() < 0.5 {
                    parent1.haplotype1().clone()
                } else {
                    parent1.haplotype2().clone()
                };

                let hap2 = if local_rng.random::<f64>() < 0.5 {
                    parent2.haplotype1().clone()
                } else {
                    parent2.haplotype2().clone()
                };

                // More efficient ID construction
                let id = format!("{gen_str}{i}");
                Individual::new(id, hap1, hap2)
            })
            .collect();

        Ok(offspring)
    }

    /// Advance simulation by one generation.
    pub fn step(&mut self) -> Result<(), String> {
        // 1. Apply mutation to current population
        self.apply_mutation()?;

        // 2. Apply recombination to current population
        self.apply_recombination()?;

        // 3. Generate offspring for next generation
        let offspring = self.generate_offspring()?;

        // 4. Replace population with offspring
        self.population.set_individuals(offspring);
        self.population.increment_generation();

        Ok(())
    }

    /// Run simulation for the configured number of generations.
    pub fn run(&mut self) -> Result<(), String> {
        for _ in 0..self.config.total_generations {
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

    fn create_test_config() -> (
        RepeatStructure,
        MutationConfig,
        RecombinationConfig,
        FitnessConfig,
        SimulationConfig,
    ) {
        let structure = RepeatStructure::new(Nucleotide::A, 10, 5, 10, 1);

        let mutation = MutationConfig::uniform(0.001).unwrap();
        let recombination = RecombinationConfig::standard(0.01, 0.7, 0.1).unwrap();
        let fitness = FitnessConfig::neutral();
        let config = SimulationConfig::new(10, 5, Some(42));

        (structure, mutation, recombination, fitness, config)
    }

    #[test]
    fn test_simulation_new() {
        let (structure, mutation, recombination, fitness, config) = create_test_config();

        let sim = Simulation::new(structure, mutation, recombination, fitness, config).unwrap();

        assert_eq!(sim.population().size(), 10);
        assert_eq!(sim.generation(), 0);
    }

    #[test]
    fn test_simulation_initial_population() {
        let (structure, mutation, recombination, fitness, config) = create_test_config();

        let sim = Simulation::new(structure, mutation, recombination, fitness, config).unwrap();

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
    fn test_simulation_step() {
        let (structure, mutation, recombination, fitness, config) = create_test_config();

        let mut sim = Simulation::new(structure, mutation, recombination, fitness, config).unwrap();

        assert_eq!(sim.generation(), 0);

        sim.step().unwrap();

        assert_eq!(sim.generation(), 1);
        assert_eq!(sim.population().size(), 10);
    }

    #[test]
    fn test_mutation_config_uniform() {
        let config = MutationConfig::uniform(0.001).unwrap();
        // Just check it was created successfully
        assert!(matches!(config.substitution, SubstitutionModel { .. }));
    }

    #[test]
    fn test_simulation_run_for() {
        let (structure, mutation, recombination, fitness, config) = create_test_config();

        let mut sim = Simulation::new(structure, mutation, recombination, fitness, config).unwrap();

        assert_eq!(sim.generation(), 0);

        sim.run_for(3).unwrap();

        assert_eq!(sim.generation(), 3);
        assert_eq!(sim.population().size(), 10);
    }

    #[test]
    fn test_simulation_run() {
        let (structure, mutation, recombination, fitness, config) = create_test_config();

        let mut sim = Simulation::new(structure, mutation, recombination, fitness, config).unwrap();

        sim.run().unwrap();

        // Should run for configured number of generations (5)
        assert_eq!(sim.generation(), 5);
    }

    #[test]
    fn test_create_uniform_individual() {
        let structure = RepeatStructure::new(Nucleotide::C, 10, 5, 10, 2);

        let ind = Simulation::create_uniform_individual("test_ind", &structure).unwrap();

        assert_eq!(ind.id(), "test_ind");
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
    fn test_apply_event_to_pair_gene_conversion_clamps() {
        // Create a small simulation with trivial configs
        let structure = RepeatStructure::new(Nucleotide::A, 1, 1, 6, 1); // chr len 6
        let mutation = MutationConfig::uniform(0.0).unwrap();
        let recombination = RecombinationConfig::standard(0.0, 0.5, 0.5).unwrap();
        let fitness = FitnessConfig::neutral();
        let config = SimulationConfig::new(1, 1, Some(42));

        let sim = Simulation::new(structure, mutation, recombination, fitness, config).unwrap();

        // Chron 1: donor (length 6), Chron 2: recipient (length 4)
        let mut donor = Chromosome::uniform("d", Nucleotide::A, 1, 1, 6);
        let mut recipient = Chromosome::uniform("r", Nucleotide::T, 1, 1, 4);

        // Create a gene conversion event where length exceeds recipient remaining length
        let event = crate::evolution::RecombinationType::GeneConversion {
            donor_start: 0, // donor start
            recipient_start: 1, // recipient start
            length: 5, // length > recipient remaining (3)
        };

        // Apply via helper — should clamp to 3 and succeed
        let params = sim.recombination.params.clone();
        Simulation::apply_event_to_pair(&params, &mut donor, &mut recipient, event).unwrap();

        assert_eq!(recipient.to_string(), "TAAA");
    }

    #[test]
    fn test_apply_event_to_pair_gene_conversion_skip_when_out_of_bounds() {
        let structure = RepeatStructure::new(Nucleotide::A, 1, 1, 4, 1);
        let mutation = MutationConfig::uniform(0.0).unwrap();
        let recombination = RecombinationConfig::standard(0.0, 0.5, 0.5).unwrap();
        let fitness = FitnessConfig::neutral();
        let config = SimulationConfig::new(1, 1, Some(42));

        let sim = Simulation::new(structure, mutation, recombination, fitness, config).unwrap();

        let mut donor = Chromosome::uniform("d", Nucleotide::A, 1, 1, 4);
        let mut recipient = Chromosome::uniform("r", Nucleotide::T, 1, 1, 4);

        // The recipient start is out of bounds
        let event = crate::evolution::RecombinationType::GeneConversion {
            donor_start: 0,
            recipient_start: 10, // invalid
            length: 2,
        };

        let params = sim.recombination.params.clone();
        Simulation::apply_event_to_pair(&params, &mut donor, &mut recipient, event).unwrap();

        // No changes to recipient - still TTTT
        assert_eq!(recipient.to_string(), "TTTT");
    }
}
