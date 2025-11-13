//! Simulation engine for evolutionary processes.
//!
//! This module provides the main simulation loop that orchestrates mutation,
//! recombination, selection, and reproduction across generations.

use crate::genome::{Chromosome, Haplotype, Individual};
use crate::simulation::{Population, RepeatStructure, MutationConfig, RecombinationConfig, FitnessConfig, SimulationConfig};
use crate::base::Sequence;
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;
use std::sync::Arc;
use rayon::prelude::*;

/// Main simulation engine.
#[derive(Debug)]
pub struct Simulation {
    /// Current population
    population: Population,
    /// Repeat structure parameters
    #[allow(dead_code)]  // Will be used in future iterations
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
        let individuals = Self::create_random_population(&structure, config.population_size, &mut rng)?;
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
        ).map_err(|e| format!("Failed to initialize from sequences: {e}"))?;

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
        let query = QueryBuilder::new(db_path)
            .map_err(|e| format!("Failed to open database: {e}"))?;
        
        // Get the latest checkpoint
        let checkpoint = query.get_latest_checkpoint(sim_id)
            .map_err(|e| format!("Failed to load checkpoint: {e}"))?;
        
        // Verify sim_id matches
        if checkpoint.sim_id != sim_id {
            return Err(format!(
                "Checkpoint sim_id mismatch: expected '{}', found '{}'",
                sim_id, checkpoint.sim_id
            ));
        }
        
        // Load complete configuration
        let snapshot = query.get_full_config(sim_id)
            .map_err(|e| format!("Failed to load configuration: {e}"))?;
        
        // Load population state at checkpoint generation
        let snapshots = query.get_generation(sim_id, checkpoint.generation)
            .map_err(|e| format!("Failed to load population: {e}"))?;
        
        if snapshots.is_empty() {
            return Err(format!(
                "No population data found for generation {}",
                checkpoint.generation
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
                "Population size mismatch: config expects {}, checkpoint has {}",
                snapshot.config.population_size,
                individuals.len()
            ));
        }
        
        // Create population with correct generation number
        let mut population = Population::new(format!("pop_{}", checkpoint.generation), individuals);
        
        // Set generation counter to checkpoint generation
        for _ in 0..checkpoint.generation {
            population.increment_generation();
        }
        
        // Restore RNG state
        let rng = bincode::deserialize(&checkpoint.rng_state)
            .map_err(|e| format!("Failed to restore RNG state: {e}"))?;
        
        // Close query builder
        query.close().map_err(|e| format!("Failed to close database: {e}"))?;
        
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
            let ind = Self::create_uniform_individual(
                format!("ind_{i}"),
                structure,
            )?;
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
            let ind = Self::create_random_individual(
                format!("ind_{i}"),
                structure,
                rng,
            )?;
            individuals.push(ind);
        }

        Ok(individuals)
    }

    /// Create a single individual with uniform sequences.
    fn create_uniform_individual(
        id: impl Into<Arc<str>>,
        structure: &RepeatStructure,
    ) -> Result<Individual, String> {
        let chr_length = structure.chr_length();

        // Create haplotypes with uniform chromosomes
        let mut hap1 = Haplotype::with_capacity(structure.chrs_per_hap);
        let mut hap2 = Haplotype::with_capacity(structure.chrs_per_hap);

        for chr_idx in 0..structure.chrs_per_hap {
            // Create uniform sequence
            let mut seq = Sequence::with_capacity(chr_length, structure.alphabet.clone());
            for _ in 0..chr_length {
                seq.push(structure.init_base);
            }

            // Create chromosomes
            let chr1 = Chromosome::new(
                format!("chr{chr_idx}"),
                seq.clone(),
                structure.ru_length,
                structure.rus_per_hor,
            );
            let chr2 = Chromosome::new(
                format!("chr{chr_idx}"),
                seq,
                structure.ru_length,
                structure.rus_per_hor,
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
        let alphabet_size = structure.alphabet.len();

        // Create haplotypes with random chromosomes
        let mut hap1 = Haplotype::with_capacity(structure.chrs_per_hap);
        let mut hap2 = Haplotype::with_capacity(structure.chrs_per_hap);

        for chr_idx in 0..structure.chrs_per_hap {
            // Create random sequence for haplotype 1
            let mut seq1 = Sequence::with_capacity(chr_length, structure.alphabet.clone());
            for _ in 0..chr_length {
                let idx = rng.random_range(0..alphabet_size);
                let ch = structure.alphabet.get_char(idx as u8)
                    .ok_or_else(|| format!("Invalid alphabet index: {idx}"))?;
                let base = Nucleotide::from_ascii(ch as u8)
                    .ok_or_else(|| format!("Invalid nucleotide character: {ch}"))?;
                seq1.push(base);
            }

            // Create random sequence for haplotype 2
            let mut seq2 = Sequence::with_capacity(chr_length, structure.alphabet.clone());
            for _ in 0..chr_length {
                let idx = rng.random_range(0..alphabet_size);
                let ch = structure.alphabet.get_char(idx as u8)
                    .ok_or_else(|| format!("Invalid alphabet index: {idx}"))?;
                let base = Nucleotide::from_ascii(ch as u8)
                    .ok_or_else(|| format!("Invalid nucleotide character: {ch}"))?;
                seq2.push(base);
            }

            // Create chromosomes
            let chr1 = Chromosome::new(
                format!("chr{chr_idx}"),
                seq1,
                structure.ru_length,
                structure.rus_per_hor,
            );
            let chr2 = Chromosome::new(
                format!("chr{chr_idx}"),
                seq2,
                structure.ru_length,
                structure.rus_per_hor,
            );

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
        bincode::serialize(&self.rng)
            .expect("Failed to serialize RNG state")
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
        let seeds: Vec<u64> = (0..pop_size)
            .map(|_| self.rng.random())
            .collect();
        
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
                    self.mutation.model.mutate_sequence_poisson(seq, &mut local_rng);
                }

                // Mutate second haplotype using Poisson pre-sampling
                for chr in individual.haplotype2_mut().chromosomes_mut() {
                    let seq = chr.sequence_mut();
                    self.mutation.model.mutate_sequence_poisson(seq, &mut local_rng);
                }
            });

        Ok(())
    }

    /// Apply recombination to all individuals in the population.
    fn apply_recombination(&mut self) -> Result<(), String> {
        let pop_size = self.population.size();
        
        // Generate seeds for each individual
        let seeds: Vec<u64> = (0..pop_size)
            .map(|_| self.rng.random())
            .collect();
        
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
                    if let (Some(chr1), Some(chr2)) = (hap1.get_mut(chr_idx), hap2.get_mut(chr_idx)) {
                        // Sample recombination event
                        let event = self.recombination.params.sample_event(
                            chr1.len(),
                            &mut local_rng,
                        );

                        // Apply the recombination event
                        match event {
                            crate::evolution::RecombinationType::None => {
                                // No recombination
                            }
                            crate::evolution::RecombinationType::Crossover { position } => {
                                // Perform crossover
                                let (new1, new2) = self.recombination.params.crossover(
                                    chr1.sequence(),
                                    chr2.sequence(),
                                    position,
                                ).map_err(|e| format!("Crossover failed: {e}"))?;
                                
                                *chr1.sequence_mut() = new1;
                                *chr2.sequence_mut() = new2;
                            }
                            crate::evolution::RecombinationType::GeneConversion { start, end } => {
                                // Perform gene conversion (chr1 -> chr2)
                                let new2 = self.recombination.params.gene_conversion(
                                    chr2.sequence(),
                                    chr1.sequence(),
                                    start,
                                    end,
                                ).map_err(|e| format!("Gene conversion failed: {e}"))?;
                                
                                *chr2.sequence_mut() = new2;
                            }
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
    use crate::base::{Nucleotide, Alphabet};

    fn create_test_config() -> (RepeatStructure, MutationConfig, RecombinationConfig, FitnessConfig, SimulationConfig) {
        let alphabet = Alphabet::dna();
        
        let structure = RepeatStructure::new(
            alphabet.clone(),
            Nucleotide::A,
            10,
            5,
            10,
            1,
        );

        let mutation = MutationConfig::uniform(alphabet, 0.001).unwrap();
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
        let alphabet = Alphabet::dna();
        let structure = RepeatStructure::new(
            alphabet,
            Nucleotide::C,
            10,
            5,
            10,
            2,
        );

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
}
