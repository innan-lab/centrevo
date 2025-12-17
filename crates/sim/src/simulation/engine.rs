//! Simulation engine for evolutionary processes.
//!
//! This module provides the main simulation loop that orchestrates mutation,
//! recombination, selection, and reproduction across generations.

use crate::errors::RecombinationError;
use crate::evolution::RecombinationType;

use crate::genome::{Chromosome, Individual};
use crate::simulation::{
    FitnessConfig, GenerationMode, MutationConfig, Population, RecombinationConfig, SequenceConfig,
    SimulationConfig, UniformRepeatStructure,
};
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;

/// Main simulation engine.
#[derive(Debug)]
pub struct Simulation {
    /// Current population
    population: Population,
    /// Repeat structure parameters
    /// Repeat structure parameters (optional, for metadata/reproducibility)
    #[allow(dead_code)]
    structure: Option<UniformRepeatStructure>,
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
    #[deprecated(
        since = "0.1.0",
        note = "Use SimulationBuilder or Simulation::from_sequences instead"
    )]
    pub fn new(
        structure: UniformRepeatStructure,
        mutation: MutationConfig,
        recombination: RecombinationConfig,
        fitness: FitnessConfig,
        config: SimulationConfig,
    ) -> Result<Self, String> {
        // Create initial population with uniform sequences
        let seq_config = SequenceConfig::Generate {
            structure: structure.clone(),
            mode: GenerationMode::Uniform,
        };
        Self::from_config(
            seq_config,
            Some(structure),
            mutation,
            recombination,
            fitness,
            config,
        )
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

    /// Create a new simulation from configuration.
    ///
    /// This allows initializing a simulation with sequences from FASTA files,
    /// JSON data, or a previous simulation database, or generating them.
    ///
    /// # Arguments
    ///
    /// * `seq_config` - Sequence initialization configuration
    /// * `structure` - Repeat structure parameters (optional, for metadata)
    /// * `mutation` - Mutation configuration
    /// * `recombination` - Recombination configuration
    /// * `fitness` - Fitness configuration
    /// * `config` - Simulation configuration
    ///
    /// # Returns
    ///
    /// A `Simulation` instance initialized with the provided sequences.
    pub fn from_config(
        seq_config: SequenceConfig,
        structure: Option<UniformRepeatStructure>,
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

        // Initialize individuals from config
        let mut individuals =
            crate::simulation::initialize(&seq_config, config.population_size, &mut rng)
                .map_err(|e| format!("Failed to initialize from sequences: {e}"))?;

        // Compute fitness for initial population
        // This ensures the invariant that all individuals have cached fitness
        individuals.par_iter_mut().for_each(|ind| {
            let mut val = crate::base::FitnessValue::default();

            if let Some(gc) = &fitness.gc_content {
                use crate::evolution::IndividualFitness;
                val *= gc.individual_fitness(ind);
            }
            if let Some(len) = &fitness.length {
                use crate::evolution::IndividualFitness;
                val *= len.individual_fitness(ind);
            }
            if let Some(sim) = &fitness.seq_similarity {
                use crate::evolution::IndividualFitness;
                val *= sim.individual_fitness(ind);
            }
            if let Some(len_sim) = &fitness.length_similarity {
                use crate::evolution::IndividualFitness;
                val *= len_sim.individual_fitness(ind);
            }

            ind.set_cached_fitness(val);
        });

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
            .map(|snap| snap.to_individual(&snapshot.config.codec))
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
        let mut population = Population::new(
            format!("pop_{gen}", gen = checkpoint.generation),
            individuals,
        );

        // Recompute fitness for restored population
        // This is critical because:
        // 1. Older checkpoints might not have cached fitness
        // 2. We want to ensure fitness matches the current code/config
        // 3. Our invariant requires all individuals to have cached fitness
        let fitness_config = &snapshot.fitness;
        population.individuals_mut().par_iter_mut().for_each(|ind| {
            let mut val = crate::base::FitnessValue::default();

            if let Some(gc) = &fitness_config.gc_content {
                use crate::evolution::IndividualFitness;
                val *= gc.individual_fitness(ind);
            }
            if let Some(len) = &fitness_config.length {
                use crate::evolution::IndividualFitness;
                val *= len.individual_fitness(ind);
            }
            if let Some(sim) = &fitness_config.seq_similarity {
                use crate::evolution::IndividualFitness;
                val *= sim.individual_fitness(ind);
            }
            if let Some(len_sim) = &fitness_config.length_similarity {
                use crate::evolution::IndividualFitness;
                val *= len_sim.individual_fitness(ind);
            }

            ind.set_cached_fitness(val);
        });

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
            structure: Some(snapshot.structure), // Assuming snapshot still has it, will fix storage later
            mutation: snapshot.mutation,
            recombination: snapshot.recombination,
            fitness: snapshot.fitness,
            config: snapshot.config,
            rng,
        })
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

    /// Get reference to repeat structure (if available).
    pub fn structure(&self) -> Option<&UniformRepeatStructure> {
        self.structure.as_ref()
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
                let mut events = self.recombination.params.sample_events(chr1, chr2, rng);

                // Process events from right to left (descending position)
                events.reverse();

                // Apply events
                for event in events {
                    Simulation::apply_event_to_pair(&self.recombination.params, chr1, chr2, event)?;
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
            self.mutation.substitution.mutate_sequence_sparse(seq, rng);

            // Apply indels if configured
            if let Some(indel_model) = &self.mutation.indel {
                indel_model.apply_indels(seq, rng);
            }
        }

        // Clear cached fitness since the sequence has changed
        gamete.set_cached_fitness(crate::base::FitnessValue::default());

        Ok(gamete)
    }

    /// Select parents and generate offspring for the next generation.
    fn generate_offspring(&mut self) -> Result<Vec<Individual>, String> {
        // Select parent pairs (using cached fitness from current population)
        let pairs = self
            .population
            .select_parents(&mut self.rng, self.config.population_size);

        // Pre-allocate generation string
        let gen_str = format!("ind_gen{}_", self.generation() + 1);

        // Generate seeds for each offspring
        let seeds: Vec<u64> = (0..self.config.population_size)
            .map(|_| self.rng.random())
            .collect();

        // Get immutable reference to population for parallel access
        let population = &self.population;
        let fitness_config = &self.fitness;

        // Generate offspring from parent pairs in parallel
        let offspring: Result<Vec<Individual>, String> = pairs
            .par_iter()
            .zip(seeds.par_iter())
            .enumerate()
            .map(|(i, ((parent1_idx, parent2_idx), &seed))| {
                let mut local_rng = Xoshiro256PlusPlus::seed_from_u64(seed);

                let parent1 = population.get(*parent1_idx).unwrap();
                let parent2 = population.get(*parent2_idx).unwrap();

                // Generate two gametes, one from each parent
                let gamete1 = self.produce_gamete(parent1, &mut local_rng)?;
                let gamete2 = self.produce_gamete(parent2, &mut local_rng)?;

                // Combine to form new individual
                let id = format!("{gen_str}{i}");
                let mut new_ind = Individual::new(id, gamete1, gamete2);

                // Compute intrinsic fitness immediately (born with it)
                let mut fitness = crate::base::FitnessValue::default();

                // 1. GC Content
                if let Some(gc) = &fitness_config.gc_content {
                    use crate::evolution::IndividualFitness;
                    fitness *= gc.individual_fitness(&new_ind);
                }
                // 2. Length
                if let Some(len_fit) = &fitness_config.length {
                    use crate::evolution::IndividualFitness;
                    fitness *= len_fit.individual_fitness(&new_ind);
                }
                // 3. Sequence Similarity
                if let Some(sim) = &fitness_config.seq_similarity {
                    use crate::evolution::IndividualFitness;
                    fitness *= sim.individual_fitness(&new_ind);
                }
                // 4. Length Similarity
                if let Some(len_sim) = &fitness_config.length_similarity {
                    use crate::evolution::IndividualFitness;
                    fitness *= len_sim.individual_fitness(&new_ind);
                }

                new_ind.set_cached_fitness(fitness);

                Ok(new_ind)
            })
            .collect();

        offspring
    }

    /// Advance simulation by one generation.
    pub fn step(&mut self) -> Result<(), String> {
        // 1. Generate offspring (Selection -> Meiosis/Mutation -> Reproduction)
        // This includes fitness calculation inside generate_offspring
        let offspring = self.generate_offspring()?;

        // 2. Replace population with offspring
        self.population.set_individuals(offspring);
        self.population.increment_generation();

        // 3. (Removed) Update fitness for new population
        // Fitness is now intrinsic and computed at birth (in generate_offspring)

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
    use crate::simulation::SimulationBuilder;

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

    #[test]
    fn test_simulation_new() {
        let sim = create_test_simulation();

        assert_eq!(sim.population().size(), 10);
        assert_eq!(sim.generation(), 0);
    }

    #[test]
    fn test_simulation_initial_population() {
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
    fn test_simulation_step() {
        let mut sim = create_test_simulation();

        assert_eq!(sim.generation(), 0);

        sim.step().unwrap();

        assert_eq!(sim.generation(), 1);
        assert_eq!(sim.population().size(), 10);
    }

    #[test]
    fn test_mutation_config_uniform() {
        let config = MutationConfig::uniform(0.001).unwrap();
        // Just check it was created successfully
        assert!(matches!(config.substitution, SubstitutionModel::Uniform(_)));
    }

    #[test]
    fn test_simulation_run_for() {
        let mut sim = create_test_simulation();

        assert_eq!(sim.generation(), 0);

        sim.run_for(3).unwrap();

        assert_eq!(sim.generation(), 3);
        assert_eq!(sim.population().size(), 10);
    }

    #[test]
    fn test_simulation_run() {
        let mut sim = create_test_simulation();

        sim.run().unwrap();

        // Should run for configured number of generations (5)
        assert_eq!(sim.generation(), 5);
    }

    #[test]
    fn test_create_uniform_individual() {
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
        let params = sim.recombination.params.clone();
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

        let params = sim.recombination.params.clone();
        Simulation::apply_event_to_pair(&params, &mut donor, &mut recipient, event).unwrap();

        // No changes to recipient - still TTTT
        assert_eq!(recipient.to_string(), "TTTT");
    }
}
