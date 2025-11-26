//! Population management and operations.
//!
//! This module provides structures and functions for managing populations
//! of individuals during evolutionary simulations.

use crate::base::{fitness, FitnessValue};
use crate::evolution::IndividualFitness;
use crate::genome::Individual;
use crate::simulation::FitnessConfig;
use rand::Rng;
use rayon::prelude::*;
use std::sync::Arc;

/// A population of diploid individuals.
#[derive(Debug, Clone)]
pub struct Population {
    /// The individuals in this population
    individuals: Vec<Individual>,
    /// Generation counter
    generation: usize,
    /// Population ID
    id: Arc<str>,
}

impl Population {
    /// Create a new population from individuals.
    pub fn new(id: impl Into<Arc<str>>, individuals: Vec<Individual>) -> Self {
        Self {
            individuals,
            generation: 0,
            id: id.into(),
        }
    }

    /// Get population ID.
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Get the current generation number.
    pub fn generation(&self) -> usize {
        self.generation
    }

    /// Increment the generation counter.
    pub fn increment_generation(&mut self) {
        self.generation += 1;
    }

    /// Get the number of individuals in the population.
    pub fn size(&self) -> usize {
        self.individuals.len()
    }

    /// Check if population is empty.
    pub fn is_empty(&self) -> bool {
        self.individuals.is_empty()
    }

    /// Get all individuals as a slice.
    pub fn individuals(&self) -> &[Individual] {
        &self.individuals
    }

    /// Get mutable access to individuals.
    pub fn individuals_mut(&mut self) -> &mut [Individual] {
        &mut self.individuals
    }

    /// Replace the entire population with new individuals.
    pub fn set_individuals(&mut self, individuals: Vec<Individual>) {
        self.individuals = individuals;
    }

    /// Get a specific individual by index.
    pub fn get(&self, index: usize) -> Option<&Individual> {
        self.individuals.get(index)
    }

    /// Get a mutable reference to a specific individual.
    pub fn get_mut(&mut self, index: usize) -> Option<&mut Individual> {
        self.individuals.get_mut(index)
    }

    /// Compute fitness values for all individuals.
    pub fn compute_fitness(&self, config: &FitnessConfig) -> Vec<FitnessValue<fitness::Normalized>> {
        self.individuals
            .par_iter()
            .map(|ind| {
                let mut fitness = FitnessValue::default();

                // Apply GC content fitness if configured
                if let Some(gc_fitness) = &config.gc_content {
                    fitness *= gc_fitness.individual_fitness(ind);
                }

                // Apply length fitness if configured
                if let Some(len_fitness) = &config.length {
                    fitness *= len_fitness.individual_fitness(ind);
                }

                // Apply sequence similarity fitness if configured
                if let Some(sim_fitness) = &config.seq_similarity {
                    fitness *= sim_fitness.individual_fitness(ind);
                }

                fitness
            })
            .collect()
    }

    /// Select parent pairs using fitness-proportional sampling.
    ///
    /// Returns `n_pairs` of parent pairs. Each pair consists of two distinct
    /// individuals (no self-pairing). Individuals can appear in multiple pairs.
    pub fn select_parents<R: Rng + ?Sized>(
        &self,
        rng: &mut R,
        fitness_values: &[FitnessValue<fitness::Normalized>],
        n_pairs: usize,
    ) -> Vec<(usize, usize)> {
        // Create cumulative distribution for weighted sampling
        let total_fitness: FitnessValue<fitness::Unnormalized> = fitness_values.iter()
            .map(|&f| f.unnormalize()).sum();
        let total_fitness_val = *total_fitness;
        if (total_fitness == FitnessValue::LETHAL_FITNESS) || ((total_fitness_val / fitness_values.len() as f64) == 1.0) {
            // Fitness values are all zeros/all ones - use uniform selection
            return (0..n_pairs)
                .map(|_| {
                    let parent1 = rng.random_range(0..self.size());
                    let mut parent2 = rng.random_range(0..self.size());
                    while parent2 == parent1 && self.size() > 1 {
                        parent2 = rng.random_range(0..self.size());
                    }
                    (parent1, parent2)
                })
                .collect();
        }

        // Weighted selection
        let cumulative: Vec<f64> = fitness_values
            .iter()
            .scan(0.0, |acc, &f| {
                *acc += *f;
                Some(*acc)
            })
            .collect();

        (0..n_pairs)
            .map(|_| {
                // Select first parent
                let r1 = rng.random_range(0.0..total_fitness_val);
                let parent1 = cumulative
                    .iter()
                    .position(|&c| c >= r1)
                    .unwrap_or(self.size() - 1);

                // Select second parent (different from first)
                let mut parent2 = parent1;
                while parent2 == parent1 && self.size() > 1 {
                    let r2 = rng.random_range(0.0..total_fitness_val);
                    parent2 = cumulative
                        .iter()
                        .position(|&c| c >= r2)
                        .unwrap_or(self.size() - 1);
                }

                (parent1, parent2)
            })
            .collect()
    }

    /// Update fitness values for all individuals in the population.
    pub fn update_fitness(&mut self, config: &FitnessConfig) {
        let fitness_values = self.compute_fitness(config);
        for (ind, fitness) in self.individuals.iter_mut().zip(fitness_values.iter()) {
            ind.set_cached_fitness(*fitness);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::base::{Nucleotide, Sequence};
    use crate::genome::{Chromosome, Haplotype, RepeatMap};
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    fn create_test_individual(id: &str, base: Nucleotide) -> Individual {
        let mut seq = Sequence::with_capacity(100);
        for _ in 0..100 {
            seq.push(base);
        }

        // 100 bp. RU=10. RUs=10. HOR=5 RUs. HORs=2.
        let map = RepeatMap::uniform(10, 5, 2);

        let chr = Chromosome::new(format!("chr_{id}"), seq, map);
        let mut hap1 = Haplotype::new();
        hap1.push(chr.clone());
        let mut hap2 = Haplotype::new();
        hap2.push(chr);

        Individual::new(id, hap1, hap2)
    }

    #[test]
    fn test_population_new() {
        let individuals = vec![
            create_test_individual("ind1", Nucleotide::A),
            create_test_individual("ind2", Nucleotide::C),
        ];

        let pop = Population::new("pop1", individuals);
        assert_eq!(pop.size(), 2);
        assert_eq!(pop.generation(), 0);
        assert_eq!(pop.id(), "pop1");
    }

    #[test]
    fn test_population_increment_generation() {
        let individuals = vec![create_test_individual("ind1", Nucleotide::A)];
        let mut pop = Population::new("pop1", individuals);

        assert_eq!(pop.generation(), 0);
        pop.increment_generation();
        assert_eq!(pop.generation(), 1);
        pop.increment_generation();
        assert_eq!(pop.generation(), 2);
    }

    #[test]
    fn test_population_size() {
        let individuals = vec![
            create_test_individual("ind1", Nucleotide::A),
            create_test_individual("ind2", Nucleotide::C),
            create_test_individual("ind3", Nucleotide::G),
        ];

        let pop = Population::new("pop1", individuals);
        assert_eq!(pop.size(), 3);
        assert!(!pop.is_empty());
    }

    #[test]
    fn test_population_empty() {
        let pop = Population::new("pop1", Vec::new());
        assert_eq!(pop.size(), 0);
        assert!(pop.is_empty());
    }

    #[test]
    fn test_population_get() {
        let individuals = vec![
            create_test_individual("ind1", Nucleotide::A),
            create_test_individual("ind2", Nucleotide::C),
        ];

        let pop = Population::new("pop1", individuals);

        assert!(pop.get(0).is_some());
        assert_eq!(pop.get(0).unwrap().id(), "ind1");
        assert!(pop.get(1).is_some());
        assert_eq!(pop.get(1).unwrap().id(), "ind2");
        assert!(pop.get(2).is_none());
    }

    #[test]
    fn test_population_compute_fitness_neutral() {
        let individuals = vec![
            create_test_individual("ind1", Nucleotide::A),
            create_test_individual("ind2", Nucleotide::C),
        ];

        let pop = Population::new("pop1", individuals);
        let config = FitnessConfig::neutral();
        let fitness_values = pop.compute_fitness(&config);

        // Neutral fitness should be 1.0 for all
        assert_eq!(fitness_values.len(), 2);
        assert_eq!(fitness_values[0], FitnessValue::default());
        assert_eq!(fitness_values[1], FitnessValue::default());
    }

    #[test]
    fn test_population_update_fitness() {
        let individuals = vec![
            create_test_individual("ind1", Nucleotide::A),
            create_test_individual("ind2", Nucleotide::C),
        ];

        let mut pop = Population::new("pop1", individuals);
        let config = FitnessConfig::neutral();

        // Initially cached fitness should be None (not yet computed)
        assert_eq!(pop.get(0).unwrap().cached_fitness(), None);

        pop.update_fitness(&config);

        // After update, should be Some(1.0) (neutral)
        assert_eq!(pop.get(0).unwrap().cached_fitness(), Some(FitnessValue::new_normalized(1.0)));
        assert_eq!(pop.get(1).unwrap().cached_fitness(), Some(FitnessValue::new_normalized(1.0)));
    }

    #[test]
    fn test_population_select_parents_uniform() {
        let individuals = vec![
            create_test_individual("ind1", Nucleotide::A),
            create_test_individual("ind2", Nucleotide::C),
            create_test_individual("ind3", Nucleotide::G),
            create_test_individual("ind4", Nucleotide::T),
        ];

        let pop = Population::new("pop1", individuals);
        let fitness_values = vec![
            FitnessValue::default(),
            FitnessValue::default(),
            FitnessValue::default(),
            FitnessValue::default(),
        ]; // All 1.0 - uniform selection

        let mut rng = StdRng::seed_from_u64(42);
        let pairs = pop.select_parents(&mut rng, &fitness_values, 5);

        assert_eq!(pairs.len(), 5);
        for (p1, p2) in pairs {
            assert!(p1 < pop.size());
            assert!(p2 < pop.size());
            assert_ne!(p1, p2); // No self-pairing
        }
    }

    #[test]
    fn test_population_select_parents_weighted() {
        let individuals = vec![
            create_test_individual("ind1", Nucleotide::A),
            create_test_individual("ind2", Nucleotide::C),
            create_test_individual("ind3", Nucleotide::G),
        ];

        let pop = Population::new("pop1", individuals);
        let total_fitness = FitnessValue::new(6.0);
        let fitness_values = vec![
            FitnessValue::new(1.0).normalize(total_fitness),
            FitnessValue::new(2.0).normalize(total_fitness),
            FitnessValue::new(3.0).normalize(total_fitness),
        ]; // Weighted selection

        let mut rng = StdRng::seed_from_u64(42);
        let pairs = pop.select_parents(&mut rng, &fitness_values, 10);

        assert_eq!(pairs.len(), 10);
        for (p1, p2) in pairs {
            assert!(p1 < pop.size());
            assert!(p2 < pop.size());
            assert_ne!(p1, p2);
        }
    }

    #[test]
    fn test_population_set_individuals() {
        let individuals1 = vec![create_test_individual("ind1", Nucleotide::A)];
        let mut pop = Population::new("pop1", individuals1);
        assert_eq!(pop.size(), 1);

        let individuals2 = vec![
            create_test_individual("ind2", Nucleotide::C),
            create_test_individual("ind3", Nucleotide::G),
        ];
        pop.set_individuals(individuals2);
        assert_eq!(pop.size(), 2);
    }
}
