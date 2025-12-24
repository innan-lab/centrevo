//! Population management and operations.
//!
//! This module provides structures and functions for managing populations
//! of individuals during evolutionary simulations.

use crate::base::FitnessValue;
use crate::genome::Individual;
use rand::Rng;
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

    /// Select parent pairs using fitness-proportional sampling.
    ///
    /// Returns `n_pairs` of parent pairs. Each pair consists of two distinct
    /// individuals (no self-pairing). Individuals can appear in multiple pairs.
    ///
    /// # Panics
    ///
    /// Panics if any individual does not have a cached fitness value (computed at birth).
    pub fn select_parents<R: Rng + ?Sized>(
        &self,
        rng: &mut R,
        n_pairs: usize,
    ) -> Vec<(usize, usize)> {
        // Collect fitness values from individuals
        // This is a transient vector used only for selection
        let fitness_values: Vec<FitnessValue> = self
            .individuals
            .iter()
            .map(|ind| {
                ind.cached_fitness()
                    .expect("Individual fitness not computed (should be intrinsic)")
            })
            .collect();

        // Create cumulative distribution for weighted sampling
        let total_fitness: FitnessValue = fitness_values.iter().copied().sum();
        let total_fitness_val = *total_fitness;

        // If either all fitnesses are lethal or all values are equal then selection
        // is uniform; the 'all equal' case corresponds to no selection pressure.
        let all_equal = if fitness_values.is_empty() {
            true
        } else {
            let first = fitness_values[0];
            fitness_values.iter().all(|&f| f == first)
        };

        if (total_fitness == FitnessValue::LETHAL_FITNESS) || all_equal {
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
}

#[cfg(test)]
mod tests {
    use super::{FitnessValue, Individual, Population};
    use crate::base::{GenomeArena, Nucleotide, Sequence};
    use crate::genome::{Chromosome, Haplotype, RepeatMap};
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    fn create_test_individual(id: &str, base: Nucleotide, arena: &mut GenomeArena) -> Individual {
        let mut seq = Sequence::with_capacity(100);
        for _ in 0..100 {
            seq.push(base);
        }

        // 100 bp. RU=10. RUs=10. HOR=5 RUs. HORs=2.
        let map = RepeatMap::uniform(10, 5, 2);

        let slice = arena.alloc(seq.as_slice());
        let chr = Chromosome::new(format!("chr_{id}"), slice, map);
        let mut hap1 = Haplotype::new();
        hap1.push(chr.clone());
        let mut hap2 = Haplotype::new();
        hap2.push(chr);

        Individual::new(id, hap1, hap2)
    }

    #[test]
    fn test_population_select_parents_selection_probabilities() {
        // Create individuals with specific fitnesses
        let mut arena = GenomeArena::new();
        let mut ind1 = create_test_individual("ind1", Nucleotide::A, &mut arena);
        ind1.set_cached_fitness(FitnessValue::new(10.0));

        let mut ind2 = create_test_individual("ind2", Nucleotide::C, &mut arena);
        ind2.set_cached_fitness(FitnessValue::new(1.0));

        let pop = Population::new("pop1", vec![ind1, ind2]);
        assert_eq!(pop.size(), 2);
        assert_eq!(pop.generation(), 0);
        assert_eq!(pop.id(), "pop1");
    }

    #[test]
    fn test_population_new() {
        let mut arena = GenomeArena::new();
        let individuals = vec![
            create_test_individual("ind1", Nucleotide::A, &mut arena),
            create_test_individual("ind2", Nucleotide::C, &mut arena),
        ];

        let pop = Population::new("pop1", individuals);
        assert_eq!(pop.size(), 2);
        assert_eq!(pop.generation(), 0);
        assert_eq!(pop.id(), "pop1");
    }

    #[test]
    fn test_population_increment_generation() {
        let mut arena = GenomeArena::new();
        let individuals = vec![create_test_individual("ind1", Nucleotide::A, &mut arena)];
        let mut pop = Population::new("pop1", individuals);

        assert_eq!(pop.generation(), 0);
        pop.increment_generation();
        assert_eq!(pop.generation(), 1);
        pop.increment_generation();
        assert_eq!(pop.generation(), 2);
    }

    #[test]
    fn test_population_size() {
        let mut arena = GenomeArena::new();
        let individuals = vec![
            create_test_individual("ind1", Nucleotide::A, &mut arena),
            create_test_individual("ind2", Nucleotide::C, &mut arena),
            create_test_individual("ind3", Nucleotide::G, &mut arena),
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
        let mut arena = GenomeArena::new();
        let individuals = vec![
            create_test_individual("ind1", Nucleotide::A, &mut arena),
            create_test_individual("ind2", Nucleotide::C, &mut arena),
        ];

        let pop = Population::new("pop1", individuals);

        assert!(pop.get(0).is_some());
        assert_eq!(pop.get(0).unwrap().id(), "ind1");
        assert!(pop.get(1).is_some());
        assert_eq!(pop.get(1).unwrap().id(), "ind2");
        assert!(pop.get(2).is_none());
    }

    #[test]
    fn test_population_select_parents_uniform() {
        let mut arena = GenomeArena::new();
        let individuals = vec![
            create_test_individual("ind1", Nucleotide::A, &mut arena),
            create_test_individual("ind2", Nucleotide::C, &mut arena),
            create_test_individual("ind3", Nucleotide::G, &mut arena),
            create_test_individual("ind4", Nucleotide::T, &mut arena),
        ];

        // Manually set fitness for testing
        let mut individuals = individuals;
        for ind in &mut individuals {
            ind.set_cached_fitness(FitnessValue::new(1.0));
        }

        let pop = Population::new("pop1", individuals);

        let mut rng = StdRng::seed_from_u64(42);
        let pairs = pop.select_parents(&mut rng, 5);

        assert_eq!(pairs.len(), 5);
        for (p1, p2) in pairs {
            assert!(p1 < pop.size());
            assert!(p2 < pop.size());
            assert_ne!(p1, p2); // No self-pairing
        }
    }

    #[test]
    fn test_population_select_parents_weighted() {
        let mut arena = GenomeArena::new();
        let individuals = vec![
            create_test_individual("ind1", Nucleotide::A, &mut arena),
            create_test_individual("ind2", Nucleotide::C, &mut arena),
            create_test_individual("ind3", Nucleotide::G, &mut arena),
        ];

        let mut individuals = individuals;
        // Manually set fitness: 1.0, 2.0, 3.0
        individuals[0].set_cached_fitness(FitnessValue::new(1.0));
        individuals[1].set_cached_fitness(FitnessValue::new(2.0));
        individuals[2].set_cached_fitness(FitnessValue::new(3.0));

        let pop = Population::new("pop1", individuals);

        let mut rng = StdRng::seed_from_u64(42);
        let pairs = pop.select_parents(&mut rng, 10);

        assert_eq!(pairs.len(), 10);
        for (p1, p2) in pairs {
            assert!(p1 < pop.size());
            assert!(p2 < pop.size());
            assert_ne!(p1, p2);
        }
    }

    // Removed tests for infinite fitness values: +infinity cases are now
    // disallowed at construction time via FitnessValue::new.

    #[test]
    fn test_population_set_individuals() {
        let mut arena = GenomeArena::new();
        let individuals1 = vec![create_test_individual("ind1", Nucleotide::A, &mut arena)];
        let mut pop = Population::new("pop1", individuals1);
        assert_eq!(pop.size(), 1);

        let individuals2 = vec![
            create_test_individual("ind2", Nucleotide::C, &mut arena),
            create_test_individual("ind3", Nucleotide::G, &mut arena),
        ];
        pop.set_individuals(individuals2);
        assert_eq!(pop.size(), 2);
    }
}
