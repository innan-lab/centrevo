//! Test reproducibility of parallel simulation with fixed seeds.

use centrevo::base::Nucleotide;
use centrevo::genome::Individual;
use centrevo::simulation::{
    FitnessConfig, MutationConfig, RecombinationConfig, RepeatStructure, Simulation,
    SimulationConfig,
};

#[test]
fn test_parallel_reproducibility() {
    // Run the same simulation twice with the same seed
    let results1 = run_simulation(42);
    let results2 = run_simulation(42);

    // Results should be identical
    assert_eq!(results1.len(), results2.len());
    for (ind1, ind2) in results1.iter().zip(results2.iter()) {
        // Compare haplotype1 sequences
        for (chr1, chr2) in ind1
            .haplotype1()
            .chromosomes()
            .iter()
            .zip(ind2.haplotype1().chromosomes().iter())
        {
            assert_eq!(
                chr1.sequence().to_string(),
                chr2.sequence().to_string(),
                "Haplotype1 sequences differ"
            );
        }

        // Compare haplotype2 sequences
        for (chr1, chr2) in ind1
            .haplotype2()
            .chromosomes()
            .iter()
            .zip(ind2.haplotype2().chromosomes().iter())
        {
            assert_eq!(
                chr1.sequence().to_string(),
                chr2.sequence().to_string(),
                "Haplotype2 sequences differ"
            );
        }
    }
}

#[test]
fn test_parallel_different_seeds() {
    // Run simulations with different seeds
    let results1 = run_simulation(42);
    let results2 = run_simulation(123);

    // Results should be different
    assert_eq!(results1.len(), results2.len());

    let mut different_found = false;
    for (ind1, ind2) in results1.iter().zip(results2.iter()) {
        for (chr1, chr2) in ind1
            .haplotype1()
            .chromosomes()
            .iter()
            .zip(ind2.haplotype1().chromosomes().iter())
        {
            if chr1.sequence().to_string() != chr2.sequence().to_string() {
                different_found = true;
                break;
            }
        }
        if different_found {
            break;
        }
    }

    assert!(
        different_found,
        "Simulations with different seeds should produce different results"
    );
}

fn run_simulation(seed: u64) -> Vec<Individual> {
    let structure = RepeatStructure::new(
        Nucleotide::A,
        10, // ru_length
        5,  // rus_per_hor
        10, // hors_per_chr
        1,  // chrs_per_hap
    );

    let mutation = MutationConfig::uniform(0.01).unwrap();
    let recombination = RecombinationConfig::standard(0.01, 0.7, 0.1).unwrap();
    let fitness = FitnessConfig::neutral();
    let config = SimulationConfig::new(20, 5, Some(seed));

    let mut sim = Simulation::new(structure, mutation, recombination, fitness, config).unwrap();
    sim.run().unwrap();

    sim.population().individuals().to_vec()
}
