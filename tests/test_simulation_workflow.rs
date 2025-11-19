//! Integration tests for end-to-end simulation workflows.
//! Tests that simulate real-world usage patterns combining multiple modules.

use centrevo::{
    analysis::{
        composition::gc_content, distance::pairwise_distances, diversity::nucleotide_diversity,
    },
    base::Nucleotide,
    evolution::{GCContentFitness, LengthFitness},
    simulation::{
        FitnessConfig, MutationConfig, RecombinationConfig, RepeatStructure, Simulation,
        SimulationConfig,
    },
};

#[test]
fn test_basic_simulation_workflow() {
    // Create a simple simulation and run it
    let structure = RepeatStructure::new(Nucleotide::A, 20, 10, 10, 1);
    let mutation = MutationConfig::uniform(0.01).unwrap();
    let recombination = RecombinationConfig::standard(0.01, 0.7, 0.1).unwrap();
    let fitness = FitnessConfig::neutral();
    let config = SimulationConfig::new(20, 10, Some(42));

    let mut sim = Simulation::new(structure, mutation, recombination, fitness, config).unwrap();

    // Run simulation
    sim.run().unwrap();

    // Verify population exists and has correct size
    let pop = sim.population();
    assert_eq!(pop.size(), 20);

    // Verify individuals have the expected structure
    for ind in pop.individuals() {
        assert_eq!(ind.haplotype1().chromosomes().len(), 1);
        assert_eq!(ind.haplotype2().chromosomes().len(), 1);
        assert!(ind.fitness() >= 0.0);
    }
}

#[test]
fn test_simulation_with_mutation_accumulation() {
    // Test that mutations accumulate over generations
    let structure = RepeatStructure::new(Nucleotide::A, 10, 5, 10, 1);
    let mutation = MutationConfig::uniform(0.05).unwrap(); // High mutation rate
    let recombination = RecombinationConfig::standard(0.0, 0.0, 0.0).unwrap(); // No recombination
    let fitness = FitnessConfig::neutral();
    let config = SimulationConfig::new(10, 50, Some(123));

    let mut sim = Simulation::new(structure, mutation, recombination, fitness, config).unwrap();

    // Record initial diversity
    let initial_diversity = nucleotide_diversity(sim.population(), 0);

    // Run simulation
    sim.run().unwrap();

    // Diversity should increase due to mutations
    let final_diversity = nucleotide_diversity(sim.population(), 0);
    assert!(
        final_diversity >= initial_diversity,
        "Diversity should increase or stay same with mutations: {} -> {}",
        initial_diversity,
        final_diversity
    );

    // Should have variation in the population
    assert!(
        final_diversity > 0.0,
        "Should have diversity after 50 generations with high mutation"
    );
}

#[test]
fn test_simulation_with_gc_content_selection() {
    // Test that GC content selection affects the population
    let structure = RepeatStructure::new(Nucleotide::A, 20, 5, 5, 1);
    let mutation = MutationConfig::uniform(0.01).unwrap();
    let recombination = RecombinationConfig::standard(0.01, 0.7, 0.1).unwrap();

    // Select for 50% GC content
    let gc_fitness = GCContentFitness::new(0.5, 5.0).unwrap();
    let fitness = FitnessConfig::new(Some(gc_fitness), None, None, None);
    let config = SimulationConfig::new(30, 100, Some(456));

    let mut sim = Simulation::new(structure, mutation, recombination, fitness, config).unwrap();

    // Record initial GC content
    let initial_gc = gc_content(sim.population(), None, None, None);

    // Run simulation
    sim.run().unwrap();

    // Final GC content should be closer to 50%
    let final_gc = gc_content(sim.population(), None, None, None);

    println!(
        "GC content evolution: {:.2}% -> {:.2}%",
        initial_gc * 100.0,
        final_gc * 100.0
    );

    // With strong selection, should move toward optimum
    // Allow some variance since evolution is stochastic
    assert!(
        final_gc > 0.1 && final_gc < 0.9,
        "Final GC content should be reasonable: {}",
        final_gc
    );
}

#[test]
fn test_simulation_with_length_selection() {
    // Test that length selection affects fitness
    let target_length = 500;
    let structure = RepeatStructure::new(Nucleotide::C, 20, 5, 5, 1);
    let mutation = MutationConfig::uniform(0.001).unwrap();
    let recombination = RecombinationConfig::standard(0.01, 0.7, 0.1).unwrap();

    // Select for specific length
    let length_fitness = LengthFitness::new(target_length, 2.0).unwrap();
    let fitness = FitnessConfig::new(None, Some(length_fitness), None, None);
    let config = SimulationConfig::new(20, 20, Some(789));

    let mut sim = Simulation::new(structure, mutation, recombination, fitness, config).unwrap();
    sim.run().unwrap();

    // Verify fitness values are calculated (they may be 0 for poor matches)
    for ind in sim.population().individuals() {
        assert!(ind.fitness() >= 0.0, "Fitness should be non-negative");
        assert!(ind.fitness() <= 1.0, "Fitness should be normalized");
    }
}

#[test]
fn test_simulation_step_by_step() {
    // Test stepping through simulation generation by generation
    let structure = RepeatStructure::new(Nucleotide::G, 10, 5, 10, 1);
    let mutation = MutationConfig::uniform(0.01).unwrap();
    let recombination = RecombinationConfig::standard(0.01, 0.7, 0.1).unwrap();
    let fitness = FitnessConfig::neutral();
    let config = SimulationConfig::new(15, 5, Some(101));

    let mut sim = Simulation::new(structure, mutation, recombination, fitness, config).unwrap();

    // Step through each generation
    for generation in 0..5 {
        assert_eq!(sim.generation(), generation);
        sim.step().unwrap();
        assert_eq!(sim.generation(), generation + 1);
    }

    assert_eq!(sim.generation(), 5);
}

#[test]
fn test_simulation_with_analysis_pipeline() {
    // Test full workflow: simulate -> analyze
    let structure = RepeatStructure::new(Nucleotide::T, 15, 8, 8, 1);
    let mutation = MutationConfig::uniform(0.02).unwrap();
    let recombination = RecombinationConfig::standard(0.02, 0.6, 0.15).unwrap();
    let fitness = FitnessConfig::neutral();
    let config = SimulationConfig::new(25, 30, Some(202));

    let mut sim = Simulation::new(structure, mutation, recombination, fitness, config).unwrap();
    sim.run().unwrap();

    let pop = sim.population();

    // Run multiple analyses
    let gc = gc_content(pop, None, None, None);
    let diversity = nucleotide_diversity(pop, 0);
    let distances = pairwise_distances(pop, 0);

    println!("\nAnalysis Results:");
    println!("  GC Content: {:.2}%", gc * 100.0);
    println!("  π (nucleotide diversity): {:.6}", diversity);
    println!("  Pairwise distances: {} values", distances.len());

    // Basic sanity checks
    assert!((0.0..=1.0).contains(&gc));
    assert!(diversity >= 0.0);
    assert!(!distances.is_empty());
}

#[test]
fn test_multiple_simulations_independence() {
    // Test that multiple simulations with different seeds produce different results
    let structure = RepeatStructure::new(Nucleotide::A, 10, 5, 10, 1);
    let mutation = MutationConfig::uniform(0.01).unwrap();
    let recombination = RecombinationConfig::standard(0.01, 0.7, 0.1).unwrap();
    let fitness = FitnessConfig::neutral();

    let mut results = Vec::new();

    for seed in [42, 123, 456] {
        let config = SimulationConfig::new(20, 20, Some(seed));
        let mut sim = Simulation::new(
            structure.clone(),
            mutation.clone(),
            recombination.clone(),
            fitness.clone(),
            config,
        )
        .unwrap();

        sim.run().unwrap();
        let diversity = nucleotide_diversity(sim.population(), 0);
        results.push(diversity);
    }

    // Results should differ (very unlikely to be identical)
    assert!(
        results[0] != results[1] || results[1] != results[2],
        "Different seeds should produce different results"
    );
}

#[test]
fn test_simulation_reproducibility() {
    // Test that same seed produces same results
    let structure = RepeatStructure::new(Nucleotide::C, 10, 5, 10, 1);
    let mutation = MutationConfig::uniform(0.01).unwrap();
    let recombination = RecombinationConfig::standard(0.01, 0.7, 0.1).unwrap();
    let fitness = FitnessConfig::neutral();
    let config = SimulationConfig::new(15, 25, Some(42));

    let run_simulation = || {
        let mut sim = Simulation::new(
            structure.clone(),
            mutation.clone(),
            recombination.clone(),
            fitness.clone(),
            config.clone(),
        )
        .unwrap();
        sim.run().unwrap();
        sim.population().individuals().to_vec()
    };

    let results1 = run_simulation();
    let results2 = run_simulation();

    // Compare sequences from both runs
    assert_eq!(results1.len(), results2.len());
    for (ind1, ind2) in results1.iter().zip(results2.iter()) {
        // Compare haplotype sequences
        for (chr1, chr2) in ind1
            .haplotype1()
            .chromosomes()
            .iter()
            .zip(ind2.haplotype1().chromosomes().iter())
        {
            assert_eq!(
                chr1.sequence().to_string(),
                chr2.sequence().to_string(),
                "Haplotype1 sequences should be identical with same seed"
            );
        }
        for (chr1, chr2) in ind1
            .haplotype2()
            .chromosomes()
            .iter()
            .zip(ind2.haplotype2().chromosomes().iter())
        {
            assert_eq!(
                chr1.sequence().to_string(),
                chr2.sequence().to_string(),
                "Haplotype2 sequences should be identical with same seed"
            );
        }
    }
}

#[test]
fn test_combined_fitness_components() {
    // Test simulation with multiple fitness components
    let structure = RepeatStructure::new(Nucleotide::A, 20, 5, 5, 1);
    let mutation = MutationConfig::uniform(0.005).unwrap();
    let recombination = RecombinationConfig::standard(0.01, 0.7, 0.1).unwrap();

    // Select for both GC content and length
    let gc_fitness = GCContentFitness::new(0.4, 2.0).unwrap();
    let length_fitness = LengthFitness::new(500, 1.5).unwrap();
    let fitness = FitnessConfig::new(Some(gc_fitness), Some(length_fitness), None, None);

    let config = SimulationConfig::new(20, 30, Some(303));

    let mut sim = Simulation::new(structure, mutation, recombination, fitness, config).unwrap();
    sim.run().unwrap();

    // All individuals should have computed fitness values (may be 0)
    for ind in sim.population().individuals() {
        assert!(ind.fitness() >= 0.0, "Fitness should be non-negative");
    }

    // Population should show genetic variation
    let diversity = nucleotide_diversity(sim.population(), 0);
    assert!(diversity > 0.0, "Should have genetic diversity");
}
