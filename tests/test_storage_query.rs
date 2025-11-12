//! Integration tests for storage and query operations.
//! Tests the interaction between recording simulation data and querying it back.

use centrevo::{
    base::{Alphabet, Nucleotide},
    simulation::{
        FitnessConfig, MutationConfig, RecombinationConfig, RepeatStructure, Simulation,
        SimulationConfig,
    },
    storage::{Recorder, RecordingStrategy, QueryBuilder},
};

#[test]
fn test_record_and_query_simulation() {
    let path = "/tmp/test_record_query.sqlite";
    let _ = std::fs::remove_file(path);

    // Create and run simulation
    let alphabet = Alphabet::dna();
    let structure = RepeatStructure::new(alphabet.clone(), Nucleotide::A, 20, 10, 10, 1);
    let mutation = MutationConfig::uniform(alphabet, 0.01).unwrap();
    let recombination = RecombinationConfig::standard(0.01, 0.7, 0.1).unwrap();
    let fitness = FitnessConfig::neutral();
    let config = SimulationConfig::new(20, 10, Some(42));

    let mut sim = Simulation::new(
        structure,
        mutation,
        recombination,
        fitness.clone(),
        config.clone(),
    )
    .unwrap();

    // Create recorder and record metadata
    let mut recorder = Recorder::new(path, "test_sim", RecordingStrategy::All).unwrap();
    recorder.record_metadata(&config).unwrap();

    // Run and record
    for generation in 0..10 {
        sim.step().unwrap();
        recorder.record_generation(sim.population(), generation).unwrap();
    }

    recorder.close().unwrap();

    // Query back the data using QueryBuilder
    let query = QueryBuilder::new(path).unwrap();

    // Test: List simulations
    let sim_list = query.list_simulations().unwrap();
    assert_eq!(sim_list.len(), 1);
    assert_eq!(sim_list[0], "test_sim");

    // Test: Get simulation info
    let info = query.get_simulation_info("test_sim").unwrap();
    assert_eq!(info.sim_id, "test_sim");
    assert_eq!(info.num_generations, 10);

    // Test: Get recorded generations
    let generations = query.get_recorded_generations("test_sim").unwrap();
    assert_eq!(generations.len(), 10);
    for (i, &generation) in generations.iter().enumerate() {
        assert_eq!(generation, i);
    }

    // Test: Get specific generation
    let gen_data = query.get_generation("test_sim", 5).unwrap();
    assert_eq!(gen_data.len(), 20); // 20 individuals

    // Test: Get fitness history
    let fitness_history = query.get_fitness_history("test_sim").unwrap();
    assert_eq!(fitness_history.len(), 10); // 10 generations

    std::fs::remove_file(path).ok();
}

#[test]
fn test_partial_recording_strategy() {
    let path = "/tmp/test_partial_recording.sqlite";
    let _ = std::fs::remove_file(path);

    // Create simulation
    let alphabet = Alphabet::dna();
    let structure = RepeatStructure::new(alphabet.clone(), Nucleotide::C, 15, 10, 10, 1);
    let mutation = MutationConfig::uniform(alphabet, 0.01).unwrap();
    let recombination = RecombinationConfig::standard(0.01, 0.7, 0.1).unwrap();
    let fitness = FitnessConfig::neutral();
    let config = SimulationConfig::new(30, 20, Some(123));

    let mut sim = Simulation::new(
        structure,
        mutation,
        recombination,
        fitness.clone(),
        config.clone(),
    )
    .unwrap();

    // Use selective recording - record every 5 generations
    let mut recorder = Recorder::new(path, "partial_sim", RecordingStrategy::EveryN(5)).unwrap();
    recorder.record_metadata(&config).unwrap();

    // Run simulation
    for generation in 0..20 {
        sim.step().unwrap();
        if generation % 5 == 0 {
            recorder
                .record_generation(sim.population(), generation)
                .unwrap();
        }
    }

    recorder.close().unwrap();

    // Query and verify only selected generations were recorded
    let query = QueryBuilder::new(path).unwrap();
    let generations = query.get_recorded_generations("partial_sim").unwrap();

    assert_eq!(generations.len(), 4); // Generations 0, 5, 10, 15
    assert_eq!(generations, vec![0, 5, 10, 15]);

    std::fs::remove_file(path).ok();
}

#[test]
fn test_multiple_simulations_in_one_database() {
    let path = "/tmp/test_multi_sim.sqlite";
    let _ = std::fs::remove_file(path);

    let alphabet = Alphabet::dna();
    let structure = RepeatStructure::new(alphabet.clone(), Nucleotide::G, 10, 10, 10, 1);
    let mutation = MutationConfig::uniform(alphabet, 0.01).unwrap();
    let recombination = RecombinationConfig::standard(0.01, 0.7, 0.1).unwrap();
    let fitness = FitnessConfig::neutral();

    // Run two simulations with different parameters
    for (sim_id, pop_size, seed) in [("sim_a", 20, 42), ("sim_b", 15, 123)] {
        let config = SimulationConfig::new(pop_size, 5, Some(seed));
        let mut sim = Simulation::new(
            structure.clone(),
            mutation.clone(),
            recombination.clone(),
            fitness.clone(),
            config.clone(),
        )
        .unwrap();

        let mut recorder = Recorder::new(path, sim_id, RecordingStrategy::All).unwrap();
        recorder.record_metadata(&config).unwrap();

        for generation in 0..5 {
            sim.step().unwrap();
            recorder
                .record_generation(sim.population(), generation)
                .unwrap();
        }

        recorder.close().unwrap();
    }

    // Query database
    let query = QueryBuilder::new(path).unwrap();

    // Should have both simulations
    let sim_list = query.list_simulations().unwrap();
    assert_eq!(sim_list.len(), 2);
    assert!(sim_list.contains(&"sim_a".to_string()));
    assert!(sim_list.contains(&"sim_b".to_string()));

    // Verify each simulation's data
    let info_a = query.get_simulation_info("sim_a").unwrap();
    let info_b = query.get_simulation_info("sim_b").unwrap();

    assert_eq!(info_a.num_generations, 5);
    assert_eq!(info_b.num_generations, 5);

    std::fs::remove_file(path).ok();
}

#[test]
fn test_fitness_tracking_over_generations() {
    let path = "/tmp/test_fitness_tracking.sqlite";
    let _ = std::fs::remove_file(path);

    // Create simulation with selection
    let alphabet = Alphabet::dna();
    let structure = RepeatStructure::new(alphabet.clone(), Nucleotide::A, 20, 10, 5, 1);
    let mutation = MutationConfig::uniform(alphabet, 0.01).unwrap();
    let recombination = RecombinationConfig::standard(0.01, 0.7, 0.1).unwrap();

    use centrevo::evolution::GCContentFitness;
    let gc_fitness = GCContentFitness::new(0.5, 3.0).unwrap();
    let fitness = FitnessConfig::new(Some(gc_fitness), None, None);

    let config = SimulationConfig::new(25, 20, Some(456));

    let mut sim = Simulation::new(
        structure,
        mutation,
        recombination,
        fitness.clone(),
        config.clone(),
    )
    .unwrap();

    let mut recorder = Recorder::new(path, "fitness_sim", RecordingStrategy::All).unwrap();
    recorder.record_metadata(&config).unwrap();

    // Run and record
    for generation in 0..20 {
        sim.step().unwrap();
        recorder
            .record_generation(sim.population(), generation)
            .unwrap();
    }

    recorder.close().unwrap();

    // Query fitness history
    let query = QueryBuilder::new(path).unwrap();
    let fitness_history = query.get_fitness_history("fitness_sim").unwrap();

    assert_eq!(fitness_history.len(), 20);

    // Verify fitness data structure (it's a tuple of (generation, FitnessStats))
    for (generation, stats) in &fitness_history {
        assert!(generation < &20);
        assert!(stats.mean >= 0.0);
        assert!(stats.std >= 0.0);
    }

    println!("\nFitness evolution:");
    println!("  Gen 0:  mean={:.4}", fitness_history[0].1.mean);
    println!("  Gen 19: mean={:.4}", fitness_history[19].1.mean);

    std::fs::remove_file(path).ok();
}

#[test]
fn test_query_nonexistent_data() {
    let path = "/tmp/test_query_errors.sqlite";
    let _ = std::fs::remove_file(path);

    // Create simulation and database
    let alphabet = Alphabet::dna();
    let structure = RepeatStructure::new(alphabet.clone(), Nucleotide::T, 10, 10, 10, 1);
    let mutation = MutationConfig::uniform(alphabet, 0.01).unwrap();
    let recombination = RecombinationConfig::standard(0.01, 0.7, 0.1).unwrap();
    let fitness = FitnessConfig::neutral();
    let config = SimulationConfig::new(10, 1, Some(789));

    let mut sim = Simulation::new(structure, mutation, recombination, fitness, config.clone()).unwrap();
    let mut recorder = Recorder::new(path, "dummy_sim", RecordingStrategy::All).unwrap();
    recorder.record_metadata(&config).unwrap();
    sim.step().unwrap();
    recorder.record_generation(sim.population(), 0).unwrap();
    recorder.close().unwrap();

    // Query using QueryBuilder
    let query = QueryBuilder::new(path).unwrap();

    // Querying nonexistent simulation should return error
    let result = query.get_simulation_info("nonexistent");
    assert!(result.is_err());

    // Should have the sim we created
    let sim_list = query.list_simulations().unwrap();
    assert_eq!(sim_list.len(), 1);

    std::fs::remove_file(path).ok();
}
