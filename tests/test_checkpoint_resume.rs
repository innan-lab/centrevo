//! Integration tests for checkpoint and resume functionality.

use centrevo::base::Nucleotide;
use centrevo::simulation::{
    FitnessConfig, MutationConfig, RecombinationConfig, RepeatStructure, Simulation,
    SimulationConfig,
};
use centrevo::storage::{Recorder, RecordingStrategy, SimulationSnapshot};
use std::path::PathBuf;

/// Helper to create a test database path
fn test_db_path(name: &str) -> PathBuf {
    PathBuf::from(format!("/tmp/centrevo_test_checkpoint_{}.sqlite", name))
}

/// Helper to remove test database
fn cleanup_db(path: &PathBuf) {
    std::fs::remove_file(path).ok();
    std::fs::remove_file(format!("{}-wal", path.display())).ok();
    std::fs::remove_file(format!("{}-shm", path.display())).ok();
}

#[test]
fn test_checkpoint_and_resume_basic() {
    let db_path = test_db_path("basic");
    cleanup_db(&db_path);

    // Setup simulation parameters
    let structure = RepeatStructure::new(
        Nucleotide::A,
        10, // Small for testing
        5,
        10,
        1,
    );

    let mutation = MutationConfig::uniform(0.001).unwrap();
    let recombination = RecombinationConfig::standard(0.01, 0.7, 0.1).unwrap();
    let fitness = FitnessConfig::neutral();
    let config = SimulationConfig::new(10, 100, Some(42)); // 100 generations

    let sim_id = "test_sim";

    // Part 1: Run simulation for first 50 generations with checkpoints
    {
        let mut sim = Simulation::new(
            structure.clone(),
            mutation.clone(),
            recombination.clone(),
            fitness.clone(),
            config.clone(),
        )
        .unwrap();

        // Setup recorder with full config
        let snapshot = SimulationSnapshot {
            structure: structure.clone(),
            mutation: mutation.clone(),
            recombination: recombination.clone(),
            fitness: fitness.clone(),
            config: config.clone(),
        };

        let mut recorder = Recorder::new(
            &db_path,
            sim_id,
            RecordingStrategy::EveryN(10), // Record every 10 generations
        )
        .unwrap();

        recorder.record_full_config(&snapshot).unwrap();

        // Run first 50 generations
        for generation in 1..=50 {
            sim.step().unwrap();

            if recorder.should_record(generation) {
                let rng_state = sim.rng_state_bytes();
                recorder
                    .record_generation(sim.population(), generation)
                    .unwrap();
                recorder.record_checkpoint(generation, &rng_state).unwrap();
            }
        }

        recorder.close().ok();

        assert_eq!(sim.generation(), 50);
        println!("✓ Completed first 50 generations");
    }

    // Part 2: Resume from checkpoint and continue
    {
        let mut sim = Simulation::from_checkpoint(&db_path, sim_id).unwrap();

        // Verify we resumed from generation 50 (last checkpoint at gen 50)
        assert_eq!(sim.generation(), 50, "Should resume from generation 50");

        println!("✓ Resumed from generation {}", sim.generation());

        // Setup recorder again
        let mut recorder = Recorder::new(&db_path, sim_id, RecordingStrategy::EveryN(10)).unwrap();

        // Continue for remaining 50 generations
        for generation in 51..=100 {
            sim.step().unwrap();

            if recorder.should_record(generation) {
                let rng_state = sim.rng_state_bytes();
                recorder
                    .record_generation(sim.population(), generation)
                    .unwrap();
                recorder.record_checkpoint(generation, &rng_state).unwrap();
            }
        }

        recorder.finalize_metadata().unwrap();
        recorder.close().ok();

        assert_eq!(sim.generation(), 100);
        println!("✓ Completed remaining 50 generations");
    }

    // Part 3: Verify final state
    {
        use centrevo::storage::QueryBuilder;

        let query = QueryBuilder::new(&db_path).unwrap();

        // Check recorded generations
        let recorded_gens = query.get_recorded_generations(sim_id).unwrap();
        assert!(!recorded_gens.is_empty());
        assert_eq!(*recorded_gens.last().unwrap(), 100);

        println!("✓ Recorded generations: {:?}", recorded_gens);

        // Verify we can load final population
        let final_pop = query.get_generation(sim_id, 100).unwrap();
        assert_eq!(final_pop.len(), 10);

        println!("✓ Final population has {} individuals", final_pop.len());

        query.close().ok();
    }

    cleanup_db(&db_path);
    println!("✅ Checkpoint and resume test passed!");
}

#[test]
fn test_resume_multiple_times() {
    let db_path = test_db_path("multiple");
    cleanup_db(&db_path);

    // Setup
    let structure = RepeatStructure::new(Nucleotide::A, 10, 5, 10, 1);

    let mutation = MutationConfig::uniform(0.001).unwrap();
    let recombination = RecombinationConfig::standard(0.01, 0.7, 0.1).unwrap();
    let fitness = FitnessConfig::neutral();
    let config = SimulationConfig::new(5, 100, Some(123));

    let sim_id = "multi_resume_sim";

    // Run simulation in 4 chunks: 0->25, 25->50, 50->75, 75->100
    let checkpoints = vec![25, 50, 75, 100];

    for (i, &target_gen) in checkpoints.iter().enumerate() {
        let mut sim = if i == 0 {
            // First run: create new simulation
            let mut sim = Simulation::new(
                structure.clone(),
                mutation.clone(),
                recombination.clone(),
                fitness.clone(),
                config.clone(),
            )
            .unwrap();

            // Record full config
            let snapshot = SimulationSnapshot {
                structure: structure.clone(),
                mutation: mutation.clone(),
                recombination: recombination.clone(),
                fitness: fitness.clone(),
                config: config.clone(),
            };

            let mut recorder =
                Recorder::new(&db_path, sim_id, RecordingStrategy::EveryN(25)).unwrap();

            recorder.record_full_config(&snapshot).unwrap();
            recorder.close().ok();

            sim
        } else {
            // Subsequent runs: resume from checkpoint
            Simulation::from_checkpoint(&db_path, sim_id).unwrap()
        };

        let start_gen = sim.generation();
        println!("Run {}: Starting from generation {}", i + 1, start_gen);

        let mut recorder = Recorder::new(&db_path, sim_id, RecordingStrategy::EveryN(25)).unwrap();

        // Run to target generation
        for generation in (start_gen + 1)..=target_gen {
            sim.step().unwrap();

            if recorder.should_record(generation) {
                let rng_state = sim.rng_state_bytes();
                recorder
                    .record_generation(sim.population(), generation)
                    .unwrap();
                recorder.record_checkpoint(generation, &rng_state).unwrap();
            }
        }

        recorder.close().ok();

        assert_eq!(sim.generation(), target_gen);
        println!("✓ Reached generation {}", target_gen);
    }

    cleanup_db(&db_path);
    println!("✅ Multiple resume test passed!");
}

#[test]
fn test_resume_preserves_rng_state() {
    let db_path1 = test_db_path("rng_state_1");
    let db_path2 = test_db_path("rng_state_2");
    cleanup_db(&db_path1);
    cleanup_db(&db_path2);

    // Setup identical simulations
    let structure = RepeatStructure::new(Nucleotide::A, 10, 5, 10, 1);

    let mutation = MutationConfig::uniform(0.001).unwrap();
    let recombination = RecombinationConfig::standard(0.01, 0.7, 0.1).unwrap();
    let fitness = FitnessConfig::neutral();
    let seed = 999u64;
    let config = SimulationConfig::new(10, 100, Some(seed));

    // Simulation 1: Run straight through to gen 100
    let final_pop_continuous = {
        let mut sim = Simulation::new(
            structure.clone(),
            mutation.clone(),
            recombination.clone(),
            fitness.clone(),
            config.clone(),
        )
        .unwrap();

        for _generation in 1..=100 {
            sim.step().unwrap();
        }

        // Get final state
        sim.population().individuals()[0]
            .haplotype1()
            .get(0)
            .unwrap()
            .sequence()
            .as_slice()
            .to_vec()
    };

    // Simulation 2: Run to gen 50, checkpoint, resume, then complete
    let final_pop_resumed = {
        // Part 1: Run to gen 50
        {
            let mut sim = Simulation::new(
                structure.clone(),
                mutation.clone(),
                recombination.clone(),
                fitness.clone(),
                config.clone(),
            )
            .unwrap();

            let snapshot = SimulationSnapshot {
                structure: structure.clone(),
                mutation: mutation.clone(),
                recombination: recombination.clone(),
                fitness: fitness.clone(),
                config: config.clone(),
            };

            let mut recorder =
                Recorder::new(&db_path2, "rng_test_sim", RecordingStrategy::EveryN(50)).unwrap();

            recorder.record_full_config(&snapshot).unwrap();

            for generation in 1..=50 {
                sim.step().unwrap();

                if recorder.should_record(generation) {
                    let rng_state = sim.rng_state_bytes();
                    recorder
                        .record_generation(sim.population(), generation)
                        .unwrap();
                    recorder.record_checkpoint(generation, &rng_state).unwrap();
                }
            }

            recorder.close().ok();
        }

        // Part 2: Resume and complete
        {
            let mut sim = Simulation::from_checkpoint(&db_path2, "rng_test_sim").unwrap();

            assert_eq!(sim.generation(), 50);

            for _generation in 51..=100 {
                sim.step().unwrap();
            }

            // Get final state
            sim.population().individuals()[0]
                .haplotype1()
                .get(0)
                .unwrap()
                .sequence()
                .as_slice()
                .to_vec()
        }
    };

    // Compare final states - they should be identical due to RNG state preservation
    assert_eq!(
        final_pop_continuous.len(),
        final_pop_resumed.len(),
        "Sequence lengths should match"
    );

    // Count differences
    let differences: usize = final_pop_continuous
        .iter()
        .zip(final_pop_resumed.iter())
        .filter(|(a, b)| a != b)
        .count();

    println!(
        "Differences between continuous and resumed: {} / {}",
        differences,
        final_pop_continuous.len()
    );

    // With proper RNG state preservation, there should be NO differences
    assert_eq!(
        differences, 0,
        "RNG state preservation failed: sequences diverged after resume"
    );

    cleanup_db(&db_path1);
    cleanup_db(&db_path2);
    println!("✅ RNG state preservation test passed!");
}

#[test]
fn test_resume_with_wrong_sim_id_fails() {
    let db_path = test_db_path("wrong_id");
    cleanup_db(&db_path);

    // Create a simulation
    let structure = RepeatStructure::new(Nucleotide::A, 10, 5, 10, 1);

    let mutation = MutationConfig::uniform(0.001).unwrap();
    let recombination = RecombinationConfig::standard(0.01, 0.7, 0.1).unwrap();
    let fitness = FitnessConfig::neutral();
    let config = SimulationConfig::new(5, 50, Some(42));

    let correct_sim_id = "correct_sim";
    let wrong_sim_id = "wrong_sim";

    // Run and record with correct_sim_id
    {
        let mut sim = Simulation::new(
            structure.clone(),
            mutation.clone(),
            recombination.clone(),
            fitness.clone(),
            config.clone(),
        )
        .unwrap();

        let snapshot = SimulationSnapshot {
            structure,
            mutation,
            recombination,
            fitness,
            config,
        };

        let mut recorder =
            Recorder::new(&db_path, correct_sim_id, RecordingStrategy::EveryN(10)).unwrap();

        recorder.record_full_config(&snapshot).unwrap();

        for generation in 1..=20 {
            sim.step().unwrap();

            if recorder.should_record(generation) {
                let rng_state = sim.rng_state_bytes();
                recorder
                    .record_generation(sim.population(), generation)
                    .unwrap();
                recorder.record_checkpoint(generation, &rng_state).unwrap();
            }
        }

        recorder.close().ok();
    }

    // Try to resume with wrong sim_id - should fail
    let result = Simulation::from_checkpoint(&db_path, wrong_sim_id);

    assert!(result.is_err(), "Should fail when using wrong sim_id");
    println!("✓ Correctly rejected wrong sim_id");

    // Resume with correct sim_id - should succeed
    let result = Simulation::from_checkpoint(&db_path, correct_sim_id);
    assert!(result.is_ok(), "Should succeed with correct sim_id");

    cleanup_db(&db_path);
    println!("✅ Wrong sim_id test passed!");
}
