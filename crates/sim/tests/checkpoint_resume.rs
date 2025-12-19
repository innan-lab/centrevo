//! Integration tests for checkpoint and resume functionality.

use centrevo_sim::base::Nucleotide;
use centrevo_sim::simulation::{
    Configuration, EvolutionConfig, ExecutionConfig, FitnessConfig, GenerationMode,
    InitializationConfig, MutationConfig, RecombinationConfig, Simulation, UniformRepeatStructure,
};
use centrevo_sim::storage::{AsyncRecorder, BufferConfig, RecordingStrategy};
use std::path::PathBuf;
use tokio::runtime::Runtime;

/// Helper to create a test database path
fn test_db_path(name: &str) -> PathBuf {
    PathBuf::from(format!("/tmp/centrevo_test_checkpoint_{name}.sqlite"))
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
    let structure = UniformRepeatStructure::new(
        Nucleotide::A,
        10, // Small for testing
        5,
        10,
        1,
    );

    let mutation = MutationConfig::uniform(0.001).unwrap();
    let recombination = RecombinationConfig::standard(0.01, 0.7, 0.1).unwrap();
    let fitness = FitnessConfig::neutral();
    let config = ExecutionConfig::new(10, 100, Some(42)); // 100 generations

    let _sim_id = "test_sim";
    let rt = Runtime::new().unwrap();

    // Part 1: Run simulation for first 50 generations with checkpoints
    rt.block_on(async {
        let sim_config = Configuration {
            execution: config.clone(),
            evolution: EvolutionConfig {
                mutation: mutation.clone(),
                recombination: recombination.clone(),
                fitness: fitness.clone(),
            },
            initialization: InitializationConfig::Generate {
                structure: structure.clone(),
                mode: GenerationMode::Uniform,
            },
        };
        let mut sim = Simulation::new(sim_config.clone()).unwrap();

        // Setup recorder with full config
        // Setup recorder with full config
        // Use sim_config directly

        let strategy = RecordingStrategy::EveryN(10);
        let buffer_config = BufferConfig {
            compression_level: 0,
            ..Default::default()
        };

        let recorder = AsyncRecorder::new(
            &db_path,
            &sim_config,
            buffer_config,
            centrevo_sim::simulation::CodecStrategy::default(),
        )
        .unwrap();

        // Run first 50 generations
        for generation in 1..=50 {
            sim.step().unwrap();

            if strategy.should_record(generation) {
                let rng_state = sim.rng_state_bytes();
                recorder
                    .record_generation(sim.population(), generation, Some(rng_state))
                    .await
                    .unwrap();
                // Checkpoint is implicit in AsyncRecorder?
                // Wait, AsyncRecorder DOES NOT have record_checkpoint method exposed directly?
                // Let's check AsyncRecorder API.
                // It has `record_checkpoint(generation, rng_state)`.
                // But `record_generation` also calls internal buffering.
                // `start` step had `recorder.record_checkpoint`.
                // AsyncRecorder DOES have `record_checkpoint`.
                // I need to confirm.
                // Assuming it has.
                // Wait, in previous task I implemented `record_metadata` and `record_full_config`.
                // `recorder.rs` has `record_checkpoint`?
                // I should reuse logic or double check.
                // If I don't have it, I can't test it.
                // `AsyncRecorder` usually sends `Snapshot` message.
                // `Snapshot` message contains `rng_state`.
                // So `record_generation` is enough?
                // `Recorder` (sync) had separate `record_checkpoint`.
                // `AsyncRecorder::record_generation` takes `rng_state`!
                // So calling `record_generation` IS recording checkpoint data (RNG state).
                // `Recorder::record_generation` took `rng_state`?
                // Sync Recorder `record_generation` signature: `(pop, gen, rng_state)`.
                // So `AsyncRecorder::record_generation` signature matches.
                // Does `recorder.record_checkpoint` exist separately?
                // Code in test calls `recorder.record_checkpoint`.
                // If `AsyncRecorder` doesn't have it, I should remove call?
                // `AsyncRecorder` likely combines them.
                // I'll assume `record_generation` handles both or checks.
                // But wait, the test calls both.
                // Sync recorder: `record_generation` records population. `record_checkpoint` writes rng state to `checkpoints` table.
                // `AsyncRecorder`: `record_generation` takes `rng_state`.
                // So it probably does both.
                // I will NOT call `record_checkpoint` separately if it doesn't exist.
                // I'll check `recorder.rs` quickly if I can...
                // But I'll assume `record_generation` covers it for now and verify.
                // So I remove explicit `record_checkpoint`.
            }
        }

        recorder.close().await.unwrap();

        assert_eq!(sim.generation(), 50);
        println!("✓ Completed first 50 generations");
    });

    // Part 2: Resume from checkpoint and continue
    rt.block_on(async {
        let mut sim = Simulation::from_checkpoint(&db_path).unwrap();

        // Verify we resumed from generation 50 (last checkpoint at gen 50)
        assert_eq!(sim.generation(), 50, "Should resume from generation 50");

        println!("✓ Resumed from generation {}", sim.generation());

        // Setup recorder again
        let buffer_config = BufferConfig {
            compression_level: 0,
            ..Default::default()
        };
        let recorder = AsyncRecorder::new(
            &db_path,
            sim.configuration(),
            buffer_config,
            centrevo_sim::simulation::CodecStrategy::default(),
        )
        .unwrap();
        let strategy = RecordingStrategy::EveryN(10);

        // Continue for remaining 50 generations
        for generation in 51..=100 {
            sim.step().unwrap();

            if strategy.should_record(generation) {
                let rng_state = sim.rng_state_bytes();
                recorder
                    .record_generation(sim.population(), generation, Some(rng_state))
                    .await
                    .unwrap();
                // recorder.record_checkpoint(generation, &rng_state).unwrap();
            }
        }

        recorder.close().await.unwrap();

        assert_eq!(sim.generation(), 100);
        println!("✓ Completed remaining 50 generations");
    });

    // Part 3: Verify final state
    {
        use centrevo_sim::storage::QueryBuilder;

        let query = QueryBuilder::new(&db_path).unwrap();

        // Check recorded generations
        let recorded_gens = query.get_recorded_generations().unwrap();
        assert!(!recorded_gens.is_empty());
        assert_eq!(*recorded_gens.last().unwrap(), 100);

        println!("✓ Recorded generations: {recorded_gens:?}");

        // Verify we can load final population
        let final_pop = query.get_generation(100).unwrap();
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
    let structure = UniformRepeatStructure::new(Nucleotide::A, 10, 5, 10, 1);

    let mutation = MutationConfig::uniform(0.001).unwrap();
    let recombination = RecombinationConfig::standard(0.01, 0.7, 0.1).unwrap();
    let fitness = FitnessConfig::neutral();
    let config = ExecutionConfig::new(5, 100, Some(123));

    let _sim_id = "multi_resume_sim";
    let rt = Runtime::new().unwrap();

    // Run simulation in 4 chunks: 0->25, 25->50, 50->75, 75->100
    let checkpoints = [25, 50, 75, 100];

    for (i, &target_gen) in checkpoints.iter().enumerate() {
        rt.block_on(async {
            let mut sim = if i == 0 {
                // First run: create new simulation
                let sim_config = Configuration {
                    execution: config.clone(),
                    evolution: EvolutionConfig {
                        mutation: mutation.clone(),
                        recombination: recombination.clone(),
                        fitness: fitness.clone(),
                    },
                    initialization: InitializationConfig::Generate {
                        structure: structure.clone(),
                        mode: GenerationMode::Uniform,
                    },
                };
                let sim = Simulation::new(sim_config.clone()).unwrap();

                // Record full config
                // Record full config
                // Use sim_config directly

                let buffer_config = BufferConfig {
                    compression_level: 0,
                    ..Default::default()
                };
                let recorder = AsyncRecorder::new(
                    &db_path,
                    &sim_config,
                    buffer_config,
                    centrevo_sim::simulation::CodecStrategy::default(),
                )
                .unwrap();

                recorder.close().await.unwrap();

                sim
            } else {
                // Subsequent runs: resume from checkpoint
                Simulation::from_checkpoint(&db_path).unwrap()
            };

            let start_gen = sim.generation();
            println!("Run {}: Starting from generation {}", i + 1, start_gen);

            let buffer_config = BufferConfig {
                compression_level: 0,
                ..Default::default()
            };
            let recorder = AsyncRecorder::new(
                &db_path,
                sim.configuration(),
                buffer_config,
                centrevo_sim::simulation::CodecStrategy::default(),
            )
            .unwrap();
            let strategy = RecordingStrategy::EveryN(25);

            // Run to target generation
            for generation in (start_gen + 1)..=target_gen {
                sim.step().unwrap();

                if strategy.should_record(generation) {
                    let rng_state = sim.rng_state_bytes();
                    recorder
                        .record_generation(sim.population(), generation, Some(rng_state))
                        .await
                        .unwrap();
                }
            }

            recorder.close().await.unwrap();

            assert_eq!(sim.generation(), target_gen);
            println!("✓ Reached generation {target_gen}");
        });
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
    let structure = UniformRepeatStructure::new(Nucleotide::A, 10, 5, 10, 1);

    let mutation = MutationConfig::uniform(0.001).unwrap();
    let recombination = RecombinationConfig::standard(0.01, 0.7, 0.1).unwrap();
    let fitness = FitnessConfig::neutral();
    let seed = 999u64;
    let config = ExecutionConfig::new(10, 100, Some(seed));

    // Simulation 1: Run straight through to gen 100
    let final_pop_continuous = {
        let sim_config = Configuration {
            execution: config.clone(),
            evolution: EvolutionConfig {
                mutation: mutation.clone(),
                recombination: recombination.clone(),
                fitness: fitness.clone(),
            },
            initialization: InitializationConfig::Generate {
                structure: structure.clone(),
                mode: GenerationMode::Uniform,
            },
        };
        let mut sim = Simulation::new(sim_config).unwrap();

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

    let rt = Runtime::new().unwrap();

    // Simulation 2: Run to gen 50, checkpoint, resume, then complete
    let final_pop_resumed = {
        // Part 1: Run to gen 50
        rt.block_on(async {
            let sim_config = Configuration {
                execution: config.clone(),
                evolution: EvolutionConfig {
                    mutation: mutation.clone(),
                    recombination: recombination.clone(),
                    fitness: fitness.clone(),
                },
                initialization: InitializationConfig::Generate {
                    structure: structure.clone(),
                    mode: GenerationMode::Uniform,
                },
            };
            let mut sim = Simulation::new(sim_config.clone()).unwrap();

            let strategy = RecordingStrategy::EveryN(50);
            let buffer_config = BufferConfig {
                compression_level: 0,
                ..Default::default()
            };

            let recorder = AsyncRecorder::new(
                &db_path2,
                &sim_config,
                buffer_config,
                centrevo_sim::simulation::CodecStrategy::default(),
            )
            .unwrap();

            for generation in 1..=50 {
                sim.step().unwrap();

                if strategy.should_record(generation) {
                    let rng_state = sim.rng_state_bytes();
                    recorder
                        .record_generation(sim.population(), generation, Some(rng_state))
                        .await
                        .unwrap();
                }
            }

            recorder.close().await.unwrap();
        });

        // Part 2: Resume and complete
        {
            let mut sim = Simulation::from_checkpoint(&db_path2).unwrap();

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
