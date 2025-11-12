//! Integration tests for async recording and compression.

use centrevo::{
    base::{Alphabet, Nucleotide},
    genome::{Chromosome, Haplotype, Individual},
    simulation::{
        FitnessConfig, MutationConfig, Population, RecombinationConfig, RepeatStructure,
        Simulation, SimulationConfig,
    },
    storage::{AsyncRecorder, BufferConfig, RecorderStats},
};

fn create_test_individual(id: &str, length: usize) -> Individual {
    let alphabet = Alphabet::dna();
    let chr1 = Chromosome::uniform(
        format!("{}_h1_chr1", id),
        Nucleotide::A,
        length,
        20,
        5,
        alphabet.clone(),
    );
    let chr2 = Chromosome::uniform(
        format!("{}_h2_chr1", id),
        Nucleotide::C,
        length,
        20,
        5,
        alphabet,
    );

    let h1 = Haplotype::from_chromosomes(vec![chr1]);
    let h2 = Haplotype::from_chromosomes(vec![chr2]);

    Individual::new(id, h1, h2)
}

fn create_test_population(size: usize, chr_length: usize) -> Population {
    let mut individuals = Vec::new();
    for i in 0..size {
        let mut ind = create_test_individual(&format!("ind_{}", i), chr_length);
        ind.set_fitness(i as f64 / size as f64);
        individuals.push(ind);
    }
    Population::new("test_pop", individuals)
}

#[tokio::test]
async fn test_full_simulation_with_async_recording() {
    let path = "/tmp/test_full_simulation.sqlite";
    let _ = std::fs::remove_file(path);

    // Create simulation
    let alphabet = Alphabet::dna();
    let structure = RepeatStructure::new(alphabet.clone(), Nucleotide::A, 20, 10, 10, 1);
    let mutation = MutationConfig::uniform(alphabet, 0.001).unwrap();
    let recombination = RecombinationConfig::standard(0.01, 0.7, 0.1).unwrap();
    let fitness = FitnessConfig::neutral();
    let config = SimulationConfig::new(50, 100, Some(42));

    let mut sim = Simulation::new(structure, mutation, recombination, fitness, config).unwrap();

    // Create async recorder
    let buffer_config = BufferConfig::medium();
    let recorder = AsyncRecorder::new(path, "full_sim_test", buffer_config).unwrap();

    // Run simulation and record periodically
    for generation in 0..100 {
        sim.step().unwrap();

        if generation % 10 == 0 {
            let dummy_rng_state = vec![0u8; 32];
            recorder
                .record_generation(sim.population(), generation, dummy_rng_state)
                .await
                .unwrap();
        }
    }

    // Close and get stats
    let stats = recorder.close().await.unwrap();

    // Verify results
    assert_eq!(stats.generations_recorded, 10);
    assert!(stats.bytes_compressed > 0);
    assert!(stats.bytes_written > 0);
    assert!(stats.bytes_written < stats.bytes_compressed);
    assert!(stats.compression_ratio < 0.5); // Should achieve at least 50% compression
    assert_eq!(stats.buffer_full_count, 0); // Buffer should never have been full

    std::fs::remove_file(path).ok();
}

#[tokio::test]
async fn test_async_recording_with_rapid_succession() {
    let path = "/tmp/test_rapid_recording.sqlite";
    let _ = std::fs::remove_file(path);

    let config = BufferConfig {
        capacity: 20, // Large buffer for rapid recording
        compression_level: 5, // Fast compression
        warn_threshold: 0.8,
    };

    let recorder = AsyncRecorder::new(path, "rapid_test", config).unwrap();
    let pop = create_test_population(30, 500);

    // Record 20 generations rapidly
    for generation in 0..20 {
        let dummy_rng_state = vec![0u8; 32];
        recorder
            .record_generation(&pop, generation, dummy_rng_state)
            .await
            .unwrap();
    }

    let stats = recorder.close().await.unwrap();

    assert_eq!(stats.generations_recorded, 20);
    assert!(stats.avg_compression_ms > 0.0);
    assert!(stats.avg_write_ms > 0.0);

    std::fs::remove_file(path).ok();
}

#[tokio::test]
async fn test_buffer_overflow_handling() {
    let path = "/tmp/test_buffer_overflow.sqlite";
    let _ = std::fs::remove_file(path);

    // Create very small buffer with slow compression
    let config = BufferConfig {
        capacity: 2,
        compression_level: 20, // Very slow
        warn_threshold: 0.5,
    };

    let recorder = AsyncRecorder::new(path, "overflow_test", config).unwrap();
    let pop = create_test_population(100, 1000); // Large population

    // Try to fill the buffer quickly
    for generation in 0..5 {
        let dummy_rng_state = vec![0u8; 32];
        recorder
            .record_generation(&pop, generation, dummy_rng_state)
            .await
            .unwrap();
    }

    let stats = recorder.close().await.unwrap();

    // Should have recorded all generations despite small buffer
    assert_eq!(stats.generations_recorded, 5);
    // Just verify the field exists (it's always >= 0 since it's unsigned)
    let _ = stats.buffer_full_count;

    std::fs::remove_file(path).ok();
}

#[tokio::test]
async fn test_compression_ratio_comparison() {
    let mut results: Vec<(i32, RecorderStats)> = Vec::new();

    for level in [1, 10, 20] {
        let path = format!("/tmp/test_compression_level_{}.sqlite", level);
        let _ = std::fs::remove_file(&path);

        let config = BufferConfig {
            capacity: 5,
            compression_level: level,
            warn_threshold: 0.8,
        };

        let recorder = AsyncRecorder::new(&path, "compression_test", config).unwrap();
        let pop = create_test_population(50, 1000);

        let dummy_rng_state = vec![0u8; 32];
        recorder.record_generation(&pop, 0, dummy_rng_state).await.unwrap();

        let stats = recorder.close().await.unwrap();
        results.push((level, stats));

        std::fs::remove_file(&path).ok();
    }

    // Verify that higher levels achieve better compression
    assert!(results[0].1.compression_ratio >= results[1].1.compression_ratio);
    assert!(results[1].1.compression_ratio >= results[2].1.compression_ratio);

    // Note: Timing can vary, so we just log it rather than assert
    // Generally lower levels should be faster, but OS/system factors can affect this

    println!("\nCompression Level Comparison:");
    for (level, stats) in results {
        println!(
            "  Level {}: {:.1}% ratio, {:.2}ms time",
            level,
            stats.compression_ratio * 100.0,
            stats.avg_compression_ms
        );
    }
}

#[tokio::test]
async fn test_memory_usage_estimation() {
    // Test memory estimation accuracy
    let configs = vec![
        (BufferConfig::small(), 50, 1000),
        (BufferConfig::medium(), 100, 5000),
        (BufferConfig::large(), 200, 10000),
    ];

    for (config, pop_size, chr_length) in configs {
        let estimated = config.estimate_memory_usage(pop_size, chr_length);
        let expected = config.capacity * pop_size * chr_length * 2;

        assert_eq!(
            estimated, expected,
            "Memory estimation incorrect for capacity={}, pop={}, chr={}",
            config.capacity, pop_size, chr_length
        );

        println!(
            "Config(cap={}): {}MB for {} ind × {}bp",
            config.capacity,
            estimated / 1_048_576,
            pop_size,
            chr_length
        );
    }
}

#[tokio::test]
async fn test_concurrent_recorders() {
    // Test multiple recorders running concurrently
    let paths: Vec<String> = (0..3)
        .map(|i| format!("/tmp/test_concurrent_{}.sqlite", i))
        .collect();

    // Clean up
    for path in &paths {
        let _ = std::fs::remove_file(path);
    }

    let mut handles = vec![];

    for (i, path) in paths.iter().enumerate() {
        let path_clone = path.clone();
        let handle = tokio::spawn(async move {
            let config = BufferConfig::medium();
            let recorder = AsyncRecorder::new(&path_clone, "concurrent_test", config).unwrap();

            let pop = create_test_population(30, 500);

            for generation in 0..5 {
                let dummy_rng_state = vec![0u8; 32];
                recorder
                    .record_generation(&pop, generation, dummy_rng_state)
                    .await
                    .unwrap();
            }

            let stats = recorder.close().await.unwrap();
            (i, stats)
        });

        handles.push(handle);
    }

    // Wait for all to complete
    for handle in handles {
        let (idx, stats) = handle.await.unwrap();
        assert_eq!(stats.generations_recorded, 5);
        println!("Recorder {}: compression={:.1}%", idx, stats.compression_ratio * 100.0);
    }

    // Clean up
    for path in &paths {
        std::fs::remove_file(path).ok();
    }
}

#[tokio::test]
async fn test_snapshot_consistency() {
    let path = "/tmp/test_snapshot_consistency.sqlite";
    let _ = std::fs::remove_file(path);

    let config = BufferConfig::small();
    let recorder = AsyncRecorder::new(path, "consistency_test", config).unwrap();

    let pop = create_test_population(20, 200);

    // Record same population twice
    let dummy_rng_state = vec![0u8; 32];
    recorder.record_generation(&pop, 0, dummy_rng_state.clone()).await.unwrap();
    recorder.record_generation(&pop, 1, dummy_rng_state).await.unwrap();

    let stats = recorder.close().await.unwrap();

    assert_eq!(stats.generations_recorded, 2);
    // Both generations should compress similarly (same data)
    assert!(stats.compression_ratio < 0.3);

    std::fs::remove_file(path).ok();
}

#[tokio::test]
async fn test_performance_vs_sync() {
    use centrevo::storage::{Recorder, RecordingStrategy};
    use std::time::Instant;

    let sync_path = "/tmp/test_sync_perf.sqlite";
    let async_path = "/tmp/test_async_perf.sqlite";

    let _ = std::fs::remove_file(sync_path);
    let _ = std::fs::remove_file(async_path);

    let pop = create_test_population(50, 1000);

    // Benchmark sync recording
    let sync_start = Instant::now();
    {
        let mut recorder = Recorder::new(sync_path, "sync_test", RecordingStrategy::All).unwrap();
        for generation in 0..10 {
            recorder.record_generation(&pop, generation).unwrap();
        }
        recorder.close().unwrap();
    }
    let sync_time = sync_start.elapsed();

    // Benchmark async recording
    let async_start = Instant::now();
    {
        let config = BufferConfig::medium();
        let recorder = AsyncRecorder::new(async_path, "async_test", config).unwrap();
        for generation in 0..10 {
            let dummy_rng_state = vec![0u8; 32];
            recorder.record_generation(&pop, generation, dummy_rng_state).await.unwrap();
        }
        recorder.close().await.unwrap();
    }
    let async_time = async_start.elapsed();

    println!("\nPerformance Comparison (10 generations, 50 ind × 1kb):");
    println!("  Sync:  {:.2}s", sync_time.as_secs_f64());
    println!("  Async: {:.2}s", async_time.as_secs_f64());
    println!(
        "  Speedup: {:.1}x",
        sync_time.as_secs_f64() / async_time.as_secs_f64()
    );

    std::fs::remove_file(sync_path).ok();
    std::fs::remove_file(async_path).ok();
}
