use centrevo::{
    base::Nucleotide,
    genome::{Chromosome, Haplotype, Individual},
    simulation::Population,
    storage::{AsyncRecorder, BufferConfig, Recorder, RecordingStrategy},
};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use std::hint::black_box;

fn create_test_individual(id: &str, length: usize) -> Individual {
    
    let chr1 = Chromosome::uniform(
        format!("{id}_h1_chr1"),
        Nucleotide::A,
        length,
        20,
        5,
    );
    let chr2 = Chromosome::uniform(
        format!("{id}_h2_chr1"),
        Nucleotide::C,
        length,
        20,
        5,
    );

    let h1 = Haplotype::from_chromosomes(vec![chr1]);
    let h2 = Haplotype::from_chromosomes(vec![chr2]);

    Individual::new(id, h1, h2)
}

fn create_test_population(size: usize, chr_length: usize) -> Population {
    let mut individuals = Vec::new();
    for i in 0..size {
        let mut ind = create_test_individual(&format!("ind_{i}"), chr_length);
        ind.set_fitness(i as f64 / size as f64);
        individuals.push(ind);
    }
    Population::new("test_pop", individuals)
}

fn dummy_rng_state() -> Vec<u8> {
    vec![0u8; 32]  // Dummy RNG state for benchmarking
}

fn bench_sync_recording(c: &mut Criterion) {
    let mut group = c.benchmark_group("sync_recording");

    for pop_size in [10, 50, 100, 200] {
        group.throughput(Throughput::Elements(pop_size as u64));

        group.bench_with_input(
            BenchmarkId::from_parameter(format!("{}ind_1kb", pop_size)),
            &pop_size,
            |b, &size| {
                b.iter(|| {
                    let path = format!("/tmp/bench_sync_{}.sqlite", size);
                    let _ = std::fs::remove_file(&path);

                    let mut recorder = Recorder::new(&path, "bench_sim", RecordingStrategy::All)
                        .expect("Failed to create recorder");

                    let pop = create_test_population(size, 1000);

                    recorder
                        .record_generation(&pop, 0)
                        .expect("Failed to record");

                    recorder.close().expect("Failed to close");
                    std::fs::remove_file(&path).ok();

                    black_box(());
                });
            },
        );
    }

    group.finish();
}

fn bench_async_recording(c: &mut Criterion) {
    let mut group = c.benchmark_group("async_recording");

    let runtime = tokio::runtime::Runtime::new().unwrap();

    for pop_size in [10, 50, 100, 200] {
        group.throughput(Throughput::Elements(pop_size as u64));

        group.bench_with_input(
            BenchmarkId::from_parameter(format!("{}ind_1kb", pop_size)),
            &pop_size,
            |b, &size| {
                b.iter(|| {
                    runtime.block_on(async {
                        let path = format!("/tmp/bench_async_{}.sqlite", size);
                        let _ = std::fs::remove_file(&path);

                        let config = BufferConfig::medium();
                        let recorder = AsyncRecorder::new(&path, "bench_sim", config)
                            .expect("Failed to create recorder");

                        let pop = create_test_population(size, 1000);

                        recorder
                            .record_generation(&pop, 0, dummy_rng_state())
                            .await
                            .expect("Failed to record");

                        recorder.close().await.expect("Failed to close");
                        std::fs::remove_file(&path).ok();

                        black_box(());
                    });
                });
            },
        );
    }

    group.finish();
}

fn bench_compression_levels(c: &mut Criterion) {
    let mut group = c.benchmark_group("compression_levels");

    let runtime = tokio::runtime::Runtime::new().unwrap();
    let pop = create_test_population(50, 1000);

    for level in [1, 3, 5, 10, 15, 20] {
        group.bench_with_input(
            BenchmarkId::from_parameter(format!("level_{}", level)),
            &level,
            |b, &lvl| {
                b.iter(|| {
                    runtime.block_on(async {
                        let path = format!("/tmp/bench_compress_{}.sqlite", lvl);
                        let _ = std::fs::remove_file(&path);

                        let config = BufferConfig {
                            capacity: 5,
                            compression_level: lvl,
                            warn_threshold: 0.8,
                        };

                        let recorder = AsyncRecorder::new(&path, "bench_sim", config)
                            .expect("Failed to create recorder");

                        recorder
                            .record_generation(&pop, 0, dummy_rng_state())
                            .await
                            .expect("Failed to record");

                        recorder.close().await.expect("Failed to close");
                        std::fs::remove_file(&path).ok();

                        black_box(());
                    });
                });
            },
        );
    }

    group.finish();
}

fn bench_population_sizes(c: &mut Criterion) {
    let mut group = c.benchmark_group("population_sizes");

    let runtime = tokio::runtime::Runtime::new().unwrap();

    for (pop_size, chr_length) in [(10, 10000), (50, 10000), (100, 10000), (200, 5000)] {
        let data_size = pop_size * chr_length * 2; // 2 haplotypes
        group.throughput(Throughput::Bytes(data_size as u64));

        group.bench_with_input(
            BenchmarkId::from_parameter(format!("{}ind_{}kb", pop_size, chr_length / 1000)),
            &(pop_size, chr_length),
            |b, &(size, length)| {
                b.iter(|| {
                    runtime.block_on(async {
                        let path = format!("/tmp/bench_size_{}_{}.sqlite", size, length);
                        let _ = std::fs::remove_file(&path);

                        let config = if size < 100 {
                            BufferConfig::small()
                        } else {
                            BufferConfig::medium()
                        };

                        let recorder = AsyncRecorder::new(&path, "bench_sim", config)
                            .expect("Failed to create recorder");

                        let pop = create_test_population(size, length);

                        recorder
                            .record_generation(&pop, 0, dummy_rng_state())
                            .await
                            .expect("Failed to record");

                        recorder.close().await.expect("Failed to close");
                        std::fs::remove_file(&path).ok();

                        black_box(());
                    });
                });
            },
        );
    }

    group.finish();
}

fn bench_multiple_generations(c: &mut Criterion) {
    let mut group = c.benchmark_group("multiple_generations");

    let runtime = tokio::runtime::Runtime::new().unwrap();

    for num_gens in [1, 5, 10, 20] {
        group.throughput(Throughput::Elements(num_gens as u64));

        group.bench_with_input(
            BenchmarkId::from_parameter(format!("{}gens", num_gens)),
            &num_gens,
            |b, &gens| {
                b.iter(|| {
                    runtime.block_on(async {
                        let path = format!("/tmp/bench_multigens_{}.sqlite", gens);
                        let _ = std::fs::remove_file(&path);

                        let config = BufferConfig::medium();
                        let recorder = AsyncRecorder::new(&path, "bench_sim", config)
                            .expect("Failed to create recorder");

                        let pop = create_test_population(50, 1000);

                        for generation in 0..gens {
                            recorder
                                .record_generation(&pop, generation, dummy_rng_state())
                                .await
                                .expect("Failed to record");
                        }

                        recorder.close().await.expect("Failed to close");
                        std::fs::remove_file(&path).ok();

                        black_box(());
                    });
                });
            },
        );
    }

    group.finish();
}

fn bench_buffer_configurations(c: &mut Criterion) {
    let mut group = c.benchmark_group("buffer_configurations");

    let runtime = tokio::runtime::Runtime::new().unwrap();
    let pop = create_test_population(100, 1000);

    for (name, config) in [
        ("small", BufferConfig::small()),
        ("medium", BufferConfig::medium()),
        ("large", BufferConfig::large()),
    ] {
        group.bench_with_input(BenchmarkId::from_parameter(name), &config, |b, cfg| {
            b.iter(|| {
                runtime.block_on(async {
                    let path = format!("/tmp/bench_buffer_{}.sqlite", name);
                    let _ = std::fs::remove_file(&path);

                    let recorder = AsyncRecorder::new(&path, "bench_sim", cfg.clone())
                        .expect("Failed to create recorder");

                    for generation in 0..5 {
                        recorder
                            .record_generation(&pop, generation, dummy_rng_state())
                            .await
                            .expect("Failed to record");
                    }

                    recorder.close().await.expect("Failed to close");
                    std::fs::remove_file(&path).ok();

                    black_box(());
                });
            });
        });
    }

    group.finish();
}

fn bench_snapshot_creation(c: &mut Criterion) {
    let mut group = c.benchmark_group("snapshot_creation");

    for pop_size in [10, 50, 100, 200] {
        group.throughput(Throughput::Elements(pop_size as u64));

        group.bench_with_input(
            BenchmarkId::from_parameter(format!("{}ind", pop_size)),
            &pop_size,
            |b, &size| {
                let pop = create_test_population(size, 1000);
                b.iter(|| {
                    use centrevo::storage::IndividualSnapshot;
                    use rayon::prelude::*;

                    let snapshots: Vec<IndividualSnapshot> = pop
                        .individuals()
                        .par_iter()
                        .map(IndividualSnapshot::from_individual)
                        .collect();

                    black_box(snapshots);
                });
            },
        );
    }

    group.finish();
}

fn bench_compression_only(c: &mut Criterion) {
    let mut group = c.benchmark_group("compression_only");

    use centrevo::storage::IndividualSnapshot;
    use rayon::prelude::*;

    for pop_size in [10, 50, 100] {
        group.throughput(Throughput::Elements(pop_size as u64));

        let pop = create_test_population(pop_size, 1000);
        let snapshots: Vec<IndividualSnapshot> = pop
            .individuals()
            .par_iter()
            .map(IndividualSnapshot::from_individual)
            .collect();

        for level in [1, 10, 20] {
            group.bench_with_input(
                BenchmarkId::from_parameter(format!("{}ind_lvl{}", pop_size, level)),
                &level,
                |b, &lvl| {
                    b.iter(|| {
                        let compressed: Vec<Vec<u8>> = snapshots
                            .par_iter()
                            .flat_map(|s| {
                                vec![
                                    zstd::bulk::compress(&s.haplotype1_seq, lvl).unwrap(),
                                    zstd::bulk::compress(&s.haplotype2_seq, lvl).unwrap(),
                                ]
                            })
                            .collect();

                        black_box(compressed);
                    });
                },
            );
        }
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_sync_recording,
    bench_async_recording,
    bench_compression_levels,
    bench_population_sizes,
    bench_multiple_generations,
    bench_buffer_configurations,
    bench_snapshot_creation,
    bench_compression_only,
);
criterion_main!(benches);
