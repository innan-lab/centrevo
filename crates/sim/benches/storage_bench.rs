use centrevo_sim::base::{FitnessValue, Nucleotide};
use centrevo_sim::genome::{Chromosome, Haplotype, Individual};
use centrevo_sim::simulation::Population;

use criterion::{Criterion, Throughput, black_box, criterion_group, criterion_main};
use tempfile::NamedTempFile;

fn create_test_population(size: usize, chr_length: usize) -> Population {
    let mut individuals = Vec::with_capacity(size);
    // Convert total length to num_hors: ru_length=20, rus_per_hor=5, so one HOR = 100 bp
    let num_hors = chr_length / 100;

    for i in 0..size {
        let chr1 = Chromosome::uniform(format!("ind_{i}_h1_chr1"), Nucleotide::A, 20, 5, num_hors);
        let chr2 = Chromosome::uniform(format!("ind_{i}_h2_chr1"), Nucleotide::C, 20, 5, num_hors);

        let h1 = Haplotype::from_chromosomes(vec![chr1]);
        let h2 = Haplotype::from_chromosomes(vec![chr2]);

        let mut ind = Individual::new(format!("ind_{i}"), h1, h2);
        ind.set_cached_fitness(FitnessValue::new(i as f64 / size as f64));
        individuals.push(ind);
    }
    Population::new("test_pop", individuals)
}

fn bench_recorder_write(c: &mut Criterion) {
    let mut group = c.benchmark_group("recorder_write");
    let pop_size = 50;
    let chr_length = 10_000;
    let rt = tokio::runtime::Runtime::new().unwrap();

    group.throughput(Throughput::Elements(pop_size as u64));

    group.bench_function("record_generation", |b| {
        b.iter_batched(
            || {
                let file = NamedTempFile::new().unwrap();
                let path = file.path().to_owned();
                // Create recorder (returns future?) No, creation is synchronous but interacting is async?
                // AsyncRecorder setup spawns a background task.
                let buffer_config = centrevo_sim::storage::BufferConfig {
                    compression_level: 0,
                    ..Default::default()
                };

                // We need to establish recorder inside the runtime context if it spawns tasks?
                // AsyncRecorder::new spawns a naked thread currently? Or tokio::spawn?
                // It uses `tokio::spawn`. So we need to be in runtime context.
                // Create dummy config for bench
                let config = centrevo_sim::simulation::Configuration {
                    execution: centrevo_sim::simulation::ExecutionConfig::new(pop_size, 100, None),
                    evolution: centrevo_sim::simulation::EvolutionConfig {
                        mutation: centrevo_sim::simulation::MutationConfig::uniform(0.0).unwrap(),
                        recombination: centrevo_sim::simulation::RecombinationConfig::standard(
                            0.0, 0.0, 0.0,
                        )
                        .unwrap(),
                        fitness: centrevo_sim::simulation::FitnessConfig::neutral(),
                    },
                    initialization: centrevo_sim::simulation::InitializationConfig::Generate {
                        structure: centrevo_sim::simulation::UniformRepeatStructure::new(
                            centrevo_sim::base::Nucleotide::A,
                            20,
                            5,
                            1,
                            1,
                        ),
                        mode: centrevo_sim::simulation::GenerationMode::Uniform,
                    },
                };

                let _guard = rt.enter();
                let recorder = centrevo_sim::storage::AsyncRecorder::new(
                    &path,
                    &config,
                    buffer_config,
                    centrevo_sim::simulation::CodecStrategy::default(),
                )
                .unwrap();
                let pop = create_test_population(pop_size, chr_length);
                (file, recorder, pop)
            },
            |(_file, recorder, pop)| {
                // Perform async recording
                rt.block_on(async {
                    let dummy_rng = vec![0u8; 32];
                    recorder
                        .record_generation(
                            black_box(&pop),
                            black_box(0),
                            black_box(Some(dummy_rng)),
                        )
                        .await
                        .unwrap();
                    // We should close to ensure flush in benchmark?
                    // Previous benchmark didn't close, just dropped.
                    // AsyncRecorder drop sends Shutdown message.
                    // But we want to measure `record_generation` speed mainly.
                    // If we don't await completion of background task, we are just measuring channel send speed.
                    // AsyncRecorder::record_generation awaits the send.
                    // So we are measuring serialization + channel send.
                    // The background task processing happens concurrently.
                    // To measure strict throughput, we might want to wait?
                    // But `record_generation` returns when message is sent.
                    // That's the API performance characteristic.
                });
            },
            criterion::BatchSize::SmallInput,
        )
    });

    group.finish();
}

criterion_group!(benches, bench_recorder_write);
criterion_main!(benches);
