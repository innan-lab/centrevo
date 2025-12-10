use centrevo_sim::base::Sequence;
use centrevo_sim::evolution::{IndelModel, RecombinationModel, SubstitutionModel};
use centrevo_sim::genome::{Chromosome, RepeatMap};
use criterion::{BenchmarkId, Criterion, Throughput, black_box, criterion_group, criterion_main};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
use std::str::FromStr;

fn bench_mutation(c: &mut Criterion) {
    let mut group = c.benchmark_group("mutation");
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

    let rates = [1e-6, 1e-3, 0.1];
    let lengths = [1_000, 100_000, 1_000_000];

    for &len in &lengths {
        let seq_str = "A".repeat(len);
        let sequence = Sequence::from_str(&seq_str).unwrap();

        group.throughput(Throughput::Elements(len as u64));

        for &rate in &rates {
            // Optimized Uniform model
            let uniform_model = SubstitutionModel::uniform(rate).unwrap();

            // General model with same uniform rates
            let r = rate / 3.0;
            let matrix = [
                [0.0, r, r, r],
                [r, 0.0, r, r],
                [r, r, 0.0, r],
                [r, r, r, 0.0],
            ];
            let general_model = SubstitutionModel::new(matrix).unwrap();

            let parameter_string = format!("len={len}/rate={rate}");

            // 1. Uniform Variant Direct
            group.bench_with_input(
                BenchmarkId::new("uniform_direct", &parameter_string),
                &(len, rate),
                |b, _| {
                    b.iter_batched(
                        || sequence.clone(),
                        |mut seq| uniform_model.mutate_sequence(&mut seq, &mut rng),
                        criterion::BatchSize::SmallInput,
                    )
                },
            );

            // 2. Uniform Variant Sparse
            group.bench_with_input(
                BenchmarkId::new("uniform_sparse", &parameter_string),
                &(len, rate),
                |b, _| {
                    b.iter_batched(
                        || sequence.clone(),
                        |mut seq| uniform_model.mutate_sequence_sparse(&mut seq, &mut rng),
                        criterion::BatchSize::SmallInput,
                    )
                },
            );

            // 3. General Variant Direct
            group.bench_with_input(
                BenchmarkId::new("general_direct", &parameter_string),
                &(len, rate),
                |b, _| {
                    b.iter_batched(
                        || sequence.clone(),
                        |mut seq| general_model.mutate_sequence(&mut seq, &mut rng),
                        criterion::BatchSize::SmallInput,
                    )
                },
            );

            // 4. General Variant Sparse
            group.bench_with_input(
                BenchmarkId::new("general_sparse", &parameter_string),
                &(len, rate),
                |b, _| {
                    b.iter_batched(
                        || sequence.clone(),
                        |mut seq| general_model.mutate_sequence_sparse(&mut seq, &mut rng),
                        criterion::BatchSize::SmallInput,
                    )
                },
            );
        }
    }

    group.finish();
}

fn bench_indels(c: &mut Criterion) {
    let mut group = c.benchmark_group("indels");
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

    let lengths = [1_000, 100_000, 1_000_000];
    let rates = [1e-6, 1e-3, 0.1];
    let length_params = [0.01, 0.1, 0.5]; // p values: 0.01 (long), 0.1 (med), 0.5 (short)

    for &len in &lengths {
        let seq_str = "A".repeat(len);
        let sequence = Sequence::from_str(&seq_str).unwrap();

        group.throughput(Throughput::Elements(len as u64));

        for &rate in &rates {
            for &p in &length_params {
                // Use same rate for insertion and deletion for benchmarking simplicity
                let model = IndelModel::new(rate, rate, p).unwrap();
                let parameter_string = format!("len={len}/rate={rate}/p={p}");

                group.bench_with_input(
                    BenchmarkId::new("apply_indels", &parameter_string),
                    &(len, rate, p),
                    |b, _| {
                        b.iter_batched(
                            || sequence.clone(),
                            |mut seq| model.apply_indels(&mut seq, &mut rng),
                            criterion::BatchSize::SmallInput,
                        )
                    },
                );
            }
        }
    }

    group.finish();
}

fn bench_recombination(c: &mut Criterion) {
    let mut group = c.benchmark_group("recombination");
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

    let break_probs = [1e-6, 1e-3, 0.1];
    let crossover_probs = [0.1, 0.5, 0.9];
    let seq_len = 100_000; // Fixed length for recombination bench

    // Setup chromosomes
    let seq_str = "A".repeat(seq_len);
    let seq = Sequence::from_str(&seq_str).unwrap();
    let map = RepeatMap::uniform(1, 1, seq_len);
    let chr1 = Chromosome::new("chr1", seq.clone(), map.clone());
    let chr2 = Chromosome::new("chr2", seq, map);

    for &break_prob in &break_probs {
        for &crossover_prob in &crossover_probs {
            let model = RecombinationModel::builder()
                .break_prob(break_prob)
                .crossover_prob(crossover_prob)
                .build()
                .unwrap();

            let parameter_string = format!("break={break_prob}/cross={crossover_prob}");

            group.bench_with_input(
                BenchmarkId::new("sample_events", &parameter_string),
                &(break_prob, crossover_prob),
                |b, _| {
                    b.iter(|| {
                        black_box(model.sample_events(
                            black_box(&chr1),
                            black_box(&chr2),
                            &mut rng,
                        ));
                    })
                },
            );
        }
    }

    group.finish();
}

criterion_group!(benches, bench_mutation, bench_indels, bench_recombination);
criterion_main!(benches);
