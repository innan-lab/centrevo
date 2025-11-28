use centrevo::analysis::{
    nucleotide_diversity, tajimas_d, wattersons_theta, haplotype_diversity,
    linkage_disequilibrium, pairwise_distances, gc_content,
};
use centrevo::base::{Nucleotide, Sequence};
use centrevo::genome::{Chromosome, Haplotype, Individual};
use centrevo::simulation::Population;
use criterion::{criterion_group, criterion_main, Criterion, BenchmarkId};
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;
use std::hint::black_box;

fn create_random_individual(id: &str, length: usize, rng: &mut Xoshiro256PlusPlus) -> Individual {

    let nucleotides = [Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T];

    let mut seq = Sequence::with_capacity(length);
    for _ in 0..length {
        let nuc = nucleotides[rng.random_range(0..4)];
        seq.push(nuc);
    }

    let chr = Chromosome::new(format!("chr_{id}"), seq.clone(), 171, 12);
    let mut hap1 = Haplotype::new();
    hap1.push(chr.clone());
    let mut hap2 = Haplotype::new();
    hap2.push(chr);

    Individual::new(id, hap1, hap2)
}

fn create_test_population(n_individuals: usize, seq_length: usize) -> Population {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let individuals: Vec<Individual> = (0..n_individuals)
        .map(|i| create_random_individual(&format!("ind{i}"), seq_length, &mut rng))
        .collect();

    Population::new("test_pop", individuals)
}

fn bench_nucleotide_diversity(c: &mut Criterion) {
    let mut group = c.benchmark_group("nucleotide_diversity");

    for &(n_ind, seq_len) in &[(10, 1000), (100, 1000), (100, 10000)] {
        let pop = create_test_population(n_ind, seq_len);
        group.bench_with_input(
            BenchmarkId::from_parameter(format!("{n_ind}ind_{seq_len}bp")),
            &pop,
            |b, pop| {
                b.iter(|| {
                    black_box(nucleotide_diversity(pop, 0))
                });
            },
        );
    }

    group.finish();
}

fn bench_tajimas_d(c: &mut Criterion) {
    let mut group = c.benchmark_group("tajimas_d");

    for &(n_ind, seq_len) in &[(10, 1000), (50, 1000), (100, 1000)] {
        let pop = create_test_population(n_ind, seq_len);
        group.bench_with_input(
            BenchmarkId::from_parameter(format!("{n_ind}ind_{seq_len}bp")),
            &pop,
            |b, pop| {
                b.iter(|| {
                    black_box(tajimas_d(pop, 0))
                });
            },
        );
    }

    group.finish();
}

fn bench_wattersons_theta(c: &mut Criterion) {
    let mut group = c.benchmark_group("wattersons_theta");

    for &(n_ind, seq_len) in &[(10, 1000), (50, 1000), (100, 1000)] {
        let pop = create_test_population(n_ind, seq_len);
        group.bench_with_input(
            BenchmarkId::from_parameter(format!("{n_ind}ind_{seq_len}bp")),
            &pop,
            |b, pop| {
                b.iter(|| {
                    black_box(wattersons_theta(pop, 0))
                });
            },
        );
    }

    group.finish();
}

fn bench_haplotype_diversity(c: &mut Criterion) {
    let mut group = c.benchmark_group("haplotype_diversity");

    for &(n_ind, seq_len) in &[(10, 1000), (50, 1000), (100, 1000)] {
        let pop = create_test_population(n_ind, seq_len);
        group.bench_with_input(
            BenchmarkId::from_parameter(format!("{n_ind}ind_{seq_len}bp")),
            &pop,
            |b, pop| {
                b.iter(|| {
                    black_box(haplotype_diversity(pop, 0))
                });
            },
        );
    }

    group.finish();
}

fn bench_linkage_disequilibrium(c: &mut Criterion) {
    let mut group = c.benchmark_group("linkage_disequilibrium");

    for &(n_ind, seq_len) in &[(10, 1000), (50, 1000), (100, 10000)] {
        let pop = create_test_population(n_ind, seq_len);
        group.bench_with_input(
            BenchmarkId::from_parameter(format!("{n_ind}ind_{seq_len}bp")),
            &pop,
            |b, pop| {
                b.iter(|| {
                    black_box(linkage_disequilibrium(pop, 100, 500, 0, 0))
                });
            },
        );
    }

    group.finish();
}

fn bench_pairwise_distances(c: &mut Criterion) {
    let mut group = c.benchmark_group("pairwise_distances");

    for &(n_ind, seq_len) in &[(10, 1000), (50, 1000), (100, 1000)] {
        let pop = create_test_population(n_ind, seq_len);
        group.bench_with_input(
            BenchmarkId::from_parameter(format!("{n_ind}ind_{seq_len}bp")),
            &pop,
            |b, pop| {
                b.iter(|| {
                    black_box(pairwise_distances(pop, 0))
                });
            },
        );
    }

    group.finish();
}

fn bench_gc_content(c: &mut Criterion) {
    let mut group = c.benchmark_group("gc_content");

    for &seq_len in &[1000, 10000, 100000] {
        let pop = create_test_population(10, seq_len);
        group.bench_with_input(
            BenchmarkId::from_parameter(format!("{seq_len}bp_chromosome")),
            &pop,
            |b, pop| {
                b.iter(|| {
                    black_box(gc_content(pop, Some(0), Some(0), Some(0)))
                });
            },
        );
        group.bench_with_input(
            BenchmarkId::from_parameter(format!("{seq_len}bp_population")),
            &pop,
            |b, pop| {
                b.iter(|| {
                    black_box(gc_content(pop, None, None, None))
                });
            },
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_nucleotide_diversity,
    bench_tajimas_d,
    bench_wattersons_theta,
    bench_haplotype_diversity,
    bench_linkage_disequilibrium,
    bench_pairwise_distances,
    bench_gc_content,
);
criterion_main!(benches);
