//! Benchmarks for genome module (chromosome, haplotype, individual operations).
use std::hint::black_box;

use criterion::{criterion_group, criterion_main, Criterion, BenchmarkId, Throughput};
use centrevo::base::{Alphabet, Nucleotide, Sequence};
use centrevo::genome::{Chromosome, Haplotype, Individual};

fn create_test_chromosome(size: usize, alphabet: Alphabet) -> Chromosome {
    let mut seq = Sequence::with_capacity(size, alphabet);
    for i in 0..size {
        seq.push(match i % 4 {
            0 => Nucleotide::A,
            1 => Nucleotide::C,
            2 => Nucleotide::G,
            _ => Nucleotide::T,
        });
    }
    Chromosome::new("chr1", seq, 171, 12)
}

fn create_test_individual(chr_size: usize, n_chromosomes: usize) -> Individual {
    let alphabet = Alphabet::dna();
    let mut hap1 = Haplotype::new();
    let mut hap2 = Haplotype::new();
    
    for _ in 0..n_chromosomes {
        let chr1 = create_test_chromosome(chr_size, alphabet.clone());
        let chr2 = create_test_chromosome(chr_size, alphabet.clone());
        hap1.push(chr1);
        hap2.push(chr2);
    }
    
    Individual::new("ind1", hap1, hap2)
}

/// Benchmark chromosome creation
fn bench_chromosome_creation(c: &mut Criterion) {
    let mut group = c.benchmark_group("chromosome_creation");
    let alphabet = Alphabet::dna();
    let sizes = [1_000, 10_000, 100_000];
    
    for size in sizes {
        group.throughput(Throughput::Elements(size as u64));
        group.bench_with_input(BenchmarkId::new("create", size), &size, |b, &s| {
            b.iter(|| black_box(create_test_chromosome(s, alphabet.clone())));
        });
    }
    
    group.finish();
}

/// Benchmark chromosome operations
fn bench_chromosome_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("chromosome_operations");
    let alphabet = Alphabet::dna();
    let sizes = [1_000, 10_000, 100_000];
    
    for size in sizes {
        let chr = create_test_chromosome(size, alphabet.clone());
        
        group.throughput(Throughput::Elements(size as u64));
        
        group.bench_with_input(BenchmarkId::new("clone", size), &chr, |b, c| {
            b.iter(|| black_box(c.clone()));
        });
        
        group.bench_with_input(BenchmarkId::new("gc_content", size), &chr, |b, c| {
            b.iter(|| black_box(c.gc_content()));
        });
        
        group.bench_with_input(BenchmarkId::new("to_string", size), &chr, |b, c| {
            b.iter(|| black_box(c.to_string()));
        });
        
        group.bench_with_input(BenchmarkId::new("to_formatted_string", size), &chr, |b, c| {
            b.iter(|| black_box(c.to_formatted_string('|', ' ')));
        });
        
        group.bench_with_input(BenchmarkId::new("to_shared", size), &chr, |b, c| {
            b.iter(|| black_box(c.to_shared()));
        });
    }
    
    group.finish();
}

/// Benchmark haplotype operations
fn bench_haplotype_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("haplotype_operations");
    let alphabet = Alphabet::dna();
    let chr_sizes = [1_000, 10_000];
    let n_chromosomes_list = [1, 10, 50];
    
    for chr_size in chr_sizes {
        for n_chromosomes in n_chromosomes_list {
            let mut hap = Haplotype::new();
            for _ in 0..n_chromosomes {
                hap.push(create_test_chromosome(chr_size, alphabet.clone()));
            }
            
            let label = format!("{}chr_x_{}", n_chromosomes, chr_size);
            
            group.bench_with_input(
                BenchmarkId::new("clone", &label),
                &hap,
                |b, h| b.iter(|| black_box(h.clone()))
            );
            
            group.bench_with_input(
                BenchmarkId::new("total_length", &label),
                &hap,
                |b, h| b.iter(|| black_box(h.total_length()))
            );
            
            group.bench_with_input(
                BenchmarkId::new("iterate", &label),
                &hap,
                |b, h| {
                    b.iter(|| {
                        let mut sum = 0;
                        for chr in h.iter() {
                            sum += chr.len();
                        }
                        black_box(sum)
                    })
                }
            );
        }
    }
    
    group.finish();
}

/// Benchmark individual operations
fn bench_individual_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("individual_operations");
    let chr_sizes = [1_000, 10_000];
    let n_chromosomes_list = [1, 10];
    
    for chr_size in chr_sizes {
        for n_chromosomes in n_chromosomes_list {
            let ind = create_test_individual(chr_size, n_chromosomes);
            let label = format!("{}chr_x_{}", n_chromosomes, chr_size);
            
            group.bench_with_input(
                BenchmarkId::new("clone", &label),
                &ind,
                |b, i| b.iter(|| black_box(i.clone()))
            );
            
            group.bench_with_input(
                BenchmarkId::new("access_haplotypes", &label),
                &ind,
                |b, i| {
                    b.iter(|| {
                        black_box(i.haplotype1());
                        black_box(i.haplotype2());
                    })
                }
            );
        }
    }
    
    group.finish();
}

/// Benchmark memory patterns for populations
fn bench_population_memory(c: &mut Criterion) {
    let mut group = c.benchmark_group("population_memory");
    
    group.bench_function("create_100_individuals", |b| {
        b.iter(|| {
            let mut population = Vec::with_capacity(100);
            for _ in 0..100 {
                let ind = create_test_individual(10_000, 1);
                population.push(ind);
            }
            black_box(population)
        });
    });
    
    group.bench_function("clone_100_individuals", |b| {
        let mut population = Vec::with_capacity(100);
        for _ in 0..100 {
            population.push(create_test_individual(10_000, 1));
        }
        
        b.iter(|| {
            let cloned: Vec<_> = population.iter().map(|ind| ind.clone()).collect();
            black_box(cloned)
        });
    });
    
    group.finish();
}

criterion_group!(
    benches,
    bench_chromosome_creation,
    bench_chromosome_operations,
    bench_haplotype_operations,
    bench_individual_operations,
    bench_population_memory,
);

criterion_main!(benches);
