//! Benchmarks for evolution module (mutation, recombination, selection operations).
use std::hint::black_box;
use criterion::{criterion_group, criterion_main, Criterion, BenchmarkId, Throughput};
use centrevo::base::{Alphabet, Nucleotide, Sequence};
use centrevo::genome::{Chromosome, Haplotype, Individual};
use centrevo::evolution::{
    SubstitutionModel, RecombinationParams, 
    GCContentFitness, LengthFitness, SequenceSimilarityFitness,
    HaplotypeFitness, IndividualFitness,
};
use rand::SeedableRng;
use rand::rngs::StdRng;

fn create_test_sequence(size: usize, alphabet: Alphabet) -> Sequence {
    let mut seq = Sequence::with_capacity(size, alphabet);
    for i in 0..size {
        seq.push(match i % 4 {
            0 => Nucleotide::A,
            1 => Nucleotide::C,
            2 => Nucleotide::G,
            _ => Nucleotide::T,
        });
    }
    seq
}

fn create_test_individual(size: usize) -> Individual {
    let alphabet = Alphabet::dna();
    let seq = create_test_sequence(size, alphabet.clone());
    let chr1 = Chromosome::new("chr1", seq.clone(), 171, 12);
    let chr2 = Chromosome::new("chr2", seq, 171, 12);
    
    let mut hap1 = Haplotype::new();
    let mut hap2 = Haplotype::new();
    hap1.push(chr1);
    hap2.push(chr2);
    
    Individual::new("ind1", hap1, hap2)
}

/// Benchmark mutation operations
fn bench_mutation(c: &mut Criterion) {
    let mut group = c.benchmark_group("mutation");
    let alphabet = Alphabet::dna();
    let mut rng = StdRng::seed_from_u64(42);
    let rates = [0.001, 0.01, 0.1];
    let sizes = [1_000, 10_000, 100_000];
    
    for rate in rates {
        let model = SubstitutionModel::uniform(alphabet.clone(), rate).unwrap();
        
        for size in sizes {
            let label = format!("rate_{}_size_{}", rate, size);
            group.throughput(Throughput::Elements(size as u64));
            
            // Standard mutation
            group.bench_with_input(
                BenchmarkId::new("standard", &label),
                &size,
                |b, &s| {
                    b.iter(|| {
                        let mut seq = create_test_sequence(s, alphabet.clone());
                        model.mutate_sequence(&mut seq, &mut rng);
                        black_box(seq)
                    });
                }
            );
            
            // Poisson mutation
            group.bench_with_input(
                BenchmarkId::new("poisson", &label),
                &size,
                |b, &s| {
                    b.iter(|| {
                        let mut seq = create_test_sequence(s, alphabet.clone());
                        model.mutate_sequence_poisson(&mut seq, &mut rng);
                        black_box(seq)
                    });
                }
            );
        }
    }
    
    group.finish();
}

/// Benchmark single base mutation
fn bench_mutate_base(c: &mut Criterion) {
    let mut group = c.benchmark_group("mutate_base");
    let alphabet = Alphabet::dna();
    let mut rng = StdRng::seed_from_u64(42);
    let rates = [0.001, 0.01, 0.1];
    
    for rate in rates {
        let model = SubstitutionModel::uniform(alphabet.clone(), rate).unwrap();
        
        group.bench_with_input(
            BenchmarkId::new("single", rate),
            &rate,
            |b, _| {
                b.iter(|| {
                    black_box(model.mutate_base(Nucleotide::A, &mut rng))
                });
            }
        );
    }
    
    group.finish();
}

/// Benchmark recombination operations
fn bench_recombination(c: &mut Criterion) {
    let mut group = c.benchmark_group("recombination");
    let alphabet = Alphabet::dna();
    let mut rng = StdRng::seed_from_u64(42);
    let sizes = [1_000, 10_000, 100_000];
    
    let params = RecombinationParams::new(0.01, 0.7, 0.1).unwrap();
    
    for size in sizes {
        group.throughput(Throughput::Elements(size as u64));
        
        // Benchmark event sampling
        group.bench_with_input(
            BenchmarkId::new("sample_event", size),
            &size,
            |b, &s| {
                b.iter(|| black_box(params.sample_event(s, &mut rng)));
            }
        );
        
        // Benchmark crossover
        let seq1 = create_test_sequence(size, alphabet.clone());
        let seq2 = create_test_sequence(size, alphabet.clone());
        let position = size / 2;
        
        group.bench_with_input(
            BenchmarkId::new("crossover", size),
            &(seq1.clone(), seq2.clone()),
            |b, (s1, s2)| {
                b.iter(|| {
                    black_box(params.crossover(s1, s2, position).unwrap())
                });
            }
        );
        
        // Benchmark gene conversion
        let tract_len = 100.min(size / 10);
        group.bench_with_input(
            BenchmarkId::new("gene_conversion", size),
            &(seq1.clone(), seq2.clone()),
            |b, (s1, s2)| {
                b.iter(|| {
                    black_box(params.gene_conversion(s1, s2, position, position + tract_len).unwrap())
                });
            }
        );
    }
    
    group.finish();
}

/// Benchmark fitness calculations
fn bench_fitness(c: &mut Criterion) {
    let mut group = c.benchmark_group("fitness");
    let alphabet = Alphabet::dna();
    let sizes = [1_000, 10_000, 100_000];
    
    // GC content fitness
    let gc_fitness = GCContentFitness::new(0.5, 2.0).unwrap();
    
    for size in sizes {
        let chr = Chromosome::new(
            "chr1",
            create_test_sequence(size, alphabet.clone()),
            171,
            12,
        );
        
        group.throughput(Throughput::Elements(size as u64));
        
        group.bench_with_input(
            BenchmarkId::new("gc_content", size),
            &chr,
            |b, c| {
                b.iter(|| black_box(gc_fitness.haplotype_fitness(c)));
            }
        );
    }
    
    // Length fitness
    let len_fitness = LengthFitness::new(10_000, 0.5).unwrap();
    
    for size in sizes {
        let chr = Chromosome::new(
            "chr1",
            create_test_sequence(size, alphabet.clone()),
            171,
            12,
        );
        
        group.bench_with_input(
            BenchmarkId::new("length", size),
            &chr,
            |b, c| {
                b.iter(|| black_box(len_fitness.haplotype_fitness(c)));
            }
        );
    }
    
    // Sequence similarity fitness
    let sim_fitness = SequenceSimilarityFitness::new(1.0).unwrap();
    
    for size in sizes {
        let ind = create_test_individual(size);
        
        group.bench_with_input(
            BenchmarkId::new("similarity", size),
            &ind,
            |b, i| {
                b.iter(|| black_box(sim_fitness.individual_fitness(i)));
            }
        );
    }
    
    group.finish();
}

/// Benchmark combined evolution operations
fn bench_combined_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("combined_operations");
    let alphabet = Alphabet::dna();
    let mut rng = StdRng::seed_from_u64(42);
    let size = 10_000;
    
    let mutation_model = SubstitutionModel::uniform(alphabet.clone(), 0.001).unwrap();
    let recomb_params = RecombinationParams::new(0.01, 0.7, 0.1).unwrap();
    let gc_fitness = GCContentFitness::new(0.5, 2.0).unwrap();
    
    group.bench_function("mutation_then_fitness", |b| {
        b.iter(|| {
            let mut seq = create_test_sequence(size, alphabet.clone());
            mutation_model.mutate_sequence(&mut seq, &mut rng);
            let chr = Chromosome::new("chr1", seq, 171, 12);
            black_box(gc_fitness.haplotype_fitness(&chr))
        });
    });
    
    group.bench_function("recombination_then_fitness", |b| {
        b.iter(|| {
            let seq1 = create_test_sequence(size, alphabet.clone());
            let seq2 = create_test_sequence(size, alphabet.clone());
            let event = recomb_params.sample_event(size, &mut rng);
            
            let result_seq = match event {
                centrevo::evolution::RecombinationType::Crossover { position } => {
                    let (s1, _) = recomb_params.crossover(&seq1, &seq2, position).unwrap();
                    s1
                }
                _ => seq1,
            };
            
            let chr = Chromosome::new("chr1", result_seq, 171, 12);
            black_box(gc_fitness.haplotype_fitness(&chr))
        });
    });
    
    group.finish();
}

/// Benchmark parallel fitness calculations (simulated)
fn bench_population_fitness(c: &mut Criterion) {
    let mut group = c.benchmark_group("population_fitness");
    let pop_sizes = [10, 100, 1000];
    let gc_fitness = GCContentFitness::new(0.5, 2.0).unwrap();
    
    for pop_size in pop_sizes {
        let individuals: Vec<_> = (0..pop_size)
            .map(|_| create_test_individual(10_000))
            .collect();
        
        group.bench_with_input(
            BenchmarkId::new("sequential", pop_size),
            &individuals,
            |b, inds| {
                b.iter(|| {
                    let fitness: Vec<_> = inds.iter()
                        .map(|ind| gc_fitness.individual_fitness(ind))
                        .collect();
                    black_box(fitness)
                });
            }
        );
    }
    
    group.finish();
}

criterion_group!(
    benches,
    bench_mutation,
    bench_mutate_base,
    bench_recombination,
    bench_fitness,
    bench_combined_operations,
    bench_population_fitness,
);

criterion_main!(benches);
