//! Benchmarks for simulation module (population, engine, full simulation runs).

use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId, Throughput};
use centrevo::base::{Alphabet, Nucleotide};
use centrevo::simulation::{
    Simulation, RepeatStructure, MutationConfig, RecombinationConfig,
    FitnessConfig, SimulationConfig,
};
use centrevo::evolution::{GCContentFitness};

fn create_test_simulation(
    pop_size: usize,
    chr_length: usize,
    generations: usize,
) -> Simulation {
    let alphabet = Alphabet::dna();
    
    let structure = RepeatStructure::new(
        alphabet.clone(),
        Nucleotide::A,
        171,
        12,
        chr_length / (171 * 12),  // Calculate HORs to get desired length
        1,
    );
    
    let mutation = MutationConfig::uniform(alphabet, 0.001).unwrap();
    let recombination = RecombinationConfig::standard(0.01, 0.7, 0.1).unwrap();
    let fitness = FitnessConfig::neutral();
    let config = SimulationConfig::new(pop_size, generations, Some(42));
    
    Simulation::new(structure, mutation, recombination, fitness, config).unwrap()
}

fn create_simulation_with_fitness(
    pop_size: usize,
    chr_length: usize,
    generations: usize,
) -> Simulation {
    let alphabet = Alphabet::dna();
    
    let structure = RepeatStructure::new(
        alphabet.clone(),
        Nucleotide::A,
        171,
        12,
        chr_length / (171 * 12),
        1,
    );
    
    let mutation = MutationConfig::uniform(alphabet, 0.001).unwrap();
    let recombination = RecombinationConfig::standard(0.01, 0.7, 0.1).unwrap();
    
    let gc_fitness = GCContentFitness::new(0.5, 2.0).unwrap();
    let fitness = FitnessConfig::new(Some(gc_fitness), None, None);
    
    let config = SimulationConfig::new(pop_size, generations, Some(42));
    
    Simulation::new(structure, mutation, recombination, fitness, config).unwrap()
}

/// Benchmark simulation initialization
fn bench_simulation_init(c: &mut Criterion) {
    let mut group = c.benchmark_group("simulation_init");
    let pop_sizes = [10, 50, 100];
    let chr_lengths = [10_000, 50_000, 100_000];
    
    for pop_size in pop_sizes {
        for chr_length in chr_lengths {
            let label = format!("pop{}_len{}", pop_size, chr_length);
            
            group.bench_with_input(
                BenchmarkId::new("create_neutral", &label),
                &(pop_size, chr_length),
                |b, &(p, l)| {
                    b.iter(|| black_box(create_test_simulation(p, l, 10)));
                }
            );
            
            group.bench_with_input(
                BenchmarkId::new("create_with_fitness", &label),
                &(pop_size, chr_length),
                |b, &(p, l)| {
                    b.iter(|| black_box(create_simulation_with_fitness(p, l, 10)));
                }
            );
        }
    }
    
    group.finish();
}

/// Benchmark single generation step
fn bench_single_generation(c: &mut Criterion) {
    let mut group = c.benchmark_group("single_generation");
    let pop_sizes = [10, 50, 100];
    let chr_length = 10_000;
    
    for pop_size in pop_sizes {
        group.throughput(Throughput::Elements(pop_size as u64));
        
        // Neutral fitness
        group.bench_with_input(
            BenchmarkId::new("neutral", pop_size),
            &pop_size,
            |b, &p| {
                b.iter_batched(
                    || create_test_simulation(p, chr_length, 1),
                    |mut sim| {
                        sim.step().unwrap();
                        black_box(sim)
                    },
                    criterion::BatchSize::SmallInput,
                );
            }
        );
        
        // With fitness selection
        group.bench_with_input(
            BenchmarkId::new("with_fitness", pop_size),
            &pop_size,
            |b, &p| {
                b.iter_batched(
                    || create_simulation_with_fitness(p, chr_length, 1),
                    |mut sim| {
                        sim.step().unwrap();
                        black_box(sim)
                    },
                    criterion::BatchSize::SmallInput,
                );
            }
        );
    }
    
    group.finish();
}

/// Benchmark multi-generation runs
fn bench_multi_generation(c: &mut Criterion) {
    let mut group = c.benchmark_group("multi_generation");
    group.sample_size(10);  // Reduce sample size for longer benchmarks
    
    let pop_size = 50;
    let chr_length = 10_000;
    let generation_counts = [5, 10, 20];
    
    for n_gens in generation_counts {
        group.throughput(Throughput::Elements((pop_size * n_gens) as u64));
        
        group.bench_with_input(
            BenchmarkId::new("run", n_gens),
            &n_gens,
            |b, &n| {
                b.iter_batched(
                    || create_test_simulation(pop_size, chr_length, n),
                    |mut sim| {
                        sim.run().unwrap();
                        black_box(sim)
                    },
                    criterion::BatchSize::SmallInput,
                );
            }
        );
    }
    
    group.finish();
}

/// Benchmark population operations
fn bench_population_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("population_operations");
    let pop_sizes = [10, 50, 100, 500];
    let chr_length = 10_000;
    
    for pop_size in pop_sizes {
        let mut sim = create_test_simulation(pop_size, chr_length, 1);
        
        group.bench_with_input(
            BenchmarkId::new("compute_fitness", pop_size),
            &sim,
            |b, s| {
                b.iter(|| {
                    let fitness_config = FitnessConfig::neutral();
                    black_box(s.population().compute_fitness(&fitness_config))
                });
            }
        );
        
        // Test with GC fitness
        let gc_fitness = GCContentFitness::new(0.5, 2.0).unwrap();
        let fitness_config = FitnessConfig::new(Some(gc_fitness), None, None);
        
        group.bench_with_input(
            BenchmarkId::new("compute_fitness_gc", pop_size),
            &sim,
            |b, s| {
                b.iter(|| {
                    black_box(s.population().compute_fitness(&fitness_config))
                });
            }
        );
    }
    
    group.finish();
}

/// Benchmark memory allocation patterns
fn bench_memory_allocation(c: &mut Criterion) {
    let mut group = c.benchmark_group("memory_allocation");
    
    group.bench_function("init_pop100_chr10k", |b| {
        b.iter(|| {
            black_box(create_test_simulation(100, 10_000, 1))
        });
    });
    
    group.bench_function("init_pop10_chr100k", |b| {
        b.iter(|| {
            black_box(create_test_simulation(10, 100_000, 1))
        });
    });
    
    group.bench_function("run_5gen_pop50", |b| {
        b.iter_batched(
            || create_test_simulation(50, 10_000, 5),
            |mut sim| {
                sim.run().unwrap();
                black_box(sim)
            },
            criterion::BatchSize::SmallInput,
        );
    });
    
    group.finish();
}

/// Benchmark different mutation rates
fn bench_mutation_rates(c: &mut Criterion) {
    let mut group = c.benchmark_group("mutation_rates");
    let pop_size = 50;
    let chr_length = 10_000;
    let rates = [0.0001, 0.001, 0.01];
    
    for rate in rates {
        let alphabet = Alphabet::dna();
        let structure = RepeatStructure::new(
            alphabet.clone(),
            Nucleotide::A,
            171,
            12,
            chr_length / (171 * 12),
            1,
        );
        
        let mutation = MutationConfig::uniform(alphabet, rate).unwrap();
        let recombination = RecombinationConfig::standard(0.01, 0.7, 0.1).unwrap();
        let fitness = FitnessConfig::neutral();
        let config = SimulationConfig::new(pop_size, 1, Some(42));
        
        group.bench_with_input(
            BenchmarkId::new("rate", rate),
            &rate,
            |b, _| {
                b.iter_batched(
                    || Simulation::new(
                        structure.clone(),
                        mutation.clone(),
                        recombination.clone(),
                        fitness.clone(),
                        config.clone(),
                    ).unwrap(),
                    |mut sim| {
                        sim.step().unwrap();
                        black_box(sim)
                    },
                    criterion::BatchSize::SmallInput,
                );
            }
        );
    }
    
    group.finish();
}

/// Benchmark scaling characteristics
fn bench_scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("scaling");
    group.sample_size(10);
    
    // Test population size scaling
    let chr_length = 10_000;
    let pop_sizes = [10, 50, 100, 200];
    
    for pop_size in pop_sizes {
        group.throughput(Throughput::Elements(pop_size as u64));
        group.bench_with_input(
            BenchmarkId::new("pop_scaling", pop_size),
            &pop_size,
            |b, &p| {
                b.iter_batched(
                    || create_test_simulation(p, chr_length, 1),
                    |mut sim| {
                        sim.step().unwrap();
                        black_box(sim)
                    },
                    criterion::BatchSize::SmallInput,
                );
            }
        );
    }
    
    // Test chromosome length scaling
    let pop_size = 50;
    let chr_lengths = [5_000, 10_000, 20_000, 50_000];
    
    for chr_length in chr_lengths {
        group.throughput(Throughput::Elements((chr_length * pop_size) as u64));
        group.bench_with_input(
            BenchmarkId::new("length_scaling", chr_length),
            &chr_length,
            |b, &l| {
                b.iter_batched(
                    || create_test_simulation(pop_size, l, 1),
                    |mut sim| {
                        sim.step().unwrap();
                        black_box(sim)
                    },
                    criterion::BatchSize::SmallInput,
                );
            }
        );
    }
    
    group.finish();
}

criterion_group!(
    benches,
    bench_simulation_init,
    bench_single_generation,
    bench_multi_generation,
    bench_population_operations,
    bench_memory_allocation,
    bench_mutation_rates,
    bench_scaling,
);

criterion_main!(benches);
