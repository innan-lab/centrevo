//! Benchmarks for base module (nucleotide, sequence operations).
use std::str::FromStr;
use std::hint::black_box;
use criterion::{criterion_group, criterion_main, Criterion, BenchmarkId, Throughput};
use centrevo::base::{Nucleotide, Sequence};

/// Benchmark nucleotide conversions
fn bench_nucleotide_conversions(c: &mut Criterion) {
    let mut group = c.benchmark_group("nucleotide_conversions");
    
    group.bench_function("from_ascii", |b| {
        b.iter(|| {
            black_box(Nucleotide::from_ascii(b'A'));
            black_box(Nucleotide::from_ascii(b'C'));
            black_box(Nucleotide::from_ascii(b'G'));
            black_box(Nucleotide::from_ascii(b'T'));
        });
    });
    
    group.bench_function("to_ascii", |b| {
        b.iter(|| {
            black_box(Nucleotide::A.to_ascii());
            black_box(Nucleotide::C.to_ascii());
            black_box(Nucleotide::G.to_ascii());
            black_box(Nucleotide::T.to_ascii());
        });
    });
    
    group.bench_function("complement", |b| {
        b.iter(|| {
            black_box(Nucleotide::A.complement());
            black_box(Nucleotide::C.complement());
            black_box(Nucleotide::G.complement());
            black_box(Nucleotide::T.complement());
        });
    });
    
    group.finish();
}

/// Benchmark sequence creation
fn bench_sequence_creation(c: &mut Criterion) {
    let mut group = c.benchmark_group("sequence_creation");
    
    let sizes = [100, 1_000, 10_000, 100_000];
    
    for size in sizes {
        let seq_str = "ACGT".repeat(size / 4);
        
        group.throughput(Throughput::Elements(size as u64));
        group.bench_with_input(BenchmarkId::new("from_str", size), &size, |b, _| {
            b.iter(|| {
                black_box(Sequence::from_str(&seq_str).unwrap())
            });
        });
        
        group.bench_with_input(BenchmarkId::new("with_capacity", size), &size, |b, &s| {
            b.iter(|| {
                black_box(Sequence::with_capacity(s))
            });
        });
    }
    
    group.finish();
}

/// Benchmark sequence operations
fn bench_sequence_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("sequence_operations");
    
    let sizes = [1_000, 10_000, 100_000];
    
    for size in sizes {
        let seq_str = "ACGT".repeat(size / 4);
        let seq = Sequence::from_str(&seq_str).unwrap();
        
        group.throughput(Throughput::Elements(size as u64));
        
        group.bench_with_input(BenchmarkId::new("clone_mutable", size), &seq, |b, s| {
            b.iter(|| black_box(s.clone()));
        });
        
        group.bench_with_input(BenchmarkId::new("to_shared", size), &seq, |b, s| {
            b.iter(|| black_box(s.to_shared()));
        });
        
        group.bench_with_input(BenchmarkId::new("to_string", size), &seq, |b, s| {
            b.iter(|| black_box(s.to_string()));
        });
        
        // Benchmark shared sequence clone
        let shared = seq.to_shared();
        group.bench_with_input(BenchmarkId::new("clone_shared", size), &shared, |b, s| {
            b.iter(|| black_box(s.clone()));
        });
        
        // Benchmark get operations
        group.bench_with_input(BenchmarkId::new("get_random_access", size), &seq, |b, s| {
            let indices = [0, s.len() / 4, s.len() / 2, 3 * s.len() / 4];
            b.iter(|| {
                for &idx in &indices {
                    black_box(s.get(idx));
                }
            });
        });
    }
    
    group.finish();
}

/// Benchmark sequence iteration
fn bench_sequence_iteration(c: &mut Criterion) {
    let mut group = c.benchmark_group("sequence_iteration");
    
    let sizes = [1_000, 10_000, 100_000];
    
    for size in sizes {
        let seq_str = "ACGT".repeat(size / 4);
        let seq = Sequence::from_str(&seq_str).unwrap();
        
        group.throughput(Throughput::Elements(size as u64));
        
        group.bench_with_input(BenchmarkId::new("iterate_indices", size), &seq, |b, s| {
            b.iter(|| {
                let mut sum = 0u64;
                for &idx in s.as_slice() {
                    sum += idx as u64;
                }
                black_box(sum)
            });
        });
        
        group.bench_with_input(BenchmarkId::new("iterate_get", size), &seq, |b, s| {
            b.iter(|| {
                let mut count = 0;
                for i in 0..s.len() {
                    if s.get(i).is_some() {
                        count += 1;
                    }
                }
                black_box(count)
            });
        });
    }
    
    group.finish();
}

/// Benchmark sequence modifications
fn bench_sequence_modifications(c: &mut Criterion) {
    let mut group = c.benchmark_group("sequence_modifications");
    
    
    group.bench_function("push_100", |b| {
        b.iter(|| {
            let mut seq = Sequence::with_capacity(100);
            for _ in 0..100 {
                seq.push(Nucleotide::A);
            }
            black_box(seq)
        });
    });
    
    group.bench_function("set_sequential", |b| {
        b.iter(|| {
            let mut seq = Sequence::from_str(&"A".repeat(1000)).unwrap();
            for i in 0..seq.len() {
                let _ = seq.set(i, Nucleotide::C);
            }
            black_box(seq)
        });
    });
    
    group.bench_function("insert_middle", |b| {
        b.iter(|| {
            let mut seq = Sequence::from_str(&"ACGT".repeat(100)).unwrap();
            seq.insert(seq.len() / 2, Nucleotide::G);
            black_box(seq)
        });
    });
    
    group.finish();
}

/// Benchmark memory usage patterns
fn bench_memory_patterns(c: &mut Criterion) {
    let mut group = c.benchmark_group("memory_patterns");
    
    let size = 10_000;
    let seq_str = "ACGT".repeat(size / 4);
    
    group.bench_function("create_1000_mutable", |b| {
        b.iter(|| {
            let mut sequences = Vec::with_capacity(1000);
            for _ in 0..1000 {
                sequences.push(Sequence::from_str(&seq_str).unwrap());
            }
            black_box(sequences)
        });
    });
    
    group.bench_function("create_1000_shared", |b| {
        let base_seq = Sequence::from_str(&seq_str).unwrap().to_shared();
        b.iter(|| {
            let mut sequences = Vec::with_capacity(1000);
            for _ in 0..1000 {
                sequences.push(base_seq.clone());
            }
            black_box(sequences)
        });
    });
    
    group.finish();
}

criterion_group!(
    benches,
    bench_nucleotide_conversions,
    bench_sequence_creation,
    bench_sequence_operations,
    bench_sequence_iteration,
    bench_sequence_modifications,
    bench_memory_patterns,
);

criterion_main!(benches);
