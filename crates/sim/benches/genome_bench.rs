use centrevo_sim::base::Nucleotide;
use centrevo_sim::genome::Chromosome;
use criterion::{BenchmarkId, Criterion, black_box, criterion_group, criterion_main};

fn bench_chromosome_ops(c: &mut Criterion) {
    let mut group = c.benchmark_group("chromosome_ops");
    let mut arena = centrevo_sim::base::GenomeArena::new();

    group.bench_function("uniform_creation", |b| {
        b.iter(|| {
            black_box(Chromosome::uniform(
                black_box("chr1"),
                black_box(Nucleotide::A),
                black_box(171),
                black_box(12),
                black_box(100),
                &mut arena,
            ));
        })
    });

    // Similarity benchmarking
    let ru_length = 171;
    let hor_lengths = [1, 10, 100]; // RUs per HOR
    let kmer_sizes = [5, 7, 9];

    for &hors in &hor_lengths {
        // Create chromosomes with different HOR lengths
        // Total length = 171 * hors * 1 (1 HOR per chromosome for simplicity of this bench)
        let chr1 = Chromosome::uniform("chr1", Nucleotide::A, ru_length, hors, 1, &mut arena);
        let chr2 = Chromosome::uniform("chr2", Nucleotide::T, ru_length, hors, 1, &mut arena);

        for &k in &kmer_sizes {
            let parameter_string = format!("hors={hors}/k={k}");

            group.bench_with_input(
                BenchmarkId::new("calculate_similarity", &parameter_string),
                &(hors, k),
                |b, &(_, k)| {
                    b.iter(|| {
                        black_box(chr1.calculate_similarity(
                            black_box(0),
                            black_box(&chr2),
                            black_box(0),
                            black_box(k),
                            &arena,
                        ));
                    })
                },
            );
        }
    }

    // Other operations with fixed size for regression testing
    let chr1 = Chromosome::uniform("chr1", Nucleotide::A, 171, 12, 10, &mut arena);
    let chr2 = Chromosome::uniform("chr2", Nucleotide::T, 171, 12, 10, &mut arena);

    group.bench_function("gc_content", |b| {
        b.iter(|| {
            black_box(chr1.gc_content(&arena));
        })
    });

    group.bench_function("crossover", |b| {
        b.iter(|| {
            black_box(
                chr1.crossover(
                    black_box(&chr2),
                    black_box(5000),
                    black_box(5000),
                    &mut arena,
                )
                .unwrap(),
            );
        })
    });

    group.bench_function("gene_conversion", |b| {
        b.iter(|| {
            black_box(
                chr1.gene_conversion(
                    black_box(&chr2),
                    black_box(1000),
                    black_box(1000),
                    black_box(500),
                    &mut arena,
                )
                .unwrap(),
            );
        })
    });

    group.finish();
}

criterion_group!(benches, bench_chromosome_ops);
criterion_main!(benches);
