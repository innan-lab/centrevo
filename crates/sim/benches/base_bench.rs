use centrevo_sim::base::Nucleotide;
use criterion::{Criterion, black_box, criterion_group, criterion_main};

fn bench_nucleotide_conversions(c: &mut Criterion) {
    let mut group = c.benchmark_group("nucleotide_conversions");

    group.bench_function("from_index", |b| {
        b.iter(|| {
            black_box(Nucleotide::from_index(black_box(0)));
            black_box(Nucleotide::from_index(black_box(1)));
            black_box(Nucleotide::from_index(black_box(2)));
            black_box(Nucleotide::from_index(black_box(3)));
        })
    });

    group.bench_function("to_index", |b| {
        b.iter(|| {
            black_box(Nucleotide::A.to_index());
            black_box(Nucleotide::C.to_index());
            black_box(Nucleotide::G.to_index());
            black_box(Nucleotide::T.to_index());
        })
    });

    group.bench_function("from_ascii", |b| {
        b.iter(|| {
            black_box(Nucleotide::from_ascii(black_box(b'A')));
            black_box(Nucleotide::from_ascii(black_box(b'C')));
            black_box(Nucleotide::from_ascii(black_box(b'G')));
            black_box(Nucleotide::from_ascii(black_box(b'T')));
        })
    });

    group.bench_function("complement", |b| {
        b.iter(|| {
            black_box(Nucleotide::A.complement());
            black_box(Nucleotide::C.complement());
            black_box(Nucleotide::G.complement());
            black_box(Nucleotide::T.complement());
        })
    });

    group.finish();
}

criterion_group!(benches, bench_nucleotide_conversions);
criterion_main!(benches);
