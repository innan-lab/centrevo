use centrevo_sim::base::{FitnessValue, Nucleotide};
use centrevo_sim::genome::{Chromosome, Haplotype, Individual};
use centrevo_sim::simulation::Population;
use centrevo_sim::storage::{Recorder, RecordingStrategy};
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

    group.throughput(Throughput::Elements(pop_size as u64));

    group.bench_function("record_generation", |b| {
        b.iter_batched(
            || {
                let file = NamedTempFile::new().unwrap();
                let path = file.path().to_owned();
                // Keep file alive
                let recorder = Recorder::new(&path, "bench_sim", RecordingStrategy::All).unwrap();
                let pop = create_test_population(pop_size, chr_length);
                (file, recorder, pop)
            },
            |(_file, mut recorder, pop)| {
                recorder
                    .record_generation(black_box(&pop), black_box(0))
                    .unwrap();
            },
            criterion::BatchSize::SmallInput,
        )
    });

    group.finish();
}

criterion_group!(benches, bench_recorder_write);
criterion_main!(benches);
