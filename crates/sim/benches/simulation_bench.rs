use centrevo_sim::simulation::SimulationBuilder;
use criterion::{BenchmarkId, Criterion, Throughput, black_box, criterion_group, criterion_main};

fn bench_simulation_init(c: &mut Criterion) {
    let mut group = c.benchmark_group("simulation_init");

    group.bench_function("default_init", |b| {
        b.iter(|| {
            black_box(
                SimulationBuilder::new()
                    .population_size(black_box(50))
                    .generations(black_box(10))
                    .repeat_structure(black_box(171), black_box(12), black_box(10))
                    .build()
                    .unwrap(),
            );
        })
    });

    group.finish();
}

fn bench_simulation_step(c: &mut Criterion) {
    let mut group = c.benchmark_group("simulation_step");
    let pop_size = 50;

    group.throughput(Throughput::Elements(pop_size as u64));

    group.bench_function("step_neutral", |b| {
        b.iter_batched(
            || {
                SimulationBuilder::new()
                    .population_size(pop_size)
                    .generations(10)
                    .repeat_structure(171, 12, 10)
                    .mutation_rate(0.001)
                    .recombination(0.01, 0.7, 0.1)
                    .build()
                    .unwrap()
            },
            |mut sim| {
                sim.step().unwrap();
                black_box(sim)
            },
            criterion::BatchSize::SmallInput,
        )
    });

    group.finish();
}

fn bench_simulation_run(c: &mut Criterion) {
    let mut group = c.benchmark_group("simulation_run");
    let pop_size = 50;
    let generations = 5;

    group.throughput(Throughput::Elements((pop_size * generations) as u64));

    group.bench_with_input(
        BenchmarkId::new("run_full", generations),
        &generations,
        |b, &gens| {
            b.iter_batched(
                || {
                    SimulationBuilder::new()
                        .population_size(pop_size)
                        .generations(gens)
                        .repeat_structure(171, 12, 10)
                        .mutation_rate(0.001)
                        .recombination(0.01, 0.7, 0.1)
                        .build()
                        .unwrap()
                },
                |mut sim| {
                    sim.run().unwrap();
                    black_box(sim)
                },
                criterion::BatchSize::SmallInput,
            )
        },
    );

    group.finish();
}

criterion_group!(
    benches,
    bench_simulation_init,
    bench_simulation_step,
    bench_simulation_run
);
criterion_main!(benches);
