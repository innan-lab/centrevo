# Centrevo Benchmarks

This directory contains comprehensive performance benchmarks for all Centrevo modules using the [Criterion](https://github.com/bheisler/criterion.rs) benchmarking framework.

## Benchmark Suites

### 1. Base Module (`base_benchmarks.rs`)
Benchmarks for fundamental data structures:
- **Nucleotide Conversions**: ASCII conversion, complement operations
- **Alphabet Operations**: Creation, cloning, character/index lookups
- **Sequence Creation**: From string, with capacity
- **Sequence Operations**: Cloning (mutable vs shared), conversion to string
- **Sequence Iteration**: Index iteration, get() method access
- **Sequence Modifications**: Push, set, insert operations
- **Memory Patterns**: Comparing mutable vs shared sequence allocation

### 2. Genome Module (`genome_benchmarks.rs`)
Benchmarks for genome structures:
- **Chromosome Creation**: Different sizes (1K, 10K, 100K bases)
- **Chromosome Operations**: Clone, GC content calculation, string formatting
- **Haplotype Operations**: Multiple chromosomes per haplotype
- **Individual Operations**: Diploid organism operations
- **Population Memory**: Creating and cloning large populations

### 3. Evolution Module (`evolution_benchmarks.rs`)
Benchmarks for evolutionary processes:
- **Mutation Operations**: Different mutation rates (0.001, 0.01, 0.1)
- **Single Base Mutation**: Individual nucleotide mutation
- **Recombination**: Event sampling, crossover, gene conversion
- **Fitness Calculations**: GC content, length-based, sequence similarity
- **Combined Operations**: Mutation+fitness, recombination+fitness
- **Population Fitness**: Sequential fitness calculation for populations

### 4. Simulation Module (`simulation_benchmarks.rs`)
Benchmarks for complete simulations:
- **Simulation Initialization**: Different population sizes and chromosome lengths
- **Single Generation**: With neutral and selective fitness
- **Multi-Generation Runs**: 5, 10, 20 generations
- **Population Operations**: Fitness computation, parent selection
- **Memory Allocation**: Different population/chromosome size combinations
- **Mutation Rate Comparison**: Testing different mutation rates
- **Scaling Tests**: Population size and chromosome length scaling

## Running Benchmarks

### Run All Benchmarks
```bash
cargo bench
```

### Run Specific Benchmark Suite
```bash
cargo bench --bench base_benchmarks
cargo bench --bench genome_benchmarks
cargo bench --bench evolution_benchmarks
cargo bench --bench simulation_benchmarks
```

### Run Specific Benchmark Group
```bash
# Run only nucleotide conversion benchmarks
cargo bench --bench base_benchmarks nucleotide_conversions

# Run only mutation benchmarks
cargo bench --bench evolution_benchmarks mutation

# Run only scaling benchmarks
cargo bench --bench simulation_benchmarks scaling
```

### Save Baseline for Comparison
```bash
# Save current performance as baseline
cargo bench -- --save-baseline my-baseline

# Compare against baseline
cargo bench -- --baseline my-baseline
```

## Interpreting Results

Criterion provides:
- **Mean time**: Average execution time
- **Standard deviation**: Variability in measurements
- **Throughput**: Operations per second (when applicable)
- **Outliers**: Measurements that deviate significantly
- **Change detection**: Statistical comparison with previous runs

### Example Output
```
sequence_creation/from_str/1000
                        time:   [1.2345 µs 1.2456 µs 1.2567 µs]
                        thrpt: [795.85 Kelem/s 802.89 Kelem/s 809.93 Kelem/s]
Found 3 outliers among 100 measurements (3.00%)
  2 (2.00%) high mild
  1 (1.00%) high severe
```

## HTML Reports

Criterion generates detailed HTML reports in `target/criterion/`:
- Statistical analysis
- Performance plots
- Comparison charts
- Regression analysis

View reports:
```bash
open target/criterion/report/index.html
```

## Performance Tips

### What to Look For

1. **Memory Allocation**:
   - Shared sequences should be ~100x faster to clone than mutable sequences
   - Arc-based alphabet sharing should show minimal overhead

2. **Scaling**:
   - Operations should scale linearly with sequence length (O(n))
   - Population operations should scale linearly with population size

3. **Hot Paths**:
   - Mutation: Should process 100K bases in <100µs
   - Fitness: GC content calculation should be <10µs for 10K bases
   - Recombination: Crossover should be <50µs for 10K bases

### Optimization Targets

Based on benchmark results, focus optimization on:
- Operations that don't scale linearly
- High-frequency operations (mutation, fitness calculation)
- Memory-intensive operations (population cloning)

## Continuous Integration

These benchmarks can be integrated into CI/CD:
```bash
# Quick smoke test (reduced sample size)
cargo bench -- --test

# Full benchmark with baseline comparison
cargo bench -- --save-baseline ci-baseline
```

## Adding New Benchmarks

When adding new functionality:

1. Add benchmarks to appropriate suite
2. Include multiple data sizes
3. Test edge cases (empty, very large)
4. Compare against baseline
5. Document expected performance

Example:
```rust
fn bench_my_operation(c: &mut Criterion) {
    let mut group = c.benchmark_group("my_operation");
    let sizes = [100, 1_000, 10_000];
    
    for size in sizes {
        group.throughput(Throughput::Elements(size as u64));
        group.bench_with_input(
            BenchmarkId::new("operation", size),
            &size,
            |b, &s| {
                b.iter(|| black_box(my_operation(s)));
            }
        );
    }
    
    group.finish();
}
```

## Performance Regression Detection

Criterion automatically detects performance regressions:
- **<5% change**: Likely noise
- **5-10% change**: Worth investigating
- **>10% change**: Significant regression/improvement

Use baselines to track performance over time:
```bash
# Before changes
cargo bench -- --save-baseline before

# After changes
cargo bench -- --baseline before
```

## Memory Profiling

For detailed memory analysis, use additional tools:

```bash
# Valgrind (Linux)
valgrind --tool=massif cargo bench --bench simulation_benchmarks -- --test

# Instruments (macOS)
cargo instruments --bench simulation_benchmarks --template Allocations

# Heaptrack (Linux)
heaptrack cargo bench --bench simulation_benchmarks -- --test
```

## Troubleshooting

### Benchmarks Taking Too Long
```bash
# Reduce sample size
cargo bench -- --sample-size 10

# Quick test mode
cargo bench -- --test
```

### Inconsistent Results
- Close other applications
- Disable CPU frequency scaling
- Run multiple times and average
- Use `--save-baseline` for comparison

### Out of Memory
- Reduce test data sizes
- Run benchmarks individually
- Increase system swap space

## Resources

- [Criterion.rs Documentation](https://bheisler.github.io/criterion.rs/book/)
- [Rust Performance Book](https://nnethercote.github.io/perf-book/)
- [Centrevo Implementation Guide](../IMPLEMENTATION_GUIDE.md)
