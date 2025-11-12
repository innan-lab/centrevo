# Centrevo Implementation Guide

This document provides in-depth documentation of the Centrevo Rust implementation. For **usage documentation**, see [README.md](../README.md), [CLI.md](../CLI.md), or [PYTHON.md](../PYTHON.md) in the root directory.

## Architecture Overview

Centrevo is organized into modular crates within a single workspace:

```
src/
├── lib.rs              # Library entry point
├── prelude.rs          # Common imports
├── base/               # Core data structures
├── genome/             # Genetic structures
├── evolution/          # Evolutionary processes
├── simulation/         # Simulation engine
├── storage/            # Database persistence
├── analysis/           # Population genetics analysis
├── python/             # Python bindings (optional)
├── utils/              # Utility functions
└── bin/
    └── centrevo.rs     # CLI entry point
```

## Module Documentation

### base/ - Core Data Structures

**Purpose**: Fundamental building blocks for sequences and alphabets.

#### `nucleotide.rs`
- `Nucleotide` enum: A, C, G, T representation
- Conversion to/from ASCII characters
- Complement operations
- Efficient encoding (single byte per nucleotide)

#### `alphabet.rs`
- `Alphabet` struct: Collection of valid nucleotides
- DNA/RNA alphabet definitions
- Arc-based sharing for memory efficiency
- Character/index lookup tables

#### `sequence.rs`
- `Sequence` struct: Vector of nucleotides with Arc-based sharing
- Efficient cloning through reference counting
- Mutation and recombination support
- String conversion utilities

**Key Design Decisions**:
- Arc wrapping for cheap cloning during reproduction
- Mutable operations return new sequences (copy-on-write semantics)
- Optimized for repeated sequence operations

### genome/ - Genetic Structures

**Purpose**: Higher-level genetic organization (chromosomes, individuals, populations).

#### `chromosome.rs`
- `Chromosome` struct: Single chromosome with metadata
- Repeat unit (RU) tracking
- Higher-order repeat (HOR) organization
- GC content calculation
- Uniform initialization from single base

#### `haplotype.rs`
- `Haplotype` struct: Collection of chromosomes (1n)
- Supports multiple chromosomes per haplotype
- Efficient cloning through Arc sharing

#### `individual.rs`
- `Individual` struct: Diploid organism (2n)
- Two haplotypes per individual
- Unique ID tracking
- Fitness caching

#### `population.rs`
- `Population` struct: Collection of individuals
- Generation tracking
- Population-level statistics
- Parallel operations via Rayon

**Key Design Decisions**:
- Separation of ploidy levels (haplotype vs diploid)
- ID-based tracking for lineage analysis
- Arc sharing at sequence level reduces memory

### evolution/ - Evolutionary Processes

**Purpose**: Mutation, recombination, and selection mechanisms.

#### `mutation.rs`
- Point mutation implementation
- Configurable per-site mutation rates
- Uniform and position-specific mutation models
- Parallel mutation across sequences
- RNG state management for reproducibility

**Mutation Models**:
- `Uniform`: Constant rate across all sites
- `PositionSpecific`: Different rates per site (future)

#### `recombination.rs`
- Crossover and gene conversion
- Break point selection
- Tract length distributions
- Non-crossover gene conversion support

**Recombination Models**:
- Uniform break probability
- Crossover vs. non-crossover resolution
- Configurable tract lengths

#### `selection.rs`
- Fitness calculation framework
- Neutral, GC-based, and custom fitness functions
- Efficient Wright-Fisher sampling with fitness
- Deterministic and stochastic selection

**Fitness Models**:
- `Neutral`: All individuals equal fitness
- `GCContent`: Fitness based on GC content
- `SequenceSimilarity`: Fitness based on similarity to target
- `Custom`: User-defined fitness functions (Python API)

**Key Design Decisions**:
- Trait-based fitness for extensibility
- RNG passed explicitly for reproducibility
- Parallel fitness calculation for large populations

### simulation/ - Simulation Engine

**Purpose**: Coordinate evolutionary processes and manage simulation state.

#### `engine.rs`
- `Simulation` struct: Main simulation driver
- Generation stepping
- Population evolution orchestration
- Recording integration
- Checkpoint/resume support

**Workflow**:
1. Initialize population
2. Calculate fitness
3. Select parents (with fitness)
4. Apply recombination
5. Apply mutation
6. Update generation
7. Record if needed
8. Repeat

#### `initialization.rs`
- Uniform population creation
- Custom sequence initialization
  - FASTA file parsing
  - JSON parsing
  - Database loading
- Validation and error reporting

**Custom Initialization Features**:
- Flexible sequence ID formats
- Length validation against structure
- Detailed error messages with line numbers
- Support for resuming from previous simulations

#### `parameters.rs`
- `SimulationConfig`: Central configuration
- Parameter validation
- Serialization for storage
- Default values

#### `population.rs`
- Population management utilities
- Parent selection algorithms
- Fitness statistics

**Key Design Decisions**:
- Separation of initialization from evolution
- Explicit RNG state for reproducibility
- Checkpoint includes full RNG state
- Validation at initialization time

### storage/ - Database Persistence

**Purpose**: SQLite-based persistence for simulation state and results.

See [storage/README.md](storage/README.md) for complete documentation.

#### Key Features:
- SQLite with WAL mode for concurrent access
- Multiple recording strategies (All, EveryN, Specific, None)
- Efficient BLOB storage for sequences
- Indexed queries for fast retrieval
- Checkpoint storage for resume capability

#### Schema:
- `simulations`: Metadata and configuration
- `population_state`: Individual snapshots
- `fitness_history`: Aggregated fitness statistics
- `checkpoints`: RNG state and generation for resume

### analysis/ - Population Genetics Analysis

**Purpose**: Calculate population genetics metrics and statistics.

#### `diversity.rs`
- Nucleotide diversity (π)
- Watterson's theta (θ_W)
- Tajima's D
- Haplotype diversity
- Parallel computation for large populations

#### `linkage.rs`
- Pairwise LD (D, D', r²)
- LD decay curves
- Haplotype block identification
- Efficient allele frequency calculation

#### `distance.rs`
- Hamming distance between sequences
- Pairwise distance matrices
- Parallel distance computation
- Normalized distances (0.0 to 1.0)

#### `composition.rs`
- GC content calculation
- Nucleotide composition
- Flexible aggregation (per-chromosome, per-individual, per-population)

#### `polymorphism.rs`
- Segregating site counting
- Site frequency spectrum (SFS)
- Polymorphism statistics

#### `temporal.rs`
- Time-series analysis utilities
- Trajectory extraction
- Change detection

#### `structure.rs`
- Population structure analysis
- Repeat unit organization
- HOR identification

#### `utils.rs`
- Common analysis utilities
- Statistical functions
- Binning and windowing

**Key Design Decisions**:
- Parallel computation via Rayon where beneficial
- Flexible aggregation levels (chromosome/individual/population)
- Efficient allele frequency caching
- Iterator-based APIs for memory efficiency

### python/ - Python Bindings

**Purpose**: Expose Rust functionality to Python via PyO3.

See [python/README.md](../python/README.md) for detailed implementation documentation.

#### `bindings.rs`
- Core type wrappers (Nucleotide, Sequence, Chromosome, etc.)
- Python-friendly constructors
- Automatic type conversion
- Error handling and exceptions

#### `analysis.rs`
- Analysis function bindings
- PyArrow export functions
- Result conversion to Python types

#### `mod.rs`
- Module registration
- Function exports
- Module initialization

**Key Design Decisions**:
- Zero-copy data exchange via PyArrow
- Minimal Python-side computation
- Type stubs for IDE support
- Error messages include context

### utils/ - Utility Functions

**Purpose**: Cross-cutting utilities and helpers.

- RNG utilities and seeding
- Math and statistics helpers
- String formatting
- Validation functions
- Testing utilities

## Design Principles

### Memory Management

1. **Arc-based Sharing**: Sequences use `Arc` for cheap cloning
   - Cloning a population is O(n) not O(n×m) where m is sequence length
   - Enables efficient parent copying during reproduction

2. **Copy-on-Write**: Mutations create new sequences
   - Original sequences remain immutable
   - Safe for concurrent access

3. **No Unsafe Code**: Pure safe Rust
   - Memory safety guaranteed by compiler
   - No manual allocation/deallocation

### Parallelism

1. **Rayon for Data Parallelism**:
   - Population-level operations parallelized
   - Fitness calculation
   - Analysis functions
   - Distance matrices

2. **Thread Safety**:
   - Sequences are Arc-based (thread-safe)
   - RNG state not shared (each thread has own)
   - No global mutable state

3. **When Not to Parallelize**:
   - Small populations (<100 individuals)
   - Overhead exceeds benefit
   - Sequential dependencies (e.g., generation stepping)

### Reproducibility

1. **Explicit RNG State**:
   - ChaCha8Rng for cryptographic-quality randomness
   - Seeded explicitly
   - Passed to all random operations
   - Checkpoint includes RNG state

2. **Deterministic Operations**:
   - No system randomness
   - No timing-based seeds (unless requested)
   - Parallel operations produce same results

3. **Checkpoint/Resume**:
   - Full state serialization
   - Bit-for-bit reproducibility
   - Validated on resume

### Error Handling

1. **Result Types**:
   - Functions return `Result<T, E>`
   - Errors propagate with `?` operator
   - Rich error context

2. **Custom Error Types**:
   - `DatabaseError` for storage issues
   - `SimulationError` for runtime issues
   - `InitializationError` for setup problems

3. **Validation**:
   - Early validation in constructors
   - Panic-free APIs (except for bugs)
   - Descriptive error messages

## Performance Optimization

### Hot Paths

The most performance-critical code paths:

1. **Mutation** (`evolution/mutation.rs`)
   - Per-site random number generation
   - Nucleotide selection
   - Sequence copying

2. **Fitness Calculation** (`evolution/selection.rs`)
   - GC content computation
   - Sequence comparison
   - Statistical calculations

3. **Parent Selection** (`simulation/population.rs`)
   - Wright-Fisher sampling
   - Fitness-weighted selection
   - RNG for random choice

4. **Analysis Functions** (`analysis/`)
   - Pairwise comparisons (O(n²))
   - Allele frequency calculation
   - LD computation

### Optimization Strategies

1. **Avoid Allocations**:
   - Reuse buffers where possible
   - Iterator-based APIs
   - In-place operations when safe

2. **Parallelism**:
   - Rayon for embarassingly parallel work
   - Parallel iterators (`par_iter()`)
   - Work stealing for load balancing

3. **Caching**:
   - Fitness values cached in Individual
   - Allele frequencies cached during analysis
   - Sequence lengths cached

4. **Algorithm Selection**:
   - Linear-time algorithms preferred
   - Avoid unnecessary O(n²) operations
   - Approximations for large datasets (future)

### Benchmarking

See [benches/README.md](../benches/README.md) for performance benchmarks.

Key benchmarks:
- Base operations (sequence creation, cloning)
- Genome operations (chromosome, individual, population)
- Evolution operations (mutation, recombination, fitness)
- Simulation (per-generation time, scaling)
- Analysis (diversity, LD, distance)
- Storage (recording, querying)

## Testing Strategy

### Unit Tests

Each module has comprehensive unit tests:
- Functions tested in isolation
- Edge cases covered
- Property-based testing where appropriate

```bash
# Run all tests
cargo test

# Run specific module
cargo test --lib genome

# Run with coverage
cargo tarpaulin
```

### Integration Tests

Located in `tests/`:
- End-to-end workflows
- Cross-module interactions
- CLI integration
- Database operations

```bash
# Run integration tests
cargo test --test simulation_workflow
cargo test --test storage_query
```

### Property-Based Testing

Using `proptest` for:
- Sequence operations maintain invariants
- Mutation preserves sequence length
- Recombination produces valid offspring

### Reproducibility Tests

Special tests verify:
- Same seed produces same results
- Checkpoint/resume is exact
- Parallel operations are deterministic

## Development Workflow

### Adding a New Feature

1. **Design**: Write design doc or RFC (for major features)
2. **Implementation**: Write code following style guide
3. **Tests**: Add unit and integration tests
4. **Benchmarks**: Add benchmarks if performance-critical
5. **Documentation**: Update relevant docs
6. **PR**: Submit pull request with clear description

### Code Style

- Follow Rust conventions (`rustfmt`)
- Use `clippy` for linting
- Write doc comments for public APIs
- Include examples in doc comments

```bash
# Format code
cargo fmt

# Run linter
cargo clippy

# Generate docs
cargo doc --open
```

### Debugging

```bash
# Debug build
cargo build

# Run with debug output
RUST_LOG=debug cargo run

# Use debugger (lldb on macOS)
rust-lldb target/debug/centrevo
```

## Future Enhancements

### Planned Features

1. **Performance**:
   - SIMD for sequence operations
   - GPU acceleration for large populations
   - Streaming analysis for huge datasets

2. **Models**:
   - Variable recombination rates
   - Complex fitness landscapes
   - Gene flow and migration
   - Non-equilibrium demography

3. **Analysis**:
   - Coalescent trees
   - Ancestral sequence reconstruction
   - Selection tests (MK test, HKA test)
   - ABC inference

4. **Storage**:
   - Compression for large sequences
   - Hierarchical storage (hot/cold)
   - Cloud storage backends

5. **API**:
   - Streaming API for real-time analysis
   - Async simulation for long runs
   - Web API for remote execution

### Contributing

See [CONTRIBUTING.md](../CONTRIBUTING.md) for guidelines.

## References

### External Documentation

- [Rust Book](https://doc.rust-lang.org/book/)
- [PyO3 Guide](https://pyo3.rs/)
- [Rayon Documentation](https://docs.rs/rayon/)
- [SQLite Documentation](https://www.sqlite.org/docs.html)

### Scientific Background

- Kimura, M. (1983). The Neutral Theory of Molecular Evolution
- Hartl & Clark (2007). Principles of Population Genetics
- Nielsen & Slatkin (2013). An Introduction to Population Genetics

## License

MIT License - see [LICENSE](../LICENSE) file for details.
