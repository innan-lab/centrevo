# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed (Breaking)
- Removed the `Unnormalized` / `Normalized` type-state pattern for `FitnessValue` and `LogFitnessValue`. Both are now single-type wrappers over `f64`:
  - Use `FitnessValue::new(f64)` and `LogFitnessValue::new(f64)` for construction.
  - `FitnessValue::LETHAL_FITNESS` = `0.0`, `FitnessValue::NEUTRAL_FITNESS` = `1.0`.
  - `LogFitnessValue::LETHAL_FITNESS` = `-∞` (negative infinity), `LogFitnessValue::NEUTRAL_FITNESS` = `0.0`.
  - The `normalize()` / `new_normalized()` / `unnormalize()` helper APIs have been removed. If you need to normalize by the sum of weights, do that manually (e.g., divide by the total or use log-space subtraction for `LogFitnessValue`).
  - Arithmetic: multiplication remains implemented via log-space addition to preserve numerical stability; addition is linear and the result may exceed `1.0`.

  Migration: update any usages such as `FitnessValue::<Normalized>::new_normalized(f)` to `FitnessValue::new(f)`, and replace `.normalize()` invocations with explicit `value / total` or log-space subtraction as appropriate.


### Planned (Phase 3 - Next Release)
- CI/CD pipeline with automated testing and code coverage tracking
- Plotting utilities module with matplotlib integration
- Example Jupyter notebooks demonstrating analysis workflows
- Documentation improvements and analysis guide

## [0.2.1] - 2025-11-13

### Added

#### Clippy Linting Improvements
- Fixed 81+ clippy warnings for improved code quality
- Updated all format strings to use inline variable syntax (`format!("{variable}")`)
- Improved loop patterns using iterators instead of range-based indexing

### Fixed
- **Code Quality**: Eliminated all clippy warnings via automatic fixes and manual refactoring
  - Replaced needless range loops with iterator-based patterns
  - Optimized format strings with inline variables (81 fixes in bin and lib)
  - Applied `#[allow(clippy::too_many_arguments)]` to complex functions requiring many parameters
  - `src/evolution/mutation.rs`: Fixed loop variable indexing pattern
  - `src/storage/recorder.rs`: Allowed too_many_arguments for database operations
  - `src/bin/centrevo.rs`: Fixed format string inlining and function signatures

### Changed
- Enhanced code maintainability through consistent formatting practices
- Improved compiler warnings hygiene (now passes `cargo clippy` with no warnings)

### Documentation
- All documentation remains consistent with v0.2.0 feature set
- Changelog updated to reflect release preparation

### Testing
- All existing tests pass without modification
- Code quality improvements do not affect functionality

## [0.2.0] - 2025-11-12

### Added

#### Complete CLI Suite
- **`centrevo init`** - Initialize simulation with comprehensive parameters
  - Population size, generations, repeat structure configuration
  - Random seed support for reproducibility
  - Creates simulation database with metadata

- **`centrevo run`** - Execute simulations with full control
  - Resume from checkpoint support for interrupted simulations
  - Configurable mutation rate, recombination parameters
  - Progress bar display for long-running simulations
  - Flexible recording intervals

- **`centrevo export`** - Export data in standard formats
  - CSV export for sequences and fitness data
  - JSON export for metadata and configurations
  - FASTA export for sequence data
  - Three data types: sequences, metadata, fitness

- **`centrevo analyze`** - Comprehensive analysis from command line
  - Diversity metrics (π, Tajima's D, θ_W, haplotype diversity)
  - Linkage disequilibrium analysis
  - Distance matrices and pairwise distances
  - GC content and nucleotide composition
  - Output in pretty-printed or JSON format

- **`centrevo validate`** - Database integrity checking
  - Verify simulation completeness
  - Check for missing generations
  - Database health validation
  - Optional fix mode for repair

- **`centrevo setup`** - Interactive configuration wizard
  - Step-by-step parameter selection
  - Intelligent defaults
  - Complete workflow from setup to execution

#### Analysis Module (Complete Implementation)

**Diversity Metrics** (`src/analysis/diversity.rs`)
- `nucleotide_diversity()` - Calculate π (average pairwise differences)
- `tajimas_d()` - Neutrality test statistic
- `wattersons_theta()` - θ_W estimator from segregating sites
- `haplotype_diversity()` - Probability two random haplotypes differ
- Parallel computation using Rayon for performance
- Comprehensive test coverage (305 tests total, 100% pass rate)

**Linkage Disequilibrium Analysis** (`src/analysis/linkage.rs`)
- `linkage_disequilibrium()` - Calculate D, D', and r² statistics
- `ld_decay()` - LD decay patterns with distance
- `haplotype_blocks()` - Identify regions of high LD
- Gabriel et al. method for block identification

**Distance and Composition** (`src/analysis/distance.rs`, `src/analysis/composition.rs`)
- `pairwise_distances()` - Hamming distances for all sequence pairs
- `distance_matrix()` - Full n×n symmetric distance matrix
- `gc_content()` - Per-sequence and population GC content
- `nucleotide_composition()` - Full ACGT composition analysis

**Polymorphism Analysis** (`src/analysis/polymorphism.rs`)
- `count_segregating_sites()` - Count polymorphic sites
- `site_frequency_spectrum()` - SFS calculation (placeholder for future enhancement)

**Population Structure** (`src/analysis/structure.rs`)
- `fst()` - FST calculation (stub for future implementation)
- `pca()` - Principal component analysis (stub for future implementation)

**Temporal Analysis** (`src/analysis/temporal.rs`)
- `allele_trajectory()` - Track allele frequencies over time (stub)
- `fitness_dynamics()` - Mean fitness tracking (stub)

#### Storage and Recording Enhancements
- **Async Recording** - Non-blocking database writes
  - Tokio-based async recorder with buffering
  - Background thread for I/O operations
  - Significant performance improvement for large populations
  - Files: `src/storage/recorder.rs`

- **Checkpoint/Resume Functionality**
  - Save and restore simulation state
  - Resume interrupted simulations seamlessly
  - Checkpoint metadata tracking
  - Comprehensive tests for checkpoint reliability

#### Performance Optimizations
- **Xoshiro256++ RNG** - Replaced StdRng throughout
  - 20-30% performance improvement in mutation operations
  - Better statistical properties
  - Faster generation of random numbers

- **Poisson Pre-sampling** - Mutation count optimization
  - Pre-compute mutation counts using Poisson distribution
  - Reduces per-site overhead in mutation operations
  - 15-25% speedup in mutation-heavy simulations

- **Parallelized Simulation** - Multi-threaded processing
  - Parallel fitness evaluation
  - Parallel distance calculations
  - Optimized memory usage with Arc-based sharing

- **Fast Sequence Operations** - Optimized core algorithms
  - Hamming distance with chunked processing (83-96% faster)
  - Direct index access patterns
  - Cache-friendly memory access
  - Nucleotide diversity: 20.8ms → 1.34ms (15.5x speedup)
  - Tajima's D: 21.8ms → 1.36ms (16x speedup)

#### Testing and Quality
- **Comprehensive Test Suite**
  - 305 total tests with 100% pass rate
  - 41 analysis module tests
  - Property-based testing with proptest
  - Integration tests for end-to-end workflows
  - Storage and serialization tests
  - ~90% code coverage for core functionality

- **Benchmark Suite** - Performance tracking
  - Criterion-based benchmarks for all modules
  - Analysis benchmarks for diversity, LD, distances
  - Storage operation benchmarks
  - Regression detection for performance

### Changed
- **BREAKING**: Diversity metrics API simplified
  - Removed `haplotype_idx` parameter from diversity functions
  - `nucleotide_diversity()`, `tajimas_d()`, `wattersons_theta()`, `haplotype_diversity()`
  - Now analyze all 2n sequences (both haplotypes) by default
  - Aligns with standard population genetics practice

- **Performance**: Significant improvements across the board
  - 15-22x speedup in diversity metric calculations
  - 20-30% improvement in mutation operations
  - Async recording for non-blocking I/O

- **Documentation**: Enhanced rustdoc with mathematical formulas
  - KaTeX-formatted equations throughout analysis module
  - References to original papers
  - Comprehensive examples in doc comments

### Dependencies Added
- `nalgebra` v0.33 - Linear algebra for future PCA implementation
- `statrs` v0.18 - Statistical functions for analysis
- `proptest` v1.5 - Property-based testing (dev dependency)
- `tokio` v1.41 - Async runtime for async recording
- `zstd` v0.13 - Compression for storage optimization

### Performance Metrics
- **Nucleotide diversity**: <2ms for 100 individuals × 1kb sequences
- **Tajima's D**: <2ms for typical datasets
- **LD calculations**: Efficient parallel implementation
- **Database writes**: Non-blocking with async recorder
- **Memory**: Arc-based sharing for efficient clones

### Documentation
- Updated CLI.md with all command examples
- Enhanced PYTHON.md with usage patterns
- Comprehensive rustdoc with mathematical formulas
- Storage README with database schema details
- Benchmark suite documentation

### Known Limitations
- Population structure analysis (FST, PCA) - stubs only, full implementation pending
- Temporal analysis functions - stubs only, full implementation pending
- Python bindings for analysis functions - not yet exposed
- Plotting utilities - planned for Phase 3
- Compression not yet implemented (zstd dependency added, implementation pending)

### Migration Notes
If upgrading from v0.1.x:
- Update diversity metric calls to remove `haplotype_idx` parameter
- Analysis now automatically considers all 2n sequences
- Database format is backward compatible
- No changes required to simulation configurations

## [0.1.1] - 2025-11-10 (Phase 1 Complete)

### Added
- Comprehensive README.md with quick start guide and examples
- MIT LICENSE file
- CONTRIBUTING.md with detailed contribution guidelines
- CHANGELOG.md following Keep a Changelog format

### Changed
- Implemented `Display` trait for `Sequence`, `SharedSequence`, and `Chromosome`
- Updated benchmarks to use `std::hint::black_box` instead of deprecated `criterion::black_box`
- Improved code quality with clippy lint fixes

### Development
- Phase 1 of roadmap completed (code quality and documentation)
- Established professional open-source project structure
- Ready for Phase 2 (analysis module implementation)

## [0.1.0] - 2025-11-10

### Added
- **Core Modules**
  - Base module with Nucleotide, Alphabet, and Sequence types
  - Genome module with Chromosome, Haplotype, and Individual representations
  - Evolution module with mutation, recombination, and selection models
  - Simulation engine with population dynamics
  - Storage module with SQLite persistence

- **CLI Interface**
  - `centrevo init` command to initialize simulations
  - `centrevo list` command to list simulations
  - `centrevo info` command to show simulation details
  - `centrevo generations` command to list recorded generations

- **Python Bindings**
  - PyO3-based Python interface
  - Full API access to core functionality
  - See PYTHON.md for documentation

- **Features**
  - Parallel computation using Rayon
  - Memory-efficient Arc-based sequence sharing
  - Flexible recording strategies (EveryN, Specific, All, None)
  - Comprehensive fitness models (GC content, length, similarity)
  - Mutation models (JC69-style substitution)
  - Recombination (crossover and gene conversion)

- **Testing & Quality**
  - 247 unit tests with 100% pass rate
  - Comprehensive benchmark suite using Criterion
  - Clean clippy lints (0 warnings)
  - Well-documented code with examples

- **Documentation**
  - README.md with quick start guide
  - CLI.md with complete command reference
  - PYTHON.md with Python API documentation
  - Storage module README with database details
  - Benchmarks README with performance guide
  - ROADMAP.md with development plans
  - CONTRIBUTING.md with contribution guidelines

### Implementation Details
- **Base Module** (~700 LOC)
  - Nucleotide operations with complement and ASCII conversion
  - Alphabet with character/index mapping
  - Mutable and shared sequence types
  - Display trait implementations for string conversion

- **Genome Module** (~1,500 LOC)
  - Chromosome with repeat structure metadata
  - Haplotype as collection of chromosomes
  - Diploid Individual representation
  - Shared (Arc-based) and mutable variants

- **Evolution Module** (~600 LOC)
  - Substitution model with configurable rates
  - Recombination with crossover and gene conversion
  - Multiple fitness functions (GC, length, similarity)
  - Selection mechanisms for reproduction

- **Simulation Module** (~800 LOC)
  - Main simulation engine
  - Population management
  - Parameter configurations (mutation, recombination, fitness)
  - Generation stepping and full simulation runs

- **Storage Module** (~1,200 LOC)
  - SQLite database with WAL mode
  - Recording strategies for flexible data capture
  - Query builder for data retrieval
  - Fitness history tracking
  - Individual snapshots at specific generations

- **Benchmarks** (~1,200 LOC)
  - Base operations (nucleotide, sequence)
  - Genome operations (chromosome, haplotype, individual)
  - Evolution operations (mutation, recombination, fitness)
  - Full simulation benchmarks
  - Scaling tests for various population sizes

### Performance Characteristics
- Sequence operations: O(n) with length
- Population fitness: Parallelized with Rayon
- Database writes: Batched transactions
- Memory: Arc-based sharing for cheap clones
- Typical throughput: 1-10 seconds per generation (100-1000 individuals)

### Known Limitations
- CLI can only initialize, not run simulations (use Python/Rust API)
- No compression for database storage
- Limited analysis capabilities (add in 0.2.0)
- No checkpoint/resume functionality
- Single-threaded simulation loop (parallelized operations)

### Breaking Changes
None (initial release)

### Security
- No known security issues
- SQLite database uses safe bindings
- No unsafe Rust code in core modules

## Release Notes Template (for future releases)

```markdown
## [X.Y.Z] - YYYY-MM-DD

### Added
- New features

### Changed
- Changes to existing functionality

### Deprecated
- Soon-to-be removed features

### Removed
- Removed features

### Fixed
- Bug fixes

### Security
- Security fixes
```

---

[Unreleased]: https://github.com/innan-lab/centrevo/compare/v0.2.1...HEAD
[0.2.1]: https://github.com/innan-lab/centrevo/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/innan-lab/centrevo/compare/v0.1.1...v0.2.0
[0.1.1]: https://github.com/innan-lab/centrevo/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/innan-lab/centrevo/releases/tag/v0.1.0
