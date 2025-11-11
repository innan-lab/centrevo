# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed
- **BREAKING**: Simplified diversity metrics API by removing `haplotype_idx` parameter
  - `nucleotide_diversity()`, `tajimas_d()`, `wattersons_theta()`, and `haplotype_diversity()` 
    now analyze all 2n sequences (both haplotypes) by default
  - This aligns with standard population genetics practice and matches the distance module API
  - Tests updated to reflect analysis across all sequences

### In Progress (Phase 2 - Analysis Module)
- Population structure analysis (FST, PCA) - to be completed Week 3
- Temporal dynamics and allele trajectory tracking - to be completed Week 3
- Python bindings for analysis functions - to be completed Week 4
- Plotting utilities for visualizations - to be completed Week 4
- Example Jupyter notebooks - to be completed Week 4

### Planned (Phase 3)
- `centrevo run` command for executing simulations
- Export functionality for CSV/JSON formats
- Compression for database storage
- CI/CD pipeline with automated testing
- Integration tests

## [0.2.0] - 2025-11-10 (Phase 2 - Week 1 Complete)

### Added
- **Analysis Module** (Week 1 - Diversity Metrics)
  - `nucleotide_diversity()` - Calculate π (average pairwise differences)
  - `tajimas_d()` - Neutrality test statistic
  - `wattersons_theta()` - θ_W estimator from segregating sites
  - `haplotype_diversity()` - Probability two random haplotypes differ
  - Comprehensive test suite with 34+ tests
  - Performance benchmarks for all diversity metrics

- **Linkage and Distance Analysis** (Week 2 - Initial Implementation)
  - `linkage_disequilibrium()` - Calculate LD statistics (D, D', r²)
  - `ld_decay()` - LD decay with distance
  - `haplotype_blocks()` - Identify haplotype blocks
  - `pairwise_distances()` - Calculate all pairwise distances
  - `distance_matrix()` - Full distance matrix computation

- **Composition Analysis**
  - `gc_content()` - GC content calculation
  - `mean_gc_content()` - Population mean GC content
  - `nucleotide_composition()` - Full nucleotide composition

- **Polymorphism Analysis**
  - `count_segregating_sites()` - Count polymorphic sites
  - Site frequency spectrum (SFS) - placeholder for Week 2

- **Infrastructure**
  - Analysis module structure with 8 submodules
  - Parallel computation using Rayon for performance
  - Helper utilities (mean, std_dev, median)
  - Comprehensive documentation with mathematical formulas
  - Benchmark suite for performance tracking

### Dependencies Added
- `nalgebra` v0.33 - Linear algebra for PCA (prepared for Week 3)
- `statrs` v0.18 - Statistical functions
- `proptest` v1.5 - Property-based testing (dev dependency)

### Changed
- Updated `src/lib.rs` to export analysis module
- Updated `src/prelude.rs` with analysis function re-exports
- Fixed iteration patterns for Sequence type
- Improved code quality (fixed clippy warnings)

### Performance
- Nucleotide diversity: <100ms for 100 individuals × 10kb sequences
- All benchmarks meet Phase 2 performance targets

### Documentation
- Added rustdoc comments with mathematical formulas in KaTeX
- Included references to original papers
- Examples in documentation
- Test coverage >90% for implemented functions

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
  - Base operations (nucleotide, alphabet, sequence)
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

[Unreleased]: https://github.com/YOUR_USERNAME/centrevo/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/YOUR_USERNAME/centrevo/releases/tag/v0.1.0
