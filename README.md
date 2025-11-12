# Centrevo

**High-performance simulator for centromeric sequence evolution**

[![Rust](https://img.shields.io/badge/rust-1.70%2B-orange.svg)](https://www.rust-lang.org/)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

Centrevo is a fast, flexible simulator for studying the evolution of tandemly repeated DNA sequences, with a focus on centromeric arrays. Built in Rust for maximum performance, it provides both command-line and Python interfaces for running population genetic simulations and comprehensive analysis.

## Key Features

- **High Performance**: Written in Rust with parallel computation support via Rayon
- **Realistic Models**: Mutation, recombination, and selection with customizable fitness functions
- **Persistent Storage**: SQLite-based recording with flexible strategies
- **Python Bindings**: Full Python API via PyO3 with PyArrow integration for efficient data export
- **Complete CLI**: Command-line interface for simulation management, execution, analysis, and export
- **Analysis Module**: Population genetics metrics including diversity, LD, distance, and composition analysis
- **Visualization**: Built-in plotting utilities with matplotlib for common analysis tasks
- **Well-Tested**: Comprehensive test coverage with unit and integration tests
- **Benchmarked**: Performance benchmarks included for optimization

## Installation

### Prerequisites

- **Rust** (1.70+): Install from [rust-lang.org](https://www.rust-lang.org/tools/install)
- **Cargo**: Comes with Rust

### Build from Source

#### macOS / Linux (Native)

```bash
git clone https://github.com/innan-lab/centrevo.git
cd centrevo
cargo build --release
```

The CLI binary will be available at `./target/release/centrevo`.

#### Linux (Cross-compilation using Docker)

If you don't have Rust installed or want to build in a containerized environment:

```bash
git clone https://github.com/innan-lab/centrevo.git
cd centrevo

# Build using the official Rust Docker image
docker run --rm -v "$(pwd)":/workspace -w /workspace rust:1.70 cargo build --release

# The binary will be in target/release/centrevo
./target/release/centrevo --version
```

For Python wheel building on Linux with Docker, see [PYTHON.md](PYTHON.md#building-for-distribution).

### Python Package

#### Local Development

To use the Python bindings locally, you'll need `maturin`:

```bash
pip install maturin
maturin develop --release
```

#### Build Python Wheels (Using Docker)

To build distributable Python wheels for Linux using Docker:

```bash
# Build for Python 3.12 on Linux x86_64
docker run --rm -v "$(pwd)":/io ghcr.io/pyo3/maturin build --release --target x86_64-unknown-linux-gnu -i python3.12

# Build for multiple Python versions
docker run --rm -v "$(pwd)":/io ghcr.io/pyo3/maturin build --release --target x86_64-unknown-linux-gnu -i python3.10 -i python3.11 -i python3.12

# Wheels will be in target/wheels/
ls target/wheels/
```

For more Python build options and usage documentation, see [PYTHON.md](PYTHON.md).

## Quick Start

### Command Line Interface

```bash
# Use the interactive setup wizard
./target/release/centrevo setup

# Or initialize and run manually
./target/release/centrevo init -N my_sim -n 100 -g 1000 --seed 42
./target/release/centrevo run -N my_sim
./target/release/centrevo analyze -N my_sim -g 1000
./target/release/centrevo export -N my_sim -g 1000 --format csv
```

See [CLI.md](CLI.md) for complete CLI documentation.

### Python API

```python
import centrevo

# Create simulation
alphabet = centrevo.Alphabet.dna()
structure = centrevo.RepeatStructure(
    alphabet=alphabet,
    init_base=centrevo.Nucleotide.A(),
    ru_length=171,
    rus_per_hor=12,
    hors_per_chr=100,
    chrs_per_hap=1
)

population = centrevo.create_initial_population(100, structure)

# Analyze
pi = centrevo.nucleotide_diversity(population, 0)
print(f"Nucleotide diversity: {pi:.6f}")
```

See [PYTHON.md](PYTHON.md) for complete Python API documentation.

## Architecture

Centrevo is organized into several modules:

- **`base`**: Core data structures (Nucleotide, Alphabet, Sequence)
- **`genome`**: Chromosome, Haplotype, and Individual representations
- **`evolution`**: Mutation, recombination, and fitness models
- **`simulation`**: Simulation engine and population dynamics
- **`storage`**: SQLite persistence with recording strategies
- **`analysis`**: Population genetics metrics (diversity, LD, distance, composition)
- **`python`**: Python bindings (optional, requires PyO3)

See [src/README.md](src/README.md) for detailed implementation documentation.

## Usage

For detailed usage documentation:

- **[CLI Guide](CLI.md)**: Complete command-line interface reference
- **[Python API](PYTHON.md)**: Python binding documentation and examples
- **[Implementation Guide](src/README.md)**: In-depth module implementation details
- **[Python Implementation](python/README.md)**: Python package implementation details
- **[Storage Module](src/storage/README.md)**: Database and persistence details
- **[Benchmarks](benches/README.md)**: Performance benchmarking guide

## Performance

Centrevo is designed for performance:

- **Parallel computation**: Uses Rayon for multi-core parallelism
- **Memory efficient**: Arc-based sequence sharing for efficient cloning
- **Optimized storage**: Binary BLOB storage with WAL mode for concurrent access
- **Benchmarked**: Comprehensive benchmark suite included

See [benches/README.md](benches/README.md) for benchmarking details and performance metrics.

## Testing

Centrevo has comprehensive test coverage:

```bash
# Run all tests
cargo test

# Run specific module tests
cargo test --lib storage
cargo test --lib analysis

# Run benchmarks
cargo bench
```

## Documentation

### User Documentation

- **[CLI Guide](CLI.md)**: Complete command-line interface reference
- **[Python API](PYTHON.md)**: Python binding documentation and examples

### Developer Documentation

- **[Implementation Guide](src/README.md)**: In-depth Rust module documentation
- **[Python Implementation](python/README.md)**: Python package internals
- **[Storage Module](src/storage/README.md)**: Database schema and persistence details
- **[Benchmarks](benches/README.md)**: Performance benchmarking guide
- **[Roadmap](ROADMAP.md)**: Development roadmap and future plans
- **[Changelog](CHANGELOG.md)**: Version history and changes
- **[Contributing](CONTRIBUTING.md)**: Contribution guidelines

Generate API documentation:
```bash
cargo doc --open
```

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

Built with:
- [Rust](https://www.rust-lang.org/) - Systems programming language
- [Rayon](https://github.com/rayon-rs/rayon) - Data parallelism
- [SQLite](https://www.sqlite.org/) via [rusqlite](https://github.com/rusqlite/rusqlite) - Database
- [PyO3](https://github.com/PyO3/pyo3) - Python bindings
- [PyArrow](https://arrow.apache.org/docs/python/) - Efficient data interchange
- [Criterion](https://github.com/bheisler/criterion.rs) - Benchmarking
- [Clap](https://github.com/clap-rs/clap) - CLI parsing
- [nalgebra](https://nalgebra.org/) - Linear algebra
- [statrs](https://github.com/statrs-dev/statrs) - Statistical functions

## Contact

- **Issues**: [GitHub Issues](https://github.com/innan-lab/centrevo/issues)
- **Discussions**: [GitHub Discussions](https://github.com/innan-lab/centrevo/discussions)

---

**Version**: 0.2.0  
**Status**: Active development  
**Rust Version**: 1.70+
