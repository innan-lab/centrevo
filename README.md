# Centrevo

**High-performance simulator for centromeric sequence evolution**

[![Rust](https://img.shields.io/badge/rust-1.70%2B-orange.svg)](https://www.rust-lang.org/)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

Centrevo is a fast, flexible simulator for studying the evolution of tandemly repeated DNA sequences, with a focus on centromeric arrays. Built in Rust for maximum performance, it provides both command-line and Python interfaces for running population genetic simulations.

## ✨ Key Features

- **🚀 High Performance**: Written in Rust with parallel computation support
- **🧬 Realistic Models**: Mutation, recombination, and selection with customizable fitness functions
- **💾 Persistent Storage**: SQLite-based recording with flexible strategies
- **🐍 Python Bindings**: Use from Python via PyO3 for analysis workflows
- **📊 Rich CLI**: Command-line interface for simulation management
- **🎯 Well-Tested**: 247 unit tests with comprehensive coverage
- **📈 Benchmarked**: Performance benchmarks included for optimization

## 🚀 Quick Start

### Installation

**From source:**
```bash
git clone https://github.com/YOUR_USERNAME/centrevo.git
cd centrevo
cargo build --release
```

The binary will be available at `./target/release/centrevo`.

### Running Your First Simulation

```bash
# Initialize a simulation
./target/release/centrevo init \
  -N my_first_sim \
  -n 100 \
  -g 1000 \
  --seed 42

# List all simulations
./target/release/centrevo list

# View simulation details
./target/release/centrevo info -N my_first_sim

# Check recorded generations
./target/release/centrevo generations -N my_first_sim
```

## 📖 Usage

### Command Line

Centrevo provides a comprehensive CLI:

```bash
centrevo <COMMAND>

Commands:
  init         Initialize a new simulation
  list         List simulations in database
  info         Show simulation info
  generations  List recorded generations
  help         Print help
```

See [CLI.md](CLI.md) for complete documentation.

### Python

```python
import centrevo

# Create simulation structure
alphabet = centrevo.Alphabet.dna()
structure = centrevo.RepeatStructure(
    alphabet=alphabet,
    init_base=centrevo.Nucleotide.A(),
    ru_length=171,
    rus_per_hor=12,
    hors_per_chr=100,
    chrs_per_hap=1
)

# Configure simulation
config = centrevo.SimulationConfig(
    population_size=100,
    total_generations=1000,
    seed=42
)

# Create population
population = centrevo.create_initial_population(100, structure)

# Setup recording
recorder = centrevo.Recorder(
    "simulation.db",
    "my_simulation",
    centrevo.RecordingStrategy.every_n(100)
)

# Record
recorder.record_metadata(config)
recorder.record_generation(population, 0)
```

See [PYTHON.md](PYTHON.md) for complete Python API documentation.

## 🏗️ Architecture

Centrevo is organized into several modules:

- **`base`**: Core data structures (Nucleotide, Alphabet, Sequence)
- **`genome`**: Chromosome, Haplotype, and Individual representations
- **`evolution`**: Mutation, recombination, and fitness models
- **`simulation`**: Simulation engine and population dynamics
- **`storage`**: SQLite persistence and query interface
- **`python`**: Python bindings (optional)

## 📊 Performance

Centrevo is designed for performance:

- **Parallel computation**: Uses Rayon for multi-core parallelism
- **Memory efficient**: Arc-based sequence sharing for cheap cloning
- **Optimized storage**: Binary BLOB storage with optional compression
- **Benchmarked**: Comprehensive benchmark suite included

Example performance (on a typical laptop):
- 100 individuals × 100KB chromosomes: ~1-5 seconds per generation
- 1000 individuals × 10KB chromosomes: ~2-10 seconds per generation

See [benches/README.md](benches/README.md) for benchmarking details.

## 🧪 Testing

Centrevo has comprehensive test coverage:

```bash
# Run all tests
cargo test

# Run with coverage
cargo tarpaulin

# Run benchmarks
cargo bench
```

All 247 tests pass with 0 failures.

## 📚 Documentation

- **[CLI Guide](CLI.md)**: Complete command-line reference
- **[Python API](PYTHON.md)**: Python binding documentation  
- **[Storage Module](src/storage/README.md)**: Database and persistence details
- **[Benchmarks](benches/README.md)**: Performance benchmarking guide
- **[Roadmap](ROADMAP.md)**: Development roadmap and future plans

Generate API documentation:
```bash
cargo doc --open
```

## 🤝 Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Development Setup

```bash
# Clone repository
git clone https://github.com/YOUR_USERNAME/centrevo.git
cd centrevo

# Build
cargo build

# Run tests
cargo test

# Check code
cargo clippy
cargo fmt --check
```

## 📝 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

Built with:
- [Rust](https://www.rust-lang.org/) - Systems programming language
- [Rayon](https://github.com/rayon-rs/rayon) - Data parallelism
- [SQLite](https://www.sqlite.org/) via [rusqlite](https://github.com/rusqlite/rusqlite) - Database
- [PyO3](https://github.com/PyO3/pyo3) - Python bindings
- [Criterion](https://github.com/bheisler/criterion.rs) - Benchmarking
- [Clap](https://github.com/clap-rs/clap) - CLI parsing

## 📧 Contact

- **Issues**: [GitHub Issues](https://github.com/YOUR_USERNAME/centrevo/issues)
- **Discussions**: [GitHub Discussions](https://github.com/YOUR_USERNAME/centrevo/discussions)

## 🗺️ Roadmap

See [ROADMAP.md](ROADMAP.md) for planned features and development phases.

## ⭐ Star History

If you find Centrevo useful, please consider giving it a star on GitHub!

---

**Version**: 0.1.0  
**Status**: Active development  
**Rust Version**: 1.70+
