# Python Bindings for Centrevo

Centrevo provides Python bindings using PyO3, allowing you to use the high-performance Rust library from Python.

## Installation

### Building from Source

```bash
# Install maturin for building Python extensions
pip install maturin

# Build and install the Python package
maturin develop --release

# Or build a wheel
maturin build --release
```

## Usage Example

```python
import centrevo

# Create DNA alphabet
alphabet = centrevo.Alphabet.dna()

# Create a nucleotide
base_a = centrevo.Nucleotide.A()

# Create repeat structure
structure = centrevo.RepeatStructure(
    alphabet=alphabet,
    init_base=base_a,
    ru_length=171,
    rus_per_hor=12,
    hors_per_chr=100,
    chrs_per_hap=1
)

print(f"Chromosome length: {structure.chr_length()} bp")

# Create simulation config
config = centrevo.SimulationConfig(
    population_size=100,
    total_generations=1000,
    seed=42
)

# Create initial population
population = centrevo.create_initial_population(100, structure)
print(f"Population size: {len(population)}")

# Setup recorder
recorder = centrevo.Recorder(
    "simulation.db",
    "my_simulation",
    centrevo.RecordingStrategy.every_n(100)
)

# Record metadata and first generation
recorder.record_metadata(config)
recorder.record_generation(population, 0)

print("✓ Simulation initialized!")

# Query results
query = centrevo.QueryBuilder("simulation.db")
simulations = query.list_simulations()
print(f"Simulations in database: {simulations}")

generations = query.get_recorded_generations("my_simulation")
print(f"Recorded generations: {generations}")
```

## API Reference

### Classes

- **`Nucleotide`**: DNA base (A, C, G, T)
  - Static methods: `A()`, `C()`, `G()`, `T()`
  
- **`Alphabet`**: Collection of nucleotide characters
  - Static method: `dna()` - Creates standard DNA alphabet
  
- **`Chromosome`**: Chromosome with repeat structure
  - Static method: `uniform(id, base, length, ru_length, rus_per_hor, alphabet)`
  
- **`Haplotype`**: Collection of chromosomes
  - Constructor: `Haplotype()`
  - Static method: `from_chromosomes(chromosomes)`
  
- **`Individual`**: Diploid individual
  - Constructor: `Individual(id, haplotype1, haplotype2)`
  
- **`Population`**: Collection of individuals
  - Constructor: `Population(id, individuals)`
  - Methods: `size()`, `generation()`
  
- **`RepeatStructure`**: Configuration for repeat structure
  - Constructor: `RepeatStructure(alphabet, init_base, ru_length, rus_per_hor, hors_per_chr, chrs_per_hap)`
  - Method: `chr_length()` - Get total chromosome length
  
- **`SimulationConfig`**: Simulation parameters
  - Constructor: `SimulationConfig(population_size, total_generations, seed)`
  
- **`RecordingStrategy`**: When to save snapshots
  - Static methods:
    - `every_n(n)` - Record every N generations
    - `specific(generations)` - Record specific generations
    - `all()` - Record all generations
    - `none()` - No recording
  
- **`Recorder`**: Save simulation state to database
  - Constructor: `Recorder(db_path, sim_id, strategy)`
  - Methods:
    - `record_metadata(config)` - Save configuration
    - `record_generation(population, generation)` - Save generation state
  
- **`QueryBuilder`**: Query saved simulations
  - Constructor: `QueryBuilder(db_path)`
  - Methods:
    - `list_simulations()` - Get all simulation names
    - `get_recorded_generations(sim_id)` - Get recorded generation numbers

### Functions

- **`create_initial_population(size, structure)`**: Create uniform initial population

## Building for Distribution

### Local Development

```bash
# Install in current environment for development
maturin develop --release
```

### Building for Current Platform

```bash
# Build for current platform (macOS ARM or Intel)
maturin build --release

# Build for multiple Python versions (requires them to be installed)
maturin build --release --find-interpreter
```

### Cross-Platform Builds

#### macOS Universal (Apple Silicon + Intel)

```bash
# Build a universal wheel that works on both ARM and Intel Macs
maturin build --release --target universal2-apple-darwin
```

First add the required targets:
```bash
rustup target add aarch64-apple-darwin
rustup target add x86_64-apple-darwin
```

#### Linux x86_64 (using Docker)

**Note:** Requires Docker to be installed and running.

```bash
# Build for Linux x86_64 with Python 3.12
docker run --rm -v "$(pwd)":/io ghcr.io/pyo3/maturin build --release --target x86_64-unknown-linux-gnu -i python3.12

# Build for multiple Python versions
docker run --rm -v "$(pwd)":/io ghcr.io/pyo3/maturin build --release --target x86_64-unknown-linux-gnu -i python3.10 -i python3.11 -i python3.12
```

First add the Linux target:
```bash
rustup target add x86_64-unknown-linux-gnu
```

#### Build All Platforms

```bash
# macOS universal
maturin build --release --target universal2-apple-darwin

# Linux x86_64 (via Docker)
docker run --rm -v "$(pwd)":/io ghcr.io/pyo3/maturin build --release --target x86_64-unknown-linux-gnu -i python3.10 -i python3.11 -i python3.12
```

All wheels will be created in `target/wheels/`.

## Performance Notes

- All core computations are done in Rust for maximum performance
- Database operations use SQLite with WAL mode for concurrent access
- Large populations can be processed in parallel using Rayon

## Requirements

- Python 3.8+
- Rust toolchain (for building from source)
- maturin (for building Python packages)
