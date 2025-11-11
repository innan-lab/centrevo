# Python Bindings for Centrevo

Centrevo provides comprehensive Python bindings using PyO3, allowing you to use the high-performance Rust library from Python with full access to simulation, analysis, and visualization capabilities.

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

### Installing Dependencies

```bash
# Core dependencies (automatically installed)
pip install cffi pyarrow matplotlib numpy

# Optional dev dependencies
pip install pytest jupyter pandas polars
```

## Quick Start

```python
import centrevo

# Create and run simulation
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
tajima_d = centrevo.tajimas_d(population, 0)
gc = centrevo.gc_content(population, None, None, None)

print(f"Nucleotide diversity: {pi:.6f}")
print(f"Tajima's D: {tajima_d:.4f}")
print(f"GC content: {gc:.2%}")
```

## API Reference

### Core Classes

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

### Analysis Functions

#### Diversity Metrics

- **`nucleotide_diversity(population, chromosome_idx)`**: Calculate π (nucleotide diversity)
  - Returns: Average pairwise differences per site across all 2n sequences
  - Example: `pi = centrevo.nucleotide_diversity(pop, 0)`

- **`tajimas_d(population, chromosome_idx)`**: Calculate Tajima's D statistic
  - Returns: Neutrality test statistic (positive = balancing selection, negative = purifying selection)
  - Example: `d = centrevo.tajimas_d(pop, 0)`

- **`wattersons_theta(population, chromosome_idx)`**: Calculate Watterson's θ
  - Returns: θ_W estimator from segregating sites
  - Example: `theta = centrevo.wattersons_theta(pop, 0)`

- **`haplotype_diversity(population, chromosome_idx)`**: Calculate haplotype diversity
  - Returns: Probability that two random haplotypes differ (0.0 to 1.0)
  - Example: `h = centrevo.haplotype_diversity(pop, 0)`

#### Linkage Disequilibrium

- **`linkage_disequilibrium(population, pos1, pos2, chromosome_idx, haplotype_idx)`**: Calculate LD statistics
  - Returns: Dict with keys `'D'`, `'D_prime'`, `'r_squared'`, or `None` if cannot calculate
  - Example: `ld = centrevo.linkage_disequilibrium(pop, 100, 200, 0, 0)`

- **`ld_decay(population, chromosome_idx, haplotype_idx, max_distance, bin_size)`**: Calculate LD decay
  - Returns: Dict with `'distances'` and `'r_squared_values'` lists
  - Example: `decay = centrevo.ld_decay(pop, 0, 0, 1000, 10)`

- **`haplotype_blocks(population, chromosome_idx, haplotype_idx, r_squared_threshold=0.8)`**: Identify haplotype blocks
  - Returns: List of tuples `(start, end)` representing blocks
  - Example: `blocks = centrevo.haplotype_blocks(pop, 0, 0, 0.8)`

#### Distance Analysis

- **`pairwise_distances(population, chromosome_idx)`**: Calculate all pairwise distances
  - Returns: List of normalized Hamming distances between all sequence pairs
  - Example: `distances = centrevo.pairwise_distances(pop, 0)`

- **`distance_matrix(population, chromosome_idx)`**: Calculate full distance matrix
  - Returns: 2D list (2n × 2n matrix) of pairwise distances
  - Example: `matrix = centrevo.distance_matrix(pop, 0)`

#### Composition Analysis

- **`gc_content(population, individual_idx, haplotype_idx, chromosome_idx)`**: Calculate GC content
  - Flexible based on which parameters are provided:
    - All None: population-level average
    - individual_idx only: average across both haplotypes
    - individual_idx + haplotype_idx: specific haplotype average
    - All three: specific chromosome
  - Returns: GC content (0.0 to 1.0)
  - Examples:
    - `gc = centrevo.gc_content(pop, None, None, None)  # Population mean`
    - `gc = centrevo.gc_content(pop, 0, None, None)     # Individual 0 mean`
    - `gc = centrevo.gc_content(pop, 0, 1, 0)           # Specific chromosome`

- **`nucleotide_composition(population, individual_idx, haplotype_idx, chromosome_idx)`**: Calculate nucleotide composition
  - Same flexibility as `gc_content()`
  - Returns: Dict mapping nucleotide (str) to frequency (float)
  - Example: `comp = centrevo.nucleotide_composition(pop, None, None, None)`

#### Polymorphism Analysis

- **`count_segregating_sites(population, chromosome_idx, haplotype_idx)`**: Count segregating sites
  - Returns: Number of polymorphic sites
  - Example: `s = centrevo.count_segregating_sites(pop, 0, 0)`

### Export Functions (PyArrow Integration)

- **`export_diversity_metrics(population, chromosome_idx)`**: Export diversity metrics as dict
  - Returns: Dict with keys: `'nucleotide_diversity'`, `'tajimas_d'`, `'wattersons_theta'`, 
    `'haplotype_diversity'`, `'generation'`, `'population_size'`
  - Example: `metrics = centrevo.export_diversity_metrics(pop, 0)`

- **`export_distance_matrix(population, chromosome_idx)`**: Export distance matrix as list of dicts
  - Returns: List of dicts with keys `'sequence_i'`, `'sequence_j'`, `'distance'`
  - Example: `dist_data = centrevo.export_distance_matrix(pop, 0)`

- **`export_ld_decay(population, chromosome_idx, haplotype_idx, max_distance, bin_size)`**: Export LD decay data
  - Returns: List of dicts with keys `'distance'`, `'r_squared'`
  - Example: `ld_data = centrevo.export_ld_decay(pop, 0, 0, 1000, 10)`

### Plotting Functions

The `centrevo.plotting` module provides visualization utilities:

- **`plot_diversity_trajectory(diversity_data, metric='nucleotide_diversity', ...)`**: Plot diversity over time
- **`plot_ld_decay(ld_data, ...)`**: Plot LD decay curve
- **`plot_distance_matrix(distance_matrix, ...)`**: Plot distance matrix heatmap
- **`plot_nucleotide_composition(composition, ...)`**: Plot nucleotide frequencies
- **`plot_multiple_diversity_metrics(diversity_data, metrics=None, ...)`**: Plot multiple metrics

### Helper Functions

- **`create_initial_population(size, structure)`**: Create uniform initial population

## Complete Example

```python
import centrevo
import centrevo.plotting as cplt

# 1. Create simulation
alphabet = centrevo.Alphabet.dna()
structure = centrevo.RepeatStructure(
    alphabet=alphabet,
    init_base=centrevo.Nucleotide.A(),
    ru_length=171,
    rus_per_hor=12,
    hors_per_chr=100,
    chrs_per_hap=1
)

config = centrevo.SimulationConfig(
    population_size=100,
    total_generations=1000,
    seed=42
)

# 2. Initialize population
population = centrevo.create_initial_population(100, structure)

# 3. Setup recording
recorder = centrevo.Recorder(
    "my_simulation.db",
    "experiment1",
    centrevo.RecordingStrategy.every_n(100)
)

recorder.record_metadata(config)
recorder.record_generation(population, 0)

# 4. Analyze population
print("Population Analysis")
print("=" * 50)

# Diversity metrics
pi = centrevo.nucleotide_diversity(population, 0)
tajima = centrevo.tajimas_d(population, 0)
theta = centrevo.wattersons_theta(population, 0)
hap_div = centrevo.haplotype_diversity(population, 0)

print(f"Nucleotide diversity (π): {pi:.6f}")
print(f"Watterson's θ: {theta:.6f}")
print(f"Tajima's D: {tajima:.4f}")
print(f"Haplotype diversity: {hap_div:.4f}")

# Composition
gc = centrevo.gc_content(population, None, None, None)
comp = centrevo.nucleotide_composition(population, None, None, None)

print(f"\nGC content: {gc:.2%}")
print("Nucleotide composition:")
for nuc, freq in sorted(comp.items()):
    print(f"  {nuc}: {freq:.4f}")

# Polymorphism
seg_sites_h1 = centrevo.count_segregating_sites(population, 0, 0)
seg_sites_h2 = centrevo.count_segregating_sites(population, 0, 1)

print(f"\nSegregating sites (haplotype 1): {seg_sites_h1}")
print(f"Segregating sites (haplotype 2): {seg_sites_h2}")

# 5. Export and visualize
print("\n" + "=" * 50)
print("Exporting data for visualization...")

# Export diversity metrics (for trajectory plots)
metrics = centrevo.export_diversity_metrics(population, 0)

# Export LD decay
ld_data = centrevo.export_ld_decay(population, 0, 0, 1000, 10)

# Export distance matrix
dist_matrix = centrevo.distance_matrix(population, 0)

# Create plots
fig_ld = cplt.plot_ld_decay(ld_data, save_path="ld_decay.png")
fig_dist = cplt.plot_distance_matrix(dist_matrix, save_path="distances.png")
fig_comp = cplt.plot_nucleotide_composition(comp, save_path="composition.png")

print("Plots saved!")

# 6. PyArrow integration for pandas/polars
import pyarrow as pa
import pandas as pd

# Convert LD data to PyArrow table
table = pa.Table.from_pylist(ld_data)

# Convert to pandas
df = table.to_pandas()
print(f"\nLD decay dataframe shape: {df.shape}")
print(df.head())

# Or use polars
try:
    import polars as pl
    df_polars = pl.from_arrow(table)
    print(f"\nPolars dataframe shape: {df_polars.shape}")
except ImportError:
    print("\nPolars not installed (optional)")

print("\nAnalysis complete!")
```

## Usage Patterns

### Pattern 1: Time-series Analysis

```python
import centrevo
import centrevo.plotting as cplt

# Collect diversity metrics over time
diversity_trajectory = []

# Load populations from database for each generation
query = centrevo.QueryBuilder("simulation.db")
generations = query.get_recorded_generations("my_sim")

for gen in generations:
    # Load population (implement your loading logic)
    population = load_population_from_db("my_sim", gen)
    
    # Export metrics
    metrics = centrevo.export_diversity_metrics(population, 0)
    diversity_trajectory.append(metrics)

# Plot trajectory
fig = cplt.plot_diversity_trajectory(
    diversity_trajectory,
    metric='nucleotide_diversity',
    save_path='diversity_trajectory.png'
)

# Plot multiple metrics
fig_multi = cplt.plot_multiple_diversity_metrics(
    diversity_trajectory,
    metrics=['nucleotide_diversity', 'tajimas_d', 'haplotype_diversity'],
    save_path='multi_metrics.png'
)
```

### Pattern 2: Comparative Analysis

```python
import centrevo
import pandas as pd

def analyze_population(pop, name):
    """Analyze a population and return summary dict."""
    return {
        'name': name,
        'pi': centrevo.nucleotide_diversity(pop, 0),
        'tajima_d': centrevo.tajimas_d(pop, 0),
        'theta_w': centrevo.wattersons_theta(pop, 0),
        'hap_div': centrevo.haplotype_diversity(pop, 0),
        'gc_content': centrevo.gc_content(pop, None, None, None),
        'seg_sites': centrevo.count_segregating_sites(pop, 0, 0),
    }

# Analyze multiple simulations
results = []
for sim_name in ['sim1', 'sim2', 'sim3']:
    pop = load_population(sim_name, generation=1000)
    results.append(analyze_population(pop, sim_name))

# Create comparison dataframe
df = pd.DataFrame(results)
print(df)

# Statistical comparison
print("\nSummary statistics:")
print(df.describe())
```

### Pattern 3: Spatial Analysis

```python
import centrevo
import numpy as np
import matplotlib.pyplot as plt

# Calculate LD across multiple distance bins
ld_results = []
for start_pos in range(0, 10000, 1000):
    ld_data = centrevo.export_ld_decay(
        population, 0, 0,
        max_distance=1000,
        bin_size=50
    )
    ld_results.append((start_pos, ld_data))

# Analyze LD patterns along chromosome
# ... custom analysis code ...

# Calculate distance distribution
distances = centrevo.pairwise_distances(population, 0)
plt.hist(distances, bins=50)
plt.xlabel('Genetic Distance')
plt.ylabel('Frequency')
plt.title('Distribution of Pairwise Distances')
plt.savefig('distance_distribution.png')
```

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
- Analysis functions use parallel computation via Rayon when beneficial
- Database operations use SQLite with WAL mode for concurrent access
- PyArrow provides zero-copy data transfer between Rust and Python
- Large populations (>1000 individuals) analyze in seconds

## Requirements

- Python 3.8+
- Rust toolchain (for building from source)
- maturin (for building Python packages)

## Python Dependencies

Core (required):
- `cffi` - Foreign function interface
- `pyarrow>=14.0.0` - Efficient data interchange
- `matplotlib>=3.5.0` - Plotting
- `numpy>=1.21.0` - Numerical operations

Optional (development):
- `pytest>=7.0.0` - Testing
- `jupyter>=1.0.0` - Notebooks
- `pandas>=1.5.0` - Data analysis
- `polars>=0.19.0` - Fast dataframes
