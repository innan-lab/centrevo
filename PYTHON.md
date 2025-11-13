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
import centrevo as cv

# Simple simulation with builder pattern (recommended)
sim = (cv.SimulationBuilder()
    .population_size(100)
    .generations(1000)
    .repeat_structure(ru_length=171, rus_per_hor=12, hors_per_chr=100)
    .mutation_rate(0.0001)
    .seed(42)
    .build())

# Run simulation
sim.run()

# Analyze final population
population = sim.population()
pi = cv.nucleotide_diversity(population, 0)
tajima_d = cv.tajimas_d(population, 0)
gc = cv.gc_content(population, None, None, None)

print(f"Nucleotide diversity: {pi:.6f}")
print(f"Tajima's D: {tajima_d:.4f}")
print(f"GC content: {gc:.2%}")
```

## API Reference

### Simulation Builder (Recommended)

The **`SimulationBuilder`** provides a fluent API for creating simulations with sensible defaults and clear validation.

#### Basic Usage

```python
import centrevo as cv

# Minimal example - defaults to uniform(A), no mutation, no recombination
sim = (cv.SimulationBuilder()
    .population_size(50)
    .generations(100)
    .repeat_structure(ru_length=171, rus_per_hor=12, hors_per_chr=10)
    .build())

# With evolutionary parameters
sim = (cv.SimulationBuilder()
    .population_size(50)
    .generations(100)
    .repeat_structure(ru_length=171, rus_per_hor=12, hors_per_chr=10)
    .mutation_rate(0.0001)
    .recombination(break_prob=0.01, crossover_prob=0.7, gc_extension_prob=0.1)
    .seed(42)
    .build())

# Run the simulation
sim.run()
```

#### Builder Methods

**Required Parameters:**
- `population_size(size)` - Number of individuals (required)
- `generations(n)` - Number of generations to run (required)

**Repeat Structure (required for random/uniform initialization):**
- `repeat_structure(ru_length, rus_per_hor, hors_per_chr)` - Define repeat structure
- `chromosomes_per_haplotype(n)` - Number of chromosomes per haplotype (default: 1)

**Initialization Methods:**
- `init_uniform(base)` - Initialize with uniform sequence (default: Nucleotide.A())
- `init_random()` - Initialize with random bases from alphabet
- `init_from_fasta(path)` - Load sequences from FASTA file
- `init_from_json(input)` - Load sequences from JSON file or string
- `init_from_checkpoint(db_path, sim_id, generation=None)` - Load sequences from checkpoint

**Evolutionary Parameters (all optional):**
- `alphabet(alphabet)` - Set alphabet (default: Alphabet.dna())
- `mutation_rate(rate)` - Set mutation rate (default: 0.0)
- `recombination(break_prob, crossover_prob, gc_extension_prob)` - Set recombination rates (default: no recombination)
- `fitness(fitness_config)` - Set fitness/selection configuration (default: neutral, see Fitness Configuration below)
- `seed(seed)` - Set random seed (default: None = random)

**Build:**
- `build()` - Validate and create the simulation

#### Initialization Examples

```python
import centrevo as cv

# 1. Uniform initialization (default)
sim = (cv.SimulationBuilder()
    .population_size(50)
    .generations(100)
    .repeat_structure(171, 12, 10)
    .build())  # All sequences initialized with 'A'

# 2. Uniform with different base
sim = (cv.SimulationBuilder()
    .population_size(50)
    .generations(100)
    .repeat_structure(171, 12, 10)
    .init_uniform(cv.Nucleotide.G())
    .build())  # All sequences initialized with 'G'

# 3. Random initialization
sim = (cv.SimulationBuilder()
    .population_size(50)
    .generations(100)
    .repeat_structure(171, 12, 10)
    .init_random()
    .seed(42)  # Reproducible random sequences
    .build())

# 4. From FASTA file
sim = (cv.SimulationBuilder()
    .population_size(50)
    .generations(100)
    .init_from_fasta("sequences.fasta")
    .mutation_rate(0.0001)
    .build())

# 5. From JSON
sim = (cv.SimulationBuilder()
    .population_size(50)
    .generations(100)
    .init_from_json("sequences.json")
    .build())

# 6. From checkpoint (sequences only, new parameters)
sim = (cv.SimulationBuilder()
    .generations(500)  # Run 500 more generations
    .init_from_checkpoint("sim.db", sim_id="exp1", generation=1000)
    .mutation_rate(0.0002)  # Different mutation rate
    .build())
```

#### Resume from Checkpoint

To resume a simulation exactly where it left off (preserving all state):

```python
import centrevo as cv

# Resume complete simulation state
sim = cv.Simulation.resume_from_checkpoint("simulation.db", sim_id="exp1")
sim.run_for(100)  # Continue for 100 more generations
```

#### Default Values

The builder uses sensible defaults that create a "null model" (no evolution):
- **Mutation rate**: 0.0 (no mutation)
- **Recombination**: All rates 0.0 (no recombination)
- **Fitness**: Neutral (no selection)
- **Alphabet**: DNA (A, C, G, T)
- **Init base**: A
- **Seed**: None (non-deterministic)
- **Chromosomes per haplotype**: 1

This allows you to easily create a baseline and then add complexity:

```python
# Start with null model
sim = (cv.SimulationBuilder()
    .population_size(50)
    .generations(100)
    .repeat_structure(171, 12, 10)
    .build())  # No mutation, no recombination

# Add mutation only
sim = (cv.SimulationBuilder()
    .population_size(50)
    .generations(100)
    .repeat_structure(171, 12, 10)
    .mutation_rate(0.0001)  # Now with mutation
    .build())

# Add recombination too
sim = (cv.SimulationBuilder()
    .population_size(50)
    .generations(100)
    .repeat_structure(171, 12, 10)
    .mutation_rate(0.0001)
    .recombination(0.01, 0.7, 0.1)  # Now with both
    .build())
```

### Fitness Configuration

Centrevo supports various fitness functions for natural selection. By default, fitness is **neutral** (no selection), but you can specify selection pressures using the `FitnessConfig` builder.

#### Available Fitness Functions

**Neutral Fitness (default)**
```python
# Explicitly set neutral fitness (though this is the default)
fitness = cv.FitnessConfig.neutral()
sim = (cv.SimulationBuilder()
    .population_size(100)
    .generations(50)
    .repeat_structure(171, 12, 10)
    .fitness(fitness)  # Optional - neutral is default
    .build())
```

**GC Content Fitness**
```python
# Selection for optimal GC content
fitness = cv.FitnessConfig.with_gc_content(
    optimum=0.5,        # Optimal GC content (50%)
    concentration=2.0   # Selection strength (higher = stronger)
).build()

sim = (cv.SimulationBuilder()
    .population_size(100)
    .generations(50)
    .repeat_structure(171, 12, 10)
    .mutation_rate(0.0001)
    .fitness(fitness)
    .build())
```

**Length-Based Fitness**
```python
# Selection for optimal sequence length
fitness = cv.FitnessConfig.with_length(
    optimum=20000,  # Optimal length in bases
    std_dev=0.5     # Standard deviation (lower = stronger selection)
).build()
```

**Sequence Similarity Fitness**
```python
# Selection favoring similarity between haplotypes
fitness = cv.FitnessConfig.with_similarity(
    shape=2.0  # Shape parameter (controls decline rate)
).build()
```

**Combined Fitness Functions**
```python
# Multiple selection pressures (fitness values are multiplied)
fitness = (cv.FitnessConfig.with_gc_content(0.5, 2.0)
    .with_length(20000, 0.5)
    .with_similarity(2.0)
    .build())

sim = (cv.SimulationBuilder()
    .population_size(100)
    .generations(100)
    .repeat_structure(171, 12, 10)
    .mutation_rate(0.0001)
    .recombination(0.01, 0.7, 0.1)
    .fitness(fitness)  # Combined selection
    .build())
```

#### Complete Example with Selection

```python
import centrevo as cv

# Create fitness configuration
fitness = cv.FitnessConfig.with_gc_content(0.5, 2.0).build()

# Build simulation with selection
sim = (cv.SimulationBuilder()
    .population_size(100)
    .generations(50)
    .repeat_structure(171, 12, 10)
    .mutation_rate(0.0001)
    .recombination(0.01, 0.7, 0.1)
    .fitness(fitness)
    .seed(42)
    .build())

# Track GC content over time
print(f"{'Gen':>4} {'Mean GC':>10}")
for gen in range(0, 51, 10):
    if gen > 0:
        sim.run_for(10)
    population = sim.population()
    gc = cv.gc_content(population, None, None, 0)
    print(f"{gen:4d} {gc:10.4f}")

# With selection, GC content approaches the optimum (0.5)
```

For more examples, see `examples/python_fitness_example.py` and `examples/python_builder_example.py`.

### Simulation Class

Once built, the `Simulation` object provides methods to run and inspect simulations:

- `run()` - Run for the configured number of generations
- `run_for(n)` - Run for N more generations
- `step()` - Advance by one generation
- `population()` - Get current population
- `generation()` - Get current generation number

```python
# Create simulation
sim = cv.SimulationBuilder().population_size(50).generations(100).repeat_structure(171, 12, 10).build()

# Run step by step
for i in range(10):
    sim.step()
    if i % 10 == 0:
        pop = sim.population()
        pi = cv.nucleotide_diversity(pop, 0)
        print(f"Gen {sim.generation()}: π = {pi:.6f}")

# Run for 90 more generations
sim.run_for(90)

# Or run to completion
sim.run()
```

### Core Classes

**Note:** While the core classes can still be used directly (legacy API), the `SimulationBuilder` is now the recommended way to create simulations.

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
    - `record_full_config(structure, mutation, recombination, fitness, config)` - Save complete config
    - `record_generation(population, generation)` - Save generation state
    - `record_checkpoint(simulation, generation)` - Save checkpoint for resume
    - `finalize_metadata()` - Finalize simulation metadata
    - `close()` - Close the recorder
  
- **`QueryBuilder`**: Query saved simulations
  - Constructor: `QueryBuilder(db_path)`
  - Methods:
    - `list_simulations()` - Get all simulation names
    - `get_recorded_generations(sim_id)` - Get recorded generation numbers
    - `get_simulation_info(sim_id)` - Get simulation metadata as dict
    - `get_generation(sim_id, generation)` - Load a specific generation as Population
    - `get_fitness_history(sim_id)` - Get fitness statistics over time
    - `close()` - Close the database connection

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

- **`export_fasta(population)`**: Export sequences to FASTA format
  - Returns: FASTA formatted string
  - Example: `fasta = centrevo.export_fasta(pop)`

- **`export_csv(population)`**: Export sequences to CSV format
  - Returns: CSV formatted string with columns: individual_id, haplotype, chromosome, sequence
  - Example: `csv = centrevo.export_csv(pop)`

- **`export_json(population)`**: Export sequences to JSON format
  - Returns: JSON formatted string
  - Example: `json_str = centrevo.export_json(pop)`

- **`export_metadata_json(info_dict)`**: Export metadata to JSON
  - Args: Dictionary from `QueryBuilder.get_simulation_info()`
  - Returns: JSON formatted metadata
  - Example: `json_str = centrevo.export_metadata_json(info)`

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
import centrevo as cv
import centrevo.plotting as cplt

# 1. Create and run simulation with builder pattern
sim = (cv.SimulationBuilder()
    .population_size(100)
    .generations(1000)
    .repeat_structure(ru_length=171, rus_per_hor=12, hors_per_chr=100)
    .mutation_rate(0.0001)
    .recombination(break_prob=0.01, crossover_prob=0.7, gc_extension_prob=0.1)
    .seed(42)
    .build())

# Run simulation
sim.run()

# Get final population
population = sim.population()

# 2. Analyze population
print("Population Analysis")
print("=" * 50)

# Diversity metrics
pi = cv.nucleotide_diversity(population, 0)
tajima = cv.tajimas_d(population, 0)
theta = cv.wattersons_theta(population, 0)
hap_div = cv.haplotype_diversity(population, 0)

print(f"Nucleotide diversity (π): {pi:.6f}")
print(f"Watterson's θ: {theta:.6f}")
print(f"Tajima's D: {tajima:.4f}")
print(f"Haplotype diversity: {hap_div:.4f}")

# Composition
gc = cv.gc_content(population, None, None, None)
comp = cv.nucleotide_composition(population, None, None, None)

print(f"\nGC content: {gc:.2%}")
print("Nucleotide composition:")
for nuc, freq in sorted(comp.items()):
    print(f"  {nuc}: {freq:.4f}")

# Polymorphism
seg_sites_h1 = cv.count_segregating_sites(population, 0, 0)
seg_sites_h2 = cv.count_segregating_sites(population, 0, 1)

print(f"\nSegregating sites (haplotype 1): {seg_sites_h1}")
print(f"Segregating sites (haplotype 2): {seg_sites_h2}")

# 3. Export and visualize
print("\n" + "=" * 50)
print("Exporting data for visualization...")

# Export diversity metrics (for trajectory plots)
metrics = cv.export_diversity_metrics(population, 0)

# Export LD decay
ld_data = cv.export_ld_decay(population, 0, 0, 1000, 10)

# Export distance matrix
dist_matrix = cv.distance_matrix(population, 0)

# Export sequences
fasta_str = cv.export_fasta(population)
with open("sequences.fasta", "w") as f:
    f.write(fasta_str)

csv_str = cv.export_csv(population)
with open("sequences.csv", "w") as f:
    f.write(csv_str)

# Create plots
fig_ld = cplt.plot_ld_decay(ld_data, save_path="ld_decay.png")
fig_dist = cplt.plot_distance_matrix(dist_matrix, save_path="distances.png")
fig_comp = cplt.plot_nucleotide_composition(comp, save_path="composition.png")

print("Plots saved!")

# 4. PyArrow integration for pandas/polars
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

## Example: Simulation with Recording

```python
import centrevo as cv

# Setup recording strategy
recorder = cv.Recorder(
    "my_simulation.db",
    "experiment1",
    cv.RecordingStrategy.every_n(100)
)

# Create simulation
sim = (cv.SimulationBuilder()
    .population_size(100)
    .generations(1000)
    .repeat_structure(171, 12, 100)
    .mutation_rate(0.0001)
    .seed(42)
    .build())

# Record initial configuration
config = cv.SimulationConfig(100, 1000, 42)
recorder.record_metadata(config)

# Record initial population
recorder.record_generation(sim.population(), 0)

# Run simulation with periodic recording
for i in range(10):
    sim.run_for(100)
    recorder.record_generation(sim.population(), sim.generation())
    print(f"Recorded generation {sim.generation()}")

# Finalize and close
recorder.finalize_metadata()
recorder.close()
print("Simulation complete and saved!")
```

## Example: Loading from Custom Sequences

```python
import centrevo as cv

# Load and evolve sequences from FASTA
sim = (cv.SimulationBuilder()
    .population_size(50)
    .generations(500)
    .init_from_fasta("initial_sequences.fasta")
    .mutation_rate(0.0001)
    .recombination(0.01, 0.7, 0.1)
    .seed(42)
    .build())

sim.run()
population = sim.population()

# Analyze evolved sequences
pi = cv.nucleotide_diversity(population, 0)
print(f"Final nucleotide diversity: {pi:.6f}")

# Export evolved sequences
evolved_fasta = cv.export_fasta(population)
with open("evolved_sequences.fasta", "w") as f:
    f.write(evolved_fasta)
```

## Example: Parameter Sweep Starting from Checkpoint

```python
import centrevo as cv

# Run initial simulation
initial_sim = (cv.SimulationBuilder()
    .population_size(100)
    .generations(1000)
    .repeat_structure(171, 12, 50)
    .mutation_rate(0.0001)
    .seed(42)
    .build())

initial_sim.run()

# Save checkpoint (implement saving logic)
# ... save initial_sim state to database ...

# Now run parameter sweep from same initial state
mutation_rates = [0.00005, 0.0001, 0.0002, 0.0004]
results = []

for rate in mutation_rates:
    # Start from checkpoint but with different mutation rate
    sim = (cv.SimulationBuilder()
        .generations(500)
        .init_from_checkpoint("checkpoint.db", "initial", generation=1000)
        .mutation_rate(rate)
        .build())
    
    sim.run()
    population = sim.population()
    pi = cv.nucleotide_diversity(population, 0)
    
    results.append({
        'mutation_rate': rate,
        'nucleotide_diversity': pi
    })
    print(f"Rate {rate}: π = {pi:.6f}")

# Analyze results
import pandas as pd
df = pd.DataFrame(results)
print(df)
```

## Loading and Querying Saved Simulations

```python
import centrevo

# Open database
query = centrevo.QueryBuilder("my_simulation.db")

# List all simulations
sims = query.list_simulations()
print(f"Simulations: {sims}")

# Get simulation info
info = query.get_simulation_info("experiment1")
print(f"Population size: {info['pop_size']}")
print(f"Generations: {info['num_generations']}")
print(f"Mutation rate: {info['mutation_rate']}")

# Get recorded generations
generations = query.get_recorded_generations("experiment1")
print(f"Recorded generations: {generations}")

# Load a specific generation
population = query.get_generation("experiment1", 1000)
print(f"Loaded population with {population.size()} individuals")

# Analyze loaded population
pi = centrevo.nucleotide_diversity(population, 0)
print(f"Nucleotide diversity at gen 1000: {pi:.6f}")

# Get fitness history
fitness_history = query.get_fitness_history("experiment1")
for record in fitness_history:
    print(f"Gen {record['generation']}: mean={record['mean']:.4f}, std={record['std']:.4f}")

# Export loaded population
fasta = centrevo.export_fasta(population)
with open("gen1000.fasta", "w") as f:
    f.write(fasta)

# Close query builder
query.close()
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

### Pattern 4: Database Workflow

```python
import centrevo

# Query and load from database
query = centrevo.QueryBuilder("simulation.db")

# Get all simulations
sims = query.list_simulations()
print(f"Found {len(sims)} simulations")

# Analyze each simulation
for sim_name in sims:
    # Get metadata
    info = query.get_simulation_info(sim_name)
    print(f"\nAnalyzing: {sim_name}")
    print(f"  Population: {info['pop_size']}")
    print(f"  Generations: {info['num_generations']}")
    
    # Get recorded generations
    gens = query.get_recorded_generations(sim_name)
    
    # Load and analyze final generation
    if gens:
        final_gen = max(gens)
        pop = query.get_generation(sim_name, final_gen)
        
        # Run analysis
        pi = centrevo.nucleotide_diversity(pop, 0)
        print(f"  Final π: {pi:.6f}")
        
        # Export results
        fasta = centrevo.export_fasta(pop)
        with open(f"{sim_name}_final.fasta", "w") as f:
            f.write(fasta)

# Get fitness trajectories
for sim_name in sims:
    history = query.get_fitness_history(sim_name)
    
    # Plot fitness over time
    import matplotlib.pyplot as plt
    generations = [h['generation'] for h in history]
    mean_fitness = [h['mean'] for h in history]
    
    plt.plot(generations, mean_fitness, label=sim_name)

plt.xlabel('Generation')
plt.ylabel('Mean Fitness')
plt.legend()
plt.savefig('fitness_comparison.png')

query.close()
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
docker run --rm -v "$(pwd)":/io ghcr.io/pyo3/maturin build --release --target x86_64-unknown-linux-gnu -i python3.12
```

All wheels will be created in `target/wheels/`.

## Performance Notes

- All core computations are done in Rust for maximum performance
- Analysis functions use parallel computation via Rayon when beneficial
- Database operations use SQLite with WAL mode for concurrent access
- PyArrow provides zero-copy data transfer between Rust and Python
- Large populations (>1000 individuals) analyze in seconds

## Requirements

- Python 3.12+
- Rust toolchain 1.88+ (for building from source)
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
