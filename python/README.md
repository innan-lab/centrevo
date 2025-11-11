# Centrevo Python Package

High-performance population genetics simulator with comprehensive analysis tools.

## Installation

```bash
# Install from source (requires Rust and maturin)
pip install maturin
maturin develop --release

# Or build wheel
maturin build --release
pip install target/wheels/centrevo-*.whl
```

## Quick Start

### Basic Simulation

```python
import centrevo

# Create a population
alphabet = centrevo.Alphabet.dna()
base_a = centrevo.Nucleotide.A()

structure = centrevo.RepeatStructure(
    alphabet=alphabet,
    init_base=base_a,
    ru_length=171,
    rus_per_hor=12,
    hors_per_chr=50,
    chrs_per_hap=1,
)

pop = centrevo.create_initial_population(size=100, structure=structure)
print(f"Population size: {pop.size()}")
```

### Diversity Analysis

```python
# Calculate diversity metrics
pi = centrevo.nucleotide_diversity(pop, chromosome_idx=0)
tajima_d = centrevo.tajimas_d(pop, chromosome_idx=0)
theta_w = centrevo.wattersons_theta(pop, chromosome_idx=0)
hap_div = centrevo.haplotype_diversity(pop, chromosome_idx=0)

print(f"π: {pi:.6f}")
print(f"Tajima's D: {tajima_d:.6f}")
print(f"θ_W: {theta_w:.6f}")
print(f"Haplotype diversity: {hap_div:.6f}")
```

### Linkage Disequilibrium

```python
# Calculate LD between two sites
ld_stats = centrevo.linkage_disequilibrium(
    pop, pos1=100, pos2=500,
    chromosome_idx=0, haplotype_idx=0
)

print(f"D: {ld_stats['D']:.6f}")
print(f"D': {ld_stats['D_prime']:.6f}")
print(f"r²: {ld_stats['r_squared']:.6f}")

# Calculate LD decay
ld_decay_dict = centrevo.ld_decay(
    pop, chromosome_idx=0, haplotype_idx=0,
    max_distance=1000, bin_size=50
)

print(f"Distances: {ld_decay_dict['distances'][:5]}")
print(f"r² values: {ld_decay_dict['r_squared_values'][:5]}")
```

### Distance Analysis

```python
# Calculate pairwise distances
distances = centrevo.pairwise_distances(pop, chromosome_idx=0)
print(f"Number of pairwise distances: {len(distances)}")

# Calculate full distance matrix
matrix = centrevo.distance_matrix(pop, chromosome_idx=0)
print(f"Distance matrix size: {len(matrix)}×{len(matrix[0])}")
```

### Composition Analysis

```python
# Population-level GC content
gc_pop = centrevo.gc_content(pop, None, None, None)
print(f"Population GC content: {gc_pop:.4f}")

# Individual-level GC content
gc_ind = centrevo.gc_content(pop, individual_idx=0, haplotype_idx=None, chromosome_idx=None)
print(f"Individual 0 GC content: {gc_ind:.4f}")

# Nucleotide composition
comp = centrevo.nucleotide_composition(pop, None, None, None)
for nuc, freq in sorted(comp.items()):
    print(f"{nuc}: {freq:.4f}")
```

## Data Export with PyArrow

Export analysis results to PyArrow format for use with pandas or polars:

```python
import centrevo
from centrevo.plotting import export_to_pyarrow_table
import pandas as pd
import polars as pl

# Export diversity metrics
metrics = centrevo.export_diversity_metrics(pop, chromosome_idx=0)

# Convert to PyArrow Table
diversity_data = [metrics]  # List of dicts from multiple generations
table = export_to_pyarrow_table(diversity_data)

# Convert to pandas or polars
df_pandas = table.to_pandas()
df_polars = pl.from_arrow(table)

print(df_pandas)
```

### Export LD Decay

```python
from centrevo.plotting import export_ld_decay_to_pyarrow

# Export LD decay data
ld_data = centrevo.export_ld_decay(
    pop, chromosome_idx=0, haplotype_idx=0,
    max_distance=1000, bin_size=50
)

# Convert to PyArrow Table
table = export_ld_decay_to_pyarrow(ld_data)

# Use with pandas/polars
df = table.to_pandas()
```

### Export Distance Matrix

```python
# Export distance matrix in long format
dist_data = centrevo.export_distance_matrix(pop, chromosome_idx=0)

# Each entry has sequence_i, sequence_j, distance
df = pd.DataFrame(dist_data)
pivot_matrix = df.pivot(index='sequence_i', columns='sequence_j', values='distance')
```

## Visualization

Simple plotting utilities for common analyses:

```python
from centrevo.plotting import (
    plot_diversity_trajectory,
    plot_ld_decay,
    plot_distance_matrix,
    plot_nucleotide_composition,
)

# Plot nucleotide composition
comp = centrevo.nucleotide_composition(pop, None, None, None)
fig = plot_nucleotide_composition(comp, save_path="composition.png")

# Plot distance matrix heatmap
matrix = centrevo.distance_matrix(pop, chromosome_idx=0)
fig = plot_distance_matrix(matrix, save_path="distances.png")

# Plot LD decay
ld_data = centrevo.export_ld_decay(pop, 0, 0, 1000, 50)
fig = plot_ld_decay(ld_data, save_path="ld_decay.png")

# Plot diversity over time (requires data from multiple generations)
# diversity_data = [centrevo.export_diversity_metrics(pop_t, 0) for pop_t in populations]
# fig = plot_diversity_trajectory(diversity_data, metric='nucleotide_diversity')
```

## Advanced Usage

### Custom Visualizations

For complex visualizations, work directly with PyArrow data:

```python
import pyarrow as pa
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Get data
ld_data = centrevo.export_ld_decay(pop, 0, 0, 1000, 50)
df = pd.DataFrame(ld_data)

# Custom plot
sns.set_style("whitegrid")
fig, ax = plt.subplots(figsize=(12, 6))
sns.lineplot(data=df, x='distance', y='r_squared', ax=ax)
ax.set_title("Custom LD Decay Plot")
plt.show()
```

### Batch Analysis

```python
# Analyze multiple populations or generations
results = []

for gen in range(0, 1000, 100):
    pop = load_population(gen)  # Your loading logic
    metrics = centrevo.export_diversity_metrics(pop, 0)
    results.append(metrics)

# Convert to DataFrame
df = pd.DataFrame(results)

# Analyze trends
print(df.describe())
```

## API Reference

### Analysis Functions

- `nucleotide_diversity(pop, chr_idx)` - Calculate π
- `tajimas_d(pop, chr_idx)` - Calculate Tajima's D
- `wattersons_theta(pop, chr_idx)` - Calculate θ_W
- `haplotype_diversity(pop, chr_idx)` - Calculate haplotype diversity
- `linkage_disequilibrium(pop, pos1, pos2, chr_idx, hap_idx)` - Calculate LD
- `ld_decay(pop, chr_idx, hap_idx, max_dist, bin_size)` - Calculate LD decay
- `haplotype_blocks(pop, chr_idx, hap_idx, threshold)` - Identify LD blocks
- `pairwise_distances(pop, chr_idx)` - Calculate all pairwise distances
- `distance_matrix(pop, chr_idx)` - Calculate distance matrix
- `gc_content(pop, ind_idx, hap_idx, chr_idx)` - Calculate GC content (flexible)
- `nucleotide_composition(pop, ind_idx, hap_idx, chr_idx)` - Calculate composition
- `count_segregating_sites(pop, chr_idx, hap_idx)` - Count polymorphic sites

### Export Functions

- `export_diversity_metrics(pop, chr_idx)` - Export all diversity metrics
- `export_distance_matrix(pop, chr_idx)` - Export distance matrix (long format)
- `export_ld_decay(pop, chr_idx, hap_idx, max_dist, bin_size)` - Export LD decay

### Plotting Functions

- `plot_diversity_trajectory(data, metric, ...)` - Plot metric over time
- `plot_ld_decay(data, ...)` - Plot LD decay
- `plot_distance_matrix(matrix, ...)` - Plot distance heatmap
- `plot_nucleotide_composition(comp, ...)` - Plot nucleotide frequencies
- `plot_multiple_diversity_metrics(data, ...)` - Multi-panel diversity plots

### PyArrow Helpers

- `export_to_pyarrow_table(data, schema=None)` - Convert to PyArrow Table
- `export_ld_decay_to_pyarrow(data)` - Convert LD decay to PyArrow
- `export_distance_matrix_to_pyarrow(data)` - Convert distances to PyArrow

## Examples

See `examples/python_analysis_example.py` for comprehensive examples.

## Dependencies

- Python >= 3.8
- pyarrow >= 14.0.0
- matplotlib >= 3.5.0
- numpy >= 1.21.0

Optional:
- pandas >= 1.5.0 (for DataFrame conversion)
- polars >= 0.19.0 (for DataFrame conversion)
- jupyter >= 1.0.0 (for notebooks)

## License

MIT License - see LICENSE file for details.
