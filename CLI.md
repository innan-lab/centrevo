# Centrevo CLI Guide

## Complete command-line interface reference for the Centrevo centromeric evolution simulator.

**Note:** For programmatic control and advanced workflows, see the [Python Bindings Guide](PYTHON.md), which provides additional features including:
- Direct simulation control (`step()`, `run_for()` methods)
- Custom sequence initialization from FASTA/JSON/database
- Programmatic data loading and querying
- Integration with PyArrow, pandas, and polars
- Flexible export to multiple formats

## Installation

Build the CLI from source:

```bash
cd /Users/kent/repos/centrevo
cargo build --release
```

The binary will be available at `./target/release/centrevo`.

## Quick Start

Use the interactive setup wizard:
```bash
./target/release/centrevo setup
```

Or set up and run manually:
```bash
# Initialize a simulation
./target/release/centrevo init -N my_sim -n 100 -g 1000

# Run the simulation
./target/release/centrevo run -N my_sim

# Analyze results
./target/release/centrevo analyze -N my_sim -g 1000

# Export data
./target/release/centrevo export -N my_sim -g 1000 --format fasta
```

## Commands

### `centrevo setup` - Interactive Setup Wizard

Launch an interactive wizard that guides you through simulation configuration and execution.

**Usage:**
```bash
centrevo setup [OPTIONS]
```

**Options:**

| Flag | Long | Type | Default | Description |
|------|------|------|---------|-------------|
| | `--defaults` | bool | `false` | Skip prompts and use default values |

**Examples:**

Interactive setup:
```bash
centrevo setup
```

Non-interactive with defaults:
```bash
centrevo setup --defaults
```

**Features:**
- Guided configuration with helpful prompts
- Parameter validation
- Configuration summary before execution
- Option to initialize and run in one step
- Calculates and displays genome size

---

### `centrevo init` - Initialize a New Simulation

Create and initialize a new simulation with specified parameters. This sets up the database and creates the initial population but does not run the simulation.

**Usage:**
```bash
centrevo init [OPTIONS]
```

**Usage:**
```bash
centrevo init [OPTIONS]
```

**Options:**

| Flag | Long | Type | Default | Description |
|------|------|------|---------|-------------|
| `-N` | `--name` | String | `"simulation"` | Simulation name/identifier |
| `-o` | `--output` | Path | `"simulation.db"` | Output database file path |
| `-n` | `--population-size` | usize | `100` | Number of diploid individuals |
| `-g` | `--generations` | usize | `1000` | Total generations to simulate |
| | `--ru-length` | usize | `171` | Repeat unit length (base pairs) |
| | `--rus-per-hor` | usize | `12` | Repeat units per higher-order repeat |
| | `--hors-per-chr` | usize | `100` | Higher-order repeats per chromosome |
| | `--seed` | u64 | None | Random seed for reproducibility |

**Examples:**

Basic initialization:
```bash
centrevo init -N test_sim -n 50 -g 500
```

Custom repeat structure:
```bash
centrevo init \
  -N alpha_sat \
  -n 200 \
  --ru-length 171 \
  --rus-per-hor 12 \
  --hors-per-chr 100 \
  -g 5000
```

With reproducible seed:
```bash
centrevo init -N reproducible -n 100 -g 1000 --seed 42
```

Custom output location:
```bash
centrevo init -N my_sim -o /path/to/custom.db -n 100
```

**Output:**
```
Centrevo - Centromeric Evolution Simulator
============================================

Initializing simulation: my_sim

Creating initial population...
Created 100 individuals

Setting up database...
Database created: simulation.db

Simulation initialized successfully!
  Name: my_sim
  Population size: 100
  Generations: 1000
  Structure: 171bp RU × 12 RUs/HOR × 100 HORs

Use 'centrevo run -N my_sim' to run the simulation
```

---

### `centrevo run` - Run a Simulation

Execute a simulation that was initialized with `centrevo init`. This runs the evolutionary process for the specified number of generations.

**Usage:**
```bash
centrevo run [OPTIONS] -N <NAME>
```

**Options:**

| Flag | Long | Type | Default | Description |
|------|------|------|---------|-------------|
| `-d` | `--database` | Path | `"simulation.db"` | Database file to use |
| `-N` | `--name` | String | **Required** | Simulation name |
| | `--mutation-rate` | f64 | `0.001` | Per-site mutation rate |
| | `--recomb-rate` | f64 | `0.01` | Recombination break probability |
| | `--crossover-prob` | f64 | `0.7` | Crossover probability (given break) |
| | `--record-every` | usize | `100` | Record every N generations |
| | `--progress` | bool | `true` | Show progress bar |

**Examples:**

Run with default parameters:
```bash
centrevo run -N my_sim
```

Run with custom evolution parameters:
```bash
centrevo run -N my_sim \
  --mutation-rate 0.002 \
  --recomb-rate 0.05 \
  --crossover-prob 0.8
```

Run with frequent recording:
```bash
centrevo run -N my_sim --record-every 10
```

**Output:**
```
Centrevo - Running Simulation
============================================

Simulation: my_sim
Population size: 100
Target generations: 1000
Mutation rate: 0.001
Recombination rate: 0.01

Starting simulation from generation 0...

Running 1000 generations...
Progress: 100.0% (1000/1000)

Simulation complete!
  Final generation: 1000
  Recorded generations: 11 snapshots

Use 'centrevo info -N my_sim' to view results
```

---

### `centrevo analyze` - Analyze Simulation Data

Calculate population genetics metrics for a specific generation. Computes diversity statistics, composition, and polymorphism measures.

**Usage:**
```bash
centrevo analyze [OPTIONS] -N <NAME> -g <GENERATION>
```

**Options:**

| Flag | Long | Type | Default | Description |
|------|------|------|---------|-------------|
| `-d` | `--database` | Path | `"simulation.db"` | Database file to query |
| `-N` | `--name` | String | **Required** | Simulation name |
| `-g` | `--generation` | usize | **Required** | Generation to analyze |
| | `--chromosome` | usize | `0` | Chromosome index to analyze |
| `-f` | `--format` | String | `"pretty"` | Output format (pretty, json) |
| `-o` | `--output` | Path | None | Output file (stdout if not specified) |

**Examples:**

Analyze generation 1000:
```bash
centrevo analyze -N my_sim -g 1000
```

Export analysis as JSON:
```bash
centrevo analyze -N my_sim -g 1000 --format json -o analysis.json
```

**Output:**
```

Population Genetics Summary
================================
Population size: 100
Sequences analyzed: 200 (2n)
Sequence length: 205200 bp

Diversity Metrics:
------------------
Nucleotide diversity (π): 0.003245
Watterson's theta (θ_W): 0.003189
Tajima's D: 0.2314
Haplotype diversity: 0.9856

Polymorphism:
-------------
Segregating sites: 1245

Composition:
------------
Mean GC content: 48.32%
```

---

### `centrevo export` - Export Simulation Data

Export sequences, metadata, or fitness data in various formats (CSV, JSON, FASTA).

**Usage:**
```bash
centrevo export [OPTIONS] -N <NAME> -g <GENERATION>
```

**Options:**

| Flag | Long | Type | Default | Description |
|------|------|------|---------|-------------|
| `-d` | `--database` | Path | `"simulation.db"` | Database file to query |
| `-N` | `--name` | String | **Required** | Simulation name |
| `-g` | `--generation` | usize | **Required** | Generation to export |
| `-f` | `--format` | String | `"csv"` | Output format (csv, json, fasta) |
| `-o` | `--output` | Path | None | Output file (stdout if not specified) |
| | `--data-type` | String | `"sequences"` | What to export (sequences, metadata, fitness) |

**Examples:**

Export sequences as FASTA:
```bash
centrevo export -N my_sim -g 1000 --format fasta -o sequences.fasta
```

Export sequences as CSV:
```bash
centrevo export -N my_sim -g 1000 --format csv -o sequences.csv
```

Export metadata as JSON:
```bash
centrevo export -N my_sim -g 0 --data-type metadata --format json
```

Export fitness history:
```bash
centrevo export -N my_sim -g 0 --data-type fitness --format csv -o fitness.csv
```

**Output (FASTA format):**
```
>ind0|h1|chr0
ACGTACGTACGTACGT...
>ind0|h2|chr0
ACGTACGTACGTACGT...
>ind1|h1|chr0
ACGTACGTACGTACGT...
...
```

**Output (CSV format):**
```
individual_id,haplotype,chromosome,sequence
ind0,h1,0,ACGTACGTACGTACGT...
ind0,h2,0,ACGTACGTACGTACGT...
ind1,h1,0,ACGTACGTACGTACGT...
...
```

---

### `centrevo validate` - Validate Database Integrity

Check database integrity and verify that simulation data is complete and consistent.

**Usage:**
```bash
centrevo validate [OPTIONS]
```

**Options:**

| Flag | Long | Type | Default | Description |
|------|------|------|---------|-------------|
| `-d` | `--database` | Path | `"simulation.db"` | Database file to validate |
| `-N` | `--name` | String | None | Simulation name (all if not specified) |
| | `--fix` | bool | `false` | Attempt to fix issues automatically |

**Examples:**

Validate entire database:
```bash
centrevo validate
```

Validate specific simulation:
```bash
centrevo validate -N my_sim
```

Validate and attempt fixes:
```bash
centrevo validate --fix
```

**Output:**
```
Validating database: simulation.db

Validating simulation: my_sim
--------------------------------------------------
Metadata: OK
Recorded generations: 11 snapshots
Generation continuity: OK
Population data: OK (100 individuals)
Fitness history: 11 entries

==================================================
Validation complete: No issues found
```

---

### `centrevo list` - List Simulations

Show all simulations stored in a database.

**Usage:**
```bash
centrevo list [OPTIONS]
```

**Options:**

| Flag | Long | Type | Default | Description |
|------|------|------|---------|-------------|
| `-d` | `--database` | Path | `"simulation.db"` | Database file to query |

**Examples:**

List simulations in default database:
```bash
centrevo list
```

List from specific database:
```bash
centrevo list --database /path/to/results.db
```

**Output:**
```
Simulations in simulation.db:
==================================================
  • test_sim
  • alpha_sat
  • my_simulation

Use 'centrevo info --name <name>' for details
```

---

### `centrevo info` - Show Simulation Information

Display detailed information about a specific simulation.

**Usage:**
```bash
centrevo info [OPTIONS] -N <NAME>
```

**Options:**

| Flag | Long | Type | Default | Description |
|------|------|------|---------|-------------|
| `-d` | `--database` | Path | `"simulation.db"` | Database file to query |
| `-N` | `--name` | String | **Required** | Simulation name |

**Examples:**

Show info for simulation:
```bash
centrevo info -N my_sim
```

From specific database:
```bash
centrevo info --database results.db -N alpha_sat
```

**Output:**
```
Simulation Information: my_sim
==================================================
Created: 1762747335
Population size: 100
Generations: 1000

Parameters:
{"population_size":100,"total_generations":1000,"seed":42}
```

---

### `centrevo generations` - List Recorded Generations

Show which generations have been recorded for a simulation.

**Usage:**
```bash
centrevo generations [OPTIONS] -N <NAME>
```

**Options:**

| Flag | Long | Type | Default | Description |
|------|------|------|---------|-------------|
| `-d` | `--database` | Path | `"simulation.db"` | Database file to query |
| `-N` | `--name` | String | **Required** | Simulation name |

**Examples:**

List recorded generations:
```bash
centrevo generations -N my_sim
```

From specific database:
```bash
centrevo generations --database results.db -N test_sim
```

**Output:**
```
Recorded Generations for 'my_sim':
==================================================
Generations: [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
Total: 11 snapshots
```

---

## Database Structure

Centrevo uses SQLite databases with the following schema:

### Tables

1. **`simulations`** - Simulation metadata
   - `sim_id` (TEXT PRIMARY KEY)
   - `start_time` (INTEGER)
   - `end_time` (INTEGER, nullable)
   - `pop_size` (INTEGER)
   - `num_generations` (INTEGER)
   - `mutation_rate` (REAL)
   - `recombination_rate` (REAL)
   - `parameters_json` (TEXT)

2. **`population_state`** - Individual snapshots
   - `sim_id` (TEXT)
   - `generation` (INTEGER)
   - `individual_id` (TEXT)
   - `haplotype1_chr_id` (TEXT)
   - `haplotype1_seq` (BLOB)
   - `haplotype2_chr_id` (TEXT)
   - `haplotype2_seq` (BLOB)
   - `fitness` (REAL)

3. **`fitness_history`** - Aggregated statistics
   - `sim_id` (TEXT)
   - `generation` (INTEGER)
   - `mean_fitness` (REAL)
   - `min_fitness` (REAL)
   - `max_fitness` (REAL)
   - `std_fitness` (REAL)

### Indices

- `idx_sim_gen` on `population_state(sim_id, generation)`
- `idx_fitness_sim` on `fitness_history(sim_id)`
- `idx_fitness_sim_gen` on `fitness_history(sim_id, generation)`

---

## Common Workflows

### 1. Complete Workflow: Initialize, Run, Analyze, Export

```bash
# Initialize simulation
centrevo init -N experiment1 -n 200 -g 5000 --seed 42

# Run simulation
centrevo run -N experiment1 --mutation-rate 0.001

# Analyze final generation
centrevo analyze -N experiment1 -g 5000 --format json -o analysis.json

# Export sequences
centrevo export -N experiment1 -g 5000 --format fasta -o final_pop.fasta

# Export fitness trajectory
centrevo export -N experiment1 -g 0 --data-type fitness -o fitness.csv

# Validate data
centrevo validate -N experiment1
```

### 2. Interactive Setup

```bash
# Launch wizard for guided setup
centrevo setup

# Or use defaults for batch processing
centrevo setup --defaults
```

### 3. Multiple Simulations in One Database

```bash
# Create first simulation
centrevo init -N sim1 -n 100 -g 1000 -o shared.db
centrevo run -N sim1 -d shared.db

# Create second simulation (same database)
centrevo init -N sim2 -n 200 -g 500 -o shared.db
centrevo run -N sim2 -d shared.db

# List both
centrevo list -d shared.db

# Analyze both
centrevo analyze -N sim1 -d shared.db -g 1000
centrevo analyze -N sim2 -d shared.db -g 500
```

### 4. Reproducible Research

```bash
# Use fixed seed for reproducibility
centrevo init -N reproducible -n 100 -g 1000 --seed 12345
centrevo run -N reproducible

# Document in paper:
# "Simulations were run with centrevo v0.2.0 using seed 12345"
```

### 5. Parameter Sweep

```bash
# Different mutation rates
for rate in 0.0001 0.001 0.01; do
  centrevo init -N "mut_${rate}" -n 100 -g 1000 -o sweep.db
  centrevo run -N "mut_${rate}" -d sweep.db --mutation-rate $rate
  centrevo analyze -N "mut_${rate}" -d sweep.db -g 1000 --format json -o "analysis_${rate}.json"
done

# Compare results
centrevo list -d sweep.db
```

---

## Tips & Best Practices

### Performance

- **Large populations**: Recording is the bottleneck - use `--record-every` to reduce frequency
- **Long simulations**: Consider using specific generation recording via the API
- **Multiple runs**: Each simulation in the same database is isolated

### Analysis Workflow

Export data for downstream analysis:
```bash
# Export sequences for external tools
centrevo export -N my_sim -g 1000 --format fasta -o seqs.fasta

# Export population genetics metrics
centrevo analyze -N my_sim -g 1000 --format json -o metrics.json

# Use with Python for visualization
python -c "
import json
import matplotlib.pyplot as plt

with open('metrics.json') as f:
    data = json.load(f)
    print(f'π = {data[\"diversity\"][\"pi\"]:.6f}')
    print(f'Tajima D = {data[\"diversity\"][\"tajima_d\"]:.4f}')
"
```

### Naming Conventions

Use descriptive simulation names:
```bash
# Good names (descriptive)
centrevo init -N "alpha_sat_n100_g5000_mu0.001_seed42"
centrevo init -N "high_recomb_crossover_0.9"

# Bad names (unclear)
centrevo init -N "sim1"
centrevo init -N "test"
```

### Data Management

Organize databases by experiment:
```bash
# Project structure
experiments/
  ├── baseline/
  │   └── simulation.db
  ├── high_mutation/
  │   └── simulation.db
  └── selection/
      └── simulation.db
```

### Recording Strategies

Configure recording frequency based on needs:
```bash
# Frequent recording (every 10 generations) - larger database
centrevo run -N my_sim --record-every 10

# Sparse recording (every 500 generations) - smaller database
centrevo run -N my_sim --record-every 500

# Default (every 100 generations) - balanced
centrevo run -N my_sim
```

For custom recording patterns, use the Python API with `RecordingStrategy.specific([0, 100, 500, 1000])`.

---

## Advanced Usage

### Scripting

Bash script for batch processing:
```bash
#!/bin/bash

# Run multiple simulations with different parameters
for seed in {1..10}; do
  for mut_rate in 0.0001 0.001 0.01; do
    sim_name="batch_seed${seed}_mut${mut_rate}"
    echo "Running: $sim_name"

    centrevo init \
      -N "$sim_name" \
      -n 100 \
      -g 1000 \
      --seed $seed \
      -o batch_results.db

    centrevo run \
      -N "$sim_name" \
      -d batch_results.db \
      --mutation-rate $mut_rate \
      --progress false

    centrevo analyze \
      -N "$sim_name" \
      -d batch_results.db \
      -g 1000 \
      --format json \
      -o "results/${sim_name}_analysis.json"

    if [ $? -eq 0 ]; then
      echo "✓ $sim_name complete"
    else
      echo "✗ $sim_name failed"
    fi
  done
done

# Validate all simulations
centrevo validate -d batch_results.db

# Summarize results
echo ""
echo "Results:"
centrevo list -d batch_results.db
```

### Integration with Python

Use CLI from Python for automation:
```python
import subprocess
import json
import sys

def run_simulation(name, pop_size, generations, seed):
    """Initialize and run a simulation."""
    # Initialize
    init_result = subprocess.run([
        "./centrevo", "init",
        "-N", name,
        "-n", str(pop_size),
        "-g", str(generations),
        "--seed", str(seed)
    ], capture_output=True, text=True)

    if init_result.returncode != 0:
        print(f"Init failed: {init_result.stderr}", file=sys.stderr)
        return False

    # Run
    run_result = subprocess.run([
        "./centrevo", "run",
        "-N", name,
        "--progress", "false"
    ], capture_output=True, text=True)

    if run_result.returncode != 0:
        print(f"Run failed: {run_result.stderr}", file=sys.stderr)
        return False

    return True

def analyze_simulation(name, generation):
    """Analyze simulation and return metrics."""
    result = subprocess.run([
        "./centrevo", "analyze",
        "-N", name,
        "-g", str(generation),
        "--format", "json"
    ], capture_output=True, text=True)

    if result.returncode != 0:
        print(f"Analysis failed: {result.stderr}", file=sys.stderr)
        return None

    return json.loads(result.stdout)

# Run pipeline
if run_simulation("python_controlled", 100, 1000, 42):
    metrics = analyze_simulation("python_controlled", 1000)
    if metrics:
        print(f"Nucleotide diversity: {metrics['diversity']['pi']:.6f}")
        print(f"Tajima's D: {metrics['diversity']['tajima_d']:.4f}")
```

---

## Troubleshooting

### Database Locked

**Problem:** Error opening database (already in use)

**Solution:**
```bash
# Check for processes using the database
lsof simulation.db

# Kill if necessary
kill <PID>

# Or use a different database
centrevo init -N my_sim -o /tmp/simulation_$(date +%s).db
```

### Simulation Not Found

**Problem:** "Failed to get simulation info"

**Solution:**
```bash
# List available simulations (case-sensitive)
centrevo list

# Use exact name
centrevo info -N exact_name
```

### Out of Disk Space

**Problem:** Database grows too large

**Solution:**
```bash
# Check database size
du -h simulation.db

# Reduce recording frequency
centrevo run -N new_sim --record-every 1000

# Or export and delete old data
centrevo export -N old_sim -g 1000 --format fasta -o archive.fasta
# Then delete simulation from database (manual SQL)
```

### Analysis Returns Unexpected Values

**Problem:** Metrics seem incorrect

**Solution:**
```bash
# Validate database first
centrevo validate -N my_sim

# Check that generation was recorded
centrevo generations -N my_sim

# Verify population size matches expectation
centrevo info -N my_sim

# Export and manually inspect sequences
centrevo export -N my_sim -g 1000 --format fasta | head -n 20
```

---

## Version Information

```bash
centrevo --version
```

Output: `centrevo 0.2.0`

---

## Support & Documentation

- **CLI Guide**: `CLI.md` (this file)
- **Python Bindings**: `PYTHON.md`
- **Storage Module**: `src/storage/README.md`
- **Changelog**: `CHANGELOG.md`
- **API Documentation**: Run `cargo doc --open`

---

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | General error (command failed) |
| 101 | Configuration error (invalid parameters) |

---

Built with Rust and Clap 4.5
