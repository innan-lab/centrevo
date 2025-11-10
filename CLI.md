# Centrevo CLI Guide

Complete command-line interface reference for the Centrevo centromeric evolution simulator.

## Installation

Build the CLI from source:

```bash
cd /Users/kent/repos/centrevo
cargo build --release
```

The binary will be available at `./target/release/centrevo`.

## Quick Start

Initialize a simulation:
```bash
./target/release/centrevo init -N my_sim -n 100 -g 1000
```

List simulations:
```bash
./target/release/centrevo list
```

View simulation info:
```bash
./target/release/centrevo info -N my_sim
```

## Commands

### `centrevo init` - Initialize a New Simulation

Create and initialize a new simulation with specified parameters.

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
🧬 Centrevo - Centromeric Evolution Simulator
============================================

Initializing simulation: my_sim

Creating initial population...
✓ Created 100 individuals

Setting up database...
✓ Database created: simulation.db

Simulation initialized successfully!
  Name: my_sim
  Population size: 100
  Generations: 1000
  Structure: 171bp RU × 12 RUs/HOR × 100 HORs

💡 Use 'centrevo list' to view simulations
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
📊 Simulations in simulation.db:
==================================================
  • test_sim
  • alpha_sat
  • my_simulation

💡 Use 'centrevo info --name <name>' for details
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
📊 Simulation Information: my_sim
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
📈 Recorded Generations for 'my_sim':
==================================================
Generations: [0, 100, 200, 300, 400, 500]
Total: 6 snapshots
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

### 1. Initialize and Inspect

```bash
# Create simulation
centrevo init -N test -n 50 -g 100

# List all simulations
centrevo list

# Check details
centrevo info -N test

# View recorded generations
centrevo generations -N test
```

### 2. Multiple Simulations in One Database

```bash
# Create first simulation
centrevo init -N sim1 -n 100 -g 1000 -o shared.db

# Create second simulation (same database)
centrevo init -N sim2 -n 200 -g 500 -o shared.db

# List both
centrevo list -d shared.db
```

### 3. Reproducible Research

```bash
# Use fixed seed for reproducibility
centrevo init -N reproducible -n 100 -g 1000 --seed 12345

# Document in paper:
# "Simulations were run with centrevo v0.1.0 using seed 12345"
```

### 4. Parameter Sweep

```bash
# Different population sizes
for n in 50 100 200 500; do
  centrevo init -N "pop_${n}" -n $n -g 1000 -o sweep.db
done

# Check all simulations
centrevo list -d sweep.db
```

---

## Tips & Best Practices

### Performance

- **Large populations**: Consider reducing recording frequency
- **Long simulations**: Use specific generation recording
- **Multiple runs**: Use separate databases to avoid locking

### Naming Conventions

Use descriptive simulation names:
```bash
# Good names (descriptive)
centrevo init -N "alpha_sat_100ind_5000gen_seed42"
centrevo init -N "test_mutation_rate_0.001"

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

Currently implemented (in code):
- `EveryN(n)` - Record every N generations (default: 100)
- `Specific(vec)` - Record only specified generations
- `All` - Record all generations (storage intensive)
- `None` - No recording (analysis mode)

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
```

### Simulation Not Found

**Problem:** "Failed to get simulation info"

**Solution:**
```bash
# List available simulations
centrevo list

# Use exact name (case-sensitive)
centrevo info -N exact_name
```

### Out of Disk Space

**Problem:** Database grows too large

**Solution:**
```bash
# Check database size
du -h simulation.db

# Reduce recording frequency in future runs
# Or use Specific strategy for important generations only
```

---

## Advanced Usage

### Scripting

Bash script for batch processing:
```bash
#!/bin/bash

# Run multiple simulations
for seed in {1..10}; do
  echo "Running simulation with seed $seed"
  centrevo init \
    -N "batch_seed_${seed}" \
    -n 100 \
    -g 1000 \
    --seed $seed \
    -o batch_results.db
  
  # Check if successful
  if [ $? -eq 0 ]; then
    echo "✓ Seed $seed complete"
  else
    echo "✗ Seed $seed failed"
  fi
done

# Summarize results
echo ""
echo "Results:"
centrevo list -d batch_results.db
```

### Integration with Python

Use CLI from Python:
```python
import subprocess
import json

# Run simulation
subprocess.run([
    "./centrevo", "init",
    "-N", "python_sim",
    "-n", "100",
    "-g", "1000",
    "--seed", "42"
])

# Query results
result = subprocess.run(
    ["./centrevo", "info", "-N", "python_sim"],
    capture_output=True,
    text=True
)
print(result.stdout)
```

---

## Version Information

```bash
centrevo --version
```

Output: `centrevo 0.1.0`

---

## Support & Documentation

- **Implementation Guide**: `IMPLEMENTATION_GUIDE.md`
- **Python Bindings**: `PYTHON_BINDINGS.md`
- **Storage Module**: `src/storage/README.md`
- **API Documentation**: Run `cargo doc --open`

---

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | General error (command failed) |
| 101 | Configuration error (invalid parameters) |

---

Built with ❤️ using Rust and Clap 4.5
