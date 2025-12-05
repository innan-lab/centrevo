# Centrevo CLI

The command-line interface for the Centrevo simulation engine.

## Usage

```bash
centrevo <COMMAND> [OPTIONS]
```

### Common Commands

- `init`: Initialize a new simulation
- `run`: Run a simulation
- `list`: List simulations in database
- `info`: Show simulation info
- `export`: Export simulation data
- `analyze`: Analyze data

## Example

Here is a quick start guide to running your first simulation.

### 1. Initialize

Create a new simulation named `test_sim`. This creates a database (`simulation.db` by default) and initializes Generation 0.

```bash
centrevo init -N test_sim --population-size 100 --generations 100
```

### 2. Run

Run the simulation. By default, this evolves the population for the configured number of generations.

```bash
centrevo run -N test_sim
```

To see progress or use a specific database:
```bash
centrevo run -N test_sim --database simulation.db
```

### 3. Verify & Export

Check the status of the simulation:

```bash
centrevo info -N test_sim
```

Export sequences from the final generation (e.g., generation 100) to FASTA:

```bash
centrevo export -N test_sim --generation 100 --format fasta --output final_pop.fasta
```

## Development

To run via Cargo:

```bash
cargo run --bin centrevo -- <COMMAND>
```
