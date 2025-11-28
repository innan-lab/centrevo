# Storage Module

This module provides SQLite-based persistence for Centrevo simulations.

## Features

- **SQLite Database**: Embedded database with no server required
- **Recording Strategies**: Flexible recording options (all generations, every N, specific generations)
- **Fitness Tracking**: Aggregated fitness statistics per generation
- **Query Interface**: Rich query API for post-simulation analysis
- **Performance**: Parallel data preparation with transactional writes
- **WAL Mode**: Write-Ahead Logging for concurrent reads during writes

## Quick Start

### Recording a Simulation

```rust
use centrevo::storage::{Recorder, RecordingStrategy};
use centrevo::simulation::{Simulation, SimulationConfig, Population};

// Create recorder
let mut recorder = Recorder::new(
    "my_simulation.db",
    "sim_001",
    RecordingStrategy::EveryN(10) // Record every 10 generations
)?;

// Record simulation metadata
let config = SimulationConfig::new(100, 1000, Some(42));
recorder.record_metadata(&config)?;

// Run simulation and record
for generation in 0..1000 {
    // ... simulation step ...

    recorder.record_if_needed(&population, generation)?;
}

// Finalize
recorder.finalize_metadata()?;
recorder.close()?;
```

### Querying Recorded Data

```rust
use centrevo::storage::QueryBuilder;

// Open database for analysis
let query = QueryBuilder::new("my_simulation.db")?;

// List all simulations
let sims = query.list_simulations()?;

// Get simulation info
let info = query.get_simulation_info("sim_001")?;
println!("Population size: {}", info.pop_size);

// Get individuals at generation 100
let individuals = query.get_generation("sim_001", 100)?;

// Get fitness trajectory
let history = query.get_fitness_history("sim_001")?;
for (gen, stats) in history {
    println!("Gen {}: mean={:.4}, std={:.4}", gen, stats.mean, stats.std);
}

// Find high-fitness individuals
let top_individuals = query.get_high_fitness_individuals("sim_001", 500, 0.9)?;

query.close()?;
```

## Recording Strategies

### EveryN(n)
Records every N generations. Good for long simulations where you want regular checkpoints.

```rust
RecordingStrategy::EveryN(10)  // Record at gen 0, 10, 20, 30, ...
```

### Specific(vec)
Records at specific generations. Useful when you know key timepoints.

```rust
RecordingStrategy::Specific(vec![0, 10, 50, 100, 500, 1000])
```

### All
Records all generations. Use for short simulations or when you need complete history.

```rust
RecordingStrategy::All
```

### None
No recording. Useful for dry runs or parameter tuning.

```rust
RecordingStrategy::None
```

## Database Schema

### Tables

#### `simulations`
Stores simulation metadata:
- `sim_id`: Unique simulation identifier
- `start_time`, `end_time`: Unix timestamps
- `pop_size`: Population size
- `num_generations`: Total generations
- `mutation_rate`, `recombination_rate`: Key parameters
- `parameters_json`: Full configuration as JSON

#### `population_state`
Stores individual-level data:
- `sim_id`, `generation`: Composite key
- `individual_id`: Individual identifier
- `haplotype1_chr_id`, `haplotype2_chr_id`: Chromosome IDs
- `haplotype1_seq`, `haplotype2_seq`: Sequence data (BLOB)
- `fitness`: Cached fitness value
- `timestamp`: Record creation time

#### `fitness_history`
Stores aggregated fitness statistics:
- `sim_id`, `generation`: Composite key (primary)
- `mean_fitness`, `min_fitness`, `max_fitness`, `std_fitness`: Statistics

### Indices
- `idx_pop_sim_gen`: Fast lookup by (sim_id, generation)
- `idx_pop_individual`: Fast lookup by individual_id
- `idx_fitness_sim_gen`: Fast lookup for fitness history

## Performance Considerations

### Memory Usage
- **Sequences stored as BLOBs**: Efficient binary storage
- **Batch inserts**: One transaction per generation reduces overhead
- **Parallel preparation**: Individual snapshots created in parallel

### Disk Usage
For a typical simulation:
- 100 individuals × 10KB chromosome × 2 haplotypes = ~2MB per generation
- Recording every 10 generations for 1000 gens = ~200MB
- Add ~10% for metadata and indices

### Optimization Tips
1. Use `RecordingStrategy::EveryN(n)` for long simulations
2. Record only first chromosome if tracking single locus
3. Compress sequences if > 100KB per individual
4. Run `VACUUM` periodically to reclaim space

## Example Workflows

### Baseline Run
```rust
// Record every 100 generations
let strategy = RecordingStrategy::EveryN(100);
let mut recorder = Recorder::new("baseline.db", "baseline_001", strategy)?;
// ... run simulation ...
```

### Detailed Analysis
```rust
// Record all generations
let strategy = RecordingStrategy::All;
let mut recorder = Recorder::new("detailed.db", "detailed_001", strategy)?;
// ... run simulation ...
```

### Critical Timepoints
```rust
// Record at specific generations of interest
let strategy = RecordingStrategy::Specific(vec![
    0,    // Initial
    100,  // Early
    1000, // Established
    5000, // Late
]);
let mut recorder = Recorder::new("critical.db", "critical_001", strategy)?;
// ... run simulation ...
```

## Error Handling

All operations return `Result<T, DatabaseError>`:

```rust
match recorder.record_generation(&pop, gen) {
    Ok(_) => println!("Recorded generation {}", gen),
    Err(e) => eprintln!("Failed to record: {}", e),
}
```

## Integration with Simulation Engine

```rust
use centrevo::simulation::Simulation;
use centrevo::storage::Recorder;

// Simulation with recording
let mut sim = Simulation::new(config)?;
let mut recorder = Recorder::new("sim.db", "sim_001", RecordingStrategy::EveryN(10))?;

recorder.record_metadata(&sim_config)?;

for gen in 0..config.total_generations {
    sim.step()?;
    recorder.record_if_needed(sim.population(), gen)?;
}

recorder.finalize_metadata()?;
recorder.close()?;
```

## Testing

The module includes comprehensive unit tests:

```bash
cargo test --lib storage
```

## Future Enhancements

Planned improvements:
- [ ] Compression for large sequences (gzip/zstd)
- [ ] Checkpoint/resume functionality
- [ ] Per-chromosome storage option
- [ ] Export to CSV/JSON for external analysis
- [ ] Real-time monitoring queries during simulation
- [ ] Lineage tracking across generations

## License

Same as the parent Centrevo project.
