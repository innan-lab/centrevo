# Centrevo Python Package

High-performance population genetics simulator with comprehensive analysis tools.

This README describes the **implementation details** of the Python package. For **usage documentation**, see [PYTHON.md](../PYTHON.md) in the root directory.

## Package Structure

```
python/
├── centrevo/           # Python package
│   ├── __init__.py    # Package initialization and re-exports
│   ├── __init__.pyi   # Type stubs for IDE support
│   ├── plotting.py    # Visualization utilities
│   ├── plotting.pyi   # Type stubs for plotting module
│   └── py.typed       # PEP 561 marker for type checking
├── tests/             # Python test suite
│   ├── test_analysis.py
│   ├── test_core.py
│   ├── test_custom_sequences.py
│   └── test_storage.py
├── pytest.ini         # pytest configuration
└── README.md          # This file
```

## Building the Package

### Prerequisites

- Python >= 3.12
- Rust toolchain >= 1.88 (for compiling the extension)
- maturin (Python build tool)

### Build Process

```bash
# Install maturin
pip install maturin

# Development build (debug)
maturin develop

# Release build (optimized)
maturin develop --release

# Build wheel for distribution
maturin build --release
```

The build process:
1. Compiles Rust code in `src/python/` to a native extension module
2. Packages the extension with Python code in `python/centrevo/`
3. Creates a wheel file in `target/wheels/`

## Implementation Architecture

### Rust Bindings (src/python/)

The Python package is implemented using [PyO3](https://pyo3.rs/), which provides:

- **Zero-copy data exchange**: Using PyArrow for efficient data transfer
- **Automatic type conversion**: Between Rust and Python types
- **Native performance**: All computations done in Rust
- **Memory safety**: Rust's ownership system prevents memory leaks

Key binding modules:
- `src/python/bindings.rs` - Core type bindings (Nucleotide, Sequence, Population, etc.)
- `src/python/analysis.rs` - Analysis function bindings (diversity, LD, distance, etc.)
- `src/python/mod.rs` - Module registration and exports

### Python Layer (python/centrevo/)

#### `__init__.py`

Main package initialization:
- Re-exports all Rust types and functions from the native module
- Provides helper functions (e.g., `create_initial_population`)
- No Python-side computation (pure Rust delegation)

#### `plotting.py`

Visualization utilities using matplotlib:
- `plot_diversity_trajectory()` - Time series plots
- `plot_ld_decay()` - LD decay curves
- `plot_distance_matrix()` - Heatmaps
- `plot_nucleotide_composition()` - Bar charts
- `plot_multiple_diversity_metrics()` - Multi-panel plots

PyArrow export helpers:
- `export_to_pyarrow_table()` - Generic dict → PyArrow Table
- `export_ld_decay_to_pyarrow()` - LD data → PyArrow Table
- `export_distance_matrix_to_pyarrow()` - Distance matrix → PyArrow Table

## Type Stubs (`.pyi` files)

The package includes type stubs for IDE support and static type checking:

- `__init__.pyi` - Complete type signatures for all exported functions and classes
- `plotting.pyi` - Type signatures for plotting functions

These enable:
- Autocomplete in IDEs (VS Code, PyCharm, etc.)
- Type checking with mypy
- Better documentation in Jupyter notebooks

## Testing

### Test Structure

```
tests/
├── test_core.py             # Core functionality tests
├── test_analysis.py         # Analysis function tests
├── test_custom_sequences.py # Custom initialization tests
└── test_storage.py          # Database/persistence tests
```

### Running Tests

```bash
# Install test dependencies
pip install pytest

# Run all tests
pytest

# Run with verbose output
pytest -v

# Run specific test file
pytest tests/test_analysis.py

# Run specific test
pytest tests/test_analysis.py::test_nucleotide_diversity
```

### Test Coverage

Tests cover:
- Core type creation and manipulation
- Analysis functions (diversity, LD, distance, composition)
- Custom sequence initialization (FASTA, JSON, database)
- Checkpoint and resume functionality
- Database operations and queries
- Error handling and edge cases

## PyArrow Integration

The package uses PyArrow for efficient data exchange:

### Why PyArrow?

- **Zero-copy**: No data copying between Rust and Python
- **Columnar format**: Efficient for analytics
- **Interoperability**: Works with pandas, polars, DuckDB, etc.
- **Type safety**: Strong typing in Arrow format

### Implementation

Export functions return Python dicts that can be converted to PyArrow:

```rust
// Rust side (simplified)
#[pyfunction]
fn export_diversity_metrics(population: &Population, chr_idx: usize) -> PyResult<PyObject> {
    // Calculate metrics
    let metrics = calculate_diversity(&population, chr_idx);

    // Return as Python dict
    Python::with_gil(|py| {
        let dict = PyDict::new(py);
        dict.set_item("nucleotide_diversity", metrics.pi)?;
        dict.set_item("tajimas_d", metrics.tajima_d)?;
        // ... more metrics
        Ok(dict.into())
    })
}
```

```python
# Python side
import pyarrow as pa

# Get data
metrics = centrevo.export_diversity_metrics(pop, 0)

# Convert to Arrow
table = pa.Table.from_pylist([metrics])

# Use with pandas/polars
df_pandas = table.to_pandas()
df_polars = pl.from_arrow(table)
```

## Custom Sequence Initialization

### Implementation Details

Located in `src/simulation/initialization.rs`:

#### FASTA Parser
- Uses basic line-by-line parsing
- Validates sequence IDs match expected format (`ind{N}_h{1|2}`)
- Checks sequence lengths against structure
- Reports detailed errors with line numbers

#### JSON Parser
- Expects array of objects with `id` and `seq` fields
- Validates same criteria as FASTA
- More compact for programmatic generation

#### Database Loader
- Queries SQLite database for generation data
- Reconstructs population from stored sequences
- Validates chromosome count and lengths
- Uses latest generation if not specified

### Checkpoint/Resume

Implemented in `src/simulation/engine.rs`:

#### Checkpoint Storage
- Saves full RNG state (rand_chacha::ChaCha8Rng)
- Stores generation number
- References population data via foreign key
- Atomic transaction for consistency

#### Resume Process
1. Load configuration from database
2. Reconstruct population from sequences
3. Restore RNG state from checkpoint
4. Validate generation continuity
5. Continue simulation from exact state

**Key benefit**: Bit-for-bit reproducibility when resuming

## Performance Considerations

### Rust → Python

- All heavy computation happens in Rust
- Python only handles:
  - Function calls and argument passing
  - Plotting with matplotlib
  - Data format conversion (optional)

### Memory Management

- Rust types wrapped in PyO3 smart pointers
- Automatic reference counting
- No manual memory management needed
- GIL released during long computations (where possible)

### Parallelism

- Rayon parallelism in Rust (transparent to Python)
- No Python GIL issues during parallel computation
- Linear scaling with core count for large populations

## Development Workflow

### Making Changes

1. **Edit Rust code** in `src/python/`
2. **Rebuild**: `maturin develop --release`
3. **Test**: `pytest`
4. **Update type stubs** if API changed

### Adding New Functions

Example: Adding a new analysis function

1. **Implement in Rust** (`src/analysis/`)
```rust
pub fn my_new_metric(population: &Population) -> f64 {
    // Implementation
}
```

2. **Create Python binding** (`src/python/analysis.rs`)
```rust
#[pyfunction]
fn my_new_metric(population: &Population) -> PyResult<f64> {
    Ok(crate::analysis::my_new_metric(population))
}
```

3. **Register function** (`src/python/mod.rs`)
```rust
m.add_function(wrap_pyfunction!(my_new_metric, m)?)?;
```

4. **Add type stub** (`python/centrevo/__init__.pyi`)
```python
def my_new_metric(population: Population) -> float: ...
```

5. **Write tests** (`tests/test_analysis.py`)
```python
def test_my_new_metric():
    pop = create_test_population()
    result = centrevo.my_new_metric(pop)
    assert isinstance(result, float)
    assert result >= 0.0
```

6. **Update documentation** (PYTHON.md)

## Dependencies

### Rust Dependencies (from Cargo.toml)

```toml
[dependencies]
pyo3 = { version = "0.20", features = ["extension-module"] }
arrow = "50.0"  # PyArrow interop
# ... other deps
```

### Python Dependencies (runtime)

- `pyarrow >= 14.0.0` - Data interchange
- `matplotlib >= 3.5.0` - Plotting
- `numpy >= 1.21.0` - Numerical operations

### Python Dependencies (development)

- `pytest >= 7.0.0` - Testing framework
- `maturin >= 1.0.0` - Build tool
- `mypy` - Type checking (optional)

## Troubleshooting

### Build Issues

**Error: "pyo3 version mismatch"**
```bash
cargo clean
maturin develop --release
```

**Error: "cannot find -lpython"**
- Ensure Python development headers are installed
- macOS: `xcode-select --install`
- Linux: `apt-get install python3-dev`

### Import Issues

**Error: "No module named 'centrevo'"**
```bash
# Rebuild and install
maturin develop --release
```

**Error: "Symbol not found"**
- Rebuild from clean state
- Check Python version matches build

### Performance Issues

- Ensure using `--release` flag (10-100x faster)
- Check that parallel features are enabled
- Profile with `py-spy` if needed

## Future Enhancements

Planned improvements:
- [ ] Streaming API for large datasets
- [ ] More plotting utilities
- [ ] Jupyter widget for interactive exploration
- [ ] Async API for long-running simulations
- [ ] Better error messages with suggestions
- [ ] Type hints for all plotting functions

## License

MIT License - same as parent Centrevo project.
