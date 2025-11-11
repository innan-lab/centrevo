# Centrevo Python Tests

Comprehensive test suite for the Centrevo Python bindings.

## Test Structure

- **test_core.py**: Tests for core classes and functionality
  - Nucleotide, Alphabet, Chromosome, Haplotype
  - Individual, Population, RepeatStructure
  - SimulationConfig, RecordingStrategy
  - Population creation

- **test_analysis.py**: Tests for analysis functions
  - Diversity metrics (π, Tajima's D, Watterson's θ, haplotype diversity)
  - Linkage disequilibrium (LD, LD decay, haplotype blocks)
  - Distance calculations (pairwise, matrix)
  - Composition analysis (GC content, nucleotide composition)
  - Polymorphism analysis (segregating sites)
  - PyArrow export functions

- **test_storage.py**: Tests for database storage and querying
  - Recorder class (metadata, generation recording)
  - QueryBuilder class (listing simulations, retrieving generations)
  - Recording strategies (EveryN, Specific, All, None)
  - Workflow integration
  - Error handling

## Running Tests

### Prerequisites

1. Build the Python package:
   ```bash
   cd /path/to/centrevo
   maturin develop --release
   ```

2. Install test dependencies:
   ```bash
   pip install pytest pytest-cov
   ```

### Run All Tests

```bash
cd python
pytest
```

### Run Specific Test Files

```bash
pytest tests/test_core.py
pytest tests/test_analysis.py
pytest tests/test_storage.py
```

### Run with Coverage

```bash
pytest --cov=centrevo --cov-report=html
```

### Run with Verbose Output

```bash
pytest -v
```

### Run Specific Test Classes or Functions

```bash
# Run specific class
pytest tests/test_core.py::TestNucleotide

# Run specific test
pytest tests/test_core.py::TestNucleotide::test_create_nucleotides
```

## Test Coverage

The test suite covers:

- ✅ All core classes and their methods
- ✅ Population creation and manipulation
- ✅ All diversity metric calculations
- ✅ Linkage disequilibrium analysis
- ✅ Distance matrix calculations
- ✅ Composition and polymorphism analysis
- ✅ Database recording and querying
- ✅ All recording strategies
- ✅ PyArrow export functions
- ✅ Error handling and edge cases

## Notes

- Tests use fixtures for common setup (populations, temporary databases)
- Linting warnings about "redefining from outer scope" are expected for pytest fixtures
- Some tests may take longer for larger populations (marked with @pytest.mark.slow)
- Temporary database files are automatically cleaned up after tests

## Adding New Tests

When adding new functionality to the Python bindings:

1. Add corresponding tests in the appropriate test file
2. Use existing fixtures when possible
3. Add new fixtures to conftest.py if needed
4. Ensure tests are independent and can run in any order
5. Add appropriate markers (@pytest.mark.slow, @pytest.mark.integration)

## Continuous Integration

These tests should be run as part of CI/CD pipeline to ensure:

- Python bindings compile correctly
- All exposed functions work as expected
- No regressions in API compatibility
- Type hints match actual behavior
