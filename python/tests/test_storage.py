"""
Test suite for Centrevo Python bindings - Storage functionality.

Tests database recording and querying capabilities.
"""

import pytest
import centrevo
import tempfile
import os
from pathlib import Path


@pytest.fixture
def temp_db():
    """Create a temporary database file."""
    fd, path = tempfile.mkstemp(suffix='.db')
    os.close(fd)
    yield path
    # Cleanup
    if os.path.exists(path):
        os.remove(path)


@pytest.fixture
def test_population():
    """Create a test population."""
    alphabet = centrevo.Alphabet.dna()
    base = centrevo.Nucleotide.A()

    structure = centrevo.RepeatStructure(
        alphabet=alphabet,
        init_base=base,
        ru_length=100,
        rus_per_hor=5,
        hors_per_chr=10,
        chrs_per_hap=1
    )

    return centrevo.create_initial_population(size=20, structure=structure)


class TestRecorder:
    """Tests for Recorder class."""

    def test_create_recorder(self, temp_db):
        """Test creating a recorder."""
        strategy = centrevo.RecordingStrategy.every_n(10)
        recorder = centrevo.Recorder(temp_db, "test_sim", strategy)

        assert "Recorder" in repr(recorder)
        assert os.path.exists(temp_db)

    def test_record_metadata(self, temp_db, test_population):
        """Test recording simulation metadata."""
        strategy = centrevo.RecordingStrategy.all()
        recorder = centrevo.Recorder(temp_db, "test_sim", strategy)

        config = centrevo.SimulationConfig(
            population_size=20,
            total_generations=100,
            seed=42
        )

        # Should not raise an error
        recorder.record_metadata(config)

    def test_record_generation(self, temp_db, test_population):
        """Test recording a generation."""
        strategy = centrevo.RecordingStrategy.all()
        recorder = centrevo.Recorder(temp_db, "test_sim", strategy)

        config = centrevo.SimulationConfig(
            population_size=test_population.size(),
            total_generations=100,
            seed=1
        )

        recorder.record_metadata(config)
        recorder.record_generation(test_population, generation=0)

        # Verify database was written
        assert os.path.getsize(temp_db) > 0

    def test_record_multiple_generations(self, temp_db, test_population):
        """Test recording multiple generations."""
        strategy = centrevo.RecordingStrategy.specific([0, 10, 20])
        recorder = centrevo.Recorder(temp_db, "test_sim", strategy)

        config = centrevo.SimulationConfig(
            population_size=test_population.size(),
            total_generations=100,
            seed=1
        )

        recorder.record_metadata(config)
        recorder.record_generation(test_population, generation=0)
        recorder.record_generation(test_population, generation=10)
        recorder.record_generation(test_population, generation=20)


class TestQueryBuilder:
    """Tests for QueryBuilder class."""

    def test_create_query_builder(self, temp_db):
        """Test creating a query builder."""
        # Create empty database first
        strategy = centrevo.RecordingStrategy.none()
        recorder = centrevo.Recorder(temp_db, "dummy", strategy)
        config = centrevo.SimulationConfig(10, 10, 1)
        recorder.record_metadata(config)

        query = centrevo.QueryBuilder(temp_db)
        assert "QueryBuilder" in repr(query)

    def test_list_simulations_empty(self, temp_db):
        """Test listing simulations in empty database."""
        # Create database with one simulation
        strategy = centrevo.RecordingStrategy.none()
        recorder = centrevo.Recorder(temp_db, "test_sim", strategy)
        config = centrevo.SimulationConfig(10, 10, 1)
        recorder.record_metadata(config)

        query = centrevo.QueryBuilder(temp_db)
        sims = query.list_simulations()

        assert isinstance(sims, list)
        assert "test_sim" in sims

    def test_list_simulations_multiple(self, temp_db, test_population):
        """Test listing multiple simulations."""
        # Create multiple simulations
        for sim_name in ["sim1", "sim2", "sim3"]:
            strategy = centrevo.RecordingStrategy.all()
            recorder = centrevo.Recorder(temp_db, sim_name, strategy)
            config = centrevo.SimulationConfig(
            population_size=test_population.size(),
            total_generations=10,
            seed=1
            )
            recorder.record_metadata(config)

        query = centrevo.QueryBuilder(temp_db)
        sims = query.list_simulations()

        assert isinstance(sims, list)
        assert len(sims) == 3
        assert "sim1" in sims
        assert "sim2" in sims
        assert "sim3" in sims

    def test_get_recorded_generations(self, temp_db, test_population):
        """Test getting recorded generations."""
        strategy = centrevo.RecordingStrategy.specific([0, 5, 10, 15])
        recorder = centrevo.Recorder(temp_db, "test_sim", strategy)

        config = centrevo.SimulationConfig(
            population_size=test_population.size(),
            total_generations=20,
            seed=1
        )

        recorder.record_metadata(config)

        # Record specific generations
        for gen in [0, 5, 10, 15]:
            recorder.record_generation(test_population, generation=gen)

        query = centrevo.QueryBuilder(temp_db)
        generations = query.get_recorded_generations("test_sim")

        assert isinstance(generations, list)
        assert len(generations) == 4
        assert 0 in generations
        assert 5 in generations
        assert 10 in generations
        assert 15 in generations

    def test_get_generations_for_nonexistent_simulation(self, temp_db):
        """Test getting generations for non-existent simulation."""
        # Create a simulation first
        strategy = centrevo.RecordingStrategy.none()
        recorder = centrevo.Recorder(temp_db, "real_sim", strategy)
        config = centrevo.SimulationConfig(10, 10, 1)
        recorder.record_metadata(config)

        query = centrevo.QueryBuilder(temp_db)

        # Try to get generations for non-existent simulation
        # May return empty list or raise exception depending on implementation
        try:
            gens = query.get_recorded_generations("nonexistent_sim")
            # If it doesn't raise, should return empty list
            assert isinstance(gens, list)
        except Exception:
            # Exception is also acceptable
            pass


class TestRecordingStrategies:
    """Tests for different recording strategies."""

    def test_every_n_strategy(self, temp_db, test_population):
        """Test EveryN recording strategy."""
        strategy = centrevo.RecordingStrategy.every_n(5)
        recorder = centrevo.Recorder(temp_db, "every_n_test", strategy)

        config = centrevo.SimulationConfig(
            population_size=test_population.size(),
            total_generations=25,
            seed=1
        )

        recorder.record_metadata(config)

        # Record generations 0, 5, 10, 15, 20
        for gen in range(0, 25, 5):
            recorder.record_generation(test_population, generation=gen)

        query = centrevo.QueryBuilder(temp_db)
        generations = query.get_recorded_generations("every_n_test")

        assert len(generations) == 5

    def test_specific_generations_strategy(self, temp_db, test_population):
        """Test Specific generations recording strategy."""
        specific_gens = [0, 1, 10, 99]
        strategy = centrevo.RecordingStrategy.specific(specific_gens)
        recorder = centrevo.Recorder(temp_db, "specific_test", strategy)

        config = centrevo.SimulationConfig(
            population_size=test_population.size(),
            total_generations=100,
            seed=1
        )

        recorder.record_metadata(config)

        for gen in specific_gens:
            recorder.record_generation(test_population, generation=gen)

        query = centrevo.QueryBuilder(temp_db)
        generations = query.get_recorded_generations("specific_test")

        assert set(generations) == set(specific_gens)

    def test_all_strategy_small(self, temp_db, test_population):
        """Test All recording strategy with small number of generations."""
        strategy = centrevo.RecordingStrategy.all()
        recorder = centrevo.Recorder(temp_db, "all_test", strategy)

        config = centrevo.SimulationConfig(
            population_size=test_population.size(),
            total_generations=5,
            seed=1
        )

        recorder.record_metadata(config)

        # Record all 5 generations
        for gen in range(5):
            recorder.record_generation(test_population, generation=gen)

        query = centrevo.QueryBuilder(temp_db)
        generations = query.get_recorded_generations("all_test")

        assert len(generations) == 5
        assert list(range(5)) == sorted(generations)

    def test_none_strategy(self, temp_db):
        """Test None recording strategy."""
        strategy = centrevo.RecordingStrategy.none()
        recorder = centrevo.Recorder(temp_db, "none_test", strategy)

        config = centrevo.SimulationConfig(
            population_size=10,
            total_generations=100,
            seed=1
        )

        # Only metadata should be recorded
        recorder.record_metadata(config)

        query = centrevo.QueryBuilder(temp_db)
        sims = query.list_simulations()

        assert "none_test" in sims


class TestWorkflowIntegration:
    """Tests for complete workflow integration."""

    def test_full_recording_workflow(self, temp_db, test_population):
        """Test complete recording and querying workflow."""
        # Setup simulation
        sim_name = "full_workflow"
        strategy = centrevo.RecordingStrategy.every_n(10)
        recorder = centrevo.Recorder(temp_db, sim_name, strategy)

        config = centrevo.SimulationConfig(
            population_size=test_population.size(),
            total_generations=50,
            seed=12345
        )

        # Record metadata
        recorder.record_metadata(config)

        # Record several generations
        generations_to_record = [0, 10, 20, 30, 40, 50]
        for gen in generations_to_record:
            recorder.record_generation(test_population, generation=gen)

        # Query back
        query = centrevo.QueryBuilder(temp_db)

        # List simulations
        sims = query.list_simulations()
        assert sim_name in sims

        # Get recorded generations
        recorded_gens = query.get_recorded_generations(sim_name)
        assert set(recorded_gens) == set(generations_to_record)

    def test_multiple_simulations_same_db(self, temp_db, test_population):
        """Test recording multiple simulations in same database."""
        sim_names = ["sim_a", "sim_b", "sim_c"]

        for sim_name in sim_names:
            strategy = centrevo.RecordingStrategy.specific([0, 5])
            recorder = centrevo.Recorder(temp_db, sim_name, strategy)

            config = centrevo.SimulationConfig(
            population_size=test_population.size(),
            total_generations=10,
            seed=1
            )

            recorder.record_metadata(config)
            recorder.record_generation(test_population, generation=0)
            recorder.record_generation(test_population, generation=5)

        # Query all simulations
        query = centrevo.QueryBuilder(temp_db)
        all_sims = query.list_simulations()

        assert len(all_sims) == 3
        for name in sim_names:
            assert name in all_sims

            # Each should have recorded generations [0, 5]
            gens = query.get_recorded_generations(name)
            assert set(gens) == {0, 5}


class TestErrorHandling:
    """Tests for error handling in storage operations."""

    def test_invalid_database_path(self):
        """Test with invalid database path."""
        with pytest.raises(Exception):
            centrevo.QueryBuilder("/nonexistent/path/database.db")

    def test_recorder_with_invalid_generation(self, temp_db, test_population):
        """Test recording with negative generation."""
        strategy = centrevo.RecordingStrategy.all()
        recorder = centrevo.Recorder(temp_db, "test", strategy)

        config = centrevo.SimulationConfig(
            population_size=test_population.size(),
            total_generations=10,
            seed=1
        )
        recorder.record_metadata(config)

        # Negative generation might be rejected
        # (depends on implementation)
        try:
            recorder.record_generation(test_population, generation=-1)
        except Exception:
            pass  # Expected to fail
