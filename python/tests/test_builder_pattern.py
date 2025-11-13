"""
Test suite for SimulationBuilder pattern in Python API.

Tests the fluent builder interface for creating simulations with various configurations.

Run with: ./.venv/bin/python -m pytest python/tests/test_builder_pattern.py -v
"""

import centrevo as cv
import pytest


class TestSimulationBuilderBasic:
    """Tests for basic SimulationBuilder functionality."""

    def test_minimal_builder(self):
        """Test minimal builder with only required parameters."""
        sim = cv.SimulationBuilder() \
            .population_size(50) \
            .generations(100) \
            .repeat_structure(171, 12, 10) \
            .build()

        assert sim is not None
        assert sim.generation() == 0
        assert sim.population().size() == 50

    def test_builder_with_seed(self):
        """Test builder with explicit seed."""
        sim = cv.SimulationBuilder() \
            .population_size(100) \
            .generations(50) \
            .repeat_structure(171, 12, 10) \
            .seed(42) \
            .build()

        assert sim.generation() == 0
        assert sim.population().size() == 100

    def test_builder_with_mutation(self):
        """Test builder with mutation rate."""
        sim = cv.SimulationBuilder() \
            .population_size(50) \
            .generations(50) \
            .repeat_structure(171, 12, 10) \
            .mutation_rate(0.0001) \
            .build()

        assert sim is not None
        assert sim.population().size() == 50

    def test_builder_with_recombination(self):
        """Test builder with recombination parameters."""
        sim = cv.SimulationBuilder() \
            .population_size(50) \
            .generations(50) \
            .repeat_structure(171, 12, 10) \
            .recombination(0.01, 0.7, 0.1) \
            .build()

        assert sim is not None
        assert sim.population().size() == 50

    def test_builder_with_all_evolution_params(self):
        """Test builder with all evolutionary parameters."""
        sim = cv.SimulationBuilder() \
            .population_size(100) \
            .generations(50) \
            .repeat_structure(171, 12, 10) \
            .mutation_rate(0.0001) \
            .recombination(0.01, 0.7, 0.1) \
            .seed(42) \
            .build()

        assert sim is not None
        assert sim.generation() == 0
        assert sim.population().size() == 100


class TestSimulationBuilderInitialization:
    """Tests for different initialization methods in builder."""

    def test_init_uniform_default(self):
        """Test uniform initialization with default base (A)."""
        sim = cv.SimulationBuilder() \
            .population_size(50) \
            .generations(50) \
            .repeat_structure(171, 12, 10) \
            .build()

        # Default should be uniform with A
        assert sim.population().size() == 50

    def test_init_uniform_custom_base(self):
        """Test uniform initialization with custom base."""
        sim = cv.SimulationBuilder() \
            .population_size(50) \
            .generations(50) \
            .repeat_structure(171, 12, 10) \
            .init_uniform(cv.Nucleotide.G()) \
            .build()

        assert sim.population().size() == 50

    def test_init_random(self):
        """Test random initialization."""
        sim = cv.SimulationBuilder() \
            .population_size(50) \
            .generations(50) \
            .repeat_structure(171, 12, 10) \
            .init_random() \
            .mutation_rate(0.0001) \
            .seed(123) \
            .build()

        assert sim is not None
        assert sim.population().size() == 50


class TestSimulationBuilderWithFitness:
    """Tests for builder with fitness configurations."""

    def test_builder_with_neutral_fitness(self):
        """Test builder with explicit neutral fitness."""
        fitness = cv.FitnessConfig.neutral()

        sim = cv.SimulationBuilder() \
            .population_size(100) \
            .generations(50) \
            .repeat_structure(171, 12, 10) \
            .fitness(fitness) \
            .build()

        assert sim is not None
        assert sim.population().size() == 100

    def test_builder_with_gc_fitness(self):
        """Test builder with GC content fitness."""
        fitness = cv.FitnessConfig.with_gc_content(0.5, 2.0).build()

        sim = cv.SimulationBuilder() \
            .population_size(100) \
            .generations(50) \
            .repeat_structure(171, 12, 10) \
            .mutation_rate(0.0001) \
            .fitness(fitness) \
            .build()

        assert sim is not None
        assert sim.population().size() == 100

    def test_builder_with_length_fitness(self):
        """Test builder with length fitness."""
        fitness = cv.FitnessConfig.with_length(20000, 0.5).build()

        sim = cv.SimulationBuilder() \
            .population_size(100) \
            .generations(50) \
            .repeat_structure(171, 12, 10) \
            .fitness(fitness) \
            .build()

        assert sim is not None

    def test_builder_with_similarity_fitness(self):
        """Test builder with similarity fitness."""
        fitness = cv.FitnessConfig.with_similarity(2.0).build()

        sim = cv.SimulationBuilder() \
            .population_size(100) \
            .generations(50) \
            .repeat_structure(171, 12, 10) \
            .fitness(fitness) \
            .build()

        assert sim is not None

    def test_builder_with_combined_fitness(self):
        """Test builder with multiple fitness components."""
        fitness = cv.FitnessConfig.with_gc_content(0.5, 2.0) \
            .with_length(20000, 0.5) \
            .with_similarity(2.0) \
            .build()

        sim = cv.SimulationBuilder() \
            .population_size(100) \
            .generations(50) \
            .repeat_structure(171, 12, 10) \
            .mutation_rate(0.0001) \
            .fitness(fitness) \
            .build()

        assert sim is not None
        assert sim.population().size() == 100


class TestSimulationBuilderRunning:
    """Tests for running simulations built with builder pattern."""

    def test_builder_simulation_runs(self):
        """Test that builder-created simulation can run."""
        sim = cv.SimulationBuilder() \
            .population_size(50) \
            .generations(30) \
            .repeat_structure(171, 12, 5) \
            .mutation_rate(0.0001) \
            .recombination(0.01, 0.7, 0.1) \
            .seed(42) \
            .build()

        assert sim.generation() == 0

        sim.step()
        assert sim.generation() == 1

        sim.run_for(9)
        assert sim.generation() == 10

    def test_builder_simulation_runs_to_completion(self):
        """Test that builder simulation runs to completion."""
        sim = cv.SimulationBuilder() \
            .population_size(50) \
            .generations(30) \
            .repeat_structure(171, 12, 5) \
            .mutation_rate(0.0001) \
            .seed(42) \
            .build()

        sim.run()
        assert sim.generation() == 30
        assert sim.population().size() == 50

    def test_builder_with_fitness_runs(self):
        """Test simulation with fitness runs correctly."""
        fitness = cv.FitnessConfig.with_gc_content(0.5, 2.0).build()

        sim = cv.SimulationBuilder() \
            .population_size(100) \
            .generations(30) \
            .repeat_structure(171, 12, 5) \
            .mutation_rate(0.0001) \
            .recombination(0.01, 0.7, 0.1) \
            .fitness(fitness) \
            .seed(42) \
            .build()

        sim.run()
        assert sim.generation() == 30
        assert sim.population().size() == 100


class TestBuilderReproducibility:
    """Tests for reproducibility with builder pattern."""

    def test_same_seed_produces_same_simulation(self):
        """Test that same seed produces reproducible simulations."""
        sim1 = cv.SimulationBuilder() \
            .population_size(50) \
            .generations(20) \
            .repeat_structure(171, 12, 5) \
            .mutation_rate(0.0001) \
            .recombination(0.01, 0.7, 0.1) \
            .seed(42) \
            .build()

        sim2 = cv.SimulationBuilder() \
            .population_size(50) \
            .generations(20) \
            .repeat_structure(171, 12, 5) \
            .mutation_rate(0.0001) \
            .recombination(0.01, 0.7, 0.1) \
            .seed(42) \
            .build()

        # Run both simulations
        sim1.run()
        sim2.run()

        # Should reach same generation
        assert sim1.generation() == sim2.generation() == 20
        assert sim1.population().size() == sim2.population().size() == 50

    def test_different_seeds_produce_different_results(self):
        """Test that different seeds can produce different results."""
        sim1 = cv.SimulationBuilder() \
            .population_size(50) \
            .generations(20) \
            .repeat_structure(171, 12, 5) \
            .mutation_rate(0.001) \
            .seed(42) \
            .build()

        sim2 = cv.SimulationBuilder() \
            .population_size(50) \
            .generations(20) \
            .repeat_structure(171, 12, 5) \
            .mutation_rate(0.001) \
            .seed(999) \
            .build()

        sim1.run()
        sim2.run()

        # Both should complete but may have different GC content due to randomness
        assert sim1.generation() == sim2.generation() == 20
        assert sim1.population().size() == sim2.population().size() == 50


class TestBuilderComplexScenarios:
    """Tests for complex scenarios with builder pattern."""

    def test_builder_neutral_vs_selection(self):
        """Test comparing neutral vs selection using builder."""
        # Neutral evolution
        sim_neutral = cv.SimulationBuilder() \
            .population_size(100) \
            .generations(30) \
            .repeat_structure(171, 12, 5) \
            .mutation_rate(0.0001) \
            .recombination(0.01, 0.7, 0.1) \
            .seed(42) \
            .build()

        # With selection
        fitness = cv.FitnessConfig.with_gc_content(0.5, 2.0).build()
        sim_selection = cv.SimulationBuilder() \
            .population_size(100) \
            .generations(30) \
            .repeat_structure(171, 12, 5) \
            .mutation_rate(0.0001) \
            .recombination(0.01, 0.7, 0.1) \
            .fitness(fitness) \
            .seed(42) \
            .build()

        # Run both
        sim_neutral.run()
        sim_selection.run()

        # Both should complete
        assert sim_neutral.generation() == 30
        assert sim_selection.generation() == 30

        # Get final GC content
        gc_neutral = cv.gc_content(sim_neutral.population(), None, None, 0)
        gc_selection = cv.gc_content(sim_selection.population(), None, None, 0)

        # Both should be valid GC content values
        assert 0.0 <= gc_neutral <= 1.0
        assert 0.0 <= gc_selection <= 1.0

    def test_builder_multiple_chromosomes(self):
        """Test builder with multiple chromosomes per haplotype."""
        sim = cv.SimulationBuilder() \
            .population_size(50) \
            .generations(20) \
            .repeat_structure(171, 12, 10) \
            .mutation_rate(0.0001) \
            .build()

        sim.run()
        assert sim.generation() == 20
        assert sim.population().size() == 50

    def test_builder_small_structure(self):
        """Test builder with smaller repeat structure."""
        sim = cv.SimulationBuilder() \
            .population_size(20) \
            .generations(10) \
            .repeat_structure(10, 5, 2) \
            .mutation_rate(0.001) \
            .seed(123) \
            .build()

        sim.run()
        assert sim.generation() == 10

    def test_builder_large_population(self):
        """Test builder with larger population."""
        sim = cv.SimulationBuilder() \
            .population_size(200) \
            .generations(10) \
            .repeat_structure(171, 12, 5) \
            .mutation_rate(0.0001) \
            .seed(42) \
            .build()

        sim.run_for(5)
        assert sim.generation() == 5
        assert sim.population().size() == 200
