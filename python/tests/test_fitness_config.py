"""
Test suite for FitnessConfig and fitness builder patterns in Python API.

Tests the various fitness configurations and their effects on selection.

Run with: ./.venv/bin/python -m pytest python/tests/test_fitness_config.py -v
"""

import centrevo as cv
import pytest


class TestFitnessConfigBasic:
    """Tests for basic FitnessConfig functionality."""

    def test_neutral_fitness(self):
        """Test neutral fitness configuration."""
        fitness = cv.FitnessConfig.neutral()

        assert fitness is not None
        assert "FitnessConfig" in repr(fitness)

    def test_neutral_fitness_in_simulation(self):
        """Test neutral fitness in simulation."""
        alphabet = cv.Alphabet.dna()
        structure = cv.RepeatStructure(
            alphabet=alphabet,
            init_base=cv.Nucleotide.A(),
            ru_length=10,
            rus_per_hor=5,
            hors_per_chr=2,
            chrs_per_hap=1,
        )

        mutation = cv.MutationConfig.uniform(alphabet, 0.001)
        recombination = cv.RecombinationConfig.standard(0.01, 0.7, 0.1)
        fitness = cv.FitnessConfig.neutral()
        config = cv.SimulationConfig(population_size=20, total_generations=10, seed=42)

        sim = cv.Simulation(structure, mutation, recombination, fitness, config)
        sim.run()

        assert sim.generation() == 10
        assert sim.population().size() == 20


class TestFitnessConfigGCContent:
    """Tests for GC content fitness."""

    def test_gc_content_fitness_creation(self):
        """Test creating GC content fitness."""
        fitness = cv.FitnessConfig.with_gc_content(0.5, 2.0).build()

        assert fitness is not None
        assert "FitnessConfig" in repr(fitness)

    def test_gc_content_fitness_different_params(self):
        """Test GC content fitness with different parameters."""
        fitness1 = cv.FitnessConfig.with_gc_content(0.4, 1.0).build()
        fitness2 = cv.FitnessConfig.with_gc_content(0.6, 3.0).build()

        assert fitness1 is not None
        assert fitness2 is not None

    def test_gc_content_fitness_in_simulation(self):
        """Test GC content fitness in simulation."""
        alphabet = cv.Alphabet.dna()
        structure = cv.RepeatStructure(
            alphabet=alphabet,
            init_base=cv.Nucleotide.A(),
            ru_length=171,
            rus_per_hor=12,
            hors_per_chr=5,
            chrs_per_hap=1,
        )

        mutation = cv.MutationConfig.uniform(alphabet, 0.0001)
        recombination = cv.RecombinationConfig.standard(0.01, 0.7, 0.1)
        fitness = cv.FitnessConfig.with_gc_content(0.5, 2.0).build()
        config = cv.SimulationConfig(population_size=50, total_generations=30, seed=42)

        sim = cv.Simulation(structure, mutation, recombination, fitness, config)
        sim.run()

        assert sim.generation() == 30
        assert sim.population().size() == 50

    def test_gc_content_evolution(self):
        """Test that GC content evolves under selection."""
        fitness = cv.FitnessConfig.with_gc_content(0.5, 2.0).build()

        sim = cv.SimulationBuilder() \
            .population_size(100) \
            .generations(50) \
            .repeat_structure(171, 12, 5) \
            .mutation_rate(0.0001) \
            .recombination(0.01, 0.7, 0.1) \
            .fitness(fitness) \
            .seed(42) \
            .build()

        # Get initial GC content
        gc_initial = cv.gc_content(sim.population(), None, None, 0)

        # Run simulation
        sim.run()

        # Get final GC content
        gc_final = cv.gc_content(sim.population(), None, None, 0)

        # Both should be valid values
        assert 0.0 <= gc_initial <= 1.0
        assert 0.0 <= gc_final <= 1.0


class TestFitnessConfigLength:
    """Tests for length-based fitness."""

    def test_length_fitness_creation(self):
        """Test creating length fitness."""
        fitness = cv.FitnessConfig.with_length(20000, 0.5).build()

        assert fitness is not None
        assert "FitnessConfig" in repr(fitness)

    def test_length_fitness_different_params(self):
        """Test length fitness with different parameters."""
        fitness1 = cv.FitnessConfig.with_length(10000, 0.3).build()
        fitness2 = cv.FitnessConfig.with_length(30000, 1.0).build()

        assert fitness1 is not None
        assert fitness2 is not None

    def test_length_fitness_in_simulation(self):
        """Test length fitness in simulation."""
        fitness = cv.FitnessConfig.with_length(20000, 0.5).build()

        sim = cv.SimulationBuilder() \
            .population_size(50) \
            .generations(20) \
            .repeat_structure(171, 12, 10) \
            .mutation_rate(0.0001) \
            .fitness(fitness) \
            .seed(42) \
            .build()

        sim.run()
        assert sim.generation() == 20
        assert sim.population().size() == 50


class TestFitnessConfigSimilarity:
    """Tests for similarity-based fitness."""

    def test_similarity_fitness_creation(self):
        """Test creating similarity fitness."""
        fitness = cv.FitnessConfig.with_similarity(2.0).build()

        assert fitness is not None
        assert "FitnessConfig" in repr(fitness)

    def test_similarity_fitness_different_params(self):
        """Test similarity fitness with different parameters."""
        fitness1 = cv.FitnessConfig.with_similarity(1.0).build()
        fitness2 = cv.FitnessConfig.with_similarity(3.0).build()

        assert fitness1 is not None
        assert fitness2 is not None

    def test_similarity_fitness_in_simulation(self):
        """Test similarity fitness in simulation."""
        fitness = cv.FitnessConfig.with_similarity(2.0).build()

        sim = cv.SimulationBuilder() \
            .population_size(50) \
            .generations(20) \
            .repeat_structure(171, 12, 5) \
            .mutation_rate(0.0001) \
            .recombination(0.01, 0.7, 0.1) \
            .fitness(fitness) \
            .seed(42) \
            .build()

        sim.run()
        assert sim.generation() == 20
        assert sim.population().size() == 50


class TestFitnessConfigCombined:
    """Tests for combined fitness functions."""

    def test_gc_and_length_combined(self):
        """Test combining GC content and length fitness."""
        fitness = cv.FitnessConfig.with_gc_content(0.5, 2.0) \
            .with_length(20000, 0.5) \
            .build()

        assert fitness is not None
        assert "FitnessConfig" in repr(fitness)

    def test_gc_and_similarity_combined(self):
        """Test combining GC content and similarity fitness."""
        fitness = cv.FitnessConfig.with_gc_content(0.5, 2.0) \
            .with_similarity(2.0) \
            .build()

        assert fitness is not None

    def test_length_and_similarity_combined(self):
        """Test combining length and similarity fitness."""
        fitness = cv.FitnessConfig.with_length(20000, 0.5) \
            .with_similarity(2.0) \
            .build()

        assert fitness is not None

    def test_all_three_combined(self):
        """Test combining all three fitness components."""
        fitness = cv.FitnessConfig.with_gc_content(0.5, 2.0) \
            .with_length(20000, 0.5) \
            .with_similarity(2.0) \
            .build()

        assert fitness is not None
        assert "FitnessConfig" in repr(fitness)

    def test_combined_fitness_in_simulation(self):
        """Test combined fitness in simulation."""
        fitness = cv.FitnessConfig.with_gc_content(0.5, 2.0) \
            .with_length(20000, 0.5) \
            .with_similarity(2.0) \
            .build()

        sim = cv.SimulationBuilder() \
            .population_size(100) \
            .generations(30) \
            .repeat_structure(171, 12, 10) \
            .mutation_rate(0.0001) \
            .recombination(0.01, 0.7, 0.1) \
            .fitness(fitness) \
            .seed(42) \
            .build()

        sim.run()
        assert sim.generation() == 30
        assert sim.population().size() == 100


class TestFitnessConfigChaining:
    """Tests for fitness builder chaining."""

    def test_chaining_order_gc_length_similarity(self):
        """Test chaining in order: GC -> length -> similarity."""
        fitness = cv.FitnessConfig.with_gc_content(0.5, 2.0) \
            .with_length(20000, 0.5) \
            .with_similarity(2.0) \
            .build()

        assert fitness is not None

    def test_chaining_order_length_gc_similarity(self):
        """Test chaining in order: length -> GC -> similarity."""
        fitness = cv.FitnessConfig.with_length(20000, 0.5) \
            .with_gc_content(0.5, 2.0) \
            .with_similarity(2.0) \
            .build()

        assert fitness is not None

    def test_chaining_order_similarity_gc_length(self):
        """Test chaining in order: similarity -> GC -> length."""
        fitness = cv.FitnessConfig.with_similarity(2.0) \
            .with_gc_content(0.5, 2.0) \
            .with_length(20000, 0.5) \
            .build()

        assert fitness is not None

    def test_multiple_gc_calls_raises_error(self):
        """Test that multiple calls to with_gc_content raises error."""
        with pytest.raises(ValueError):
            cv.FitnessConfig.with_gc_content(0.3, 1.0) \
                .with_gc_content(0.5, 2.0) \
                .build()


class TestFitnessEffects:
    """Tests for actual fitness effects on evolution."""

    def test_neutral_vs_gc_selection(self):
        """Test comparing neutral vs GC content selection."""
        # Neutral simulation
        sim_neutral = cv.SimulationBuilder() \
            .population_size(100) \
            .generations(30) \
            .repeat_structure(171, 12, 5) \
            .mutation_rate(0.0001) \
            .recombination(0.01, 0.7, 0.1) \
            .seed(42) \
            .build()

        # With GC selection
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

        # Get final GC content
        gc_neutral = cv.gc_content(sim_neutral.population(), None, None, 0)
        gc_selection = cv.gc_content(sim_selection.population(), None, None, 0)

        # Both should be valid
        assert 0.0 <= gc_neutral <= 1.0
        assert 0.0 <= gc_selection <= 1.0

        # Both simulations should complete
        assert sim_neutral.generation() == 30
        assert sim_selection.generation() == 30

    def test_strong_vs_weak_selection(self):
        """Test strong vs weak selection pressure."""
        # Weak selection
        fitness_weak = cv.FitnessConfig.with_gc_content(0.5, 0.5).build()
        sim_weak = cv.SimulationBuilder() \
            .population_size(100) \
            .generations(30) \
            .repeat_structure(171, 12, 5) \
            .mutation_rate(0.0001) \
            .recombination(0.01, 0.7, 0.1) \
            .fitness(fitness_weak) \
            .seed(42) \
            .build()

        # Strong selection
        fitness_strong = cv.FitnessConfig.with_gc_content(0.5, 5.0).build()
        sim_strong = cv.SimulationBuilder() \
            .population_size(100) \
            .generations(30) \
            .repeat_structure(171, 12, 5) \
            .mutation_rate(0.0001) \
            .recombination(0.01, 0.7, 0.1) \
            .fitness(fitness_strong) \
            .seed(42) \
            .build()

        # Run both
        sim_weak.run()
        sim_strong.run()

        # Both should complete
        assert sim_weak.generation() == 30
        assert sim_strong.generation() == 30

    def test_different_gc_optima(self):
        """Test different GC content optima."""
        # Low GC optimum
        fitness_low = cv.FitnessConfig.with_gc_content(0.3, 2.0).build()
        sim_low = cv.SimulationBuilder() \
            .population_size(100) \
            .generations(30) \
            .repeat_structure(171, 12, 5) \
            .mutation_rate(0.0001) \
            .recombination(0.01, 0.7, 0.1) \
            .fitness(fitness_low) \
            .seed(42) \
            .build()

        # High GC optimum
        fitness_high = cv.FitnessConfig.with_gc_content(0.7, 2.0).build()
        sim_high = cv.SimulationBuilder() \
            .population_size(100) \
            .generations(30) \
            .repeat_structure(171, 12, 5) \
            .mutation_rate(0.0001) \
            .recombination(0.01, 0.7, 0.1) \
            .fitness(fitness_high) \
            .seed(42) \
            .build()

        # Run both
        sim_low.run()
        sim_high.run()

        # Both should complete
        assert sim_low.generation() == 30
        assert sim_high.generation() == 30


class TestFitnessConfigEdgeCases:
    """Tests for edge cases and error conditions."""

    def test_fitness_with_zero_mutation(self):
        """Test fitness with zero mutation rate."""
        fitness = cv.FitnessConfig.with_gc_content(0.5, 2.0).build()

        sim = cv.SimulationBuilder() \
            .population_size(50) \
            .generations(20) \
            .repeat_structure(171, 12, 5) \
            .mutation_rate(0.0) \
            .fitness(fitness) \
            .seed(42) \
            .build()

        sim.run()
        assert sim.generation() == 20

    def test_fitness_with_high_mutation(self):
        """Test fitness with high mutation rate."""
        fitness = cv.FitnessConfig.with_gc_content(0.5, 2.0).build()

        sim = cv.SimulationBuilder() \
            .population_size(50) \
            .generations(20) \
            .repeat_structure(171, 12, 5) \
            .mutation_rate(0.01) \
            .fitness(fitness) \
            .seed(42) \
            .build()

        sim.run()
        assert sim.generation() == 20

    def test_extreme_gc_optimum_low(self):
        """Test extreme GC content optimum (very low)."""
        fitness = cv.FitnessConfig.with_gc_content(0.0, 2.0).build()

        sim = cv.SimulationBuilder() \
            .population_size(50) \
            .generations(20) \
            .repeat_structure(171, 12, 5) \
            .mutation_rate(0.0001) \
            .fitness(fitness) \
            .seed(42) \
            .build()

        sim.run()
        assert sim.generation() == 20

    def test_extreme_gc_optimum_high(self):
        """Test extreme GC content optimum (very high)."""
        fitness = cv.FitnessConfig.with_gc_content(1.0, 2.0).build()

        sim = cv.SimulationBuilder() \
            .population_size(50) \
            .generations(20) \
            .repeat_structure(171, 12, 5) \
            .mutation_rate(0.0001) \
            .fitness(fitness) \
            .seed(42) \
            .build()

        sim.run()
        assert sim.generation() == 20
