"""
Test suite for running simulations in Python.

Tests the ability to run simulations from start to finish using various methods:
- step() - Run one generation at a time
- run_for() - Run for N generations
- run() - Run until completion

Run with: python -m pytest python/tests/test_simulation_run.py
"""

import centrevo as cv


class TestSimulationRunning:
    """Tests for simulation running functionality."""
    
    def test_step_increments_generation(self):
        """Test that step() correctly increments generation counter."""
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
        config = cv.SimulationConfig(population_size=10, total_generations=100, seed=42)
        
        sim = cv.Simulation(structure, mutation, recombination, fitness, config)
        
        # Initial state
        assert sim.generation() == 0
        
        # After one step
        sim.step()
        assert sim.generation() == 1
        
        # After multiple steps
        for _ in range(9):
            sim.step()
        assert sim.generation() == 10
    
    def test_run_for_multiple_generations(self):
        """Test that run_for() runs exactly the specified number of generations."""
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
        config = cv.SimulationConfig(population_size=10, total_generations=100, seed=42)
        
        sim = cv.Simulation(structure, mutation, recombination, fitness, config)
        
        # Run for 25 generations
        sim.run_for(25)
        assert sim.generation() == 25
        
        # Run for 30 more generations
        sim.run_for(30)
        assert sim.generation() == 55
        
        # Run for 45 more generations
        sim.run_for(45)
        assert sim.generation() == 100
    
    def test_run_completes_simulation(self):
        """Test that run() completes the simulation to configured total."""
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
        config = cv.SimulationConfig(population_size=10, total_generations=50, seed=42)
        
        sim = cv.Simulation(structure, mutation, recombination, fitness, config)
        
        # Run to completion
        sim.run()
        
        # Should reach configured total
        assert sim.generation() == 50
    
    def test_run_from_partial_state(self):
        """Test that run() continues for the configured total_generations from start."""
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
        config = cv.SimulationConfig(population_size=10, total_generations=60, seed=42)
        
        sim = cv.Simulation(structure, mutation, recombination, fitness, config)
        
        # Run partially
        sim.run_for(20)
        assert sim.generation() == 20
        
        # Continue with run() - runs for the full total_generations from start
        sim.run()
        assert sim.generation() == 80  # 20 + 60 generations
    
    def test_population_size_maintained(self):
        """Test that population size is maintained throughout simulation."""
        population_size = 25
        
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
        config = cv.SimulationConfig(
            population_size=population_size,
            total_generations=40,
            seed=42
        )
        
        sim = cv.Simulation(structure, mutation, recombination, fitness, config)
        
        # Check population size at various points
        assert sim.population().size() == population_size
        
        sim.run_for(10)
        assert sim.population().size() == population_size
        
        sim.run_for(10)
        assert sim.population().size() == population_size
        
        sim.run()
        assert sim.population().size() == population_size
    
    def test_simulation_run_start_to_finish(self):
        """Test complete simulation run from start to finish."""
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
        fitness = cv.FitnessConfig.neutral()
        
        total_generations = 100
        population_size = 50
        config = cv.SimulationConfig(
            population_size=population_size,
            total_generations=total_generations,
            seed=12345
        )
        
        # Create simulation
        sim = cv.Simulation(structure, mutation, recombination, fitness, config)
        
        # Verify initial state
        assert sim.generation() == 0
        assert sim.population().size() == population_size
        
        # Run to completion
        sim.run()
        
        # Verify final state
        assert sim.generation() == total_generations
        assert sim.population().size() == population_size
        
        # Verify population is not empty and has valid data
        pop = sim.population()
        assert pop.size() > 0
        assert pop.generation() == total_generations
    
    def test_mixed_running_methods(self):
        """Test mixing different running methods (step, run_for, run)."""
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
        config = cv.SimulationConfig(population_size=10, total_generations=60, seed=42)
        
        sim = cv.Simulation(structure, mutation, recombination, fitness, config)
        
        # Mix different methods
        sim.step()  # 1 generation
        assert sim.generation() == 1
        
        sim.run_for(19)  # +19 = 20 generations
        assert sim.generation() == 20
        
        sim.step()  # +1 = 21 generations
        assert sim.generation() == 21
        
        sim.run_for(9)  # +9 = 30 generations
        assert sim.generation() == 30
        
        sim.run()  # Runs total_generations (60) more from current state
        assert sim.generation() == 90  # 30 + 60
    
    def test_large_population_simulation(self):
        """Test simulation with larger population and more generations."""
        alphabet = cv.Alphabet.dna()
        structure = cv.RepeatStructure(
            alphabet=alphabet,
            init_base=cv.Nucleotide.A(),
            ru_length=171,
            rus_per_hor=12,
            hors_per_chr=10,
            chrs_per_hap=1,
        )
        
        mutation = cv.MutationConfig.uniform(alphabet, 0.00005)
        recombination = cv.RecombinationConfig.standard(0.001, 0.7, 0.1)
        fitness = cv.FitnessConfig.neutral()
        config = cv.SimulationConfig(
            population_size=200,
            total_generations=200,
            seed=999
        )
        
        sim = cv.Simulation(structure, mutation, recombination, fitness, config)
        
        # Track progress
        checkpoints = [50, 100, 150, 200]
        for checkpoint in checkpoints:
            remaining = checkpoint - sim.generation()
            sim.run_for(remaining)
            assert sim.generation() == checkpoint
            assert sim.population().size() == 200
    
    def test_reproducibility_with_seed(self):
        """Test that simulations with same seed are reproducible."""
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
        
        seed = 54321
        config1 = cv.SimulationConfig(population_size=10, total_generations=30, seed=seed)
        config2 = cv.SimulationConfig(population_size=10, total_generations=30, seed=seed)
        
        # Run first simulation
        sim1 = cv.Simulation(structure, mutation, recombination, fitness, config1)
        sim1.run()
        
        # Run second simulation with same seed
        sim2 = cv.Simulation(structure, mutation, recombination, fitness, config2)
        sim2.run()
        
        # Both should reach same final generation
        assert sim1.generation() == sim2.generation() == 30
        # Both should have same final population size
        assert sim1.population().size() == sim2.population().size() == 10
    
    def test_get_population_during_run(self):
        """Test retrieving population state during simulation."""
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
        config = cv.SimulationConfig(population_size=15, total_generations=50, seed=42)
        
        sim = cv.Simulation(structure, mutation, recombination, fitness, config)
        
        populations_at_generations = {}
        
        # Collect population snapshots
        for _ in range(5):
            sim.run_for(10)
            gen = sim.generation()
            pop = sim.population()
            populations_at_generations[gen] = pop.size()
            
            # Verify generation and population match
            assert pop.generation() == gen
        
        # Verify we captured all checkpoints
        assert populations_at_generations == {10: 15, 20: 15, 30: 15, 40: 15, 50: 15}
