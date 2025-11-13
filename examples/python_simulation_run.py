#!/usr/bin/env python3
"""
Example: Running Simulations in Python

This example demonstrates the main ways to run evolutionary simulations:
1. Create and run a basic simulation step-by-step
2. Run a simulation for a specific number of generations
3. Run a complete simulation to completion
4. Monitor population changes during simulation
"""

import centrevo as cv

def example_1_basic_simulation():
    """Example 1: Create and run a basic simulation."""
    print("=" * 70)
    print("Example 1: Basic Simulation Setup and Execution")
    print("=" * 70)

    # Set up simulation parameters
    alphabet = cv.Alphabet.dna()
    structure = cv.RepeatStructure(
        alphabet=alphabet,
        init_base=cv.Nucleotide.A(),
        ru_length=171,      # Repeat unit length
        rus_per_hor=12,     # Repeat units per HOR
        hors_per_chr=10,    # HORs per chromosome
        chrs_per_hap=1,     # Chromosomes per haplotype
    )

    # Configure evolutionary processes
    mutation = cv.MutationConfig.uniform(alphabet, 0.0001)  # mutation rate
    recombination = cv.RecombinationConfig.standard(
        break_prob=0.01,
        crossover_prob=0.7,
        gc_extension_prob=0.1
    )
    fitness = cv.FitnessConfig.neutral()  # No selection pressure

    # Configure simulation parameters
    config = cv.SimulationConfig(
        population_size=50,
        total_generations=100,
        seed=42
    )

    # Create simulation
    sim = cv.Simulation(structure, mutation, recombination, fitness, config)
    print(f"Initial state: {sim}")
    print(f"  Generation: {sim.generation()}")
    print(f"  Population size: {sim.population().size()}\n")

    return sim, config


def example_2_step_by_step(sim):
    """Example 2: Run simulation step by step and monitor progress."""
    print("\n" + "=" * 70)
    print("Example 2: Step-by-Step Simulation")
    print("=" * 70)
    print("Running one generation at a time with monitoring:\n")

    for _ in range(10):
        # Advance simulation by one generation
        sim.step()

        # Get current state
        generation = sim.generation()
        population = sim.population()

        # Print progress
        print(f"  Gen {generation:3d}: Population size = {population.size()}")

    print(f"\nAfter 10 steps: {sim}")


def example_3_run_for_generations(sim):
    """Example 3: Run simulation for a batch of generations."""
    print("\n" + "=" * 70)
    print("Example 3: Run for Specific Number of Generations")
    print("=" * 70)

    current_gen = sim.generation()
    print(f"Current generation: {current_gen}\n")

    # Run for 25 more generations
    print("Running for 25 generations...")
    sim.run_for(25)

    new_gen = sim.generation()
    print(f"New generation: {new_gen}")
    print(f"Generations completed: {new_gen - current_gen}\n")


def example_4_run_to_completion(sim, config):
    """Example 4: Run simulation until configured total generations."""
    print("\n" + "=" * 70)
    print("Example 4: Run to Completion")
    print("=" * 70)

    current_gen = sim.generation()
    total_gens = config.total_generations
    remaining = total_gens - current_gen

    print(f"Current generation: {current_gen}")
    print(f"Total configured generations: {total_gens}")
    print(f"Remaining to complete: {remaining}\n")

    print(f"Running {remaining} more generations to completion...")
    sim.run()

    final_gen = sim.generation()
    print(f"Final generation: {final_gen}")
    print(f"Simulation completed: {final_gen >= total_gens}")


def example_5_population_analysis():
    """Example 5: Run simulation with population analysis."""
    print("\n" + "=" * 70)
    print("Example 5: Simulation with Population Analysis")
    print("=" * 70)

    # Create a smaller simulation for demonstration
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

    config = cv.SimulationConfig(
        population_size=100,
        total_generations=50,
        seed=123
    )

    sim = cv.Simulation(structure, mutation, recombination, fitness, config)

    print("Running simulation with periodic analysis:\n")
    print(f"{'Gen':>4} {'Pop Size':>10} {'ND (π)':>12}")
    print("-" * 28)

    # Run and analyze every 10 generations
    for _ in range(5):
        sim.run_for(10)
        generation = sim.generation()
        population = sim.population()

        # Calculate nucleotide diversity
        nd = cv.nucleotide_diversity(population, chromosome_idx=0)

        print(f"{generation:4d} {population.size():10d} {nd:12.6f}")

    print("\nSimulation completed!")


def example_6_compare_runs():
    """Example 6: Compare two simulations with different seeds."""
    print("\n" + "=" * 70)
    print("Example 6: Comparing Different Simulation Runs")
    print("=" * 70)

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

    # Run two simulations with different seeds
    seeds = [42, 999]
    results = []

    for seed in seeds:
        print(f"\nRun with seed={seed}:")
        config = cv.SimulationConfig(
            population_size=50,
            total_generations=30,
            seed=seed
        )

        sim = cv.Simulation(structure, mutation, recombination, fitness, config)
        sim.run()

        population = sim.population()
        nd = cv.nucleotide_diversity(population, 0)
        hd = cv.haplotype_diversity(population, 0)

        results.append({
            'seed': seed,
            'nucleotide_diversity': nd,
            'haplotype_diversity': hd
        })

        print(f"  Nucleotide diversity (π): {nd:.6f}")
        print(f"  Haplotype diversity: {hd:.6f}")

    print("\n" + "-" * 70)
    print("Comparison of runs:")
    for result in results:
        print(f"  Seed {result['seed']:4d}: π={result['nucleotide_diversity']:.6f}, "
              f"H={result['haplotype_diversity']:.6f}")


def main():
    """Run all examples."""
    print("\n" + "=" * 70)
    print("Centrevo Python API: Running Simulations")
    print("=" * 70)

    # Example 1: Basic setup
    sim, config = example_1_basic_simulation()

    # Example 2: Step by step
    example_2_step_by_step(sim)

    # Example 3: Run for N generations
    example_3_run_for_generations(sim)

    # Example 4: Run to completion
    example_4_run_to_completion(sim, config)

    # Example 5: With analysis
    example_5_population_analysis()

    # Example 6: Compare different runs
    example_6_compare_runs()

    print("\n" + "=" * 70)
    print("All simulation examples completed successfully!")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    main()
