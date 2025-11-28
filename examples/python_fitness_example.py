#!/usr/bin/env python3
"""
Example: Using Fitness Functions with Selection

This example demonstrates how to configure fitness functions for selection
in evolutionary simulations.
"""

import centrevo as cv


def example_1_neutral_fitness():
    """Example 1: Neutral evolution (no selection)."""
    print("=" * 70)
    print("Example 1: Neutral Fitness (No Selection)")
    print("=" * 70)

    # Simple neutral fitness - no selection pressure
    fitness = cv.FitnessConfig.neutral()
    print(f"Fitness config: {fitness}")
    print("All individuals have equal fitness (no selection)\n")

    return fitness


def example_2_gc_content_fitness():
    """Example 2: GC content-based selection."""
    print("=" * 70)
    print("Example 2: GC Content Fitness")
    print("=" * 70)

    # Selection favoring sequences with ~50% GC content
    # optimum: 0.5 (50% GC)
    # concentration: 2.0 (controls sharpness of selection)
    fitness = cv.FitnessConfig.with_gc_content(optimum=0.5, concentration=2.0).build()
    print(f"Fitness config: {fitness}")
    print("Selection favors sequences with 50% GC content")
    print("Higher concentration = stronger selection\n")

    return fitness


def example_3_length_fitness():
    """Example 3: Length-based selection."""
    print("=" * 70)
    print("Example 3: Length-Based Fitness")
    print("=" * 70)

    # Selection favoring sequences near optimal length
    # optimum: 20000 bases
    # std_dev: 0.5 (controls width of fitness curve)
    fitness = cv.FitnessConfig.with_length(optimum=20000, std_dev=0.5).build()
    print(f"Fitness config: {fitness}")
    print("Selection favors sequences around 20,000 bases")
    print("Lower std_dev = stronger selection for exact length\n")

    return fitness


def example_4_similarity_fitness():
    """Example 4: Sequence similarity-based selection."""
    print("=" * 70)
    print("Example 4: Sequence Similarity Fitness")
    print("=" * 70)

    # Selection favoring similarity between haplotypes
    # shape: 2.0 (controls decline rate)
    fitness = cv.FitnessConfig.with_similarity(shape=2.0).build()
    print(f"Fitness config: {fitness}")
    print("Selection favors individuals with similar haplotypes")
    print("Higher shape = steeper decline with dissimilarity\n")

    return fitness


def example_5_combined_fitness():
    """Example 5: Multiple fitness components combined."""
    print("=" * 70)
    print("Example 5: Combined Fitness Functions")
    print("=" * 70)

    # Combine multiple selection pressures
    # All fitness values are multiplied together
    fitness = cv.FitnessConfig.with_gc_content(optimum=0.5, concentration=2.0) \
        .with_length(optimum=20000, std_dev=0.5) \
        .with_similarity(shape=2.0) \
        .build()

    print(f"Fitness config: {fitness}")
    print("Selection combines:")
    print("  - GC content near 50%")
    print("  - Length near 20,000 bases")
    print("  - Similarity between haplotypes")
    print("\nFitness = gc_fitness × length_fitness × similarity_fitness\n")

    return fitness


def example_6_run_simulation_with_selection():
    """Example 6: Run a simulation with selection."""
    print("=" * 70)
    print("Example 6: Running Simulation with Selection")
    print("=" * 70)

    # Set up simulation parameters
    alphabet = cv.Alphabet.dna()
    structure = cv.RepeatStructure(
        alphabet=alphabet,
        init_base=cv.Nucleotide.A(),
        ru_length=171,
        rus_per_hor=12,
        hors_per_chr=10,
        chrs_per_hap=1,
    )

    mutation = cv.MutationConfig.uniform(alphabet, 0.0001)
    recombination = cv.RecombinationConfig.standard(0.01, 0.7, 0.1)

    # Configure selection for GC content
    fitness = cv.FitnessConfig.with_gc_content(optimum=0.5, concentration=2.0).build()

    config = cv.SimulationConfig(
        population_size=100,
        total_generations=50,
        seed=42
    )

    print("Running simulation with GC content selection...")
    print(f"  Population size: {config.population_size}")
    print(f"  Generations: {config.total_generations}")
    print("  Selection: GC content (optimum=0.5, concentration=2.0)\n")

    # Create and run simulation
    sim = cv.Simulation(structure, mutation, recombination, fitness, config)

    # Track GC content over time
    print(f"{'Gen':>4} {'Mean GC':>10}")
    print("-" * 16)

    for gen in range(0, 51, 10):
        if gen > 0:
            sim.run_for(10)

        population = sim.population()
        gc = cv.gc_content(population, individual_idx=None, haplotype_idx=None, chromosome_idx=0)
        print(f"{gen:4d} {gc:10.4f}")

    print("\nWith selection, GC content should approach the optimum (0.5)\n")


def example_7_compare_neutral_vs_selection():
    """Example 7: Compare neutral vs selection."""
    print("=" * 70)
    print("Example 7: Neutral vs Selection Comparison")
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

    # Run with neutral fitness
    config_neutral = cv.SimulationConfig(population_size=50, total_generations=30, seed=42)
    fitness_neutral = cv.FitnessConfig.neutral()
    sim_neutral = cv.Simulation(structure, mutation, recombination, fitness_neutral, config_neutral)
    sim_neutral.run()

    # Run with GC selection
    config_selection = cv.SimulationConfig(population_size=50, total_generations=30, seed=42)
    fitness_selection = cv.FitnessConfig.with_gc_content(0.5, 2.0).build()
    sim_selection = cv.Simulation(
        structure,
        mutation,
        recombination,
        fitness_selection,
        config_selection
    )
    sim_selection.run()

    # Compare final GC content
    gc_neutral = cv.gc_content(sim_neutral.population(), None, None, 0)
    gc_selection = cv.gc_content(sim_selection.population(), None, None, 0)

    print(f"After {config_neutral.total_generations} generations:")
    print(f"  Neutral evolution:  GC = {gc_neutral:.4f}")
    print(f"  With selection:     GC = {gc_selection:.4f} (target: 0.5)")
    print("\nDifference from target:")
    print(f"  Neutral:   {abs(gc_neutral - 0.5):.4f}")
    print(f"  Selection: {abs(gc_selection - 0.5):.4f}")
    print("\nSelection drives GC content toward the optimum!\n")


def main():
    """Run all examples."""
    print("\n" + "=" * 70)
    print("Centrevo Python API: Fitness Functions and Selection")
    print("=" * 70 + "\n")

    # Basic fitness configurations
    example_1_neutral_fitness()
    example_2_gc_content_fitness()
    example_3_length_fitness()
    example_4_similarity_fitness()
    example_5_combined_fitness()

    # Running simulations with selection
    example_6_run_simulation_with_selection()
    example_7_compare_neutral_vs_selection()

    print("=" * 70)
    print("All fitness examples completed successfully!")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    main()
