#!/usr/bin/env python3
"""
Example: Using SimulationBuilder Pattern

This example demonstrates how to use the SimulationBuilder to create
simulations with a fluent, ergonomic API.
"""

import centrevo as cv


def example_1_minimal_builder():
    """Example 1: Minimal simulation using builder pattern."""
    print("=" * 70)
    print("Example 1: Minimal Builder Pattern")
    print("=" * 70)
    
    # Simplest possible simulation - only required parameters
    sim = cv.SimulationBuilder() \
        .population_size(50) \
        .generations(100) \
        .repeat_structure(171, 12, 10) \
        .build()
    
    print(f"Created simulation: {sim}")
    print(f"  Population size: {sim.population().size()}")
    print(f"  Initial generation: {sim.generation()}")
    print("\nDefaults applied:")
    print("  - Mutation rate: 0.0 (no mutation)")
    print("  - Recombination: None")
    print("  - Fitness: Neutral (no selection)")
    print("  - Initialization: Uniform with base A")
    print("  - Chromosomes per haplotype: 1\n")
    
    return sim


def example_2_with_evolution():
    """Example 2: Builder with mutation and recombination."""
    print("=" * 70)
    print("Example 2: Builder with Evolutionary Parameters")
    print("=" * 70)
    
    sim = cv.SimulationBuilder() \
        .population_size(100) \
        .generations(50) \
        .repeat_structure(171, 12, 10) \
        .mutation_rate(0.0001) \
        .recombination(0.01, 0.7, 0.1) \
        .seed(42) \
        .build()
    
    print(f"Created simulation with evolution: {sim}")
    print("  Mutation rate: 0.0001")
    print("  Recombination: break=0.01, crossover=0.7, gc_ext=0.1")
    print("  Random seed: 42 (reproducible)\n")
    
    return sim


def example_3_with_fitness():
    """Example 3: Builder with fitness/selection."""
    print("=" * 70)
    print("Example 3: Builder with Fitness (Selection)")
    print("=" * 70)
    
    # Create fitness configuration
    fitness = cv.FitnessConfig.with_gc_content(0.5, 2.0).build()
    
    # Build simulation with selection
    sim = cv.SimulationBuilder() \
        .population_size(100) \
        .generations(50) \
        .repeat_structure(171, 12, 10) \
        .mutation_rate(0.0001) \
        .recombination(0.01, 0.7, 0.1) \
        .fitness(fitness) \
        .seed(42) \
        .build()
    
    print(f"Created simulation with selection: {sim}")
    print(f"  Fitness: {fitness}")
    print("  Selection favors GC content near 50%\n")
    
    return sim


def example_4_combined_fitness():
    """Example 4: Builder with multiple fitness components."""
    print("=" * 70)
    print("Example 4: Multiple Fitness Components")
    print("=" * 70)
    
    # Build complex fitness configuration
    fitness = cv.FitnessConfig.with_gc_content(0.5, 2.0) \
        .with_length(20000, 0.5) \
        .with_similarity(2.0) \
        .build()
    
    sim = cv.SimulationBuilder() \
        .population_size(100) \
        .generations(50) \
        .repeat_structure(171, 12, 10) \
        .mutation_rate(0.0001) \
        .recombination(0.01, 0.7, 0.1) \
        .fitness(fitness) \
        .build()
    
    print(f"Created simulation: {sim}")
    print(f"  Fitness: {fitness}")
    print("  Selection combines:")
    print("    - GC content (optimum=0.5, concentration=2.0)")
    print("    - Length (optimum=20000, std_dev=0.5)")
    print("    - Similarity (shape=2.0)\n")
    
    return sim


def example_5_random_initialization():
    """Example 5: Random initialization instead of uniform."""
    print("=" * 70)
    print("Example 5: Random Initialization")
    print("=" * 70)
    
    sim = cv.SimulationBuilder() \
        .population_size(50) \
        .generations(50) \
        .repeat_structure(171, 12, 10) \
        .init_random() \
        .mutation_rate(0.0001) \
        .seed(123) \
        .build()
    
    print(f"Created simulation: {sim}")
    print("  Initialization: Random (not uniform)")
    print("  Each base position gets a random nucleotide\n")
    
    return sim


def example_6_custom_base():
    """Example 6: Uniform initialization with custom base."""
    print("=" * 70)
    print("Example 6: Custom Initial Base")
    print("=" * 70)
    
    sim = cv.SimulationBuilder() \
        .population_size(50) \
        .generations(50) \
        .repeat_structure(171, 12, 10) \
        .init_uniform(cv.Nucleotide.G()) \
        .mutation_rate(0.0001) \
        .build()
    
    print(f"Created simulation: {sim}")
    print("  Initialization: Uniform with base G")
    print("  All positions start as G instead of A\n")
    
    return sim


def example_7_from_fasta():
    """Example 7: Initialize from FASTA file."""
    print("=" * 70)
    print("Example 7: Initialize from FASTA (Example)")
    print("=" * 70)
    
    # This is how you would do it (requires actual FASTA file)
    print("Example code (requires actual FASTA file):")
    print("""
    sim = cv.SimulationBuilder() \\
        .population_size(100) \\
        .generations(50) \\
        .init_from_fasta("sequences.fasta") \\
        .mutation_rate(0.0001) \\
        .build()
    """)
    print("This would:")
    print("  - Load sequences from FASTA file")
    print("  - Infer repeat structure from sequence length")
    print("  - Use provided population size and evolutionary parameters\n")


def example_8_compare_approaches():
    """Example 8: Compare builder vs direct construction."""
    print("=" * 70)
    print("Example 8: Builder vs Direct Construction")
    print("=" * 70)
    
    print("Old approach (verbose):")
    print("""
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
    fitness = cv.FitnessConfig.neutral()
    config = cv.SimulationConfig(
        population_size=100,
        total_generations=50,
        seed=42
    )
    sim = cv.Simulation(structure, mutation, recombination, fitness, config)
    """)
    
    print("\nNew approach (builder pattern):")
    print("""
    sim = cv.SimulationBuilder() \\
        .population_size(100) \\
        .generations(50) \\
        .repeat_structure(171, 12, 10) \\
        .mutation_rate(0.0001) \\
        .recombination(0.01, 0.7, 0.1) \\
        .seed(42) \\
        .build()
    """)
    
    print("\nBoth approaches work! Builder is more concise.\n")


def example_9_run_simulation():
    """Example 9: Build and run a simulation with fitness."""
    print("=" * 70)
    print("Example 9: Build, Run, and Analyze")
    print("=" * 70)
    
    # Build with GC content selection
    fitness = cv.FitnessConfig.with_gc_content(0.5, 2.0).build()
    
    sim = cv.SimulationBuilder() \
        .population_size(100) \
        .generations(50) \
        .repeat_structure(171, 12, 10) \
        .mutation_rate(0.0001) \
        .recombination(0.01, 0.7, 0.1) \
        .fitness(fitness) \
        .seed(42) \
        .build()
    
    print(f"Built simulation: {sim}")
    print(f"Fitness: {fitness}\n")
    
    print("Running simulation with GC content selection...")
    print(f"{'Gen':>4} {'Mean GC':>10}")
    print("-" * 16)
    
    for gen in range(0, 51, 10):
        if gen > 0:
            sim.run_for(10)
        
        population = sim.population()
        gc = cv.gc_content(population, individual_idx=None, haplotype_idx=None, chromosome_idx=0)
        print(f"{gen:4d} {gc:10.4f}")
    
    print("\nGC content evolves toward the optimum (0.5)!\n")


def example_10_neutral_vs_selection():
    """Example 10: Compare neutral vs selection using builder."""
    print("=" * 70)
    print("Example 10: Neutral vs Selection Comparison")
    print("=" * 70)
    
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
    
    print("Running both simulations...")
    sim_neutral.run()
    sim_selection.run()
    
    # Compare results
    gc_neutral = cv.gc_content(sim_neutral.population(), None, None, 0)
    gc_selection = cv.gc_content(sim_selection.population(), None, None, 0)
    
    print(f"\nAfter 30 generations:")
    print(f"  Neutral:   GC = {gc_neutral:.4f}")
    print(f"  Selection: GC = {gc_selection:.4f} (target: 0.5)")
    print(f"\nDifference from target:")
    print(f"  Neutral:   {abs(gc_neutral - 0.5):.4f}")
    print(f"  Selection: {abs(gc_selection - 0.5):.4f}")
    print("\nSelection makes a difference!\n")


def main():
    """Run all examples."""
    print("\n" + "=" * 70)
    print("Centrevo Python API: SimulationBuilder Pattern")
    print("=" * 70 + "\n")
    
    # Basic usage
    example_1_minimal_builder()
    example_2_with_evolution()
    
    # With fitness
    example_3_with_fitness()
    example_4_combined_fitness()
    
    # Initialization options
    example_5_random_initialization()
    example_6_custom_base()
    example_7_from_fasta()
    
    # Comparisons
    example_8_compare_approaches()
    example_9_run_simulation()
    example_10_neutral_vs_selection()
    
    print("=" * 70)
    print("All builder examples completed successfully!")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    main()
