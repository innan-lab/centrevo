#!/usr/bin/env python3
"""
Example: Using Custom Sequences and Resume Functionality in Python

This example demonstrates how to:
1. Create simulations with custom sequences from FASTA/JSON
2. Resume simulations from checkpoints
3. Continue simulations from previous runs
"""

import centrevo as cv
import json
import tempfile
import os

def example_1_fasta_initialization():
    """Example 1: Initialize simulation from FASTA file."""
    print("=" * 60)
    print("Example 1: Initialize from FASTA file")
    print("=" * 60)

    # Create a temporary FASTA file
    fasta_content = """
>ind0_h1_chr0
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>ind0_h2_chr0
TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA
>ind1_h1_chr0
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>ind1_h2_chr0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
""".strip()

    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(fasta_content)
        fasta_path = f.name

    try:
        # Define simulation parameters
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
        config = cv.SimulationConfig(population_size=2, total_generations=10, seed=42)

        # Initialize from FASTA
        sim = cv.Simulation.from_sequences(
            source_type="fasta",
            source_path=fasta_path,
            structure=structure,
            mutation=mutation,
            recombination=recombination,
            fitness=fitness,
            config=config,
        )

        print(f"Created simulation: {sim}")
        print(f"Population size: {sim.population().size()}")
        print(f"Generation: {sim.generation()}")

        # Run for a few generations
        sim.run_for(5)
        print(f"After running 5 generations: {sim}")

    finally:
        os.unlink(fasta_path)

def example_2_json_initialization():
    """Example 2: Initialize simulation from JSON."""
    print("\n" + "=" * 60)
    print("Example 2: Initialize from JSON")
    print("=" * 60)

    # Create JSON with sequences (100 bases each)
    sequences = [
        {"id": "ind0_h1", "seq": "ACGT" * 25},
        {"id": "ind0_h2", "seq": "TGCA" * 25},
        {"id": "ind1_h1", "seq": "AAAA" * 25},
        {"id": "ind1_h2", "seq": "CCCC" * 25},
    ]
    json_str = json.dumps(sequences)

    # Define simulation parameters
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
    config = cv.SimulationConfig(population_size=2, total_generations=10, seed=123)

    # Initialize from JSON string
    sim = cv.Simulation.from_sequences(
        source_type="json",
        source_path=json_str,
        structure=structure,
        mutation=mutation,
        recombination=recombination,
        fitness=fitness,
        config=config,
    )

    print(f"Created simulation: {sim}")
    print(f"Population size: {sim.population().size()}")

    # Run simulation
    sim.run_for(3)
    print(f"After running 3 generations: {sim}")

def example_3_checkpoint_resume():
    """Example 3: Resume simulation from checkpoint."""
    print("\n" + "=" * 60)
    print("Example 3: Resume from checkpoint")
    print("=" * 60)

    with tempfile.NamedTemporaryFile(suffix='.db', delete=False) as f:
        db_path = f.name

    try:
        # Create and run initial simulation
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
        config = cv.SimulationConfig(population_size=5, total_generations=100, seed=42)

        # Create simulation
        sim = cv.Simulation(structure, mutation, recombination, fitness, config)
        print(f"Created initial simulation: {sim}")

        # Create recorder with full config
        recorder = cv.Recorder(db_path, "checkpoint_example", cv.RecordingStrategy.every_n(5))
        recorder.record_full_config(structure, mutation, recombination, fitness, config)

        # Run and record checkpoints
        for gen_batch in range(5):
            sim.step()
            generation = sim.generation()

            # Record population state
            recorder.record_generation(sim.population(), generation)

            # Record checkpoint every 5 generations
            if generation % 5 == 0:
                recorder.record_checkpoint(sim, generation)
                print(f"Checkpoint saved at generation {generation}")

        print(f"Simulation state: {sim}")

        # Now resume from checkpoint
        print("\nResuming from checkpoint...")
        resumed_sim = cv.Simulation.from_checkpoint(db_path, "checkpoint_example")
        print(f"Resumed simulation: {resumed_sim}")

        # Continue running
        resumed_sim.run_for(5)
        print(f"After continuing 5 more generations: {resumed_sim}")

    finally:
        if os.path.exists(db_path):
            os.unlink(db_path)

def example_4_database_initialization():
    """Example 4: Start new simulation from previous run's state."""
    print("\n" + "=" * 60)
    print("Example 4: Initialize from database")
    print("=" * 60)

    with tempfile.NamedTemporaryFile(suffix='.db', delete=False) as f:
        db_path = f.name

    try:
        # Run initial simulation
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
        config = cv.SimulationConfig(population_size=3, total_generations=20, seed=100)

        # Create and run initial simulation
        sim1 = cv.Simulation(structure, mutation, recombination, fitness, config)
        print(f"Initial simulation: {sim1}")

        # Record generations
        recorder = cv.Recorder(db_path, "original", cv.RecordingStrategy.all())
        recorder.record_full_config(structure, mutation, recombination, fitness, config)

        for i in range(10):
            sim1.step()
            recorder.record_generation(sim1.population(), sim1.generation())

        print(f"After 10 generations: {sim1}")

        # Start a new simulation from generation 5 of the previous run
        print("\nStarting new simulation from generation 5...")
        new_config = cv.SimulationConfig(population_size=3, total_generations=10, seed=999)

        sim2 = cv.Simulation.from_sequences(
            source_type="database",
            source_path=db_path,
            structure=structure,
            mutation=mutation,
            recombination=recombination,
            fitness=fitness,
            config=new_config,
            sim_id="original",
            generation=5,
        )

        print(f"New simulation from generation 5: {sim2}")
        print(f"This starts at generation 0 with the population state from gen 5 of original")

        # Run the new simulation
        sim2.run_for(5)
        print(f"After running 5 more generations: {sim2}")

    finally:
        if os.path.exists(db_path):
            os.unlink(db_path)

def main():
    """Run all examples."""
    print("\n" + "=" * 60)
    print("Custom Sequences and Resume Examples")
    print("=" * 60 + "\n")

    example_1_fasta_initialization()
    example_2_json_initialization()
    example_3_checkpoint_resume()
    example_4_database_initialization()

    print("\n" + "=" * 60)
    print("All examples completed successfully!")
    print("=" * 60)

if __name__ == "__main__":
    main()
