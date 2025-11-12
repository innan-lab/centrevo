"""
Tests for custom sequence initialization and checkpoint resume in Python API.

Run with: python -m pytest python/tests/test_custom_sequences.py
"""

import centrevo as cv
import tempfile
import os
import json
import pytest


def test_uniform_simulation():
    """Test creating a basic simulation with uniform sequences."""
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
    
    sim = cv.Simulation(structure, mutation, recombination, fitness, config)
    
    assert sim.generation() == 0
    assert sim.population().size() == 2
    
    sim.step()
    assert sim.generation() == 1


def test_fasta_initialization():
    """Test initializing simulation from FASTA file."""
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
        
        sim = cv.Simulation.from_sequences(
            source_type="fasta",
            source_path=fasta_path,
            structure=structure,
            mutation=mutation,
            recombination=recombination,
            fitness=fitness,
            config=config,
        )
        
        assert sim.generation() == 0
        assert sim.population().size() == 2
        
        sim.run_for(3)
        assert sim.generation() == 3
        
    finally:
        os.unlink(fasta_path)


def test_json_initialization():
    """Test initializing simulation from JSON."""
    sequences = [
        {"id": "ind0_h1", "seq": "ACGT" * 25},
        {"id": "ind0_h2", "seq": "TGCA" * 25},
        {"id": "ind1_h1", "seq": "AAAA" * 25},
        {"id": "ind1_h2", "seq": "CCCC" * 25},
    ]
    json_str = json.dumps(sequences)
    
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
    
    sim = cv.Simulation.from_sequences(
        source_type="json",
        source_path=json_str,
        structure=structure,
        mutation=mutation,
        recombination=recombination,
        fitness=fitness,
        config=config,
    )
    
    assert sim.generation() == 0
    assert sim.population().size() == 2


def test_checkpoint_resume():
    """Test resuming simulation from checkpoint."""
    with tempfile.NamedTemporaryFile(suffix='.db', delete=False) as f:
        db_path = f.name
    
    try:
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
        config = cv.SimulationConfig(population_size=3, total_generations=100, seed=42)
        
        # Create and run simulation
        sim = cv.Simulation(structure, mutation, recombination, fitness, config)
        
        # Record with checkpoints
        recorder = cv.Recorder(db_path, "test_sim", cv.RecordingStrategy.all())
        recorder.record_full_config(structure, mutation, recombination, fitness, config)
        
        for _ in range(5):
            sim.step()
            recorder.record_generation(sim.population(), sim.generation())
        
        recorder.record_checkpoint(sim, sim.generation())
        
        assert sim.generation() == 5
        
        # Resume from checkpoint
        resumed_sim = cv.Simulation.from_checkpoint(db_path, "test_sim")
        
        assert resumed_sim.generation() == 5
        assert resumed_sim.population().size() == 3
        
        # Continue running
        resumed_sim.run_for(3)
        assert resumed_sim.generation() == 8
        
    finally:
        if os.path.exists(db_path):
            os.unlink(db_path)


def test_database_initialization():
    """Test initializing from previous simulation database."""
    with tempfile.NamedTemporaryFile(suffix='.db', delete=False) as f:
        db_path = f.name
    
    try:
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
        config = cv.SimulationConfig(population_size=2, total_generations=20, seed=100)
        
        # Run initial simulation
        sim1 = cv.Simulation(structure, mutation, recombination, fitness, config)
        
        recorder = cv.Recorder(db_path, "original", cv.RecordingStrategy.all())
        recorder.record_full_config(structure, mutation, recombination, fitness, config)
        
        for _ in range(10):
            sim1.step()
            recorder.record_generation(sim1.population(), sim1.generation())
        
        assert sim1.generation() == 10
        
        # Start new simulation from generation 5
        new_config = cv.SimulationConfig(population_size=2, total_generations=10, seed=999)
        
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
        
        assert sim2.generation() == 0  # New simulation starts at gen 0
        assert sim2.population().size() == 2
        
        sim2.run_for(5)
        assert sim2.generation() == 5
        
    finally:
        if os.path.exists(db_path):
            os.unlink(db_path)


def test_invalid_source_type():
    """Test that invalid source type raises error."""
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
    
    with pytest.raises(Exception):  # Should raise ValueError
        cv.Simulation.from_sequences(
            source_type="invalid",
            source_path="dummy",
            structure=structure,
            mutation=mutation,
            recombination=recombination,
            fitness=fitness,
            config=config,
        )
