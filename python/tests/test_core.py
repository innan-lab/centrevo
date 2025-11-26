"""
Test suite for Centrevo Python bindings - Core functionality.

Tests the basic classes and population creation.
"""

import pytest
import centrevo


class TestNucleotide:
    """Tests for Nucleotide class."""

    def test_create_nucleotides(self):
        """Test creating nucleotides from string."""
        a = centrevo.Nucleotide("A")
        c = centrevo.Nucleotide("C")
        g = centrevo.Nucleotide("G")
        t = centrevo.Nucleotide("T")

        assert str(a) == "A"
        assert str(c) == "C"
        assert str(g) == "G"
        assert str(t) == "T"

    def test_static_constructors(self):
        """Test static constructor methods."""
        a = centrevo.Nucleotide.A()
        c = centrevo.Nucleotide.C()
        g = centrevo.Nucleotide.G()
        t = centrevo.Nucleotide.T()

        assert str(a) == "A"
        assert str(c) == "C"
        assert str(g) == "G"
        assert str(t) == "T"

    def test_repr(self):
        """Test repr method."""
        a = centrevo.Nucleotide.A()
        assert "Nucleotide" in repr(a)
        assert "A" in repr(a)

    def test_invalid_nucleotide(self):
        """Test error handling for invalid nucleotides."""
        with pytest.raises(ValueError):
            centrevo.Nucleotide("X")

        with pytest.raises(ValueError):
            centrevo.Nucleotide("AA")

        with pytest.raises(ValueError):
            centrevo.Nucleotide("")


class TestAlphabet:
    """Tests for Alphabet class."""

    def test_dna_alphabet(self):
        """Test DNA alphabet creation."""
        alphabet = centrevo.Alphabet.dna()
        assert len(alphabet) == 4
        assert "Alphabet" in repr(alphabet)

    def test_custom_alphabet(self):
        """Test custom alphabet creation."""
        alphabet = centrevo.Alphabet(["A", "C"])
        assert len(alphabet) == 2


class TestChromosome:
    """Tests for Chromosome class."""

    def test_uniform_chromosome(self):
        """Test creating uniform chromosome."""
        alphabet = centrevo.Alphabet.dna()
        base = centrevo.Nucleotide.A()

        chr = centrevo.Chromosome.uniform(
            id="chr1",
            base=base,
            length=1000,
            ru_length=100,
            rus_per_hor=5,
            alphabet=alphabet
        )

        assert chr.id == "chr1"
        assert len(chr) == 1000
        assert "Chromosome" in repr(chr)

    def test_gc_content(self):
        """Test GC content calculation."""
        alphabet = centrevo.Alphabet.dna()
        base = centrevo.Nucleotide.A()

        chr = centrevo.Chromosome.uniform(
            id="chr1",
            base=base,
            length=100,
            ru_length=10,
            rus_per_hor=5,
            alphabet=alphabet
        )

        gc = chr.gc_content()
        assert 0.0 <= gc <= 1.0

    def test_to_string(self):
        """Test chromosome sequence as string."""
        alphabet = centrevo.Alphabet.dna()
        base = centrevo.Nucleotide.A()

        chr = centrevo.Chromosome.uniform(
            id="chr1",
            base=base,
            length=50,
            ru_length=10,
            rus_per_hor=5,
            alphabet=alphabet
        )

        seq = chr.to_string()
        assert len(seq) == 50
        assert all(c in "ACGT" for c in seq)


class TestHaplotype:
    """Tests for Haplotype class."""

    def test_create_empty_haplotype(self):
        """Test creating empty haplotype."""
        hap = centrevo.Haplotype()
        assert len(hap) == 0

    def test_from_chromosomes(self):
        """Test creating haplotype from chromosomes."""
        base = centrevo.Nucleotide.A()

        chr1 = centrevo.Chromosome.uniform(
            id="chr1", base=base, length=100,
            ru_length=10, rus_per_hor=5
        )
        chr2 = centrevo.Chromosome.uniform(
            id="chr2", base=base, length=100,
            ru_length=10, rus_per_hor=5
        )

        hap = centrevo.Haplotype.from_chromosomes([chr1, chr2])
        assert len(hap) == 2
        assert "Haplotype" in repr(hap)

    def test_haplotype_fitness_property(self):
        base = centrevo.Nucleotide.A()

        chrom = centrevo.Chromosome.uniform(
            id="chr1", base=base, length=100,
            ru_length=10, rus_per_hor=5
        )

        hap = centrevo.Haplotype.from_chromosomes([chrom])
        # Default is None
        assert hap.cached_fitness is None
        # Set and get cached fitness
        hap.cached_fitness = 0.25
        assert hap.cached_fitness == 0.25
        # Clear cached fitness
        hap.clear_cached_fitness()
        assert hap.cached_fitness is None


class TestIndividual:
    """Tests for Individual class."""

    def test_create_individual(self):
        """Test creating an individual."""
        alphabet = centrevo.Alphabet.dna()
        base = centrevo.Nucleotide.A()

        chr = centrevo.Chromosome.uniform(
            id="chr1", base=base, length=100,
            ru_length=10, rus_per_hor=5
        )

        hap1 = centrevo.Haplotype.from_chromosomes([chr])
        hap2 = centrevo.Haplotype.from_chromosomes([chr])

        ind = centrevo.Individual("ind1", hap1, hap2)
        assert ind.id == "ind1"
        assert "Individual" in repr(ind)

    def test_fitness_property(self):
        """Test fitness getter and setter."""
        alphabet = centrevo.Alphabet.dna()
        base = centrevo.Nucleotide.A()

        chrom = centrevo.Chromosome.uniform(
            id="chr1", base=base, length=100,
            ru_length=10, rus_per_hor=5
        )

        hap = centrevo.Haplotype.from_chromosomes([chrom])
        ind = centrevo.Individual("ind1", hap, hap)

        # Default fitness is None (not yet computed)
        assert ind.cached_fitness is None

        # Set and get cached fitness
        ind.cached_fitness = 0.5
        assert ind.cached_fitness == 0.5

class TestPopulation:
    """Tests for Population class."""

    def test_create_population(self):
        """Test creating a population."""
        alphabet = centrevo.Alphabet.dna()
        base = centrevo.Nucleotide.A()

        chr = centrevo.Chromosome.uniform(
            id="chr1", base=base, length=100,
            ru_length=10, rus_per_hor=5
        )

        hap = centrevo.Haplotype.from_chromosomes([chr])

        individuals = [
            centrevo.Individual(f"ind{i}", hap, hap)
            for i in range(5)
        ]

        pop = centrevo.Population("pop1", individuals)
        assert pop.id == "pop1"
        assert pop.size() == 5
        assert len(pop) == 5
        assert pop.generation() == 0
        assert "Population" in repr(pop)


class TestRepeatStructure:
    """Tests for RepeatStructure class."""

    def test_create_structure(self):
        """Test creating repeat structure."""
        alphabet = centrevo.Alphabet.dna()
        base = centrevo.Nucleotide.A()

        structure = centrevo.RepeatStructure(
            alphabet=alphabet,
            init_base=base,
            ru_length=171,
            rus_per_hor=12,
            hors_per_chr=100,
            chrs_per_hap=1
        )

        assert structure.ru_length == 171
        assert structure.rus_per_hor == 12
        assert structure.hors_per_chr == 100
        assert structure.chr_length() == 171 * 12 * 100
        assert "RepeatStructure" in repr(structure)

    def test_structure_properties(self):
        """Test structure property access."""
        alphabet = centrevo.Alphabet.dna()
        base = centrevo.Nucleotide.A()

        structure = centrevo.RepeatStructure(
            alphabet=alphabet,
            init_base=base,
            ru_length=100,
            rus_per_hor=5,
            hors_per_chr=10,
            chrs_per_hap=2
        )

        assert structure.ru_length == 100
        assert structure.rus_per_hor == 5
        assert structure.hors_per_chr == 10


class TestSimulationConfig:
    """Tests for SimulationConfig class."""

    def test_create_config_with_seed(self):
        """Test creating config with seed (required parameter)."""
        config = centrevo.SimulationConfig(
            population_size=100,
            total_generations=1000,
            seed=42
        )

        assert config.population_size == 100
        assert config.total_generations == 1000
        assert config.seed == 42
        assert "SimulationConfig" in repr(config)

    def test_create_config_different_seed(self):
        """Test creating config with different seed."""
        config = centrevo.SimulationConfig(
            population_size=50,
            total_generations=500,
            seed=12345
        )

        assert config.population_size == 50
        assert config.total_generations == 500
        assert config.seed == 12345


class TestCreatePopulation:
    """Tests for create_initial_population function."""

    def test_create_small_population(self):
        """Test creating a small population."""
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

        pop = centrevo.create_initial_population(size=10, structure=structure)

        assert pop.size() == 10
        assert len(pop) == 10
        assert pop.generation() == 0

    def test_create_larger_population(self):
        """Test creating a larger population."""
        alphabet = centrevo.Alphabet.dna()
        base = centrevo.Nucleotide.A()

        structure = centrevo.RepeatStructure(
            alphabet=alphabet,
            init_base=base,
            ru_length=171,
            rus_per_hor=12,
            hors_per_chr=50,
            chrs_per_hap=1
        )

        pop = centrevo.create_initial_population(size=100, structure=structure)

        assert pop.size() == 100
        assert pop.generation() == 0

    def test_population_reproducibility(self):
        """Test that populations are created consistently."""
        alphabet = centrevo.Alphabet.dna()
        base = centrevo.Nucleotide.A()

        structure = centrevo.RepeatStructure(
            alphabet=alphabet,
            init_base=base,
            ru_length=50,
            rus_per_hor=5,
            hors_per_chr=10,
            chrs_per_hap=1
        )

        pop1 = centrevo.create_initial_population(size=5, structure=structure)
        pop2 = centrevo.create_initial_population(size=5, structure=structure)

        # Both should have same size
        assert pop1.size() == pop2.size() == 5


class TestRecordingStrategy:
    """Tests for RecordingStrategy class."""

    def test_every_n_strategy(self):
        """Test every N generations strategy."""
        strategy = centrevo.RecordingStrategy.every_n(100)
        assert "RecordingStrategy" in repr(strategy)
        assert "100" in repr(strategy)

    def test_specific_generations_strategy(self):
        """Test specific generations strategy."""
        strategy = centrevo.RecordingStrategy.specific([0, 10, 50, 100])
        assert "RecordingStrategy" in repr(strategy)

    def test_all_strategy(self):
        """Test record all generations strategy."""
        strategy = centrevo.RecordingStrategy.all()
        assert "RecordingStrategy" in repr(strategy)

    def test_none_strategy(self):
        """Test record no generations strategy."""
        strategy = centrevo.RecordingStrategy.none()
        assert "RecordingStrategy" in repr(strategy)
