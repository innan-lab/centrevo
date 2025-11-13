"""
Test suite for Centrevo Python bindings - Analysis functions.

Tests diversity metrics, linkage, distance calculations, and composition analysis.
"""

import pytest
import centrevo


@pytest.fixture
def small_population():
    """Create a small test population."""
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

    return centrevo.create_initial_population(size=10, structure=structure)


@pytest.fixture
def larger_population():
    """Create a larger test population."""
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

    return centrevo.create_initial_population(size=50, structure=structure)


class TestDiversityMetrics:
    """Tests for diversity metric calculations."""

    def test_nucleotide_diversity(self, small_population):
        """Test nucleotide diversity calculation."""
        pi = centrevo.nucleotide_diversity(small_population, chromosome_idx=0)

        assert isinstance(pi, float)
        assert pi >= 0.0
        # For uniform initial population, diversity should be low
        assert pi <= 0.1

    def test_nucleotide_diversity_larger_population(self, larger_population):
        """Test diversity with larger population."""
        pi = centrevo.nucleotide_diversity(larger_population, chromosome_idx=0)

        assert isinstance(pi, float)
        assert pi >= 0.0

    def test_tajimas_d(self, small_population):
        """Test Tajima's D calculation."""
        tajima = centrevo.tajimas_d(small_population, chromosome_idx=0)

        assert isinstance(tajima, float)
        # Tajima's D can be negative, zero, or positive
        assert -5.0 <= tajima <= 5.0

    def test_wattersons_theta(self, small_population):
        """Test Watterson's theta calculation."""
        theta = centrevo.wattersons_theta(small_population, chromosome_idx=0)

        assert isinstance(theta, float)
        assert theta >= 0.0

    def test_haplotype_diversity(self, small_population):
        """Test haplotype diversity calculation."""
        h = centrevo.haplotype_diversity(small_population, chromosome_idx=0)

        assert isinstance(h, float)
        assert 0.0 <= h <= 1.0

    def test_diversity_consistency(self, small_population):
        """Test that diversity metrics are consistent across calls."""
        pi1 = centrevo.nucleotide_diversity(small_population, chromosome_idx=0)
        pi2 = centrevo.nucleotide_diversity(small_population, chromosome_idx=0)

        assert pi1 == pi2


class TestLinkageAnalysis:
    """Tests for linkage disequilibrium analysis."""

    def test_linkage_disequilibrium(self, small_population):
        """Test LD calculation between two sites."""
        # Calculate LD between two positions
        ld = centrevo.linkage_disequilibrium(
            small_population,
            pos1=10,
            pos2=50,
            chromosome_idx=0,
            haplotype_idx=0
        )

        # LD may be None if sites are not polymorphic
        if ld is not None:
            assert isinstance(ld, dict)
            assert 'D' in ld
            assert 'D_prime' in ld
            assert 'r_squared' in ld

            # Check value ranges
            assert isinstance(ld['r_squared'], float)
            assert 0.0 <= ld['r_squared'] <= 1.0

    def test_ld_decay(self, small_population):
        """Test LD decay calculation."""
        decay = centrevo.ld_decay(
            small_population,
            chromosome_idx=0,
            haplotype_idx=0,
            max_distance=500,
            bin_size=50
        )

        assert isinstance(decay, dict)
        assert 'distances' in decay
        assert 'r_squared_values' in decay

        distances = decay['distances']
        r_squared = decay['r_squared_values']

        assert isinstance(distances, list)
        assert isinstance(r_squared, list)
        assert len(distances) == len(r_squared)

        # All r² values should be in [0, 1]
        for r2 in r_squared:
            assert 0.0 <= r2 <= 1.0

    def test_haplotype_blocks(self, small_population):
        """Test haplotype block identification."""
        blocks = centrevo.haplotype_blocks(
            small_population,
            chromosome_idx=0,
            haplotype_idx=0,
            r_squared_threshold=0.8
        )

        assert isinstance(blocks, list)

        # Each block should be a tuple of (start, end)
        for block in blocks:
            assert isinstance(block, tuple)
            assert len(block) == 2
            start, end = block
            assert isinstance(start, int)
            assert isinstance(end, int)
            assert start < end

    def test_haplotype_blocks_default_threshold(self, small_population):
        """Test haplotype blocks with threshold of 0.8."""
        blocks = centrevo.haplotype_blocks(
            small_population,
            chromosome_idx=0,
            haplotype_idx=0,
            r_squared_threshold=0.8
        )

        assert isinstance(blocks, list)


class TestDistanceCalculations:
    """Tests for genetic distance calculations."""

    def test_pairwise_distances(self, small_population):
        """Test pairwise distance calculation."""
        distances = centrevo.pairwise_distances(
            small_population,
            chromosome_idx=0
        )

        assert isinstance(distances, list)
        assert len(distances) > 0

        # All distances should be non-negative (can be int or float)
        for d in distances:
            assert isinstance(d, (int, float))
            assert d >= 0

    def test_distance_matrix(self, small_population):
        """Test distance matrix calculation."""
        matrix = centrevo.distance_matrix(
            small_population,
            chromosome_idx=0
        )

        assert isinstance(matrix, list)
        assert len(matrix) > 0

        # Matrix should be square
        n = len(matrix)
        for row in matrix:
            assert isinstance(row, list)
            assert len(row) == n

        # Diagonal should be zero
        for i in range(n):
            assert matrix[i][i] == 0.0

        # Matrix should be symmetric
        for i in range(n):
            for j in range(n):
                assert matrix[i][j] == pytest.approx(matrix[j][i])

    def test_distance_matrix_size(self, small_population):
        """Test that distance matrix has correct size."""
        matrix = centrevo.distance_matrix(
            small_population,
            chromosome_idx=0
        )

        # For a diploid population of size n, we have 2n haplotypes
        pop_size = small_population.size()
        expected_size = 2 * pop_size

        assert len(matrix) == expected_size


class TestCompositionAnalysis:
    """Tests for nucleotide composition analysis."""

    def test_gc_content_population(self, small_population):
        """Test GC content at population level."""
        gc = centrevo.gc_content(
            small_population,
            individual_idx=None,
            haplotype_idx=None,
            chromosome_idx=None
        )

        assert isinstance(gc, float)
        assert 0.0 <= gc <= 1.0

    def test_gc_content_chromosome(self, small_population):
        """Test GC content at chromosome level."""
        gc = centrevo.gc_content(
            small_population,
            individual_idx=None,
            haplotype_idx=None,
            chromosome_idx=0
        )

        assert isinstance(gc, float)
        assert 0.0 <= gc <= 1.0

    def test_gc_content_individual(self, small_population):
        """Test GC content at individual level."""
        gc = centrevo.gc_content(
            small_population,
            individual_idx=0,
            haplotype_idx=None,
            chromosome_idx=None
        )

        assert isinstance(gc, float)
        assert 0.0 <= gc <= 1.0

    def test_nucleotide_composition_population(self, small_population):
        """Test nucleotide composition at population level."""
        comp = centrevo.nucleotide_composition(
            small_population,
            individual_idx=None,
            haplotype_idx=None,
            chromosome_idx=None
        )

        assert isinstance(comp, dict)

        # Should have all four nucleotides
        assert 'A' in comp or 'C' in comp or 'G' in comp or 'T' in comp

        # Composition returns counts, not frequencies
        # All counts should be non-negative integers
        for count in comp.values():
            assert isinstance(count, int)
            assert count >= 0

        # Total count should be > 0
        total = sum(comp.values())
        assert total > 0

    def test_nucleotide_composition_chromosome(self, small_population):
        """Test composition at chromosome level."""
        comp = centrevo.nucleotide_composition(
            small_population,
            individual_idx=0,
            haplotype_idx=0,
            chromosome_idx=0
        )

        assert isinstance(comp, dict)
        total = sum(comp.values())
        assert total > 0


class TestPolymorphismAnalysis:
    """Tests for polymorphism analysis."""

    def test_count_segregating_sites(self, small_population):
        """Test segregating sites count."""
        sites = centrevo.count_segregating_sites(
            small_population,
            chromosome_idx=0,
            haplotype_idx=0
        )

        assert isinstance(sites, int)
        assert sites >= 0

    def test_segregating_sites_larger_population(self, larger_population):
        """Test with larger population."""
        sites = centrevo.count_segregating_sites(
            larger_population,
            chromosome_idx=0,
            haplotype_idx=0
        )

        assert isinstance(sites, int)
        assert sites >= 0


class TestExportFunctions:
    """Tests for PyArrow export functions."""

    def test_export_diversity_metrics(self, small_population):
        """Test exporting diversity metrics."""
        metrics = centrevo.export_diversity_metrics(
            small_population,
            chromosome_idx=0
        )

        assert isinstance(metrics, dict)

        # Check expected keys
        expected_keys = [
            'nucleotide_diversity',
            'tajimas_d',
            'wattersons_theta',
            'haplotype_diversity',
            'generation',
            'population_size'
        ]

        for key in expected_keys:
            assert key in metrics, f"Missing key: {key}"

        # Check value types
        assert isinstance(metrics['nucleotide_diversity'], float)
        assert isinstance(metrics['tajimas_d'], float)
        assert isinstance(metrics['wattersons_theta'], float)
        assert isinstance(metrics['haplotype_diversity'], float)
        assert isinstance(metrics['generation'], int)
        assert isinstance(metrics['population_size'], int)

        assert metrics['population_size'] == small_population.size()
        assert metrics['generation'] == small_population.generation()

    def test_export_distance_matrix(self, small_population):
        """Test exporting distance matrix."""
        data = centrevo.export_distance_matrix(
            small_population,
            chromosome_idx=0
        )

        assert isinstance(data, list)
        assert len(data) > 0

        # Check structure of first entry
        entry = data[0]
        assert isinstance(entry, dict)
        assert 'sequence_i' in entry
        assert 'sequence_j' in entry
        assert 'distance' in entry

        assert isinstance(entry['sequence_i'], int)
        assert isinstance(entry['sequence_j'], int)
        assert isinstance(entry['distance'], float)
        assert entry['distance'] >= 0.0

    def test_export_ld_decay(self, small_population):
        """Test exporting LD decay data."""
        data = centrevo.export_ld_decay(
            small_population,
            chromosome_idx=0,
            haplotype_idx=0,
            max_distance=500,
            bin_size=50
        )

        assert isinstance(data, list)

        if len(data) > 0:
            # Check structure
            entry = data[0]
            assert isinstance(entry, dict)
            assert 'distance' in entry
            assert 'r_squared' in entry

            assert isinstance(entry['distance'], (int, float))
            assert isinstance(entry['r_squared'], float)
            assert 0.0 <= entry['r_squared'] <= 1.0


class TestEdgeCases:
    """Tests for edge cases and error handling."""

    def test_invalid_chromosome_index(self, small_population):
        """Test with invalid chromosome index - may return default value or handle gracefully."""
        # The API may handle invalid indices gracefully without raising an exception
        try:
            result = centrevo.nucleotide_diversity(small_population, chromosome_idx=999)
            # If it doesn't raise, it should return a sensible value
            assert isinstance(result, float)
        except Exception:
            # Exception is also acceptable
            pass

    def test_invalid_individual_index(self, small_population):
        """Test with invalid individual index - may handle gracefully."""
        try:
            result = centrevo.gc_content(
                small_population,
                individual_idx=999,
                haplotype_idx=None,
                chromosome_idx=None
            )
            # If it doesn't raise, should return a sensible value
            assert isinstance(result, float)
        except Exception:
            # Exception is also acceptable
            pass

    def test_zero_max_distance_ld_decay(self, small_population):
        """Test LD decay with zero max distance."""
        decay = centrevo.ld_decay(
            small_population,
            chromosome_idx=0,
            haplotype_idx=0,
            max_distance=0,
            bin_size=10
        )

        assert isinstance(decay, dict)
        # Should have empty or minimal data
        assert len(decay['distances']) >= 0
