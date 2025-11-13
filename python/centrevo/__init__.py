"""
Centrevo Python Package

High-performance simulator for centromeric sequence evolution.
"""

__version__ = "0.2.1"

# Import the compiled Rust module
# The binary module is named 'centrevo' via PyO3's #[pymodule]
try:
    from .centrevo import (
        # Core classes
        Nucleotide,
        Alphabet,
        Chromosome,
        Haplotype,
        Individual,
        Population,
        RepeatStructure,
        SimulationConfig,
        Simulation,
        SimulationBuilder,
        MutationConfig,
        RecombinationConfig,
        FitnessConfig,
        FitnessConfigBuilder,
        Recorder,
        QueryBuilder,
        RecordingStrategy,

        # Core functions
        create_initial_population,

        # Analysis functions - Diversity metrics
        nucleotide_diversity,
        tajimas_d,
        wattersons_theta,
        haplotype_diversity,

        # Analysis functions - Linkage
        linkage_disequilibrium,
        ld_decay,
        haplotype_blocks,

        # Analysis functions - Distance
        pairwise_distances,
        distance_matrix,

        # Analysis functions - Composition
        gc_content,
        nucleotide_composition,

        # Analysis functions - Polymorphism
        count_segregating_sites,

        # Export functions for PyArrow integration
        export_diversity_metrics,
        export_distance_matrix,
        export_ld_decay,
    )
except ImportError as e:
    raise ImportError(
        f"Failed to import centrevo binary module. "
        f"Make sure you've built the package with maturin: {e}"
    ) from e

# Import plotting utilities
from . import plotting

__all__ = [
    # Core classes
    'Nucleotide',
    'Alphabet',
    'Chromosome',
    'Haplotype',
    'Individual',
    'Population',
    'RepeatStructure',
    'SimulationConfig',
    'Simulation',
    'SimulationBuilder',
    'MutationConfig',
    'RecombinationConfig',
    'FitnessConfig',
    'FitnessConfigBuilder',
    'Recorder',
    'QueryBuilder',
    'RecordingStrategy',

    # Core functions
    'create_initial_population',

    # Analysis functions - Diversity metrics
    'nucleotide_diversity',
    'tajimas_d',
    'wattersons_theta',
    'haplotype_diversity',

    # Analysis functions - Linkage
    'linkage_disequilibrium',
    'ld_decay',
    'haplotype_blocks',

    # Analysis functions - Distance
    'pairwise_distances',
    'distance_matrix',

    # Analysis functions - Composition
    'gc_content',
    'nucleotide_composition',

    # Analysis functions - Polymorphism
    'count_segregating_sites',

    # Export functions for PyArrow integration
    'export_diversity_metrics',
    'export_distance_matrix',
    'export_ld_decay',

    # Submodules
    'plotting',
]
