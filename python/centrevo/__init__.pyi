"""
Type stubs for Centrevo Python Package

Provides type hints for the compiled Rust extension module.
"""

from __future__ import annotations
from typing import Dict, List, Optional, Any, Tuple

__version__: str

# Core classes
class Nucleotide:
    """Nucleotide base (A, C, G, T)."""
    def __init__(self, base: str) -> None: ...
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...
    @staticmethod
    def A() -> Nucleotide: ...
    @staticmethod
    def C() -> Nucleotide: ...
    @staticmethod
    def G() -> Nucleotide: ...
    @staticmethod
    def T() -> Nucleotide: ...

class Alphabet:
    """Alphabet of nucleotide bases."""
    def __init__(self, chars: List[str]) -> None: ...
    @staticmethod
    def dna() -> Alphabet: ...
    def __len__(self) -> int: ...
    def __repr__(self) -> str: ...

class Chromosome:
    """Chromosome with repeat structure."""
    @staticmethod
    def uniform(
        id: str,
        base: Nucleotide,
        length: int,
        ru_length: int,
        rus_per_hor: int,
        alphabet: Alphabet,
    ) -> Chromosome: ...
    @property
    def id(self) -> str: ...
    def __len__(self) -> int: ...
    def to_string(self) -> str: ...
    def gc_content(self) -> float: ...
    def __repr__(self) -> str: ...

class Haplotype:
    """Haplotype (collection of chromosomes)."""
    def __init__(self) -> None: ...
    @staticmethod
    def from_chromosomes(chromosomes: List[Chromosome]) -> Haplotype: ...
    def __len__(self) -> int: ...
    def __repr__(self) -> str: ...

class Individual:
    """Individual organism with diploid genome."""
    def __init__(self, id: str, haplotype1: Haplotype, haplotype2: Haplotype) -> None: ...
    @property
    def id(self) -> str: ...
    @property
    def fitness(self) -> float: ...
    @fitness.setter
    def fitness(self, value: float) -> None: ...
    def __repr__(self) -> str: ...

class Population:
    """Population of individuals."""
    def __init__(self, id: str, individuals: List[Individual]) -> None: ...
    @property
    def id(self) -> str: ...
    def size(self) -> int: ...
    def generation(self) -> int: ...
    def __len__(self) -> int: ...
    def __repr__(self) -> str: ...

class RepeatStructure:
    """Repeat structure configuration."""
    def __init__(
        self,
        alphabet: Alphabet,
        init_base: Nucleotide,
        ru_length: int,
        rus_per_hor: int,
        hors_per_chr: int,
        chrs_per_hap: int,
    ) -> None: ...
    @property
    def ru_length(self) -> int: ...
    @property
    def rus_per_hor(self) -> int: ...
    @property
    def hors_per_chr(self) -> int: ...
    def chr_length(self) -> int: ...
    def __repr__(self) -> str: ...

class SimulationConfig:
    """Simulation configuration."""
    def __init__(self, population_size: int, total_generations: int, seed: Optional[int] = None) -> None: ...
    @property
    def population_size(self) -> int: ...
    @property
    def total_generations(self) -> int: ...
    @property
    def seed(self) -> Optional[int]: ...
    def __repr__(self) -> str: ...

class RecordingStrategy:
    """Recording strategy for simulation snapshots."""
    @staticmethod
    def every_n(n: int) -> RecordingStrategy: ...
    @staticmethod
    def specific(generations: List[int]) -> RecordingStrategy: ...
    @staticmethod
    def all() -> RecordingStrategy: ...
    @staticmethod
    def none() -> RecordingStrategy: ...
    def __repr__(self) -> str: ...

class Recorder:
    """Simulation recorder."""
    def __init__(self, db_path: str, sim_id: str, strategy: RecordingStrategy) -> None: ...
    def record_metadata(self, config: SimulationConfig) -> None: ...
    def record_generation(self, population: Population, generation: int) -> None: ...
    def __repr__(self) -> str: ...

class QueryBuilder:
    """Query builder for simulation results."""
    def __init__(self, db_path: str) -> None: ...
    def list_simulations(self) -> List[str]: ...
    def get_recorded_generations(self, sim_id: str) -> List[int]: ...
    def __repr__(self) -> str: ...

# Core functions
def create_initial_population(size: int, structure: RepeatStructure) -> Population:
    """Create an initial population with uniform sequences."""
    ...

# Analysis functions - Diversity metrics
def nucleotide_diversity(population: Population, chromosome_idx: int) -> float:
    """Calculate nucleotide diversity (π) for a population."""
    ...

def tajimas_d(population: Population, chromosome_idx: int) -> float:
    """Calculate Tajima's D statistic."""
    ...

def wattersons_theta(population: Population, chromosome_idx: int) -> float:
    """Calculate Watterson's estimator (θ_W)."""
    ...

def haplotype_diversity(population: Population, chromosome_idx: int) -> float:
    """Calculate haplotype diversity."""
    ...

# Analysis functions - Linkage
def linkage_disequilibrium(
    population: Population,
    pos1: int,
    pos2: int,
    chromosome_idx: int,
    haplotype_idx: int,
) -> Optional[Dict[str, float]]:
    """Calculate linkage disequilibrium statistics between two sites.
    
    Returns:
        Dictionary with keys 'D', 'D_prime', 'r_squared', or None if calculation fails
    """
    ...

def ld_decay(
    population: Population,
    chromosome_idx: int,
    haplotype_idx: int,
    max_distance: int,
    bin_size: int,
) -> Dict[str, List[float]]:
    """Calculate LD decay across a range of distances.
    
    Returns:
        Dictionary with 'distances' and 'r_squared_values' as lists
    """
    ...

def haplotype_blocks(
    population: Population,
    chromosome_idx: int,
    haplotype_idx: int,
    r_squared_threshold: Optional[float] = None,
) -> List[Tuple[int, int]]:
    """Identify haplotype blocks based on LD threshold.
    
    Returns:
        List of blocks, each block is a tuple (start, end)
    """
    ...

# Analysis functions - Distance
def pairwise_distances(population: Population, chromosome_idx: int) -> List[float]:
    """Calculate pairwise genetic distances between all sequences."""
    ...

def distance_matrix(population: Population, chromosome_idx: int) -> List[List[float]]:
    """Calculate full distance matrix between all sequences.
    
    Returns:
        2D list representing the distance matrix (2n × 2n)
    """
    ...

# Analysis functions - Composition
def gc_content(
    population: Population,
    individual_idx: Optional[int],
    haplotype_idx: Optional[int],
    chromosome_idx: Optional[int],
) -> float:
    """Calculate GC content at different levels.
    
    Args:
        population: Population to analyze
        individual_idx: Optional individual index
        haplotype_idx: Optional haplotype index (0 or 1)
        chromosome_idx: Optional chromosome index
        
    Returns:
        GC content (0.0 to 1.0)
    """
    ...

def nucleotide_composition(
    population: Population,
    individual_idx: Optional[int],
    haplotype_idx: Optional[int],
    chromosome_idx: Optional[int],
) -> Dict[str, float]:
    """Calculate nucleotide composition at different levels.
    
    Returns:
        Dictionary mapping nucleotide (str) to frequency (float)
    """
    ...

# Analysis functions - Polymorphism
def count_segregating_sites(
    population: Population,
    chromosome_idx: int,
    haplotype_idx: int,
) -> int:
    """Count number of segregating sites in a population."""
    ...

# Export functions for PyArrow integration
def export_diversity_metrics(
    population: Population,
    chromosome_idx: int,
) -> Dict[str, Any]:
    """Export diversity metrics to PyArrow table format.
    
    Returns:
        Dictionary with metric names as keys and values, including:
        - nucleotide_diversity
        - tajimas_d
        - wattersons_theta
        - haplotype_diversity
        - generation
        - population_size
    """
    ...

def export_distance_matrix(
    population: Population,
    chromosome_idx: int,
) -> List[Dict[str, Any]]:
    """Export distance matrix to flat format for PyArrow.
    
    Returns:
        List of dicts with keys 'sequence_i', 'sequence_j', 'distance'
    """
    ...

def export_ld_decay(
    population: Population,
    chromosome_idx: int,
    haplotype_idx: int,
    max_distance: int,
    bin_size: int,
) -> List[Dict[str, float]]:
    """Export LD decay data to PyArrow-friendly format.
    
    Returns:
        List of dicts with keys 'distance', 'r_squared'
    """
    ...
