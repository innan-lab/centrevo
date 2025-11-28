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
    def __init__(self, population_size: int, total_generations: int, seed: int) -> None: ...
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
    def record_full_config(
        self,
        structure: RepeatStructure,
        mutation: MutationConfig,
        recombination: RecombinationConfig,
        fitness: FitnessConfig,
        config: SimulationConfig,
    ) -> None:
        """Record complete simulation configuration for resumability."""
        ...
    def record_generation(self, population: Population, generation: int) -> None: ...
    def record_checkpoint(self, simulation: Simulation, generation: int) -> None:
        """Record a checkpoint for resuming the simulation."""
        ...
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

class MutationConfig:
    """Mutation configuration for simulations."""
    @staticmethod
    def uniform(alphabet: Alphabet, rate: float) -> MutationConfig:
        """Create a uniform mutation rate configuration."""
        ...
    def __repr__(self) -> str: ...

class RecombinationConfig:
    """Recombination configuration for simulations."""
    @staticmethod
    def standard(
        break_prob: float,
        crossover_prob: float,
        gc_extension_prob: float,
    ) -> RecombinationConfig:
        """Create a standard recombination configuration."""
        ...
    def __repr__(self) -> str: ...

class FitnessConfig:
    """Fitness configuration for simulations."""
    @staticmethod
    def neutral() -> FitnessConfig:
        """Create a neutral fitness configuration (no selection)."""
        ...

    @staticmethod
    def with_gc_content(optimum: float, concentration: float) -> FitnessConfigBuilder:
        """Start building with GC content fitness.

        Args:
            optimum: Optimal GC content (0.0 to 1.0)
            concentration: Concentration parameter (> 0.0), controls sharpness

        Returns:
            FitnessConfigBuilder for chaining

        Example:
            >>> # Single fitness component
            >>> fitness = cv.FitnessConfig.with_gc_content(0.5, 2.0).build()
            >>>
            >>> # Multiple components
            >>> fitness = cv.FitnessConfig.with_gc_content(0.5, 2.0) \\
            ...     .with_length(20000, 0.5) \\
            ...     .build()
        """
        ...

    @staticmethod
    def with_length(optimum: int, std_dev: float) -> FitnessConfigBuilder:
        """Start building with length-based fitness.

        Args:
            optimum: Optimal sequence length in bases (> 0)
            std_dev: Standard deviation (> 0.0)

        Returns:
            FitnessConfigBuilder for chaining

        Example:
            >>> fitness = cv.FitnessConfig.with_length(20000, 0.5).build()
        """
        ...

    @staticmethod
    def with_similarity(shape: float) -> FitnessConfigBuilder:
        """Start building with sequence similarity fitness.

        Args:
            shape: Shape parameter (> 0.0), controls decline rate

        Returns:
            FitnessConfigBuilder for chaining

        Example:
            >>> fitness = cv.FitnessConfig.with_similarity(2.0).build()
        """
        ...

    def __repr__(self) -> str: ...

class FitnessConfigBuilder:
    """Builder for constructing fitness configurations with multiple components.

    Use static methods on FitnessConfig to create a builder, then chain
    additional fitness components with .with_*() methods.

    Example:
        >>> fitness = cv.FitnessConfig.with_gc_content(0.5, 2.0) \\
        ...     .with_length(20000, 0.5) \\
        ...     .with_similarity(2.0) \\
        ...     .build()
    """

    def with_gc_content(self, optimum: float, concentration: float) -> FitnessConfigBuilder:
        """Add GC content fitness to the configuration.

        Args:
            optimum: Optimal GC content (0.0 to 1.0)
            concentration: Concentration parameter (> 0.0)

        Returns:
            Self for chaining

        Raises:
            ValueError: If GC content fitness is already set
        """
        ...

    def with_length(self, optimum: int, std_dev: float) -> FitnessConfigBuilder:
        """Add length-based fitness to the configuration.

        Args:
            optimum: Optimal sequence length in bases (> 0)
            std_dev: Standard deviation (> 0.0)

        Returns:
            Self for chaining

        Raises:
            ValueError: If length fitness is already set
        """
        ...

    def with_similarity(self, shape: float) -> FitnessConfigBuilder:
        """Add sequence similarity fitness to the configuration.

        Args:
            shape: Shape parameter (> 0.0)

        Returns:
            Self for chaining

        Raises:
            ValueError: If similarity fitness is already set
        """
        ...

    def build(self) -> FitnessConfig:
        """Build the final fitness configuration.

        Returns:
            FitnessConfig with all specified components
        """
        ...

    def __repr__(self) -> str: ...

class Simulation:
    """Simulation engine for evolutionary processes."""
    def __init__(
        self,
        structure: RepeatStructure,
        mutation: MutationConfig,
        recombination: RecombinationConfig,
        fitness: FitnessConfig,
        config: SimulationConfig,
    ) -> None:
        """Create a new simulation with uniform initial sequences."""
        ...

    @staticmethod
    def from_sequences(
        source_type: str,
        source_path: str,
        structure: RepeatStructure,
        mutation: MutationConfig,
        recombination: RecombinationConfig,
        fitness: FitnessConfig,
        config: SimulationConfig,
        sim_id: Optional[str] = None,
        generation: Optional[int] = None,
    ) -> Simulation:
        """Create a simulation from custom sequences.

        Args:
            source_type: Type of input - "fasta", "json", or "database"
            source_path: Path to FASTA/JSON file, JSON string, or database path
            structure: Repeat structure configuration
            mutation: Mutation configuration
            recombination: Recombination configuration
            fitness: Fitness configuration
            config: Simulation configuration
            sim_id: Simulation ID (required for "database" source)
            generation: Generation to load (optional for "database" source, defaults to last)

        Returns:
            Simulation instance initialized with custom sequences

        Examples:
            >>> # From FASTA file
            >>> sim = Simulation.from_sequences(
            ...     "fasta", "sequences.fasta", structure, mutation,
            ...     recombination, fitness, config
            ... )

            >>> # From JSON string
            >>> json_data = '[{"id": "seq1", "seq": "ACGT..."}]'
            >>> sim = Simulation.from_sequences(
            ...     "json", json_data, structure, mutation,
            ...     recombination, fitness, config
            ... )

            >>> # From database
            >>> sim = Simulation.from_sequences(
            ...     "database", "simulation.db", structure, mutation,
            ...     recombination, fitness, config,
            ...     sim_id="exp1", generation=1000
            ... )
        """
        ...

    @staticmethod
    def from_checkpoint(db_path: str, sim_id: str) -> Simulation:
        """Resume a simulation from a checkpoint.

        Args:
            db_path: Path to the database file
            sim_id: Simulation ID to resume

        Returns:
            Simulation instance ready to continue from checkpoint

        Example:
            >>> sim = Simulation.from_checkpoint("simulation.db", "exp1")
            >>> sim.run_for(100)  # Continue for 100 more generations
        """
        ...

    def population(self) -> Population:
        """Get the current population."""
        ...

    def generation(self) -> int:
        """Get the current generation number."""
        ...

    def step(self) -> None:
        """Advance simulation by one generation."""
        ...

    def run_for(self, generations: int) -> None:
        """Run simulation for a specific number of generations."""
        ...

    def run(self) -> None:
        """Run simulation for the configured number of generations."""
        ...

    def __repr__(self) -> str: ...

class SimulationBuilder:
    """Builder for constructing simulations with a fluent API.

    Provides an ergonomic way to configure and create simulations with
    sensible defaults and comprehensive validation.

    Examples:
        >>> # Simple simulation with defaults
        >>> sim = cv.SimulationBuilder() \\
        ...     .population_size(100) \\
        ...     .generations(50) \\
        ...     .repeat_structure(171, 12, 10) \\
        ...     .build()

        >>> # With mutation and recombination
        >>> sim = cv.SimulationBuilder() \\
        ...     .population_size(100) \\
        ...     .generations(50) \\
        ...     .repeat_structure(171, 12, 10) \\
        ...     .mutation_rate(0.0001) \\
        ...     .recombination(0.01, 0.7, 0.1) \\
        ...     .seed(42) \\
        ...     .build()

        >>> # With fitness/selection
        >>> fitness = cv.FitnessConfig.with_gc_content(0.5, 2.0).build()
        >>> sim = cv.SimulationBuilder() \\
        ...     .population_size(100) \\
        ...     .generations(50) \\
        ...     .repeat_structure(171, 12, 10) \\
        ...     .fitness(fitness) \\
        ...     .build()

        >>> # From FASTA file
        >>> sim = cv.SimulationBuilder() \\
        ...     .population_size(100) \\
        ...     .generations(50) \\
        ...     .init_from_fasta("sequences.fasta") \\
        ...     .mutation_rate(0.0001) \\
        ...     .build()
    """

    def __init__(self) -> None:
        """Create a new simulation builder with default values."""
        ...

    def population_size(self, size: int) -> SimulationBuilder:
        """Set the population size (required).

        Args:
            size: Number of diploid individuals in the population

        Returns:
            Self for method chaining
        """
        ...

    def generations(self, generations: int) -> SimulationBuilder:
        """Set the number of generations to run (required).

        Args:
            generations: Total number of generations to simulate

        Returns:
            Self for method chaining
        """
        ...

    def repeat_structure(
        self,
        ru_length: int,
        rus_per_hor: int,
        hors_per_chr: int,
    ) -> SimulationBuilder:
        """Set the repeat structure parameters (required for uniform/random init).

        Args:
            ru_length: Repeat unit length in bases
            rus_per_hor: Repeat units per higher-order repeat
            hors_per_chr: Higher-order repeats per chromosome

        Returns:
            Self for method chaining
        """
        ...

    def chromosomes_per_haplotype(self, chrs_per_hap: int) -> SimulationBuilder:
        """Set the number of chromosomes per haplotype (default: 1).

        Args:
            chrs_per_hap: Number of chromosomes in each haplotype

        Returns:
            Self for method chaining
        """
        ...

    def init_uniform(self, base: Nucleotide) -> SimulationBuilder:
        """Initialize with uniform sequences using the specified base.

        Args:
            base: Nucleotide to fill all positions (e.g., Nucleotide.A())

        Returns:
            Self for method chaining
        """
        ...

    def init_random(self) -> SimulationBuilder:
        """Initialize with random sequences.

        Each position gets a random base from the alphabet.

        Returns:
            Self for method chaining
        """
        ...

    def init_from_fasta(self, path: str) -> SimulationBuilder:
        """Initialize from a FASTA file.

        Args:
            path: Path to FASTA file containing sequences

        Returns:
            Self for method chaining
        """
        ...

    def init_from_json(self, input: str) -> SimulationBuilder:
        """Initialize from JSON (file path or JSON string).

        Args:
            input: Path to JSON file or JSON string

        Returns:
            Self for method chaining
        """
        ...

    def init_from_checkpoint(
        self,
        db_path: str,
        sim_id: str,
        generation: Optional[int] = None,
    ) -> SimulationBuilder:
        """Initialize from a checkpoint database.

        Args:
            db_path: Path to the checkpoint database
            sim_id: Simulation ID to load
            generation: Optional generation to load (defaults to last)

        Returns:
            Self for method chaining
        """
        ...

    def alphabet(self: Alphabet) -> SimulationBuilder:
        """Set the alphabet (default: DNA).

        Args:
            alphabet: Alphabet to use for sequences

        Returns:
            Self for method chaining
        """
        ...

    def mutation_rate(self, rate: float) -> SimulationBuilder:
        """Set the mutation rate (default: 0.0).

        Args:
            rate: Per-base mutation rate (0.0 to 1.0)

        Returns:
            Self for method chaining
        """
        ...

    def recombination(
        self,
        break_prob: float,
        crossover_prob: float,
        gc_extension_prob: float,
    ) -> SimulationBuilder:
        """Set recombination parameters.

        Args:
            break_prob: Probability of DNA strand break
            crossover_prob: Probability of crossover
            gc_extension_prob: Probability of gene conversion extension

        Returns:
            Self for method chaining
        """
        ...

    def fitness(self, fitness: FitnessConfig) -> SimulationBuilder:
        """Set the fitness configuration (default: neutral).

        If not specified, neutral fitness (no selection) is used.

        Args:
            fitness: FitnessConfig to use for selection

        Returns:
            Self for method chaining

        Example:
            >>> fitness = cv.FitnessConfig.with_gc_content(0.5, 2.0).build()
            >>> sim = cv.SimulationBuilder() \\
            ...     .population_size(100) \\
            ...     .generations(50) \\
            ...     .repeat_structure(171, 12, 10) \\
            ...     .fitness(fitness) \\
            ...     .build()
        """
        ...

    def seed(self, seed: int) -> SimulationBuilder:
        """Set the random seed for reproducibility (default: random).

        Args:
            seed: Random seed value

        Returns:
            Self for method chaining
        """
        ...

    def build(self) -> Simulation:
        """Build and validate the simulation.

        Returns:
            Simulation ready to run

        Raises:
            ValueError: If required parameters are missing or invalid
            RuntimeError: If simulation creation fails
        """
        ...

    def __repr__(self) -> str: ...

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
