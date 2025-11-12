"""Type stubs for centrevo.plotting module."""

from typing import Dict, List, Optional, Tuple, Any
from matplotlib.figure import Figure
import pyarrow as pa

def plot_diversity_trajectory(
    diversity_data: List[Dict[str, Any]],
    metric: str = "nucleotide_diversity",
    figsize: Tuple[int, int] = (10, 6),
    title: Optional[str] = None,
    xlabel: str = "Generation",
    ylabel: Optional[str] = None,
    save_path: Optional[str] = None,
) -> Figure: ...

def plot_ld_decay(
    ld_data: List[Dict[str, float]],
    figsize: Tuple[int, int] = (10, 6),
    title: str = "Linkage Disequilibrium Decay",
    xlabel: str = "Distance (bp)",
    ylabel: str = "r²",
    save_path: Optional[str] = None,
) -> Figure: ...

def plot_distance_matrix(
    distance_matrix: List[List[float]],
    figsize: Tuple[int, int] = (10, 8),
    title: str = "Genetic Distance Matrix",
    cmap: str = "viridis",
    save_path: Optional[str] = None,
) -> Figure: ...

def plot_nucleotide_composition(
    composition: Dict[str, float],
    figsize: Tuple[int, int] = (8, 6),
    title: str = "Nucleotide Composition",
    save_path: Optional[str] = None,
) -> Figure: ...

def plot_multiple_diversity_metrics(
    diversity_data: List[Dict[str, Any]],
    metrics: Optional[List[str]] = None,
    figsize: Tuple[int, int] = (12, 8),
    save_path: Optional[str] = None,
) -> Figure: ...

def export_to_pyarrow_table(
    data: List[Dict[str, Any]],
    schema: Optional[pa.Schema] = None,
) -> pa.Table: ...

def export_distance_matrix_to_pyarrow(
    distance_data: List[Dict[str, Any]]
) -> pa.Table: ...

def export_ld_decay_to_pyarrow(
    ld_data: List[Dict[str, float]]
) -> pa.Table: ...
