"""
Plotting utilities for Centrevo analysis results.

Provides simple, common visualization functions using matplotlib.
For complex visualizations, users can work directly with the PyArrow data.
"""

from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np
from typing import List, Dict, Optional, Tuple, Any
import pyarrow as pa


def plot_diversity_trajectory(
    diversity_data: List[Dict[str, Any]],
    metric: str = "nucleotide_diversity",
    figsize: Tuple[int, int] = (10, 6),
    title: Optional[str] = None,
    xlabel: str = "Generation",
    ylabel: Optional[str] = None,
    save_path: Optional[str] = None,
) -> Figure:
    """
    Plot diversity metric over generations.
    
    Args:
        diversity_data: List of dicts from export_diversity_metrics(),
                       each containing generation and diversity metrics
        metric: Which metric to plot ('nucleotide_diversity', 'tajimas_d', 
                'wattersons_theta', 'haplotype_diversity')
        figsize: Figure size as (width, height)
        title: Plot title (auto-generated if None)
        xlabel: X-axis label
        ylabel: Y-axis label (auto-generated if None)
        save_path: Path to save figure (displays if None)
    
    Returns:
        matplotlib Figure object
    
    Example:
        ```python
        import centrevo
        import pyarrow as pa
        
        # Collect data over generations
        diversity_data = []
        for gen in range(0, 100, 10):
            pop = load_population(gen)  # Your loading logic
            metrics = centrevo.export_diversity_metrics(pop, 0)
            diversity_data.append(metrics)
        
        # Plot
        fig = plot_diversity_trajectory(diversity_data, metric='nucleotide_diversity')
        ```
    """
    generations = [d["generation"] for d in diversity_data]
    values = [d[metric] for d in diversity_data]
    
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(generations, values, marker='o', linewidth=2, markersize=6)
    
    ax.set_xlabel(xlabel, fontsize=12)
    
    if ylabel is None:
        ylabel_map = {
            "nucleotide_diversity": "Nucleotide Diversity (π)",
            "tajimas_d": "Tajima's D",
            "wattersons_theta": "Watterson's θ",
            "haplotype_diversity": "Haplotype Diversity",
        }
        ylabel = ylabel_map.get(metric, metric)
    ax.set_ylabel(ylabel, fontsize=12)
    
    if title is None:
        title = f"{ylabel} Over Time"
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig


def plot_ld_decay(
    ld_data: List[Dict[str, float]],
    figsize: Tuple[int, int] = (10, 6),
    title: str = "Linkage Disequilibrium Decay",
    xlabel: str = "Distance (bp)",
    ylabel: str = "r²",
    save_path: Optional[str] = None,
) -> Figure:
    """
    Plot LD decay (r² vs distance).
    
    Args:
        ld_data: List of dicts from export_ld_decay() with 'distance' and 'r_squared' keys
        figsize: Figure size as (width, height)
        title: Plot title
        xlabel: X-axis label
        ylabel: Y-axis label
        save_path: Path to save figure (displays if None)
    
    Returns:
        matplotlib Figure object
    
    Example:
        ```python
        import centrevo
        
        pop = load_population()
        ld_data = centrevo.export_ld_decay(pop, max_distance=1000, step_size=10, 
                                           chromosome_idx=0, haplotype_idx=0)
        
        fig = plot_ld_decay(ld_data)
        ```
    """
    distances = [d["distance"] for d in ld_data]
    r_squared = [d["r_squared"] for d in ld_data]
    
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(distances, r_squared, marker='o', linewidth=2, markersize=4, alpha=0.7)
    
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    ax.set_ylim(-0.05, 1.05)
    ax.axhline(y=0.5, color='red', linestyle='--', alpha=0.5, label='r² = 0.5')
    ax.axhline(y=0.8, color='orange', linestyle='--', alpha=0.5, label='r² = 0.8')
    
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig


def plot_distance_matrix(
    distance_matrix: List[List[float]],
    figsize: Tuple[int, int] = (10, 8),
    title: str = "Genetic Distance Matrix",
    cmap: str = "viridis",
    save_path: Optional[str] = None,
) -> Figure:
    """
    Plot distance matrix as a heatmap.
    
    Args:
        distance_matrix: 2D list from distance_matrix() function
        figsize: Figure size as (width, height)
        title: Plot title
        cmap: Matplotlib colormap name
        save_path: Path to save figure (displays if None)
    
    Returns:
        matplotlib Figure object
    
    Example:
        ```python
        import centrevo
        
        pop = load_population()
        matrix = centrevo.distance_matrix(pop, chromosome_idx=0)
        
        fig = plot_distance_matrix(matrix)
        ```
    """
    matrix_array = np.array(distance_matrix)
    
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(matrix_array, cmap=cmap, aspect='auto')
    
    ax.set_xlabel("Sequence Index", fontsize=12)
    ax.set_ylabel("Sequence Index", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Genetic Distance", rotation=270, labelpad=20, fontsize=12)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig


def plot_nucleotide_composition(
    composition: Dict[str, float],
    figsize: Tuple[int, int] = (8, 6),
    title: str = "Nucleotide Composition",
    save_path: Optional[str] = None,
) -> Figure:
    """
    Plot nucleotide composition as a bar chart.
    
    Args:
        composition: Dict from nucleotide_composition() function
        figsize: Figure size as (width, height)
        title: Plot title
        save_path: Path to save figure (displays if None)
    
    Returns:
        matplotlib Figure object
    
    Example:
        ```python
        import centrevo
        
        pop = load_population()
        comp = centrevo.nucleotide_composition(pop, None, None, None)
        
        fig = plot_nucleotide_composition(comp)
        ```
    """
    # Extract nucleotide letters (A, C, G, T)
    nucleotides = sorted(composition.keys())
    frequencies = [composition[nuc] for nuc in nucleotides]
    
    # Color mapping
    colors = {'A': '#FF6B6B', 'C': '#4ECDC4', 'G': '#FFE66D', 'T': '#95E1D3'}
    bar_colors = [colors.get(nuc, 'gray') for nuc in nucleotides]
    
    fig, ax = plt.subplots(figsize=figsize)
    bars = ax.bar(nucleotides, frequencies, color=bar_colors, alpha=0.8, edgecolor='black')
    
    ax.set_xlabel("Nucleotide", fontsize=12)
    ax.set_ylabel("Frequency", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_ylim(0, max(frequencies) * 1.1)
    
    # Add value labels on bars
    for bar, freq in zip(bars, frequencies):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{freq:.3f}',
                ha='center', va='bottom', fontsize=10)
    
    ax.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig


def plot_multiple_diversity_metrics(
    diversity_data: List[Dict[str, Any]],
    metrics: Optional[List[str]] = None,
    figsize: Tuple[int, int] = (12, 8),
    save_path: Optional[str] = None,
) -> Figure:
    """
    Plot multiple diversity metrics in subplots.
    
    Args:
        diversity_data: List of dicts from export_diversity_metrics()
        metrics: List of metrics to plot (all if None)
        figsize: Figure size as (width, height)
        save_path: Path to save figure (displays if None)
    
    Returns:
        matplotlib Figure object
    
    Example:
        ```python
        diversity_data = [...]  # List of diversity metrics over time
        
        fig = plot_multiple_diversity_metrics(
            diversity_data,
            metrics=['nucleotide_diversity', 'tajimas_d', 'haplotype_diversity']
        )
        ```
    """
    if metrics is None:
        metrics = ['nucleotide_diversity', 'tajimas_d', 
                   'wattersons_theta', 'haplotype_diversity']
    
    generations = [d["generation"] for d in diversity_data]
    
    n_metrics = len(metrics)
    fig, axes = plt.subplots(n_metrics, 1, figsize=figsize, sharex=True)
    
    if n_metrics == 1:
        axes = [axes]
    
    metric_labels = {
        "nucleotide_diversity": "π",
        "tajimas_d": "Tajima's D",
        "wattersons_theta": "θ_W",
        "haplotype_diversity": "Haplotype Diversity",
    }
    
    for ax, metric in zip(axes, metrics):
        values = [d[metric] for d in diversity_data]
        ax.plot(generations, values, marker='o', linewidth=2, markersize=6)
        ax.set_ylabel(metric_labels.get(metric, metric), fontsize=11)
        ax.grid(True, alpha=0.3)
    
    axes[-1].set_xlabel("Generation", fontsize=12)
    fig.suptitle("Diversity Metrics Over Time", fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig


def export_to_pyarrow_table(
    data: List[Dict[str, Any]],
    schema: Optional[pa.Schema] = None,
) -> pa.Table:
    """
    Convert list of dictionaries to PyArrow Table.
    
    This helper function makes it easy to convert Centrevo analysis results
    to PyArrow format, which can then be converted to pandas or polars.
    
    Args:
        data: List of dictionaries (e.g., from export_diversity_metrics())
        schema: Optional PyArrow schema (inferred if None)
    
    Returns:
        PyArrow Table
    
    Example:
        ```python
        import centrevo
        import pyarrow as pa
        import pandas as pd
        import polars as pl
        
        # Collect diversity metrics
        diversity_data = [...]
        
        # Convert to PyArrow Table
        table = export_to_pyarrow_table(diversity_data)
        
        # Convert to pandas or polars
        df_pandas = table.to_pandas()
        df_polars = pl.from_arrow(table)
        ```
    """
    if schema is None:
        return pa.Table.from_pylist(data)
    else:
        return pa.Table.from_pylist(data, schema=schema)


def export_distance_matrix_to_pyarrow(
    distance_data: List[Dict[str, Any]]
) -> pa.Table:
    """
    Convert distance matrix export to PyArrow Table.
    
    Args:
        distance_data: List of dicts from export_distance_matrix()
    
    Returns:
        PyArrow Table with columns: sequence_i, sequence_j, distance
    
    Example:
        ```python
        import centrevo
        
        pop = load_population()
        dist_data = centrevo.export_distance_matrix(pop, 0)
        table = export_distance_matrix_to_pyarrow(dist_data)
        
        # Now use with pandas/polars
        import pandas as pd
        df = table.to_pandas()
        ```
    """
    return pa.Table.from_pylist(distance_data)


def export_ld_decay_to_pyarrow(
    ld_data: List[Dict[str, float]]
) -> pa.Table:
    """
    Convert LD decay export to PyArrow Table.
    
    Args:
        ld_data: List of dicts from export_ld_decay()
    
    Returns:
        PyArrow Table with columns: distance, r_squared
    
    Example:
        ```python
        import centrevo
        
        pop = load_population()
        ld_data = centrevo.export_ld_decay(pop, 1000, 10, 0, 0)
        table = export_ld_decay_to_pyarrow(ld_data)
        
        # Now use with pandas/polars
        import polars as pl
        df = pl.from_arrow(table)
        ```
    """
    return pa.Table.from_pylist(ld_data)
