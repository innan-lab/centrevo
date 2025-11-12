#!/usr/bin/env python3
"""
Example: Python Analysis with Centrevo

Demonstrates how to:
1. Use analysis functions from Python
2. Export data to PyArrow format
3. Convert to pandas/polars
4. Create visualizations

Requirements:
    pip install centrevo pyarrow matplotlib numpy pandas polars
"""

import pandas as pd
import polars as pl

import centrevo
from centrevo.plotting import (export_ld_decay_to_pyarrow,
                               export_to_pyarrow_table, plot_distance_matrix,
                               plot_nucleotide_composition)


def example_diversity_analysis():
    """Demonstrate diversity analysis."""
    print("=" * 60)
    print("Diversity Analysis Example")
    print("=" * 60)
    
    # Create a simple population
    alphabet = centrevo.Alphabet.dna()
    base_a = centrevo.Nucleotide.A()
    
    structure = centrevo.RepeatStructure(
        alphabet=alphabet,
        init_base=base_a,
        ru_length=171,
        rus_per_hor=12,
        hors_per_chr=50,
        chrs_per_hap=1,
    )
    
    pop = centrevo.create_initial_population(size=20, structure=structure)
    
    # Calculate diversity metrics
    pi = centrevo.nucleotide_diversity(pop, chromosome_idx=0)
    tajima_d = centrevo.tajimas_d(pop, chromosome_idx=0)
    theta_w = centrevo.wattersons_theta(pop, chromosome_idx=0)
    hap_div = centrevo.haplotype_diversity(pop, chromosome_idx=0)
    
    print("\nDiversity Metrics:")
    print(f"  Nucleotide diversity (π): {pi:.6f}")
    print(f"  Tajima's D: {tajima_d:.6f}")
    print(f"  Watterson's θ: {theta_w:.6f}")
    print(f"  Haplotype diversity: {hap_div:.6f}")


def example_export_to_pyarrow():
    """Demonstrate PyArrow export and conversion."""
    print("\n" + "=" * 60)
    print("PyArrow Export Example")
    print("=" * 60)
    
    # Create population
    alphabet = centrevo.Alphabet.dna()
    base_a = centrevo.Nucleotide.A()
    structure = centrevo.RepeatStructure(
        alphabet=alphabet,
        init_base=base_a,
        ru_length=171,
        rus_per_hor=12,
        hors_per_chr=50,
        chrs_per_hap=1,
    )
    
    pop = centrevo.create_initial_population(size=10, structure=structure)
    
    # Export diversity metrics
    metrics = centrevo.export_diversity_metrics(pop, chromosome_idx=0)
    print(f"\nExported metrics: {metrics}")
    
    # Convert to PyArrow Table
    # In a real scenario, you'd collect metrics over multiple generations
    diversity_data = [metrics]  # List of dicts
    table = export_to_pyarrow_table(diversity_data)
    
    print("\nPyArrow Table Schema:")
    print(table.schema)
    
    # Convert to pandas
    df_pandas = table.to_pandas()
    print("\nPandas DataFrame:")
    print(df_pandas)
    
    # Convert to polars
    df_polars = pl.from_arrow(table)
    print("\nPolars DataFrame:")
    print(df_polars)


def example_ld_analysis():
    """Demonstrate LD analysis."""
    print("\n" + "=" * 60)
    print("LD Analysis Example")
    print("=" * 60)
    
    # Create population
    alphabet = centrevo.Alphabet.dna()
    base_a = centrevo.Nucleotide.A()
    structure = centrevo.RepeatStructure(
        alphabet=alphabet,
        init_base=base_a,
        ru_length=171,
        rus_per_hor=12,
        hors_per_chr=50,
        chrs_per_hap=1,
    )
    
    pop = centrevo.create_initial_population(size=20, structure=structure)
    
    # Calculate LD between two positions
    ld_stats = centrevo.linkage_disequilibrium(
        pop,
        pos1=100,
        pos2=500,
        chromosome_idx=0,
        haplotype_idx=0,
    )
    
    print("\nLD Statistics (pos 100 vs 500):")
    if ld_stats is not None:
        print(f"  D: {ld_stats['D']:.6f}")
        print(f"  D': {ld_stats['D_prime']:.6f}")
        print(f"  r²: {ld_stats['r_squared']:.6f}")
    else:
        print("  LD calculation returned None (not enough variation)")
    
    # Export LD decay data
    ld_decay_data = centrevo.export_ld_decay(
        pop,
        chromosome_idx=0,
        haplotype_idx=0,
        max_distance=1000,
        bin_size=50,
    )
    
    print(f"\nLD Decay data points: {len(ld_decay_data)}")
    
    # Convert to PyArrow and then polars
    table = export_ld_decay_to_pyarrow(ld_decay_data)
    df = pl.from_arrow(table)
    print("\nLD Decay DataFrame (first 5 rows):")
    print(df.head())


def example_composition_analysis():
    """Demonstrate composition analysis."""
    print("\n" + "=" * 60)
    print("Composition Analysis Example")
    print("=" * 60)
    
    # Create population
    alphabet = centrevo.Alphabet.dna()
    base_a = centrevo.Nucleotide.A()
    structure = centrevo.RepeatStructure(
        alphabet=alphabet,
        init_base=base_a,
        ru_length=171,
        rus_per_hor=12,
        hors_per_chr=50,
        chrs_per_hap=1,
    )
    
    pop = centrevo.create_initial_population(size=10, structure=structure)
    
    # Population-level GC content
    gc_pop = centrevo.gc_content(pop, None, None, None)
    print(f"\nPopulation GC content: {gc_pop:.4f}")
    
    # Individual-level GC content
    gc_ind = centrevo.gc_content(pop, individual_idx=0, haplotype_idx=None, chromosome_idx=None)
    print(f"Individual 0 GC content: {gc_ind:.4f}")
    
    # Nucleotide composition
    comp = centrevo.nucleotide_composition(pop, None, None, None)
    print("\nNucleotide composition:")
    for nuc, freq in sorted(comp.items()):
        print(f"  {nuc}: {freq:.4f}")


def example_distance_analysis():
    """Demonstrate distance analysis."""
    print("\n" + "=" * 60)
    print("Distance Analysis Example")
    print("=" * 60)
    
    # Create population
    alphabet = centrevo.Alphabet.dna()
    base_a = centrevo.Nucleotide.A()
    structure = centrevo.RepeatStructure(
        alphabet=alphabet,
        init_base=base_a,
        ru_length=171,
        rus_per_hor=12,
        hors_per_chr=50,
        chrs_per_hap=1,
    )
    
    pop = centrevo.create_initial_population(size=5, structure=structure)
    
    # Calculate distance matrix
    matrix = centrevo.distance_matrix(pop, chromosome_idx=0)
    
    print(f"\nDistance matrix shape: {len(matrix)}x{len(matrix[0])}")
    print(f"Number of sequences: {len(matrix)} (5 individuals × 2 haplotypes)")
    
    # Export to structured format
    dist_data = centrevo.export_distance_matrix(pop, chromosome_idx=0)
    print(f"Distance data points: {len(dist_data)}")
    
    # Convert to pandas
    df = pd.DataFrame(dist_data)
    print("\nDistance DataFrame (first 10 rows):")
    print(df.head(10))


def example_visualization():
    """Demonstrate simple visualizations (requires matplotlib)."""
    print("\n" + "=" * 60)
    print("Visualization Example")
    print("=" * 60)
    
    try:
        import matplotlib
        matplotlib.use('Agg')  # Non-interactive backend
    except ImportError:
        print("matplotlib not installed, skipping visualization")
        return
    
    # Create population
    alphabet = centrevo.Alphabet.dna()
    base_a = centrevo.Nucleotide.A()
    structure = centrevo.RepeatStructure(
        alphabet=alphabet,
        init_base=base_a,
        ru_length=171,
        rus_per_hor=12,
        hors_per_chr=50,
        chrs_per_hap=1,
    )
    
    pop = centrevo.create_initial_population(size=10, structure=structure)
    
    # Plot nucleotide composition
    comp = centrevo.nucleotide_composition(pop, None, None, None)
    _ = plot_nucleotide_composition(comp)
    print("\nNucleotide composition plot created")
    
    # Plot distance matrix
    matrix = centrevo.distance_matrix(pop, chromosome_idx=0)
    _ = plot_distance_matrix(matrix, figsize=(8, 6))
    print("Distance matrix heatmap created")
    
    # Note: For trajectory plots, you'd need data from multiple generations
    print("\nTo create trajectory plots, collect diversity metrics over multiple generations")
    print("and use plot_diversity_trajectory() or plot_multiple_diversity_metrics()")


def main():
    """Run all examples."""
    print("\n" + "=" * 60)
    print("Centrevo Python Analysis Examples")
    print("=" * 60)
    
    example_diversity_analysis()
    example_export_to_pyarrow()
    example_ld_analysis()
    example_composition_analysis()
    example_distance_analysis()
    example_visualization()
    
    print("\n" + "=" * 60)
    print("All examples completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()
