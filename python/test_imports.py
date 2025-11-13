#!/usr/bin/env python3
"""
Test script to verify Centrevo Python package imports and basic functionality.

Run this after building with: maturin develop --release
"""

import sys

def test_imports():
    """Test that all expected symbols can be imported."""
    print("Testing imports from centrevo package...")

    try:
        import centrevo
        print(f"✓ Successfully imported centrevo (version {centrevo.__version__})")
    except ImportError as e:
        print(f"✗ Failed to import centrevo: {e}")
        return False

    # Test core classes
    core_classes = [
        'Nucleotide', 'Alphabet', 'Chromosome', 'Haplotype',
        'Individual', 'Population', 'RepeatStructure', 'SimulationConfig',
        'Recorder', 'QueryBuilder', 'RecordingStrategy',
    ]

    for cls_name in core_classes:
        if hasattr(centrevo, cls_name):
            print(f"✓ {cls_name} is available")
        else:
            print(f"✗ {cls_name} is NOT available")
            return False

    # Test core functions
    if hasattr(centrevo, 'create_initial_population'):
        print("✓ create_initial_population is available")
    else:
        print("✗ create_initial_population is NOT available")
        return False

    # Test analysis functions
    analysis_functions = [
        'nucleotide_diversity', 'tajimas_d', 'wattersons_theta',
        'haplotype_diversity', 'linkage_disequilibrium', 'ld_decay',
        'haplotype_blocks', 'pairwise_distances', 'distance_matrix',
        'gc_content', 'nucleotide_composition', 'count_segregating_sites',
    ]

    for func_name in analysis_functions:
        if hasattr(centrevo, func_name):
            print(f"✓ {func_name} is available")
        else:
            print(f"✗ {func_name} is NOT available")
            return False

    # Test export functions
    export_functions = [
        'export_diversity_metrics', 'export_distance_matrix', 'export_ld_decay',
    ]

    for func_name in export_functions:
        if hasattr(centrevo, func_name):
            print(f"✓ {func_name} is available")
        else:
            print(f"✗ {func_name} is NOT available")
            return False

    # Test plotting module
    try:
        from centrevo import plotting
        print("✓ plotting module is available")
    except ImportError as e:
        print(f"✗ Failed to import plotting module: {e}")
        return False

    plotting_functions = [
        'plot_diversity_trajectory', 'plot_ld_decay',
        'plot_distance_matrix', 'plot_nucleotide_composition',
        'plot_multiple_diversity_metrics',
        'export_to_pyarrow_table', 'export_distance_matrix_to_pyarrow',
        'export_ld_decay_to_pyarrow',
    ]

    for func_name in plotting_functions:
        if hasattr(plotting, func_name):
            print(f"✓ plotting.{func_name} is available")
        else:
            print(f"✗ plotting.{func_name} is NOT available")
            return False

    return True


def test_basic_functionality():
    """Test basic functionality."""
    print("\nTesting basic functionality...")

    import centrevo

    # Create basic objects
    alphabet = centrevo.Alphabet.dna()
    print("✓ Created DNA alphabet")

    base_a = centrevo.Nucleotide.A()
    print("✓ Created Nucleotide A")

    structure = centrevo.RepeatStructure(
        alphabet=alphabet,
        init_base=base_a,
        ru_length=171,
        rus_per_hor=12,
        hors_per_chr=10,
        chrs_per_hap=1,
    )
    print(f"✓ Created RepeatStructure: {structure}")

    # Create small population
    pop = centrevo.create_initial_population(size=10, structure=structure)
    print(f"✓ Created initial population: {pop}")

    # Calculate diversity metrics
    pi = centrevo.nucleotide_diversity(pop, chromosome_idx=0)
    print(f"✓ Calculated nucleotide diversity: π = {pi:.6f}")

    tajima = centrevo.tajimas_d(pop, chromosome_idx=0)
    print(f"✓ Calculated Tajima's D: D = {tajima:.6f}")

    # Export metrics
    metrics = centrevo.export_diversity_metrics(pop, chromosome_idx=0)
    print(f"✓ Exported diversity metrics: {list(metrics.keys())}")

    # Test composition
    comp = centrevo.nucleotide_composition(pop, None, None, None)
    print(f"✓ Calculated nucleotide composition: {comp}")

    print("\n✓ All basic functionality tests passed!")
    return True


def main():
    """Run all tests."""
    print("=" * 70)
    print("Centrevo Python Package Import Test")
    print("=" * 70)

    if not test_imports():
        print("\n✗ Import tests FAILED")
        return 1

    print("\n" + "=" * 70)

    try:
        if not test_basic_functionality():
            print("\n✗ Functionality tests FAILED")
            return 1
    except Exception as e:
        print(f"\n✗ Functionality tests FAILED with exception: {e}")
        import traceback
        traceback.print_exc()
        return 1

    print("\n" + "=" * 70)
    print("✓ ALL TESTS PASSED")
    print("=" * 70)
    return 0


if __name__ == "__main__":
    sys.exit(main())
