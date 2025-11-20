//! Python bindings for analysis functions.
//!
//! This module provides Python interfaces to all Centrevo analysis functionality,
//! with PyArrow integration for efficient data export.

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use centrevo_analysis::analysis::{
    nucleotide_diversity, tajimas_d, wattersons_theta, haplotype_diversity,
    linkage_disequilibrium, ld_decay, haplotype_blocks,
    pairwise_distances, distance_matrix,
    gc_content, nucleotide_composition,
};
use centrevo_analysis::analysis::polymorphism::count_segregating_sites;
use super::bindings::PyPopulation;

/// Calculate nucleotide diversity (π) for a population.
///
/// Args:
///     population: Population to analyze
///     chromosome_idx: Index of chromosome to analyze (default: 0)
///
/// Returns:
///     float: Average pairwise nucleotide diversity per site across all 2n sequences
///
/// Example:
///     ```python
///     pi = nucleotide_diversity_py(population, 0)
///     print(f"Nucleotide diversity: {pi:.6f}")
///     ```
#[pyfunction]
#[pyo3(name = "nucleotide_diversity")]
fn nucleotide_diversity_py(
    population: &PyPopulation,
    chromosome_idx: usize,
) -> PyResult<f64> {
    Ok(nucleotide_diversity(&population.inner, chromosome_idx))
}

/// Calculate Tajima's D statistic.
///
/// Tests the hypothesis of neutral evolution by comparing two estimates of θ.
/// Positive D suggests balancing selection or population contraction.
/// Negative D suggests purifying selection or population expansion.
///
/// Args:
///     population: Population to analyze
///     chromosome_idx: Index of chromosome to analyze
///
/// Returns:
///     float: Tajima's D statistic
#[pyfunction]
#[pyo3(name = "tajimas_d")]
fn tajimas_d_py(
    population: &PyPopulation,
    chromosome_idx: usize,
) -> PyResult<f64> {
    Ok(tajimas_d(&population.inner, chromosome_idx))
}

/// Calculate Watterson's estimator (θ_W).
///
/// Estimates θ = 4Nμ from the number of segregating sites.
///
/// Args:
///     population: Population to analyze
///     chromosome_idx: Index of chromosome to analyze
///
/// Returns:
///     float: Watterson's theta per site
#[pyfunction]
#[pyo3(name = "wattersons_theta")]
fn wattersons_theta_py(
    population: &PyPopulation,
    chromosome_idx: usize,
) -> PyResult<f64> {
    Ok(wattersons_theta(&population.inner, chromosome_idx))
}

/// Calculate haplotype diversity.
///
/// Returns the probability that two randomly chosen haplotypes are different.
///
/// Args:
///     population: Population to analyze
///     chromosome_idx: Index of chromosome to analyze
///
/// Returns:
///     float: Haplotype diversity (0.0 to 1.0)
#[pyfunction]
#[pyo3(name = "haplotype_diversity")]
fn haplotype_diversity_py(
    population: &PyPopulation,
    chromosome_idx: usize,
) -> PyResult<f64> {
    Ok(haplotype_diversity(&population.inner, chromosome_idx))
}

/// Calculate linkage disequilibrium statistics between two sites.
///
/// Args:
///     population: Population to analyze
///     pos1: Position of first site
///     pos2: Position of second site
///     chromosome_idx: Index of chromosome to analyze
///     haplotype_idx: Index of haplotype to analyze
///
/// Returns:
///     dict: Dictionary with keys 'D', 'D_prime', 'r_squared', or None if calculation fails
#[pyfunction]
#[pyo3(name = "linkage_disequilibrium")]
fn linkage_disequilibrium_py(
    py: Python,
    population: &PyPopulation,
    pos1: usize,
    pos2: usize,
    chromosome_idx: usize,
    haplotype_idx: usize,
) -> PyResult<Option<Py<PyDict>>> {
    let ld_stats = linkage_disequilibrium(
        &population.inner,
        pos1,
        pos2,
        chromosome_idx,
        haplotype_idx,
    );
    
    match ld_stats {
        Some(stats) => {
            let dict = PyDict::new(py);
            dict.set_item("D", stats.d)?;
            dict.set_item("D_prime", stats.d_prime)?;
            dict.set_item("r_squared", stats.r_squared)?;
            Ok(Some(dict.into()))
        }
        None => Ok(None),
    }
}

/// Calculate LD decay across a range of distances.
///
/// Args:
///     population: Population to analyze
///     chromosome_idx: Index of chromosome to analyze
///     haplotype_idx: Index of haplotype to analyze
///     max_distance: Maximum distance to analyze
///     bin_size: Size of distance bins
///
/// Returns:
///     dict: Dictionary with 'distances' and 'r_squared_values' as lists
#[pyfunction]
#[pyo3(name = "ld_decay")]
fn ld_decay_py(
    py: Python,
    population: &PyPopulation,
    chromosome_idx: usize,
    haplotype_idx: usize,
    max_distance: usize,
    bin_size: usize,
) -> PyResult<Py<PyDict>> {
    let decay = ld_decay(
        &population.inner,
        chromosome_idx,
        haplotype_idx,
        max_distance,
        bin_size,
    );
    
    // Separate distances and r_squared values
    let (distances, r_squared_values): (Vec<usize>, Vec<f64>) = decay.into_iter().unzip();
    
    let dict = PyDict::new(py);
    dict.set_item("distances", PyList::new(py, &distances)?)?;
    dict.set_item("r_squared_values", PyList::new(py, &r_squared_values)?)?;
    
    Ok(dict.into())
}

/// Identify haplotype blocks based on LD threshold.
///
/// Args:
///     population: Population to analyze
///     chromosome_idx: Index of chromosome to analyze
///     haplotype_idx: Index of haplotype to analyze
///     r_squared_threshold: Minimum r² to consider sites in same block (default: 0.8)
///
/// Returns:
///     list: List of blocks, each block is a tuple (start, end)
#[pyfunction]
#[pyo3(name = "haplotype_blocks")]
fn haplotype_blocks_py(
    py: Python,
    population: &PyPopulation,
    chromosome_idx: usize,
    haplotype_idx: usize,
    r_squared_threshold: Option<f64>,
) -> PyResult<Py<PyList>> {
    let threshold = r_squared_threshold.unwrap_or(0.8);
    let blocks = haplotype_blocks(
        &population.inner,
        chromosome_idx,
        haplotype_idx,
        threshold,
    );
    
    // Convert Vec<(usize, usize)> to Python list of tuples
    let py_blocks: Vec<(usize, usize)> = blocks;
    
    Ok(PyList::new(py, py_blocks)?.into())
}

/// Calculate pairwise genetic distances between all sequences.
///
/// Args:
///     population: Population to analyze
///     chromosome_idx: Index of chromosome to analyze
///
/// Returns:
///     list: List of distances (normalized Hamming distances)
#[pyfunction]
#[pyo3(name = "pairwise_distances")]
fn pairwise_distances_py(
    py: Python,
    population: &PyPopulation,
    chromosome_idx: usize,
) -> PyResult<Py<PyList>> {
    let distances = pairwise_distances(&population.inner, chromosome_idx);
    Ok(PyList::new(py, distances)?.into())
}

/// Calculate full distance matrix between all sequences.
///
/// Args:
///     population: Population to analyze
///     chromosome_idx: Index of chromosome to analyze
///
/// Returns:
///     list: List of lists representing the distance matrix (2n × 2n)
#[pyfunction]
#[pyo3(name = "distance_matrix")]
fn distance_matrix_py(
    py: Python,
    population: &PyPopulation,
    chromosome_idx: usize,
) -> PyResult<Py<PyList>> {
    let matrix = distance_matrix(&population.inner, chromosome_idx);
    
    // Convert Vec<Vec<f64>> to Python list of lists
    let py_matrix: Vec<Py<PyList>> = matrix
        .iter()
        .map(|row| PyList::new(py, row).unwrap().into())
        .collect();
    
    Ok(PyList::new(py, py_matrix)?.into())
}

/// Calculate GC content.
///
/// Flexible function that calculates GC content at different levels based on
/// which parameters are provided:
/// - All None: population-level average
/// - individual_idx only: average across both haplotypes
/// - individual_idx + haplotype_idx: specific haplotype average
/// - All three: specific chromosome
///
/// Args:
///     population: Population to analyze
///     individual_idx: Optional individual index
///     haplotype_idx: Optional haplotype index (0 or 1)
///     chromosome_idx: Optional chromosome index
///
/// Returns:
///     float: GC content (0.0 to 1.0)
#[pyfunction]
#[pyo3(name = "gc_content")]
fn gc_content_py(
    population: &PyPopulation,
    individual_idx: Option<usize>,
    haplotype_idx: Option<usize>,
    chromosome_idx: Option<usize>,
) -> PyResult<f64> {
    Ok(gc_content(
        &population.inner,
        individual_idx,
        haplotype_idx,
        chromosome_idx,
    ))
}

/// Calculate nucleotide composition.
///
/// Flexible function that calculates nucleotide composition at different levels.
///
/// Args:
///     population: Population to analyze
///     individual_idx: Optional individual index
///     haplotype_idx: Optional haplotype index (0 or 1)
///     chromosome_idx: Optional chromosome index
///
/// Returns:
///     dict: Dictionary mapping nucleotide (str) to frequency (float)
#[pyfunction]
#[pyo3(name = "nucleotide_composition")]
fn nucleotide_composition_py(
    py: Python,
    population: &PyPopulation,
    individual_idx: Option<usize>,
    haplotype_idx: Option<usize>,
    chromosome_idx: Option<usize>,
) -> PyResult<Py<PyDict>> {
    let composition = nucleotide_composition(
        &population.inner,
        individual_idx,
        haplotype_idx,
        chromosome_idx,
    );
    
    let dict = PyDict::new(py);
    for (nuc, freq) in composition {
        dict.set_item(format!("{:?}", nuc), freq)?;
    }
    
    Ok(dict.into())
}

/// Count number of segregating sites in a population.
///
/// Args:
///     population: Population to analyze
///     chromosome_idx: Index of chromosome to analyze
///     haplotype_idx: Index of haplotype to analyze
///
/// Returns:
///     int: Number of polymorphic sites
#[pyfunction]
#[pyo3(name = "count_segregating_sites")]
fn count_segregating_sites_py(
    population: &PyPopulation,
    chromosome_idx: usize,
    haplotype_idx: usize,
) -> PyResult<usize> {
    Ok(count_segregating_sites(
        &population.inner,
        chromosome_idx,
        haplotype_idx,
    ))
}

/// Export diversity metrics to PyArrow table format.
///
/// Calculates multiple diversity metrics and returns them as a structured
/// dictionary that can be easily converted to PyArrow Table.
///
/// Args:
///     population: Population to analyze
///     chromosome_idx: Index of chromosome to analyze
///
/// Returns:
///     dict: Dictionary with metric names as keys and values
#[pyfunction]
#[pyo3(name = "export_diversity_metrics")]
fn export_diversity_metrics_py(
    py: Python,
    population: &PyPopulation,
    chromosome_idx: usize,
) -> PyResult<Py<PyDict>> {
    let dict = PyDict::new(py);
    
    dict.set_item(
        "nucleotide_diversity",
        nucleotide_diversity(&population.inner, chromosome_idx),
    )?;
    dict.set_item(
        "tajimas_d",
        tajimas_d(&population.inner, chromosome_idx),
    )?;
    dict.set_item(
        "wattersons_theta",
        wattersons_theta(&population.inner, chromosome_idx),
    )?;
    dict.set_item(
        "haplotype_diversity",
        haplotype_diversity(&population.inner, chromosome_idx),
    )?;
    dict.set_item("generation", population.inner.generation())?;
    dict.set_item("population_size", population.inner.size())?;
    
    Ok(dict.into())
}

/// Export distance matrix to flat format for PyArrow.
///
/// Converts distance matrix to list of dictionaries with (i, j, distance) format,
/// suitable for PyArrow Table creation.
///
/// Args:
///     population: Population to analyze
///     chromosome_idx: Index of chromosome to analyze
///
/// Returns:
///     list: List of dicts with keys 'sequence_i', 'sequence_j', 'distance'
#[pyfunction]
#[pyo3(name = "export_distance_matrix")]
fn export_distance_matrix_py(
    py: Python,
    population: &PyPopulation,
    chromosome_idx: usize,
) -> PyResult<Py<PyList>> {
    let matrix = distance_matrix(&population.inner, chromosome_idx);
    let n = matrix.len();
    
    let mut records: Vec<Py<PyDict>> = Vec::new();
    for i in 0..n {
        for j in 0..n {
            let dict = PyDict::new(py);
            dict.set_item("sequence_i", i)?;
            dict.set_item("sequence_j", j)?;
            dict.set_item("distance", matrix[i][j])?;
            records.push(dict.into());
        }
    }
    
    Ok(PyList::new(py, records)?.into())
}

/// Export LD decay data to PyArrow-friendly format.
///
/// Args:
///     population: Population to analyze
///     chromosome_idx: Index of chromosome to analyze
///     haplotype_idx: Index of haplotype to analyze
///     max_distance: Maximum distance to analyze
///     bin_size: Size of distance bins
///
/// Returns:
///     list: List of dicts with keys 'distance', 'r_squared'
#[pyfunction]
#[pyo3(name = "export_ld_decay")]
fn export_ld_decay_py(
    py: Python,
    population: &PyPopulation,
    chromosome_idx: usize,
    haplotype_idx: usize,
    max_distance: usize,
    bin_size: usize,
) -> PyResult<Py<PyList>> {
    let decay = ld_decay(
        &population.inner,
        chromosome_idx,
        haplotype_idx,
        max_distance,
        bin_size,
    );
    
    let mut records: Vec<Py<PyDict>> = Vec::new();
    for (dist, r2) in decay {
        let dict = PyDict::new(py);
        dict.set_item("distance", dist)?;
        dict.set_item("r_squared", r2)?;
        records.push(dict.into());
    }
    
    Ok(PyList::new(py, records)?.into())
}

/// Register analysis functions with the Python module.
pub fn register_analysis_functions(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Diversity metrics
    m.add_function(wrap_pyfunction!(nucleotide_diversity_py, m)?)?;
    m.add_function(wrap_pyfunction!(tajimas_d_py, m)?)?;
    m.add_function(wrap_pyfunction!(wattersons_theta_py, m)?)?;
    m.add_function(wrap_pyfunction!(haplotype_diversity_py, m)?)?;
    
    // Linkage analysis
    m.add_function(wrap_pyfunction!(linkage_disequilibrium_py, m)?)?;
    m.add_function(wrap_pyfunction!(ld_decay_py, m)?)?;
    m.add_function(wrap_pyfunction!(haplotype_blocks_py, m)?)?;
    
    // Distance analysis
    m.add_function(wrap_pyfunction!(pairwise_distances_py, m)?)?;
    m.add_function(wrap_pyfunction!(distance_matrix_py, m)?)?;
    
    // Composition analysis
    m.add_function(wrap_pyfunction!(gc_content_py, m)?)?;
    m.add_function(wrap_pyfunction!(nucleotide_composition_py, m)?)?;
    
    // Polymorphism analysis
    m.add_function(wrap_pyfunction!(count_segregating_sites_py, m)?)?;
    
    // Export functions for PyArrow
    m.add_function(wrap_pyfunction!(export_diversity_metrics_py, m)?)?;
    m.add_function(wrap_pyfunction!(export_distance_matrix_py, m)?)?;
    m.add_function(wrap_pyfunction!(export_ld_decay_py, m)?)?;
    
    Ok(())
}
