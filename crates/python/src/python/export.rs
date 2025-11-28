//! Python bindings for data export functionality.

use pyo3::prelude::*;
use pyo3::types::PyDict;

use super::bindings::PyPopulation;

/// Export population sequences to FASTA format.
///
/// Args:
///     population: Population to export
///
/// Returns:
///     str: FASTA formatted sequences
#[pyfunction]
#[pyo3(name = "export_fasta")]
fn export_fasta_py(population: &PyPopulation) -> PyResult<String> {
    let mut fasta = String::new();

    for (i, ind) in population.inner.individuals().iter().enumerate() {
        let id = ind.id();

        // Export haplotype 1
        for (chr_idx, chr) in ind.haplotype1().chromosomes().iter().enumerate() {
            fasta.push_str(&format!(">{}|h1|chr{}\n", id, chr_idx));
            fasta.push_str(&chr.to_string());
            fasta.push('\n');
        }

        // Export haplotype 2
        for (chr_idx, chr) in ind.haplotype2().chromosomes().iter().enumerate() {
            fasta.push_str(&format!(">{}|h2|chr{}\n", id, chr_idx));
            fasta.push_str(&chr.to_string());
            fasta.push('\n');
        }
    }

    Ok(fasta)
}

/// Export population sequences to CSV format.
///
/// Args:
///     population: Population to export
///
/// Returns:
///     str: CSV formatted data with columns: individual_id, haplotype, chromosome, sequence
#[pyfunction]
#[pyo3(name = "export_csv")]
fn export_csv_py(population: &PyPopulation) -> PyResult<String> {
    let mut csv = String::from("individual_id,haplotype,chromosome,sequence\n");

    for ind in population.inner.individuals() {
        let id = ind.id();

        // Export haplotype 1
        for (chr_idx, chr) in ind.haplotype1().chromosomes().iter().enumerate() {
            csv.push_str(&format!("{},h1,{},{}\n", id, chr_idx, chr.to_string()));
        }

        // Export haplotype 2
        for (chr_idx, chr) in ind.haplotype2().chromosomes().iter().enumerate() {
            csv.push_str(&format!("{},h2,{},{}\n", id, chr_idx, chr.to_string()));
        }
    }

    Ok(csv)
}

/// Export population sequences to JSON format.
///
/// Args:
///     population: Population to export
///
/// Returns:
///     str: JSON formatted data
#[pyfunction]
#[pyo3(name = "export_json")]
fn export_json_py(population: &PyPopulation) -> PyResult<String> {
    use serde_json::json;

    let data: Vec<_> = population
        .inner
        .individuals()
        .iter()
        .map(|ind| {
            let h1_seqs: Vec<String> = ind
                .haplotype1()
                .chromosomes()
                .iter()
                .map(|chr| chr.to_string())
                .collect();

            let h2_seqs: Vec<String> = ind
                .haplotype2()
                .chromosomes()
                .iter()
                .map(|chr| chr.to_string())
                .collect();

            json!({
                "id": ind.id(),
                "fitness": ind.cached_fitness(),
                "haplotype1": h1_seqs,
                "haplotype2": h2_seqs,
            })
        })
        .collect();

    serde_json::to_string_pretty(&data)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("JSON serialization failed: {}", e)))
}

/// Export simulation metadata to JSON format.
///
/// Args:
///     info_dict: Dictionary from QueryBuilder.get_simulation_info()
///
/// Returns:
///     str: JSON formatted metadata
#[pyfunction]
#[pyo3(name = "export_metadata_json")]
fn export_metadata_json_py(info_dict: &Bound<'_, PyDict>) -> PyResult<String> {
    use serde_json::json;

    let obj = json!({
        "name": info_dict.get_item("sim_id")?.and_then(|v| v.extract::<String>().ok()),
        "population_size": info_dict.get_item("pop_size")?.and_then(|v| v.extract::<usize>().ok()),
        "generations": info_dict.get_item("num_generations")?.and_then(|v| v.extract::<usize>().ok()),
        "start_time": info_dict.get_item("start_time")?.and_then(|v| v.extract::<i64>().ok()),
        "end_time": info_dict.get_item("end_time")?.and_then(|v| v.extract::<Option<i64>>().ok()).flatten(),
        "mutation_rate": info_dict.get_item("mutation_rate")?.and_then(|v| v.extract::<f64>().ok()),
        "recombination_rate": info_dict.get_item("recombination_rate")?.and_then(|v| v.extract::<f64>().ok()),
    });

    serde_json::to_string_pretty(&obj)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("JSON serialization failed: {}", e)))
}

/// Register export functions with the Python module.
pub fn register_export_functions(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(export_fasta_py, m)?)?;
    m.add_function(wrap_pyfunction!(export_csv_py, m)?)?;
    m.add_function(wrap_pyfunction!(export_json_py, m)?)?;
    m.add_function(wrap_pyfunction!(export_metadata_json_py, m)?)?;

    Ok(())
}
