use anyhow::{Context, Result};
use centrevo_sim::storage::QueryBuilder;
use std::path::PathBuf;

use crate::utils::sequence_from_indices;

pub fn export_data(
    database: &PathBuf,
    generation: usize,
    format: &str,
    output: Option<&PathBuf>,
    data_type: &str,
) -> Result<()> {
    let name = database
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("simulation");

    println!("📤 Exporting data from simulation '{name}'");
    println!("Generation: {generation}, Format: {format}, Type: {data_type}");

    let query = QueryBuilder::new(database).context("Failed to open database")?;

    match data_type {
        "sequences" => export_sequences(&query, generation, format, output)?,
        "metadata" => export_metadata(&query, name, format, output)?,
        "fitness" => export_fitness(&query, format, output)?,
        _ => anyhow::bail!("Unknown data type '{data_type}'. Use: sequences, metadata, or fitness"),
    }

    if let Some(path) = output {
        println!("✓ Data exported to: {}", path.display());
    } else {
        println!("\n✓ Export complete");
    }

    Ok(())
}

fn export_sequences(
    query: &QueryBuilder,
    generation: usize,
    format: &str,
    output: Option<&PathBuf>,
) -> Result<()> {
    let snapshots = query
        .get_generation(generation)
        .context("Failed to load generation")?;

    if snapshots.is_empty() {
        anyhow::bail!("No data found for generation {generation}");
    }

    let mut content = String::new();

    match format {
        "csv" => {
            content.push_str("individual_id,haplotype,chromosome,sequence\n");
            for (id, snap) in &snapshots {
                let ind_id = format!("ind_{id}");
                // Haplotype 1
                let seq1 = sequence_from_indices(snap.haplotype1_seq.clone());
                content.push_str(&format!("{},h1,0,{}\n", ind_id, seq1));

                // Haplotype 2
                let seq2 = sequence_from_indices(snap.haplotype2_seq.clone());
                content.push_str(&format!("{},h2,0,{}\n", ind_id, seq2));
            }
        }
        "fasta" => {
            for (id, snap) in &snapshots {
                let ind_id = format!("ind_{id}");
                // Haplotype 1
                let seq1 = sequence_from_indices(snap.haplotype1_seq.clone());
                content.push_str(&format!(">{}|h1|chr0\n{}\n", ind_id, seq1));

                // Haplotype 2
                let seq2 = sequence_from_indices(snap.haplotype2_seq.clone());
                content.push_str(&format!(">{}|h2|chr0\n{}\n", ind_id, seq2));
            }
        }
        "json" => {
            use serde_json::json;
            let data: Vec<_> = snapshots
                .iter()
                .map(|(id, snap)| {
                    let ind_id = format!("ind_{id}");
                    let seq1 = sequence_from_indices(snap.haplotype1_seq.clone());
                    let seq2 = sequence_from_indices(snap.haplotype2_seq.clone());
                    json!({
                        "id": ind_id,
                        "fitness": snap.fitness,
                        "haplotype1": [seq1.to_string()],
                        "haplotype2": [seq2.to_string()],
                    })
                })
                .collect();
            content = serde_json::to_string_pretty(&data)?;
        }
        _ => anyhow::bail!("Unknown format '{format}'. Use: csv, fasta, or json"),
    }

    if let Some(path) = output {
        std::fs::write(path, content)?;
    } else {
        println!("{content}");
    }

    Ok(())
}

fn export_metadata(
    query: &QueryBuilder,
    name: &str,
    format: &str,
    output: Option<&PathBuf>,
) -> Result<()> {
    let config = query.get_full_config()?;
    let pop_size = config.execution.population_size;
    let generations = config.execution.total_generations;
    // We don't have start_time in config atm, use current or omit

    let content = match format {
        "json" => {
            use serde_json::json;
            serde_json::to_string_pretty(&json!({
                "name": name,
                "population_size": pop_size,
                "generations": generations,
                "parameters": config
            }))?
        }
        "csv" => {
            format!(
                "key,value\nname,{}\npopulation_size,{}\ngenerations,{}\n",
                name, pop_size, generations
            )
        }
        _ => anyhow::bail!("Format '{format}' not supported for metadata. Use: json or csv"),
    };

    if let Some(path) = output {
        std::fs::write(path, content)?;
    } else {
        println!("{content}");
    }

    Ok(())
}

fn export_fitness(query: &QueryBuilder, format: &str, output: Option<&PathBuf>) -> Result<()> {
    let history = query.get_fitness_history()?;

    if history.is_empty() {
        anyhow::bail!("No fitness data found");
    }

    let content = match format {
        "csv" => {
            let mut csv =
                String::from("generation,mean_fitness,std_fitness,min_fitness,max_fitness\n");
            for (generation, stats) in &history {
                csv.push_str(&format!(
                    "{},{},{},{},{}\n",
                    generation, stats.mean, stats.std, stats.min, stats.max
                ));
            }
            csv
        }
        "json" => {
            use serde_json::json;
            let data: Vec<_> = history
                .iter()
                .map(|(generation, stats)| {
                    json!({
                        "generation": generation,
                        "mean": stats.mean,
                        "std_dev": stats.std,
                        "min": stats.min,
                        "max": stats.max,
                    })
                })
                .collect();
            serde_json::to_string_pretty(&data)?
        }
        _ => anyhow::bail!("Format '{format}' not supported for fitness. Use: csv or json"),
    };

    if let Some(path) = output {
        std::fs::write(path, content)?;
    } else {
        println!("{content}");
    }

    Ok(())
}
