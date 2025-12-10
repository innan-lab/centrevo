use anyhow::{Context, Result};
use centrevo_sim::storage::QueryBuilder;
use std::path::PathBuf;

use crate::utils::sequence_from_indices;

pub fn export_data(
    database: &PathBuf,
    name: &str,
    generation: usize,
    format: &str,
    output: Option<&PathBuf>,
    data_type: &str,
) -> Result<()> {
    println!("📤 Exporting data from simulation '{name}'");
    println!("Generation: {generation}, Format: {format}, Type: {data_type}");

    let query = QueryBuilder::new(database).context("Failed to open database")?;

    match data_type {
        "sequences" => export_sequences(&query, name, generation, format, output)?,
        "metadata" => export_metadata(&query, name, format, output)?,
        "fitness" => export_fitness(&query, name, format, output)?,
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
    name: &str,
    generation: usize,
    format: &str,
    output: Option<&PathBuf>,
) -> Result<()> {
    let snapshots = query
        .get_generation(name, generation)
        .context("Failed to load generation")?;

    if snapshots.is_empty() {
        anyhow::bail!("No data found for generation {generation}");
    }

    let mut content = String::new();

    match format {
        "csv" => {
            content.push_str("individual_id,haplotype,chromosome,sequence\n");
            for snap in &snapshots {
                // Haplotype 1
                let seq1 = sequence_from_indices(snap.haplotype1_seq.clone());
                content.push_str(&format!("{},h1,0,{}\n", snap.individual_id, seq1));

                // Haplotype 2
                let seq2 = sequence_from_indices(snap.haplotype2_seq.clone());
                content.push_str(&format!("{},h2,0,{}\n", snap.individual_id, seq2));
            }
        }
        "fasta" => {
            for snap in &snapshots {
                // Haplotype 1
                let seq1 = sequence_from_indices(snap.haplotype1_seq.clone());
                content.push_str(&format!(">{}|h1|chr0\n{}\n", snap.individual_id, seq1));

                // Haplotype 2
                let seq2 = sequence_from_indices(snap.haplotype2_seq.clone());
                content.push_str(&format!(">{}|h2|chr0\n{}\n", snap.individual_id, seq2));
            }
        }
        "json" => {
            use serde_json::json;
            let data: Vec<_> = snapshots
                .iter()
                .map(|snap| {
                    let seq1 = sequence_from_indices(snap.haplotype1_seq.clone());
                    let seq2 = sequence_from_indices(snap.haplotype2_seq.clone());
                    json!({
                        "id": snap.individual_id,
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
    let info = query.get_simulation_info(name)?;

    let content = match format {
        "json" => {
            use serde_json::json;
            serde_json::to_string_pretty(&json!({
                "name": name,
                "population_size": info.pop_size,
                "generations": info.num_generations,
                "start_time": info.start_time,
                "parameters": serde_json::from_str::<serde_json::Value>(&info.parameters_json)?
            }))?
        }
        "csv" => {
            format!(
                "key,value\nname,{}\npopulation_size,{}\ngenerations,{}\nstart_time,{}\n",
                name, info.pop_size, info.num_generations, info.start_time
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

fn export_fitness(
    query: &QueryBuilder,
    name: &str,
    format: &str,
    output: Option<&PathBuf>,
) -> Result<()> {
    let history = query.get_fitness_history(name)?;

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
