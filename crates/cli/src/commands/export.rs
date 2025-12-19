use anyhow::{Context, Result};
use centrevo_sim::storage::QueryBuilder;
use std::path::PathBuf;

use crate::utils::sequence_from_indices;

pub fn export_data(
    database: &PathBuf,
    generations_arg: Option<String>,
    format: &str,
    output: Option<&PathBuf>,
    data_type: &str,
    ru_delimiter: &str,
    hor_delimiter: &str,
    chr_delimiter: &str,
) -> Result<()> {
    let name = database
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("simulation");

    let query = QueryBuilder::new(database).context("Failed to open database")?;

    // --- Validation & Parsing ---

    // 1. Output required for Sequence + FASTA (due to sidecar BED)
    if data_type == "sequence" && format == "fasta" && output.is_none() {
        anyhow::bail!(
            "--output is required for 'sequence' data in 'fasta' format (to generate sidecar BED file)"
        );
    }

    // 2. Generation required for Sequence and Fitness
    if (data_type == "sequence" || data_type == "fitness") && generations_arg.is_none() {
        anyhow::bail!("--generations is required for '{data_type}' data export");
    }

    // 3. Metadata doesn't need generations, handling it first
    if data_type == "metadata" {
        return export_metadata(&query, name, format, output);
    }

    // Parse generations
    let available_gens = get_available_generations(&query)?;
    let selected_gens = parse_generations(generations_arg.as_deref().unwrap(), &available_gens)?;

    if selected_gens.is_empty() {
        anyhow::bail!("No matching generations found in the database.");
    }

    // 4. Directory Output Validation
    let is_multi = selected_gens.len() > 1;
    if is_multi {
        if let Some(path) = output {
            if path.exists() && !path.is_dir() {
                anyhow::bail!(
                    "Output path must be a directory when exporting multiple generations."
                );
            }
            if !path.exists() {
                std::fs::create_dir_all(path).context("Failed to create output directory")?;
            }
        } else {
            // Stdout not allowed for multi-generation (unless maybe we format it with headers? But user agreed to directory)
            anyhow::bail!("--output directory is required when exporting multiple generations.");
        }
    }

    println!("📤 Exporting data from simulation '{name}'");
    println!(
        "Generations: {} selected, Format: {format}, Type: {data_type}",
        selected_gens.len()
    );

    // Loop through generations
    for g in selected_gens {
        // determine output path for this generation
        let gen_output = if is_multi {
            let dir = output.unwrap(); // Validated above
            let ext = match format {
                "fasta" => "fasta",
                "csv" => "csv",
                "json" => "json",
                _ => "txt",
            };
            Some(dir.join(format!("gen_{g}.{ext}")))
        } else {
            output.map(|p| p.to_path_buf())
        };

        match data_type {
            "sequences" | "sequence" => export_sequences(
                &query,
                g,
                format,
                gen_output.as_ref(),
                ru_delimiter,
                hor_delimiter,
                chr_delimiter,
            )?,
            "fitness" => export_fitness(&query, g, format, gen_output.as_ref())?,
            _ => anyhow::bail!(
                "Unknown data type '{data_type}'. Use: sequences, metadata, or fitness"
            ),
        }

        if is_multi {
            println!("✓ Exported generation {g}");
        }
    }

    if !is_multi {
        if let Some(path) = output {
            println!("✓ Data exported to: {}", path.display());
        } else {
            println!("\n✓ Export complete");
        }
    } else {
        println!(
            "\n✓ All requested generations exported to: {}",
            output.unwrap().display()
        );
    }

    Ok(())
}

fn get_available_generations(query: &QueryBuilder) -> Result<Vec<usize>> {
    let history = query.get_fitness_history()?;
    let mut gens: Vec<usize> = history.iter().map(|(g, _)| *g).collect();
    if gens.is_empty() {
        // Fallback: check if gen 0 exists if history is empty (e.g. no stats recorded yet?)
        // But usually history has entries. If empty, maybe only gen 0 created and not finalized?
        // query.get_generation(0) check?
        if !query.get_generation(0)?.is_empty() {
            gens.push(0);
        }
    }
    gens.sort();
    Ok(gens)
}

fn parse_generations(arg: &str, available: &[usize]) -> Result<Vec<usize>> {
    if arg == "all" {
        return Ok(available.to_vec());
    }

    let mut selected = std::collections::HashSet::new();

    for part in arg.split(',') {
        let part = part.trim();
        if part.contains("..") {
            let ranges: Vec<&str> = part.split("..").collect();
            if ranges.len() != 2 {
                anyhow::bail!("Invalid range syntax: {part}");
            }
            let start: usize = ranges[0].parse().context("Invalid range start")?;
            let end: usize = ranges[1].parse().context("Invalid range end")?;

            for g in available {
                if *g >= start && *g <= end {
                    selected.insert(*g);
                }
            }
        } else {
            let g: usize = part.parse().context("Invalid generation number")?;
            if available.contains(&g) {
                selected.insert(g);
            } else {
                // Determine behaviour: strict or lax?
                // Plan implied lax filtering "Filters against actual", but user might expect error if specific missing.
                // Let's warn but not fail? Or fail?
                // For now, strict: if user explicitly requested '100' and it's missing, maybe warn?
                // Logic above: "Filters against actual recorded".
                // I will skip if missing, and check empty at end.
                eprintln!("Warning: Generation {g} not found in database.");
            }
        }
    }

    let mut result: Vec<usize> = selected.into_iter().collect();
    result.sort();
    Ok(result)
}

fn export_sequences(
    query: &QueryBuilder,
    generation: usize,
    format: &str,
    output: Option<&PathBuf>,
    ru_delimiter: &str,
    hor_delimiter: &str,
    chr_delimiter: &str,
) -> Result<()> {
    let snapshots = query
        .get_generation(generation)
        .context("Failed to load generation")?;

    if snapshots.is_empty() {
        anyhow::bail!("No data found for generation {generation}");
    }

    let config = query.get_full_config().context("Failed to load config")?;
    let codec = config.execution.codec;

    match format {
        "csv" => {
            let mut content = String::from("individual_id,haplotype,chromosome,sequence\n");
            for (id, snap) in &snapshots {
                let ind_id = format!("ind_{id}");
                let seq1 = sequence_from_indices(snap.haplotype1_seq.clone());
                content.push_str(&format!("{ind_id},h1,0,{seq1}\n"));

                let seq2 = sequence_from_indices(snap.haplotype2_seq.clone());
                content.push_str(&format!("{ind_id},h2,0,{seq2}\n"));
            }
            if let Some(path) = output {
                std::fs::write(path, content)?;
            } else {
                print!("{content}"); // print! avoids extra newline
            }
        }
        "fasta" => {
            let path = output.unwrap(); // Validated
            let bed_path = path.with_extension("bed");

            let mut fasta_content = String::new();
            let mut bed_content = String::new();

            for (id, snap) in &snapshots {
                let ind = snap
                    .to_individual(&id.to_string(), &codec)
                    .map_err(|e| anyhow::anyhow!(e))?;

                // Haplotype 1
                if let Some(chr) = ind.haplotype1().get(0) {
                    let header = format!("ind_{id}|h1|chr0");
                    fasta_content.push_str(&format!(">{header}\n{}\n", chr.sequence()));

                    let map = chr.map();
                    for i in 0..map.num_rus() {
                        if let Some((start, end)) = map.get_ru_interval(i) {
                            bed_content
                                .push_str(&format!("{header}\t{start}\t{end}\tRU_{i}\t0\t+\n"));
                        }
                    }
                }

                // Haplotype 2
                if let Some(chr) = ind.haplotype2().get(0) {
                    let header = format!("ind_{id}|h2|chr0");
                    fasta_content.push_str(&format!(">{header}\n{}\n", chr.sequence()));

                    let map = chr.map();
                    for i in 0..map.num_rus() {
                        if let Some((start, end)) = map.get_ru_interval(i) {
                            bed_content
                                .push_str(&format!("{header}\t{start}\t{end}\tRU_{i}\t0\t+\n"));
                        }
                    }
                }
            }

            std::fs::write(path, fasta_content)?;
            std::fs::write(&bed_path, bed_content)?;
            // Only print sidecar msg if not checking logs too verbosely
            // println!("✓ Sidecar BED file created: {}", bed_path.display());
        }
        "formatted-text" => {
            let mut content = String::new();

            for (id, snap) in &snapshots {
                let ind = snap
                    .to_individual(&id.to_string(), &codec)
                    .map_err(|e| anyhow::anyhow!(e))?;

                // Helper to format a chromosome
                let format_chr = |chr: &centrevo_sim::genome::Chromosome| -> String {
                    let seq_str = chr.sequence().to_string();
                    let map = chr.map();
                    let mut formatted = String::new();

                    // Iterate RUs
                    for i in 0..map.num_rus() {
                        if i > 0 {
                            // Check if we crossed an HOR boundary
                            let is_hor_start = map
                                .find_hor_index(i)
                                .map(|hor_idx| {
                                    map.get_hor_interval(hor_idx)
                                        .map(|(s, _)| s == map.get_ru_interval(i).unwrap().0)
                                        .unwrap_or(false)
                                })
                                .unwrap_or(false);

                            if is_hor_start {
                                formatted.push_str(hor_delimiter);
                            } else {
                                formatted.push_str(ru_delimiter);
                            }
                        }

                        if let Some((start, end)) = map.get_ru_interval(i) {
                            formatted.push_str(&seq_str[start..end]);
                        }
                    }
                    formatted
                };

                let h1_str = ind
                    .haplotype1()
                    .iter()
                    .map(format_chr)
                    .collect::<Vec<_>>()
                    .join(hor_delimiter);

                let h2_str = ind
                    .haplotype2()
                    .iter()
                    .map(format_chr)
                    .collect::<Vec<_>>()
                    .join(hor_delimiter);

                content.push_str(&format!("{h1_str}\t{h2_str}{chr_delimiter}"));
            }

            if let Some(path) = output {
                std::fs::write(path, content)?;
            } else {
                print!("{content}");
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
            let content = serde_json::to_string_pretty(&data)?;
            if let Some(path) = output {
                std::fs::write(path, content)?;
            } else {
                println!("{content}");
            }
        }
        _ => anyhow::bail!("Unknown format '{format}'. Use: csv, fasta, formatted-text, or json"),
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
                "key,value\nname,{name}\npopulation_size,{pop_size}\ngenerations,{generations}\n",
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
    generation: usize,
    format: &str,
    output: Option<&PathBuf>,
) -> Result<()> {
    // Export individual fitness for specific generation
    let snapshots = query
        .get_generation(generation)
        .context("Failed to load generation")?;
    if snapshots.is_empty() {
        anyhow::bail!("No data for generation {generation}");
    }

    match format {
        "csv" => {
            let mut content = String::from("individual_id,fitness_log,fitness_linear\n");
            for (id, snap) in snapshots {
                if let Some(log_fitness) = snap.fitness {
                    let linear = log_fitness.exp();
                    content.push_str(&format!("ind_{id},{log_fitness},{linear}\n"));
                }
            }
            if let Some(path) = output {
                std::fs::write(path, content)?;
            } else {
                println!("{content}");
            }
        }
        _ => anyhow::bail!("Format '{format}' not supported for individual fitness. Use: csv"),
    }

    Ok(())
}
