use anyhow::{Context, Result};
use centrevo_sim::base::FitnessValue;
use centrevo_sim::genome::{Chromosome, Haplotype, Individual};
use centrevo_sim::simulation::Population;
use centrevo_sim::storage::QueryBuilder;
use std::path::PathBuf;

use crate::utils::sequence_from_indices;

pub fn analyze_data(
    database: &PathBuf,
    name: &str,
    generation: usize,
    chromosome: usize,
    format: &str,
    output: Option<&PathBuf>,
) -> Result<()> {
    use centrevo_analysis::analysis::{
        composition::gc_content, haplotype_diversity, nucleotide_diversity,
        polymorphism::count_segregating_sites, tajimas_d, wattersons_theta,
    };

    println!("🔬 Analyzing simulation '{name}'");
    println!("Generation: {generation}, Chromosome: {chromosome}");

    let query = QueryBuilder::new(database).context("Failed to open database")?;
    let snapshots = query
        .get_generation(name, generation)
        .context("Failed to load generation")?;

    if snapshots.is_empty() {
        anyhow::bail!("No data found for generation {generation}");
    }

    // Convert snapshots to individuals for analysis
    let mut individuals = Vec::new();

    for snap in snapshots {
        let seq1 = sequence_from_indices(snap.haplotype1_seq.clone());
        let seq2 = sequence_from_indices(snap.haplotype2_seq.clone());

        // Assume uniform structure for analysis reconstruction
        // TODO: Load actual structure from database
        let ru_len = 171;
        let rus_per_hor = 12;
        let hor_len = ru_len * rus_per_hor;
        // Avoid division by zero if length is 0 (though unlikely for valid sim)
        let hors_per_chr = if hor_len > 0 { seq1.len() / hor_len } else { 0 };

        let map =
            centrevo_sim::genome::repeat_map::RepeatMap::uniform(ru_len, rus_per_hor, hors_per_chr);

        let chr1 = Chromosome::new(snap.haplotype1_chr_id, seq1, map.clone());
        let chr2 = Chromosome::new(snap.haplotype2_chr_id, seq2, map);

        let h1 = Haplotype::from_chromosomes(vec![chr1]);
        let h2 = Haplotype::from_chromosomes(vec![chr2]);

        let mut ind = Individual::new(snap.individual_id, h1, h2);
        if let Some(f) = snap.fitness {
            ind.set_cached_fitness(FitnessValue::new(f));
        }
        individuals.push(ind);
    }

    let population = Population::new(format!("{name}_gen{generation}"), individuals);
    let seq_len = population.individuals()[0]
        .haplotype1()
        .get(chromosome)
        .map(|chr| chr.sequence().len())
        .unwrap_or(0);

    // Calculate metrics
    let pi = nucleotide_diversity(&population, chromosome);
    let tajima = tajimas_d(&population, chromosome);
    let theta_w = wattersons_theta(&population, chromosome);
    let hap_div = haplotype_diversity(&population, chromosome);
    let seg_sites = count_segregating_sites(&population, chromosome, 0)
        + count_segregating_sites(&population, chromosome, 1);

    // Calculate GC content at population level
    let gc = gc_content(&population, None, None, None);

    let content = match format {
        "pretty" => {
            format!(
                "\n📊 Population Genetics Summary\n\
                 ================================\n\
                 Population size: {}\n\
                 Sequences analyzed: {} (2n)\n\
                 Sequence length: {} bp\n\
                 \n\
                 Diversity Metrics:\n\
                 ------------------\n\
                 Nucleotide diversity (π): {:.6}\n\
                 Watterson's theta (θ_W): {:.6}\n\
                 Tajima's D: {:.4}\n\
                 Haplotype diversity: {:.4}\n\
                 \n\
                 Polymorphism:\n\
                 -------------\n\
                 Segregating sites: {}\n\
                 \n\
                 Composition:\n\
                 ------------\n\
                 Mean GC content: {:.2}%\n",
                population.size(),
                population.size() * 2,
                seq_len,
                pi,
                theta_w,
                tajima,
                hap_div,
                seg_sites,
                gc * 100.0,
            )
        }
        "json" => {
            use serde_json::json;
            serde_json::to_string_pretty(&json!({
                "simulation": name,
                "generation": generation,
                "chromosome": chromosome,
                "population_size": population.size(),
                "sequence_count": population.size() * 2,
                "sequence_length": seq_len,
                "diversity": {
                    "pi": pi,
                    "theta_w": theta_w,
                    "tajima_d": tajima,
                    "haplotype_diversity": hap_div,
                },
                "polymorphism": {
                    "segregating_sites": seg_sites,
                },
                "composition": {
                    "mean_gc_content": gc,
                }
            }))?
        }
        _ => anyhow::bail!("Unknown format '{format}'. Use: pretty or json"),
    };

    if let Some(path) = output {
        std::fs::write(path, content)?;
    } else {
        println!("{content}");
    }

    Ok(())
}
