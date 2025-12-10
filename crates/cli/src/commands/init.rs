use anyhow::{Context, Result};
use centrevo_sim::base::Nucleotide;
use centrevo_sim::genome::{Chromosome, Haplotype, Individual};
use centrevo_sim::simulation::{
    FitnessConfig, MutationConfig, Population, RecombinationConfig, SimulationConfig,
    UniformRepeatStructure,
};
use centrevo_sim::storage::{Recorder, RecordingStrategy, SimulationSnapshot};
use std::path::PathBuf;

use crate::printing::print_parameters;

#[allow(clippy::too_many_arguments)]
pub fn init_simulation(
    name: &str,
    output: &PathBuf,
    population_size: usize,
    generations: usize,
    ru_length: usize,
    rus_per_hor: usize,
    hors_per_chr: usize,
    mutation_rate: f64,
    recomb_rate: f64,
    crossover_prob: f64,
    gc_extension_prob: f64,
    homology_strength: f64,
    search_window: usize,
    record_every: usize,
    seed: Option<u64>,
) -> Result<()> {
    println!("🧬 Centrevo - Centromeric Evolution Simulator");
    println!("============================================\n");
    println!("Initializing simulation: {name}");

    // Create configurations
    let structure =
        UniformRepeatStructure::new(Nucleotide::A, ru_length, rus_per_hor, hors_per_chr, 1);

    let config = SimulationConfig::new(population_size, generations, seed);

    // Create mutation and recombination configs
    let mutation = MutationConfig::uniform(mutation_rate)
        .map_err(|e| anyhow::anyhow!("Failed to create mutation configuration: {e}"))?;

    let recomb_params = centrevo_sim::evolution::RecombinationModel::builder()
        .break_prob(recomb_rate)
        .crossover_prob(crossover_prob)
        .gc_extension_prob(gc_extension_prob)
        .homology_strength(homology_strength)
        .search_window(search_window)
        .build()
        .map_err(|e| anyhow::anyhow!("Failed to create recombination model: {e}"))?;
    let recombination = RecombinationConfig::new(recomb_params);

    let fitness = FitnessConfig::neutral();

    println!("\nConfiguration:");
    print_parameters(
        &config,
        Some(&structure),
        &mutation,
        &recombination,
        &fitness,
    );

    // Create initial population (Generation 0)
    // We use the initial seed (if provided) or a random one to generate Gen 0.
    // Note: Since Gen 0 is uniform, the seed actually doesn't matter for the content,
    // but we use it for consistency if we add random initialization later.
    println!("\nCreating initial population (Generation 0)...");
    let population = create_initial_population(population_size, &structure);
    println!("✓ Created {} individuals", population.size());

    // Setup database recorder
    println!("\nSetting up database...");
    let mut recorder = Recorder::new(output, name, RecordingStrategy::EveryN(record_every))
        .context("Failed to create recorder")?;

    // Record full configuration
    let snapshot = SimulationSnapshot {
        structure: structure.clone(),
        mutation: mutation.clone(),
        recombination: recombination.clone(),
        fitness: fitness.clone(),
        config: config.clone(),
    };

    recorder
        .record_full_config(&snapshot)
        .context("Failed to record configuration")?;

    // Record initial generation
    recorder
        .record_generation(&population, 0)
        .context("Failed to record initial generation")?;

    println!("✓ Database created: {}", output.display());
    println!("\nSimulation initialized successfully!");
    println!("  Name: {name}");
    println!("  Population size: {population_size}");
    println!("  Generations: {generations}");
    println!("\n💡 Use 'centrevo run -N {name}' to start the simulation");

    Ok(())
}

fn create_initial_population(size: usize, structure: &UniformRepeatStructure) -> Population {
    let mut individuals = Vec::with_capacity(size);

    for i in 0..size {
        let chr1 = Chromosome::uniform(
            format!("ind{i}_h1_chr1"),
            structure.init_base,
            structure.ru_length,
            structure.rus_per_hor,
            structure.hors_per_chr,
        );

        let chr2 = Chromosome::uniform(
            format!("ind{i}_h2_chr1"),
            structure.init_base,
            structure.ru_length,
            structure.rus_per_hor,
            structure.hors_per_chr,
        );

        let h1 = Haplotype::from_chromosomes(vec![chr1]);
        let h2 = Haplotype::from_chromosomes(vec![chr2]);

        individuals.push(Individual::new(format!("ind{i}"), h1, h2));
    }

    Population::new("initial_pop", individuals)
}

#[cfg(test)]
mod tests {
    use super::*;
    use centrevo_sim::base::Nucleotide;

    #[test]
    fn test_create_initial_population_regression() {
        // Use small structure for testing
        let structure = UniformRepeatStructure::new(Nucleotide::A, 10, 5, 2, 1);
        let pop = create_initial_population(10, &structure);

        assert_eq!(pop.size(), 10);
        let ind = pop.get(0).unwrap();
        // Check chromosome length: 10 * 5 * 2 = 100
        assert_eq!(ind.haplotype1().len(), 1);
        assert_eq!(ind.haplotype1().total_length(), 100);
    }
}
