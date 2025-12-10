use anyhow::{Context, Result};
use centrevo_sim::base::Nucleotide;
use centrevo_sim::evolution::{IndelModel, SubstitutionModel};
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
    chrs_per_hap: usize,
    mutation_rate: Option<f64>,
    rate_ac: Option<f64>,
    rate_ag: Option<f64>,
    rate_at: Option<f64>,
    rate_cg: Option<f64>,
    rate_ct: Option<f64>,
    rate_gt: Option<f64>,
    indel_ins_rate: f64,
    indel_del_rate: f64,
    indel_length_p: f64,
    recomb_rate: f64,
    crossover_prob: f64,
    gc_extension_prob: f64,
    homology_strength: f64,
    search_window: usize,
    fit_gc_opt: Option<f64>,
    fit_gc_conc: Option<f64>,
    fit_len_opt: Option<usize>,
    fit_len_std: Option<f64>,
    fit_seq_sim: Option<f64>,
    fit_len_sim: Option<f64>,
    record_every: usize,
    seed: Option<u64>,
) -> Result<()> {
    println!("🧬 Centrevo - Centromeric Evolution Simulator");
    println!("============================================\n");
    println!("Initializing simulation: {name}");

    // Create configurations
    let structure = UniformRepeatStructure::new(
        Nucleotide::A,
        ru_length,
        rus_per_hor,
        hors_per_chr,
        chrs_per_hap,
    );

    let config = SimulationConfig::new(population_size, generations, seed);

    // Create mutation configuration
    // Check for advanced mutation rates
    let substitution = if let (Some(ac), Some(ag), Some(at), Some(cg), Some(ct), Some(gt)) =
        (rate_ac, rate_ag, rate_at, rate_cg, rate_ct, rate_gt)
    {
        // Construct matrix (upper triangle provided)
        // A=0, C=1, G=2, T=3
        let matrix = [
            [0.0, ac, ag, at],
            [ac, 0.0, cg, ct],
            [ag, cg, 0.0, gt],
            [at, ct, gt, 0.0],
        ];
        SubstitutionModel::new(matrix)
            .map_err(|e| anyhow::anyhow!("Failed to create substitution model: {e}"))?
    } else if let Some(rate) = mutation_rate {
        // Uniform
        SubstitutionModel::uniform(rate)
            .map_err(|e| anyhow::anyhow!("Failed to create substitution model: {e}"))?
    } else {
        // Fallback or Error?
        // Since main.rs defines mutation_rate default, this branch might be unreachable if logic is tight,
        // BUT main.rs default implies mutation_rate is strictly set unless user overrides?
        // Wait, I changed mutation_rate to Option<f64> but I kept default_value="1e-5".
        // Clap will parse the default value into Some(1e-5) if not provided.
        // So this else should strictly be unreachable unless I messed up Clap config.
        // However, if the user provides explicit rates, 'mutation_rate' encounters conflict.
        // So mutation_rate being None means explicit rates MUST be provided.
        // But what if the user provides 5 explicit rates and forgets 1?
        // My main.rs didn't enforce "All or Nothing" for explicit rates group.
        // So I should check that here.
        if rate_ac.is_some()
            || rate_ag.is_some()
            || rate_at.is_some()
            || rate_cg.is_some()
            || rate_ct.is_some()
            || rate_gt.is_some()
        {
            anyhow::bail!(
                "When using specific mutation rates, ALL 6 rates (ac, ag, at, cg, ct, gt) must be provided."
            );
        }
        // If we really get here with everything None, it's weird. defaulting to 1e-5.
        SubstitutionModel::uniform(1e-5)
            .map_err(|e| anyhow::anyhow!("Failed to create substitution model: {e}"))?
    };

    let mut mutation = MutationConfig::new(substitution);

    // Add Indels if enabled
    if indel_ins_rate > 0.0 || indel_del_rate > 0.0 {
        let indel_model = IndelModel::new(indel_ins_rate, indel_del_rate, indel_length_p)
            .map_err(|e| anyhow::anyhow!("Failed to create indel model: {e}"))?;
        mutation.indel = Some(indel_model);
    }

    let recomb_params = centrevo_sim::evolution::RecombinationModel::builder()
        .break_prob(recomb_rate)
        .crossover_prob(crossover_prob)
        .gc_extension_prob(gc_extension_prob)
        .homology_strength(homology_strength)
        .search_window(search_window)
        .build()
        .map_err(|e| anyhow::anyhow!("Failed to create recombination model: {e}"))?;
    let recombination = RecombinationConfig::new(recomb_params);

    // Build Fitness Config
    // Since neutral() returns FitnessConfig, not a builder... wait.
    // The builder API: FitnessConfigBuilder::<BuilderEmpty>::neutral() returns FitnessConfig directly.
    // But FitnessConfigBuilder::<BuilderEmpty>::with_gc_content(...) returns a Builder.
    // I need to chain conditional `with_` calls.
    // But the builder types change! This is a type-state builder.
    // FitnessConfigBuilder<BuilderEmpty>.with_gc() -> FitnessConfigBuilder<BuilderInitialized>
    // Ah, this makes conditional chaining hard in Rust without boxing or fancy traits.
    // Actually, looking at configs.rs:
    // BuilderInitialized has `with_gc_content` too.
    // So I can transition to Initialized state once, then chain.
    // But how to start?
    // BuilderEmpty only has `neutral` (returns Config) or `with_*` (returns BuilderInitialized).
    // I can't start with an "empty initialized builder".
    // I might need to painstakingly check each option.

    // Better approach: Reconstruct the FitnessConfig manually using its `new` constructor since I have all options.
    // The Builder is useful for fluent API but maybe overkill here if I have all parts.
    // FitnessConfig::new(gc, length, sim, len_sim)
    // Let's check FitnessConfig::new signature in configs.rs.
    // pub fn new(gc, length, seq_sim, length_sim) -> Self.
    // Yes! That's much easier than fighting the type-state builder in dynamic code.

    let gc_fitness = if let (Some(opt), Some(conc)) = (fit_gc_opt, fit_gc_conc) {
        Some(
            centrevo_sim::evolution::GCContentFitness::new(opt, conc)
                .map_err(|e| anyhow::anyhow!("Invalid GC Fitness: {e}"))?,
        )
    } else {
        None
    };

    let len_fitness = if let (Some(opt), Some(std)) = (fit_len_opt, fit_len_std) {
        Some(
            centrevo_sim::evolution::LengthFitness::new(opt, std)
                .map_err(|e| anyhow::anyhow!("Invalid Length Fitness: {e}"))?,
        )
    } else {
        None
    };

    let seq_sim_fitness = if let Some(shape) = fit_seq_sim {
        Some(
            centrevo_sim::evolution::SequenceSimilarityFitness::new(shape)
                .map_err(|e| anyhow::anyhow!("Invalid Seq Sim Fitness: {e}"))?,
        )
    } else {
        None
    };

    let len_sim_fitness = if let Some(shape) = fit_len_sim {
        Some(
            centrevo_sim::evolution::LengthSimilarityFitness::new(shape)
                .map_err(|e| anyhow::anyhow!("Invalid Len Sim Fitness: {e}"))?,
        )
    } else {
        None
    };

    let fitness = FitnessConfig::new(gc_fitness, len_fitness, seq_sim_fitness, len_sim_fitness);

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
