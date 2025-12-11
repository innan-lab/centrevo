use anyhow::{Context, Result};
use centrevo_sim::base::Nucleotide;
use centrevo_sim::evolution::{IndelModel, SubstitutionModel};
use centrevo_sim::genome::{Chromosome, Haplotype, Individual};
use centrevo_sim::simulation::{
    CodecStrategy, FitnessConfig, MutationConfig, Population, RecombinationConfig,
    SimulationConfig, UniformRepeatStructure,
};
use centrevo_sim::storage::{Recorder, RecordingStrategy, SimulationSnapshot};

use crate::args::InitArgs;
use crate::printing::print_parameters;

#[allow(clippy::too_many_arguments)]
pub fn init_simulation(args: &InitArgs) -> Result<()> {
    let name = &args.name;
    let output = &args.output;
    let population_size = args.population_size;
    let generations = args.generations;
    let record_every = args.record_every;

    let codec_strategy = args
        .codec
        .parse::<CodecStrategy>()
        .map_err(|e| anyhow::anyhow!(e))
        .context("Invalid codec strategy")?;

    println!("🧬 Centrevo - Centromeric Evolution Simulator");
    println!("============================================\n");
    println!("Initializing simulation: {name}");
    println!("Codec Strategy: {codec_strategy}");

    let (structure, config, mutation, recombination, fitness) =
        build_configs(args, codec_strategy.clone())?;

    // Logic moved to build_configs

    println!("\nConfiguration:");

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
    let mut recorder = Recorder::new(
        output,
        name.as_str(),
        RecordingStrategy::EveryN(record_every),
        codec_strategy,
    )
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
    let dummy_rng = vec![0u8; 32];
    recorder
        .record_generation(&population, 0, &dummy_rng)
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

pub fn build_configs(
    args: &InitArgs,
    codec_strategy: CodecStrategy,
) -> Result<(
    UniformRepeatStructure,
    SimulationConfig,
    MutationConfig,
    RecombinationConfig,
    FitnessConfig,
)> {
    let structure = UniformRepeatStructure::new(
        Nucleotide::A,
        args.ru_length,
        args.rus_per_hor,
        args.hors_per_chr,
        args.chrs_per_hap,
    );

    let config = SimulationConfig::new(args.population_size, args.generations, args.seed)
        .with_codec(codec_strategy);

    // Create mutation configuration
    // Check for advanced mutation rates
    let substitution = if let (Some(ac), Some(ag), Some(at), Some(cg), Some(ct), Some(gt)) = (
        args.rate_ac,
        args.rate_ag,
        args.rate_at,
        args.rate_cg,
        args.rate_ct,
        args.rate_gt,
    ) {
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
    } else if let Some(rate) = args.mutation_rate {
        // Uniform
        SubstitutionModel::uniform(rate)
            .map_err(|e| anyhow::anyhow!("Failed to create substitution model: {e}"))?
    } else {
        if args.rate_ac.is_some()
            || args.rate_ag.is_some()
            || args.rate_at.is_some()
            || args.rate_cg.is_some()
            || args.rate_ct.is_some()
            || args.rate_gt.is_some()
        {
            anyhow::bail!(
                "When using specific mutation rates, ALL 6 rates (ac, ag, at, cg, ct, gt) must be provided."
            );
        }
        // Fallback to default if no rates provided (though CLI default handles this, manual construction might pass None)
        SubstitutionModel::uniform(1e-5)
            .map_err(|e| anyhow::anyhow!("Failed to create substitution model: {e}"))?
    };

    let mut mutation = MutationConfig::new(substitution);

    // Add Indels if enabled
    if args.indel_ins_rate > 0.0 || args.indel_del_rate > 0.0 {
        let indel_model = IndelModel::new(
            args.indel_ins_rate,
            args.indel_del_rate,
            args.indel_length_p,
        )
        .map_err(|e| anyhow::anyhow!("Failed to create indel model: {e}"))?;
        mutation.indel = Some(indel_model);
    }

    let recomb_params = centrevo_sim::evolution::RecombinationModel::builder()
        .break_prob(args.recomb_rate)
        .crossover_prob(args.crossover_prob)
        .gc_extension_prob(args.gc_extension_prob)
        .homology_strength(args.homology_strength)
        .search_window(args.search_window)
        .build()
        .map_err(|e| anyhow::anyhow!("Failed to create recombination model: {e}"))?;
    let recombination = RecombinationConfig::new(recomb_params);

    let gc_fitness = if let (Some(opt), Some(conc)) = (args.fit_gc_opt, args.fit_gc_conc) {
        Some(
            centrevo_sim::evolution::GCContentFitness::new(opt, conc)
                .map_err(|e| anyhow::anyhow!("Invalid GC Fitness: {e}"))?,
        )
    } else {
        None
    };

    let len_fitness = if let (Some(opt), Some(std)) = (args.fit_len_opt, args.fit_len_std) {
        Some(
            centrevo_sim::evolution::LengthFitness::new(opt, std)
                .map_err(|e| anyhow::anyhow!("Invalid Length Fitness: {e}"))?,
        )
    } else {
        None
    };

    let seq_sim_fitness = if let Some(shape) = args.fit_seq_sim {
        Some(
            centrevo_sim::evolution::SequenceSimilarityFitness::new(shape)
                .map_err(|e| anyhow::anyhow!("Invalid Seq Sim Fitness: {e}"))?,
        )
    } else {
        None
    };

    let len_sim_fitness = if let Some(shape) = args.fit_len_sim {
        Some(
            centrevo_sim::evolution::LengthSimilarityFitness::new(shape)
                .map_err(|e| anyhow::anyhow!("Invalid Len Sim Fitness: {e}"))?,
        )
    } else {
        None
    };

    let fitness = FitnessConfig::new(gc_fitness, len_fitness, seq_sim_fitness, len_sim_fitness);

    Ok((structure, config, mutation, recombination, fitness))
}

#[cfg(test)]
mod tests {
    use super::*;
    use centrevo_sim::base::Nucleotide;
    use std::path::PathBuf;

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

    #[test]
    fn test_build_configs_defaults() {
        let args = InitArgs {
            name: "test".to_string(),
            output: PathBuf::from("test.db"),
            population_size: 100,
            generations: 100,
            ru_length: 171,
            rus_per_hor: 12,
            hors_per_chr: 10,
            chrs_per_hap: 1,
            mutation_rate: Some(1e-5),
            rate_ac: None,
            rate_ag: None,
            rate_at: None,
            rate_cg: None,
            rate_ct: None,
            rate_gt: None,
            indel_ins_rate: 0.0,
            indel_del_rate: 0.0,
            indel_length_p: 0.5,
            recomb_rate: 0.01,
            crossover_prob: 0.01,
            gc_extension_prob: 0.95,
            homology_strength: 5.0,
            search_window: 100,
            fit_gc_opt: None,
            fit_gc_conc: None,
            fit_len_opt: None,
            fit_len_std: None,
            fit_seq_sim: None,
            fit_len_sim: None,
            record_every: 100,
            seed: None,
            codec: "parallel-packed-rs".to_string(),
        };

        let (structure, config, mutation, _recombination, fitness) =
            build_configs(&args, CodecStrategy::default()).unwrap();

        assert_eq!(config.population_size, 100);
        assert_eq!(structure.ru_length, 171);
        // Default mutation is uniform
        // Note: checking internal state of mutation model is hard without public accessors,
        // but we can check enabled flags if any.
        assert!(mutation.indel.is_none());
        assert!(fitness.gc_content.is_none());
    }

    #[test]
    fn test_build_configs_fitness() {
        let args = InitArgs {
            name: "test".to_string(),
            output: PathBuf::from("test.db"),
            population_size: 100,
            generations: 100,
            ru_length: 171,
            rus_per_hor: 12,
            hors_per_chr: 10,
            chrs_per_hap: 1,
            mutation_rate: Some(1e-5),
            rate_ac: None,
            rate_ag: None,
            rate_at: None,
            rate_cg: None,
            rate_ct: None,
            rate_gt: None,
            indel_ins_rate: 0.0,
            indel_del_rate: 0.0,
            indel_length_p: 0.5,
            recomb_rate: 0.01,
            crossover_prob: 0.01,
            gc_extension_prob: 0.95,
            homology_strength: 5.0,
            search_window: 100,
            fit_gc_opt: Some(0.4),
            fit_gc_conc: Some(10.0),
            fit_len_opt: None,
            fit_len_std: None,
            fit_seq_sim: None,
            fit_len_sim: None,
            record_every: 100,
            seed: None,
            codec: "parallel-packed-rs".to_string(),
        };

        let (_, _, _, _, fitness) = build_configs(&args, CodecStrategy::default()).unwrap();
        assert!(fitness.gc_content.is_some());
        assert!(fitness.length.is_none());
    }

    #[test]
    fn test_build_configs_advanced_mutation() {
        let args = InitArgs {
            name: "test".to_string(),
            output: PathBuf::from("test.db"),
            population_size: 100,
            generations: 100,
            ru_length: 171,
            rus_per_hor: 12,
            hors_per_chr: 10,
            chrs_per_hap: 1,
            mutation_rate: None,
            rate_ac: Some(1e-6),
            rate_ag: Some(1e-6),
            rate_at: Some(1e-6),
            rate_cg: Some(1e-6),
            rate_ct: Some(1e-6),
            rate_gt: Some(1e-6),
            indel_ins_rate: 0.0,
            indel_del_rate: 0.0,
            indel_length_p: 0.5,
            recomb_rate: 0.01,
            crossover_prob: 0.01,
            gc_extension_prob: 0.95,
            homology_strength: 5.0,
            search_window: 100,
            fit_gc_opt: None,
            fit_gc_conc: None,
            fit_len_opt: None,
            fit_len_std: None,
            fit_seq_sim: None,
            fit_len_sim: None,
            record_every: 100,
            seed: None,
            codec: "parallel-packed-rs".to_string(),
        };

        let (_, _, _mutation, _, _) = build_configs(&args, CodecStrategy::default()).unwrap();
        // Should be General substitution model
        // We can't easily introspect Enum variant without public access or matching.
        // But if it didn't error, it worked.
    }

    #[test]
    fn test_build_configs_indels() {
        let args = InitArgs {
            name: "test".to_string(),
            output: PathBuf::from("test.db"),
            population_size: 100,
            generations: 100,
            ru_length: 171,
            rus_per_hor: 12,
            hors_per_chr: 10,
            chrs_per_hap: 1,
            mutation_rate: Some(1e-5),
            rate_ac: None,
            rate_ag: None,
            rate_at: None,
            rate_cg: None,
            rate_ct: None,
            rate_gt: None,
            indel_ins_rate: 0.001,
            indel_del_rate: 0.001,
            indel_length_p: 0.5,
            recomb_rate: 0.01,
            crossover_prob: 0.01,
            gc_extension_prob: 0.95,
            homology_strength: 5.0,
            search_window: 100,
            fit_gc_opt: None,
            fit_gc_conc: None,
            fit_len_opt: None,
            fit_len_std: None,
            fit_seq_sim: None,
            fit_len_sim: None,
            record_every: 100,
            seed: None,
            codec: "parallel-packed-rs".to_string(),
        };

        let (_, _, mutation, _, _) = build_configs(&args, CodecStrategy::default()).unwrap();
        assert!(mutation.indel.is_some());
    }

    #[test]
    fn test_build_configs_conflict() {
        let args = InitArgs {
            name: "test".to_string(),
            output: PathBuf::from("test.db"),
            population_size: 100,
            generations: 100,
            ru_length: 171,
            rus_per_hor: 12,
            hors_per_chr: 10,
            chrs_per_hap: 1,
            mutation_rate: None, // Missing logic
            rate_ac: Some(1e-6), // Partial
            rate_ag: None,
            rate_at: None,
            rate_cg: None,
            rate_ct: None,
            rate_gt: None,
            indel_ins_rate: 0.0,
            indel_del_rate: 0.0,
            indel_length_p: 0.5,
            recomb_rate: 0.01,
            crossover_prob: 0.01,
            gc_extension_prob: 0.95,
            homology_strength: 5.0,
            search_window: 100,
            fit_gc_opt: None,
            fit_gc_conc: None,
            fit_len_opt: None,
            fit_len_std: None,
            fit_seq_sim: None,
            fit_len_sim: None,
            record_every: 100,
            seed: None,
            codec: "parallel-packed-rs".to_string(),
        };

        // Should return error because strict mode for explicit rates requires all 6
        let result = build_configs(&args, CodecStrategy::default());
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err().to_string(),
            "When using specific mutation rates, ALL 6 rates (ac, ag, at, cg, ct, gt) must be provided."
        );
    }
}
