use centrevo_sim::simulation::Configuration;
use centrevo_sim::simulation::InitializationConfig;
use centrevo_sim::simulation::Simulation;

pub fn print_simulation_parameters(sim: &Simulation) {
    print_parameters(sim.configuration());
}

pub fn print_parameters(config: &Configuration) {
    let sim_config = &config.execution;
    let mutation = &config.evolution.mutation;
    let recombination = &config.evolution.recombination;
    let fitness = &config.evolution.fitness;
    let structure = match &config.initialization {
        InitializationConfig::Generate { structure, .. } => Some(structure),
        InitializationConfig::Load { .. } => None,
    };
    println!("\n📋 Simulation Configuration");
    println!(
        "  • Population Size: {} [-n, --population-size]",
        sim_config.population_size
    );
    println!(
        "  • Generations: {} [-g, --generations]",
        sim_config.total_generations
    );
    if let Some(seed) = sim_config.seed {
        println!("  • Random Seed: {seed} [--seed]");
    } else {
        println!("  • Random Seed: Random [--seed]");
    }

    if let Some(s) = structure {
        println!("\n🧬 Genome Structure");
        println!("  • Initial Base: {:?}", s.init_base);
        println!("  • RU Length: {} bp [--ru-length]", s.ru_length);
        println!(
            "  • Structure: {} RUs/HOR [--rus-per-hor] × {} HORs/Chr [--hors-per-chr]",
            s.rus_per_hor, s.hors_per_chr
        );
        println!("  • Chromosome Length: {} bp", s.chr_length());
        println!("  • Ploidy: Diploid (2 haplotypes)");
    }

    println!("\n⚡ Mutation Parameters");
    let rates = mutation.substitution.rates();
    let total_sub = rates.iter().sum::<f64>() / 2.0; // Approximation for summary
    println!("  • Substitution Rate: {total_sub:.2e} (approx per base) [--mutation-rate]",);
    println!(
        "    - A↔C: {:.2e}, A↔G: {:.2e}, A↔T: {:.2e}",
        rates[0], rates[1], rates[2]
    );
    println!(
        "    - C↔G: {:.2e}, C↔T: {:.2e}, G↔T: {:.2e}",
        rates[3], rates[4], rates[5]
    );

    if let Some(indel) = &mutation.indel {
        println!("  • Indel Model: Enabled");
        println!("    - Insertion Rate: {:.2e}", indel.insertion_rate());
        println!("    - Deletion Rate: {:.2e}", indel.deletion_rate());
        println!(
            "    - Length Param (p): {:.2e} (mean len ~{:.1} bp)",
            indel.length_p(),
            1.0 / indel.length_p()
        );
    } else {
        println!("  • Indel Model: Disabled");
    }

    println!("\n🔀 Recombination Parameters");
    println!(
        "  • Break Probability: {:.2e} /bp/gen [--recomb-rate]",
        recombination.params.break_prob()
    );
    println!(
        "  • Crossover Probability: {:.2} [--crossover-prob]",
        recombination.params.crossover_prob()
    );
    println!(
        "  • GC Extension Prob: {:.2} [--gc-extension-prob]",
        recombination.params.gc_extension_prob()
    );
    println!(
        "  • Homology Strength: {:.2} [--homology-strength]",
        recombination.params.homology_strength()
    );
    println!(
        "  • Search Window: {} RUs [--search-window]",
        recombination.params.search_window()
    );

    println!("\n🎯 Fitness & Selection");
    if fitness.is_neutral() {
        println!("  • Regime: Neutral Evolution (No Selection)");
    } else {
        println!("  • Regime: Selection Enabled");
        if let Some(gc) = &fitness.gc_content {
            println!(
                "    - GC Content: Optimum={:.2}, Concentration={:.2}",
                gc.optimum, gc.concentration
            );
        }
        if let Some(len) = &fitness.length {
            println!(
                "    - Length: Optimum={}, StdDev={:.2}",
                len.optimum, len.std_dev
            );
        }
        if let Some(sim) = &fitness.seq_similarity {
            println!("    - Seq Similarity: Shape={:.2}", sim.shape);
        }
        if let Some(len_sim) = &fitness.length_similarity {
            println!("    - Length Similarity: Shape={:.2}", len_sim.shape);
        }
    }
    println!();
}
