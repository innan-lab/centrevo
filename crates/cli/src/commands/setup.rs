use anyhow::{Context, Result};
use std::io::{self, Write};
use std::path::PathBuf;

use crate::commands::init::init_simulation;
use crate::commands::run::run_simulation;

pub fn setup_wizard(defaults: bool) -> Result<()> {
    println!("\n🧙 Centrevo Setup Wizard");
    println!("========================\n");
    println!("This wizard will guide you through setting up and running a simulation.\n");

    use crate::defaults;

    // Simulation configuration
    let name = if defaults {
        defaults::SIMULATION_NAME.to_string()
    } else {
        prompt_string("Simulation name", Some(defaults::SIMULATION_NAME))?
    };

    let database = if defaults {
        PathBuf::from(defaults::OUTPUT_DB)
    } else {
        PathBuf::from(prompt_string("Database file", Some(defaults::OUTPUT_DB))?)
    };

    println!("\n📐 Population Parameters");
    println!("-----------------------");

    let population_size = if defaults {
        defaults::POPULATION_SIZE
    } else {
        prompt_usize("Population size", Some(defaults::POPULATION_SIZE))?
    };

    let generations = if defaults {
        defaults::GENERATIONS
    } else {
        prompt_usize("Number of generations", Some(defaults::GENERATIONS))?
    };

    println!("\n🧬 Genome Structure");
    println!("------------------");

    let ru_length = if defaults {
        defaults::RU_LENGTH
    } else {
        prompt_usize("Repeat unit (RU) length (bp)", Some(defaults::RU_LENGTH))?
    };

    let rus_per_hor = if defaults {
        defaults::RUS_PER_HOR
    } else {
        prompt_usize(
            "RUs per Higher-Order Repeat (HOR)",
            Some(defaults::RUS_PER_HOR),
        )?
    };

    let hors_per_chr = if defaults {
        defaults::HORS_PER_CHR
    } else {
        prompt_usize("HORs per chromosome", Some(defaults::HORS_PER_CHR))?
    };

    println!("\n🧪 Evolution Parameters");
    println!("----------------------");

    let (mutation_rate, rate_ac, rate_ag, rate_at, rate_cg, rate_ct, rate_gt) = if defaults {
        (
            Some(defaults::MUTATION_RATE),
            None,
            None,
            None,
            None,
            None,
            None,
        )
    } else if !prompt_confirm(
        "Configure advanced mutation model (e.g. specific rates)?",
        false,
    )? {
        let rate = prompt_f64("Mutation rate (uniform)", Some(defaults::MUTATION_RATE))?;
        (Some(rate), None, None, None, None, None, None)
    } else {
        println!("  Enter substitution rates (probability per base per generation):");
        let ac = prompt_f64("  Rate A<->C", Some(defaults::SPECIFIC_RATE_DEFAULT))?;
        let ag = prompt_f64("  Rate A<->G", Some(defaults::SPECIFIC_RATE_DEFAULT))?;
        let at = prompt_f64("  Rate A<->T", Some(defaults::SPECIFIC_RATE_DEFAULT))?;
        let cg = prompt_f64("  Rate C<->G", Some(defaults::SPECIFIC_RATE_DEFAULT))?;
        let ct = prompt_f64("  Rate C<->T", Some(defaults::SPECIFIC_RATE_DEFAULT))?;
        let gt = prompt_f64("  Rate G<->T", Some(defaults::SPECIFIC_RATE_DEFAULT))?;
        (
            None,
            Some(ac),
            Some(ag),
            Some(at),
            Some(cg),
            Some(ct),
            Some(gt),
        )
    };

    let (indel_ins_rate, indel_del_rate, indel_length_p) = if defaults {
        (
            defaults::INDEL_INS_RATE,
            defaults::INDEL_DEL_RATE,
            defaults::INDEL_LENGTH_P,
        )
    } else if prompt_confirm("Enable indel mutations?", false)? {
        let ins = prompt_f64("  Insertion rate", Some(1e-5))?;
        let del = prompt_f64("  Deletion rate", Some(1e-5))?;
        let p = prompt_f64(
            "  Indel length parameter p (higher = shorter indels)",
            Some(defaults::INDEL_LENGTH_P),
        )?;
        (ins, del, p)
    } else {
        (
            defaults::INDEL_INS_RATE,
            defaults::INDEL_DEL_RATE,
            defaults::INDEL_LENGTH_P,
        )
    };

    let recomb_rate = if defaults {
        defaults::RECOMB_RATE
    } else {
        prompt_f64(
            "Recombination break probability (per base)",
            Some(defaults::RECOMB_RATE),
        )?
    };

    let crossover_prob = if defaults {
        defaults::CROSSOVER_PROB
    } else {
        prompt_f64(
            "Crossover probability (vs. gene conversion)",
            Some(defaults::CROSSOVER_PROB),
        )?
    };

    let (gc_extension_prob, homology_strength, search_window) = if defaults {
        (
            defaults::GC_EXTENSION_PROB,
            defaults::HOMOLOGY_STRENGTH,
            defaults::SEARCH_WINDOW,
        )
    } else if prompt_confirm("Configure advanced recombination parameters?", false)? {
        let gc = prompt_f64(
            "  Gene conversion extension prob (higher = longer tracts)",
            Some(defaults::GC_EXTENSION_PROB),
        )?;
        let hom = prompt_f64(
            "  Homology strength (higher = stricter matching)",
            Some(defaults::HOMOLOGY_STRENGTH),
        )?;
        let win = prompt_usize(
            "  Search window (in Repeat Units)",
            Some(defaults::SEARCH_WINDOW),
        )?;
        (gc, hom, win)
    } else {
        (
            defaults::GC_EXTENSION_PROB,
            defaults::HOMOLOGY_STRENGTH,
            defaults::SEARCH_WINDOW,
        )
    };

    println!("\n🎯 Fitness & Selection");
    println!("---------------------");

    let (fit_gc_opt, fit_gc_conc) = if defaults {
        (None, None)
    } else if prompt_confirm("Enable GC content selection?", false)? {
        let opt = prompt_f64("  Optimal GC content (0.0-1.0)", Some(defaults::FIT_GC_OPT))?;
        let conc = prompt_f64("  Concentration parameter", Some(defaults::FIT_GC_CONC))?;
        (Some(opt), Some(conc))
    } else {
        (None, None)
    };

    let (fit_len_opt, fit_len_std) = if defaults {
        (None, None)
    } else if prompt_confirm("Enable Chromosome Length selection?", false)? {
        let total_bp_est = ru_length * rus_per_hor * hors_per_chr;
        let opt = prompt_usize("  Optimal length (bp)", Some(total_bp_est))?;
        let std = prompt_f64("  Length std dev (bp)", Some(defaults::FIT_LEN_STD))?;
        (Some(opt), Some(std))
    } else {
        (None, None)
    };

    let fit_seq_sim = if defaults {
        None
    } else if prompt_confirm("Enable Sequence Similarity selection?", false)? {
        Some(prompt_f64(
            "  Similarity shape parameter",
            Some(defaults::FIT_SEQ_SIM_SHAPE),
        )?)
    } else {
        None
    };

    let fit_len_sim = if defaults {
        None
    } else if prompt_confirm("Enable Length Similarity selection?", false)? {
        Some(prompt_f64(
            "  Similarity shape parameter",
            Some(defaults::FIT_LEN_SIM_SHAPE),
        )?)
    } else {
        None
    };

    println!("\n💾 Recording Options");
    println!("-------------------");

    let record_every = if defaults {
        defaults::RECORD_EVERY
    } else {
        prompt_usize("Record every N generations", Some(defaults::RECORD_EVERY))?
    };

    let seed = if defaults {
        None
    } else {
        match prompt_string("Random seed (press Enter for random)", None)? {
            s if s.is_empty() => None,
            s => Some(s.parse::<u64>().context("Invalid seed value")?),
        }
    };

    // Summary
    let total_bp = ru_length * rus_per_hor * hors_per_chr;

    println!("\n📋 Configuration Summary");
    println!("========================");
    println!("Simulation name: {name}");
    println!("Database: {}", database.display());
    println!();
    println!("Population:");
    println!("  Size: {population_size} individuals");
    println!("  Generations: {generations}");
    println!();
    println!("Genome:");
    println!("  RU length: {ru_length} bp");
    println!("  RUs per HOR: {rus_per_hor}");
    println!("  HORs per chromosome: {hors_per_chr}");
    println!(
        "  Total chromosome length: {} bp ({:.1} kb)",
        total_bp,
        total_bp as f64 / 1000.0
    );
    println!();
    println!("Evolution:");
    println!("Evolution:");
    if let Some(rate) = mutation_rate {
        println!("  Mutation rate: {rate} (uniform)");
    } else {
        println!("  Mutation rate: Custom (Advanced)");
    }
    if indel_ins_rate > 0.0 {
        println!("  Indels: Enabled (Ins={indel_ins_rate}, Del={indel_del_rate})");
    } else {
        println!("  Indels: Disabled");
    }
    println!("  Recombination rate: {recomb_rate}");
    println!("  Crossover probability: {crossover_prob}");

    let mut fitness_modes = Vec::new();
    if fit_gc_opt.is_some() {
        fitness_modes.push("GC");
    }
    if fit_len_opt.is_some() {
        fitness_modes.push("Length");
    }
    if fit_seq_sim.is_some() {
        fitness_modes.push("SeqSim");
    }
    if fit_len_sim.is_some() {
        fitness_modes.push("LenSim");
    }

    if fitness_modes.is_empty() {
        println!("  Fitness: Neutral");
    } else {
        println!("  Fitness: Selection ({})", fitness_modes.join(", "));
    }
    println!();
    println!("Recording:");
    println!("  Every {record_every} generations");
    println!("  Total snapshots: ~{}", (generations / record_every) + 1);
    if let Some(s) = seed {
        println!("  Seed: {s}");
    } else {
        println!("  Seed: random");
    }

    // Confirm
    if !defaults {
        println!();
        if !prompt_confirm("Proceed with simulation?", true)? {
            println!("\n❌ Setup cancelled.");
            return Ok(());
        }
    }

    use crate::args::InitArgs;

    // Initialize simulation
    println!("\n🚀 Starting simulation setup...\n");

    let args = InitArgs {
        name: name.clone(),
        output: database.clone(),
        population_size,
        generations,
        ru_length,
        rus_per_hor,
        hors_per_chr,
        chrs_per_hap: 1, // default
        mutation_rate,
        rate_ac,
        rate_ag,
        rate_at,
        rate_cg,
        rate_ct,
        rate_gt,
        indel_ins_rate,
        indel_del_rate,
        indel_length_p,
        recomb_rate,
        crossover_prob,
        gc_extension_prob,
        homology_strength,
        search_window,
        fit_gc_opt,
        fit_gc_conc,
        fit_len_opt,
        fit_len_std,
        fit_seq_sim,
        fit_len_sim,
        record_every,
        seed,
        codec: "parallel-packed-rs".to_string(),
    };

    init_simulation(&args)?;

    // Ask if user wants to run now
    let should_run = if defaults {
        true
    } else {
        println!();
        prompt_confirm("Run simulation now?", true)?
    };

    if should_run {
        println!();
        run_simulation(
            &database, &name, false, // not resume
            None,  // no seed override (use configured)
            None,  // no record override (use configured)
            true,  // show progress
        )?;
    } else {
        println!("\n💡 To run later, use: centrevo run -N {name}");
    }

    Ok(())
}

fn prompt_string(prompt: &str, default: Option<&str>) -> Result<String> {
    print!("{prompt}");
    if let Some(def) = default {
        print!(" [{def}]");
    }
    print!(": ");
    io::stdout().flush()?;

    let mut input = String::new();
    io::stdin().read_line(&mut input)?;
    let input = input.trim();

    if input.is_empty() {
        if let Some(def) = default {
            Ok(def.to_string())
        } else {
            anyhow::bail!("Input required");
        }
    } else {
        Ok(input.to_string())
    }
}

fn prompt_usize(prompt: &str, default: Option<usize>) -> Result<usize> {
    print!("{prompt}");
    if let Some(def) = default {
        print!(" [{def}]");
    }
    print!(": ");
    io::stdout().flush()?;

    let mut input = String::new();
    io::stdin().read_line(&mut input)?;
    let input = input.trim();

    if input.is_empty() {
        default.ok_or_else(|| anyhow::anyhow!("Input required"))
    } else {
        input
            .parse::<usize>()
            .with_context(|| format!("Invalid number: {input}"))
    }
}

fn prompt_f64(prompt: &str, default: Option<f64>) -> Result<f64> {
    print!("{prompt}");
    if let Some(def) = default {
        print!(" [{def}]");
    }
    print!(": ");
    io::stdout().flush()?;

    let mut input = String::new();
    io::stdin().read_line(&mut input)?;
    let input = input.trim();

    if input.is_empty() {
        default.ok_or_else(|| anyhow::anyhow!("Input required"))
    } else {
        input
            .parse::<f64>()
            .with_context(|| format!("Invalid number: {input}"))
    }
}

fn prompt_confirm(prompt: &str, default: bool) -> Result<bool> {
    let default_str = if default { "Y/n" } else { "y/N" };
    print!("{prompt} [{default_str}]: ");
    io::stdout().flush()?;

    let mut input = String::new();
    io::stdin().read_line(&mut input)?;
    let input = input.trim().to_lowercase();

    if input.is_empty() {
        Ok(default)
    } else if input == "y" || input == "yes" {
        Ok(true)
    } else if input == "n" || input == "no" {
        Ok(false)
    } else {
        anyhow::bail!("Please answer 'y' or 'n'")
    }
}
