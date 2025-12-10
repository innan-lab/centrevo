use anyhow::{Context, Result};
use std::io::{self, Write};
use std::path::PathBuf;

use crate::commands::init::init_simulation;
use crate::commands::run::run_simulation;

pub fn setup_wizard(defaults: bool) -> Result<()> {
    println!("\n🧙 Centrevo Setup Wizard");
    println!("========================\n");
    println!("This wizard will guide you through setting up and running a simulation.\n");

    // Simulation configuration
    let name = if defaults {
        "simulation".to_string()
    } else {
        prompt_string("Simulation name", Some("simulation"))?
    };

    let database = if defaults {
        PathBuf::from("simulation.db")
    } else {
        PathBuf::from(prompt_string("Database file", Some("simulation.db"))?)
    };

    println!("\n📐 Population Parameters");
    println!("-----------------------");

    let population_size = if defaults {
        100
    } else {
        prompt_usize("Population size", Some(100))?
    };

    let generations = if defaults {
        1000
    } else {
        prompt_usize("Number of generations", Some(1000))?
    };

    println!("\n🧬 Genome Structure");
    println!("------------------");

    let ru_length = if defaults {
        171
    } else {
        prompt_usize("Repeat unit (RU) length (bp)", Some(171))?
    };

    let rus_per_hor = if defaults {
        12
    } else {
        prompt_usize("RUs per Higher-Order Repeat (HOR)", Some(12))?
    };

    let hors_per_chr = if defaults {
        100
    } else {
        prompt_usize("HORs per chromosome", Some(100))?
    };

    println!("\n🧪 Evolution Parameters");
    println!("----------------------");

    let mutation_rate = if defaults {
        0.001
    } else {
        prompt_f64("Mutation rate", Some(0.001))?
    };

    let recomb_rate = if defaults {
        0.01
    } else {
        prompt_f64("Recombination break probability", Some(0.01))?
    };

    let crossover_prob = if defaults {
        0.7
    } else {
        prompt_f64("Crossover probability (given break)", Some(0.7))?
    };

    println!("\n💾 Recording Options");
    println!("-------------------");

    let record_every = if defaults {
        100
    } else {
        prompt_usize("Record every N generations", Some(100))?
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
    println!("  Mutation rate: {mutation_rate}");
    println!("  Recombination rate: {recomb_rate}");
    println!("  Crossover probability: {crossover_prob}");
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
        mutation_rate: Some(mutation_rate),
        rate_ac: None,
        rate_ag: None,
        rate_at: None,
        rate_cg: None,
        rate_ct: None,
        rate_gt: None,
        indel_ins_rate: 0.0,
        indel_del_rate: 0.0,
        indel_length_p: 0.5,
        recomb_rate,
        crossover_prob,
        gc_extension_prob: 0.95,
        homology_strength: 5.0,
        search_window: 100,
        fit_gc_opt: None,
        fit_gc_conc: None,
        fit_len_opt: None,
        fit_len_std: None,
        fit_seq_sim: None,
        fit_len_sim: None,
        record_every,
        seed,
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
