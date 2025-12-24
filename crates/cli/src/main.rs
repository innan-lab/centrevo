mod args;
mod commands;
pub mod defaults;
mod printing;
mod utils;

use anyhow::Result;
use clap::{Parser, Subcommand};
use std::path::PathBuf;

use args::InitArgs;
use commands::{export, init, inspect, run, setup, validate};

/// Centrevo: A Centromere Evolution Simulator
///
/// This tool simulates how centromeres (complex DNA regions) change over time
/// due to mutation, recombination, and natural selection.
#[derive(Parser, Debug)]
#[command(name = "centrevo")]
#[command(author, version, about = "Simulates the evolution of centromeres over time", long_about = None)]
struct Cli {
    /// Number of threads to use for parallel processing
    ///
    /// If not specified, defaults to the number of logical CPUs.
    #[arg(short = 't', long, global = true)]
    threads: Option<usize>,

    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Initialize a new simulation configuration.
    ///
    /// Sets up the parameters for a new experiment (population size, mutation rates, etc.)
    /// but does not run it yet.
    Init(Box<InitArgs>),

    /// Run an existing or new simulation.
    ///
    /// Executes the simulation generation by generation.
    Run {
        /// Database path (where to save data)
        #[arg(short, long, default_value = "simulation.db")]
        database: PathBuf,

        /// Resume from the last saved checkpoint
        ///
        /// Use this if a previous run was interrupted or if you want to extend it.
        #[arg(long)]
        resume: bool,

        /// Override random seed (default: use configured seed)
        #[arg(long)]
        seed: Option<u64>,

        /// Override recording interval (default: use configured interval)
        #[arg(long)]
        record_every: Option<usize>,

        /// Show progress bar
        #[arg(long, default_value = "true")]
        progress: bool,
    },

    /// Info: Show detailed configuration of a simulation.
    Info {
        /// Database path
        #[arg(short, long, default_value = "simulation.db")]
        database: PathBuf,
    },

    /// Generations: List all recorded timepoints for a simulation.
    Generations {
        /// Database path
        #[arg(short, long, default_value = "simulation.db")]
        database: PathBuf,
    },

    /// Export data to other formats (CSV, FASTA, JSON).
    ///
    /// Use this to get data out for analysis in Python, R, or other tools.
    Export {
        /// Database path
        #[arg(short, long, default_value = "simulation.db")]
        database: PathBuf,

        /// Generation(s) to export (e.g., "100", "0..1000", "10,20", "all")
        #[arg(short = 'g', long)]
        generations: Option<String>,

        /// Output format (csv, json, fasta, formatted-text)
        #[arg(short, long, default_value = "csv")]
        format: String,

        /// Output file
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// What to export (sequences, metadata, fitness)
        #[arg(long, default_value = "sequences")]
        data_type: String,

        /// Delimiter for Repeat Units (formatted-text only)
        #[arg(long, default_value = ":")]
        ru_delimiter: String,

        /// Delimiter for Higher Order Repeats (formatted-text only)
        #[arg(long, default_value = "|")]
        hor_delimiter: String,

        /// Delimiter for Chromosomes (formatted-text only)
        #[arg(long, default_value = "\n")]
        chr_delimiter: String,
    },

    // /// Analyze simulation data (statistics).
    // ///
    // /// Calculates stats like sequence length, GC content, and entropy.
    // Analyze {
    //     /// Database path
    //     #[arg(short, long, default_value = "simulation.db")]
    //     database: PathBuf,

    //     /// Generation to analyze
    //     #[arg(short, long)]
    //     generation: usize,

    //     /// Chromosome index (which chromosome to analyze)
    //     #[arg(long, default_value = "0")]
    //     chromosome: usize,

    //     /// Output format (pretty, json)
    //     #[arg(short, long, default_value = "pretty")]
    //     format: String,

    //     /// Output file (stdout if not specified)
    //     #[arg(short, long)]
    //     output: Option<PathBuf>,
    // },
    /// Validate database integrity.
    ///
    /// Checks if the database file is healthy and fixes simple issues.
    Validate {
        /// Database path
        #[arg(short, long, default_value = "simulation.db")]
        database: PathBuf,

        /// Attempt to fix found issues
        #[arg(long)]
        fix: bool,
    },

    /// Setup Wizard: Interactive guide to create a simulation.
    ///
    /// The easiest way to get started. Just follow the prompts!
    Setup {
        /// Skip interactive prompts and use defaults
        #[arg(long)]
        defaults: bool,
    },
}

// InitArgs moved to args.rs

fn main() -> Result<()> {
    let cli = Cli::parse();

    if let Some(threads) = cli.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()?;
    }

    match cli.command {
        Commands::Init(args) => {
            init::init_simulation(&args)?;
        }
        Commands::Run {
            database,
            resume,
            seed,
            record_every,
            progress,
        } => {
            run::run_simulation(&database, resume, seed, record_every, progress)?;
        }
        Commands::Info { database } => {
            inspect::show_info(&database)?;
        }
        Commands::Generations { database } => {
            inspect::show_generations(&database)?;
        }
        Commands::Export {
            database,
            generations,
            format,
            output,
            data_type,
            ru_delimiter,
            hor_delimiter,
            chr_delimiter,
        } => {
            export::export_data(
                &database,
                generations,
                &format,
                output.as_ref(),
                &data_type,
                &ru_delimiter,
                &hor_delimiter,
                &chr_delimiter,
            )?;
        }
        // Commands::Analyze {
        //     database,
        //     generation,
        //     chromosome,
        //     format,
        //     output,
        // } => {
        //     analyze::analyze_data(&database, generation, chromosome, &format, output.as_ref())?;
        // }
        Commands::Validate { database, fix } => {
            validate::validate_database(&database, fix)?;
        }
        Commands::Setup { defaults } => {
            setup::setup_wizard(defaults)?;
        }
    }

    Ok(())
}
