mod args;
mod commands;
mod printing;
mod utils;

use anyhow::Result;
use clap::{Parser, Subcommand};
use std::path::PathBuf;

use args::InitArgs;
use commands::{analyze, export, init, list, run, setup, validate};

/// Centrevo - Centromeric evolution simulator
#[derive(Parser, Debug)]
#[command(name = "centrevo")]
#[command(author, version, about = "Centromeric evolution simulator", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Initialize a new simulation
    Init(Box<InitArgs>),

    /// Run a simulation
    Run {
        /// Database path
        #[arg(short, long, default_value = "simulation.db")]
        database: PathBuf,

        /// Simulation name
        #[arg(short = 'N', long)]
        name: String,

        /// Resume from checkpoint (ignores other parameters)
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

    /// List simulations in database
    List {
        /// Database path
        #[arg(short, long, default_value = "simulation.db")]
        database: PathBuf,
    },

    /// Show simulation info
    Info {
        /// Database path
        #[arg(short, long, default_value = "simulation.db")]
        database: PathBuf,

        /// Simulation name
        #[arg(short = 'N', long)]
        name: String,
    },

    /// List recorded generations
    Generations {
        /// Database path
        #[arg(short, long, default_value = "simulation.db")]
        database: PathBuf,

        /// Simulation name
        #[arg(short = 'N', long)]
        name: String,
    },

    /// Export simulation data
    Export {
        /// Database path
        #[arg(short, long, default_value = "simulation.db")]
        database: PathBuf,

        /// Simulation name
        #[arg(short = 'N', long)]
        name: String,

        /// Generation to export
        #[arg(short, long)]
        generation: usize,

        /// Output format (csv, json, fasta)
        #[arg(short, long, default_value = "csv")]
        format: String,

        /// Output file (stdout if not specified)
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// What to export (sequences, metadata, fitness)
        #[arg(long, default_value = "sequences")]
        data_type: String,
    },

    /// Analyze simulation data
    Analyze {
        /// Database path
        #[arg(short, long, default_value = "simulation.db")]
        database: PathBuf,

        /// Simulation name
        #[arg(short = 'N', long)]
        name: String,

        /// Generation to analyze
        #[arg(short, long)]
        generation: usize,

        /// Chromosome index
        #[arg(long, default_value = "0")]
        chromosome: usize,

        /// Output format (pretty, json)
        #[arg(short, long, default_value = "pretty")]
        format: String,

        /// Output file (stdout if not specified)
        #[arg(short, long)]
        output: Option<PathBuf>,
    },

    /// Validate database integrity
    Validate {
        /// Database path
        #[arg(short, long, default_value = "simulation.db")]
        database: PathBuf,

        /// Simulation name (all if not specified)
        #[arg(short = 'N', long)]
        name: Option<String>,

        /// Fix issues if possible
        #[arg(long)]
        fix: bool,
    },

    /// Interactive wizard to setup and run a simulation
    Setup {
        /// Skip interactive prompts and use defaults
        #[arg(long)]
        defaults: bool,
    },
}

// InitArgs moved to args.rs

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Init(args) => {
            init::init_simulation(&args)?;
        }
        Commands::Run {
            database,
            name,
            resume,
            seed,
            record_every,
            progress,
        } => {
            run::run_simulation(&database, &name, resume, seed, record_every, progress)?;
        }
        Commands::List { database } => {
            list::list_simulations(&database)?;
        }
        Commands::Info { database, name } => {
            list::show_info(&database, &name)?;
        }
        Commands::Generations { database, name } => {
            list::show_generations(&database, &name)?;
        }
        Commands::Export {
            database,
            name,
            generation,
            format,
            output,
            data_type,
        } => {
            export::export_data(
                &database,
                &name,
                generation,
                &format,
                output.as_ref(),
                &data_type,
            )?;
        }
        Commands::Analyze {
            database,
            name,
            generation,
            chromosome,
            format,
            output,
        } => {
            analyze::analyze_data(
                &database,
                &name,
                generation,
                chromosome,
                &format,
                output.as_ref(),
            )?;
        }
        Commands::Validate {
            database,
            name,
            fix,
        } => {
            validate::validate_database(&database, name.as_deref(), fix)?;
        }
        Commands::Setup { defaults } => {
            setup::setup_wizard(defaults)?;
        }
    }

    Ok(())
}
