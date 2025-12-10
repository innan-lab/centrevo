mod commands;
mod printing;
mod utils;

use anyhow::Result;
use clap::{Parser, Subcommand};
use std::path::PathBuf;

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
    Init {
        /// Simulation name
        #[arg(short = 'N', long, default_value = "simulation")]
        name: String,

        /// Output database path
        #[arg(short, long, default_value = "simulation.db")]
        output: PathBuf,

        /// Population size
        #[arg(short = 'n', long, default_value = "100")]
        population_size: usize,

        /// Number of generations
        #[arg(short = 'g', long, default_value = "1000")]
        generations: usize,

        /// Repeat unit length
        #[arg(long, default_value = "171")]
        ru_length: usize,

        /// Repeat units per HOR
        #[arg(long, default_value = "12")]
        rus_per_hor: usize,

        /// HORs per chromosome
        #[arg(long, default_value = "100")]
        hors_per_chr: usize,

        /// Mutation rate
        ///
        /// Defaults to 1e-5. Although typical per-base mutation rate is ~1e-8, 1e-5 is useful for simulation timescales.
        #[arg(long, default_value = "1e-5")]
        mutation_rate: f64,

        /// Recombination break probability
        ///
        /// Defaults to 1e-6 or approx. 1 DSB per Mb per generation.
        #[arg(long, default_value = "1e-6")]
        recomb_rate: f64,

        /// Crossover probability (given break)
        ///
        /// Defaults to 0.01. Crossovers are heavily suppressed in centromeres; most events are gene conversions.
        #[arg(long, default_value = "0.01")]
        crossover_prob: f64,

        /// Gene conversion extension probability
        ///
        /// Defaults to 0.95. Results in average tract length ~20bp (1/(1-0.95)).
        #[arg(long, default_value = "0.95")]
        gc_extension_prob: f64,

        /// Homology strength (0.0 = random, >0.0 = preference for similarity)
        ///
        /// Defaults to 5.0. (Strong preference for homologous sequences to drive homogenization)
        #[arg(long, default_value = "5.0")]
        homology_strength: f64,

        /// Search window for homology (in RUs)
        ///
        /// Defaults to 100. (Allows interaction with neighboring ~10-20kb of sequence)
        #[arg(long, default_value = "100")]
        search_window: usize,

        /// Default recording interval (record every N generations)
        #[arg(long, default_value = "100")]
        record_every: usize,

        /// Random seed
        #[arg(long)]
        seed: Option<u64>,
    },

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

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Init {
            name,
            output,
            population_size,
            generations,
            ru_length,
            rus_per_hor,
            hors_per_chr,
            mutation_rate,
            recomb_rate,
            crossover_prob,
            gc_extension_prob,
            homology_strength,
            search_window,
            record_every,
            seed,
        } => {
            init::init_simulation(
                &name,
                &output,
                population_size,
                generations,
                ru_length,
                rus_per_hor,
                hors_per_chr,
                mutation_rate,
                recomb_rate,
                crossover_prob,
                gc_extension_prob,
                homology_strength,
                search_window,
                record_every,
                seed,
            )?;
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
