//! Centrevo CLI - Command-line interface for centromeric evolution simulations.

use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use centrevo::base::{Alphabet, Nucleotide};
use centrevo::genome::{Chromosome, Haplotype, Individual};
use centrevo::simulation::{
    Population, RepeatStructure, SimulationConfig,
    Simulation, MutationConfig, RecombinationConfig, FitnessConfig,
};
use centrevo::storage::{QueryBuilder, Recorder, RecordingStrategy};
use std::path::PathBuf;
use std::io::{self, Write};

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

        /// Mutation rate
        #[arg(long, default_value = "0.001")]
        mutation_rate: f64,

        /// Recombination break probability
        #[arg(long, default_value = "0.01")]
        recomb_rate: f64,

        /// Crossover probability (given break)
        #[arg(long, default_value = "0.7")]
        crossover_prob: f64,

        /// Recording interval (record every N generations)
        #[arg(long, default_value = "100")]
        record_every: usize,

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
            seed,
        } => {
            init_simulation(
                &name,
                &output,
                population_size,
                generations,
                ru_length,
                rus_per_hor,
                hors_per_chr,
                seed,
            )?;
        }
        Commands::Run {
            database,
            name,
            mutation_rate,
            recomb_rate,
            crossover_prob,
            record_every,
            progress,
        } => {
            run_simulation(
                &database,
                &name,
                mutation_rate,
                recomb_rate,
                crossover_prob,
                record_every,
                progress,
            )?;
        }
        Commands::List { database } => {
            list_simulations(&database)?;
        }
        Commands::Info { database, name } => {
            show_info(&database, &name)?;
        }
        Commands::Generations { database, name } => {
            show_generations(&database, &name)?;
        }
    }

    Ok(())
}

fn init_simulation(
    name: &str,
    output: &PathBuf,
    population_size: usize,
    generations: usize,
    ru_length: usize,
    rus_per_hor: usize,
    hors_per_chr: usize,
    seed: Option<u64>,
) -> Result<()> {
    println!("🧬 Centrevo - Centromeric Evolution Simulator");
    println!("============================================\n");
    println!("Initializing simulation: {}", name);

    // Create configurations
    let alphabet = Alphabet::dna();
    let structure = RepeatStructure::new(
        alphabet.clone(),
        Nucleotide::A,
        ru_length,
        rus_per_hor,
        hors_per_chr,
        1,
    );

    let config = SimulationConfig::new(population_size, generations, seed);

    // Create initial population
    println!("\nCreating initial population...");
    let population = create_initial_population(population_size, &structure);
    println!("✓ Created {} individuals", population.size());

    // Setup database recorder
    println!("\nSetting up database...");
    let mut recorder = Recorder::new(output, name, RecordingStrategy::EveryN(100))
        .context("Failed to create recorder")?;

    // Record metadata and initial generation
    recorder
        .record_metadata(&config)
        .context("Failed to record metadata")?;

    recorder
        .record_generation(&population, 0)
        .context("Failed to record initial generation")?;

    println!("✓ Database created: {}", output.display());
    println!("\nSimulation initialized successfully!");
    println!("  Name: {}", name);
    println!("  Population size: {}", population_size);
    println!("  Generations: {}", generations);
    println!("  Structure: {}bp RU × {} RUs/HOR × {} HORs", 
             ru_length, rus_per_hor, hors_per_chr);
    println!("\n💡 Use 'centrevo run -N {}' to run the simulation", name);

    Ok(())
}

fn run_simulation(
    database: &PathBuf,
    name: &str,
    mutation_rate: f64,
    recomb_rate: f64,
    crossover_prob: f64,
    record_every: usize,
    show_progress: bool,
) -> Result<()> {
    println!("🧬 Centrevo - Running Simulation");
    println!("============================================\n");
    
    // Open database and get simulation info
    let query = QueryBuilder::new(database).context("Failed to open database")?;
    let info = query
        .get_simulation_info(name)
        .context("Failed to get simulation info")?;
    
    println!("Simulation: {}", name);
    println!("Population size: {}", info.pop_size);
    println!("Target generations: {}", info.num_generations);
    println!("Mutation rate: {}", mutation_rate);
    println!("Recombination rate: {}", recomb_rate);
    println!();
    
    // Check if simulation already has data
    let recorded_gens = query
        .get_recorded_generations(name)
        .context("Failed to get recorded generations")?;
    
    if !recorded_gens.is_empty() && recorded_gens.iter().max().copied().unwrap_or(0) >= info.num_generations {
        println!("✓ Simulation already complete!");
        println!("  Recorded {} generations", recorded_gens.len());
        return Ok(());
    }
    
    println!("Starting simulation from generation 0...\n");
    
    // Load initial population
    let initial_individuals = query.get_generation(name, 0)
        .context("Failed to load initial population. Did you run 'centrevo init' first?")?;
    
    if initial_individuals.is_empty() {
        anyhow::bail!("No initial population found. Please run 'centrevo init' first.");
    }
    
    // Setup simulation components
    let alphabet = Alphabet::dna();
    let structure = RepeatStructure::new(
        alphabet.clone(),
        Nucleotide::A,
        171,  // TODO: These should match init params from database
        12,
        100,
        1,
    );
    
    let mutation = MutationConfig::uniform(alphabet, mutation_rate)
        .map_err(|e| anyhow::anyhow!(e))?;
    let recombination = RecombinationConfig::standard(recomb_rate, crossover_prob, 0.1)
        .map_err(|e| anyhow::anyhow!(e))?;
    let fitness = FitnessConfig::neutral();
    let config = SimulationConfig::new(info.pop_size, info.num_generations, None);
    
    // Create simulation engine
    let mut sim = Simulation::new(structure, mutation, recombination, fitness, config)
        .map_err(|e| anyhow::anyhow!(e))?;
    
    // Setup recorder
    let mut recorder = Recorder::new(
        database,
        name,
        RecordingStrategy::EveryN(record_every),
    ).context("Failed to create recorder")?;
    
    // Run simulation
    println!("Running {} generations...", info.num_generations);
    
    for generation in 1..=info.num_generations {
        sim.step().map_err(|e| anyhow::anyhow!("Generation {}: {}", generation, e))?;
        
        // Record if needed
        recorder.record_if_needed(sim.population(), generation)
            .with_context(|| format!("Failed to record generation {}", generation))?;
        
        // Show progress
        if show_progress && (generation % 10 == 0 || generation == info.num_generations) {
            let progress = generation as f64 / info.num_generations as f64 * 100.0;
            print!("\rProgress: {:.1}% ({}/{})", progress, generation, info.num_generations);
            io::stdout().flush().ok();
        }
    }
    
    if show_progress {
        println!();
    }
    
    // Finalize
    recorder.finalize_metadata()
        .context("Failed to finalize metadata")?;
    
    println!("\n✓ Simulation complete!");
    println!("  Final generation: {}", info.num_generations);
    println!("  Recorded generations: {} snapshots", 
             (info.num_generations / record_every) + 1);
    println!("\n💡 Use 'centrevo info -N {}' to view results", name);
    
    Ok(())
}

fn list_simulations(database: &PathBuf) -> Result<()> {
    let query = QueryBuilder::new(database).context("Failed to open database")?;
    let simulations = query.list_simulations().context("Failed to list simulations")?;

    if simulations.is_empty() {
        println!("No simulations found in database.");
        return Ok(());
    }

    println!("\n📊 Simulations in {}:", database.display());
    println!("{}", "=".repeat(50));

    for name in simulations {
        println!("  • {}", name);
    }

    println!("\n💡 Use 'centrevo info --name <name>' for details");

    Ok(())
}

fn show_info(database: &PathBuf, name: &str) -> Result<()> {
    let query = QueryBuilder::new(database).context("Failed to open database")?;
    let info = query
        .get_simulation_info(name)
        .context("Failed to get simulation info")?;

    println!("\n📊 Simulation Information: {}", name);
    println!("{}", "=".repeat(50));
    println!("Created: {}", info.start_time);
    println!("Population size: {}", info.pop_size);
    println!("Generations: {}", info.num_generations);
    println!("\nParameters:");
    println!("{}", info.parameters_json);

    Ok(())
}

fn show_generations(database: &PathBuf, name: &str) -> Result<()> {
    let query = QueryBuilder::new(database).context("Failed to open database")?;
    let generations = query
        .get_recorded_generations(name)
        .context("Failed to get generations")?;

    if generations.is_empty() {
        println!("No recorded generations found for '{}'.", name);
        return Ok(());
    }

    println!("\n📈 Recorded Generations for '{}':", name);
    println!("{}", "=".repeat(50));
    println!("Generations: {:?}", generations);
    println!("Total: {} snapshots", generations.len());

    Ok(())
}

fn create_initial_population(size: usize, structure: &RepeatStructure) -> Population {
    let mut individuals = Vec::with_capacity(size);

    for i in 0..size {
        let chr1 = Chromosome::uniform(
            format!("ind{}_h1_chr1", i),
            structure.init_base,
            structure.chr_length(),
            structure.ru_length,
            structure.rus_per_hor,
            structure.alphabet.clone(),
        );

        let chr2 = Chromosome::uniform(
            format!("ind{}_h2_chr1", i),
            structure.init_base,
            structure.chr_length(),
            structure.ru_length,
            structure.rus_per_hor,
            structure.alphabet.clone(),
        );

        let h1 = Haplotype::from_chromosomes(vec![chr1]);
        let h2 = Haplotype::from_chromosomes(vec![chr2]);

        individuals.push(Individual::new(format!("ind{}", i), h1, h2));
    }

    Population::new("initial_pop", individuals)
}
