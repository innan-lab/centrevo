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

        /// Resume from checkpoint (ignores other parameters)
        #[arg(long)]
        resume: bool,

        /// Mutation rate (ignored if --resume is set)
        #[arg(long, default_value = "0.001")]
        mutation_rate: f64,

        /// Recombination break probability (ignored if --resume is set)
        #[arg(long, default_value = "0.01")]
        recomb_rate: f64,

        /// Crossover probability (given break) (ignored if --resume is set)
        #[arg(long, default_value = "0.7")]
        crossover_prob: f64,

        /// Recording interval (record every N generations) (ignored if --resume is set)
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
            resume,
            mutation_rate,
            recomb_rate,
            crossover_prob,
            record_every,
            progress,
        } => {
            run_simulation(
                &database,
                &name,
                resume,
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
        Commands::Export {
            database,
            name,
            generation,
            format,
            output,
            data_type,
        } => {
            export_data(&database, &name, generation, &format, output.as_ref(), &data_type)?;
        }
        Commands::Analyze {
            database,
            name,
            generation,
            chromosome,
            format,
            output,
        } => {
            analyze_data(&database, &name, generation, chromosome, &format, output.as_ref())?;
        }
        Commands::Validate {
            database,
            name,
            fix,
        } => {
            validate_database(&database, name.as_deref(), fix)?;
        }
        Commands::Setup { defaults } => {
            setup_wizard(defaults)?;
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
    resume: bool,
    mutation_rate: f64,
    recomb_rate: f64,
    crossover_prob: f64,
    record_every: usize,
    show_progress: bool,
) -> Result<()> {
    println!("🧬 Centrevo - Running Simulation");
    println!("============================================\n");
    
    // If resuming, load from checkpoint
    if resume {
        println!("📂 Resuming simulation from checkpoint...");
        
        // Load simulation from checkpoint
        let mut sim = Simulation::from_checkpoint(database, name)
            .map_err(|e| anyhow::anyhow!("Failed to resume: {}", e))?;
        
        let start_generation = sim.generation();
        let total_generations = sim.config().total_generations;
        
        println!("✓ Loaded checkpoint from generation {}", start_generation);
        println!("  Population size: {}", sim.config().population_size);
        println!("  Target generations: {}", total_generations);
        println!();
        
        if start_generation >= total_generations {
            println!("✓ Simulation already complete!");
            return Ok(());
        }
        
        // Setup recorder with full config
        use centrevo::storage::SimulationSnapshot;
        let snapshot = SimulationSnapshot {
            structure: sim.structure().clone(),
            mutation: sim.mutation().clone(),
            recombination: sim.recombination().clone(),
            fitness: sim.fitness().clone(),
            config: sim.config().clone(),
        };
        
        let mut recorder = Recorder::new(
            database,
            name,
            RecordingStrategy::EveryN(record_every),
        ).context("Failed to create recorder")?;
        
        // Record full config if not already recorded
        recorder.record_full_config(&snapshot)
            .context("Failed to record configuration")?;
        
        // Run simulation from checkpoint
        let remaining_generations = total_generations - start_generation;
        println!("Running {} remaining generations...", remaining_generations);
        
        for i in 1..=remaining_generations {
            let generation = start_generation + i;
            sim.step().map_err(|e| anyhow::anyhow!("Generation {}: {}", generation, e))?;
            
            // Record if needed (with RNG state for checkpoint)
            if recorder.should_record(generation) {
                let rng_state = sim.rng_state_bytes();
                recorder.record_generation(sim.population(), generation)
                    .with_context(|| format!("Failed to record generation {}", generation))?;
                recorder.record_checkpoint(generation, &rng_state)
                    .with_context(|| format!("Failed to record checkpoint {}", generation))?;
            }
            
            // Show progress
            if show_progress && (i % 10 == 0 || generation == total_generations) {
                let progress = generation as f64 / total_generations as f64 * 100.0;
                print!("\rProgress: {:.1}% ({}/{})", progress, generation, total_generations);
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
        println!("  Final generation: {}", total_generations);
        
    } else {
        // Original run logic (from generation 0)
        
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
            println!("\n💡 Use '--resume' flag to continue from checkpoint");
            return Ok(());
        }
        
        println!("Starting simulation from generation 0...\n");
        
        // Load initial population
        let initial_individuals = query.get_generation(name, 0)
            .context("Failed to load initial population. Did you run 'centrevo init' first?")?;
        
        if initial_individuals.is_empty() {
            anyhow::bail!("No initial population found. Please run 'centrevo init' first.");
        }
        
        query.close().ok();
        
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
        let mut sim = Simulation::new(structure.clone(), mutation.clone(), recombination.clone(), fitness.clone(), config.clone())
            .map_err(|e| anyhow::anyhow!(e))?;
        
        // Setup recorder
        use centrevo::storage::SimulationSnapshot;
        let snapshot = SimulationSnapshot {
            structure,
            mutation,
            recombination,
            fitness,
            config: config.clone(),
        };
        
        let mut recorder = Recorder::new(
            database,
            name,
            RecordingStrategy::EveryN(record_every),
        ).context("Failed to create recorder")?;
        
        recorder.record_full_config(&snapshot)
            .context("Failed to record configuration")?;
        
        // Run simulation
        println!("Running {} generations...", info.num_generations);
        
        for generation in 1..=info.num_generations {
            sim.step().map_err(|e| anyhow::anyhow!("Generation {}: {}", generation, e))?;
            
            // Record if needed (with RNG state for checkpoint)
            if recorder.should_record(generation) {
                let rng_state = sim.rng_state_bytes();
                recorder.record_generation(sim.population(), generation)
                    .with_context(|| format!("Failed to record generation {}", generation))?;
                recorder.record_checkpoint(generation, &rng_state)
                    .with_context(|| format!("Failed to record checkpoint {}", generation))?;
            }
            
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
    }
    
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

fn export_data(
    database: &PathBuf,
    name: &str,
    generation: usize,
    format: &str,
    output: Option<&PathBuf>,
    data_type: &str,
) -> Result<()> {
    println!("📤 Exporting data from simulation '{}'", name);
    println!("Generation: {}, Format: {}, Type: {}", generation, format, data_type);
    
    let query = QueryBuilder::new(database).context("Failed to open database")?;
    
    match data_type {
        "sequences" => export_sequences(&query, name, generation, format, output)?,
        "metadata" => export_metadata(&query, name, format, output)?,
        "fitness" => export_fitness(&query, name, format, output)?,
        _ => anyhow::bail!("Unknown data type '{}'. Use: sequences, metadata, or fitness", data_type),
    }
    
    if let Some(path) = output {
        println!("✓ Data exported to: {}", path.display());
    } else {
        println!("\n✓ Export complete");
    }
    
    Ok(())
}

fn export_sequences(
    query: &QueryBuilder,
    name: &str,
    generation: usize,
    format: &str,
    output: Option<&PathBuf>,
) -> Result<()> {
    let snapshots = query.get_generation(name, generation)
        .context("Failed to load generation")?;
    
    if snapshots.is_empty() {
        anyhow::bail!("No data found for generation {}", generation);
    }
    
    let alphabet = Alphabet::dna();
    let mut content = String::new();
    
    match format {
        "csv" => {
            content.push_str("individual_id,haplotype,chromosome,sequence\n");
            for snap in &snapshots {
                // Haplotype 1
                let seq1 = centrevo::base::Sequence::from_indices(snap.haplotype1_seq.clone(), alphabet.clone());
                content.push_str(&format!("{},h1,0,{}\n", snap.individual_id, seq1));
                
                // Haplotype 2
                let seq2 = centrevo::base::Sequence::from_indices(snap.haplotype2_seq.clone(), alphabet.clone());
                content.push_str(&format!("{},h2,0,{}\n", snap.individual_id, seq2));
            }
        }
        "fasta" => {
            for snap in &snapshots {
                // Haplotype 1
                let seq1 = centrevo::base::Sequence::from_indices(snap.haplotype1_seq.clone(), alphabet.clone());
                content.push_str(&format!(">{}|h1|chr0\n{}\n", snap.individual_id, seq1));
                
                // Haplotype 2
                let seq2 = centrevo::base::Sequence::from_indices(snap.haplotype2_seq.clone(), alphabet.clone());
                content.push_str(&format!(">{}|h2|chr0\n{}\n", snap.individual_id, seq2));
            }
        }
        "json" => {
            use serde_json::json;
            let data: Vec<_> = snapshots.iter().map(|snap| {
                let seq1 = centrevo::base::Sequence::from_indices(snap.haplotype1_seq.clone(), alphabet.clone());
                let seq2 = centrevo::base::Sequence::from_indices(snap.haplotype2_seq.clone(), alphabet.clone());
                json!({
                    "id": snap.individual_id,
                    "fitness": snap.fitness,
                    "haplotype1": [seq1.to_string()],
                    "haplotype2": [seq2.to_string()],
                })
            }).collect();
            content = serde_json::to_string_pretty(&data)?;
        }
        _ => anyhow::bail!("Unknown format '{}'. Use: csv, fasta, or json", format),
    }
    
    if let Some(path) = output {
        std::fs::write(path, content)?;
    } else {
        println!("{}", content);
    }
    
    Ok(())
}

fn export_metadata(
    query: &QueryBuilder,
    name: &str,
    format: &str,
    output: Option<&PathBuf>,
) -> Result<()> {
    let info = query.get_simulation_info(name)?;
    
    let content = match format {
        "json" => {
            use serde_json::json;
            serde_json::to_string_pretty(&json!({
                "name": name,
                "population_size": info.pop_size,
                "generations": info.num_generations,
                "start_time": info.start_time,
                "parameters": serde_json::from_str::<serde_json::Value>(&info.parameters_json)?
            }))?
        }
        "csv" => {
            format!("key,value\nname,{}\npopulation_size,{}\ngenerations,{}\nstart_time,{}\n",
                name, info.pop_size, info.num_generations, info.start_time)
        }
        _ => anyhow::bail!("Format '{}' not supported for metadata. Use: json or csv", format),
    };
    
    if let Some(path) = output {
        std::fs::write(path, content)?;
    } else {
        println!("{}", content);
    }
    
    Ok(())
}

fn export_fitness(
    query: &QueryBuilder,
    name: &str,
    format: &str,
    output: Option<&PathBuf>,
) -> Result<()> {
    let history = query.get_fitness_history(name)?;
    
    if history.is_empty() {
        anyhow::bail!("No fitness data found");
    }
    
    let content = match format {
        "csv" => {
            let mut csv = String::from("generation,mean_fitness,std_fitness,min_fitness,max_fitness\n");
            for (generation, stats) in &history {
                csv.push_str(&format!("{},{},{},{},{}\n",
                    generation, stats.mean, stats.std, stats.min, stats.max));
            }
            csv
        }
        "json" => {
            use serde_json::json;
            let data: Vec<_> = history.iter().map(|(generation, stats)| {
                json!({
                    "generation": generation,
                    "mean": stats.mean,
                    "std_dev": stats.std,
                    "min": stats.min,
                    "max": stats.max,
                })
            }).collect();
            serde_json::to_string_pretty(&data)?
        }
        _ => anyhow::bail!("Format '{}' not supported for fitness. Use: csv or json", format),
    };
    
    if let Some(path) = output {
        std::fs::write(path, content)?;
    } else {
        println!("{}", content);
    }
    
    Ok(())
}

fn analyze_data(
    database: &PathBuf,
    name: &str,
    generation: usize,
    chromosome: usize,
    format: &str,
    output: Option<&PathBuf>,
) -> Result<()> {
    use centrevo::analysis::{
        nucleotide_diversity, tajimas_d, wattersons_theta, haplotype_diversity,
        composition::gc_content, polymorphism::count_segregating_sites,
    };
    
    println!("🔬 Analyzing simulation '{}'", name);
    println!("Generation: {}, Chromosome: {}", generation, chromosome);
    
    let query = QueryBuilder::new(database).context("Failed to open database")?;
    let snapshots = query.get_generation(name, generation)
        .context("Failed to load generation")?;
    
    if snapshots.is_empty() {
        anyhow::bail!("No data found for generation {}", generation);
    }
    
    // Convert snapshots to individuals for analysis
    let alphabet = Alphabet::dna();
    let mut individuals = Vec::new();
    
    for snap in snapshots {
        let seq1 = centrevo::base::Sequence::from_indices(snap.haplotype1_seq.clone(), alphabet.clone());
        let seq2 = centrevo::base::Sequence::from_indices(snap.haplotype2_seq.clone(), alphabet.clone());
        
        let chr1 = Chromosome::new(snap.haplotype1_chr_id, seq1, 171, 12);
        let chr2 = Chromosome::new(snap.haplotype2_chr_id, seq2, 171, 12);
        
        let h1 = Haplotype::from_chromosomes(vec![chr1]);
        let h2 = Haplotype::from_chromosomes(vec![chr2]);
        
        let mut ind = Individual::new(snap.individual_id, h1, h2);
        ind.set_fitness(snap.fitness);
        individuals.push(ind);
    }
    
    let population = Population::new(format!("{}_gen{}", name, generation), individuals);
    let seq_len = population.individuals()[0].haplotype1().get(chromosome)
        .map(|chr| chr.sequence().len())
        .unwrap_or(0);
    
    // Calculate metrics
    let pi = nucleotide_diversity(&population, chromosome);
    let tajima = tajimas_d(&population, chromosome);
    let theta_w = wattersons_theta(&population, chromosome);
    let hap_div = haplotype_diversity(&population, chromosome);
    let seg_sites = count_segregating_sites(&population, chromosome, 0)
        + count_segregating_sites(&population, chromosome, 1);
    
    // Calculate GC content at population level
    let gc = gc_content(&population, None, None, None);
    
    let content = match format {
        "pretty" => {
            format!(
                "\n📊 Population Genetics Summary\n\
                 ================================\n\
                 Population size: {}\n\
                 Sequences analyzed: {} (2n)\n\
                 Sequence length: {} bp\n\
                 \n\
                 Diversity Metrics:\n\
                 ------------------\n\
                 Nucleotide diversity (π): {:.6}\n\
                 Watterson's theta (θ_W): {:.6}\n\
                 Tajima's D: {:.4}\n\
                 Haplotype diversity: {:.4}\n\
                 \n\
                 Polymorphism:\n\
                 -------------\n\
                 Segregating sites: {}\n\
                 \n\
                 Composition:\n\
                 ------------\n\
                 Mean GC content: {:.2}%\n",
                population.size(),
                population.size() * 2,
                seq_len,
                pi,
                theta_w,
                tajima,
                hap_div,
                seg_sites,
                gc * 100.0,
            )
        }
        "json" => {
            use serde_json::json;
            serde_json::to_string_pretty(&json!({
                "simulation": name,
                "generation": generation,
                "chromosome": chromosome,
                "population_size": population.size(),
                "sequence_count": population.size() * 2,
                "sequence_length": seq_len,
                "diversity": {
                    "pi": pi,
                    "theta_w": theta_w,
                    "tajima_d": tajima,
                    "haplotype_diversity": hap_div,
                },
                "polymorphism": {
                    "segregating_sites": seg_sites,
                },
                "composition": {
                    "mean_gc_content": gc,
                }
            }))?
        }
        _ => anyhow::bail!("Unknown format '{}'. Use: pretty or json", format),
    };
    
    if let Some(path) = output {
        std::fs::write(path, content)?;
    } else {
        println!("{}", content);
    }
    
    Ok(())
}

fn validate_database(
    database: &PathBuf,
    name: Option<&str>,
    fix: bool,
) -> Result<()> {
    println!("🔍 Validating database: {}", database.display());
    
    if !database.exists() {
        anyhow::bail!("Database file does not exist");
    }
    
    let query = QueryBuilder::new(database).context("Failed to open database")?;
    
    let simulations = if let Some(sim_name) = name {
        vec![sim_name.to_string()]
    } else {
        query.list_simulations().context("Failed to list simulations")?
    };
    
    if simulations.is_empty() {
        println!("⚠️  No simulations found in database");
        return Ok(());
    }
    
    let mut total_issues = 0;
    
    for sim_name in &simulations {
        println!("\nValidating simulation: {}", sim_name);
        println!("{}", "-".repeat(50));
        
        // Check metadata
        let info = match query.get_simulation_info(sim_name) {
            Ok(info) => {
                println!("✓ Metadata: OK");
                info
            }
            Err(e) => {
                println!("✗ Metadata: FAILED - {}", e);
                total_issues += 1;
                continue;
            }
        };
        
        // Check recorded generations
        let recorded_gens = match query.get_recorded_generations(sim_name) {
            Ok(gens) => {
                println!("✓ Recorded generations: {} snapshots", gens.len());
                gens
            }
            Err(e) => {
                println!("✗ Failed to query generations: {}", e);
                total_issues += 1;
                continue;
            }
        };
        
        if recorded_gens.is_empty() {
            println!("⚠️  No generation data recorded");
            total_issues += 1;
            continue;
        }
        
        // Check for gaps in generations
        let expected_gens: Vec<usize> = (0..=info.num_generations).step_by(100).collect();
        let missing: Vec<usize> = expected_gens.iter()
            .filter(|g| !recorded_gens.contains(g))
            .copied()
            .collect();
        
        if !missing.is_empty() {
            println!("⚠️  Missing generations (expected every 100): {:?}", 
                if missing.len() > 10 {
                    format!("{} generations missing", missing.len())
                } else {
                    format!("{:?}", missing)
                });
            total_issues += 1;
        } else {
            println!("✓ Generation continuity: OK");
        }
        
        // Check if we can load a sample generation
        if let Some(&generation_num) = recorded_gens.first() {
            match query.get_generation(sim_name, generation_num) {
                Ok(snapshots) => {
                    if snapshots.is_empty() {
                        println!("⚠️  Generation {} has no individuals", generation_num);
                        total_issues += 1;
                    } else if snapshots.len() != info.pop_size {
                        println!("⚠️  Generation {} has {} individuals (expected {})", 
                            generation_num, snapshots.len(), info.pop_size);
                        total_issues += 1;
                    } else {
                        println!("✓ Population data: OK ({} individuals)", snapshots.len());
                    }
                }
                Err(e) => {
                    println!("✗ Failed to load generation {}: {}", generation_num, e);
                    total_issues += 1;
                }
            }
        }
        
        // Check fitness history
        match query.get_fitness_history(sim_name) {
            Ok(history) => {
                if history.is_empty() {
                    println!("⚠️  No fitness history recorded");
                } else {
                    println!("✓ Fitness history: {} entries", history.len());
                }
            }
            Err(e) => {
                println!("⚠️  Fitness history query failed: {}", e);
            }
        }
    }
    
    println!("\n{}", "=".repeat(50));
    if total_issues == 0 {
        println!("✓ Validation complete: No issues found");
    } else {
        println!("⚠️  Validation complete: {} issue(s) found", total_issues);
        if !fix {
            println!("💡 Use --fix flag to attempt automatic repairs");
        }
    }
    
    if fix && total_issues > 0 {
        println!("\n⚠️  Automatic repair not yet implemented");
        println!("Please check the issues manually or re-run the simulation");
    }
    
    Ok(())
}

fn setup_wizard(use_defaults: bool) -> Result<()> {
    println!("\n🧙 Centrevo Setup Wizard");
    println!("========================\n");
    println!("This wizard will guide you through setting up and running a simulation.\n");

    // Simulation configuration
    let name = if use_defaults {
        "simulation".to_string()
    } else {
        prompt_string("Simulation name", Some("simulation"))?
    };

    let database = if use_defaults {
        PathBuf::from("simulation.db")
    } else {
        PathBuf::from(prompt_string("Database file", Some("simulation.db"))?)
    };

    println!("\n📐 Population Parameters");
    println!("-----------------------");
    
    let population_size = if use_defaults {
        100
    } else {
        prompt_usize("Population size", Some(100))?
    };

    let generations = if use_defaults {
        1000
    } else {
        prompt_usize("Number of generations", Some(1000))?
    };

    println!("\n🧬 Genome Structure");
    println!("------------------");
    
    let ru_length = if use_defaults {
        171
    } else {
        prompt_usize("Repeat unit (RU) length (bp)", Some(171))?
    };

    let rus_per_hor = if use_defaults {
        12
    } else {
        prompt_usize("RUs per Higher-Order Repeat (HOR)", Some(12))?
    };

    let hors_per_chr = if use_defaults {
        100
    } else {
        prompt_usize("HORs per chromosome", Some(100))?
    };

    println!("\n🧪 Evolution Parameters");
    println!("----------------------");
    
    let mutation_rate = if use_defaults {
        0.001
    } else {
        prompt_f64("Mutation rate", Some(0.001))?
    };

    let recomb_rate = if use_defaults {
        0.01
    } else {
        prompt_f64("Recombination break probability", Some(0.01))?
    };

    let crossover_prob = if use_defaults {
        0.7
    } else {
        prompt_f64("Crossover probability (given break)", Some(0.7))?
    };

    println!("\n💾 Recording Options");
    println!("-------------------");
    
    let record_every = if use_defaults {
        100
    } else {
        prompt_usize("Record every N generations", Some(100))?
    };

    let seed = if use_defaults {
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
    println!("Simulation name: {}", name);
    println!("Database: {}", database.display());
    println!();
    println!("Population:");
    println!("  Size: {} individuals", population_size);
    println!("  Generations: {}", generations);
    println!();
    println!("Genome:");
    println!("  RU length: {} bp", ru_length);
    println!("  RUs per HOR: {}", rus_per_hor);
    println!("  HORs per chromosome: {}", hors_per_chr);
    println!("  Total chromosome length: {} bp ({:.1} kb)", total_bp, total_bp as f64 / 1000.0);
    println!();
    println!("Evolution:");
    println!("  Mutation rate: {}", mutation_rate);
    println!("  Recombination rate: {}", recomb_rate);
    println!("  Crossover probability: {}", crossover_prob);
    println!();
    println!("Recording:");
    println!("  Every {} generations", record_every);
    println!("  Total snapshots: ~{}", (generations / record_every) + 1);
    if let Some(s) = seed {
        println!("  Seed: {}", s);
    } else {
        println!("  Seed: random");
    }

    // Confirm
    if !use_defaults {
        println!();
        if !prompt_confirm("Proceed with simulation?", true)? {
            println!("\n❌ Setup cancelled.");
            return Ok(());
        }
    }

    // Initialize simulation
    println!("\n🚀 Starting simulation setup...\n");
    
    init_simulation(
        &name,
        &database,
        population_size,
        generations,
        ru_length,
        rus_per_hor,
        hors_per_chr,
        seed,
    )?;

    // Ask if user wants to run now
    let should_run = if use_defaults {
        true
    } else {
        println!();
        prompt_confirm("Run simulation now?", true)?
    };

    if should_run {
        println!();
        run_simulation(
            &database,
            &name,
            false, // not resume
            mutation_rate,
            recomb_rate,
            crossover_prob,
            record_every,
            true, // show progress
        )?;
    } else {
        println!("\n💡 To run later, use: centrevo run -N {}", name);
    }

    Ok(())
}

fn prompt_string(prompt: &str, default: Option<&str>) -> Result<String> {
    print!("{}", prompt);
    if let Some(def) = default {
        print!(" [{}]", def);
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
    print!("{}", prompt);
    if let Some(def) = default {
        print!(" [{}]", def);
    }
    print!(": ");
    io::stdout().flush()?;

    let mut input = String::new();
    io::stdin().read_line(&mut input)?;
    let input = input.trim();

    if input.is_empty() {
        default.ok_or_else(|| anyhow::anyhow!("Input required"))
    } else {
        input.parse::<usize>()
            .with_context(|| format!("Invalid number: {}", input))
    }
}

fn prompt_f64(prompt: &str, default: Option<f64>) -> Result<f64> {
    print!("{}", prompt);
    if let Some(def) = default {
        print!(" [{}]", def);
    }
    print!(": ");
    io::stdout().flush()?;

    let mut input = String::new();
    io::stdin().read_line(&mut input)?;
    let input = input.trim();

    if input.is_empty() {
        default.ok_or_else(|| anyhow::anyhow!("Input required"))
    } else {
        input.parse::<f64>()
            .with_context(|| format!("Invalid number: {}", input))
    }
}

fn prompt_confirm(prompt: &str, default: bool) -> Result<bool> {
    let default_str = if default { "Y/n" } else { "y/N" };
    print!("{} [{}]: ", prompt, default_str);
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
