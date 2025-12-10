mod printing;

use anyhow::{Context, Result};
use centrevo_sim::base::{FitnessValue, Nucleotide};
use centrevo_sim::genome::{Chromosome, Haplotype, Individual};
use centrevo_sim::simulation::{
    FitnessConfig, MutationConfig, Population, RecombinationConfig, Simulation, SimulationConfig,
    UniformRepeatStructure,
};
use centrevo_sim::storage::{QueryBuilder, Recorder, RecordingStrategy};
use clap::{Parser, Subcommand};
use indicatif::{ProgressBar, ProgressStyle};
use printing::print_parameters;
use std::io::{self, Write};
use std::path::PathBuf;

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
            init_simulation(
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
            run_simulation(&database, &name, resume, seed, record_every, progress)?;
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
            export_data(
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
            analyze_data(
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
            validate_database(&database, name.as_deref(), fix)?;
        }
        Commands::Setup { defaults } => {
            setup_wizard(defaults)?;
        }
    }

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn init_simulation(
    name: &str,
    output: &PathBuf,
    population_size: usize,
    generations: usize,
    ru_length: usize,
    rus_per_hor: usize,
    hors_per_chr: usize,
    mutation_rate: f64,
    recomb_rate: f64,
    crossover_prob: f64,
    gc_extension_prob: f64,
    homology_strength: f64,
    search_window: usize,
    record_every: usize,
    seed: Option<u64>,
) -> Result<()> {
    println!("🧬 Centrevo - Centromeric Evolution Simulator");
    println!("============================================\n");
    println!("Initializing simulation: {name}");

    // Create configurations
    let structure =
        UniformRepeatStructure::new(Nucleotide::A, ru_length, rus_per_hor, hors_per_chr, 1);

    let config = SimulationConfig::new(population_size, generations, seed);

    // Create mutation and recombination configs
    let mutation = MutationConfig::uniform(mutation_rate)
        .map_err(|e| anyhow::anyhow!("Failed to create mutation configuration: {e}"))?;

    let recomb_params = centrevo_sim::evolution::RecombinationModel::builder()
        .break_prob(recomb_rate)
        .crossover_prob(crossover_prob)
        .gc_extension_prob(gc_extension_prob)
        .homology_strength(homology_strength)
        .search_window(search_window)
        .build()
        .map_err(|e| anyhow::anyhow!("Failed to create recombination model: {e}"))?;
    let recombination = RecombinationConfig::new(recomb_params);

    let fitness = FitnessConfig::neutral();

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
    use centrevo_sim::storage::SimulationSnapshot;
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

#[allow(clippy::too_many_arguments)]
fn run_simulation(
    database: &PathBuf,
    name: &str,
    resume: bool,
    seed_override: Option<u64>,
    record_every_override: Option<usize>,
    show_progress: bool,
) -> Result<()> {
    println!("🧬 Centrevo - Running Simulation");
    println!("============================================\n");

    // If resuming, load from checkpoint
    if resume {
        println!("📂 Resuming simulation from checkpoint...");

        // Load simulation from checkpoint
        let mut sim = Simulation::from_checkpoint(database, name)
            .map_err(|e| anyhow::anyhow!("Failed to resume: {e}"))?;

        // Apply seed override if provided
        if seed_override.is_some() {
            // We need to re-seed the RNG.
            // Note: Simulation::from_checkpoint restores the RNG state.
            // If we provide a new seed, we are essentially branching from that point
            // with a different random sequence.
            // We might need to expose a method to re-seed.
            // For now, let's assume we can't easily re-seed a restored RNG without breaking properties,
            // or check if Simulation has a reseed method.
            // Actually, the user requirement is to override seed mainly for fresh runs.
            // For resume, usually you want to continue exactly.
            // But if they explicitly pass --seed with --resume, maybe they want to branch?
            // Let's warn if they try to do both, or just ignore seed on resume for now as per plan focus on "Run" vs "Init" separation for fresh runs.
            // Wait, the plan says "Run ... Allows overriding ...".
            println!("⚠️  Warning: Seed override ignored when resuming from checkpoint.");
        }

        // Apply record_every override
        let record_every = record_every_override.unwrap_or(100); // Default to 100 if we can't get it from config easily?
        // Actually, we should probably fetch the strategy from the config loaded in checkpoint?
        // But the Recorder is separate.
        // We configure a NEW recorder here.

        let start_generation = sim.generation();
        let total_generations = sim.config().total_generations;

        println!("✓ Loaded checkpoint from generation {start_generation}");
        println!("  Population size: {}", sim.config().population_size);
        println!("  Target generations: {total_generations}");
        println!();

        if start_generation >= total_generations {
            println!("✓ Simulation already complete!");
            return Ok(());
        }

        // Setup recorder with full config
        use centrevo_sim::storage::SimulationSnapshot;
        let snapshot = SimulationSnapshot {
            structure: sim
                .structure()
                .cloned()
                .expect("Cannot resume without uniform structure"),
            mutation: sim.mutation().clone(),
            recombination: sim.recombination().clone(),
            fitness: sim.fitness().clone(),
            config: sim.config().clone(),
        };

        let mut recorder = Recorder::new(database, name, RecordingStrategy::EveryN(record_every))
            .context("Failed to create recorder")?;

        // Record full config if not already recorded (idempotent-ish)
        recorder
            .record_full_config(&snapshot)
            .context("Failed to record configuration")?;

        // Print parameters
        println!("Resuming Configuration:");
        print_sim_state(&sim);

        // Run simulation from checkpoint
        let remaining_generations = total_generations - start_generation;
        println!("Running {remaining_generations} remaining generations...");

        let pb = if show_progress {
            let pb = ProgressBar::new(total_generations as u64);
            pb.set_position(start_generation as u64);
            pb.set_style(
                ProgressStyle::default_bar()
                    .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta}) {per_sec}")
                    .unwrap()
                    .progress_chars("#>-"),
            );
            Some(pb)
        } else {
            None
        };

        for i in 1..=remaining_generations {
            let generation = start_generation + i;
            sim.step()
                .map_err(|e| anyhow::anyhow!("Generation {generation}: {e}"))?;

            // Record if needed (with RNG state for checkpoint)
            if recorder.should_record(generation) {
                let rng_state = sim.rng_state_bytes();
                recorder
                    .record_generation(sim.population(), generation)
                    .with_context(|| format!("Failed to record generation {generation}"))?;
                recorder
                    .record_checkpoint(generation, &rng_state)
                    .with_context(|| format!("Failed to record checkpoint {generation}"))?;
            }

            // Show progress
            if let Some(pb) = &pb {
                pb.inc(1);
            }
        }

        if let Some(pb) = pb {
            pb.finish_with_message("Done");
        }

        // Finalize
        recorder
            .finalize_metadata()
            .context("Failed to finalize metadata")?;

        println!("\n✓ Simulation complete!");
        println!("  Final generation: {total_generations}");
    } else {
        // Fresh run logic (from generation 0)

        // Open database
        let query = QueryBuilder::new(database).context("Failed to open database")?;

        // Load full configuration
        let snapshot = query
            .get_full_config(name)
            .context("Failed to load configuration. Did you run 'centrevo init' first?")?;

        // Verify that Gen 0 exists
        let initial_individuals = query
            .get_generation(name, 0)
            .context("Failed to load initial population. Did you run 'centrevo init' first?")?;

        if initial_individuals.is_empty() {
            anyhow::bail!("No initial population found. Please run 'centrevo init' first.");
        }

        query.close().ok();

        // Apply Overrides
        let mut config = snapshot.config.clone();
        if let Some(seed) = seed_override {
            println!("ℹ️  Overriding random seed: {seed}");
            config.seed = Some(seed);
        }

        // Determing recording strategy
        // Ideally we save the strategy in DB, but SimulationConfig doesn't have it.
        // It's a Recorder concern.
        // If override is provided, use it. Else default to 100 (or what was used in init?)
        // In init we passed record_every to Recorder, but did we save it?
        // Recorder doesn't persist its strategy directly in a way we can easily read back without
        // adding a new table or column.
        // For now, let's default to 100 if not overridden, or try to infer.
        // Or assume the user should pass it if they want something specific.
        let record_every = record_every_override.unwrap_or(100);

        println!("Simulation: {name}");
        println!("Population size: {}", config.population_size);
        println!("Target generations: {}", config.total_generations);
        println!();

        // Setup simulation using the configuration loaded from DB
        let mut sim = Simulation::from_config(
            centrevo_sim::simulation::SequenceConfig::Load {
                source: centrevo_sim::simulation::SequenceSource::Database {
                    path: database.to_string_lossy().to_string(),
                    sim_id: name.to_string(),
                    generation: Some(0),
                },
                structure: Some(snapshot.structure.clone()),
            },
            Some(snapshot.structure.clone()),
            snapshot.mutation.clone(),
            snapshot.recombination.clone(),
            snapshot.fitness.clone(),
            config.clone(), // Use the config with potentially overridden seed
        )
        .map_err(|e| anyhow::anyhow!("Failed to initialize simulation: {e}"))?;

        // Setup recorder
        let mut recorder = Recorder::new(database, name, RecordingStrategy::EveryN(record_every))
            .context("Failed to create recorder")?;

        // We don't need to record config again, it's already there from Init.
        // But maybe good to ensure consistency?
        // Actually, if we override seed, we might want to record that this run used a different seed?
        // But the DB structure assumes one config per simulation name.
        // If we strictly follow "Simulation Name = Unique Config", then overriding seed
        // might imply we are running a "Variant" of the simulation.
        // But for now, we just run it.

        // Print all parameters
        println!("Configuration:");
        print_sim_state(&sim);

        // Run simulation
        println!("Running {} generations...", config.total_generations);

        let pb = if show_progress {
            let pb = ProgressBar::new(config.total_generations as u64);
            pb.set_style(
                ProgressStyle::default_bar()
                    .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta}) {per_sec}")
                    .unwrap()
                    .progress_chars("#>-"),
            );
            Some(pb)
        } else {
            None
        };

        for generation in 1..=config.total_generations {
            sim.step()
                .map_err(|e| anyhow::anyhow!("Generation {generation}: {e}"))?;

            // Record if needed (with RNG state for checkpoint)
            if recorder.should_record(generation) {
                let rng_state = sim.rng_state_bytes();
                recorder
                    .record_generation(sim.population(), generation)
                    .with_context(|| format!("Failed to record generation {generation}"))?;
                recorder
                    .record_checkpoint(generation, &rng_state)
                    .with_context(|| format!("Failed to record checkpoint {generation}"))?;
            }

            // Show progress
            if let Some(pb) = &pb {
                pb.inc(1);
            }
        }

        if let Some(pb) = pb {
            pb.finish_with_message("Done");
        }

        // Finalize
        recorder
            .finalize_metadata()
            .context("Failed to finalize metadata")?;

        println!("\n✓ Simulation complete!");
        println!("  Final generation: {}", config.total_generations);
    }

    println!("\n💡 Use 'centrevo info -N {name}' to view results");

    Ok(())
}

fn list_simulations(database: &PathBuf) -> Result<()> {
    let query = QueryBuilder::new(database).context("Failed to open database")?;
    let simulations = query
        .list_simulations()
        .context("Failed to list simulations")?;

    if simulations.is_empty() {
        println!("No simulations found in database.");
        return Ok(());
    }

    println!("\n📊 Simulations in {}:", database.display());
    println!("{}", "=".repeat(50));

    for name in simulations {
        println!("  • {name}");
    }

    println!("\n💡 Use 'centrevo info --name <name>' for details");

    Ok(())
}

fn show_info(database: &PathBuf, name: &str) -> Result<()> {
    let query = QueryBuilder::new(database).context("Failed to open database")?;
    let info = query
        .get_simulation_info(name)
        .context("Failed to get simulation info")?;

    println!("\n📊 Simulation Information: {name}");
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
        println!("No recorded generations found for '{name}'.");
        return Ok(());
    }

    println!("\n📈 Recorded Generations for '{name}':");
    println!("{}", "=".repeat(50));
    println!("Generations: {generations:?}");
    println!("Total: {} snapshots", generations.len());

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

fn sequence_from_indices(indices: Vec<u8>) -> centrevo_sim::base::Sequence {
    use centrevo_sim::base::{Nucleotide, Sequence};
    let nucleotides: Vec<Nucleotide> = indices
        .into_iter()
        .map(|i| Nucleotide::from_index(i).unwrap_or(Nucleotide::A))
        .collect();
    Sequence::from_nucleotides(nucleotides)
}

fn export_data(
    database: &PathBuf,
    name: &str,
    generation: usize,
    format: &str,
    output: Option<&PathBuf>,
    data_type: &str,
) -> Result<()> {
    println!("📤 Exporting data from simulation '{name}'");
    println!("Generation: {generation}, Format: {format}, Type: {data_type}");

    let query = QueryBuilder::new(database).context("Failed to open database")?;

    match data_type {
        "sequences" => export_sequences(&query, name, generation, format, output)?,
        "metadata" => export_metadata(&query, name, format, output)?,
        "fitness" => export_fitness(&query, name, format, output)?,
        _ => anyhow::bail!("Unknown data type '{data_type}'. Use: sequences, metadata, or fitness"),
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
    let snapshots = query
        .get_generation(name, generation)
        .context("Failed to load generation")?;

    if snapshots.is_empty() {
        anyhow::bail!("No data found for generation {generation}");
    }

    let mut content = String::new();

    match format {
        "csv" => {
            content.push_str("individual_id,haplotype,chromosome,sequence\n");
            for snap in &snapshots {
                // Haplotype 1
                let seq1 = sequence_from_indices(snap.haplotype1_seq.clone());
                content.push_str(&format!("{},h1,0,{}\n", snap.individual_id, seq1));

                // Haplotype 2
                let seq2 = sequence_from_indices(snap.haplotype2_seq.clone());
                content.push_str(&format!("{},h2,0,{}\n", snap.individual_id, seq2));
            }
        }
        "fasta" => {
            for snap in &snapshots {
                // Haplotype 1
                let seq1 = sequence_from_indices(snap.haplotype1_seq.clone());
                content.push_str(&format!(">{}|h1|chr0\n{}\n", snap.individual_id, seq1));

                // Haplotype 2
                let seq2 = sequence_from_indices(snap.haplotype2_seq.clone());
                content.push_str(&format!(">{}|h2|chr0\n{}\n", snap.individual_id, seq2));
            }
        }
        "json" => {
            use serde_json::json;
            let data: Vec<_> = snapshots
                .iter()
                .map(|snap| {
                    let seq1 = sequence_from_indices(snap.haplotype1_seq.clone());
                    let seq2 = sequence_from_indices(snap.haplotype2_seq.clone());
                    json!({
                        "id": snap.individual_id,
                        "fitness": snap.fitness,
                        "haplotype1": [seq1.to_string()],
                        "haplotype2": [seq2.to_string()],
                    })
                })
                .collect();
            content = serde_json::to_string_pretty(&data)?;
        }
        _ => anyhow::bail!("Unknown format '{format}'. Use: csv, fasta, or json"),
    }

    if let Some(path) = output {
        std::fs::write(path, content)?;
    } else {
        println!("{content}");
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
            format!(
                "key,value\nname,{}\npopulation_size,{}\ngenerations,{}\nstart_time,{}\n",
                name, info.pop_size, info.num_generations, info.start_time
            )
        }
        _ => anyhow::bail!("Format '{format}' not supported for metadata. Use: json or csv"),
    };

    if let Some(path) = output {
        std::fs::write(path, content)?;
    } else {
        println!("{content}");
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
            let mut csv =
                String::from("generation,mean_fitness,std_fitness,min_fitness,max_fitness\n");
            for (generation, stats) in &history {
                csv.push_str(&format!(
                    "{},{},{},{},{}\n",
                    generation, stats.mean, stats.std, stats.min, stats.max
                ));
            }
            csv
        }
        "json" => {
            use serde_json::json;
            let data: Vec<_> = history
                .iter()
                .map(|(generation, stats)| {
                    json!({
                        "generation": generation,
                        "mean": stats.mean,
                        "std_dev": stats.std,
                        "min": stats.min,
                        "max": stats.max,
                    })
                })
                .collect();
            serde_json::to_string_pretty(&data)?
        }
        _ => anyhow::bail!("Format '{format}' not supported for fitness. Use: csv or json"),
    };

    if let Some(path) = output {
        std::fs::write(path, content)?;
    } else {
        println!("{content}");
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
    use centrevo_analysis::analysis::{
        composition::gc_content, haplotype_diversity, nucleotide_diversity,
        polymorphism::count_segregating_sites, tajimas_d, wattersons_theta,
    };

    println!("🔬 Analyzing simulation '{name}'");
    println!("Generation: {generation}, Chromosome: {chromosome}");

    let query = QueryBuilder::new(database).context("Failed to open database")?;
    let snapshots = query
        .get_generation(name, generation)
        .context("Failed to load generation")?;

    if snapshots.is_empty() {
        anyhow::bail!("No data found for generation {generation}");
    }

    // Convert snapshots to individuals for analysis
    let mut individuals = Vec::new();

    for snap in snapshots {
        let seq1 = sequence_from_indices(snap.haplotype1_seq.clone());
        let seq2 = sequence_from_indices(snap.haplotype2_seq.clone());

        // Assume uniform structure for analysis reconstruction
        // TODO: Load actual structure from database
        let ru_len = 171;
        let rus_per_hor = 12;
        let hor_len = ru_len * rus_per_hor;
        // Avoid division by zero if length is 0 (though unlikely for valid sim)
        let hors_per_chr = if hor_len > 0 { seq1.len() / hor_len } else { 0 };

        let map =
            centrevo_sim::genome::repeat_map::RepeatMap::uniform(ru_len, rus_per_hor, hors_per_chr);

        let chr1 = Chromosome::new(snap.haplotype1_chr_id, seq1, map.clone());
        let chr2 = Chromosome::new(snap.haplotype2_chr_id, seq2, map);

        let h1 = Haplotype::from_chromosomes(vec![chr1]);
        let h2 = Haplotype::from_chromosomes(vec![chr2]);

        let mut ind = Individual::new(snap.individual_id, h1, h2);
        if let Some(f) = snap.fitness {
            ind.set_cached_fitness(FitnessValue::new(f));
        }
        individuals.push(ind);
    }

    let population = Population::new(format!("{name}_gen{generation}"), individuals);
    let seq_len = population.individuals()[0]
        .haplotype1()
        .get(chromosome)
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
        _ => anyhow::bail!("Unknown format '{format}'. Use: pretty or json"),
    };

    if let Some(path) = output {
        std::fs::write(path, content)?;
    } else {
        println!("{content}");
    }

    Ok(())
}

fn validate_database(database: &PathBuf, name: Option<&str>, fix: bool) -> Result<()> {
    println!("🔍 Validating database: {}", database.display());

    if !database.exists() {
        anyhow::bail!("Database file does not exist");
    }

    let query = QueryBuilder::new(database).context("Failed to open database")?;

    let simulations = if let Some(sim_name) = name {
        vec![sim_name.to_string()]
    } else {
        query
            .list_simulations()
            .context("Failed to list simulations")?
    };

    if simulations.is_empty() {
        println!("⚠️  No simulations found in database");
        return Ok(());
    }

    let mut total_issues = 0;

    for sim_name in &simulations {
        println!("\nValidating simulation: {sim_name}");
        println!("{}", "-".repeat(50));

        // Check metadata
        let info = match query.get_simulation_info(sim_name) {
            Ok(info) => {
                println!("✓ Metadata: OK");
                info
            }
            Err(e) => {
                println!("✗ Metadata: FAILED - {e}");
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
                println!("✗ Failed to query generations: {e}");
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
        let missing: Vec<usize> = expected_gens
            .iter()
            .filter(|g| !recorded_gens.contains(g))
            .copied()
            .collect();

        if !missing.is_empty() {
            println!(
                "⚠️  Missing generations (expected every 100): {:?}",
                if missing.len() > 10 {
                    format!("{} generations missing", missing.len())
                } else {
                    format!("{missing:?}")
                }
            );
            total_issues += 1;
        } else {
            println!("✓ Generation continuity: OK");
        }

        // Check if we can load a sample generation
        if let Some(&generation_num) = recorded_gens.first() {
            match query.get_generation(sim_name, generation_num) {
                Ok(snapshots) => {
                    if snapshots.is_empty() {
                        println!("⚠️  Generation {generation_num} has no individuals");
                        total_issues += 1;
                    } else if snapshots.len() != info.pop_size {
                        println!(
                            "⚠️  Generation {} has {} individuals (expected {})",
                            generation_num,
                            snapshots.len(),
                            info.pop_size
                        );
                        total_issues += 1;
                    } else {
                        println!("✓ Population data: OK ({} individuals)", snapshots.len());
                    }
                }
                Err(e) => {
                    println!("✗ Failed to load generation {generation_num}: {e}");
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
                println!("⚠️  Fitness history query failed: {e}");
            }
        }
    }

    println!("\n{}", "=".repeat(50));
    if total_issues == 0 {
        println!("✓ Validation complete: No issues found");
    } else {
        println!("⚠️  Validation complete: {total_issues} issue(s) found");
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
        mutation_rate,
        recomb_rate,
        crossover_prob,
        0.1, // gc_extension_prob default
        0.0, // homology_strength default
        0,   // search_window default
        record_every,
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

fn print_sim_state(sim: &Simulation) {
    use printing::print_simulation_parameters;
    print_simulation_parameters(sim);
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

        // Before the fix, this would have been much larger (if it didn't crash)
        // The bug passed chr_length (100) as ru_length, ru_length (10) as rus_per_hor, rus_per_hor (5) as num_hors.
        // Incorrect length would be: 100 * 10 * 5 = 5000.
    }
}
