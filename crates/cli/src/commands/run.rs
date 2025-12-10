use anyhow::{Context, Result};
use centrevo_sim::simulation::Simulation;
use centrevo_sim::storage::{QueryBuilder, Recorder, RecordingStrategy, SimulationSnapshot};
use indicatif::{ProgressBar, ProgressStyle};
use std::path::PathBuf;

use crate::printing::print_simulation_parameters;

#[allow(clippy::too_many_arguments)]
pub fn run_simulation(
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
            println!("⚠️  Warning: Seed override ignored when resuming from checkpoint.");
        }

        // Apply record_every override
        let record_every = record_every_override.unwrap_or(100);

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
        print_simulation_parameters(&sim);

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

        // Print all parameters
        println!("Configuration:");
        print_simulation_parameters(&sim);

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
