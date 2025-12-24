use anyhow::{Context, Result};
use centrevo_sim::simulation::{Configuration, Simulation};
use centrevo_sim::storage::{QueryBuilder, Recorder, RecordingStrategy};
use indicatif::{ProgressBar, ProgressStyle};
use std::path::PathBuf;

use crate::printing::print_simulation_parameters;

#[allow(clippy::too_many_arguments)]
pub fn run_simulation(
    database: &PathBuf,
    resume: bool,
    seed_override: Option<u64>,
    record_every_override: Option<usize>,
    show_progress: bool,
) -> Result<()> {
    println!("🧬 Centrevo - Running Simulation");
    println!("============================================\n");

    let rt = tokio::runtime::Runtime::new().context("Failed to create Tokio runtime")?;

    // If resuming, load from checkpoint
    if resume {
        println!("📂 Resuming simulation from checkpoint...");

        // Load simulation from checkpoint
        let mut sim = Simulation::from_checkpoint(database)
            .map_err(|e| anyhow::anyhow!("Failed to resume: {e}"))?;

        // Apply seed override if provided
        if seed_override.is_some() {
            println!("⚠️  Warning: Seed override ignored when resuming from checkpoint.");
        }

        // Apply record_every override
        let record_every = record_every_override.unwrap_or(100);

        let start_generation = sim.generation();
        let total_generations = sim.simulation_config().total_generations;

        println!("✓ Loaded checkpoint from generation {start_generation}");
        println!("  Population size: {}", sim.simulation_config().population_size);
        println!("  Target generations: {total_generations}");
        println!();

        if start_generation >= total_generations {
            println!("✓ Simulation already complete!");
            return Ok(());
        }

        rt.block_on(async {
            // Setup recorder with full config
            // Use the simulation's configuration directly
            let config = sim.configuration();

            let buffer_config = centrevo_sim::storage::BufferConfig {
                compression_level: 0,
                ..Default::default()
            };

            let recorder = Recorder::new(
                database,
                config,
                buffer_config,
                config.execution.codec,
            )
            .context("Failed to create recorder")?;
            // Metadata is recorded in Recorder::new

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
                        .template(
                            "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta}) {per_sec}",
                        )
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

                // Determine if we need to record
                let should_record =
                    RecordingStrategy::EveryN(record_every).should_record(generation);

                // Record if needed (with RNG state for checkpoint)
                if should_record {
                    let rng_state = sim.rng_state_bytes();
                    recorder
                        .record_generation(sim.population(), generation, Some(rng_state), sim.arena())
                        .await
                        .context(format!("Failed to record generation {generation}"))?;
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
            // Metadata finalized implicitly on close or not needed


            recorder.close().await.context("Failed to close recorder")?;

            Ok::<(), anyhow::Error>(())
        })?;

        println!("\n✓ Simulation complete!");
        println!("  Final generation: {total_generations}");
    } else {
        // Fresh run logic (from generation 0)

        // Open database
        let query = QueryBuilder::new(database).context("Failed to open database")?;

        // Load full configuration (now returns Configuration directly)
        let config = query
            .get_full_config()
            .context("Failed to load configuration. Did you run 'centrevo init' first?")?;

        // Verify that Gen 0 exists
        let initial_individuals = query
            .get_generation(0)
            .context("Failed to load initial population. Did you run 'centrevo init' first?")?;

        if initial_individuals.is_empty() {
            anyhow::bail!("No initial population found. Please run 'centrevo init' first.");
        }

        query.close().ok();

        // Apply Overrides
        let record_every = record_every_override.unwrap_or(100);

        let initialization_config = centrevo_sim::simulation::InitializationConfig::Load {
            source: centrevo_sim::simulation::SequenceSource::Database {
                path: database.to_string_lossy().to_string(),
                // sim_id is unused in single-db model, passing empty
                sim_id: String::new(),
                generation: Some(0),
            },
        };

        // Construct simulation config based on loaded config but overriding initialization
        let sim_config = Configuration {
            execution: config.execution.clone(),
            evolution: config.evolution.clone(), // Use loaded evolution config
            initialization: initialization_config,
        };

        // If seed override is present, we must update execution config inside sim_config
        let mut sim_config = sim_config;
        if let Some(seed) = seed_override {
            sim_config.execution.seed = Some(seed);
        }

        let mut sim = Simulation::new(sim_config)
            .map_err(|e| anyhow::anyhow!("Failed to initialize simulation: {e}"))?;

        rt.block_on(async {
            // Setup recorder
            let buffer_config = centrevo_sim::storage::BufferConfig {
                compression_level: 0,
                ..Default::default()
            };
            let recorder = Recorder::new(
                database,
                &config,
                buffer_config,
                config.execution.codec,
            )
            .context("Failed to create recorder")?;

            // Print all parameters
            println!("Configuration:");
            print_simulation_parameters(&sim);

            // Run simulation
            println!("Running {} generations...", config.execution.total_generations);

            let pb = if show_progress {
                let pb = ProgressBar::new(config.execution.total_generations as u64);
                pb.set_style(
                    ProgressStyle::default_bar()
                        .template(
                            "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta}) {per_sec}",
                        )
                        .unwrap()
                        .progress_chars("#>-"),
                );
                Some(pb)
            } else {
                None
            };

            for generation in 1..=config.execution.total_generations {
                // Determine if we need to record BEFORE stepping, or just check generation index
                // We check if we SHOULD record this generation
                let should_record =
                    RecordingStrategy::EveryN(record_every).should_record(generation);

                // Run step (CPU intensive, synchronous)
                sim.step()
                    .map_err(|e| anyhow::anyhow!("Generation {generation}: {e}"))?;

                // Record if needed (with RNG state for checkpoint)
                if should_record {
                    let rng_state = sim.rng_state_bytes();
                    recorder
                        .record_generation(sim.population(), generation, Some(rng_state), sim.arena())
                        .await
                        .context(format!("Failed to record generation {generation}"))?;
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
            // Finalize not needed

            
            recorder.close().await.context("Failed to close recorder")?;

            Ok::<(), anyhow::Error>(())
        })?;

        println!("\n✓ Simulation complete!");
        println!("  Final generation: {}", config.execution.total_generations);
    }

    println!("\n💡 Use 'centrevo info -d {}' to view results", database.display());

    Ok(())
}
