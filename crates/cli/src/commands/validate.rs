use anyhow::{Context, Result};
use centrevo_sim::storage::QueryBuilder;
use std::path::PathBuf;

pub fn validate_database(database: &PathBuf, fix: bool) -> Result<()> {
    println!("🔍 Validating database: {}", database.display());

    if !database.exists() {
        anyhow::bail!("Database file does not exist");
    }

    let query = QueryBuilder::new(database).context("Failed to open database")?;

    // Single simulation per database
    let sim_name = database
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("simulation");

    println!("\nValidating simulation: {sim_name}");
    println!("{}", "-".repeat(50));

    // Check metadata and config
    let config = match query.get_full_config() {
        Ok(cfg) => {
            println!("✓ Metadata: OK");
            cfg
        }
        Err(e) => {
            println!("✗ Metadata: FAILED - {e}");
            return Ok(());
        }
    };

    let info_num_generations = config.execution.total_generations;
    let info_pop_size = config.execution.population_size;

    // Check recorded generations
    let recorded_gens = match query.get_recorded_generations() {
        Ok(gens) => {
            println!("✓ Recorded generations: {} snapshots", gens.len());
            gens
        }
        Err(e) => {
            println!("✗ Failed to query generations: {e}");
            return Ok(());
        }
    };

    let mut total_issues = 0;

    if recorded_gens.is_empty() {
        println!("⚠️  No generation data recorded");
        total_issues += 1;
    }

    // Check for gaps in generations
    let expected_gens: Vec<usize> = (0..=info_num_generations).step_by(100).collect();
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
        match query.get_generation(generation_num) {
            Ok(snapshots) => {
                if snapshots.is_empty() {
                    println!("⚠️  Generation {generation_num} has no individuals");
                    total_issues += 1;
                } else if snapshots.len() != info_pop_size {
                    println!(
                        "⚠️  Generation {} has {} individuals (expected {})",
                        generation_num,
                        snapshots.len(),
                        info_pop_size
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
    match query.get_fitness_history() {
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
