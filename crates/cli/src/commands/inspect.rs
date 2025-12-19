use anyhow::{Context, Result};
use centrevo_sim::storage::QueryBuilder;
use std::path::PathBuf;

pub fn show_info(database: &PathBuf) -> Result<()> {
    let query = QueryBuilder::new(database).context("Failed to open database")?;
    let config = query
        .get_full_config()
        .context("Failed to get simulation info")?;

    println!("\n📊 Simulation Information");
    println!("{}", "=".repeat(50));
    println!("Population size: {}", config.execution.population_size);
    println!("Generations: {}", config.execution.total_generations);

    let params_json = serde_json::to_string_pretty(&config).unwrap_or_default();
    println!("\nParameters:");
    println!("{params_json}");

    Ok(())
}

pub fn show_generations(database: &PathBuf) -> Result<()> {
    let query = QueryBuilder::new(database).context("Failed to open database")?;
    let generations = query
        .get_recorded_generations()
        .context("Failed to get generations")?;

    if generations.is_empty() {
        println!("No recorded generations found.");
        return Ok(());
    }

    println!("\n📈 Recorded Generations:");
    println!("{}", "=".repeat(50));
    println!("Generations: {generations:?}");
    println!("Total: {} snapshots", generations.len());

    Ok(())
}
