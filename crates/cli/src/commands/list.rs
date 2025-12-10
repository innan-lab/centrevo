use anyhow::{Context, Result};
use centrevo_sim::storage::QueryBuilder;
use std::path::PathBuf;

pub fn list_simulations(database: &PathBuf) -> Result<()> {
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

pub fn show_info(database: &PathBuf, name: &str) -> Result<()> {
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

pub fn show_generations(database: &PathBuf, name: &str) -> Result<()> {
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
