//! Example demonstrating async recording with compression.
//!
//! This example shows how to use AsyncRecorder for high-performance simulations
//! with automatic background compression and buffering.

use centrevo::{
    base::{Alphabet, Nucleotide},
    simulation::{
        SimulationConfig, RepeatStructure, MutationConfig, 
        RecombinationConfig, FitnessConfig,
    },
    storage::{AsyncRecorder, BufferConfig},
};

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== Async Recording Example ===\n");
    
    // Step 1: Determine buffer size based on simulation parameters
    let pop_size = 1000;
    let chr_length = 10_000;
    
    println!("Simulation parameters:");
    println!("  Population size: {}", pop_size);
    println!("  Chromosome length: {} bp", chr_length);
    
    // Use built-in configuration presets
    let buffer_config = if pop_size < 100 {
        BufferConfig::small()
    } else if pop_size < 500 {
        BufferConfig::medium()
    } else {
        BufferConfig::large()
    };
    
    let memory_usage = buffer_config.estimate_memory_usage(pop_size, chr_length);
    println!("\nBuffer configuration:");
    println!("  Capacity: {} snapshots", buffer_config.capacity);
    println!("  Compression level: {}", buffer_config.compression_level);
    println!("  Estimated memory: {:.1} MB", memory_usage as f64 / 1_048_576.0);
    println!("  Warning threshold: {:.0}%", buffer_config.warn_threshold * 100.0);
    
    // Step 2: Create async recorder
    let db_path = "/tmp/async_simulation.sqlite";
    let sim_id = "async_example";
    
    println!("\nCreating async recorder at: {}", db_path);
    let recorder = AsyncRecorder::new(db_path, sim_id, buffer_config)?;
    
    // Step 3: Create and run simulation
    let alphabet = Alphabet::dna();
    let structure = RepeatStructure::new(
        alphabet.clone(),
        Nucleotide::A,
        20,      // RU length
        50,      // RUs per HOR
        10,      // HORs per chromosome = 10kb total
        1,       // 1 chromosome per haplotype
    );
    
    let mutation = MutationConfig::uniform(alphabet, 0.001)?;
    let recombination = RecombinationConfig::standard(0.01, 0.7, 0.1)?;
    let fitness = FitnessConfig::neutral();
    let config = SimulationConfig::new(pop_size, 100, Some(42));
    
    let mut sim = centrevo::simulation::Simulation::new(
        structure,
        mutation,
        recombination,
        fitness,
        config,
    )?;
    
    println!("\nRunning simulation with async recording...");
    let start = std::time::Instant::now();
    
    // Record every 10th generation
    for generation in 0..100 {
        sim.step()?;
        
        if generation % 10 == 0 {
            // Record generation (non-blocking, unless buffer is full)
            recorder.record_generation(sim.population(), generation).await?;
            
            // Check buffer status
            let fill_ratio = recorder.buffer_fill_ratio();
            let status = if recorder.is_buffer_high() { "⚠️ HIGH" } else { "✓ OK" };
            
            println!(
                "  Gen {}: Buffer {:.0}% full {}",
                generation,
                fill_ratio * 100.0,
                status
            );
        }
    }
    
    let sim_time = start.elapsed();
    println!("\nSimulation completed in {:.2}s", sim_time.as_secs_f64());
    
    // Step 4: Close recorder and get statistics
    println!("\nFlushing buffer and closing recorder...");
    let stats = recorder.close().await?;
    
    let total_time = start.elapsed();
    
    // Step 5: Display results
    println!("\n=== Recording Statistics ===");
    println!("Generations recorded: {}", stats.generations_recorded);
    println!("Data compressed: {:.2} MB", stats.bytes_compressed as f64 / 1_048_576.0);
    println!("Data written: {:.2} MB", stats.bytes_written as f64 / 1_048_576.0);
    println!("Compression ratio: {:.1}%", stats.compression_ratio * 100.0);
    println!("Space saved: {:.2} MB ({:.1}%)", 
        (stats.bytes_compressed - stats.bytes_written) as f64 / 1_048_576.0,
        (1.0 - stats.compression_ratio) * 100.0
    );
    println!("\n=== Performance ===");
    println!("Avg compression time: {:.2} ms", stats.avg_compression_ms);
    println!("Avg write time: {:.2} ms", stats.avg_write_ms);
    println!("Buffer full events: {}", stats.buffer_full_count);
    println!("\n=== Total Time ===");
    println!("Simulation time: {:.2}s", sim_time.as_secs_f64());
    println!("Total time: {:.2}s", total_time.as_secs_f64());
    println!("Overhead: {:.1}%", 
        (total_time.as_secs_f64() - sim_time.as_secs_f64()) / sim_time.as_secs_f64() * 100.0
    );
    
    println!("\n✓ Example completed successfully!");
    println!("Database saved to: {}", db_path);
    
    Ok(())
}
