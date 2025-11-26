//! Asynchronous recorder with compression and buffering.
//!
//! This module provides an async recorder that compresses and writes data in the background,
//! allowing simulations to continue without waiting for I/O. The system uses a bounded channel
//! to buffer snapshots, only blocking the simulation when the buffer is full.

use crate::base::FitnessValue;
use crate::errors::DatabaseError;
use crate::simulation::Population;
use crate::storage::IndividualSnapshot;
use rusqlite::{Connection, params};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use tokio::sync::mpsc;
use tokio::task::JoinHandle;

/// Message sent from simulation to async recorder.
#[derive(Debug)]
enum RecorderMessage {
    /// Record a generation snapshot.
    Snapshot {
        generation: usize,
        snapshots: Vec<IndividualSnapshot>,
        mean_fitness: f64,
        min_fitness: f64,
        max_fitness: f64,
        std_fitness: f64,
        rng_state: Vec<u8>,
    },
    /// Flush and shutdown the recorder.
    Shutdown,
}

/// Statistics about the async recorder's performance.
#[derive(Debug, Clone, Default)]
pub struct RecorderStats {
    /// Total generations recorded.
    pub generations_recorded: usize,
    /// Total bytes compressed.
    pub bytes_compressed: usize,
    /// Total bytes written to database.
    pub bytes_written: usize,
    /// Compression ratio (compressed / original).
    pub compression_ratio: f64,
    /// Average compression time in milliseconds.
    pub avg_compression_ms: f64,
    /// Average write time in milliseconds.
    pub avg_write_ms: f64,
    /// Number of times the buffer was full (causing simulation to wait).
    pub buffer_full_count: usize,
}

/// Configuration for async recorder buffering.
#[derive(Debug, Clone)]
pub struct BufferConfig {
    /// Maximum number of snapshots to buffer before blocking.
    ///
    /// Rule of thumb for determining buffer size:
    /// - Small simulations (10-50 individuals): 5-10 snapshots
    /// - Medium simulations (100-500 individuals): 10-20 snapshots  
    /// - Large simulations (1000+ individuals): 20-50 snapshots
    ///
    /// Buffer memory usage ≈ buffer_capacity × pop_size × chromosome_size × 2
    /// Example: 10 snapshots × 100 individuals × 10KB × 2 = 20MB
    pub capacity: usize,

    /// Compression level (1-22, higher = better compression but slower).
    /// - Level 1-3: Fast compression (~200-300 MB/s)
    /// - Level 10-15: Balanced (default, ~50-100 MB/s)
    /// - Level 20-22: Maximum compression (~10-20 MB/s)
    pub compression_level: i32,

    /// Warn if buffer fill exceeds this percentage (0.0-1.0).
    pub warn_threshold: f64,
}

impl Default for BufferConfig {
    fn default() -> Self {
        Self {
            capacity: 10,          // 10 snapshots default
            compression_level: 10, // Balanced compression
            warn_threshold: 0.8,   // Warn at 80% full
        }
    }
}

impl BufferConfig {
    /// Create config optimized for small simulations (<100 individuals).
    pub fn small() -> Self {
        Self {
            capacity: 5,
            compression_level: 15, // Better compression for smaller data
            warn_threshold: 0.8,
        }
    }

    /// Create config optimized for medium simulations (100-500 individuals).
    pub fn medium() -> Self {
        Self {
            capacity: 10,
            compression_level: 10,
            warn_threshold: 0.8,
        }
    }

    /// Create config optimized for large simulations (>500 individuals).
    pub fn large() -> Self {
        Self {
            capacity: 20,
            compression_level: 5, // Faster compression for larger data
            warn_threshold: 0.85,
        }
    }

    /// Estimate memory usage in bytes for this buffer configuration.
    pub fn estimate_memory_usage(&self, pop_size: usize, chromosome_size: usize) -> usize {
        // Each snapshot has 2 haplotypes per individual
        self.capacity * pop_size * chromosome_size * 2
    }
}

/// Asynchronous recorder that compresses and writes in the background.
pub struct AsyncRecorder {
    /// Simulation ID.
    sim_id: Arc<str>,
    /// Sender for messages to background task.
    tx: mpsc::Sender<RecorderMessage>,
    /// Handle to background task.
    handle: Option<JoinHandle<Result<RecorderStats, DatabaseError>>>,
    /// Buffer configuration.
    config: BufferConfig,
    /// Current buffer fill level (approximate).
    buffer_fill: Arc<std::sync::atomic::AtomicUsize>,
}

impl AsyncRecorder {
    /// Create a new async recorder.
    ///
    /// The recorder spawns a background tokio task that handles compression and writing.
    /// The simulation can continue without waiting, unless the buffer fills up.
    pub fn new(
        db_path: impl AsRef<Path>,
        sim_id: impl Into<Arc<str>>,
        config: BufferConfig,
    ) -> Result<Self, DatabaseError> {
        let sim_id: Arc<str> = sim_id.into();
        let db_path = db_path.as_ref().to_path_buf();

        // Create bounded channel with configured capacity
        let (tx, rx) = mpsc::channel(config.capacity);

        // Shared buffer fill counter
        let buffer_fill = Arc::new(std::sync::atomic::AtomicUsize::new(0));
        let buffer_fill_clone = buffer_fill.clone();

        // Spawn background task
        let sim_id_clone = sim_id.clone();
        let compression_level = config.compression_level;
        let handle = tokio::spawn(async move {
            background_recorder_task(
                db_path,
                sim_id_clone,
                rx,
                compression_level,
                buffer_fill_clone,
            )
            .await
        });

        Ok(Self {
            sim_id,
            tx,
            handle: Some(handle),
            config,
            buffer_fill,
        })
    }

    /// Get simulation ID.
    pub fn sim_id(&self) -> &str {
        &self.sim_id
    }

    /// Get buffer configuration.
    pub fn config(&self) -> &BufferConfig {
        &self.config
    }

    /// Get current buffer fill level (0 to capacity).
    pub fn buffer_fill(&self) -> usize {
        self.buffer_fill.load(std::sync::atomic::Ordering::Relaxed)
    }

    /// Get buffer fill percentage (0.0 to 1.0).
    pub fn buffer_fill_ratio(&self) -> f64 {
        self.buffer_fill() as f64 / self.config.capacity as f64
    }

    /// Check if buffer is above warning threshold.
    pub fn is_buffer_high(&self) -> bool {
        self.buffer_fill_ratio() >= self.config.warn_threshold
    }

    /// Record a generation (non-blocking, unless buffer is full).
    ///
    /// This function will:
    /// 1. Create snapshots of the population (in caller's thread)
    /// 2. Send to background task (blocks only if buffer is full)
    /// 3. Return immediately (background task handles compression/write)
    ///
    /// # Parameters
    /// - `population`: The current population to record
    /// - `generation`: The generation number
    /// - `rng_state`: The RNG state bytes for checkpointing
    ///
    /// # Performance
    /// - Snapshot creation: ~1-5ms (parallel)
    /// - Send to channel: <1μs (unless buffer full)
    /// - Total blocking time: ~1-5ms (or more if buffer full)
    pub async fn record_generation(
        &self,
        population: &Population,
        generation: usize,
        rng_state: Vec<u8>,
    ) -> Result<(), DatabaseError> {
        // Create snapshots in parallel (happens in caller's context)
        use rayon::prelude::*;
        let snapshots: Vec<IndividualSnapshot> = population
            .individuals()
            .par_iter()
            .map(IndividualSnapshot::from_individual)
            .collect();

        // Calculate fitness stats
        let fitnesses: Vec<f64> = population
            .individuals()
            .iter()
            .filter_map(|ind| ind.cached_fitness())
            .map(|f| f.get())
            .collect();

        let (mean_fitness, min_fitness, max_fitness, std_fitness) = if fitnesses.is_empty() {
            (0.0, 0.0, 0.0, 0.0)
        } else {
            let sum: f64 = fitnesses.iter().sum();
            let mean = sum / fitnesses.len() as f64;
            let min = fitnesses.iter().copied().fold(f64::INFINITY, f64::min);
            let max = fitnesses.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let variance: f64 =
                fitnesses.iter().map(|&f| (f - mean).powi(2)).sum::<f64>() / fitnesses.len() as f64;
            let std = variance.sqrt();
            (mean, min, max, std)
        };

        // Check if buffer is getting full
        if self.is_buffer_high() {
            eprintln!(
                "Warning: Recorder buffer at {:.1}% capacity (generation {}). \
                 Simulation may slow down if compression can't keep up.",
                self.buffer_fill_ratio() * 100.0,
                generation
            );
        }

        // Increment buffer fill counter before sending
        self.buffer_fill
            .fetch_add(1, std::sync::atomic::Ordering::Relaxed);

        // Send to background task (blocks only if buffer is full)
        let result = self
            .tx
            .send(RecorderMessage::Snapshot {
                generation,
                snapshots,
                mean_fitness,
                min_fitness,
                max_fitness,
                std_fitness,
                rng_state,
            })
            .await;

        // If send failed, decrement the counter
        if result.is_err() {
            self.buffer_fill
                .fetch_sub(1, std::sync::atomic::Ordering::Relaxed);
            return Err(DatabaseError::Insert("Recorder task died".to_string()));
        }

        Ok(())
    }

    /// Flush all pending writes and shut down the recorder.
    ///
    /// This method waits for all buffered snapshots to be compressed and written,
    /// then returns statistics about the recording session.
    pub async fn close(mut self) -> Result<RecorderStats, DatabaseError> {
        // Send shutdown message
        if let Err(e) = self.tx.send(RecorderMessage::Shutdown).await {
            eprintln!("Warning: Failed to send shutdown message: {e}");
        }

        // Wait for background task to finish
        if let Some(handle) = self.handle.take() {
            handle
                .await
                .map_err(|e| DatabaseError::Close(format!("Background task panicked: {e}")))?
        } else {
            Ok(RecorderStats::default())
        }
    }
}

/// Background task that handles compression and database writes.
async fn background_recorder_task(
    db_path: PathBuf,
    sim_id: Arc<str>,
    mut rx: mpsc::Receiver<RecorderMessage>,
    compression_level: i32,
    buffer_fill: Arc<std::sync::atomic::AtomicUsize>,
) -> Result<RecorderStats, DatabaseError> {
    // Open database connection (must be done in background thread)
    let mut conn =
        Connection::open(&db_path).map_err(|e| DatabaseError::Connection(e.to_string()))?;

    // Performance pragmas
    conn.execute_batch(
        "PRAGMA synchronous = NORMAL;
         PRAGMA journal_mode = WAL;
         PRAGMA temp_store = MEMORY;
         PRAGMA cache_size = -64000;",
    )
    .map_err(|e| DatabaseError::Initialization(e.to_string()))?;

    // Initialize database schema
    conn.execute_batch(
        "-- Main population state table
        CREATE TABLE IF NOT EXISTS population_state (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            sim_id TEXT NOT NULL,
            generation INTEGER NOT NULL,
            individual_id TEXT NOT NULL,
            haplotype1_chr_id TEXT NOT NULL,
            haplotype1_seq BLOB NOT NULL,
            haplotype2_chr_id TEXT NOT NULL,
            haplotype2_seq BLOB NOT NULL,
            fitness REAL,
            timestamp INTEGER DEFAULT (strftime('%s', 'now'))
        );

        -- Simulation metadata table (expanded for full config)
        CREATE TABLE IF NOT EXISTS simulations (
            sim_id TEXT PRIMARY KEY,
            start_time INTEGER NOT NULL,
            end_time INTEGER,
            pop_size INTEGER NOT NULL,
            num_generations INTEGER NOT NULL,
            mutation_rate REAL NOT NULL,
            recombination_rate REAL NOT NULL,
            parameters_json TEXT NOT NULL,
            config_json TEXT  -- Complete simulation configuration
        );

        -- Fitness history table for aggregated statistics
        CREATE TABLE IF NOT EXISTS fitness_history (
            sim_id TEXT NOT NULL,
            generation INTEGER NOT NULL,
            mean_fitness REAL NOT NULL,
            min_fitness REAL NOT NULL,
            max_fitness REAL NOT NULL,
            std_fitness REAL NOT NULL,
            PRIMARY KEY (sim_id, generation)
        );

        -- Checkpoints table for resumability
        CREATE TABLE IF NOT EXISTS checkpoints (
            sim_id TEXT NOT NULL,
            generation INTEGER NOT NULL,
            rng_state BLOB NOT NULL,
            timestamp INTEGER NOT NULL,
            PRIMARY KEY (sim_id, generation)
        );

        -- Indices for fast queries
        CREATE INDEX IF NOT EXISTS idx_pop_sim_gen 
            ON population_state(sim_id, generation);
        CREATE INDEX IF NOT EXISTS idx_pop_individual 
            ON population_state(individual_id);
        CREATE INDEX IF NOT EXISTS idx_fitness_sim_gen 
            ON fitness_history(sim_id, generation);
        CREATE INDEX IF NOT EXISTS idx_checkpoints_sim_gen
            ON checkpoints(sim_id, generation);",
    )
    .map_err(|e| DatabaseError::Initialization(e.to_string()))?;

    let mut stats = RecorderStats::default();
    let mut total_compression_time = std::time::Duration::ZERO;
    let mut total_write_time = std::time::Duration::ZERO;

    // Process messages until shutdown
    while let Some(msg) = rx.recv().await {
        match msg {
            RecorderMessage::Snapshot {
                generation,
                snapshots,
                mean_fitness,
                min_fitness,
                max_fitness,
                std_fitness,
                rng_state,
            } => {
                // Decrement buffer fill counter
                buffer_fill.fetch_sub(1, std::sync::atomic::Ordering::Relaxed);

                // Compress and write
                let compress_start = std::time::Instant::now();
                let compressed = compress_snapshots(&snapshots, compression_level)?;
                let compress_time = compress_start.elapsed();
                total_compression_time += compress_time;

                let write_start = std::time::Instant::now();
                write_compressed_snapshots(
                    &mut conn,
                    &sim_id,
                    generation,
                    &compressed,
                    mean_fitness,
                    min_fitness,
                    max_fitness,
                    std_fitness,
                    &rng_state,
                )?;
                let write_time = write_start.elapsed();
                total_write_time += write_time;

                // Update stats
                stats.generations_recorded += 1;
                let original_size: usize = snapshots
                    .iter()
                    .map(|s| s.haplotype1_seq.len() + s.haplotype2_seq.len())
                    .sum();
                let compressed_size: usize = compressed
                    .iter()
                    .map(|c| c.haplotype1_compressed.len() + c.haplotype2_compressed.len())
                    .sum();

                stats.bytes_compressed += original_size;
                stats.bytes_written += compressed_size;
            }
            RecorderMessage::Shutdown => {
                break;
            }
        }
    }

    // Calculate final statistics
    if stats.generations_recorded > 0 {
        stats.compression_ratio = stats.bytes_written as f64 / stats.bytes_compressed as f64;
        stats.avg_compression_ms =
            total_compression_time.as_secs_f64() * 1000.0 / stats.generations_recorded as f64;
        stats.avg_write_ms =
            total_write_time.as_secs_f64() * 1000.0 / stats.generations_recorded as f64;
    }

    // Close connection
    conn.close()
        .map_err(|(_, e)| DatabaseError::Close(e.to_string()))?;

    Ok(stats)
}

/// Compressed snapshot data.
struct CompressedSnapshot {
    individual_id: String,
    haplotype1_chr_id: String,
    haplotype1_compressed: Vec<u8>,
    haplotype2_chr_id: String,
    haplotype2_compressed: Vec<u8>,
    fitness: Option<f64>,
}

/// Compress snapshots using zstd.
fn compress_snapshots(
    snapshots: &[IndividualSnapshot],
    level: i32,
) -> Result<Vec<CompressedSnapshot>, DatabaseError> {
    use rayon::prelude::*;

    let compressed = snapshots
        .par_iter()
        .map(|snapshot| {
            let h1_compressed = zstd::bulk::compress(&snapshot.haplotype1_seq, level)
                .map_err(|e| DatabaseError::Insert(format!("Compression failed: {e}")))?;
            let h2_compressed = zstd::bulk::compress(&snapshot.haplotype2_seq, level)
                .map_err(|e| DatabaseError::Insert(format!("Compression failed: {e}")))?;

            Ok(CompressedSnapshot {
                individual_id: snapshot.individual_id.clone(),
                haplotype1_chr_id: snapshot.haplotype1_chr_id.clone(),
                haplotype1_compressed: h1_compressed,
                haplotype2_chr_id: snapshot.haplotype2_chr_id.clone(),
                haplotype2_compressed: h2_compressed,
                fitness: snapshot.fitness,
            })
        })
        .collect::<Result<Vec<_>, _>>()?;

    Ok(compressed)
}

/// Write compressed snapshots to database.
#[allow(clippy::too_many_arguments)]
fn write_compressed_snapshots(
    conn: &mut Connection,
    sim_id: &str,
    generation: usize,
    snapshots: &[CompressedSnapshot],
    mean_fitness: f64,
    min_fitness: f64,
    max_fitness: f64,
    std_fitness: f64,
    rng_state: &[u8],
) -> Result<(), DatabaseError> {
    let tx = conn
        .transaction()
        .map_err(|e| DatabaseError::Transaction(e.to_string()))?;

    // Insert population state
    {
        let mut stmt = tx
            .prepare_cached(
                "INSERT INTO population_state 
                (sim_id, generation, individual_id, haplotype1_chr_id, haplotype1_seq, 
                 haplotype2_chr_id, haplotype2_seq, fitness)
                VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8)",
            )
            .map_err(|e| DatabaseError::Insert(e.to_string()))?;

            for snapshot in snapshots {
            stmt.execute(params![
                sim_id,
                generation as i64,
                snapshot.individual_id,
                snapshot.haplotype1_chr_id,
                &snapshot.haplotype1_compressed,
                snapshot.haplotype2_chr_id,
                &snapshot.haplotype2_compressed,
                    snapshot.fitness,
            ])
            .map_err(|e| DatabaseError::Insert(e.to_string()))?;
        }
    }

    // Insert fitness history
    tx.execute(
        "INSERT OR REPLACE INTO fitness_history 
        (sim_id, generation, mean_fitness, min_fitness, max_fitness, std_fitness)
        VALUES (?1, ?2, ?3, ?4, ?5, ?6)",
        params![
            sim_id,
            generation as i64,
            mean_fitness,
            min_fitness,
            max_fitness,
            std_fitness
        ],
    )
    .map_err(|e| DatabaseError::Insert(e.to_string()))?;

    // Insert checkpoint with RNG state
    let timestamp = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_secs() as i64;

    tx.execute(
        "INSERT OR REPLACE INTO checkpoints 
        (sim_id, generation, rng_state, timestamp)
        VALUES (?1, ?2, ?3, ?4)",
        params![sim_id, generation as i64, rng_state, timestamp],
    )
    .map_err(|e| DatabaseError::Insert(e.to_string()))?;

    tx.commit()
        .map_err(|e| DatabaseError::Transaction(e.to_string()))?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::base::Nucleotide;
    use crate::genome::{Chromosome, Haplotype, Individual};
    use crate::simulation::Population;

    fn create_test_individual(id: &str, length: usize) -> Individual {
        // Convert total length to num_hors: ru_length=20, rus_per_hor=5, so one HOR = 100 bp
        let num_hors = length / 100;
        let chr1 = Chromosome::uniform(format!("{}_h1_chr1", id), Nucleotide::A, 20, 5, num_hors);
        let chr2 = Chromosome::uniform(format!("{}_h2_chr1", id), Nucleotide::C, 20, 5, num_hors);

        let h1 = Haplotype::from_chromosomes(vec![chr1]);
        let h2 = Haplotype::from_chromosomes(vec![chr2]);

        Individual::new(id, h1, h2)
    }

    fn create_test_population(size: usize, chr_length: usize) -> Population {
        let mut individuals = Vec::new();
        for i in 0..size {
            let mut ind = create_test_individual(&format!("ind_{}", i), chr_length);
            ind.set_cached_fitness(FitnessValue::new(i as f64 / size as f64));
            individuals.push(ind);
        }
        Population::new("test_pop", individuals)
    }

    fn dummy_rng_state() -> Vec<u8> {
        vec![0u8; 32] // Dummy RNG state for testing
    }

    #[tokio::test]
    async fn test_buffer_config() {
        let config = BufferConfig::default();
        assert_eq!(config.capacity, 10);

        let small = BufferConfig::small();
        assert_eq!(small.capacity, 5);

        let medium = BufferConfig::medium();
        assert_eq!(medium.capacity, 10);

        let large = BufferConfig::large();
        assert_eq!(large.capacity, 20);
    }

    #[tokio::test]
    async fn test_buffer_memory_estimation() {
        let config = BufferConfig::default();
        let memory = config.estimate_memory_usage(100, 10_000);
        // 10 snapshots × 100 individuals × 10KB × 2 haplotypes = 20MB
        assert_eq!(memory, 20_000_000);
    }

    #[tokio::test]
    async fn test_async_recorder_creation() {
        let path = "/tmp/test_async_recorder.sqlite";
        let _ = std::fs::remove_file(path);

        let config = BufferConfig::small();
        let recorder =
            AsyncRecorder::new(path, "test_sim", config).expect("Failed to create async recorder");

        assert_eq!(recorder.sim_id(), "test_sim");
        assert_eq!(recorder.buffer_fill(), 0);

        let stats = recorder.close().await.expect("Failed to close recorder");
        assert_eq!(stats.generations_recorded, 0);

        std::fs::remove_file(path).ok();
    }

    #[tokio::test]
    async fn test_record_generation() {
        let path = "/tmp/test_async_record.sqlite";
        let _ = std::fs::remove_file(path);

        let config = BufferConfig::small();
        let recorder =
            AsyncRecorder::new(path, "test_sim", config).expect("Failed to create async recorder");

        let pop = create_test_population(5, 100);

        recorder
            .record_generation(&pop, 0, dummy_rng_state())
            .await
            .expect("Failed to record generation");

        let stats = recorder.close().await.expect("Failed to close recorder");
        assert_eq!(stats.generations_recorded, 1);
        assert!(stats.compression_ratio < 1.0); // Should achieve compression

        std::fs::remove_file(path).ok();
    }

    #[tokio::test]
    async fn test_multiple_generations() {
        let path = "/tmp/test_async_multiple.sqlite";
        let _ = std::fs::remove_file(path);

        let config = BufferConfig::medium();
        let recorder =
            AsyncRecorder::new(path, "test_sim", config).expect("Failed to create async recorder");

        let pop = create_test_population(10, 200);

        // Record 5 generations
        for generation in 0..5 {
            recorder
                .record_generation(&pop, generation, dummy_rng_state())
                .await
                .expect("Failed to record generation");
        }

        let stats = recorder.close().await.expect("Failed to close recorder");
        assert_eq!(stats.generations_recorded, 5);
        assert!(stats.bytes_written < stats.bytes_compressed); // Compression working

        println!("Compression ratio: {:.2}%", stats.compression_ratio * 100.0);
        println!("Avg compression time: {:.2}ms", stats.avg_compression_ms);
        println!("Avg write time: {:.2}ms", stats.avg_write_ms);

        std::fs::remove_file(path).ok();
    }

    #[tokio::test]
    async fn test_buffer_fill_tracking() {
        let path = "/tmp/test_buffer_fill.sqlite";
        let _ = std::fs::remove_file(path);

        let config = BufferConfig::small(); // capacity = 5
        let recorder =
            AsyncRecorder::new(path, "test_sim", config).expect("Failed to create async recorder");

        let pop = create_test_population(5, 50);

        // Initially buffer should be empty
        assert_eq!(recorder.buffer_fill(), 0);
        assert_eq!(recorder.buffer_fill_ratio(), 0.0);

        // Record a generation
        recorder
            .record_generation(&pop, 0, dummy_rng_state())
            .await
            .expect("Failed to record");

        // Buffer should have 1 item (or be processed already)
        let fill = recorder.buffer_fill();
        assert!(fill <= 1, "Buffer fill should be 0 or 1, got {}", fill);

        let stats = recorder.close().await.expect("Failed to close");
        assert_eq!(stats.generations_recorded, 1);

        std::fs::remove_file(path).ok();
    }

    #[tokio::test]
    async fn test_compression_levels() {
        let path = "/tmp/test_compression_levels.sqlite";

        for level in [1, 5, 10, 15, 20] {
            let _ = std::fs::remove_file(path);

            let config = BufferConfig {
                capacity: 5,
                compression_level: level,
                warn_threshold: 0.8,
            };

            let recorder =
                AsyncRecorder::new(path, "test_sim", config).expect("Failed to create recorder");

            let pop = create_test_population(10, 100);

            recorder
                .record_generation(&pop, 0, dummy_rng_state())
                .await
                .expect("Failed to record");

            let stats = recorder.close().await.expect("Failed to close");

            println!(
                "Level {}: ratio={:.2}%, time={:.2}ms",
                level,
                stats.compression_ratio * 100.0,
                stats.avg_compression_ms
            );

            assert_eq!(stats.generations_recorded, 1);
            assert!(stats.compression_ratio < 1.0);
        }

        std::fs::remove_file(path).ok();
    }

    #[tokio::test]
    async fn test_large_population() {
        let path = "/tmp/test_large_pop.sqlite";
        let _ = std::fs::remove_file(path);

        let config = BufferConfig::large();
        let recorder =
            AsyncRecorder::new(path, "test_sim", config).expect("Failed to create recorder");

        // Test with larger population
        let pop = create_test_population(100, 500);

        recorder
            .record_generation(&pop, 0, dummy_rng_state())
            .await
            .expect("Failed to record");

        let stats = recorder.close().await.expect("Failed to close");
        assert_eq!(stats.generations_recorded, 1);
        assert!(stats.bytes_compressed > 0);
        assert!(stats.bytes_written > 0);
        assert!(stats.bytes_written < stats.bytes_compressed);

        std::fs::remove_file(path).ok();
    }
}
