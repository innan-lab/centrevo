//! Asynchronous recorder with compression and buffering.

use crate::errors::DatabaseError;
use crate::simulation::{Configuration, Population};
use crate::storage::types::IndividualSnapshot;
use centrevo_codec::CodecStrategy;
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
        rng_state: Option<Vec<u8>>, // Checkpoint if present
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
    pub capacity: usize,
    /// Compression level (1-22).
    pub compression_level: i32,
    /// Warn if buffer fill exceeds this percentage (0.0-1.0).
    pub warn_threshold: f64,
}

impl Default for BufferConfig {
    fn default() -> Self {
        Self {
            capacity: 10,
            compression_level: 10,
            warn_threshold: 0.8,
        }
    }
}

impl BufferConfig {
    pub fn small() -> Self {
        Self {
            capacity: 5,
            compression_level: 15,
            warn_threshold: 0.8,
        }
    }
}

/// Asynchronous recorder that compresses and writes in the background.
pub struct AsyncRecorder {
    tx: mpsc::Sender<RecorderMessage>,
    handle: Option<JoinHandle<Result<RecorderStats, DatabaseError>>>,
    config: BufferConfig,
    buffer_fill: Arc<std::sync::atomic::AtomicUsize>,
    codec: CodecStrategy,
}

impl AsyncRecorder {
    /// Create a new async recorder.
    pub fn new(
        db_path: impl AsRef<Path>,
        sim_config: &Configuration,
        buffer_config: BufferConfig,
        codec: CodecStrategy,
    ) -> Result<Self, DatabaseError> {
        let db_path = db_path.as_ref().to_path_buf();

        let (tx, rx) = mpsc::channel(buffer_config.capacity);
        let buffer_fill = Arc::new(std::sync::atomic::AtomicUsize::new(0));
        let buffer_fill_clone = buffer_fill.clone();
        let compression_level = buffer_config.compression_level;

        // Initialize Metadata (Blocking write for safety at start)
        {
            let mut db = crate::storage::database::Database::open(&db_path)?;

            // Write metadata
            let exec = &sim_config.execution;
            let meta_pairs = vec![
                ("population_size", exec.population_size.to_string()),
                ("total_generations", exec.total_generations.to_string()),
                (
                    "mutation_rate",
                    format!("{:?}", sim_config.evolution.mutation),
                ),
                (
                    "seed",
                    exec.seed
                        .map(|s| s.to_string())
                        .unwrap_or_else(|| "None".to_string()),
                ),
                (
                    "full_config_json",
                    serde_json::to_string(sim_config).unwrap_or_default(),
                ),
                ("codec", format!("{:?}", exec.codec)),
                (
                    "created_at",
                    std::time::SystemTime::now()
                        .duration_since(std::time::UNIX_EPOCH)
                        .unwrap()
                        .as_secs()
                        .to_string(),
                ),
            ];

            let tx = db.transaction()?;
            {
                let mut stmt = tx
                    .prepare("INSERT OR REPLACE INTO metadata (key, value) VALUES (?1, ?2)")
                    .map_err(|e| DatabaseError::Insert(e.to_string()))?;
                for (k, v) in meta_pairs {
                    stmt.execute(params![k, v])
                        .map_err(|e| DatabaseError::Insert(e.to_string()))?;
                }
            }
            tx.commit()
                .map_err(|e| DatabaseError::Transaction(e.to_string()))?;
            db.close()?;
        }

        let handle = tokio::spawn(async move {
            background_recorder_task(db_path, rx, compression_level, buffer_fill_clone).await
        });

        Ok(Self {
            tx,
            handle: Some(handle),
            config: buffer_config,
            buffer_fill,
            codec,
        })
    }

    pub fn buffer_fill(&self) -> usize {
        self.buffer_fill.load(std::sync::atomic::Ordering::Relaxed)
    }

    pub fn is_buffer_high(&self) -> bool {
        (self.buffer_fill() as f64 / self.config.capacity as f64) >= self.config.warn_threshold
    }

    pub async fn record_generation(
        &self,
        population: &Population,
        generation: usize,
        rng_state: Option<Vec<u8>>,
    ) -> Result<(), DatabaseError> {
        use rayon::prelude::*;
        let codec = &self.codec;

        let snapshots: Vec<IndividualSnapshot> = population
            .individuals()
            .par_iter()
            .map(|ind| IndividualSnapshot::from_individual(ind, codec))
            .collect();

        if self.is_buffer_high() {
            eprintln!(
                "Warning: Recorder buffer high at generation {}.",
                generation
            );
        }

        self.buffer_fill
            .fetch_add(1, std::sync::atomic::Ordering::Relaxed);

        let result = self
            .tx
            .send(RecorderMessage::Snapshot {
                generation,
                snapshots,
                rng_state,
            })
            .await;

        if result.is_err() {
            self.buffer_fill
                .fetch_sub(1, std::sync::atomic::Ordering::Relaxed);
            return Err(DatabaseError::Insert("Recorder task died".to_string()));
        }

        Ok(())
    }

    pub async fn close(mut self) -> Result<RecorderStats, DatabaseError> {
        if let Err(e) = self.tx.send(RecorderMessage::Shutdown).await {
            eprintln!("Warning: Failed to send shutdown message: {e}");
        }
        if let Some(handle) = self.handle.take() {
            handle
                .await
                .map_err(|e| DatabaseError::Close(format!("Background task panicked: {e}")))?
        } else {
            Ok(RecorderStats::default())
        }
    }
}

async fn background_recorder_task(
    db_path: PathBuf,
    mut rx: mpsc::Receiver<RecorderMessage>,
    compression_level: i32,
    buffer_fill: Arc<std::sync::atomic::AtomicUsize>,
) -> Result<RecorderStats, DatabaseError> {
    // Open DB (Assuming already initialized in new())
    let mut db = crate::storage::database::Database::open(&db_path)?;
    let mut stats = RecorderStats::default();

    let mut total_compression_time = std::time::Duration::ZERO;
    let mut total_write_time = std::time::Duration::ZERO;

    while let Some(msg) = rx.recv().await {
        match msg {
            RecorderMessage::Snapshot {
                generation,
                snapshots,
                rng_state,
            } => {
                buffer_fill.fetch_sub(1, std::sync::atomic::Ordering::Relaxed);

                let compress_start = std::time::Instant::now();
                let compressed = compress_snapshots(&snapshots, compression_level)?;
                total_compression_time += compress_start.elapsed();

                let write_start = std::time::Instant::now();
                write_compressed_snapshots(
                    db.connection_mut(),
                    generation,
                    &compressed,
                    &snapshots, // Need original for fitness values? Or we can extract fitness separately or put it in compressed struct.
                    // Actually, CompressedSnapshot should carry fitness.
                    rng_state.as_deref(),
                )?;
                total_write_time += write_start.elapsed();

                stats.generations_recorded += 1;
                // Rough stat calc
                let original_size: usize = snapshots.len() * 200; // Placeholder average
                let compressed_size: usize = compressed.len() * 50; // Placeholder average
                stats.bytes_compressed += original_size;
                stats.bytes_written += compressed_size;
            }
            RecorderMessage::Shutdown => {
                break;
            }
        }
    }

    if stats.generations_recorded > 0 {
        stats.compression_ratio = stats.bytes_written as f64 / stats.bytes_compressed.max(1) as f64;
        stats.avg_compression_ms =
            total_compression_time.as_secs_f64() * 1000.0 / stats.generations_recorded as f64;
        stats.avg_write_ms =
            total_write_time.as_secs_f64() * 1000.0 / stats.generations_recorded as f64;
    }

    db.close()?;
    Ok(stats)
}

struct CompressedSnapshot {
    haplotype1_compressed: Vec<u8>,
    haplotype1_map_compressed: Option<Vec<u8>>,
    haplotype2_compressed: Vec<u8>,
    haplotype2_map_compressed: Option<Vec<u8>>,
    fitness: Option<f64>,
}

fn compress_snapshots(
    snapshots: &[IndividualSnapshot],
    level: i32,
) -> Result<Vec<CompressedSnapshot>, DatabaseError> {
    use rayon::prelude::*;
    snapshots
        .par_iter()
        .map(|snapshot| {
            let (h1_compressed, h2_compressed, h1_map, h2_map) = if level == 0 {
                (
                    snapshot.haplotype1_seq.clone(),
                    snapshot.haplotype2_seq.clone(),
                    snapshot.haplotype1_map.clone(),
                    snapshot.haplotype2_map.clone(),
                )
            } else {
                let h1 = zstd::bulk::compress(&snapshot.haplotype1_seq, level)
                    .map_err(|e| DatabaseError::Insert(format!("Compression failed: {e}")))?;
                let h2 = zstd::bulk::compress(&snapshot.haplotype2_seq, level)
                    .map_err(|e| DatabaseError::Insert(format!("Compression failed: {e}")))?;

                let m1 = snapshot
                    .haplotype1_map
                    .as_ref()
                    .map(|m| zstd::bulk::compress(m, level))
                    .transpose()
                    .map_err(|e| DatabaseError::Insert(format!("Map fail: {e}")))?;
                let m2 = snapshot
                    .haplotype2_map
                    .as_ref()
                    .map(|m| zstd::bulk::compress(m, level))
                    .transpose()
                    .map_err(|e| DatabaseError::Insert(format!("Map fail: {e}")))?;
                (h1, h2, m1, m2)
            };

            Ok(CompressedSnapshot {
                haplotype1_compressed: h1_compressed,
                haplotype1_map_compressed: h1_map,
                haplotype2_compressed: h2_compressed,
                haplotype2_map_compressed: h2_map,
                fitness: snapshot.fitness,
            })
        })
        .collect()
}

fn write_compressed_snapshots(
    conn: &mut Connection,
    generation: usize,
    compressed: &[CompressedSnapshot],
    _original: &[IndividualSnapshot], // Unused but signature match from caller
    rng_state: Option<&[u8]>,
) -> Result<(), DatabaseError> {
    let tx = conn
        .transaction()
        .map_err(|e| DatabaseError::Transaction(e.to_string()))?;
    {
        // 1. Write State (Normalized: 1 row per chromosome)
        let mut stmt_state = tx
            .prepare_cached(
                "INSERT INTO state
            (generation, individual_id, haplotype_id, chromosome_id, seq, map)
            VALUES (?1, ?2, ?3, ?4, ?5, ?6)",
            )
            .map_err(|e| DatabaseError::Insert(e.to_string()))?;

        // 2. Write Individual Fitness
        let mut stmt_fit = tx
            .prepare_cached(
                "INSERT INTO individual_fitness
            (generation, individual_id, fitness, timestamp)
            VALUES (?1, ?2, ?3, ?4)",
            )
            .map_err(|e| DatabaseError::Insert(e.to_string()))?;

        let timestamp = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_secs() as i64;

        for (i, snap) in compressed.iter().enumerate() {
            let individual_id = i + 1; // 1-based ID

            // Haplotype 1 / Chromosome 1
            stmt_state
                .execute(params![
                    generation,
                    individual_id,
                    1, // Hap 1
                    1, // Chr 1 (Always 1 for now)
                    &snap.haplotype1_compressed,
                    snap.haplotype1_map_compressed
                ])
                .map_err(|e| DatabaseError::Insert(e.to_string()))?;

            // Haplotype 2 / Chromosome 1
            stmt_state
                .execute(params![
                    generation,
                    individual_id,
                    2, // Hap 2
                    1, // Chr 1
                    &snap.haplotype2_compressed,
                    snap.haplotype2_map_compressed
                ])
                .map_err(|e| DatabaseError::Insert(e.to_string()))?;

            // Fitness
            if let Some(f) = snap.fitness {
                stmt_fit
                    .execute(params![generation, individual_id, f, timestamp])
                    .map_err(|e| DatabaseError::Insert(e.to_string()))?;
            }
        }
    }

    // 3. Write Checkpoint if RNG state provided
    if let Some(rng) = rng_state {
        let timestamp = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_secs() as i64;

        tx.execute(
            "INSERT OR REPLACE INTO checkpoints (generation, rng_state, timestamp) VALUES (?1, ?2, ?3)",
            params![generation, rng, timestamp]
        ).map_err(|e| DatabaseError::Insert(e.to_string()))?;
    }

    tx.commit()
        .map_err(|e| DatabaseError::Transaction(e.to_string()))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::base::{FitnessValue, Nucleotide};
    use crate::genome::{Chromosome, Haplotype, Individual};
    use crate::simulation::{
        EvolutionConfig, ExecutionConfig, FitnessConfig, GenerationMode, InitializationConfig,
        MutationConfig, Population, RecombinationConfig, UniformRepeatStructure,
    };

    // Mock Helpers
    fn mock_config() -> Configuration {
        Configuration {
            execution: ExecutionConfig::new(10, 10, Some(1)),
            evolution: EvolutionConfig {
                mutation: MutationConfig::uniform(0.01).unwrap(),
                recombination: RecombinationConfig::standard(0.01, 0.5, 0.1).unwrap(),
                fitness: FitnessConfig::neutral(),
            },
            initialization: InitializationConfig::Generate {
                structure: UniformRepeatStructure::new(Nucleotide::A, 10, 2, 2, 1),
                mode: GenerationMode::Uniform,
            },
        }
    }

    fn create_test_population(size: usize, chr_length: usize) -> Population {
        // Minimal mock for population
        let mut individuals = Vec::new();
        for i in 0..size {
            let h1 = Haplotype::from_chromosomes(vec![Chromosome::uniform(
                "c1",
                Nucleotide::A,
                10,
                1,
                chr_length / 10,
            )]);
            let h2 = Haplotype::from_chromosomes(vec![Chromosome::uniform(
                "c1",
                Nucleotide::T,
                10,
                1,
                chr_length / 10,
            )]);
            let mut ind = Individual::new(format!("ind{}", i), h1, h2);
            ind.set_cached_fitness(FitnessValue::new(0.5));
            individuals.push(ind);
        }
        Population::new("test", individuals)
    }

    #[tokio::test]
    async fn test_recorder_flow() {
        let path = "/tmp/test_recorder_new_schema.sqlite";
        let _ = std::fs::remove_file(path);

        let config = mock_config();
        let recorder = AsyncRecorder::new(
            path,
            &config,
            BufferConfig::small(),
            CodecStrategy::default(),
        )
        .expect("Failed to create recorder");

        let pop = create_test_population(5, 100);
        recorder
            .record_generation(&pop, 0, Some(vec![1, 2, 3]))
            .await
            .expect("Failed to record");

        let stats = recorder.close().await.expect("Failed to close");
        assert_eq!(stats.generations_recorded, 1);

        // Verify DB content
        let conn = Connection::open(path).unwrap();
        let count: i64 = conn
            .query_row("SELECT COUNT(*) FROM state", [], |r| r.get(0))
            .unwrap();
        // 5 individuals * 2 haplotypes = 10 rows
        assert_eq!(count, 10);

        let count_meta: i64 = conn
            .query_row("SELECT COUNT(*) FROM metadata", [], |r| r.get(0))
            .unwrap();
        assert!(count_meta > 0);

        std::fs::remove_file(path).ok();
    }
}
