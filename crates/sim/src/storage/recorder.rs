//! High-level recording of simulation state.

use crate::base::FitnessValue;
use crate::errors::DatabaseError;
use crate::genome::Individual;
use crate::simulation::{
    FitnessConfig, MutationConfig, Population, RecombinationConfig, SimulationConfig,
    UniformRepeatStructure,
};
use crate::storage::Database;
use bincode;
use centrevo_codec::CodecStrategy;

use rusqlite::params;
use serde::{Deserialize, Serialize};
use std::sync::Arc;

/// Complete simulation configuration for resumability.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationSnapshot {
    pub structure: UniformRepeatStructure,
    pub mutation: MutationConfig,
    pub recombination: RecombinationConfig,
    pub fitness: FitnessConfig,
    pub config: SimulationConfig,
}

/// Recording strategy for when to persist simulation state.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum RecordingStrategy {
    /// Record every N generations.
    EveryN(usize),

    /// Record at specific generations.
    Specific(Vec<usize>),

    /// Record all generations.
    All,

    /// No recording.
    None,
}

impl RecordingStrategy {
    /// Check if generation should be recorded
    pub fn should_record(&self, generation: usize) -> bool {
        match self {
            Self::EveryN(n) => generation.is_multiple_of(*n),
            Self::Specific(gens) => gens.contains(&generation),
            Self::All => true,
            Self::None => false,
        }
    }
}

/// Snapshot of an individual for database storage.
#[derive(Debug, Clone)]
pub struct IndividualSnapshot {
    pub individual_id: String,
    pub haplotype1_chr_id: String,
    pub haplotype1_seq: Vec<u8>,         // Encoded
    pub haplotype1_map: Option<Vec<u8>>, // Encoded serialized map
    pub haplotype2_chr_id: String,
    pub haplotype2_seq: Vec<u8>,         // Encoded
    pub haplotype2_map: Option<Vec<u8>>, // Encoded serialized map
    pub fitness: Option<f64>,
}

impl IndividualSnapshot {
    /// Create a snapshot from an individual.
    /// Encodes sequence using the provided strategy.
    /// Encodes structure using UnpackedRS (always).
    pub fn from_individual(ind: &Individual, codec: &CodecStrategy) -> Self {
        let h1 = ind.haplotype1();
        let h2 = ind.haplotype2();

        // Helper to process chromosome
        let process_chr =
            |chr_opt: Option<&crate::genome::Chromosome>| -> (String, Vec<u8>, Option<Vec<u8>>) {
                if let Some(chr) = chr_opt {
                    // 1. Sequence -> Configured Codec
                    // Convert Nucleotide -> u8 indices
                    let raw_seq: Vec<u8> = chr
                        .sequence()
                        .as_slice()
                        .iter()
                        .map(|n| n.to_index())
                        .collect();

                    let encoded_seq = codec.encode(&raw_seq).expect("Failed to encode sequence");

                    // 2. Map -> Serialize -> UnpackedRS (Fixed)
                    let map_bytes = bincode::serialize(chr.map()).unwrap_or_default();
                    let encoded_map = CodecStrategy::UnpackedRS.encode(&map_bytes).ok(); // Option

                    (chr.id().to_string(), encoded_seq, encoded_map)
                } else {
                    (String::new(), Vec::new(), None)
                }
            };

        let (h1_chr_id, h1_seq, h1_map) = process_chr(h1.get(0));
        let (h2_chr_id, h2_seq, h2_map) = process_chr(h2.get(0));

        Self {
            individual_id: ind.id().to_string(),
            haplotype1_chr_id: h1_chr_id,
            haplotype1_seq: h1_seq,
            haplotype1_map: h1_map,
            haplotype2_chr_id: h2_chr_id,
            haplotype2_seq: h2_seq,
            haplotype2_map: h2_map,
            fitness: ind.cached_fitness().map(|f| *f),
        }
    }

    /// Reconstruct an Individual from a snapshot.
    /// Uses the provided codec for sequence, and UnpackedRS for map.
    pub fn to_individual(&self, codec: &CodecStrategy) -> Result<Individual, String> {
        use crate::base::{Nucleotide, Sequence};
        use crate::genome::repeat_map::RepeatMap;
        use crate::genome::{Chromosome, Haplotype};

        // Helper to reconstruct chromosome
        let reconstruct_chr =
            |id: &str, seq_data: &[u8], map_data: Option<&[u8]>| -> Result<Chromosome, String> {
                // 1. Decode Sequence
                let raw_seq = codec
                    .decode(seq_data)
                    .map_err(|e| format!("Seq Decode: {e}"))?;

                let nucs: Vec<Nucleotide> = raw_seq
                    .iter()
                    .map(|&i| Nucleotide::from_index(i).unwrap_or(Nucleotide::A))
                    .collect();
                let seq = Sequence::from_nucleotides(nucs);

                // 2. Decode Map
                let map = if let Some(bytes) = map_data {
                    let raw_map_bytes = CodecStrategy::UnpackedRS
                        .decode(bytes)
                        .map_err(|e| format!("Map Decode: {e}"))?;

                    bincode::deserialize::<RepeatMap>(&raw_map_bytes)
                        .map_err(|e| format!("Map Deser: {e}"))?
                } else {
                    return Err("Missing RepeatMap in snapshot".to_string());
                };

                Ok(Chromosome::new(id.to_string(), seq, map))
            };

        let chr1 = reconstruct_chr(
            &self.haplotype1_chr_id,
            &self.haplotype1_seq,
            self.haplotype1_map.as_deref(),
        )?;
        let chr2 = reconstruct_chr(
            &self.haplotype2_chr_id,
            &self.haplotype2_seq,
            self.haplotype2_map.as_deref(),
        )?;

        // Create haplotypes
        let mut hap1 = Haplotype::new();
        hap1.push(chr1);

        let mut hap2 = Haplotype::new();
        hap2.push(chr2);

        // Create individual
        let mut individual = Individual::new(self.individual_id.as_str(), hap1, hap2);
        if let Some(f) = self.fitness {
            individual.set_cached_fitness(FitnessValue::new(f));
        }

        Ok(individual)
    }
}

/// Aggregated fitness statistics for a generation.
#[derive(Debug, Clone, Copy)]
pub struct FitnessStats {
    pub mean: f64,
    pub min: f64,
    pub max: f64,
    pub std: f64,
}

impl FitnessStats {
    /// Calculate fitness statistics from a population.
    pub fn from_population(pop: &Population) -> Self {
        let fitnesses: Vec<f64> = pop
            .individuals()
            .iter()
            .filter_map(|ind| ind.cached_fitness())
            .map(|f| *f)
            .collect();

        if fitnesses.is_empty() {
            return Self {
                mean: 0.0,
                min: 0.0,
                max: 0.0,
                std: 0.0,
            };
        }

        let sum: f64 = fitnesses.iter().sum();
        let mean = sum / fitnesses.len() as f64;

        let min = fitnesses
            .iter()
            .copied()
            .min_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap_or(0.0);

        let max = fitnesses
            .iter()
            .copied()
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap_or(0.0);

        // Calculate standard deviation
        let variance: f64 =
            fitnesses.iter().map(|&f| (f - mean).powi(2)).sum::<f64>() / fitnesses.len() as f64;
        let std = variance.sqrt();

        Self {
            mean,
            min,
            max,
            std,
        }
    }
}

/// High-level recorder for simulation state.
pub struct Recorder {
    db: Database,
    sim_id: Arc<str>,
    strategy: RecordingStrategy,
    generations_recorded: usize,
    codec: CodecStrategy,
}

impl Recorder {
    /// Create a new recorder.
    pub fn new(
        db_path: impl AsRef<std::path::Path>,
        sim_id: impl Into<Arc<str>>,
        strategy: RecordingStrategy,
        codec: CodecStrategy,
    ) -> Result<Self, DatabaseError> {
        let db = Database::open(db_path)?;
        Ok(Self {
            db,
            sim_id: sim_id.into(),
            strategy,
            generations_recorded: 0,
            codec,
        })
    }

    /// Get simulation ID.
    pub fn sim_id(&self) -> &str {
        &self.sim_id
    }

    /// Get recording strategy.
    pub fn strategy(&self) -> &RecordingStrategy {
        &self.strategy
    }

    /// Get number of generations recorded.
    pub fn generations_recorded(&self) -> usize {
        self.generations_recorded
    }

    /// Record simulation metadata.
    pub fn record_metadata(&mut self, config: &SimulationConfig) -> Result<(), DatabaseError> {
        // Serialize configuration to JSON
        let params_json = serde_json::to_string(config).unwrap_or_else(|_| "{}".to_string());

        let start_time = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_secs() as i64;

        self.db
            .connection()
            .execute(
                "INSERT OR REPLACE INTO simulations
                (sim_id, start_time, pop_size, num_generations, mutation_rate, recombination_rate, parameters_json)
                VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7)",
                params![
                    self.sim_id.as_ref(),
                    start_time,
                    config.population_size,
                    config.total_generations,
                    0.0,  // placeholder for mutation rate
                    0.0,  // placeholder for recombination rate
                    params_json,
                ],
            )
            .map_err(|e| DatabaseError::Insert(e.to_string()))?;

        Ok(())
    }

    /// Record complete simulation configuration for resumability.
    pub fn record_full_config(
        &mut self,
        snapshot: &SimulationSnapshot,
    ) -> Result<(), DatabaseError> {
        let config_json = serde_json::to_string(snapshot)
            .map_err(|e| DatabaseError::Insert(format!("Failed to serialize config: {e}")))?;

        let start_time = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_secs() as i64;

        self.db
            .connection()
            .execute(
                "INSERT OR REPLACE INTO simulations
                (sim_id, start_time, pop_size, num_generations, mutation_rate, recombination_rate, parameters_json, config_json)
                VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8)",
                params![
                    self.sim_id.as_ref(),
                    start_time,
                    snapshot.config.population_size,
                    snapshot.config.total_generations,
                    0.0,  // Can extract from mutation config if needed
                    0.0,  // Can extract from recombination config if needed
                    serde_json::to_string(&snapshot.config).unwrap_or_else(|_| "{}".to_string()),
                    config_json,
                ],
            )
            .map_err(|e| DatabaseError::Insert(e.to_string()))?;

        Ok(())
    }

    /// Update simulation end time.
    pub fn finalize_metadata(&mut self) -> Result<(), DatabaseError> {
        let end_time = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_secs() as i64;

        self.db
            .connection()
            .execute(
                "UPDATE simulations SET end_time = ?1 WHERE sim_id = ?2",
                params![end_time, self.sim_id.as_ref()],
            )
            .map_err(|e| DatabaseError::Insert(e.to_string()))?;

        Ok(())
    }

    /// Check if we should record the given generation based on strategy.
    pub fn should_record(&self, generation: usize) -> bool {
        self.strategy.should_record(generation)
    }

    /// Record generation if strategy permits.
    pub fn record_if_needed(
        &mut self,
        population: &Population,
        generation: usize,
        rng_state: &[u8],
    ) -> Result<bool, DatabaseError> {
        if self.should_record(generation) {
            self.record_generation(population, generation, rng_state)?;
            Ok(true)
        } else {
            Ok(false)
        }
    }

    /// Record a generation (forced).
    pub fn record_generation(
        &mut self,
        population: &Population,
        generation: usize,
        rng_state: &[u8],
    ) -> Result<(), DatabaseError> {
        // Prepare snapshots in parallel
        use rayon::prelude::*;
        let codec = &self.codec;
        let snapshots: Vec<IndividualSnapshot> = population
            .individuals()
            .par_iter()
            .map(|ind| IndividualSnapshot::from_individual(ind, codec))
            .collect();

        // Calculate fitness statistics
        let stats = FitnessStats::from_population(population);

        // Write in transaction
        let tx = self.db.transaction()?;

        // Insert population state
        {
            let mut stmt = tx
                .prepare_cached(
                    "INSERT INTO population_state
                    (sim_id, generation, individual_id, haplotype1_chr_id, haplotype1_map, haplotype1_seq,
                     haplotype2_chr_id, haplotype2_map, haplotype2_seq, fitness)
                    VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10)",

                )
                .map_err(|e| DatabaseError::Insert(e.to_string()))?;

            // Clear any existing data for this generation to prevent duplicates
            tx.execute(
                "DELETE FROM population_state WHERE sim_id = ?1 AND generation = ?2",
                params![self.sim_id.as_ref(), generation as i64],
            )
            .map_err(|e| {
                DatabaseError::Insert(format!("Failed to clear existing generation: {e}"))
            })?;

            for snapshot in snapshots {
                stmt.execute(params![
                    self.sim_id.as_ref(),
                    generation as i64,
                    snapshot.individual_id,
                    snapshot.haplotype1_chr_id,
                    snapshot.haplotype1_map,
                    snapshot.haplotype1_seq,
                    snapshot.haplotype2_chr_id,
                    snapshot.haplotype2_map,
                    snapshot.haplotype2_seq,
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
                self.sim_id.as_ref(),
                generation as i64,
                stats.mean,
                stats.min,
                stats.max,
                stats.std,
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
            params![
                self.sim_id.as_ref(),
                generation as i64,
                rng_state,
                timestamp,
            ],
        )
        .map_err(|e| DatabaseError::Insert(e.to_string()))?;

        tx.commit()
            .map_err(|e| DatabaseError::Transaction(e.to_string()))?;

        self.generations_recorded += 1;
        Ok(())
    }

    /// Record a checkpoint with RNG state for resumability.
    /// This should be called alongside record_generation to enable resume functionality.
    pub fn record_checkpoint(
        &mut self,
        generation: usize,
        rng_state: &[u8],
    ) -> Result<(), DatabaseError> {
        let timestamp = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_secs() as i64;

        self.db
            .connection()
            .execute(
                "INSERT OR REPLACE INTO checkpoints
                (sim_id, generation, rng_state, timestamp)
                VALUES (?1, ?2, ?3, ?4)",
                params![
                    self.sim_id.as_ref(),
                    generation as i64,
                    rng_state,
                    timestamp,
                ],
            )
            .map_err(|e| DatabaseError::Insert(e.to_string()))?;

        Ok(())
    }

    /// Close the recorder and clean up.
    pub fn close(self) -> Result<(), DatabaseError> {
        self.db.close()
    }

    /// Get database statistics.
    pub fn stats(&self) -> Result<crate::storage::database::DatabaseStats, DatabaseError> {
        self.db.stats()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::base::Nucleotide;
    use crate::genome::{Chromosome, Haplotype, Individual};

    fn create_test_individual(id: &str, length: usize) -> Individual {
        // Convert total length to num_hors: ru_length=20, rus_per_hor=5, so one HOR = 100 bp
        let num_hors = length / 100;
        let chr1 = Chromosome::uniform(format!("{id}_h1_chr1"), Nucleotide::A, 20, 5, num_hors);
        let chr2 = Chromosome::uniform(format!("{id}_h2_chr1"), Nucleotide::C, 20, 5, num_hors);

        let h1 = Haplotype::from_chromosomes(vec![chr1]);
        let h2 = Haplotype::from_chromosomes(vec![chr2]);

        Individual::new(id, h1, h2)
    }

    fn create_test_population(size: usize, chr_length: usize) -> Population {
        let mut individuals = Vec::new();
        for i in 0..size {
            let mut ind = create_test_individual(&format!("ind_{i}"), chr_length);
            ind.set_cached_fitness(FitnessValue::new(i as f64 / size as f64));
            individuals.push(ind);
        }
        Population::new("test_pop", individuals)
    }

    fn create_test_config() -> SimulationConfig {
        SimulationConfig::new(10, 100, Some(42))
    }

    #[test]
    fn test_recording_strategy() {
        let every_n = RecordingStrategy::EveryN(10);
        assert!(every_n.should_record(0));
        assert!(!every_n.should_record(5));
        assert!(every_n.should_record(10));
        assert!(every_n.should_record(20));

        let specific = RecordingStrategy::Specific(vec![0, 10, 50, 100]);
        assert!(specific.should_record(0));
        assert!(!specific.should_record(5));
        assert!(specific.should_record(10));
        assert!(specific.should_record(100));

        let all = RecordingStrategy::All;
        assert!(all.should_record(0));
        assert!(all.should_record(5));
        assert!(all.should_record(999));

        let none = RecordingStrategy::None;
        assert!(!none.should_record(0));
        assert!(!none.should_record(10));
    }

    #[test]
    fn test_individual_snapshot() {
        let ind = create_test_individual("test_ind", 100);
        let codec = CodecStrategy::BitPackedRS;
        let snapshot = IndividualSnapshot::from_individual(&ind, &codec);

        assert_eq!(snapshot.individual_id, "test_ind");
        assert_eq!(snapshot.individual_id, "test_ind");
        // Length check is harder now as it is compressed.
        assert!(snapshot.haplotype1_seq.len() > 0);
        assert!(snapshot.haplotype2_seq.len() > 0);
        assert!(snapshot.haplotype1_map.is_some());
    }

    #[test]
    fn test_fitness_stats() {
        let pop = create_test_population(10, 100);
        let stats = FitnessStats::from_population(&pop);

        assert!(stats.mean > 0.0);
        assert!(stats.min >= 0.0);
        assert!(stats.max <= 1.0);
        assert!(stats.std >= 0.0);
    }

    #[test]
    fn test_recorder_creation() {
        let path = "/tmp/test_recorder.sqlite";
        let _ = std::fs::remove_file(path);

        let recorder = Recorder::new(
            path,
            "test_sim",
            RecordingStrategy::EveryN(10),
            CodecStrategy::default(),
        )
        .expect("Failed to create recorder");

        assert_eq!(recorder.sim_id(), "test_sim");
        assert_eq!(recorder.generations_recorded(), 0);

        recorder.close().expect("Failed to close recorder");
        std::fs::remove_file(path).ok();
    }

    #[test]
    fn test_record_metadata() {
        let path = "/tmp/test_metadata.sqlite";
        let _ = std::fs::remove_file(path);

        let mut recorder = Recorder::new(
            path,
            "test_sim",
            RecordingStrategy::All,
            CodecStrategy::default(),
        )
        .expect("Failed to create recorder");

        let config = create_test_config();
        recorder
            .record_metadata(&config)
            .expect("Failed to record metadata");

        recorder.close().expect("Failed to close recorder");
        std::fs::remove_file(path).ok();
    }

    #[test]
    fn test_record_generation() {
        let path = "/tmp/test_generation.sqlite";
        let _ = std::fs::remove_file(path);

        let mut recorder = Recorder::new(
            path,
            "test_sim",
            RecordingStrategy::All,
            CodecStrategy::default(),
        )
        .expect("Failed to create recorder");

        let config = create_test_config();
        recorder
            .record_metadata(&config)
            .expect("Failed to record metadata");

        let pop = create_test_population(5, 100);
        let dummy_rng = vec![0u8; 32];
        recorder
            .record_generation(&pop, 0, &dummy_rng)
            .expect("Failed to record generation");

        assert_eq!(recorder.generations_recorded(), 1);

        let stats = recorder.stats().expect("Failed to get stats");
        assert_eq!(stats.population_records, 5);

        recorder.close().expect("Failed to close recorder");
        std::fs::remove_file(path).ok();
    }
}
