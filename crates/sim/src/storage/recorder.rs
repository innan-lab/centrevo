//! High-level recording of simulation state.

use crate::errors::DatabaseError;
use crate::genome::Individual;
use crate::simulation::{
    FitnessConfig, MutationConfig, Population, RecombinationConfig, RepeatStructure,
    SimulationConfig,
};
use crate::storage::Database;
use rusqlite::params;
use serde::{Deserialize, Serialize};
use std::sync::Arc;

/// Complete simulation configuration for resumability.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationSnapshot {
    pub structure: RepeatStructure,
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
    pub haplotype1_seq: Vec<u8>,
    pub haplotype2_chr_id: String,
    pub haplotype2_seq: Vec<u8>,
    pub fitness: f64,
}

impl IndividualSnapshot {
    /// Create a snapshot from an individual.
    /// For simplicity, we record the first chromosome of each haplotype.
    /// Can be extended to record all chromosomes if needed.
    pub fn from_individual(ind: &Individual) -> Self {
        let h1 = ind.haplotype1();
        let h2 = ind.haplotype2();

        // Get first chromosome from each haplotype
        let (h1_chr_id, h1_seq) = if let Some(chr) = h1.get(0) {
            (
                chr.id().to_string(),
                chr.sequence()
                    .as_slice()
                    .iter()
                    .map(|n| n.to_index())
                    .collect(),
            )
        } else {
            (String::new(), Vec::new())
        };

        let (h2_chr_id, h2_seq) = if let Some(chr) = h2.get(0) {
            (
                chr.id().to_string(),
                chr.sequence()
                    .as_slice()
                    .iter()
                    .map(|n| n.to_index())
                    .collect(),
            )
        } else {
            (String::new(), Vec::new())
        };

        Self {
            individual_id: ind.id().to_string(),
            haplotype1_chr_id: h1_chr_id,
            haplotype1_seq: h1_seq,
            haplotype2_chr_id: h2_chr_id,
            haplotype2_seq: h2_seq,
            fitness: ind.fitness(),
        }
    }

    /// Reconstruct an Individual from a snapshot.
    /// Requires the RepeatStructure to rebuild the proper chromosome structure.
    pub fn to_individual(
        &self,
        structure: &crate::simulation::RepeatStructure,
    ) -> Result<Individual, String> {
        use crate::base::{Nucleotide, Sequence};
        use crate::genome::{Chromosome, Haplotype};
        use crate::genome::repeat_map::RepeatMap;

        // Reconstruct sequence from indices
        let seq1_nucs: Vec<Nucleotide> = self
            .haplotype1_seq
            .iter()
            .map(|&i| Nucleotide::from_index(i).unwrap_or(Nucleotide::A))
            .collect();
        let seq1 = Sequence::from_nucleotides(seq1_nucs);

        let seq2_nucs: Vec<Nucleotide> = self
            .haplotype2_seq
            .iter()
            .map(|&i| Nucleotide::from_index(i).unwrap_or(Nucleotide::A))
            .collect();
        let seq2 = Sequence::from_nucleotides(seq2_nucs);

        // Create uniform map
        let map = RepeatMap::uniform(
            structure.ru_length,
            structure.rus_per_hor,
            structure.hors_per_chr,
        );

        // Create chromosomes
        let chr1 = Chromosome::new(
            self.haplotype1_chr_id.clone(),
            seq1,
            map.clone(),
        );
        let chr2 = Chromosome::new(
            self.haplotype2_chr_id.clone(),
            seq2,
            map,
        );

        // Create haplotypes
        let mut hap1 = Haplotype::new();
        hap1.push(chr1);

        let mut hap2 = Haplotype::new();
        hap2.push(chr2);

        // Create individual
        let mut individual = Individual::new(self.individual_id.as_str(), hap1, hap2);
        individual.set_fitness(self.fitness);

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
        let fitnesses: Vec<f64> = pop.individuals().iter().map(|ind| ind.fitness()).collect();

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
}

impl Recorder {
    /// Create a new recorder.
    pub fn new(
        db_path: impl AsRef<std::path::Path>,
        sim_id: impl Into<Arc<str>>,
        strategy: RecordingStrategy,
    ) -> Result<Self, DatabaseError> {
        let db = Database::open(db_path)?;
        Ok(Self {
            db,
            sim_id: sim_id.into(),
            strategy,
            generations_recorded: 0,
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
    ) -> Result<bool, DatabaseError> {
        if self.should_record(generation) {
            self.record_generation(population, generation)?;
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
    ) -> Result<(), DatabaseError> {
        // Prepare snapshots in parallel
        use rayon::prelude::*;
        let snapshots: Vec<IndividualSnapshot> = population
            .individuals()
            .par_iter()
            .map(IndividualSnapshot::from_individual)
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
                    (sim_id, generation, individual_id, haplotype1_chr_id, haplotype1_seq, 
                     haplotype2_chr_id, haplotype2_seq, fitness)
                    VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8)",
                )
                .map_err(|e| DatabaseError::Insert(e.to_string()))?;

            for snapshot in snapshots {
                stmt.execute(params![
                    self.sim_id.as_ref(),
                    generation as i64,
                    snapshot.individual_id,
                    snapshot.haplotype1_chr_id,
                    snapshot.haplotype1_seq,
                    snapshot.haplotype2_chr_id,
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
        let chr1 = Chromosome::uniform(format!("{}_h1_chr1", id), Nucleotide::A, length, 20, 5);
        let chr2 = Chromosome::uniform(format!("{}_h2_chr1", id), Nucleotide::C, length, 20, 5);

        let h1 = Haplotype::from_chromosomes(vec![chr1]);
        let h2 = Haplotype::from_chromosomes(vec![chr2]);

        Individual::new(id, h1, h2)
    }

    fn create_test_population(size: usize, chr_length: usize) -> Population {
        let mut individuals = Vec::new();
        for i in 0..size {
            let mut ind = create_test_individual(&format!("ind_{}", i), chr_length);
            ind.set_fitness(i as f64 / size as f64);
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
        let snapshot = IndividualSnapshot::from_individual(&ind);

        assert_eq!(snapshot.individual_id, "test_ind");
        assert_eq!(snapshot.haplotype1_seq.len(), 100);
        assert_eq!(snapshot.haplotype2_seq.len(), 100);
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

        let recorder = Recorder::new(path, "test_sim", RecordingStrategy::EveryN(10))
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

        let mut recorder = Recorder::new(path, "test_sim", RecordingStrategy::All)
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

        let mut recorder = Recorder::new(path, "test_sim", RecordingStrategy::All)
            .expect("Failed to create recorder");

        let config = create_test_config();
        recorder
            .record_metadata(&config)
            .expect("Failed to record metadata");

        let pop = create_test_population(5, 100);
        recorder
            .record_generation(&pop, 0)
            .expect("Failed to record generation");

        assert_eq!(recorder.generations_recorded(), 1);

        let stats = recorder.stats().expect("Failed to get stats");
        assert_eq!(stats.population_records, 5);

        recorder.close().expect("Failed to close recorder");
        std::fs::remove_file(path).ok();
    }
}
