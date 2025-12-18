//! Query interface for analyzing recorded simulation data.

use crate::errors::DatabaseError;
use crate::storage::{Database, FitnessStats, IndividualSnapshot};
use rusqlite::params;

/// Query builder for analyzing simulation data.
pub struct QueryBuilder {
    db: Database,
}

impl QueryBuilder {
    /// Open a database for querying.
    pub fn new(db_path: impl AsRef<std::path::Path>) -> Result<Self, DatabaseError> {
        let db = Database::open(db_path)?;
        Ok(Self { db })
    }

    /// List all simulation IDs in the database.
    pub fn list_simulations(&self) -> Result<Vec<String>, DatabaseError> {
        let mut stmt = self
            .db
            .connection()
            .prepare("SELECT sim_id FROM simulations ORDER BY start_time DESC")
            .map_err(|e| DatabaseError::Query(e.to_string()))?;

        let rows = stmt
            .query_map([], |row| row.get::<_, String>(0))
            .map_err(|e| DatabaseError::Query(e.to_string()))?;

        let mut sims = Vec::new();
        for row in rows {
            sims.push(row.map_err(|e| DatabaseError::Query(e.to_string()))?);
        }

        Ok(sims)
    }

    /// Get simulation metadata.
    pub fn get_simulation_info(&self, sim_id: &str) -> Result<SimulationInfo, DatabaseError> {
        let mut stmt = self
            .db
            .connection()
            .prepare(
                "SELECT start_time, end_time, pop_size, num_generations,
                       mutation_rate, recombination_rate, parameters_json
                FROM simulations WHERE sim_id = ?1",
            )
            .map_err(|e| DatabaseError::Query(e.to_string()))?;

        let info = stmt
            .query_row(params![sim_id], |row| {
                Ok(SimulationInfo {
                    sim_id: sim_id.to_string(),
                    start_time: row.get(0)?,
                    end_time: row.get(1)?,
                    pop_size: row.get(2)?,
                    num_generations: row.get(3)?,
                    mutation_rate: row.get(4)?,
                    recombination_rate: row.get(5)?,
                    parameters_json: row.get(6)?,
                })
            })
            .map_err(|e| DatabaseError::Query(e.to_string()))?;

        Ok(info)
    }

    /// Get all individuals at a specific generation.
    pub fn get_generation(
        &self,
        sim_id: &str,
        generation: usize,
    ) -> Result<Vec<IndividualSnapshot>, DatabaseError> {
        let mut stmt = self
            .db
            .connection()
            .prepare(
                "SELECT individual_id, haplotype1_chr_id, haplotype1_seq, haplotype1_map,
                       haplotype2_chr_id, haplotype2_seq, haplotype2_map, haplotype1_fitness, haplotype2_fitness, fitness
                FROM population_state
                WHERE sim_id = ?1 AND generation = ?2
                ORDER BY individual_id",
            )
            .map_err(|e| DatabaseError::Query(e.to_string()))?;

        let rows = stmt
            .query_map(params![sim_id, generation as i64], |row| {
                Ok(IndividualSnapshot {
                    individual_id: row.get(0)?,
                    haplotype1_chr_id: row.get(1)?,
                    haplotype1_seq: row.get(2)?,
                    haplotype1_map: row.get(3)?,
                    haplotype2_chr_id: row.get(4)?,
                    haplotype2_seq: row.get(5)?,
                    haplotype2_map: row.get(6)?,
                    haplotype1_fitness: row.get(7)?,
                    haplotype2_fitness: row.get(8)?,
                    fitness: row.get(9)?,
                })
            })
            .map_err(|e| DatabaseError::Query(e.to_string()))?;

        let mut snapshots = Vec::new();
        for row in rows {
            snapshots.push(row.map_err(|e| DatabaseError::Query(e.to_string()))?);
        }

        Ok(snapshots)
    }

    /// Get fitness history for a simulation.
    pub fn get_fitness_history(
        &self,
        sim_id: &str,
    ) -> Result<Vec<(usize, FitnessStats)>, DatabaseError> {
        let mut stmt = self
            .db
            .connection()
            .prepare(
                "SELECT generation, mean_fitness, min_fitness, max_fitness, std_fitness
                FROM fitness_history
                WHERE sim_id = ?1
                ORDER BY generation",
            )
            .map_err(|e| DatabaseError::Query(e.to_string()))?;

        let rows = stmt
            .query_map(params![sim_id], |row| {
                Ok((
                    row.get::<_, i64>(0)? as usize,
                    FitnessStats {
                        mean: row.get(1)?,
                        min: row.get(2)?,
                        max: row.get(3)?,
                        std: row.get(4)?,
                    },
                ))
            })
            .map_err(|e| DatabaseError::Query(e.to_string()))?;

        let mut history = Vec::new();
        for row in rows {
            history.push(row.map_err(|e| DatabaseError::Query(e.to_string()))?);
        }

        Ok(history)
    }

    /// Get all recorded generations for a simulation.
    pub fn get_recorded_generations(&self, sim_id: &str) -> Result<Vec<usize>, DatabaseError> {
        let mut stmt = self
            .db
            .connection()
            .prepare(
                "SELECT DISTINCT generation FROM population_state
                WHERE sim_id = ?1
                ORDER BY generation",
            )
            .map_err(|e| DatabaseError::Query(e.to_string()))?;

        let rows = stmt
            .query_map(params![sim_id], |row| row.get::<_, i64>(0))
            .map_err(|e| DatabaseError::Query(e.to_string()))?;

        let mut generations = Vec::new();
        for row in rows {
            generations.push(row.map_err(|e| DatabaseError::Query(e.to_string()))? as usize);
        }

        Ok(generations)
    }

    /// Trace an individual across generations (if same ID is maintained).
    pub fn trace_individual(
        &self,
        sim_id: &str,
        individual_id: &str,
    ) -> Result<Vec<(usize, IndividualSnapshot)>, DatabaseError> {
        let mut stmt = self
            .db
            .connection()
            .prepare(
                "SELECT generation, individual_id, haplotype1_chr_id, haplotype1_seq, haplotype1_map,
                       haplotype2_chr_id, haplotype2_seq, haplotype2_map, haplotype1_fitness, haplotype2_fitness, fitness
                FROM population_state
                WHERE sim_id = ?1 AND individual_id = ?2
                ORDER BY generation",
            )
            .map_err(|e| DatabaseError::Query(e.to_string()))?;

        let rows = stmt
            .query_map(params![sim_id, individual_id], |row| {
                Ok((
                    row.get::<_, i64>(0)? as usize,
                    IndividualSnapshot {
                        individual_id: row.get(1)?,
                        haplotype1_chr_id: row.get(2)?,
                        haplotype1_seq: row.get(3)?,
                        haplotype1_map: row.get(4)?,
                        haplotype2_chr_id: row.get(5)?,
                        haplotype2_seq: row.get(6)?,
                        haplotype2_map: row.get(7)?,
                        haplotype1_fitness: row.get(8)?,
                        haplotype2_fitness: row.get(9)?,
                        fitness: row.get(10)?,
                    },
                ))
            })
            .map_err(|e| DatabaseError::Query(e.to_string()))?;

        let mut lineage = Vec::new();
        for row in rows {
            lineage.push(row.map_err(|e| DatabaseError::Query(e.to_string()))?);
        }

        Ok(lineage)
    }

    /// Get fitness statistics for a specific generation.
    pub fn get_generation_fitness(
        &self,
        sim_id: &str,
        generation: usize,
    ) -> Result<FitnessStats, DatabaseError> {
        let mut stmt = self
            .db
            .connection()
            .prepare(
                "SELECT mean_fitness, min_fitness, max_fitness, std_fitness
                FROM fitness_history
                WHERE sim_id = ?1 AND generation = ?2",
            )
            .map_err(|e| DatabaseError::Query(e.to_string()))?;

        let stats = stmt
            .query_row(params![sim_id, generation as i64], |row| {
                Ok(FitnessStats {
                    mean: row.get(0)?,
                    min: row.get(1)?,
                    max: row.get(2)?,
                    std: row.get(3)?,
                })
            })
            .map_err(|e| DatabaseError::Query(e.to_string()))?;

        Ok(stats)
    }

    /// Get individuals with fitness above a threshold at a generation.
    pub fn get_high_fitness_individuals(
        &self,
        sim_id: &str,
        generation: usize,
        min_fitness: f64,
    ) -> Result<Vec<IndividualSnapshot>, DatabaseError> {
        let mut stmt = self
            .db
            .connection()
            .prepare(
                "SELECT individual_id, haplotype1_chr_id, haplotype1_seq, haplotype1_map,
                       haplotype2_chr_id, haplotype2_seq, haplotype2_map, haplotype1_fitness, haplotype2_fitness, fitness
                FROM population_state
                WHERE sim_id = ?1 AND generation = ?2 AND fitness >= ?3
                ORDER BY fitness DESC",
            )
            .map_err(|e| DatabaseError::Query(e.to_string()))?;

        let rows = stmt
            .query_map(params![sim_id, generation as i64, min_fitness], |row| {
                Ok(IndividualSnapshot {
                    individual_id: row.get(0)?,
                    haplotype1_chr_id: row.get(1)?,
                    haplotype1_seq: row.get(2)?,
                    haplotype1_map: row.get(3)?,
                    haplotype2_chr_id: row.get(4)?,
                    haplotype2_seq: row.get(5)?,
                    haplotype2_map: row.get(6)?,
                    haplotype1_fitness: row.get(7)?,
                    haplotype2_fitness: row.get(8)?,
                    fitness: row.get(9)?,
                })
            })
            .map_err(|e| DatabaseError::Query(e.to_string()))?;

        let mut snapshots = Vec::new();
        for row in rows {
            snapshots.push(row.map_err(|e| DatabaseError::Query(e.to_string()))?);
        }

        Ok(snapshots)
    }

    /// Get the latest checkpoint for a simulation.
    pub fn get_latest_checkpoint(&self, sim_id: &str) -> Result<CheckpointInfo, DatabaseError> {
        let mut stmt = self
            .db
            .connection()
            .prepare(
                "SELECT generation, rng_state, timestamp
                FROM checkpoints
                WHERE sim_id = ?1
                ORDER BY generation DESC
                LIMIT 1",
            )
            .map_err(|e| DatabaseError::Query(e.to_string()))?;

        let checkpoint = stmt
            .query_row(params![sim_id], |row| {
                Ok(CheckpointInfo {
                    sim_id: sim_id.to_string(),
                    generation: row.get::<_, i64>(0)? as usize,
                    rng_state: row.get(1)?,
                    timestamp: row.get(2)?,
                })
            })
            .map_err(|e| DatabaseError::Query(e.to_string()))?;

        Ok(checkpoint)
    }

    /// Get a specific checkpoint by generation.
    pub fn get_checkpoint(
        &self,
        sim_id: &str,
        generation: usize,
    ) -> Result<CheckpointInfo, DatabaseError> {
        let mut stmt = self
            .db
            .connection()
            .prepare(
                "SELECT generation, rng_state, timestamp
                FROM checkpoints
                WHERE sim_id = ?1 AND generation = ?2",
            )
            .map_err(|e| DatabaseError::Query(e.to_string()))?;

        let checkpoint = stmt
            .query_row(params![sim_id, generation as i64], |row| {
                Ok(CheckpointInfo {
                    sim_id: sim_id.to_string(),
                    generation: row.get::<_, i64>(0)? as usize,
                    rng_state: row.get(1)?,
                    timestamp: row.get(2)?,
                })
            })
            .map_err(|e| DatabaseError::Query(e.to_string()))?;

        Ok(checkpoint)
    }

    /// Get complete simulation configuration from database.
    pub fn get_full_config(
        &self,
        sim_id: &str,
    ) -> Result<crate::simulation::Configuration, DatabaseError> {
        let mut stmt = self
            .db
            .connection()
            .prepare("SELECT config_json FROM simulations WHERE sim_id = ?1")
            .map_err(|e| DatabaseError::Query(e.to_string()))?;

        let config_json: String = stmt
            .query_row(params![sim_id], |row| row.get(0))
            .map_err(|e| DatabaseError::Query(e.to_string()))?;

        let config: crate::simulation::Configuration = serde_json::from_str(&config_json)
            .map_err(|e| DatabaseError::Query(format!("Failed to parse config: {e}")))?;

        Ok(config)
    }

    /// Close the query builder.
    pub fn close(self) -> Result<(), DatabaseError> {
        self.db.close()
    }
}

/// Simulation metadata information.
#[derive(Debug, Clone)]
pub struct SimulationInfo {
    pub sim_id: String,
    pub start_time: i64,
    pub end_time: Option<i64>,
    pub pop_size: usize,
    pub num_generations: usize,
    pub mutation_rate: f64,
    pub recombination_rate: f64,
    pub parameters_json: String,
}

/// Checkpoint information for resumability.
#[derive(Debug, Clone)]
pub struct CheckpointInfo {
    pub sim_id: String,
    pub generation: usize,
    pub rng_state: Vec<u8>,
    pub timestamp: i64,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::base::{FitnessValue, Nucleotide};
    use crate::genome::{Chromosome, Haplotype, Individual};
    use crate::simulation::{ExecutionConfig, Population};
    use crate::storage::Recorder;

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

    fn create_test_config() -> ExecutionConfig {
        ExecutionConfig::new(10, 100, Some(42))
    }

    async fn setup_test_db(path: &str) -> (Recorder, ExecutionConfig) {
        let _ = std::fs::remove_file(path);
        let buffer_config = crate::storage::BufferConfig {
            compression_level: 0,
            ..Default::default()
        };
        let recorder = Recorder::new(
            path,
            "test_sim",
            buffer_config,
            crate::simulation::CodecStrategy::default(),
        )
        .expect("Failed to create recorder");

        let config = create_test_config();
        recorder
            .record_metadata(&config)
            .await
            .expect("Failed to record metadata");

        // Record a few generations
        for generation in 0..3 {
            let pop = create_test_population(5, 100);
            let dummy_rng = vec![0u8; 32];
            recorder
                .record_generation(&pop, generation, dummy_rng)
                .await
                .expect("Failed to record generation");
        }

        (recorder, config)
    }

    #[tokio::test]
    async fn test_list_simulations() {
        let path = "/tmp/test_query_list.sqlite";
        let (recorder, _config) = setup_test_db(path).await;
        recorder.close().await.ok();

        let query = QueryBuilder::new(path).expect("Failed to create query builder");
        let sims = query
            .list_simulations()
            .expect("Failed to list simulations");

        assert_eq!(sims.len(), 1);
        assert_eq!(sims[0], "test_sim");

        query.close().ok();
        std::fs::remove_file(path).ok();
    }

    #[tokio::test]
    async fn test_get_simulation_info() {
        let path = "/tmp/test_query_info.sqlite";
        let (recorder, _config) = setup_test_db(path).await;
        recorder.close().await.ok();

        let query = QueryBuilder::new(path).expect("Failed to create query builder");
        let info = query
            .get_simulation_info("test_sim")
            .expect("Failed to get simulation info");

        assert_eq!(info.sim_id, "test_sim");
        assert_eq!(info.pop_size, 10);
        assert_eq!(info.num_generations, 100);

        query.close().ok();
        std::fs::remove_file(path).ok();
    }

    #[tokio::test]
    async fn test_get_generation() {
        let path = "/tmp/test_query_gen.sqlite";
        let (recorder, _config) = setup_test_db(path).await;
        recorder.close().await.ok();

        let query = QueryBuilder::new(path).expect("Failed to create query builder");
        let individuals = query
            .get_generation("test_sim", 0)
            .expect("Failed to get generation");

        assert_eq!(individuals.len(), 5);

        query.close().ok();
        std::fs::remove_file(path).ok();
    }

    #[tokio::test]
    async fn test_get_fitness_history() {
        let path = "/tmp/test_query_fitness.sqlite";
        let (recorder, _config) = setup_test_db(path).await;
        recorder.close().await.ok();

        let query = QueryBuilder::new(path).expect("Failed to create query builder");
        let history = query
            .get_fitness_history("test_sim")
            .expect("Failed to get fitness history");

        assert_eq!(history.len(), 3); // 3 generations recorded

        query.close().ok();
        std::fs::remove_file(path).ok();
    }

    #[tokio::test]
    async fn test_get_recorded_generations() {
        let path = "/tmp/test_query_gens.sqlite";
        let (recorder, _config) = setup_test_db(path).await;
        recorder.close().await.ok();

        let query = QueryBuilder::new(path).expect("Failed to create query builder");
        let gens = query
            .get_recorded_generations("test_sim")
            .expect("Failed to get generations");

        assert_eq!(gens, vec![0, 1, 2]);

        query.close().ok();
        std::fs::remove_file(path).ok();
    }
}
