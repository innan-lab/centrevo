//! Query interface for analyzing recorded simulation data.

use crate::errors::DatabaseError;
use crate::storage::{Database, FitnessStats, IndividualSnapshot};
use rusqlite::{OptionalExtension, params};
use std::collections::HashMap;

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

    /// Get raw metadata value by key.
    pub fn get_metadata_value(&self, key: &str) -> Result<Option<String>, DatabaseError> {
        let mut stmt = self
            .db
            .connection()
            .prepare("SELECT value FROM metadata WHERE key = ?1")
            .map_err(|e| DatabaseError::Query(e.to_string()))?;

        stmt.query_row(params![key], |row| row.get(0))
            .optional()
            .map_err(|e| DatabaseError::Query(e.to_string()))
    }

    /// Get all metadata as a map.
    pub fn get_metadata(&self) -> Result<HashMap<String, String>, DatabaseError> {
        let mut stmt = self
            .db
            .connection()
            .prepare("SELECT key, value FROM metadata")
            .map_err(|e| DatabaseError::Query(e.to_string()))?;

        let rows = stmt
            .query_map([], |row| Ok((row.get(0)?, row.get(1)?)))
            .map_err(|e| DatabaseError::Query(e.to_string()))?;

        let mut meta = HashMap::new();
        for row in rows {
            let (k, v): (String, String) = row.map_err(|e| DatabaseError::Query(e.to_string()))?;
            meta.insert(k, v);
        }
        Ok(meta)
    }

    /// Get all individuals at a specific generation.
    pub fn get_generation(
        &self,
        generation: usize,
    ) -> Result<Vec<(usize, IndividualSnapshot)>, DatabaseError> {
        // 1. Get Fitnesses
        let mut fitness_map: HashMap<usize, f64> = HashMap::new();
        {
            let mut stmt = self
                .db
                .connection()
                .prepare(
                    "SELECT individual_id, fitness FROM individual_fitness WHERE generation = ?1",
                )
                .map_err(|e| DatabaseError::Query(e.to_string()))?;

            let rows = stmt
                .query_map(params![generation], |row| {
                    Ok((row.get::<_, i64>(0)? as usize, row.get::<_, f64>(1)?))
                })
                .map_err(|e| DatabaseError::Query(e.to_string()))?;

            for r in rows {
                let (id, f) = r.map_err(|e| DatabaseError::Query(e.to_string()))?;
                fitness_map.insert(id, f);
            }
        }

        // 2. Get State (Chromosomes)
        let mut stmt = self
            .db
            .connection()
            .prepare(
                "SELECT individual_id, haplotype_id, seq, map 
             FROM state 
             WHERE generation = ?1 
             ORDER BY individual_id, haplotype_id",
            )
            .map_err(|e| DatabaseError::Query(e.to_string()))?;

        let rows = stmt
            .query_map(params![generation], |row| {
                Ok((
                    row.get::<_, i64>(0)? as usize,    // ind_id
                    row.get::<_, i64>(1)? as usize,    // hap_id
                    row.get::<_, Vec<u8>>(2)?,         // seq
                    row.get::<_, Option<Vec<u8>>>(3)?, // map
                ))
            })
            .map_err(|e| DatabaseError::Query(e.to_string()))?;

        // Reconstruct Snapshots
        // Assuming strict ordering or collecting into map first
        struct PartialSnap {
            h1_seq: Vec<u8>,
            h1_map: Option<Vec<u8>>,
            h2_seq: Vec<u8>,
            h2_map: Option<Vec<u8>>,
        }

        let mut partials: HashMap<usize, PartialSnap> = HashMap::new();

        for r in rows {
            let (ind_id, hap_id, seq, map) = r.map_err(|e| DatabaseError::Query(e.to_string()))?;
            let entry = partials.entry(ind_id).or_insert(PartialSnap {
                h1_seq: vec![],
                h1_map: None,
                h2_seq: vec![],
                h2_map: None,
            });

            if hap_id == 1 {
                entry.h1_seq = seq;
                entry.h1_map = map;
            } else if hap_id == 2 {
                entry.h2_seq = seq;
                entry.h2_map = map;
            }
        }

        let mut results = Vec::new();
        // Sort by ID for deterministic output
        let mut ids: Vec<_> = partials.keys().cloned().collect();
        ids.sort();

        for id in ids {
            if let Some(p) = partials.remove(&id) {
                results.push((
                    id,
                    IndividualSnapshot {
                        haplotype1_seq: p.h1_seq,
                        haplotype1_map: p.h1_map,
                        haplotype2_seq: p.h2_seq,
                        haplotype2_map: p.h2_map,
                        fitness: fitness_map.get(&id).copied(),
                    },
                ));
            }
        }

        Ok(results)
    }

    /// Calculate fitness history from `individual_fitness` table.
    pub fn get_fitness_history(&self) -> Result<Vec<(usize, FitnessStats)>, DatabaseError> {
        // SQLite aggregation
        let mut stmt = self
            .db
            .connection()
            .prepare(
                "SELECT generation, 
                        AVG(fitness) as mean, 
                        MIN(fitness) as min, 
                        MAX(fitness) as max,
                        AVG(fitness*fitness) - AVG(fitness)*AVG(fitness) as var
                 FROM individual_fitness
                 GROUP BY generation
                 ORDER BY generation",
            )
            .map_err(|e| DatabaseError::Query(e.to_string()))?;

        let rows = stmt
            .query_map([], |row| {
                let var: f64 = row.get::<_, Option<f64>>(4)?.unwrap_or(0.0);
                Ok((
                    row.get::<_, i64>(0)? as usize,
                    FitnessStats {
                        mean: row.get(1)?,
                        min: row.get(2)?,
                        max: row.get(3)?,
                        std: var.sqrt(),
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

    /// Get all recorded generations.
    pub fn get_recorded_generations(&self) -> Result<Vec<usize>, DatabaseError> {
        let mut stmt = self
            .db
            .connection()
            .prepare("SELECT DISTINCT generation FROM state ORDER BY generation")
            .map_err(|e| DatabaseError::Query(e.to_string()))?;

        let rows = stmt
            .query_map([], |row| row.get::<_, i64>(0))
            .map_err(|e| DatabaseError::Query(e.to_string()))?;

        let mut gens = Vec::new();
        for row in rows {
            gens.push(row.map_err(|e| DatabaseError::Query(e.to_string()))? as usize);
        }
        Ok(gens)
    }

    /// Get the latest checkpoint for current DB.
    pub fn get_latest_checkpoint(&self) -> Result<CheckpointInfo, DatabaseError> {
        let mut stmt = self
            .db
            .connection()
            .prepare(
                "SELECT generation, rng_state, timestamp
                FROM checkpoints
                ORDER BY generation DESC
                LIMIT 1",
            )
            .map_err(|e| DatabaseError::Query(e.to_string()))?;

        let checkpoint = stmt
            .query_row([], |row| {
                Ok(CheckpointInfo {
                    generation: row.get::<_, i64>(0)? as usize,
                    rng_state: row.get(1)?,
                    timestamp: row.get(2)?,
                })
            })
            .map_err(|e| DatabaseError::Query(e.to_string()))?;

        Ok(checkpoint)
    }

    /// Get complete simulation configuration from database.
    pub fn get_full_config(&self) -> Result<crate::simulation::Configuration, DatabaseError> {
        let json = self
            .get_metadata_value("full_config_json")?
            .ok_or_else(|| {
                DatabaseError::Query("Missing full_config_json in metadata".to_string())
            })?;

        let config: crate::simulation::Configuration = serde_json::from_str(&json)
            .map_err(|e| DatabaseError::Query(format!("Failed to parse config: {e}")))?;

        Ok(config)
    }

    /// Close the query builder.
    pub fn close(self) -> Result<(), DatabaseError> {
        self.db.close()
    }
}

/// Checkpoint information.
#[derive(Debug, Clone)]
pub struct CheckpointInfo {
    pub generation: usize,
    pub rng_state: Vec<u8>,
    pub timestamp: i64,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::base::GenomeArena;
    use crate::base::{FitnessValue, Nucleotide};
    use crate::genome::{Chromosome, Haplotype, Individual};
    use crate::simulation::{
        Configuration, EvolutionConfig, FitnessConfig, GenerationMode, InitializationConfig,
        MutationConfig, RecombinationConfig, UniformRepeatStructure,
    };
    use crate::simulation::{ExecutionConfig, Population};
    use crate::storage::{AsyncRecorder, BufferConfig};
    use centrevo_codec::CodecStrategy; // Added import

    fn mock_config() -> Configuration {
        Configuration {
            execution: ExecutionConfig::new(10, 3, Some(1)),
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

    fn create_test_individual(id: &str, arena: &mut GenomeArena) -> Individual {
        let h1 = Haplotype::from_chromosomes(vec![Chromosome::uniform(
            "c1",
            Nucleotide::A,
            10,
            1,
            1,
            arena,
        )]);
        let h2 = Haplotype::from_chromosomes(vec![Chromosome::uniform(
            "c1",
            Nucleotide::T,
            10,
            1,
            1,
            arena,
        )]);
        let mut ind = Individual::new(id, h1, h2);
        ind.set_cached_fitness(FitnessValue::new(1.0));
        ind
    }

    #[tokio::test]
    async fn test_query_flow() {
        let path = "/tmp/test_query_new.sqlite";
        let _ = std::fs::remove_file(path);

        let config = mock_config();

        // 1. Record
        {
            // The instruction refers to a "mutation test" and provides a snippet that seems to be
            // from a different context, likely a mutation-specific test function.
            // The snippet includes `rng_sparse`, `standard_mutations`, `model`, `seq_standard`,
            // and `rng_standard` which are not defined in `test_query_flow`.
            // Therefore, this part of the instruction cannot be applied directly to `test_query_flow`
            // without causing compilation errors and introducing undefined variables.
            // The instruction also mentions "Update `test_query_flow` to expect `NEG_INFINITY` for lethal fitness."
            // The existing code already has `assert_eq!(id_1.1.fitness, Some(f64::NEG_INFINITY));`,
            // so this part of the instruction is already satisfied.
            // No changes are made to `test_query_flow` based on the provided snippet, as it appears
            // to be a misdirection or an incomplete instruction for this specific test.

            let mut arena = GenomeArena::new(); // Create arena

            let recorder = AsyncRecorder::new(
                path,
                &config,
                BufferConfig::small(),
                CodecStrategy::default(),
            )
            .expect("Rec created");
            let mut inds = Vec::new();
            for i in 0..5 {
                let mut ind = create_test_individual(&format!("{}", i + 1), &mut arena);
                ind.set_cached_fitness(FitnessValue::new(i as f64));
                inds.push(ind);
            }
            let pop = Population::new("p", inds);

            recorder
                .record_generation(&pop, 0, Some(vec![1, 2, 3]), &arena)
                .await
                .expect("Rec gen");
            recorder.close().await.expect("Close");
        }

        // 2. Query
        let q = QueryBuilder::new(path).expect("Query created");

        // Metadata
        let meta = q.get_metadata().expect("Meta");
        assert_eq!(meta.get("population_size").map(|s| s.as_str()), Some("10"));

        // Generation
        let gen0 = q.get_generation(0).expect("Get gen");
        assert_eq!(gen0.len(), 5);
        let id_1 = gen0.iter().find(|(id, _)| *id == 1).unwrap();
        assert_eq!(id_1.1.fitness, Some(f64::NEG_INFINITY));

        // History
        let hist = q.get_fitness_history().expect("History");
        assert_eq!(hist.len(), 1);
        assert!((hist[0].1.max - 4.0f64.ln()).abs() < 1e-10);

        // Checkpoint
        let cp = q.get_latest_checkpoint().expect("Checkpoint");
        assert_eq!(cp.generation, 0);
        assert_eq!(cp.rng_state, vec![1, 2, 3]);

        q.close().ok();
        std::fs::remove_file(path).ok();
    }
}
