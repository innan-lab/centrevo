//! Low-level database operations and schema management.

pub use crate::errors::DatabaseError;
use rusqlite::{Connection, Transaction};
use std::path::Path;

/// Database connection wrapper with schema management.
#[derive(Debug)]
pub struct Database {
    conn: Connection,
    db_path: String,
}

impl Database {
    /// Open (or create) a database at the specified path.
    pub fn open(path: impl AsRef<Path>) -> Result<Self, DatabaseError> {
        let path_str = path.as_ref().to_string_lossy().to_string();
        let conn =
            Connection::open(&path_str).map_err(|e| DatabaseError::Connection(e.to_string()))?;

        // Performance pragmas for faster bulk inserts
        conn.execute_batch(
            "PRAGMA synchronous = NORMAL;
             PRAGMA journal_mode = WAL;
             PRAGMA temp_store = MEMORY;
             PRAGMA cache_size = -64000;",
        )
        .map_err(|e| DatabaseError::Initialization(e.to_string()))?;

        let mut db = Self {
            conn,
            db_path: path_str,
        };

        db.initialize_schema()?;
        Ok(db)
    }

    /// Initialize database schema.
    fn initialize_schema(&mut self) -> Result<(), DatabaseError> {
        self.conn
            .execute_batch(
                "-- Main population state table
                CREATE TABLE IF NOT EXISTS population_state (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    sim_id TEXT NOT NULL,
                    generation INTEGER NOT NULL,
                    individual_id TEXT NOT NULL,
                    haplotype1_chr_id TEXT NOT NULL,
                    haplotype1_map BLOB,
                    haplotype1_seq BLOB NOT NULL,
                    haplotype2_chr_id TEXT NOT NULL,
                    haplotype2_map BLOB,
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
                    -- Full configuration for resumability
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
                    rng_state BLOB NOT NULL,  -- 32 bytes: 4x u64 for Xoshiro256PlusPlus
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

        Ok(())
    }

    /// Begin a transaction for batched operations.
    pub fn transaction(&mut self) -> Result<Transaction<'_>, DatabaseError> {
        self.conn
            .transaction()
            .map_err(|e| DatabaseError::Transaction(e.to_string()))
    }

    /// Get reference to underlying connection.
    pub fn connection(&self) -> &Connection {
        &self.conn
    }

    /// Get mutable reference to underlying connection.
    pub fn connection_mut(&mut self) -> &mut Connection {
        &mut self.conn
    }

    /// Get database path.
    pub fn path(&self) -> &str {
        &self.db_path
    }

    /// Close the database and clean up WAL files.
    pub fn close(self) -> Result<(), DatabaseError> {
        // Checkpoint and truncate WAL
        if let Err(e) = self.conn.execute_batch(
            "PRAGMA wal_checkpoint(TRUNCATE);
             PRAGMA journal_mode = DELETE;",
        ) {
            eprintln!("Warning: failed to checkpoint/truncate WAL: {e}");
        }

        // Close connection
        self.conn
            .close()
            .map_err(|(_conn, e)| DatabaseError::Close(e.to_string()))?;

        // Remove WAL/SHM files
        for suffix in &["-wal", "-shm"] {
            let fname = format!("{}{}", self.db_path, suffix);
            if let Err(e) = std::fs::remove_file(&fname) {
                if e.kind() == std::io::ErrorKind::NotFound {
                    // File doesn't exist, that's fine
                } else {
                    eprintln!("Warning: failed to remove {fname}: {e}");
                }
            }
        }

        Ok(())
    }

    /// Vacuum the database to reclaim space.
    pub fn vacuum(&self) -> Result<(), DatabaseError> {
        self.conn
            .execute_batch("VACUUM;")
            .map_err(|e| DatabaseError::Vacuum(e.to_string()))
    }

    /// Get database statistics.
    pub fn stats(&self) -> Result<DatabaseStats, DatabaseError> {
        let mut stmt = self
            .conn
            .prepare("SELECT name, COUNT(*) FROM sqlite_master WHERE type='table' GROUP BY name")
            .map_err(|e| DatabaseError::Query(e.to_string()))?;

        let mut tables = Vec::new();
        let rows = stmt
            .query_map([], |row| {
                Ok((row.get::<_, String>(0)?, row.get::<_, i64>(1)?))
            })
            .map_err(|e| DatabaseError::Query(e.to_string()))?;

        for row in rows {
            let (name, _count) = row.map_err(|e| DatabaseError::Query(e.to_string()))?;
            tables.push(name);
        }

        // Get row counts for main tables
        let pop_count: i64 = self
            .conn
            .query_row("SELECT COUNT(*) FROM population_state", [], |row| {
                row.get(0)
            })
            .unwrap_or(0);

        let sim_count: i64 = self
            .conn
            .query_row("SELECT COUNT(*) FROM simulations", [], |row| row.get(0))
            .unwrap_or(0);

        Ok(DatabaseStats {
            population_records: pop_count as usize,
            simulation_count: sim_count as usize,
            tables,
        })
    }
}

/// Database statistics.
#[derive(Debug, Clone)]
pub struct DatabaseStats {
    pub population_records: usize,
    pub simulation_count: usize,
    pub tables: Vec<String>,
}

// Removed DatabaseError definition, imported from crate::errors

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn test_database_creation() {
        let path = "/tmp/test_centrevo_db.sqlite";
        let _ = fs::remove_file(path); // Clean up if exists

        let db = Database::open(path).expect("Failed to create database");
        assert_eq!(db.path(), path);

        db.close().expect("Failed to close database");
        fs::remove_file(path).expect("Failed to remove test database");
    }

    #[test]
    fn test_schema_initialization() {
        let path = "/tmp/test_centrevo_schema.sqlite";
        let _ = fs::remove_file(path);

        let db = Database::open(path).expect("Failed to create database");
        let stats = db.stats().expect("Failed to get stats");

        assert!(stats.tables.contains(&"population_state".to_string()));
        assert!(stats.tables.contains(&"simulations".to_string()));
        assert!(stats.tables.contains(&"fitness_history".to_string()));

        db.close().expect("Failed to close database");
        fs::remove_file(path).expect("Failed to remove test database");
    }

    #[test]
    fn test_transaction() {
        let path = "/tmp/test_centrevo_tx.sqlite";
        let _ = fs::remove_file(path);

        let mut db = Database::open(path).expect("Failed to create database");
        let tx = db.transaction().expect("Failed to begin transaction");
        tx.commit().expect("Failed to commit transaction");

        db.close().expect("Failed to close database");
        fs::remove_file(path).expect("Failed to remove test database");
    }
}
