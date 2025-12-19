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
                "-- Metadata table (Configuration items)
                CREATE TABLE IF NOT EXISTS metadata (
                    key TEXT PRIMARY KEY,
                    value TEXT NOT NULL
                );

                -- Normalized Population State
                -- ONE ROW PER CHROMOSOME
                CREATE TABLE IF NOT EXISTS state (
                    generation INTEGER NOT NULL,
                    individual_id INTEGER NOT NULL,
                    haplotype_id INTEGER NOT NULL, -- 1 or 2
                    chromosome_id INTEGER NOT NULL, -- 1-based index
                    seq BLOB,       -- Encoded sequence
                    map BLOB,       -- Encoded map structure
                    PRIMARY KEY (generation, individual_id, haplotype_id, chromosome_id)
                );

                -- Individual Fitness
                CREATE TABLE IF NOT EXISTS individual_fitness (
                    generation INTEGER NOT NULL,
                    individual_id INTEGER NOT NULL,
                    fitness REAL,
                    timestamp INTEGER,
                    PRIMARY KEY (generation, individual_id)
                );

                -- Checkpoints for resumability (RNG state)
                CREATE TABLE IF NOT EXISTS checkpoints (
                    generation INTEGER PRIMARY KEY,
                    rng_state BLOB NOT NULL,
                    timestamp INTEGER NOT NULL
                );

                -- Indices
                CREATE INDEX IF NOT EXISTS idx_state_gen ON state(generation);
                CREATE INDEX IF NOT EXISTS idx_fitness_gen ON individual_fitness(generation);",
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
        let state_count: i64 = self
            .conn
            .query_row("SELECT COUNT(*) FROM state", [], |row| row.get(0))
            .unwrap_or(0);

        let fitness_count: i64 = self
            .conn
            .query_row("SELECT COUNT(*) FROM individual_fitness", [], |row| {
                row.get(0)
            })
            .unwrap_or(0);

        Ok(DatabaseStats {
            state_records: state_count as usize,
            fitness_records: fitness_count as usize,
            tables,
        })
    }
}

/// Database statistics.
#[derive(Debug, Clone)]
pub struct DatabaseStats {
    pub state_records: usize,
    pub fitness_records: usize,
    pub tables: Vec<String>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn test_database_creation() {
        let path = "/tmp/test_centrevo_db_v2.sqlite";
        let _ = fs::remove_file(path); // Clean up if exists

        let db = Database::open(path).expect("Failed to create database");
        assert_eq!(db.path(), path);

        db.close().expect("Failed to close database");
        fs::remove_file(path).expect("Failed to remove test database");
    }

    #[test]
    fn test_schema_initialization() {
        let path = "/tmp/test_centrevo_schema_v2.sqlite";
        let _ = fs::remove_file(path);

        let db = Database::open(path).expect("Failed to create database");
        let stats = db.stats().expect("Failed to get stats");

        assert!(stats.tables.contains(&"metadata".to_string()));
        assert!(stats.tables.contains(&"state".to_string()));
        assert!(stats.tables.contains(&"individual_fitness".to_string()));
        assert!(stats.tables.contains(&"checkpoints".to_string()));

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
