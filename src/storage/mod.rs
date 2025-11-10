//! Storage module for persisting simulation data.
//!
//! This module provides SQLite-based recording of simulation state,
//! allowing for reproducibility and post-simulation analysis.

mod database;
mod recorder;
mod query;

pub use database::{Database, DatabaseError};
pub use recorder::{Recorder, RecordingStrategy, IndividualSnapshot, FitnessStats};
pub use query::QueryBuilder;
