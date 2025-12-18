//! Storage module for persisting simulation data.
//!
//! This module provides SQLite-based recording of simulation state,
//! allowing for reproducibility and post-simulation analysis.

mod database;
mod query;
mod recorder;
pub mod types;

pub use database::Database;
pub use query::{CheckpointInfo, QueryBuilder};
pub use recorder::{AsyncRecorder, AsyncRecorder as Recorder, BufferConfig, RecorderStats};
pub use types::{FitnessStats, IndividualSnapshot, RecordingStrategy};
