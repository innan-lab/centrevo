//! Storage module for persisting simulation data.
//!
//! This module provides SQLite-based recording of simulation state,
//! allowing for reproducibility and post-simulation analysis.

mod async_recorder;
mod database;
mod query;
pub mod types;

pub use async_recorder::{AsyncRecorder, AsyncRecorder as Recorder, BufferConfig, RecorderStats};
pub use database::Database;
pub use query::{CheckpointInfo, QueryBuilder};
pub use types::{FitnessStats, IndividualSnapshot, RecordingStrategy, SimulationSnapshot};
