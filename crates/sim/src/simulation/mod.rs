//! Simulation engine and population management.
//!
//! This module provides the core simulation loop and population management
//! for evolutionary simulations.

//! Re-exports
//!
//! The most commonly used simulation types are re-exported here for
//! convenience so consumers can import them from `centrevo::simulation`.
//!
//! - `Simulation`: the simulation engine that runs generations and orchestrates
//!   mutation/recombination/selection.
//! - `Population`: in-memory container for individuals used during simulation.
//! - `SimulationBuilder`: fluent builder for constructing `Simulation` instances
//!   with sensible defaults and validation.

pub mod builder;
pub mod configs;
pub mod engine;
pub mod population;
pub mod sequence;

pub use builder::SimulationBuilder;
pub use configs::{
    BuilderEmpty, BuilderInitialized, FitnessConfig, FitnessConfigBuilder, GenerationMode,
    MutationConfig, RecombinationConfig, SequenceConfig, SequenceSource, SimulationConfig,
    UniformRepeatStructure,
};
pub use engine::Simulation;
pub use population::Population;
pub use sequence::{
    SequenceEntry, SequenceEntryWithIndices, create_individuals_from_sequences, initialize,
    load_from_database, parse_fasta, parse_json, validate_sequences,
};
