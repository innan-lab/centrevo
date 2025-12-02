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
pub mod engine;
pub mod inputs;
pub mod parameters;
pub mod population;

pub use builder::SimulationBuilder;
pub use engine::Simulation;
pub use inputs::{
    SequenceEntry, SequenceInput, create_individuals_from_sequences, initialize_from_source,
    load_from_database, parse_fasta, parse_json, validate_sequences,
};
pub use parameters::{
    BuilderEmpty, BuilderInitialized, FitnessConfig, FitnessConfigBuilder, MutationConfig,
    RecombinationConfig, SimulationConfig, UniformRepeatStructure,
};
pub use population::Population;
