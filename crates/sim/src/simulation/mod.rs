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

pub mod engine;
pub mod population;
pub mod parameters;
pub mod initialization;
pub mod builder;

pub use engine::Simulation;
pub use population::Population;
pub use parameters::{
    RepeatStructure, MutationConfig, RecombinationConfig,
    FitnessConfig, FitnessConfigBuilder, BuilderEmpty, BuilderInitialized,
    SimulationConfig,
};
pub use initialization::{
    SequenceInput, SequenceEntry,
    parse_fasta, parse_json, load_from_database,
    validate_sequences, create_individuals_from_sequences,
    initialize_from_source,
};
pub use builder::SimulationBuilder;
