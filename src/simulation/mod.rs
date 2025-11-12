//! Simulation engine and population management.
//!
//! This module provides the core simulation loop and population management
//! for evolutionary simulations.

pub mod engine;
pub mod population;
pub mod parameters;
pub mod initialization;

pub use engine::Simulation;
pub use population::Population;
pub use parameters::{
    RepeatStructure, MutationConfig, RecombinationConfig,
    FitnessConfig, SimulationConfig,
};
pub use initialization::{
    SequenceInput, SequenceEntry, InitializationError,
    parse_fasta, parse_json, load_from_database,
    validate_sequences, create_individuals_from_sequences,
    initialize_from_source,
};
