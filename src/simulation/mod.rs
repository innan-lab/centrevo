//! Simulation engine and population management.
//!
//! This module provides the core simulation loop and population management
//! for evolutionary simulations.

pub mod engine;
pub mod population;
pub mod parameters;

pub use engine::Simulation;
pub use population::Population;
pub use parameters::{
    RepeatStructure, MutationConfig, RecombinationConfig,
    FitnessConfig, SimulationConfig,
};
