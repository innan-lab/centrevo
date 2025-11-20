//! # Simulation Crate
//!
//! The `sim` crate provides the core logic for the evolutionary simulation.
//! It includes modules for defining genomes, managing populations,
//! executing evolutionary operators (mutation, recombination, selection),
//! and running the simulation engine.

pub mod base;
pub mod errors;
pub mod evolution;
pub mod genome;
pub mod simulation;
pub mod storage;
pub mod prelude;

pub use base::{Nucleotide, Sequence, SharedSequence};
