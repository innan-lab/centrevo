//! Centrevo: A high-performance library for simulating evolution of centromeric sequences.
//!
//! This library provides efficient data structures and algorithms for modeling
//! the evolution of tandemly repeated DNA sequences, particularly centromeric arrays.

pub mod base;
pub mod genome;
pub mod evolution;
pub mod simulation;
pub mod storage;
pub mod analysis;

// Python bindings module (conditionally compiled)
// Use `maturin develop` or `maturin build` to build as Python extension
#[cfg(feature = "python")]
pub mod python;

// TODO: Uncomment as modules are implemented
// pub mod utils;

// Re-export commonly used types
pub use base::{Nucleotide, Alphabet, Sequence, SharedSequence};
