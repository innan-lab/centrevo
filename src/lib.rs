//! Centrevo: A high-performance library for simulating evolution of centromeric sequences.
//!
//! This library provides efficient data structures and algorithms for modeling
//! the evolution of tandemly repeated DNA sequences, particularly centromeric arrays.

pub mod analysis;
pub mod base;
pub mod evolution;
pub mod genome;
pub mod simulation;
pub mod storage;

// Python bindings module (conditionally compiled)
// Use `maturin develop` or `maturin build` to build as Python extension
#[cfg(feature = "python")]
pub mod python;

// TODO: Uncomment as modules are implemented
// pub mod utils;

// Re-export commonly used types for convenient external access.
//
// These types form the public, stable surface that most consumers of the
// library will use when building analyses or running simulations. Re-exporting
// them here makes them available as `centrevo::Sequence`, `centrevo::Nucleotide`,
// etc.
pub use base::{Nucleotide, Sequence, SharedSequence};
