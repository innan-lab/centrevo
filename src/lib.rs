//! Centrevo: A high-performance library for simulating evolution of centromeric sequences.
//!
//! This library provides efficient data structures and algorithms for modeling
//! the evolution of tandemly repeated DNA sequences, particularly centromeric arrays.

pub mod base;
pub mod genome;
pub mod evolution;
// TODO: Uncomment as modules are implemented
// pub mod simulation;
// pub mod utils;
// pub mod storage;

// Re-export commonly used types
pub use base::{Nucleotide, Alphabet, Sequence, SharedSequence};
