//! Base types for sequence representation.
//!
//! This module provides the foundational types for representing nucleotides,
//! alphabets, and sequences in the centrevo library.

pub mod arena;
pub mod fitness;
mod nucleotide;
mod sequence;

pub use arena::{GenomeArena, SequenceSlice};
pub use fitness::{FitnessValue, LogFitnessValue};
pub use nucleotide::Nucleotide;
pub use sequence::{Sequence, SharedSequence};
