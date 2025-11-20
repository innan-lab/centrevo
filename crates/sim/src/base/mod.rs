//! Base types for sequence representation.
//!
//! This module provides the foundational types for representing nucleotides,
//! alphabets, and sequences in the centrevo library.

mod nucleotide;
mod sequence;

pub use nucleotide::Nucleotide;
pub use sequence::{Sequence, SharedSequence};
