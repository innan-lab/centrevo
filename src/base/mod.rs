//! Base types for sequence representation.
//!
//! This module provides the foundational types for representing nucleotides,
//! alphabets, and sequences in the centrevo library.

mod nucleotide;
mod alphabet;
mod sequence;

pub use nucleotide::{Nucleotide, InvalidNucleotide};
pub use alphabet::Alphabet;
pub use sequence::{Sequence, SharedSequence, InvalidSequence, OutOfBounds};
