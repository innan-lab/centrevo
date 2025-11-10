//! Commonly used imports for convenience.
//!
//! This prelude module provides a convenient way to import the most commonly
//! used types and traits in the centrevo library.
//!
//! # Example
//!
//! ```
//! use centrevo::prelude::*;
//! 
//! let alphabet = Alphabet::dna();
//! let seq = Sequence::from_str("ACGT", alphabet).unwrap();
//! ```

pub use crate::base::{Alphabet, InvalidNucleotide, Nucleotide, Sequence, SharedSequence};

// TODO: Add more re-exports as other modules are implemented
// pub use crate::genome::{Chromosome, Haplotype, Individual};
// pub use crate::evolution::{Mutation, Recombination, Selection};
