//! Commonly used imports for convenience.
//!
//! This prelude module provides a convenient way to import the most commonly
//! used types and traits in the centrevo library.
//!
//! # Example
//!
//! ```
//! use centrevo_sim::prelude::*;
//! use std::str::FromStr;
//!
//! let seq = Sequence::from_str("ACGT").unwrap();
//! ```

pub use crate::errors;
pub use crate::base::fitness::{self, FitnessValue, LogFitnessValue};
pub use crate::base::{Nucleotide, Sequence, SharedSequence};
pub use crate::genome::{Chromosome, Haplotype, Individual};
pub use crate::simulation::Population;

// TODO: Add more re-exports as other modules are implemented
// pub use crate::evolution::{Mutation, Recombination, Selection};
