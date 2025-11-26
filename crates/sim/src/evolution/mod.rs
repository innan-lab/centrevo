//! Evolution module providing mutation, recombination, and selection.
//!
//! This module implements the core evolutionary processes:
//! - **Mutation**: Point mutations using substitution models (JC69, uniform)
//! - **Recombination**: Crossover and gene conversion
//! - **Selection**: Fitness functions (GC content, length, sequence similarity)

pub mod mutation;
pub mod recombination;
pub mod selection;

pub use mutation::SubstitutionModel;
pub use recombination::{RecombinationModel, RecombinationType};
pub use selection::{
    GCContentFitness, HaplotypeFitness, IndividualFitness, LengthFitness, LengthSimilarityFitness,
    SequenceSimilarityFitness,
};
