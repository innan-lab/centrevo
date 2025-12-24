//! Genome structures for representing chromosomes, haplotypes, and individuals.

mod chromosome;
mod haplotype;
mod individual;
pub mod repeat_map;

pub use chromosome::Chromosome;
pub use haplotype::Haplotype;
pub use individual::Individual;
pub use repeat_map::RepeatMap;
