//! Genome structures for representing chromosomes, haplotypes, and individuals.

mod chromosome;
mod haplotype;
mod individual;

pub use chromosome::{Chromosome, SharedChromosome};
pub use haplotype::Haplotype;
pub use individual::Individual;
