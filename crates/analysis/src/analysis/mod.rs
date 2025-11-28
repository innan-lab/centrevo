//! Population genetics analysis tools for Centrevo
//!
//! This module provides comprehensive analysis capabilities including:
//! - Diversity metrics (π, Tajima's D, θ_W)
//! - Linkage disequilibrium analysis
//! - Population structure (FST, PCA)
//! - Temporal dynamics
//! - Sequence composition and distance

pub mod diversity;
pub mod linkage;
pub mod distance;
pub mod composition;
pub mod polymorphism;
pub mod structure;
pub mod temporal;
pub mod utils;

// Re-export commonly used functions
pub use diversity::{nucleotide_diversity, tajimas_d, wattersons_theta, haplotype_diversity};
pub use linkage::{linkage_disequilibrium, ld_decay, haplotype_blocks};
pub use distance::{pairwise_distances, distance_matrix};
pub use structure::{fst, pca};
pub use composition::{gc_content, nucleotide_composition};
pub use temporal::{allele_trajectory, fitness_dynamics};
