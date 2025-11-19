//! Distance calculations for sequences
//!
//! Provides functions for calculating genetic distances between sequences.
//!
//! Distance calculations compare all sequences of the same chromosome type
//! across all individuals and both haplotypes (2n sequences total).

use crate::analysis::utils::hamming_distance_fast;
use crate::base::Sequence;
use crate::simulation::Population;
use rayon::prelude::*;

/// Calculate pairwise distances between all sequences
///
/// Compares all sequences of the specified chromosome across all individuals
/// and both haplotypes. For a population of n individuals, this compares 2n sequences.
///
/// # Arguments
///
/// * `population` - The population to analyze
/// * `chromosome_idx` - Chromosome index to analyze
///
/// # Returns
///
/// Vector of pairwise Hamming distances (unnormalized) for all sequence pairs
///
/// # Examples
///
/// ```no_run
/// use centrevo::analysis::distance::pairwise_distances;
/// use centrevo::simulation::Population;
/// # // Note: This example uses a real population with individuals
/// # // let population = ...; // Create your population here
///
/// // Get all pairwise distances for chromosome 0
/// // let distances = pairwise_distances(&population, 0);
/// ```
pub fn pairwise_distances(population: &Population, chromosome_idx: usize) -> Vec<usize> {
    // Collect all sequences from both haplotypes of all individuals
    let sequences: Vec<&Sequence> = population
        .individuals()
        .iter()
        .flat_map(|ind| {
            [
                ind.haplotype1()
                    .get(chromosome_idx)
                    .map(|chr| chr.sequence()),
                ind.haplotype2()
                    .get(chromosome_idx)
                    .map(|chr| chr.sequence()),
            ]
            .into_iter()
            .flatten()
        })
        .collect();

    let n = sequences.len();
    let mut distances = Vec::with_capacity(n * (n - 1) / 2);

    // Use optimized hamming distance calculation
    for i in 0..n {
        for j in (i + 1)..n {
            let dist = hamming_distance_fast(sequences[i], sequences[j]);
            distances.push(dist);
        }
    }

    distances
}

/// Calculate full distance matrix
///
/// Returns an n×n matrix of pairwise distances (normalized by sequence length).
/// Compares all sequences of the specified chromosome across all individuals
/// and both haplotypes.
///
/// # Arguments
///
/// * `population` - The population to analyze
/// * `chromosome_idx` - Chromosome index to analyze
///
/// # Returns
///
/// 2D vector representing symmetric distance matrix where element [i][j]
/// is the normalized distance between sequence i and sequence j
///
/// # Examples
///
/// ```no_run
/// use centrevo::analysis::distance::distance_matrix;
/// use centrevo::simulation::Population;
/// # // Note: This example uses a real population with individuals
/// # // let population = ...; // Create your population here
///
/// // Get distance matrix for chromosome 0
/// // let matrix = distance_matrix(&population, 0);
/// // matrix[i][j] is the normalized distance between sequences i and j
/// ```
pub fn distance_matrix(population: &Population, chromosome_idx: usize) -> Vec<Vec<f64>> {
    // Collect all sequences from both haplotypes of all individuals
    let sequences: Vec<&Sequence> = population
        .individuals()
        .iter()
        .flat_map(|ind| {
            [
                ind.haplotype1()
                    .get(chromosome_idx)
                    .map(|chr| chr.sequence()),
                ind.haplotype2()
                    .get(chromosome_idx)
                    .map(|chr| chr.sequence()),
            ]
            .into_iter()
            .flatten()
        })
        .collect();

    let n = sequences.len();
    if n == 0 {
        return Vec::new();
    }

    let length = sequences[0].len() as f64;
    if length == 0.0 {
        return vec![vec![0.0; n]; n];
    }

    // Compute matrix in parallel by rows
    let matrix: Vec<Vec<f64>> = (0..n)
        .into_par_iter()
        .map(|i| {
            let mut row = vec![0.0; n];
            for j in 0..n {
                if i == j {
                    row[j] = 0.0;
                } else if i < j {
                    let dist = hamming_distance_fast(sequences[i], sequences[j]) as f64 / length;
                    row[j] = dist;
                } else {
                    // Fill in symmetric value (will be set by row j)
                    let dist = hamming_distance_fast(sequences[j], sequences[i]) as f64 / length;
                    row[j] = dist;
                }
            }
            row
        })
        .collect();

    matrix
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::base::Nucleotide;
    use crate::genome::{Chromosome, Haplotype, Individual};

    fn create_test_individual(id: &str, sequence: &[Nucleotide]) -> Individual {
        let mut seq = Sequence::with_capacity(sequence.len());
        for &nuc in sequence {
            seq.push(nuc);
        }

        let chr = Chromosome::new(format!("chr_{}", id), seq, 10, 5);
        let mut hap1 = Haplotype::new();
        hap1.push(chr);
        let hap2 = Haplotype::new();

        Individual::new(id, hap1, hap2)
    }

    #[test]
    fn test_pairwise_distances() {
        let seq1 = vec![Nucleotide::A; 10];
        let seq2 = vec![Nucleotide::T; 10];
        let seq3 = vec![Nucleotide::C; 10];

        let individuals = vec![
            create_test_individual("ind1", &seq1),
            create_test_individual("ind2", &seq2),
            create_test_individual("ind3", &seq3),
        ];

        let pop = Population::new("pop1", individuals);
        let distances = pairwise_distances(&pop, 0);

        // 3 individuals = 3 sequences (only haplotype1 has data)
        // = 3 pairs: (0,1), (0,2), (1,2)
        assert_eq!(distances.len(), 3);
        // All completely different
        assert!(distances.iter().all(|&d| d == 10));
    }

    #[test]
    fn test_pairwise_distances_both_haplotypes() {
        // Create individuals with data in both haplotypes
        let mut seq1 = Sequence::with_capacity(5);
        let mut seq2 = Sequence::with_capacity(5);

        for _ in 0..5 {
            seq1.push(Nucleotide::A);
            seq2.push(Nucleotide::T);
        }

        let chr1 = Chromosome::new("chr1", seq1, 5, 1);
        let chr2 = Chromosome::new("chr2", seq2, 5, 1);

        let mut hap1 = Haplotype::new();
        hap1.push(chr1);
        let mut hap2 = Haplotype::new();
        hap2.push(chr2);

        let ind = Individual::new("ind1", hap1, hap2);
        let pop = Population::new("pop1", vec![ind]);

        let distances = pairwise_distances(&pop, 0);

        // 1 individual with 2 haplotypes = 2 sequences = 1 pair
        assert_eq!(distances.len(), 1);
        // Both sequences completely different
        assert_eq!(distances[0], 5);
    }

    #[test]
    fn test_distance_matrix() {
        let seq1 = vec![Nucleotide::A; 10];
        let mut seq2 = vec![Nucleotide::A; 10];
        seq2[0] = Nucleotide::T; // 1 difference

        let individuals = vec![
            create_test_individual("ind1", &seq1),
            create_test_individual("ind2", &seq2),
        ];

        let pop = Population::new("pop1", individuals);
        let matrix = distance_matrix(&pop, 0);

        // 2 individuals with only haplotype1 = 2 sequences
        assert_eq!(matrix.len(), 2);
        assert_eq!(matrix[0].len(), 2);
        assert_eq!(matrix[0][0], 0.0); // Self distance
        assert!((matrix[0][1] - 0.1).abs() < 1e-10); // 1/10 = 0.1
        assert!((matrix[1][0] - 0.1).abs() < 1e-10); // Symmetric
        assert_eq!(matrix[1][1], 0.0);
    }

    #[test]
    fn test_distance_matrix_both_haplotypes() {
        // Create two individuals, each with different sequences in both haplotypes
        let mut seq1 = Sequence::with_capacity(4);
        let mut seq2 = Sequence::with_capacity(4);
        let mut seq3 = Sequence::with_capacity(4);
        let mut seq4 = Sequence::with_capacity(4);

        for _ in 0..4 {
            seq1.push(Nucleotide::A); // All A
            seq2.push(Nucleotide::T); // All T
            seq3.push(Nucleotide::C); // All C
            seq4.push(Nucleotide::G); // All G
        }

        let chr1 = Chromosome::new("chr1", seq1, 4, 1);
        let chr2 = Chromosome::new("chr2", seq2, 4, 1);
        let chr3 = Chromosome::new("chr3", seq3, 4, 1);
        let chr4 = Chromosome::new("chr4", seq4, 4, 1);

        let mut hap1_ind1 = Haplotype::new();
        hap1_ind1.push(chr1);
        let mut hap2_ind1 = Haplotype::new();
        hap2_ind1.push(chr2);

        let mut hap1_ind2 = Haplotype::new();
        hap1_ind2.push(chr3);
        let mut hap2_ind2 = Haplotype::new();
        hap2_ind2.push(chr4);

        let ind1 = Individual::new("ind1", hap1_ind1, hap2_ind1);
        let ind2 = Individual::new("ind2", hap1_ind2, hap2_ind2);

        let pop = Population::new("pop1", vec![ind1, ind2]);
        let matrix = distance_matrix(&pop, 0);

        // 2 individuals × 2 haplotypes = 4 sequences
        assert_eq!(matrix.len(), 4);
        assert_eq!(matrix[0].len(), 4);

        // All different sequences should have distance = 1.0
        for (i, row) in matrix.iter().enumerate() {
            for (j, &dist) in row.iter().enumerate() {
                if i == j {
                    assert_eq!(dist, 0.0);
                } else {
                    assert_eq!(dist, 1.0); // All 4 bases different
                }
            }
        }
    }
}
