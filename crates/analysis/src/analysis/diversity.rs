//! Diversity metrics for population genetics
//!
//! Implements standard measures of genetic diversity including:
//! - Nucleotide diversity (π)
//! - Tajima's D
//! - Watterson's estimator (θ_W)
//! - Haplotype diversity

use crate::analysis::utils::{hamming_distance_fast, harmonic_number};
use centrevo_sim::base::{GenomeArena, Nucleotide};
use centrevo_sim::simulation::Population;
use rayon::prelude::*;
use std::collections::HashMap;

/// Calculate nucleotide diversity (π) for a population
///
/// Nucleotide diversity is the average number of nucleotide differences
/// per site between two sequences.
///
/// # Formula
///
/// $$\pi = \frac{\sum_{i<j} d_{ij}}{n(n-1)/2 \cdot L}$$
///
/// where $d_{ij}$ is the number of differences between sequences $i$ and $j$,
/// $n$ is the number of sequences, and $L$ is the sequence length.
///
/// # Arguments
///
/// * `population` - The population to analyze
/// * `chromosome_idx` - Index of chromosome to analyze (default: 0)
///
/// # Returns
///
/// Average pairwise nucleotide diversity per site across all 2n sequences
///
/// # Examples
///
/// ```
/// use centrevo::analysis::diversity::nucleotide_diversity;
/// use centrevo::simulation::Population;
///
/// // Assuming you have a population
/// // let population = Population::new("pop1", individuals);
/// // let pi = nucleotide_diversity(&population, 0);
/// // println!("Nucleotide diversity: {:.6}", pi);
/// ```
///
/// # References
///
/// Nei, M., & Li, W. H. (1979). Mathematical model for studying genetic
/// variation in terms of restriction endonucleases. PNAS, 76(10), 5269-5273.
pub fn nucleotide_diversity(
    population: &Population,
    chromosome_idx: usize,
    arena: &GenomeArena,
) -> f64 {
    let sequences: Vec<&[Nucleotide]> = population
        .individuals()
        .iter()
        .flat_map(|ind| {
            [ind.haplotype1(), ind.haplotype2()]
                .into_iter()
                .filter_map(|hap| hap.get(chromosome_idx).map(|chr| chr.sequence(arena)))
        })
        .collect();

    let n = sequences.len();
    if n < 2 {
        return 0.0;
    }

    let length = sequences[0].len();
    if length == 0 {
        return 0.0;
    }

    // Calculate pairwise differences in parallel
    let total_differences: usize = (0..n)
        .into_par_iter()
        .map(|i| {
            (i + 1..n)
                .map(|j| hamming_distance_fast(sequences[i], sequences[j]))
                .sum::<usize>()
        })
        .sum();

    let num_comparisons = n * (n - 1) / 2;
    total_differences as f64 / (num_comparisons * length) as f64
}

/// Calculate Tajima's D statistic
///
/// Tajima's D tests the hypothesis of neutral evolution by comparing
/// two estimates of θ: one based on the number of segregating sites
/// and one based on nucleotide diversity.
///
/// # Formula
///
/// $$D = \frac{\pi - \theta_W}{\sqrt{Var(\pi - \theta_W)}}$$
///
/// Positive D suggests balancing selection or population contraction.
/// Negative D suggests purifying selection or population expansion.
///
/// # Arguments
///
/// * `population` - The population to analyze
/// * `chromosome_idx` - Index of chromosome to analyze
///
/// # Returns
///
/// Tajima's D statistic across all 2n sequences. Returns 0.0 if variance is zero or calculation is undefined.
///
/// # References
///
/// Tajima, F. (1989). Statistical method for testing the neutral mutation
/// hypothesis by DNA polymorphism. Genetics, 123(3), 585-595.
pub fn tajimas_d(population: &Population, chromosome_idx: usize, arena: &GenomeArena) -> f64 {
    let sequences: Vec<&[Nucleotide]> = population
        .individuals()
        .iter()
        .flat_map(|ind| {
            [ind.haplotype1(), ind.haplotype2()]
                .into_iter()
                .filter_map(|hap| hap.get(chromosome_idx).map(|chr| chr.sequence(arena)))
        })
        .collect();

    let n = sequences.len();
    if n < 2 {
        return 0.0;
    }

    let length = sequences[0].len();
    if length == 0 {
        return 0.0;
    }

    // Calculate π (nucleotide diversity)
    let pi = nucleotide_diversity(population, chromosome_idx, arena);

    // Calculate θ_W (Watterson's estimator)
    let theta_w = wattersons_theta(population, chromosome_idx, arena);

    // If both are zero, no variation
    if pi == 0.0 && theta_w == 0.0 {
        return 0.0;
    }

    // Calculate variance components
    let n_f64 = n as f64;
    let a1 = harmonic_number(n);
    let a2: f64 = (1..n).map(|i| 1.0 / (i * i) as f64).sum();

    let b1 = (n_f64 + 1.0) / (3.0 * (n_f64 - 1.0));
    let b2 = 2.0 * (n_f64 * n_f64 + n_f64 + 3.0) / (9.0 * n_f64 * (n_f64 - 1.0));

    let c1 = b1 - 1.0 / a1;
    let c2 = b2 - (n_f64 + 2.0) / (a1 * n_f64) + a2 / (a1 * a1);

    let e1 = c1 / a1;
    let e2 = c2 / (a1 * a1 + a2);

    // Number of segregating sites
    let s = segregating_sites(&sequences) as f64;

    if s == 0.0 {
        return 0.0;
    }

    // Variance of (π - θ_W)
    let var = e1 * s + e2 * s * (s - 1.0);

    if var <= 0.0 {
        return 0.0;
    }

    (pi * length as f64 - theta_w * length as f64) / var.sqrt()
}

/// Calculate Watterson's estimator (θ_W)
///
/// Estimates θ = 4Nμ from the number of segregating sites.
///
/// # Formula
///
/// $$\theta_W = \frac{S}{a_n}$$
///
/// where $S$ is the number of segregating sites and
/// $a_n = \sum_{i=1}^{n-1} \frac{1}{i}$
///
/// # Arguments
///
/// * `population` - The population to analyze
/// * `chromosome_idx` - Index of chromosome to analyze
///
/// # Returns
///
/// Watterson's theta per site across all 2n sequences
///
/// # References
///
/// Watterson, G. A. (1975). On the number of segregating sites in genetical
/// models without recombination. Theoretical Population Biology, 7(2), 256-276.
pub fn wattersons_theta(
    population: &Population,
    chromosome_idx: usize,
    arena: &GenomeArena,
) -> f64 {
    let sequences: Vec<&[Nucleotide]> = population
        .individuals()
        .iter()
        .flat_map(|ind| {
            [ind.haplotype1(), ind.haplotype2()]
                .into_iter()
                .filter_map(|hap| hap.get(chromosome_idx).map(|chr| chr.sequence(arena)))
        })
        .collect();

    let n = sequences.len();
    if n < 2 {
        return 0.0;
    }

    let length = sequences[0].len();
    if length == 0 {
        return 0.0;
    }

    let s = segregating_sites(&sequences) as f64;
    let a_n = harmonic_number(n);

    s / (a_n * length as f64)
}

/// Calculate haplotype diversity
///
/// Returns the probability that two randomly chosen haplotypes
/// are different. This is equivalent to expected heterozygosity
/// for haploid data.
///
/// # Formula
///
/// $$H = 1 - \sum_{i} p_i^2$$
///
/// where $p_i$ is the frequency of haplotype $i$
///
/// # Arguments
///
/// * `population` - The population to analyze
/// * `chromosome_idx` - Index of chromosome to analyze
///
/// # Returns
///
/// Haplotype diversity (0.0 to 1.0) across all 2n sequences
pub fn haplotype_diversity(
    population: &Population,
    chromosome_idx: usize,
    arena: &GenomeArena,
) -> f64 {
    let sequences: Vec<String> = population
        .individuals()
        .iter()
        .flat_map(|ind| {
            [ind.haplotype1(), ind.haplotype2()]
                .into_iter()
                .filter_map(|hap| {
                    hap.get(chromosome_idx).map(|chr| {
                        // Convert sequence to string for hashing
                        let seq = chr.sequence(arena);
                        (0..seq.len())
                            .filter_map(|i| seq.get(i))
                            .map(|n| format!("{n:?}"))
                            .collect::<String>()
                    })
                })
        })
        .collect();

    let n = sequences.len();
    if n == 0 {
        return 0.0;
    }

    // Count haplotype frequencies
    let mut counts: HashMap<String, usize> = HashMap::new();
    for seq in sequences {
        *counts.entry(seq).or_insert(0) += 1;
    }

    // Calculate diversity
    let sum_squared_freqs: f64 = counts
        .values()
        .map(|&count| {
            let freq = count as f64 / n as f64;
            freq * freq
        })
        .sum();

    1.0 - sum_squared_freqs
}

// ===== Helper Functions =====

/// Calculate number of segregating sites
fn segregating_sites(sequences: &[&[Nucleotide]]) -> usize {
    let n_seq = sequences.len();
    if n_seq == 0 {
        return 0;
    }

    let length = sequences[0].len();
    (0..length)
        .filter(|&pos| {
            let first = sequences[0].get(pos);
            sequences.iter().any(|seq| seq.get(pos) != first)
        })
        .count()
}

#[cfg(test)]
mod tests {
    use super::*;
    use centrevo_sim::base::Nucleotide;
    use centrevo_sim::genome::{Chromosome, Haplotype, Individual};

    fn create_test_individual(
        id: &str,
        sequence: &[Nucleotide],
        arena: &mut GenomeArena,
    ) -> Individual {
        let mut seq = centrevo_sim::base::Sequence::with_capacity(sequence.len());
        for &nuc in sequence {
            seq.push(nuc);
        }
        let data = arena.alloc(seq.as_slice());

        // Assume uniform structure for tests
        let ru_len = 10;
        let rus_per_hor = 5;
        let hor_len = ru_len * rus_per_hor;
        let hors_per_chr = if hor_len > 0 { seq.len() / hor_len } else { 0 };

        let map =
            centrevo_sim::genome::repeat_map::RepeatMap::uniform(ru_len, rus_per_hor, hors_per_chr);

        let chr = Chromosome::new(format!("chr_{id}"), data, map);
        let mut hap1 = Haplotype::new();
        hap1.push(chr.clone());
        let mut hap2 = Haplotype::new();
        hap2.push(chr);

        Individual::new(id, hap1, hap2)
    }

    #[test]
    fn test_nucleotide_diversity_identical_sequences() {
        let mut arena = GenomeArena::new();
        // All identical sequences should have π = 0
        let sequence = vec![Nucleotide::A; 100];
        let individuals = vec![
            create_test_individual("ind1", &sequence, &mut arena),
            create_test_individual("ind2", &sequence, &mut arena),
            create_test_individual("ind3", &sequence, &mut arena),
        ];

        let pop = Population::new("pop1", individuals);
        let pi = nucleotide_diversity(&pop, 0, &arena);

        assert_eq!(pi, 0.0);
    }

    #[test]
    fn test_nucleotide_diversity_completely_different() {
        let mut arena = GenomeArena::new();
        // Two individuals with completely different sequences
        // Total: 4 sequences (2 × 2 haplotypes)
        // Since both haplotypes are identical within each individual,
        // we have: 2× seq1 (all A) and 2× seq2 (all T)
        let seq1 = vec![Nucleotide::A; 100];
        let seq2 = vec![Nucleotide::T; 100];

        let individuals = vec![
            create_test_individual("ind1", &seq1, &mut arena),
            create_test_individual("ind2", &seq2, &mut arena),
        ];

        let pop = Population::new("pop1", individuals);
        let pi = nucleotide_diversity(&pop, 0, &arena);

        // 6 pairwise comparisons:
        // - 1 comparison within seq1 (0 differences)
        // - 1 comparison within seq2 (0 differences)
        // - 4 comparisons between seq1 and seq2 (100 differences each)
        // Total: 400 differences across 6 comparisons × 100 sites = 600
        // π = 400/600 = 2/3 ≈ 0.667
        assert!((pi - 2.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_nucleotide_diversity_known_value() {
        let mut arena = GenomeArena::new();
        // Simple test case: 3 individuals with 6 total sequences (3 × 2 haplotypes)
        // Each individual has identical haplotypes for simplicity
        // Seq1: AAAA (ind1, both haplotypes)
        // Seq2: AAAT (ind2, both haplotypes)
        // Seq3: AATT (ind3, both haplotypes)
        // With 6 sequences total: pairs from same haplotype contribute same as original
        // 15 pairwise comparisons total, differences scale accordingly

        let seq1 = vec![Nucleotide::A, Nucleotide::A, Nucleotide::A, Nucleotide::A];
        let seq2 = vec![Nucleotide::A, Nucleotide::A, Nucleotide::A, Nucleotide::T];
        let seq3 = vec![Nucleotide::A, Nucleotide::A, Nucleotide::T, Nucleotide::T];

        let individuals = vec![
            create_test_individual("ind1", &seq1, &mut arena),
            create_test_individual("ind2", &seq2, &mut arena),
            create_test_individual("ind3", &seq3, &mut arena),
        ];

        let pop = Population::new("pop1", individuals);
        let pi = nucleotide_diversity(&pop, 0, &arena);

        // With 6 identical pairs (doubled from 3 individuals):
        // Comparisons: seq1-seq2 (×4 ways) = 4 diff, seq1-seq3 (×4 ways) = 8 diff, seq2-seq3 (×4 ways) = 4 diff
        // Total: 4 + 8 + 4 = 16 differences across 15 comparisons × 4 sites = 60
        // π = 16/60 = 4/15 ≈ 0.267
        assert!((pi - 4.0 / 15.0).abs() < 1e-10);
    }

    #[test]
    fn test_nucleotide_diversity_empty_population() {
        let arena = GenomeArena::new();
        let pop = Population::new("pop1", Vec::new());
        let pi = nucleotide_diversity(&pop, 0, &arena);
        assert_eq!(pi, 0.0);
    }

    #[test]
    fn test_nucleotide_diversity_single_individual() {
        let mut arena = GenomeArena::new();
        let sequence = vec![Nucleotide::A; 100];
        let individuals = vec![create_test_individual("ind1", &sequence, &mut arena)];

        let pop = Population::new("pop1", individuals);
        let pi = nucleotide_diversity(&pop, 0, &arena);

        // Single individual has 2 identical haplotypes (by construction)
        // π between identical sequences = 0
        assert_eq!(pi, 0.0);
    }

    #[test]
    fn test_segregating_sites_no_variation() {
        let mut seq = Vec::with_capacity(10);
        for _ in 0..10 {
            seq.push(Nucleotide::A);
        }

        let slice = seq.as_slice();
        let sequences = vec![slice, slice, slice];
        let s = segregating_sites(&sequences);

        assert_eq!(s, 0);
    }

    #[test]
    fn test_segregating_sites_all_different() {
        let seq1 = vec![Nucleotide::A; 4];
        let seq2 = vec![Nucleotide::T; 4];

        let slice1 = seq1.as_slice();
        let slice2 = seq2.as_slice();
        let sequences = vec![slice1, slice2];
        let s = segregating_sites(&sequences);

        assert_eq!(s, 4);
    }

    #[test]
    fn test_harmonic_number() {
        assert_eq!(harmonic_number(1), 0.0);
        assert!((harmonic_number(2) - 1.0).abs() < 1e-10);
        assert!((harmonic_number(3) - 1.5).abs() < 1e-10);
        assert!((harmonic_number(4) - (1.0 + 0.5 + 1.0 / 3.0)).abs() < 1e-10);
    }

    #[test]
    fn test_wattersons_theta_no_variation() {
        let mut arena = GenomeArena::new();
        let sequence = vec![Nucleotide::A; 100];
        let individuals = vec![
            create_test_individual("ind1", &sequence, &mut arena),
            create_test_individual("ind2", &sequence, &mut arena),
        ];

        let pop = Population::new("pop1", individuals);
        let theta = wattersons_theta(&pop, 0, &arena);

        assert_eq!(theta, 0.0);
    }

    #[test]
    fn test_wattersons_theta_with_variation() {
        let mut arena = GenomeArena::new();
        // Two individuals, 4 total sequences
        let seq1 = vec![Nucleotide::A; 100];
        let mut seq2 = vec![Nucleotide::A; 100];
        seq2[50] = Nucleotide::T;

        let individuals = vec![
            create_test_individual("ind1", &seq1, &mut arena),
            create_test_individual("ind2", &seq2, &mut arena),
        ];

        let pop = Population::new("pop1", individuals);
        let theta = wattersons_theta(&pop, 0, &arena);

        // With 4 sequences: a_n = 1 + 1/2 + 1/3 ≈ 1.833
        // S = 1 segregating site
        // θ_W = S / (a_n * L) = 1 / (1.833 * 100) ≈ 0.00545
        let a_4 = 1.0 + 0.5 + 1.0 / 3.0;
        assert!((theta - 1.0 / (a_4 * 100.0)).abs() < 1e-10);
    }

    #[test]
    fn test_haplotype_diversity_all_unique() {
        let mut arena = GenomeArena::new();
        let seq1 = vec![Nucleotide::A; 10];
        let seq2 = vec![Nucleotide::T; 10];
        let seq3 = vec![Nucleotide::C; 10];
        let seq4 = vec![Nucleotide::G; 10];

        let individuals = vec![
            create_test_individual("ind1", &seq1, &mut arena),
            create_test_individual("ind2", &seq2, &mut arena),
            create_test_individual("ind3", &seq3, &mut arena),
            create_test_individual("ind4", &seq4, &mut arena),
        ];

        let pop = Population::new("pop1", individuals);
        let h = haplotype_diversity(&pop, 0, &arena);

        // 4 individuals × 2 haplotypes = 8 total sequences
        // Each individual has both haplotypes identical, so we have:
        // 2× seq1, 2× seq2, 2× seq3, 2× seq4
        // H = 1 - Σ(p_i^2) = 1 - 4×(2/8)^2 = 1 - 4×(1/16) = 1 - 0.25 = 0.75
        assert!((h - 0.75).abs() < 1e-10);
    }

    #[test]
    fn test_haplotype_diversity_all_identical() {
        let mut arena = GenomeArena::new();
        let sequence = vec![Nucleotide::A; 10];
        let individuals = vec![
            create_test_individual("ind1", &sequence, &mut arena),
            create_test_individual("ind2", &sequence, &mut arena),
            create_test_individual("ind3", &sequence, &mut arena),
        ];

        let pop = Population::new("pop1", individuals);
        let h = haplotype_diversity(&pop, 0, &arena);

        // All 6 sequences (3 individuals × 2 haplotypes) are identical
        // H = 1 - 1 = 0
        assert_eq!(h, 0.0);
    }

    #[test]
    fn test_tajimas_d_neutral() {
        let mut arena = GenomeArena::new();
        // For neutral population, expect D ≈ 0
        // This is a simple test with limited data
        let mut seq1 = vec![Nucleotide::A; 100];
        let mut seq2 = vec![Nucleotide::A; 100];
        seq1[10] = Nucleotide::T;
        seq2[20] = Nucleotide::C;

        let individuals = vec![
            create_test_individual("ind1", &seq1, &mut arena),
            create_test_individual("ind2", &seq2, &mut arena),
        ];

        let pop = Population::new("pop1", individuals);
        let d = tajimas_d(&pop, 0, &arena);

        // With 4 sequences (2 individuals × 2 haplotypes), the result should be defined
        // The exact value depends on the implementation details
        assert!(d.is_finite());
    }

    #[test]
    fn test_hamming_distance() {
        let mut seq1 = Vec::with_capacity(10);
        let mut seq2 = Vec::with_capacity(10);

        for i in 0..10 {
            seq1.push(if i % 2 == 0 {
                Nucleotide::A
            } else {
                Nucleotide::T
            });
            seq2.push(if i % 2 == 0 {
                Nucleotide::A
            } else {
                Nucleotide::C
            });
        }

        let dist = hamming_distance_fast(&seq1, &seq2);
        assert_eq!(dist, 5); // Half are different
    }
}
