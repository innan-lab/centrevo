//! Linkage disequilibrium analysis
//!
//! Measures non-random association of alleles at different loci.

use centrevo_sim::base::Nucleotide;
use centrevo_sim::simulation::Population;

/// Linkage disequilibrium statistics
#[derive(Debug, Clone, Copy)]
pub struct LDStatistics {
    /// D statistic (raw disequilibrium)
    pub d: f64,
    /// D' (normalized D)
    pub d_prime: f64,
    /// r² (correlation coefficient squared)
    pub r_squared: f64,
}

/// Calculate linkage disequilibrium between two positions
///
/// # Arguments
///
/// * `population` - The population to analyze
/// * `pos1` - First position
/// * `pos2` - Second position
/// * `chromosome_idx` - Chromosome index
/// * `haplotype_idx` - Haplotype index
///
/// # Returns
///
/// LD statistics (D, D', r²)
///
/// # References
///
/// Lewontin, R. C. (1964). The Interaction of Selection and Linkage. I. General
/// Considerations; Heterotic Models. Genetics, 49(1), 49-67.
pub fn linkage_disequilibrium(
    population: &Population,
    pos1: usize,
    pos2: usize,
    chromosome_idx: usize,
    haplotype_idx: usize,
) -> Option<LDStatistics> {
    // Extract alleles at both positions
    let alleles: Vec<(Nucleotide, Nucleotide)> = population
        .individuals()
        .iter()
        .filter_map(|ind| {
            let haplotype = if haplotype_idx == 0 {
                ind.haplotype1()
            } else {
                ind.haplotype2()
            };
            haplotype.get(chromosome_idx).and_then(|chr| {
                let seq = chr.sequence();
                match (seq.get(pos1), seq.get(pos2)) {
                    (Some(n1), Some(n2)) => Some((n1, n2)),
                    _ => None,
                }
            })
        })
        .collect();

    if alleles.is_empty() {
        return None;
    }

    let n = alleles.len() as f64;

    // Count haplotype frequencies
    let mut counts: std::collections::HashMap<(Nucleotide, Nucleotide), usize> =
        std::collections::HashMap::new();
    for &allele_pair in &alleles {
        *counts.entry(allele_pair).or_insert(0) += 1;
    }

    // Get the two most common alleles at each position
    let mut pos1_alleles: Vec<(Nucleotide, usize)> = Vec::new();
    let mut pos2_alleles: Vec<(Nucleotide, usize)> = Vec::new();

    for (&(a1, a2), &count) in &counts {
        if let Some(entry) = pos1_alleles.iter_mut().find(|(nuc, _)| *nuc == a1) {
            entry.1 += count;
        } else {
            pos1_alleles.push((a1, count));
        }

        if let Some(entry) = pos2_alleles.iter_mut().find(|(nuc, _)| *nuc == a2) {
            entry.1 += count;
        } else {
            pos2_alleles.push((a2, count));
        }
    }

    if pos1_alleles.len() < 2 || pos2_alleles.len() < 2 {
        // No variation at one or both sites
        return Some(LDStatistics {
            d: 0.0,
            d_prime: 0.0,
            r_squared: 0.0,
        });
    }

    pos1_alleles.sort_by_key(|(_, count)| std::cmp::Reverse(*count));
    pos2_alleles.sort_by_key(|(_, count)| std::cmp::Reverse(*count));

    let allele_a = pos1_alleles[0].0;
    let allele_b = pos2_alleles[0].0;

    // Calculate frequencies
    let p_a = pos1_alleles[0].1 as f64 / n;
    let p_b = pos2_alleles[0].1 as f64 / n;
    let p_ab = *counts.get(&(allele_a, allele_b)).unwrap_or(&0) as f64 / n;

    // Calculate D
    let d = p_ab - p_a * p_b;

    // Calculate D'
    let d_prime = if d >= 0.0 {
        let d_max = (p_a * (1.0 - p_b)).min((1.0 - p_a) * p_b);
        if d_max > 0.0 { d / d_max } else { 0.0 }
    } else {
        let d_max = (p_a * p_b).min((1.0 - p_a) * (1.0 - p_b));
        if d_max > 0.0 { -d / d_max } else { 0.0 }
    };

    // Calculate r²
    let denominator = p_a * (1.0 - p_a) * p_b * (1.0 - p_b);
    let r_squared = if denominator > 0.0 {
        (d * d) / denominator
    } else {
        0.0
    };

    Some(LDStatistics {
        d,
        d_prime,
        r_squared,
    })
}

/// Calculate LD decay with distance
///
/// Computes average r² for pairs of sites at different distances.
///
/// # Arguments
///
/// * `population` - The population to analyze
/// * `chromosome_idx` - Chromosome index
/// * `haplotype_idx` - Haplotype index
/// * `max_distance` - Maximum distance to consider
/// * `bin_size` - Size of distance bins
///
/// # Returns
///
/// Vector of (distance, average_r²) pairs
pub fn ld_decay(
    population: &Population,
    chromosome_idx: usize,
    haplotype_idx: usize,
    max_distance: usize,
    bin_size: usize,
) -> Vec<(usize, f64)> {
    // Get sequence length
    let seq_len = population
        .individuals()
        .iter()
        .filter_map(|ind| {
            let haplotype = if haplotype_idx == 0 {
                ind.haplotype1()
            } else {
                ind.haplotype2()
            };
            haplotype.get(chromosome_idx)
        })
        .next()
        .map(|chr| chr.sequence().len())
        .unwrap_or(0);

    if seq_len == 0 {
        return Vec::new();
    }

    let num_bins = max_distance.div_ceil(bin_size);
    let mut bin_sums = vec![0.0; num_bins];
    let mut bin_counts = vec![0usize; num_bins];

    // Sample pairs of positions
    let sample_step = (seq_len / 100).max(1); // Sample ~100 positions

    for pos1 in (0..seq_len).step_by(sample_step) {
        for pos2 in (pos1 + 1..seq_len.min(pos1 + max_distance)).step_by(sample_step) {
            let distance = pos2 - pos1;
            let bin_idx = (distance / bin_size).min(num_bins - 1);

            if let Some(ld) =
                linkage_disequilibrium(population, pos1, pos2, chromosome_idx, haplotype_idx)
            {
                bin_sums[bin_idx] += ld.r_squared;
                bin_counts[bin_idx] += 1;
            }
        }
    }

    // Calculate averages
    (0..num_bins)
        .filter(|&i| bin_counts[i] > 0)
        .map(|i| {
            let distance = (i * bin_size) + (bin_size / 2);
            let avg_r2 = bin_sums[i] / bin_counts[i] as f64;
            (distance, avg_r2)
        })
        .collect()
}

/// Identify haplotype blocks using simple LD threshold
///
/// Returns list of (start, end) positions for blocks with high LD.
/// This is a simplified version - full implementation would use Gabriel et al. method.
///
/// # Arguments
///
/// * `population` - The population to analyze
/// * `chromosome_idx` - Chromosome index
/// * `haplotype_idx` - Haplotype index
/// * `ld_threshold` - Minimum r² to consider sites in same block (default: 0.8)
///
/// # Returns
///
/// Vector of (start_position, end_position) for each block
pub fn haplotype_blocks(
    population: &Population,
    chromosome_idx: usize,
    haplotype_idx: usize,
    ld_threshold: f64,
) -> Vec<(usize, usize)> {
    let seq_len = population
        .individuals()
        .iter()
        .filter_map(|ind| {
            let haplotype = if haplotype_idx == 0 {
                ind.haplotype1()
            } else {
                ind.haplotype2()
            };
            haplotype.get(chromosome_idx)
        })
        .next()
        .map(|chr| chr.sequence().len())
        .unwrap_or(0);

    if seq_len == 0 {
        return Vec::new();
    }

    let mut blocks = Vec::new();
    let mut block_start = 0;
    let mut in_block = false;

    let sample_step = (seq_len / 50).max(1);

    for pos in (0..seq_len - 1).step_by(sample_step) {
        if let Some(ld) = linkage_disequilibrium(
            population,
            pos,
            pos + sample_step,
            chromosome_idx,
            haplotype_idx,
        ) {
            if ld.r_squared >= ld_threshold {
                if !in_block {
                    block_start = pos;
                    in_block = true;
                }
            } else if in_block {
                blocks.push((block_start, pos));
                in_block = false;
            }
        }
    }

    if in_block {
        blocks.push((block_start, seq_len - 1));
    }

    blocks
}

#[cfg(test)]
mod tests {
    use super::*;
    use centrevo_sim::base::Sequence;
    use centrevo_sim::genome::{Chromosome, Haplotype, Individual};

    fn create_test_individual(id: &str, sequence: &[Nucleotide]) -> Individual {
        let mut seq = Sequence::with_capacity(sequence.len());
        for &nuc in sequence {
            seq.push(nuc);
        }

        // Assume uniform structure for tests
        let ru_len = 10;
        let rus_per_hor = 5;
        let hor_len = ru_len * rus_per_hor;
        let hors_per_chr = if hor_len > 0 { seq.len() / hor_len } else { 0 };
        
        let map = centrevo_sim::genome::repeat_map::RepeatMap::uniform(ru_len, rus_per_hor, hors_per_chr);

        let chr = Chromosome::new(format!("chr_{}", id), seq, map);
        let mut hap1 = Haplotype::new();
        hap1.push(chr);
        let hap2 = Haplotype::new();

        Individual::new(id, hap1, hap2)
    }

    #[test]
    fn test_linkage_disequilibrium_perfect() {
        // Perfect LD: if pos1=A then pos2=T, if pos1=C then pos2=G
        let seq1 = vec![Nucleotide::A, Nucleotide::T];
        let seq2 = vec![Nucleotide::C, Nucleotide::G];

        let individuals = vec![
            create_test_individual("ind1", &seq1),
            create_test_individual("ind2", &seq2),
        ];

        let pop = Population::new("pop1", individuals);
        let ld = linkage_disequilibrium(&pop, 0, 1, 0, 0);

        assert!(ld.is_some());
        let ld = ld.unwrap();
        assert!(ld.r_squared > 0.9); // Should be close to 1.0
    }

    #[test]
    fn test_linkage_disequilibrium_no_variation() {
        // No variation at one site
        let seq1 = vec![Nucleotide::A, Nucleotide::A];
        let seq2 = vec![Nucleotide::A, Nucleotide::T];

        let individuals = vec![
            create_test_individual("ind1", &seq1),
            create_test_individual("ind2", &seq2),
        ];

        let pop = Population::new("pop1", individuals);
        let ld = linkage_disequilibrium(&pop, 0, 1, 0, 0);

        assert!(ld.is_some());
        let ld = ld.unwrap();
        assert_eq!(ld.r_squared, 0.0);
    }

    #[test]
    fn test_linkage_disequilibrium_independent() {
        // Independent sites
        let seq1 = vec![Nucleotide::A, Nucleotide::T];
        let seq2 = vec![Nucleotide::A, Nucleotide::C];
        let seq3 = vec![Nucleotide::C, Nucleotide::T];
        let seq4 = vec![Nucleotide::C, Nucleotide::C];

        let individuals = vec![
            create_test_individual("ind1", &seq1),
            create_test_individual("ind2", &seq2),
            create_test_individual("ind3", &seq3),
            create_test_individual("ind4", &seq4),
        ];

        let pop = Population::new("pop1", individuals);
        let ld = linkage_disequilibrium(&pop, 0, 1, 0, 0);

        assert!(ld.is_some());
        let ld = ld.unwrap();
        // With perfect independence, r² should be close to 0
        assert!(ld.r_squared < 0.5);
    }

    #[test]
    fn test_ld_decay() {
        // Create population with some LD structure
        let seq = vec![Nucleotide::A; 1000];
        let individuals = vec![
            create_test_individual("ind1", &seq),
            create_test_individual("ind2", &seq),
        ];

        let pop = Population::new("pop1", individuals);
        let decay = ld_decay(&pop, 0, 0, 500, 100);

        // Should return some results
        assert!(!decay.is_empty());
    }

    #[test]
    fn test_haplotype_blocks_simple() {
        let seq = vec![Nucleotide::A; 100];
        let individuals = vec![
            create_test_individual("ind1", &seq),
            create_test_individual("ind2", &seq),
        ];

        let pop = Population::new("pop1", individuals);
        let blocks = haplotype_blocks(&pop, 0, 0, 0.8);

        // With no variation, might get one large block or none
        assert!(blocks.len() <= 1);
    }
}
