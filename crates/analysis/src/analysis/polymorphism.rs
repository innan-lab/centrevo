//! Polymorphism analysis
//!
//! Functions for analyzing segregating sites and site frequency spectra.

use centrevo_sim::base::Sequence;
use centrevo_sim::simulation::Population;

/// Calculate site frequency spectrum
///
/// Returns the distribution of allele frequencies across polymorphic sites.
///
/// # Arguments
///
/// * `population` - The population to analyze
/// * `chromosome_idx` - Chromosome index
/// * `haplotype_idx` - Haplotype index
///
/// # Returns
///
/// Vector where index i contains count of sites with i derived alleles
pub fn site_frequency_spectrum(
    _population: &Population,
    _chromosome_idx: usize,
    _haplotype_idx: usize,
) -> Vec<usize> {
    // TODO: Implement in Week 2
    vec![]
}

/// Count segregating sites
///
/// # Arguments
///
/// * `population` - The population to analyze
/// * `chromosome_idx` - Chromosome index
/// * `haplotype_idx` - Haplotype index
///
/// # Returns
///
/// Number of segregating (polymorphic) sites
pub fn count_segregating_sites(
    population: &Population,
    chromosome_idx: usize,
    haplotype_idx: usize,
) -> usize {
    let sequences: Vec<&Sequence> = population
        .individuals()
        .iter()
        .filter_map(|ind| {
            let haplotype = if haplotype_idx == 0 {
                ind.haplotype1()
            } else {
                ind.haplotype2()
            };
            haplotype.get(chromosome_idx).map(|chr| chr.sequence())
        })
        .collect();

    if sequences.is_empty() {
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
    use centrevo_sim::base::{Nucleotide, Sequence};
    use centrevo_sim::genome::{Chromosome, Haplotype, Individual, RepeatMap};

    fn create_test_individual(id: &str, sequence: &[Nucleotide]) -> Individual {
        let mut seq = Sequence::with_capacity(sequence.len());
        for &nuc in sequence {
            seq.push(nuc);
        }

        let total_len = seq.len();
        let ru_len = if total_len > 0 { total_len } else { 10 };
        let rus_per_hor = 1;
        let num_hors = if total_len > 0 { 1 } else { 0 };

        let map = RepeatMap::uniform(ru_len, rus_per_hor, num_hors);

        let chr = Chromosome::new(format!("chr_{id}"), seq, map);
        let mut hap1 = Haplotype::new();
        hap1.push(chr);
        let hap2 = Haplotype::new();

        Individual::new(id, hap1, hap2)
    }

    #[test]
    fn test_count_segregating_sites() {
        let seq1 = vec![Nucleotide::A, Nucleotide::T, Nucleotide::C];
        let seq2 = vec![Nucleotide::A, Nucleotide::G, Nucleotide::C];

        let individuals = vec![
            create_test_individual("ind1", &seq1),
            create_test_individual("ind2", &seq2),
        ];
        let pop = Population::new("pop1", individuals);

        let count = count_segregating_sites(&pop, 0, 0);
        assert_eq!(count, 1); // Only position 1 differs
    }
}
