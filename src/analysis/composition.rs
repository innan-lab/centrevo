//! Sequence composition analysis
//! 
//! Functions for analyzing nucleotide composition and content.
//! 
//! These functions are flexible and work at multiple levels:
//! - Population level: No indices specified
//! - Individual level: Only individual_idx specified
//! - Haplotype level: individual_idx and haplotype_idx specified
//! - Chromosome level: All indices specified

use crate::base::Nucleotide;
use crate::simulation::Population;
use std::collections::HashMap;

/// Calculate GC content flexibly based on provided indices
/// 
/// This function adapts to different levels of analysis:
/// - **Population level**: `gc_content(pop, None, None, None)` - Mean across all sequences
/// - **Individual level**: `gc_content(pop, Some(i), None, None)` - Mean across both haplotypes
/// - **Haplotype level**: `gc_content(pop, Some(i), Some(h), None)` - Mean across all chromosomes in haplotype
/// - **Chromosome level**: `gc_content(pop, Some(i), Some(h), Some(c))` - Single chromosome
/// 
/// # Arguments
/// 
/// * `population` - The population to analyze
/// * `individual_idx` - Optional index of individual
/// * `haplotype_idx` - Optional haplotype index (0 or 1)
/// * `chromosome_idx` - Optional chromosome index
/// 
/// # Returns
/// 
/// GC content as proportion (0.0 to 1.0)
/// 
/// # Examples
/// 
/// ```
/// use centrevo::analysis::composition::gc_content;
/// # use centrevo::simulation::Population;
/// # let population = Population::new("test", vec![]);
/// 
/// // Population-wide GC content
/// let pop_gc = gc_content(&population, None, None, None);
/// 
/// // Individual-level GC content (both haplotypes)
/// let ind_gc = gc_content(&population, Some(0), None, None);
/// 
/// // Haplotype-level GC content (all chromosomes in haplotype)
/// let hap_gc = gc_content(&population, Some(0), Some(0), None);
/// 
/// // Chromosome-level GC content
/// let chr_gc = gc_content(&population, Some(0), Some(0), Some(0));
/// ```
pub fn gc_content(
    population: &Population,
    individual_idx: Option<usize>,
    haplotype_idx: Option<usize>,
    chromosome_idx: Option<usize>,
) -> f64 {
    match (individual_idx, haplotype_idx, chromosome_idx) {
        // Chromosome level: specific chromosome
        (Some(ind_idx), Some(hap_idx), Some(chr_idx)) => {
            gc_content_chromosome(population, ind_idx, hap_idx, chr_idx)
        }
        
        // Haplotype level: all chromosomes in haplotype
        (Some(ind_idx), Some(hap_idx), None) => {
            gc_content_haplotype(population, ind_idx, hap_idx)
        }
        
        // Individual level: both haplotypes
        (Some(ind_idx), None, None) => {
            gc_content_individual(population, ind_idx)
        }
        
        // Population level: all individuals, haplotypes, chromosomes
        (None, None, None) => {
            gc_content_population(population)
        }
        
        // Invalid combinations
        _ => {
            eprintln!("Warning: Invalid index combination for gc_content. Use None for all higher levels.");
            0.0
        }
    }
}

/// Calculate nucleotide composition flexibly based on provided indices
/// 
/// This function adapts to different levels of analysis:
/// - **Population level**: Sum across all sequences
/// - **Individual level**: Sum across both haplotypes
/// - **Haplotype level**: Sum across all chromosomes in haplotype
/// - **Chromosome level**: Single chromosome counts
/// 
/// # Arguments
/// 
/// * `population` - The population to analyze
/// * `individual_idx` - Optional index of individual
/// * `haplotype_idx` - Optional haplotype index (0 or 1)
/// * `chromosome_idx` - Optional chromosome index
/// 
/// # Returns
/// 
/// HashMap with nucleotide counts
/// 
/// # Examples
/// 
/// ```
/// use centrevo::analysis::composition::nucleotide_composition;
/// # use centrevo::simulation::Population;
/// # let population = Population::new("test", vec![]);
/// 
/// // Population-wide composition
/// let pop_comp = nucleotide_composition(&population, None, None, None);
/// 
/// // Individual-level composition
/// let ind_comp = nucleotide_composition(&population, Some(0), None, None);
/// 
/// // Chromosome-level composition
/// let chr_comp = nucleotide_composition(&population, Some(0), Some(0), Some(0));
/// ```
pub fn nucleotide_composition(
    population: &Population,
    individual_idx: Option<usize>,
    haplotype_idx: Option<usize>,
    chromosome_idx: Option<usize>,
) -> HashMap<Nucleotide, usize> {
    match (individual_idx, haplotype_idx, chromosome_idx) {
        // Chromosome level: specific chromosome
        (Some(ind_idx), Some(hap_idx), Some(chr_idx)) => {
            composition_chromosome(population, ind_idx, hap_idx, chr_idx)
        }
        
        // Haplotype level: all chromosomes in haplotype
        (Some(ind_idx), Some(hap_idx), None) => {
            composition_haplotype(population, ind_idx, hap_idx)
        }
        
        // Individual level: both haplotypes
        (Some(ind_idx), None, None) => {
            composition_individual(population, ind_idx)
        }
        
        // Population level: all individuals, haplotypes, chromosomes
        (None, None, None) => {
            composition_population(population)
        }
        
        // Invalid combinations
        _ => {
            eprintln!("Warning: Invalid index combination for nucleotide_composition. Use None for all higher levels.");
            HashMap::new()
        }
    }
}

// ===== Helper functions for GC content at different levels =====

fn gc_content_chromosome(
    population: &Population,
    individual_idx: usize,
    haplotype_idx: usize,
    chromosome_idx: usize,
) -> f64 {
    if let Some(ind) = population.individuals().get(individual_idx) {
        let haplotype = if haplotype_idx == 0 {
            ind.haplotype1()
        } else {
            ind.haplotype2()
        };

        if let Some(chr) = haplotype.get(chromosome_idx) {
            let seq = chr.sequence();
            let total = seq.len();
            if total == 0 {
                return 0.0;
            }

            let gc_count = seq
                .indices()
                .iter()
                .filter_map(|&idx| crate::base::Nucleotide::from_index(idx))
                .filter(|&n| n == Nucleotide::G || n == Nucleotide::C)
                .count();

            return gc_count as f64 / total as f64;
        }
    }
    0.0
}

fn gc_content_haplotype(
    population: &Population,
    individual_idx: usize,
    haplotype_idx: usize,
) -> f64 {
    if let Some(ind) = population.individuals().get(individual_idx) {
        let haplotype = if haplotype_idx == 0 {
            ind.haplotype1()
        } else {
            ind.haplotype2()
        };

        let mut total_gc = 0;
        let mut total_bases = 0;

        for chr in haplotype.chromosomes() {
            let seq = chr.sequence();
            total_bases += seq.len();
            total_gc += seq
                .indices()
                .iter()
                .filter_map(|&idx| crate::base::Nucleotide::from_index(idx))
                .filter(|&n| n == Nucleotide::G || n == Nucleotide::C)
                .count();
        }

        if total_bases > 0 {
            return total_gc as f64 / total_bases as f64;
        }
    }
    0.0
}

fn gc_content_individual(population: &Population, individual_idx: usize) -> f64 {
    if let Some(ind) = population.individuals().get(individual_idx) {
        let mut total_gc = 0;
        let mut total_bases = 0;

        // Process both haplotypes
        for haplotype in [ind.haplotype1(), ind.haplotype2()] {
            for chr in haplotype.chromosomes() {
                let seq = chr.sequence();
                total_bases += seq.len();
                total_gc += seq
                    .indices()
                    .iter()
                    .filter_map(|&idx| crate::base::Nucleotide::from_index(idx))
                    .filter(|&n| n == Nucleotide::G || n == Nucleotide::C)
                    .count();
            }
        }

        if total_bases > 0 {
            return total_gc as f64 / total_bases as f64;
        }
    }
    0.0
}

fn gc_content_population(population: &Population) -> f64 {
    let mut total_gc = 0;
    let mut total_bases = 0;

    for ind in population.individuals() {
        for haplotype in [ind.haplotype1(), ind.haplotype2()] {
            for chr in haplotype.chromosomes() {
                let seq = chr.sequence();
                total_bases += seq.len();
                total_gc += seq
                    .indices()
                    .iter()
                    .filter_map(|&idx| crate::base::Nucleotide::from_index(idx))
                    .filter(|&n| n == Nucleotide::G || n == Nucleotide::C)
                    .count();
            }
        }
    }

    if total_bases > 0 {
        total_gc as f64 / total_bases as f64
    } else {
        0.0
    }
}

// ===== Helper functions for nucleotide composition at different levels =====

fn composition_chromosome(
    population: &Population,
    individual_idx: usize,
    haplotype_idx: usize,
    chromosome_idx: usize,
) -> HashMap<Nucleotide, usize> {
    let mut counts = HashMap::new();
    counts.insert(Nucleotide::A, 0);
    counts.insert(Nucleotide::C, 0);
    counts.insert(Nucleotide::G, 0);
    counts.insert(Nucleotide::T, 0);

    if let Some(ind) = population.individuals().get(individual_idx) {
        let haplotype = if haplotype_idx == 0 {
            ind.haplotype1()
        } else {
            ind.haplotype2()
        };

        if let Some(chr) = haplotype.get(chromosome_idx) {
            let seq = chr.sequence();
            for i in 0..seq.len() {
                if let Some(nuc) = seq.get(i) {
                    *counts.entry(nuc).or_insert(0) += 1;
                }
            }
        }
    }

    counts
}

fn composition_haplotype(
    population: &Population,
    individual_idx: usize,
    haplotype_idx: usize,
) -> HashMap<Nucleotide, usize> {
    let mut counts = HashMap::new();
    counts.insert(Nucleotide::A, 0);
    counts.insert(Nucleotide::C, 0);
    counts.insert(Nucleotide::G, 0);
    counts.insert(Nucleotide::T, 0);

    if let Some(ind) = population.individuals().get(individual_idx) {
        let haplotype = if haplotype_idx == 0 {
            ind.haplotype1()
        } else {
            ind.haplotype2()
        };

        for chr in haplotype.chromosomes() {
            let seq = chr.sequence();
            for i in 0..seq.len() {
                if let Some(nuc) = seq.get(i) {
                    *counts.entry(nuc).or_insert(0) += 1;
                }
            }
        }
    }

    counts
}

fn composition_individual(
    population: &Population,
    individual_idx: usize,
) -> HashMap<Nucleotide, usize> {
    let mut counts = HashMap::new();
    counts.insert(Nucleotide::A, 0);
    counts.insert(Nucleotide::C, 0);
    counts.insert(Nucleotide::G, 0);
    counts.insert(Nucleotide::T, 0);

    if let Some(ind) = population.individuals().get(individual_idx) {
        for haplotype in [ind.haplotype1(), ind.haplotype2()] {
            for chr in haplotype.chromosomes() {
                let seq = chr.sequence();
                for i in 0..seq.len() {
                    if let Some(nuc) = seq.get(i) {
                        *counts.entry(nuc).or_insert(0) += 1;
                    }
                }
            }
        }
    }

    counts
}

fn composition_population(population: &Population) -> HashMap<Nucleotide, usize> {
    let mut counts = HashMap::new();
    counts.insert(Nucleotide::A, 0);
    counts.insert(Nucleotide::C, 0);
    counts.insert(Nucleotide::G, 0);
    counts.insert(Nucleotide::T, 0);

    for ind in population.individuals() {
        for haplotype in [ind.haplotype1(), ind.haplotype2()] {
            for chr in haplotype.chromosomes() {
                let seq = chr.sequence();
                for i in 0..seq.len() {
                    if let Some(nuc) = seq.get(i) {
                        *counts.entry(nuc).or_insert(0) += 1;
                    }
                }
            }
        }
    }

    counts
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::base::{Alphabet, Sequence};
    use crate::genome::{Chromosome, Haplotype, Individual};

    fn create_test_individual(id: &str, sequence: &[Nucleotide]) -> Individual {
        let alphabet = Alphabet::dna();
        let mut seq = Sequence::with_capacity(sequence.len(), alphabet);
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
    fn test_gc_content_all_gc() {
        let seq = vec![Nucleotide::G, Nucleotide::C, Nucleotide::G, Nucleotide::C];
        let individuals = vec![create_test_individual("ind1", &seq)];
        let pop = Population::new("pop1", individuals);

        // Chromosome level
        let gc = gc_content(&pop, Some(0), Some(0), Some(0));
        assert_eq!(gc, 1.0);
    }

    #[test]
    fn test_gc_content_all_at() {
        let seq = vec![Nucleotide::A, Nucleotide::T, Nucleotide::A, Nucleotide::T];
        let individuals = vec![create_test_individual("ind1", &seq)];
        let pop = Population::new("pop1", individuals);

        // Chromosome level
        let gc = gc_content(&pop, Some(0), Some(0), Some(0));
        assert_eq!(gc, 0.0);
    }

    #[test]
    fn test_gc_content_half() {
        let seq = vec![Nucleotide::G, Nucleotide::C, Nucleotide::A, Nucleotide::T];
        let individuals = vec![create_test_individual("ind1", &seq)];
        let pop = Population::new("pop1", individuals);

        // Chromosome level
        let gc = gc_content(&pop, Some(0), Some(0), Some(0));
        assert_eq!(gc, 0.5);
    }

    #[test]
    fn test_gc_content_population_level() {
        let seq1 = vec![Nucleotide::G, Nucleotide::C]; // GC = 1.0
        let seq2 = vec![Nucleotide::A, Nucleotide::T]; // GC = 0.0

        let individuals = vec![
            create_test_individual("ind1", &seq1),
            create_test_individual("ind2", &seq2),
        ];
        let pop = Population::new("pop1", individuals);

        // Population level - should average across all sequences
        let pop_gc = gc_content(&pop, None, None, None);
        assert_eq!(pop_gc, 0.5);
    }

    #[test]
    fn test_gc_content_individual_level() {
        // Create individual with different GC content in haplotype
        let seq_gc = vec![Nucleotide::G, Nucleotide::C];
        let individuals = vec![create_test_individual("ind1", &seq_gc)];
        let pop = Population::new("pop1", individuals);

        // Individual level (only haplotype1 has data, haplotype2 is empty)
        let ind_gc = gc_content(&pop, Some(0), None, None);
        assert_eq!(ind_gc, 1.0);
    }

    #[test]
    fn test_gc_content_haplotype_level() {
        let seq = vec![Nucleotide::G, Nucleotide::C];
        let individuals = vec![create_test_individual("ind1", &seq)];
        let pop = Population::new("pop1", individuals);

        // Haplotype level
        let hap_gc = gc_content(&pop, Some(0), Some(0), None);
        assert_eq!(hap_gc, 1.0);
    }

    #[test]
    fn test_nucleotide_composition_chromosome() {
        let seq = vec![
            Nucleotide::A,
            Nucleotide::A,
            Nucleotide::T,
            Nucleotide::G,
            Nucleotide::C,
        ];
        let individuals = vec![create_test_individual("ind1", &seq)];
        let pop = Population::new("pop1", individuals);

        // Chromosome level
        let comp = nucleotide_composition(&pop, Some(0), Some(0), Some(0));
        assert_eq!(*comp.get(&Nucleotide::A).unwrap(), 2);
        assert_eq!(*comp.get(&Nucleotide::T).unwrap(), 1);
        assert_eq!(*comp.get(&Nucleotide::G).unwrap(), 1);
        assert_eq!(*comp.get(&Nucleotide::C).unwrap(), 1);
    }

    #[test]
    fn test_nucleotide_composition_population() {
        let seq1 = vec![Nucleotide::A, Nucleotide::A];
        let seq2 = vec![Nucleotide::T, Nucleotide::T];
        let individuals = vec![
            create_test_individual("ind1", &seq1),
            create_test_individual("ind2", &seq2),
        ];
        let pop = Population::new("pop1", individuals);

        // Population level
        let comp = nucleotide_composition(&pop, None, None, None);
        assert_eq!(*comp.get(&Nucleotide::A).unwrap(), 2);
        assert_eq!(*comp.get(&Nucleotide::T).unwrap(), 2);
        assert_eq!(*comp.get(&Nucleotide::G).unwrap(), 0);
        assert_eq!(*comp.get(&Nucleotide::C).unwrap(), 0);
    }

    #[test]
    fn test_nucleotide_composition_individual() {
        let seq = vec![Nucleotide::A, Nucleotide::C];
        let individuals = vec![create_test_individual("ind1", &seq)];
        let pop = Population::new("pop1", individuals);

        // Individual level (haplotype1 has data, haplotype2 is empty)
        let comp = nucleotide_composition(&pop, Some(0), None, None);
        assert_eq!(*comp.get(&Nucleotide::A).unwrap(), 1);
        assert_eq!(*comp.get(&Nucleotide::C).unwrap(), 1);
        assert_eq!(*comp.get(&Nucleotide::G).unwrap(), 0);
        assert_eq!(*comp.get(&Nucleotide::T).unwrap(), 0);
    }

    #[test]
    fn test_invalid_index_combinations() {
        let seq = vec![Nucleotide::A];
        let individuals = vec![create_test_individual("ind1", &seq)];
        let pop = Population::new("pop1", individuals);

        // Invalid: haplotype specified but not individual
        let gc = gc_content(&pop, None, Some(0), None);
        assert_eq!(gc, 0.0);

        // Invalid: chromosome specified but not haplotype
        let gc = gc_content(&pop, Some(0), None, Some(0));
        assert_eq!(gc, 0.0);
    }
}
