//! Example demonstrating flexible composition analysis at different levels

use centrevo::base::{Alphabet, Nucleotide, Sequence};
use centrevo::genome::{Chromosome, Haplotype, Individual};
use centrevo::simulation::Population;
use centrevo::analysis::composition::{gc_content, nucleotide_composition};

fn main() {
    // Create a small test population
    let population = create_example_population();

    println!("=== Flexible Composition Analysis Demo ===\n");

    // Population-level analysis
    println!("1. Population Level (all sequences):");
    let pop_gc = gc_content(&population, None, None, None);
    println!("   GC content: {:.2}%", pop_gc * 100.0);
    let pop_comp = nucleotide_composition(&population, None, None, None);
    println!("   Composition: {:?}\n", pop_comp);

    // Individual-level analysis
    println!("2. Individual Level (both haplotypes):");
    for i in 0..population.size() {
        let ind_gc = gc_content(&population, Some(i), None, None);
        println!("   Individual {}: GC = {:.2}%", i, ind_gc * 100.0);
    }
    println!();

    // Haplotype-level analysis
    println!("3. Haplotype Level (all chromosomes in haplotype):");
    let hap0_gc = gc_content(&population, Some(0), Some(0), None);
    let hap1_gc = gc_content(&population, Some(0), Some(1), None);
    println!("   Individual 0, Haplotype 0: GC = {:.2}%", hap0_gc * 100.0);
    println!("   Individual 0, Haplotype 1: GC = {:.2}%", hap1_gc * 100.0);
    println!();

    // Chromosome-level analysis
    println!("4. Chromosome Level (single chromosome):");
    let chr_gc = gc_content(&population, Some(0), Some(0), Some(0));
    println!("   Individual 0, Haplotype 0, Chromosome 0: GC = {:.2}%", chr_gc * 100.0);
    let chr_comp = nucleotide_composition(&population, Some(0), Some(0), Some(0));
    println!("   Nucleotide counts: {:?}", chr_comp);
}

fn create_example_population() -> Population {
    let alphabet = Alphabet::dna();
    
    // Create sequences with different GC content
    let high_gc = create_sequence(&alphabet, &[
        Nucleotide::G, Nucleotide::C, Nucleotide::G, Nucleotide::C,
        Nucleotide::G, Nucleotide::C, Nucleotide::G, Nucleotide::C,
    ]);
    
    let low_gc = create_sequence(&alphabet, &[
        Nucleotide::A, Nucleotide::T, Nucleotide::A, Nucleotide::T,
        Nucleotide::A, Nucleotide::T, Nucleotide::A, Nucleotide::T,
    ]);
    
    let mid_gc = create_sequence(&alphabet, &[
        Nucleotide::A, Nucleotide::T, Nucleotide::G, Nucleotide::C,
        Nucleotide::A, Nucleotide::T, Nucleotide::G, Nucleotide::C,
    ]);

    // Create individuals
    let ind1 = create_individual("ind1", high_gc.clone(), low_gc.clone());
    let ind2 = create_individual("ind2", mid_gc.clone(), mid_gc.clone());
    
    Population::new("example_pop", vec![ind1, ind2])
}

fn create_sequence(alphabet: &Alphabet, nucleotides: &[Nucleotide]) -> Sequence {
    let mut seq = Sequence::with_capacity(nucleotides.len(), alphabet.clone());
    for &nuc in nucleotides {
        seq.push(nuc);
    }
    seq
}

fn create_individual(id: &str, seq1: Sequence, seq2: Sequence) -> Individual {
    let chr1 = Chromosome::new(format!("{}_chr1", id), seq1, 4, 2);
    let chr2 = Chromosome::new(format!("{}_chr2", id), seq2, 4, 2);
    
    let mut hap1 = Haplotype::new();
    hap1.push(chr1);
    
    let mut hap2 = Haplotype::new();
    hap2.push(chr2);
    
    Individual::new(id, hap1, hap2)
}
