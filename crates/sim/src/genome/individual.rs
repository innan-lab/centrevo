use crate::base::FitnessValue;
use crate::genome::Haplotype;
use std::sync::Arc;

/// An individual organism with a diploid genome.
///
/// `Individual` contains two `Haplotype`s (representing the two chromosome
/// sets) and a cached fitness value. The `id` is stored in an `Arc<str>` so
/// cloning individuals is cheap for the identifier field. Use the provided
/// accessors to read or mutate haplotypes and fitness as needed.
#[derive(Debug, Clone)]
pub struct Individual {
    /// Unique identifier
    id: Arc<str>,
    /// First haplotype
    haplotype1: Haplotype,
    /// Second haplotype
    haplotype2: Haplotype,
    /// Cached fitness value. `None` indicates that the fitness has not
    /// been computed/memoized yet.
    fitness: Option<FitnessValue>,
}

impl Individual {
    /// Create a new `Individual` from two haplotypes.
    ///
    /// The `id` can be any type convertible into `Arc<str>` (for example `&str`
    /// or `String`). The initial cached fitness is `None` and may be updated
    /// later with `set_fitness`.
    pub fn new(id: impl Into<Arc<str>>, haplotype1: Haplotype, haplotype2: Haplotype) -> Self {
        Self {
            id: id.into(),
            haplotype1,
            haplotype2,
            fitness: None,
        }
    }

    /// Return the individual's identifier as a `&str`.
    #[inline]
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Borrow the first haplotype (read-only).
    #[inline]
    pub fn haplotype1(&self) -> &Haplotype {
        &self.haplotype1
    }

    /// Borrow the first haplotype mutably to perform in-place modifications.
    #[inline]
    pub fn haplotype1_mut(&mut self) -> &mut Haplotype {
        &mut self.haplotype1
    }

    /// Borrow the second haplotype (read-only).
    #[inline]
    pub fn haplotype2(&self) -> &Haplotype {
        &self.haplotype2
    }

    /// Borrow the second haplotype mutably to perform in-place modifications.
    #[inline]
    pub fn haplotype2_mut(&mut self) -> &mut Haplotype {
        &mut self.haplotype2
    }

    /// Return the cached fitness value for this individual.
    ///
    /// Returns `None` if the fitness has not yet been computed.
    #[inline]
    pub fn cached_fitness(&self) -> Option<FitnessValue> {
        self.fitness
    }

    /// Set the cached fitness value for this individual.
    #[inline]
    pub fn set_cached_fitness(&mut self, fitness: impl Into<FitnessValue>) {
        self.fitness = Some(fitness.into());
    }

    /// Clear the cached fitness value, indicating it needs to be recomputed.
    #[inline]
    pub fn clear_cached_fitness(&mut self) {
        self.fitness = None;
    }

    /// Borrow both haplotypes as a pair of references.
    pub fn haplotypes(&self) -> (&Haplotype, &Haplotype) {
        (&self.haplotype1, &self.haplotype2)
    }

    /// Borrow both haplotypes mutably as a pair of mutable references.
    pub fn haplotypes_mut(&mut self) -> (&mut Haplotype, &mut Haplotype) {
        (&mut self.haplotype1, &mut self.haplotype2)
    }

    /// Offset the page IDs of all chromosomes.
    /// Used when merging arenas.
    pub fn offset_page_ids(&mut self, offset: u32) {
        self.haplotype1.offset_page_ids(offset);
        self.haplotype2.offset_page_ids(offset);
    }

    /// Localize all chromosomes in this individual.
    pub fn localize(
        &mut self,
        source_arena: &crate::base::GenomeArena,
        dest_arena: &mut crate::base::GenomeArena,
    ) {
        self.haplotype1.localize(source_arena, dest_arena);
        self.haplotype2.localize(source_arena, dest_arena);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::base::{GenomeArena, Nucleotide};
    use crate::genome::Chromosome;

    fn test_chromosome(id: &str, length: usize, arena: &mut GenomeArena) -> Chromosome {
        let hor_length = 100;
        let num_hors = length / hor_length;
        let num_hors = if num_hors == 0 { 1 } else { num_hors };
        Chromosome::uniform(id, Nucleotide::A, 10, 10, num_hors, arena)
    }

    fn test_haplotype(num_chrs: usize, arena: &mut GenomeArena) -> Haplotype {
        let mut hap = Haplotype::new();
        for i in 1..=num_chrs {
            hap.push(test_chromosome(&format!("chr{i}"), i * 100, arena));
        }
        hap
    }

    // ===== Individual Tests =====

    #[test]
    fn test_individual_new() {
        let mut arena = GenomeArena::new();
        let hap1 = test_haplotype(2, &mut arena);
        let hap2 = test_haplotype(2, &mut arena);

        let ind = Individual::new("ind1", hap1, hap2);

        assert_eq!(ind.id(), "ind1");
        assert_eq!(ind.haplotype1().len(), 2);
        assert_eq!(ind.haplotype2().len(), 2);
        assert_eq!(ind.cached_fitness(), None);
    }

    #[test]
    fn test_individual_id() {
        let mut arena = GenomeArena::new();
        let ind = Individual::new(
            "test_id",
            test_haplotype(1, &mut arena),
            test_haplotype(1, &mut arena),
        );
        assert_eq!(ind.id(), "test_id");
    }

    #[test]
    fn test_individual_id_shared() {
        let mut arena = GenomeArena::new();
        let ind1 = Individual::new(
            "ind1",
            test_haplotype(1, &mut arena),
            test_haplotype(1, &mut arena),
        );
        let ind2 = ind1.clone();

        // Both should share the same ID Arc
        assert_eq!(ind1.id(), ind2.id());
    }

    #[test]
    fn test_individual_haplotype1() {
        let mut arena = GenomeArena::new();
        let hap1 = test_haplotype(3, &mut arena);
        let hap2 = test_haplotype(2, &mut arena);

        let ind = Individual::new("ind1", hap1, hap2);

        assert_eq!(ind.haplotype1().len(), 3);
        assert_eq!(ind.haplotype1().total_length(), 600); // 100 + 200 + 300
    }

    #[test]
    fn test_individual_haplotype2() {
        let mut arena = GenomeArena::new();
        let hap1 = test_haplotype(2, &mut arena);
        let hap2 = test_haplotype(3, &mut arena);

        let ind = Individual::new("ind1", hap1, hap2);

        assert_eq!(ind.haplotype2().len(), 3);
        assert_eq!(ind.haplotype2().total_length(), 600);
    }

    #[test]
    fn test_individual_haplotype1_mut() {
        let mut arena = GenomeArena::new();
        let hap1 = test_haplotype(2, &mut arena);
        let hap2 = test_haplotype(2, &mut arena);

        let mut ind = Individual::new("ind1", hap1, hap2);

        // Mutate first haplotype
        ind.haplotype1_mut()
            .push(test_chromosome("chr3", 300, &mut arena));

        assert_eq!(ind.haplotype1().len(), 3);
        assert_eq!(ind.haplotype2().len(), 2); // Second haplotype unchanged
    }

    #[test]
    fn test_individual_haplotype2_mut() {
        let mut arena = GenomeArena::new();
        let hap1 = test_haplotype(2, &mut arena);
        let hap2 = test_haplotype(2, &mut arena);

        let mut ind = Individual::new("ind1", hap1, hap2);

        // Mutate second haplotype
        ind.haplotype2_mut()
            .push(test_chromosome("chr3", 300, &mut arena));

        assert_eq!(ind.haplotype1().len(), 2); // First haplotype unchanged
        assert_eq!(ind.haplotype2().len(), 3);
    }

    #[test]
    fn test_individual_cached_fitness() {
        let mut arena = GenomeArena::new();
        let ind = Individual::new(
            "ind1",
            test_haplotype(2, &mut arena),
            test_haplotype(2, &mut arena),
        );
        assert_eq!(ind.cached_fitness(), None);
    }

    #[test]
    fn test_individual_set_cached_fitness() {
        let mut arena = GenomeArena::new();
        let mut ind = Individual::new(
            "ind1",
            test_haplotype(2, &mut arena),
            test_haplotype(2, &mut arena),
        );

        ind.set_cached_fitness(FitnessValue::new(0.75));
        assert_eq!(ind.cached_fitness(), Some(FitnessValue::new(0.75)));

        ind.set_cached_fitness(FitnessValue::new(1.0));
        assert_eq!(ind.cached_fitness(), Some(FitnessValue::new(1.0)));
    }

    #[test]
    fn test_individual_clear_cached_fitness() {
        let mut arena = GenomeArena::new();
        let mut ind = Individual::new(
            "ind1",
            test_haplotype(2, &mut arena),
            test_haplotype(2, &mut arena),
        );
        ind.set_cached_fitness(FitnessValue::new(0.9));
        assert_eq!(ind.cached_fitness(), Some(FitnessValue::new(0.9)));
        ind.clear_cached_fitness();
        assert_eq!(ind.cached_fitness(), None);
    }

    #[test]
    fn test_individual_haplotypes() {
        let mut arena = GenomeArena::new();
        let hap1 = test_haplotype(2, &mut arena);
        let hap2 = test_haplotype(3, &mut arena);

        let ind = Individual::new("ind1", hap1, hap2);

        let (h1, h2) = ind.haplotypes();
        assert_eq!(h1.len(), 2);
        assert_eq!(h2.len(), 3);
    }

    #[test]
    fn test_individual_haplotypes_mut() {
        let mut arena = GenomeArena::new();
        let hap1 = test_haplotype(2, &mut arena);
        let hap2 = test_haplotype(2, &mut arena);

        let mut ind = Individual::new("ind1", hap1, hap2);

        let (h1, h2) = ind.haplotypes_mut();
        h1.push(test_chromosome("chr3", 300, &mut arena));
        h2.push(test_chromosome("chr3", 300, &mut arena));

        assert_eq!(ind.haplotype1().len(), 3);
        assert_eq!(ind.haplotype2().len(), 3);
    }

    #[test]
    fn test_individual_clone() {
        let mut arena = GenomeArena::new();
        let mut ind1 = Individual::new(
            "ind1",
            test_haplotype(2, &mut arena),
            test_haplotype(2, &mut arena),
        );
        ind1.set_cached_fitness(FitnessValue::new(0.8));

        let ind2 = ind1.clone();

        assert_eq!(ind1.id(), ind2.id());
        assert_eq!(ind1.cached_fitness(), ind2.cached_fitness());
        assert_eq!(ind1.haplotype1().len(), ind2.haplotype1().len());
        assert_eq!(ind1.haplotype2().len(), ind2.haplotype2().len());
    }

    #[test]
    fn test_individual_diploid() {
        // Test that individual properly maintains diploid state
        let mut arena = GenomeArena::new();
        let hap1 = test_haplotype(3, &mut arena);
        let hap2 = test_haplotype(3, &mut arena);

        let ind = Individual::new("ind1", hap1, hap2);

        // Both haplotypes should exist
        assert_eq!(ind.haplotype1().len(), 3);
        assert_eq!(ind.haplotype2().len(), 3);

        // Should be independent
        let chr1_0 = ind.haplotype1().get(0).unwrap();
        let chr2_0 = ind.haplotype2().get(0).unwrap();

        // Both exist but can be different
        assert_eq!(chr1_0.id(), "chr1");
        assert_eq!(chr2_0.id(), "chr1");
    }

    #[test]
    fn test_individual_heterozygous() {
        let mut arena = GenomeArena::new();
        // Create individual with different haplotypes (heterozygous)
        let mut hap1 = Haplotype::new();
        hap1.push(Chromosome::uniform(
            "chr1",
            Nucleotide::A,
            100,
            10,
            5,
            &mut arena,
        ));

        let mut hap2 = Haplotype::new();
        hap2.push(Chromosome::uniform(
            "chr1",
            Nucleotide::T,
            100,
            10,
            5,
            &mut arena,
        ));

        let ind = Individual::new("ind1", hap1, hap2);

        // Check that haplotypes are different
        let seq1 = ind.haplotype1().get(0).unwrap().sequence(&arena)[0];
        let seq2 = ind.haplotype2().get(0).unwrap().sequence(&arena)[0];

        assert_eq!(seq1, Nucleotide::A);
        assert_eq!(seq2, Nucleotide::T);
    }

    #[test]
    fn test_individual_homozygous() {
        let mut arena = GenomeArena::new();
        // Create individual with identical haplotypes (homozygous)
        let hap1 = test_haplotype(2, &mut arena);
        let hap2 = hap1.clone();

        let ind = Individual::new("ind1", hap1, hap2);

        // Both haplotypes should be the same
        assert_eq!(
            ind.haplotype1().total_length(),
            ind.haplotype2().total_length()
        );
    }

    #[test]
    fn test_individual_mutation_isolation() {
        let mut arena = GenomeArena::new();
        let ind1 = Individual::new(
            "ind1",
            test_haplotype(2, &mut arena),
            test_haplotype(2, &mut arena),
        );
        let mut ind2 = ind1.clone();

        // Mutate ind2's fitness
        ind2.set_cached_fitness(FitnessValue::new(0.5));

        // ind1 should be unchanged
        assert_eq!(ind1.cached_fitness(), None);
        assert_eq!(ind2.cached_fitness(), Some(FitnessValue::new(0.5)));

        // Mutate ind2's haplotype (by adding a chromosome, since in-place sequence mutation isn't direct anymore)
        ind2.haplotype1_mut()
            .push(test_chromosome("chr3", 300, &mut arena));

        // ind1's haplotype should be unchanged
        assert_eq!(ind1.haplotype1().len(), 2);
        assert_eq!(ind2.haplotype1().len(), 3);
    }

    #[test]
    fn test_individual_complex_genome() {
        let mut arena = GenomeArena::new();
        // Create individual with complex genome structure
        let mut hap1 = Haplotype::new();
        let mut hap2 = Haplotype::new();

        // Add multiple chromosomes with different structures
        for i in 1..=5 {
            let total_length = i * 1000;
            let ru_length = 10;
            let rus_per_hor = 10;
            let num_hors = total_length / (ru_length * rus_per_hor);

            hap1.push(Chromosome::uniform(
                format!("chr{i}"),
                Nucleotide::A,
                ru_length,
                rus_per_hor,
                num_hors,
                &mut arena,
            ));
            hap2.push(Chromosome::uniform(
                format!("chr{i}"),
                Nucleotide::C,
                ru_length,
                rus_per_hor,
                num_hors,
                &mut arena,
            ));
        }

        let ind = Individual::new("ind1", hap1, hap2);

        assert_eq!(ind.haplotype1().len(), 5);
        assert_eq!(ind.haplotype2().len(), 5);
        assert_eq!(ind.haplotype1().total_length(), 15_000); // 1k + 2k + 3k + 4k + 5k
        assert_eq!(ind.haplotype2().total_length(), 15_000);
    }

    #[test]
    fn test_individual_large_genome() {
        let mut arena = GenomeArena::new();
        let mut hap1 = Haplotype::new();
        let mut hap2 = Haplotype::new();

        for i in 1..=23 {
            hap1.push(test_chromosome(&format!("chr{i}"), i * 10_000, &mut arena));
            hap2.push(test_chromosome(&format!("chr{i}"), i * 10_000, &mut arena));
        }

        let ind = Individual::new("ind1", hap1, hap2);

        assert_eq!(ind.haplotype1().len(), 23);
        assert_eq!(ind.haplotype2().len(), 23);
        assert!(ind.haplotype1().total_length() > 1_000_000);
    }
}
