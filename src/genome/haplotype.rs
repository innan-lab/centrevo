use crate::genome::Chromosome;

/// A haplotype: an ordered collection of chromosomes representing one set of
/// genetic material for an individual.
///
/// `Haplotype` owns its `Chromosome`s and provides simple accessors to query
/// chromosomes by index, iterate over them, and obtain aggregate properties
/// (for example total length). For read-only, parallel scenarios prefer
/// converting to `SharedHaplotype` (via `to_shared`) which contains
/// `SharedChromosome`s and is cheap to clone.
#[derive(Debug, Clone)]
pub struct Haplotype {
    /// Chromosomes in this haplotype
    chromosomes: Vec<Chromosome>,
}

impl Haplotype {
    /// Create a new, empty `Haplotype`.
    pub fn new() -> Self {
        Self {
            chromosomes: Vec::new(),
        }
    }

    /// Create a `Haplotype` with reserved capacity for `capacity` chromosomes.
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            chromosomes: Vec::with_capacity(capacity),
        }
    }

    /// Create a `Haplotype` from an existing vector of `Chromosome`s.
    pub fn from_chromosomes(chromosomes: Vec<Chromosome>) -> Self {
        Self { chromosomes }
    }

    /// Return the number of chromosomes in this haplotype.
    #[inline]
    pub fn len(&self) -> usize {
        self.chromosomes.len()
    }

    /// Return `true` if this haplotype contains no chromosomes.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.chromosomes.is_empty()
    }

    /// Get a reference to the chromosome at `index`, or `None` if out of
    /// bounds.
    #[inline]
    pub fn get(&self, index: usize) -> Option<&Chromosome> {
        self.chromosomes.get(index)
    }

    /// Get a mutable reference to the chromosome at `index`, or `None` if out
    /// of bounds.
    #[inline]
    pub fn get_mut(&mut self, index: usize) -> Option<&mut Chromosome> {
        self.chromosomes.get_mut(index)
    }

    /// Append a `Chromosome` to this haplotype.
    pub fn push(&mut self, chromosome: Chromosome) {
        self.chromosomes.push(chromosome);
    }

    /// Borrow the slice of chromosomes.
    #[inline]
    pub fn chromosomes(&self) -> &[Chromosome] {
        &self.chromosomes
    }

    /// Borrow the mutable slice of chromosomes.
    #[inline]
    pub fn chromosomes_mut(&mut self) -> &mut [Chromosome] {
        &mut self.chromosomes
    }

    /// Iterate over chromosomes by reference.
    pub fn iter(&self) -> impl Iterator<Item = &Chromosome> {
        self.chromosomes.iter()
    }

    /// Iterate mutably over chromosomes.
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut Chromosome> {
        self.chromosomes.iter_mut()
    }

    /// Return the total number of bases across all chromosomes in this
    /// haplotype.
    pub fn total_length(&self) -> usize {
        self.chromosomes.iter().map(|chr| chr.len()).sum()
    }
}

impl Default for Haplotype {
    fn default() -> Self {
        Self::new()
    }
}

impl IntoIterator for Haplotype {
    type Item = Chromosome;
    type IntoIter = std::vec::IntoIter<Chromosome>;

    fn into_iter(self) -> Self::IntoIter {
        self.chromosomes.into_iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::base::{Alphabet, Nucleotide};

    fn test_alphabet() -> Alphabet {
        Alphabet::dna()
    }

    fn test_chromosome(id: &str, length: usize) -> Chromosome {
        Chromosome::uniform(id, Nucleotide::A, length, 10, 5, test_alphabet())
    }

    // ===== Haplotype Tests =====

    #[test]
    fn test_haplotype_new() {
        let hap = Haplotype::new();
        assert_eq!(hap.len(), 0);
        assert!(hap.is_empty());
    }

    #[test]
    fn test_haplotype_default() {
        let hap = Haplotype::default();
        assert_eq!(hap.len(), 0);
        assert!(hap.is_empty());
    }

    #[test]
    fn test_haplotype_with_capacity() {
        let hap = Haplotype::with_capacity(10);
        assert_eq!(hap.len(), 0);
        assert!(hap.is_empty());
    }

    #[test]
    fn test_haplotype_from_chromosomes() {
        let chr1 = test_chromosome("chr1", 100);
        let chr2 = test_chromosome("chr2", 200);
        
        let hap = Haplotype::from_chromosomes(vec![chr1, chr2]);
        assert_eq!(hap.len(), 2);
        assert!(!hap.is_empty());
    }

    #[test]
    fn test_haplotype_push() {
        let mut hap = Haplotype::new();
        
        hap.push(test_chromosome("chr1", 100));
        assert_eq!(hap.len(), 1);
        
        hap.push(test_chromosome("chr2", 200));
        assert_eq!(hap.len(), 2);
        
        hap.push(test_chromosome("chr3", 300));
        assert_eq!(hap.len(), 3);
    }

    #[test]
    fn test_haplotype_get() {
        let mut hap = Haplotype::new();
        hap.push(test_chromosome("chr1", 100));
        hap.push(test_chromosome("chr2", 200));
        
        let chr1 = hap.get(0);
        assert!(chr1.is_some());
        assert_eq!(chr1.unwrap().id(), "chr1");
        assert_eq!(chr1.unwrap().len(), 100);
        
        let chr2 = hap.get(1);
        assert!(chr2.is_some());
        assert_eq!(chr2.unwrap().id(), "chr2");
        assert_eq!(chr2.unwrap().len(), 200);
        
        let chr3 = hap.get(2);
        assert!(chr3.is_none());
    }

    #[test]
    fn test_haplotype_get_mut() {
        let mut hap = Haplotype::new();
        hap.push(test_chromosome("chr1", 100));
        
        if let Some(chr) = hap.get_mut(0) {
            chr.sequence_mut().set(0, Nucleotide::T).unwrap();
        }
        
        let chr = hap.get(0).unwrap();
        assert_eq!(chr.sequence().get(0), Some(Nucleotide::T));
    }

    #[test]
    fn test_haplotype_chromosomes() {
        let mut hap = Haplotype::new();
        hap.push(test_chromosome("chr1", 100));
        hap.push(test_chromosome("chr2", 200));
        
        let chrs = hap.chromosomes();
        assert_eq!(chrs.len(), 2);
        assert_eq!(chrs[0].id(), "chr1");
        assert_eq!(chrs[1].id(), "chr2");
    }

    #[test]
    fn test_haplotype_chromosomes_mut() {
        let mut hap = Haplotype::new();
        hap.push(test_chromosome("chr1", 100));
        
        let chrs = hap.chromosomes_mut();
        chrs[0].sequence_mut().set(0, Nucleotide::G).unwrap();
        
        assert_eq!(hap.get(0).unwrap().sequence().get(0), Some(Nucleotide::G));
    }

    #[test]
    fn test_haplotype_iter() {
        let mut hap = Haplotype::new();
        hap.push(test_chromosome("chr1", 100));
        hap.push(test_chromosome("chr2", 200));
        hap.push(test_chromosome("chr3", 300));
        
        let ids: Vec<&str> = hap.iter().map(|chr| chr.id()).collect();
        assert_eq!(ids, vec!["chr1", "chr2", "chr3"]);
    }

    #[test]
    fn test_haplotype_iter_mut() {
        let mut hap = Haplotype::new();
        hap.push(test_chromosome("chr1", 100));
        hap.push(test_chromosome("chr2", 200));
        
        // Mutate all chromosomes
        for chr in hap.iter_mut() {
            chr.sequence_mut().set(0, Nucleotide::C).unwrap();
        }
        
        // Verify mutations
        for chr in hap.iter() {
            assert_eq!(chr.sequence().get(0), Some(Nucleotide::C));
        }
    }

    #[test]
    fn test_haplotype_into_iter() {
        let mut hap = Haplotype::new();
        hap.push(test_chromosome("chr1", 100));
        hap.push(test_chromosome("chr2", 200));
        
        let ids: Vec<String> = hap.into_iter()
            .map(|chr| chr.id().to_string())
            .collect();
        
        assert_eq!(ids, vec!["chr1", "chr2"]);
    }

    #[test]
    fn test_haplotype_total_length() {
        let mut hap = Haplotype::new();
        hap.push(test_chromosome("chr1", 100));
        hap.push(test_chromosome("chr2", 200));
        hap.push(test_chromosome("chr3", 300));
        
        assert_eq!(hap.total_length(), 600);
    }

    #[test]
    fn test_haplotype_total_length_empty() {
        let hap = Haplotype::new();
        assert_eq!(hap.total_length(), 0);
    }

    #[test]
    fn test_haplotype_total_length_single() {
        let mut hap = Haplotype::new();
        hap.push(test_chromosome("chr1", 100));
        assert_eq!(hap.total_length(), 100);
    }

    #[test]
    fn test_haplotype_clone() {
        let mut hap1 = Haplotype::new();
        hap1.push(test_chromosome("chr1", 100));
        hap1.push(test_chromosome("chr2", 200));
        
        let hap2 = hap1.clone();
        
        assert_eq!(hap1.len(), hap2.len());
        assert_eq!(hap1.total_length(), hap2.total_length());
        
        // Verify chromosome IDs match
        for (chr1, chr2) in hap1.iter().zip(hap2.iter()) {
            assert_eq!(chr1.id(), chr2.id());
        }
    }

    #[test]
    fn test_haplotype_many_chromosomes() {
        let mut hap = Haplotype::new();
        
        // Add 23 chromosomes (human genome)
        for i in 1..=23 {
            hap.push(test_chromosome(&format!("chr{}", i), i * 100));
        }
        
        assert_eq!(hap.len(), 23);
        assert!(hap.total_length() > 0);
    }

    #[test]
    fn test_haplotype_empty_after_creation() {
        let hap = Haplotype::from_chromosomes(vec![]);
        assert!(hap.is_empty());
        assert_eq!(hap.len(), 0);
    }

    #[test]
    fn test_haplotype_access_pattern() {
        let mut hap = Haplotype::new();
        
        // Add some chromosomes
        for i in 1..=5 {
            hap.push(test_chromosome(&format!("chr{}", i), i * 100));
        }
        
        // Random access
        assert_eq!(hap.get(2).unwrap().id(), "chr3");
        assert_eq!(hap.get(2).unwrap().len(), 300);
        
        // Sequential access
        let mut lengths = Vec::new();
        for chr in hap.iter() {
            lengths.push(chr.len());
        }
        assert_eq!(lengths, vec![100, 200, 300, 400, 500]);
    }

    #[test]
    fn test_haplotype_mutation_isolation() {
        let mut hap1 = Haplotype::new();
        hap1.push(test_chromosome("chr1", 100));
        
        let mut hap2 = hap1.clone();
        
        // Mutate hap2
        hap2.get_mut(0).unwrap()
            .sequence_mut()
            .set(0, Nucleotide::T)
            .unwrap();
        
        // hap1 should be unchanged
        assert_eq!(hap1.get(0).unwrap().sequence().get(0), Some(Nucleotide::A));
        assert_eq!(hap2.get(0).unwrap().sequence().get(0), Some(Nucleotide::T));
    }

    #[test]
    fn test_haplotype_large() {
        let mut hap = Haplotype::new();
        
        // Create haplotype with large chromosomes
        hap.push(test_chromosome("chr1", 1_000_000));
        hap.push(test_chromosome("chr2", 500_000));
        
        assert_eq!(hap.len(), 2);
        assert_eq!(hap.total_length(), 1_500_000);
    }

    #[test]
    fn test_haplotype_realistic_human() {
        // Simplified human genome (just a few chromosomes for testing)
        let mut hap = Haplotype::new();
        
        // Roughly proportional sizes (in kb, scaled down)
        hap.push(test_chromosome("chr1", 248_000));
        hap.push(test_chromosome("chr2", 242_000));
        hap.push(test_chromosome("chrX", 156_000));
        
        assert_eq!(hap.len(), 3);
        assert!(hap.total_length() > 600_000);
    }
}
