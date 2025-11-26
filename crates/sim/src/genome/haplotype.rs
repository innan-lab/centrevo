use crate::base::FitnessValue;
use crate::genome::Chromosome;

// We store ids as `String` for locality and fast, non-atomic cloning.

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
    /// Cached chromosome ids in the same order as `chromosomes`
    ids: Vec<String>,
    /// Cached fitness for this haplotype. `None` indicates not memoized.
    fitness: Option<FitnessValue>,
}

impl Haplotype {
    /// Create a new, empty `Haplotype`.
    pub fn new() -> Self {
        Self {
            chromosomes: Vec::new(),
            ids: Vec::new(),
            fitness: None,
        }
    }

    /// Create a `Haplotype` with reserved capacity for `capacity` chromosomes.
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            chromosomes: Vec::with_capacity(capacity),
            ids: Vec::with_capacity(capacity),
            fitness: None,
        }
    }

    /// Create a `Haplotype` from an existing vector of `Chromosome`s.
    pub fn from_chromosomes(chromosomes: Vec<Chromosome>) -> Self {
        // Build ids in the same order as the chromosomes
        let ids = chromosomes.iter().map(|c| c.id().to_string()).collect();
        Self { chromosomes, ids, fitness: None }
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
        // Invalidate cached fitness when providing mutable access
        self.fitness = None;
        self.chromosomes.get_mut(index)
    }

    /// Append a `Chromosome` to this haplotype.
    pub fn push(&mut self, chromosome: Chromosome) {
        // Capture the id before moving the chromosome into the vector
        let id = chromosome.id().to_string();
        self.ids.push(id);
        self.chromosomes.push(chromosome);
        // Invalidate cached fitness when haplotype changes
        self.fitness = None;
    }

    /// Return the cached fitness value for this haplotype, or `None` if unset.
    #[inline]
    pub fn cached_fitness(&self) -> Option<FitnessValue> {
        self.fitness
    }

    /// Set the cached fitness value for this haplotype.
    #[inline]
    pub fn set_cached_fitness(&mut self, fitness: impl Into<FitnessValue>) {
        self.fitness = Some(fitness.into());
    }

    /// Clear the cached fitness value.
    #[inline]
    pub fn clear_cached_fitness(&mut self) {
        self.fitness = None;
    }

    /// Borrow the slice of chromosomes.
    #[inline]
    pub fn chromosomes(&self) -> &[Chromosome] {
        &self.chromosomes
    }

    /// Borrow the mutable slice of chromosomes.
    #[inline]
    pub fn chromosomes_mut(&mut self) -> &mut [Chromosome] {
        // Invalidate cached fitness whenever mutating access is requested
        self.fitness = None;
        &mut self.chromosomes
    }

    /// Iterate over chromosomes by reference.
    pub fn iter(&self) -> impl Iterator<Item = &Chromosome> {
        self.chromosomes.iter()
    }

    /// Iterate mutably over chromosomes.
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut Chromosome> {
        // Invalidate cached fitness when providing mutable iterators
        self.fitness = None;
        self.chromosomes.iter_mut()
    }

    /// Get a reference to a chromosome by its id, or `None` if not found.
    ///
    /// Uses the internal `ids` cache to find the index and returns the
    /// chromosome at the same position.
    #[inline]
    pub fn get_by_id(&self, id: &str) -> Option<&Chromosome> {
        self.ids.iter().position(|s| s == id).map(|i| &self.chromosomes[i])
    }

    /// Get a mutable reference to a chromosome by its id, or `None` if not
    /// found.
    #[inline]
    pub fn get_by_id_mut(&mut self, id: &str) -> Option<&mut Chromosome> {
        if let Some(pos) = self.ids.iter().position(|s| s == id) {
            self.chromosomes.get_mut(pos)
        } else {
            None
        }
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
    use crate::base::Nucleotide;

    fn test_chromosome(id: &str, length: usize) -> Chromosome {
        // Convert total length to num_hors
        // ru_length=10, rus_per_hor=10, so one HOR = 100 bp
        let hor_length = 100;
        let num_hors = length / hor_length;
        Chromosome::uniform(id, Nucleotide::A, 10, 10, num_hors)
    }

    // ===== Haplotype Tests =====

    #[test]
    fn test_haplotype_new() {
        let hap = Haplotype::new();
        assert_eq!(hap.len(), 0);
        assert!(hap.is_empty());
        assert_eq!(hap.cached_fitness(), None);
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
        assert_eq!(hap.cached_fitness(), None);
    }

    #[test]
    fn test_haplotype_push() {
        let mut hap = Haplotype::new();

        // Set a fitness and ensure it's invalidated by push
        hap.set_cached_fitness(FitnessValue::new(0.42));
        assert_eq!(hap.cached_fitness(), Some(FitnessValue::new(0.42)));

        hap.push(test_chromosome("chr1", 100));
        assert_eq!(hap.len(), 1);

        hap.push(test_chromosome("chr2", 200));
        assert_eq!(hap.len(), 2);

        hap.push(test_chromosome("chr3", 300));
        assert_eq!(hap.len(), 3);
        assert_eq!(hap.cached_fitness(), None);
    }

    #[test]
    fn test_haplotype_cached_fitness() {
        let hap = Haplotype::new();

        assert_eq!(hap.cached_fitness(), None);
    }

    #[test]
    fn test_haplotype_set_cached_fitness() {
        let mut hap = Haplotype::new();

        hap.set_cached_fitness(FitnessValue::new(0.85));
        assert_eq!(hap.cached_fitness(), Some(FitnessValue::new(0.85)));

        hap.set_cached_fitness(FitnessValue::new(0.95));
        assert_eq!(hap.cached_fitness(), Some(FitnessValue::new(0.95)));
    }

    #[test]
    fn test_haplotype_set_clear_cached_fitness() {
        let mut hap = Haplotype::new();
        
        hap.set_cached_fitness(FitnessValue::new(0.75));
        assert_eq!(hap.cached_fitness(), Some(FitnessValue::new(0.75)));

        hap.clear_cached_fitness();
        assert_eq!(hap.cached_fitness(), None);
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
        // Mutable access should invalidate cached fitness
        assert_eq!(hap.cached_fitness(), None);
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
        // Mutable access should invalidate cached fitness
        assert_eq!(hap.cached_fitness(), None);
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
        // Mutable iterator should clear cached fitness
        assert_eq!(hap.cached_fitness(), None);
    }

    #[test]
    fn test_haplotype_into_iter() {
        let mut hap = Haplotype::new();
        hap.push(test_chromosome("chr1", 100));
        hap.push(test_chromosome("chr2", 200));

        let ids: Vec<String> = hap.into_iter().map(|chr| chr.id().to_string()).collect();

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
        // Cached fitness should be cloned as well
        hap1.set_cached_fitness(FitnessValue::new(0.33));
        let hap2 = hap1.clone();
        assert_eq!(hap1.cached_fitness(), hap2.cached_fitness());
    }

    #[test]
    fn test_haplotype_many_chromosomes() {
        let mut hap = Haplotype::new();

        // Add 23 chromosomes (human genome)
        for i in 1..=23 {
            hap.push(test_chromosome(&format!("chr{i}"), i * 100));
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
            hap.push(test_chromosome(&format!("chr{i}"), i * 100));
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
        hap2.get_mut(0)
            .unwrap()
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
