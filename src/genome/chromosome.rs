use crate::base::{Nucleotide, Sequence, SharedSequence};
use std::sync::Arc;

/// A chromosome: a named sequence organized into repeat units (RUs) and
/// higher-order repeats (HORs).
///
/// A `Chromosome` owns a mutable `Sequence` and provides helpers that expose
/// the repeat structure (RU length, RUs per HOR) and convenience methods such
/// as GC content and formatted string output. For parallel/read-only
/// operations prefer `SharedChromosome` which shares sequence storage.
#[derive(Debug, Clone)]
pub struct Chromosome {
    /// Type identifier. Unique in a haplotype, but not unique in a population
    /// (shared across clones)
    id: Arc<str>,
    /// The sequence data (mutable during simulation)
    sequence: Sequence,
    /// Repeat unit length in bases
    ru_length: usize,
    /// Number of repeat units per higher-order repeat
    rus_per_hor: usize,
}

impl Chromosome {
    /// Create a new chromosome
    pub fn new(
        id: impl Into<Arc<str>>,
        sequence: Sequence,
        ru_length: usize,
        rus_per_hor: usize,
    ) -> Self {
        Self {
            id: id.into(),
            sequence,
            ru_length,
            rus_per_hor,
        }
    }

    /// Create a uniform chromosome where every base is the same `base`.
    ///
    /// This is a convenience constructor often used in tests and for creating
    /// deterministic initial populations.
    pub fn uniform(
        id: impl Into<Arc<str>>,
        base: Nucleotide,
        length: usize,
        ru_length: usize,
        rus_per_hor: usize,
    ) -> Self {
        let mut sequence = Sequence::with_capacity(length);
        for _ in 0..length {
            sequence.push(base);
        }

        Self::new(id, sequence, ru_length, rus_per_hor)
    }

    /// Return the chromosome identifier.
    ///
    /// The ID is stored in an `Arc<str>` so cloning the `Chromosome` cheaply
    /// shares the identifier.
    #[inline]
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Borrow the underlying mutable `Sequence` for read-only operations.
    #[inline]
    pub fn sequence(&self) -> &Sequence {
        &self.sequence
    }

    /// Borrow the underlying `Sequence` mutably to apply in-place modifications
    /// such as mutation or recombination.
    #[inline]
    pub fn sequence_mut(&mut self) -> &mut Sequence {
        &mut self.sequence
    }

    /// Return the length of the chromosome in bases.
    #[inline]
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// Return true if the chromosome contains no bases.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// Get the repeat unit (RU) length in bases.
    #[inline]
    pub fn ru_length(&self) -> usize {
        self.ru_length
    }

    /// Get the number of repeat units per higher-order repeat (HOR).
    #[inline]
    pub fn rus_per_hor(&self) -> usize {
        self.rus_per_hor
    }

    /// Return the HOR length in bases (RU length × RUs per HOR).
    #[inline]
    pub fn hor_length(&self) -> usize {
        self.ru_length * self.rus_per_hor
    }

    /// Return the number of complete HORs contained in this chromosome
    /// (integer division).
    #[inline]
    pub fn num_hors(&self) -> usize {
        self.len() / self.hor_length()
    }

    /// Calculate GC content of the chromosome as a proportion in [0.0, 1.0].
    pub fn gc_content(&self) -> f64 {
        let mut gc_count = 0;
        let total = self.sequence.len();

        if total == 0 {
            return 0.0;
        }

        for &nuc in self.sequence.as_slice() {
            if matches!(nuc, Nucleotide::G | Nucleotide::C) {
                gc_count += 1;
            }
        }

        gc_count as f64 / total as f64
    }

    /// Convert the chromosome sequence to a formatted string inserting
    /// delimiters between repeat units and higher-order repeats.
    ///
    /// `ru_delim` is used between repeat units, `hor_delim` between HORs.
    pub fn to_formatted_string(&self, ru_delim: char, hor_delim: char) -> String {
        let chars: Vec<char> = self
            .sequence
            .as_slice()
            .iter()
            .map(|&nuc| nuc.to_char())
            .collect();

        let mut result = String::with_capacity(chars.len() * 2);
        let hor_len = self.hor_length();

        for (i, c) in chars.iter().enumerate() {
            if i > 0 {
                if i % hor_len == 0 {
                    result.push(hor_delim);
                } else if i % self.ru_length == 0 {
                    result.push(ru_delim);
                }
            }
            result.push(*c);
        }

        result
    }

    /// Create a `SharedChromosome` which provides an immutable, cheaply
    /// clonable view of this chromosome suitable for parallel read-only use.
    pub fn to_shared(&self) -> SharedChromosome {
        SharedChromosome {
            id: self.id.clone(),
            sequence: self.sequence.to_shared(),
            ru_length: self.ru_length,
            rus_per_hor: self.rus_per_hor,
        }
    }
}

impl std::fmt::Display for Chromosome {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.sequence)
    }
}

/// Immutable, shared view of a chromosome.
///
/// `SharedChromosome` contains a reference-counted view of the sequence data
/// and the repeat-structure metadata. It is intended for read-only access and
/// cheap cloning across threads.
#[derive(Debug, Clone)]
pub struct SharedChromosome {
    id: Arc<str>,
    sequence: SharedSequence,
    ru_length: usize,
    rus_per_hor: usize,
}

impl SharedChromosome {
    /// Return the chromosome identifier.
    #[inline]
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Borrow the shared, immutable sequence.
    #[inline]
    pub fn sequence(&self) -> &SharedSequence {
        &self.sequence
    }

    /// Return the length of the shared chromosome in bases.
    #[inline]
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// Return true if the shared chromosome contains no bases.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// Get the repeat unit (RU) length in bases.
    #[inline]
    pub fn ru_length(&self) -> usize {
        self.ru_length
    }

    /// Get the number of RUs per HOR.
    #[inline]
    pub fn rus_per_hor(&self) -> usize {
        self.rus_per_hor
    }

    /// Calculate GC content for the shared chromosome.
    pub fn gc_content(&self) -> f64 {
        let mut gc_count = 0;
        let total = self.sequence.len();

        if total == 0 {
            return 0.0;
        }

        for &nuc in self.sequence.as_slice() {
            if matches!(nuc, Nucleotide::G | Nucleotide::C) {
                gc_count += 1;
            }
        }

        gc_count as f64 / total as f64
    }

    /// Convert the shared chromosome into an owned `Chromosome` with
    /// mutable sequence data (clones the underlying indices if necessary).
    pub fn to_mutable(&self) -> Chromosome {
        Chromosome {
            id: self.id.clone(),
            sequence: self.sequence.to_mutable(),
            ru_length: self.ru_length,
            rus_per_hor: self.rus_per_hor,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::str::FromStr;

    fn test_sequence() -> Sequence {
        Sequence::from_str("ACGTACGTACGTACGT").unwrap()
    }

    // ===== Chromosome Tests =====

    #[test]
    fn test_chromosome_new() {
        let seq = test_sequence();
        let chr = Chromosome::new("chr1", seq.clone(), 4, 2);

        assert_eq!(chr.id(), "chr1");
        assert_eq!(chr.len(), 16);
        assert_eq!(chr.ru_length(), 4);
        assert_eq!(chr.rus_per_hor(), 2);
    }

    #[test]
    fn test_chromosome_uniform() {
        let chr = Chromosome::uniform("chr1", Nucleotide::A, 100, 10, 5);

        assert_eq!(chr.len(), 100);
        assert_eq!(chr.to_string(), "A".repeat(100));
        assert_eq!(chr.ru_length(), 10);
        assert_eq!(chr.rus_per_hor(), 5);
    }

    #[test]
    fn test_chromosome_id() {
        let chr = Chromosome::uniform("test_id", Nucleotide::A, 10, 5, 2);
        assert_eq!(chr.id(), "test_id");
    }

    #[test]
    fn test_chromosome_id_shared() {
        let chr1 = Chromosome::uniform("chr1", Nucleotide::A, 10, 5, 2);
        let chr2 = chr1.clone();

        // Both should share the same ID Arc
        assert_eq!(chr1.id(), chr2.id());
    }

    #[test]
    fn test_chromosome_sequence() {
        let seq = test_sequence();
        let chr = Chromosome::new("chr1", seq.clone(), 4, 2);

        assert_eq!(chr.sequence().to_string(), "ACGTACGTACGTACGT");
    }

    #[test]
    fn test_chromosome_sequence_mut() {
        let seq = test_sequence();
        let mut chr = Chromosome::new("chr1", seq, 4, 2);

        chr.sequence_mut().set(0, Nucleotide::T).unwrap();
        assert_eq!(chr.sequence().to_string().chars().next(), Some('T'));
    }

    #[test]
    fn test_chromosome_len() {
        let chr = Chromosome::uniform("chr1", Nucleotide::A, 100, 10, 5);
        assert_eq!(chr.len(), 100);
    }

    #[test]
    fn test_chromosome_is_empty() {
        let empty_seq = Sequence::new();
        let chr = Chromosome::new("chr1", empty_seq, 4, 2);
        assert!(chr.is_empty());

        let non_empty = Chromosome::uniform("chr2", Nucleotide::A, 10, 5, 2);
        assert!(!non_empty.is_empty());
    }

    #[test]
    fn test_chromosome_hor_length() {
        let chr = Chromosome::uniform("chr1", Nucleotide::A, 100, 10, 5);
        assert_eq!(chr.hor_length(), 50); // 10 * 5
    }

    #[test]
    fn test_chromosome_num_hors() {
        let chr = Chromosome::uniform("chr1", Nucleotide::A, 100, 10, 5);
        assert_eq!(chr.num_hors(), 2); // 100 / 50
    }

    #[test]
    fn test_chromosome_num_hors_incomplete() {
        let chr = Chromosome::uniform("chr1", Nucleotide::A, 75, 10, 5);
        assert_eq!(chr.num_hors(), 1); // 75 / 50 = 1 (integer division)
    }

    #[test]
    fn test_chromosome_gc_content_all_gc() {
        let seq = Sequence::from_str("GCGCGCGC").unwrap();
        let chr = Chromosome::new("chr1", seq, 2, 2);
        assert_eq!(chr.gc_content(), 1.0);
    }

    #[test]
    fn test_chromosome_gc_content_all_at() {
        let seq = Sequence::from_str("ATATATATAT").unwrap();
        let chr = Chromosome::new("chr1", seq, 2, 2);
        assert_eq!(chr.gc_content(), 0.0);
    }

    #[test]
    fn test_chromosome_gc_content_half() {
        let seq = Sequence::from_str("ACGT").unwrap();
        let chr = Chromosome::new("chr1", seq, 2, 2);
        assert_eq!(chr.gc_content(), 0.5);
    }

    #[test]
    fn test_chromosome_gc_content_empty() {
        let seq = Sequence::new();
        let chr = Chromosome::new("chr1", seq, 2, 2);
        assert_eq!(chr.gc_content(), 0.0);
    }

    #[test]
    fn test_chromosome_to_string() {
        let seq = Sequence::from_str("ACGT").unwrap();
        let chr = Chromosome::new("chr1", seq, 2, 2);
        assert_eq!(chr.to_string(), "ACGT");
    }

    #[test]
    fn test_chromosome_to_formatted_string_basic() {
        let seq = Sequence::from_str("ACGTACGTACGTACGT").unwrap();
        let chr = Chromosome::new("chr1", seq, 4, 2);

        let formatted = chr.to_formatted_string('|', '#');
        // RU length = 4, HOR length = 8
        // Expected: ACGT|ACGT#ACGT|ACGT
        assert_eq!(formatted, "ACGT|ACGT#ACGT|ACGT");
    }

    #[test]
    fn test_chromosome_to_formatted_string_single_ru() {
        let seq = Sequence::from_str("ACGT").unwrap();
        let chr = Chromosome::new("chr1", seq, 4, 1);

        let formatted = chr.to_formatted_string('|', '#');
        assert_eq!(formatted, "ACGT"); // No delimiters for single RU
    }

    #[test]
    fn test_chromosome_to_formatted_string_custom_delimiters() {
        let seq = Sequence::from_str("ACGTACGT").unwrap();
        let chr = Chromosome::new("chr1", seq, 2, 2);

        let formatted = chr.to_formatted_string('-', '=');
        // RU length = 2, HOR length = 4
        // Expected: AC-GT=AC-GT
        assert_eq!(formatted, "AC-GT=AC-GT");
    }

    #[test]
    fn test_chromosome_clone() {
        let chr1 = Chromosome::uniform("chr1", Nucleotide::A, 100, 10, 5);
        let chr2 = chr1.clone();

        assert_eq!(chr1.id(), chr2.id());
        assert_eq!(chr1.len(), chr2.len());
        assert_eq!(chr1.ru_length(), chr2.ru_length());
        assert_eq!(chr1.to_string(), chr2.to_string());
    }

    // ===== SharedChromosome Tests =====

    #[test]
    fn test_chromosome_to_shared() {
        let chr = Chromosome::uniform("chr1", Nucleotide::A, 100, 10, 5);
        let shared = chr.to_shared();

        assert_eq!(shared.id(), "chr1");
        assert_eq!(shared.len(), 100);
        assert_eq!(shared.ru_length(), 10);
        assert_eq!(shared.rus_per_hor(), 5);
    }

    #[test]
    fn test_shared_chromosome_id() {
        let chr = Chromosome::uniform("test_id", Nucleotide::A, 10, 5, 2);
        let shared = chr.to_shared();
        assert_eq!(shared.id(), "test_id");
    }

    #[test]
    fn test_shared_chromosome_sequence() {
        let seq = test_sequence();
        let chr = Chromosome::new("chr1", seq, 4, 2);
        let shared = chr.to_shared();

        assert_eq!(shared.sequence().len(), 16);
    }

    #[test]
    fn test_shared_chromosome_len() {
        let chr = Chromosome::uniform("chr1", Nucleotide::A, 100, 10, 5);
        let shared = chr.to_shared();
        assert_eq!(shared.len(), 100);
    }

    #[test]
    fn test_shared_chromosome_is_empty() {
        let empty_seq = Sequence::new();
        let chr = Chromosome::new("chr1", empty_seq, 4, 2);
        let shared = chr.to_shared();
        assert!(shared.is_empty());
    }

    #[test]
    fn test_shared_chromosome_gc_content() {
        let seq = Sequence::from_str("GCGCGCGC").unwrap();
        let chr = Chromosome::new("chr1", seq, 2, 2);
        let shared = chr.to_shared();
        assert_eq!(shared.gc_content(), 1.0);
    }

    #[test]
    fn test_shared_chromosome_clone_is_cheap() {
        let chr = Chromosome::uniform("chr1", Nucleotide::A, 1000, 10, 5);
        let shared1 = chr.to_shared();
        let shared2 = shared1.clone();

        // Both should share the same sequence data
        assert_eq!(shared1.sequence().strong_count(), 2);
        assert_eq!(shared2.sequence().strong_count(), 2);
    }

    #[test]
    fn test_shared_chromosome_to_mutable() {
        let chr1 = Chromosome::uniform("chr1", Nucleotide::A, 100, 10, 5);
        let shared = chr1.to_shared();
        let chr2 = shared.to_mutable();

        assert_eq!(chr2.id(), "chr1");
        assert_eq!(chr2.len(), 100);
        assert_eq!(chr2.ru_length(), 10);
        assert_eq!(chr2.to_string(), chr1.to_string());
    }

    #[test]
    fn test_roundtrip_mutable_to_shared_to_mutable() {
        let chr1 = Chromosome::uniform("chr1", Nucleotide::A, 100, 10, 5);
        let original_str = chr1.to_string();

        let shared = chr1.to_shared();
        let chr2 = shared.to_mutable();

        assert_eq!(chr2.to_string(), original_str);
        assert_eq!(chr2.id(), "chr1");
    }

    #[test]
    fn test_shared_chromosome_immutability() {
        let chr = Chromosome::uniform("chr1", Nucleotide::A, 100, 10, 5);
        let shared = chr.to_shared();

        // Clone it
        let _cloned = shared.clone();

        // Original shared should be unchanged (immutable)
        assert_eq!(shared.len(), 100);
    }

    // ===== Integration Tests =====

    #[test]
    fn test_chromosome_complex_structure() {
        // Create a chromosome with realistic parameters
        // RU = 171 bp, 12 RUs per HOR = 2052 bp HOR
        let chr = Chromosome::uniform(
            "chr1",
            Nucleotide::A,
            20520, // 10 HORs
            171,
            12,
        );

        assert_eq!(chr.hor_length(), 2052);
        assert_eq!(chr.num_hors(), 10);
    }

    #[test]
    fn test_chromosome_mutation_scenario() {
        let mut chr = Chromosome::uniform("chr1", Nucleotide::A, 100, 10, 5);

        // Mutate some bases
        chr.sequence_mut().set(0, Nucleotide::T).unwrap();
        chr.sequence_mut().set(50, Nucleotide::G).unwrap();

        // Check that mutations were applied
        assert_eq!(chr.sequence().get(0), Some(Nucleotide::T));
        assert_eq!(chr.sequence().get(50), Some(Nucleotide::G));
    }

    #[test]
    fn test_chromosome_different_ru_configurations() {
        let chr1 = Chromosome::uniform("chr1", Nucleotide::A, 100, 10, 5);
        let chr2 = Chromosome::uniform("chr2", Nucleotide::A, 100, 5, 10);
        let chr3 = Chromosome::uniform("chr3", Nucleotide::A, 100, 20, 2);

        // All same length but different HOR structures
        assert_eq!(chr1.len(), chr2.len());
        assert_eq!(chr2.len(), chr3.len());

        // But different HOR lengths
        assert_eq!(chr1.hor_length(), 50);
        assert_eq!(chr2.hor_length(), 50);
        assert_eq!(chr3.hor_length(), 40);
    }

    #[test]
    fn test_chromosome_large() {
        // Test with larger, more realistic sizes
        let chr = Chromosome::uniform(
            "chr1",
            Nucleotide::A,
            1_000_000, // 1 Mbp
            171,
            12,
        );

        assert_eq!(chr.len(), 1_000_000);
        assert!(chr.num_hors() > 0);
    }
}
