use std::sync::Arc;
use crate::base::{Alphabet, Sequence, SharedSequence, Nucleotide};

/// A chromosome is a sequence organized into repeat units and higher-order repeats.
#[derive(Debug, Clone)]
pub struct Chromosome {
    /// Unique identifier (shared across clones)
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

    /// Create a uniform chromosome (all same base)
    pub fn uniform(
        id: impl Into<Arc<str>>,
        base: Nucleotide,
        length: usize,
        ru_length: usize,
        rus_per_hor: usize,
        alphabet: Alphabet,
    ) -> Self {
        let mut sequence = Sequence::with_capacity(length, alphabet);
        for _ in 0..length {
            sequence.push(base);
        }

        Self::new(id, sequence, ru_length, rus_per_hor)
    }

    /// Get chromosome ID
    #[inline]
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Get sequence reference
    #[inline]
    pub fn sequence(&self) -> &Sequence {
        &self.sequence
    }

    /// Get mutable sequence reference
    #[inline]
    pub fn sequence_mut(&mut self) -> &mut Sequence {
        &mut self.sequence
    }

    /// Get length
    #[inline]
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// Check if empty
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// Get repeat unit length
    #[inline]
    pub fn ru_length(&self) -> usize {
        self.ru_length
    }

    /// Get RUs per HOR
    #[inline]
    pub fn rus_per_hor(&self) -> usize {
        self.rus_per_hor
    }

    /// Get HOR length in bases
    #[inline]
    pub fn hor_length(&self) -> usize {
        self.ru_length * self.rus_per_hor
    }

    /// Get number of complete HORs
    #[inline]
    pub fn num_hors(&self) -> usize {
        self.len() / self.hor_length()
    }

    /// Get alphabet
    #[inline]
    pub fn alphabet(&self) -> &Alphabet {
        self.sequence.alphabet()
    }

    /// Calculate GC content
    pub fn gc_content(&self) -> f64 {
        let mut gc_count = 0;
        let mut total = 0;

        for &idx in self.sequence.indices() {
            if let Some(nuc) = Nucleotide::from_index(idx) {
                total += 1;
                if matches!(nuc, Nucleotide::G | Nucleotide::C) {
                    gc_count += 1;
                }
            }
        }

        if total == 0 {
            0.0
        } else {
            gc_count as f64 / total as f64
        }
    }

    /// Convert to formatted string with RU and HOR delimiters
    pub fn to_formatted_string(&self, ru_delim: char, hor_delim: char) -> String {
        let chars: Vec<char> = self.sequence
            .indices()
            .iter()
            .filter_map(|&idx| self.alphabet().get_char(idx))
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

    /// Create a shared immutable view
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

/// Immutable shared chromosome - use for parallel operations
#[derive(Debug, Clone)]
pub struct SharedChromosome {
    id: Arc<str>,
    sequence: SharedSequence,
    ru_length: usize,
    rus_per_hor: usize,
}

impl SharedChromosome {
    /// Get chromosome ID
    #[inline]
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Get sequence reference
    #[inline]
    pub fn sequence(&self) -> &SharedSequence {
        &self.sequence
    }

    /// Get length
    #[inline]
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// Check if empty
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// Get repeat unit length
    #[inline]
    pub fn ru_length(&self) -> usize {
        self.ru_length
    }

    /// Get RUs per HOR
    #[inline]
    pub fn rus_per_hor(&self) -> usize {
        self.rus_per_hor
    }

    /// Calculate GC content
    pub fn gc_content(&self) -> f64 {
        let mut gc_count = 0;
        let mut total = 0;

        for &idx in self.sequence.indices() {
            if let Some(nuc) = Nucleotide::from_index(idx) {
                total += 1;
                if matches!(nuc, Nucleotide::G | Nucleotide::C) {
                    gc_count += 1;
                }
            }
        }

        if total == 0 {
            0.0
        } else {
            gc_count as f64 / total as f64
        }
    }

    /// Clone data into mutable chromosome
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

    fn test_alphabet() -> Alphabet {
        Alphabet::dna()
    }

    fn test_sequence() -> Sequence {
        Sequence::from_str("ACGTACGTACGTACGT", test_alphabet()).unwrap()
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
        let chr = Chromosome::uniform(
            "chr1",
            Nucleotide::A,
            100,
            10,
            5,
            test_alphabet(),
        );
        
        assert_eq!(chr.len(), 100);
        assert_eq!(chr.to_string(), "A".repeat(100));
        assert_eq!(chr.ru_length(), 10);
        assert_eq!(chr.rus_per_hor(), 5);
    }

    #[test]
    fn test_chromosome_id() {
        let chr = Chromosome::uniform("test_id", Nucleotide::A, 10, 5, 2, test_alphabet());
        assert_eq!(chr.id(), "test_id");
    }

    #[test]
    fn test_chromosome_id_shared() {
        let chr1 = Chromosome::uniform("chr1", Nucleotide::A, 10, 5, 2, test_alphabet());
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
        let chr = Chromosome::uniform("chr1", Nucleotide::A, 100, 10, 5, test_alphabet());
        assert_eq!(chr.len(), 100);
    }

    #[test]
    fn test_chromosome_is_empty() {
        let empty_seq = Sequence::new(test_alphabet());
        let chr = Chromosome::new("chr1", empty_seq, 4, 2);
        assert!(chr.is_empty());
        
        let non_empty = Chromosome::uniform("chr2", Nucleotide::A, 10, 5, 2, test_alphabet());
        assert!(!non_empty.is_empty());
    }

    #[test]
    fn test_chromosome_hor_length() {
        let chr = Chromosome::uniform("chr1", Nucleotide::A, 100, 10, 5, test_alphabet());
        assert_eq!(chr.hor_length(), 50); // 10 * 5
    }

    #[test]
    fn test_chromosome_num_hors() {
        let chr = Chromosome::uniform("chr1", Nucleotide::A, 100, 10, 5, test_alphabet());
        assert_eq!(chr.num_hors(), 2); // 100 / 50
    }

    #[test]
    fn test_chromosome_num_hors_incomplete() {
        let chr = Chromosome::uniform("chr1", Nucleotide::A, 75, 10, 5, test_alphabet());
        assert_eq!(chr.num_hors(), 1); // 75 / 50 = 1 (integer division)
    }

    #[test]
    fn test_chromosome_alphabet() {
        let alphabet = test_alphabet();
        let chr = Chromosome::uniform("chr1", Nucleotide::A, 10, 5, 2, alphabet.clone());
        assert_eq!(chr.alphabet(), &alphabet);
    }

    #[test]
    fn test_chromosome_gc_content_all_gc() {
        let seq = Sequence::from_str("GCGCGCGC", test_alphabet()).unwrap();
        let chr = Chromosome::new("chr1", seq, 2, 2);
        assert_eq!(chr.gc_content(), 1.0);
    }

    #[test]
    fn test_chromosome_gc_content_all_at() {
        let seq = Sequence::from_str("ATATATATAT", test_alphabet()).unwrap();
        let chr = Chromosome::new("chr1", seq, 2, 2);
        assert_eq!(chr.gc_content(), 0.0);
    }

    #[test]
    fn test_chromosome_gc_content_half() {
        let seq = Sequence::from_str("ACGT", test_alphabet()).unwrap();
        let chr = Chromosome::new("chr1", seq, 2, 2);
        assert_eq!(chr.gc_content(), 0.5);
    }

    #[test]
    fn test_chromosome_gc_content_empty() {
        let seq = Sequence::new(test_alphabet());
        let chr = Chromosome::new("chr1", seq, 2, 2);
        assert_eq!(chr.gc_content(), 0.0);
    }

    #[test]
    fn test_chromosome_to_string() {
        let seq = Sequence::from_str("ACGT", test_alphabet()).unwrap();
        let chr = Chromosome::new("chr1", seq, 2, 2);
        assert_eq!(chr.to_string(), "ACGT");
    }

    #[test]
    fn test_chromosome_to_formatted_string_basic() {
        let seq = Sequence::from_str("ACGTACGTACGTACGT", test_alphabet()).unwrap();
        let chr = Chromosome::new("chr1", seq, 4, 2);
        
        let formatted = chr.to_formatted_string('|', '#');
        // RU length = 4, HOR length = 8
        // Expected: ACGT|ACGT#ACGT|ACGT
        assert_eq!(formatted, "ACGT|ACGT#ACGT|ACGT");
    }

    #[test]
    fn test_chromosome_to_formatted_string_single_ru() {
        let seq = Sequence::from_str("ACGT", test_alphabet()).unwrap();
        let chr = Chromosome::new("chr1", seq, 4, 1);
        
        let formatted = chr.to_formatted_string('|', '#');
        assert_eq!(formatted, "ACGT"); // No delimiters for single RU
    }

    #[test]
    fn test_chromosome_to_formatted_string_custom_delimiters() {
        let seq = Sequence::from_str("ACGTACGT", test_alphabet()).unwrap();
        let chr = Chromosome::new("chr1", seq, 2, 2);
        
        let formatted = chr.to_formatted_string('-', '=');
        // RU length = 2, HOR length = 4
        // Expected: AC-GT=AC-GT
        assert_eq!(formatted, "AC-GT=AC-GT");
    }

    #[test]
    fn test_chromosome_clone() {
        let chr1 = Chromosome::uniform("chr1", Nucleotide::A, 100, 10, 5, test_alphabet());
        let chr2 = chr1.clone();
        
        assert_eq!(chr1.id(), chr2.id());
        assert_eq!(chr1.len(), chr2.len());
        assert_eq!(chr1.ru_length(), chr2.ru_length());
        assert_eq!(chr1.to_string(), chr2.to_string());
    }

    // ===== SharedChromosome Tests =====

    #[test]
    fn test_chromosome_to_shared() {
        let chr = Chromosome::uniform("chr1", Nucleotide::A, 100, 10, 5, test_alphabet());
        let shared = chr.to_shared();
        
        assert_eq!(shared.id(), "chr1");
        assert_eq!(shared.len(), 100);
        assert_eq!(shared.ru_length(), 10);
        assert_eq!(shared.rus_per_hor(), 5);
    }

    #[test]
    fn test_shared_chromosome_id() {
        let chr = Chromosome::uniform("test_id", Nucleotide::A, 10, 5, 2, test_alphabet());
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
        let chr = Chromosome::uniform("chr1", Nucleotide::A, 100, 10, 5, test_alphabet());
        let shared = chr.to_shared();
        assert_eq!(shared.len(), 100);
    }

    #[test]
    fn test_shared_chromosome_is_empty() {
        let empty_seq = Sequence::new(test_alphabet());
        let chr = Chromosome::new("chr1", empty_seq, 4, 2);
        let shared = chr.to_shared();
        assert!(shared.is_empty());
    }

    #[test]
    fn test_shared_chromosome_gc_content() {
        let seq = Sequence::from_str("GCGCGCGC", test_alphabet()).unwrap();
        let chr = Chromosome::new("chr1", seq, 2, 2);
        let shared = chr.to_shared();
        assert_eq!(shared.gc_content(), 1.0);
    }

    #[test]
    fn test_shared_chromosome_clone_is_cheap() {
        let chr = Chromosome::uniform("chr1", Nucleotide::A, 1000, 10, 5, test_alphabet());
        let shared1 = chr.to_shared();
        let shared2 = shared1.clone();
        
        // Both should share the same sequence data
        assert_eq!(shared1.sequence().strong_count(), 2);
        assert_eq!(shared2.sequence().strong_count(), 2);
    }

    #[test]
    fn test_shared_chromosome_to_mutable() {
        let chr1 = Chromosome::uniform("chr1", Nucleotide::A, 100, 10, 5, test_alphabet());
        let shared = chr1.to_shared();
        let chr2 = shared.to_mutable();
        
        assert_eq!(chr2.id(), "chr1");
        assert_eq!(chr2.len(), 100);
        assert_eq!(chr2.ru_length(), 10);
        assert_eq!(chr2.to_string(), chr1.to_string());
    }

    #[test]
    fn test_roundtrip_mutable_to_shared_to_mutable() {
        let chr1 = Chromosome::uniform("chr1", Nucleotide::A, 100, 10, 5, test_alphabet());
        let original_str = chr1.to_string();
        
        let shared = chr1.to_shared();
        let chr2 = shared.to_mutable();
        
        assert_eq!(chr2.to_string(), original_str);
        assert_eq!(chr2.id(), "chr1");
    }

    #[test]
    fn test_shared_chromosome_immutability() {
        let chr = Chromosome::uniform("chr1", Nucleotide::A, 100, 10, 5, test_alphabet());
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
            test_alphabet(),
        );
        
        assert_eq!(chr.hor_length(), 2052);
        assert_eq!(chr.num_hors(), 10);
    }

    #[test]
    fn test_chromosome_mutation_scenario() {
        let mut chr = Chromosome::uniform("chr1", Nucleotide::A, 100, 10, 5, test_alphabet());
        
        // Mutate some bases
        chr.sequence_mut().set(0, Nucleotide::T).unwrap();
        chr.sequence_mut().set(50, Nucleotide::G).unwrap();
        
        // Check that mutations were applied
        assert_eq!(chr.sequence().get(0), Some(Nucleotide::T));
        assert_eq!(chr.sequence().get(50), Some(Nucleotide::G));
    }

    #[test]
    fn test_chromosome_different_ru_configurations() {
        let chr1 = Chromosome::uniform("chr1", Nucleotide::A, 100, 10, 5, test_alphabet());
        let chr2 = Chromosome::uniform("chr2", Nucleotide::A, 100, 5, 10, test_alphabet());
        let chr3 = Chromosome::uniform("chr3", Nucleotide::A, 100, 20, 2, test_alphabet());
        
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
            test_alphabet(),
        );
        
        assert_eq!(chr.len(), 1_000_000);
        assert!(chr.num_hors() > 0);
    }
}
