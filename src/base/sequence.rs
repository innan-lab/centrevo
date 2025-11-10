use super::{Nucleotide, Alphabet};
use std::fmt;
use std::sync::Arc;

/// Mutable sequence - use for active operations
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Sequence {
    /// Indices into alphabet (0-3 for standard DNA)
    data: Vec<u8>,
    /// Shared reference to alphabet
    alphabet: Alphabet,
}

impl Sequence {
    /// Create a new empty sequence
    pub fn new(alphabet: Alphabet) -> Self {
        Self {
            data: Vec::new(),
            alphabet,
        }
    }

    /// Create with capacity
    pub fn with_capacity(capacity: usize, alphabet: Alphabet) -> Self {
        Self {
            data: Vec::with_capacity(capacity),
            alphabet,
        }
    }

    /// Create from raw indices
    pub fn from_indices(indices: Vec<u8>, alphabet: Alphabet) -> Self {
        Self {
            data: indices,
            alphabet,
        }
    }

    /// Create from string
    pub fn from_str(s: &str, alphabet: Alphabet) -> Result<Self, InvalidSequence> {
        let data: Result<Vec<u8>, _> = s
            .chars()
            .map(|c| {
                alphabet
                    .get_index(c.to_ascii_uppercase())
                    .ok_or(InvalidSequence::InvalidChar(c))
            })
            .collect();

        Ok(Self {
            data: data?,
            alphabet,
        })
    }

    /// Get length
    #[inline(always)]
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Check if empty
    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Get base at position
    #[inline]
    pub fn get(&self, index: usize) -> Option<Nucleotide> {
        self.data.get(index).and_then(|&idx| Nucleotide::from_index(idx))
    }

    /// Set base at position
    #[inline]
    pub fn set(&mut self, index: usize, base: Nucleotide) -> Result<(), OutOfBounds> {
        self.data
            .get_mut(index)
            .map(|slot| *slot = base.to_index())
            .ok_or(OutOfBounds { index, len: self.len() })
    }

    /// Get raw indices
    #[inline]
    pub fn indices(&self) -> &[u8] {
        &self.data
    }

    /// Get mutable raw indices
    #[inline]
    pub fn indices_mut(&mut self) -> &mut [u8] {
        &mut self.data
    }

    /// Get alphabet
    #[inline]
    pub fn alphabet(&self) -> &Alphabet {
        &self.alphabet
    }

    /// Push a base
    #[inline]
    pub fn push(&mut self, base: Nucleotide) {
        self.data.push(base.to_index());
    }

    /// Pop a base
    #[inline]
    pub fn pop(&mut self) -> Option<Nucleotide> {
        self.data.pop().and_then(Nucleotide::from_index)
    }

    /// Insert a base
    #[inline]
    pub fn insert(&mut self, index: usize, base: Nucleotide) {
        self.data.insert(index, base.to_index());
    }

    /// Remove a base
    #[inline]
    pub fn remove(&mut self, index: usize) -> Nucleotide {
        Nucleotide::from_index(self.data.remove(index))
            .expect("Invalid base in sequence")
    }

    /// Convert to immutable shared sequence
    pub fn into_shared(self) -> SharedSequence {
        SharedSequence {
            data: self.data.into(),
            alphabet: self.alphabet,
        }
    }

    /// Convert to immutable shared sequence (cloning data)
    pub fn to_shared(&self) -> SharedSequence {
        SharedSequence {
            data: self.data.clone().into(),
            alphabet: self.alphabet.clone(),
        }
    }
}

impl fmt::Display for Sequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for &idx in &self.data {
            if let Some(ch) = self.alphabet.get_char(idx) {
                write!(f, "{}", ch)?;
            }
        }
        Ok(())
    }
}

/// Immutable shared sequence - use for read-only operations
#[derive(Debug, Clone)]
pub struct SharedSequence {
    /// Shared immutable indices
    data: Arc<[u8]>,
    /// Shared reference to alphabet
    alphabet: Alphabet,
}

impl SharedSequence {
    /// Create from Arc
    pub fn new(data: Arc<[u8]>, alphabet: Alphabet) -> Self {
        Self { data, alphabet }
    }

    /// Get length
    #[inline(always)]
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Check if empty
    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Get base at position
    #[inline]
    pub fn get(&self, index: usize) -> Option<Nucleotide> {
        self.data.get(index).and_then(|&idx| Nucleotide::from_index(idx))
    }

    /// Get raw indices
    #[inline]
    pub fn indices(&self) -> &[u8] {
        &self.data
    }

    /// Get alphabet
    #[inline]
    pub fn alphabet(&self) -> &Alphabet {
        &self.alphabet
    }

    /// Clone data into mutable sequence
    pub fn to_mutable(&self) -> Sequence {
        Sequence {
            data: self.data.to_vec(),
            alphabet: self.alphabet.clone(),
        }
    }

    /// Get strong reference count (for debugging)
    pub fn strong_count(&self) -> usize {
        Arc::strong_count(&self.data)
    }
}

impl PartialEq for SharedSequence {
    fn eq(&self, other: &Self) -> bool {
        // Fast path: check if same Arc
        Arc::ptr_eq(&self.data, &other.data)
            || (self.data == other.data && self.alphabet == other.alphabet)
    }
}

impl Eq for SharedSequence {}

impl fmt::Display for SharedSequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for &idx in self.data.iter() {
            if let Some(ch) = self.alphabet.get_char(idx) {
                write!(f, "{}", ch)?;
            }
        }
        Ok(())
    }
}

// Errors
#[derive(Debug, Clone)]
pub enum InvalidSequence {
    InvalidChar(char),
    EmptySequence,
}

impl std::fmt::Display for InvalidSequence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidChar(c) => write!(f, "Invalid character in sequence: '{}'", c),
            Self::EmptySequence => write!(f, "Empty sequence not allowed"),
        }
    }
}

impl std::error::Error for InvalidSequence {}

#[derive(Debug, Clone, Copy)]
pub struct OutOfBounds {
    pub index: usize,
    pub len: usize,
}

impl std::fmt::Display for OutOfBounds {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Index {} out of bounds (len = {})", self.index, self.len)
    }
}

impl std::error::Error for OutOfBounds {}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_alphabet() -> Alphabet {
        Alphabet::dna()
    }

    // ===== Sequence Tests =====

    #[test]
    fn test_sequence_new() {
        let seq = Sequence::new(test_alphabet());
        assert_eq!(seq.len(), 0);
        assert!(seq.is_empty());
    }

    #[test]
    fn test_sequence_with_capacity() {
        let seq = Sequence::with_capacity(100, test_alphabet());
        assert_eq!(seq.len(), 0);
        assert!(seq.is_empty());
        // Capacity is at least 100 (actual capacity might be larger)
        assert!(seq.data.capacity() >= 100);
    }

    #[test]
    fn test_sequence_from_indices() {
        let indices = vec![0, 1, 2, 3]; // A, C, G, T
        let seq = Sequence::from_indices(indices, test_alphabet());
        assert_eq!(seq.len(), 4);
        assert_eq!(seq.get(0), Some(Nucleotide::A));
        assert_eq!(seq.get(1), Some(Nucleotide::C));
        assert_eq!(seq.get(2), Some(Nucleotide::G));
        assert_eq!(seq.get(3), Some(Nucleotide::T));
    }

    #[test]
    fn test_sequence_from_str_valid() {
        let seq = Sequence::from_str("ACGT", test_alphabet()).unwrap();
        assert_eq!(seq.len(), 4);
        assert_eq!(seq.to_string(), "ACGT");
    }

    #[test]
    fn test_sequence_from_str_lowercase() {
        let seq = Sequence::from_str("acgt", test_alphabet()).unwrap();
        assert_eq!(seq.len(), 4);
        assert_eq!(seq.to_string(), "ACGT");
    }

    #[test]
    fn test_sequence_from_str_mixed_case() {
        let seq = Sequence::from_str("AcGt", test_alphabet()).unwrap();
        assert_eq!(seq.len(), 4);
        assert_eq!(seq.to_string(), "ACGT");
    }

    #[test]
    fn test_sequence_from_str_invalid() {
        let result = Sequence::from_str("ACGN", test_alphabet());
        assert!(result.is_err());
        
        match result.unwrap_err() {
            InvalidSequence::InvalidChar(c) => assert_eq!(c, 'N'),
            _ => panic!("Expected InvalidChar error"),
        }
    }

    #[test]
    fn test_sequence_from_str_empty() {
        let seq = Sequence::from_str("", test_alphabet()).unwrap();
        assert_eq!(seq.len(), 0);
        assert!(seq.is_empty());
    }

    #[test]
    fn test_sequence_get() {
        let seq = Sequence::from_str("ACGT", test_alphabet()).unwrap();
        assert_eq!(seq.get(0), Some(Nucleotide::A));
        assert_eq!(seq.get(1), Some(Nucleotide::C));
        assert_eq!(seq.get(2), Some(Nucleotide::G));
        assert_eq!(seq.get(3), Some(Nucleotide::T));
        assert_eq!(seq.get(4), None);
    }

    #[test]
    fn test_sequence_set() {
        let mut seq = Sequence::from_str("ACGT", test_alphabet()).unwrap();
        seq.set(1, Nucleotide::T).unwrap();
        assert_eq!(seq.to_string(), "ATGT");
    }

    #[test]
    fn test_sequence_set_out_of_bounds() {
        let mut seq = Sequence::from_str("ACGT", test_alphabet()).unwrap();
        let result = seq.set(10, Nucleotide::A);
        assert!(result.is_err());
        
        let err = result.unwrap_err();
        assert_eq!(err.index, 10);
        assert_eq!(err.len, 4);
    }

    #[test]
    fn test_sequence_push() {
        let mut seq = Sequence::new(test_alphabet());
        seq.push(Nucleotide::A);
        seq.push(Nucleotide::C);
        seq.push(Nucleotide::G);
        seq.push(Nucleotide::T);
        
        assert_eq!(seq.len(), 4);
        assert_eq!(seq.to_string(), "ACGT");
    }

    #[test]
    fn test_sequence_pop() {
        let mut seq = Sequence::from_str("ACGT", test_alphabet()).unwrap();
        
        assert_eq!(seq.pop(), Some(Nucleotide::T));
        assert_eq!(seq.len(), 3);
        assert_eq!(seq.to_string(), "ACG");
        
        assert_eq!(seq.pop(), Some(Nucleotide::G));
        assert_eq!(seq.pop(), Some(Nucleotide::C));
        assert_eq!(seq.pop(), Some(Nucleotide::A));
        assert_eq!(seq.pop(), None);
        
        assert!(seq.is_empty());
    }

    #[test]
    fn test_sequence_insert() {
        let mut seq = Sequence::from_str("ACT", test_alphabet()).unwrap();
        seq.insert(2, Nucleotide::G);
        assert_eq!(seq.to_string(), "ACGT");
    }

    #[test]
    fn test_sequence_insert_at_beginning() {
        let mut seq = Sequence::from_str("CGT", test_alphabet()).unwrap();
        seq.insert(0, Nucleotide::A);
        assert_eq!(seq.to_string(), "ACGT");
    }

    #[test]
    fn test_sequence_insert_at_end() {
        let mut seq = Sequence::from_str("ACG", test_alphabet()).unwrap();
        seq.insert(3, Nucleotide::T);
        assert_eq!(seq.to_string(), "ACGT");
    }

    #[test]
    fn test_sequence_remove() {
        let mut seq = Sequence::from_str("ACGT", test_alphabet()).unwrap();
        let removed = seq.remove(1);
        assert_eq!(removed, Nucleotide::C);
        assert_eq!(seq.to_string(), "AGT");
    }

    #[test]
    #[should_panic]
    fn test_sequence_remove_out_of_bounds() {
        let mut seq = Sequence::from_str("ACGT", test_alphabet()).unwrap();
        seq.remove(10);
    }

    #[test]
    fn test_sequence_indices() {
        let seq = Sequence::from_str("ACGT", test_alphabet()).unwrap();
        let indices = seq.indices();
        assert_eq!(indices, &[0, 1, 2, 3]);
    }

    #[test]
    fn test_sequence_indices_mut() {
        let mut seq = Sequence::from_str("ACGT", test_alphabet()).unwrap();
        let indices = seq.indices_mut();
        indices[1] = 3; // Change C to T
        assert_eq!(seq.to_string(), "ATGT");
    }

    #[test]
    fn test_sequence_alphabet() {
        let alphabet = test_alphabet();
        let seq = Sequence::new(alphabet.clone());
        assert_eq!(seq.alphabet(), &alphabet);
    }

    #[test]
    fn test_sequence_clone() {
        let seq1 = Sequence::from_str("ACGT", test_alphabet()).unwrap();
        let seq2 = seq1.clone();
        
        assert_eq!(seq1, seq2);
        assert_eq!(seq1.to_string(), seq2.to_string());
    }

    #[test]
    fn test_sequence_equality() {
        let seq1 = Sequence::from_str("ACGT", test_alphabet()).unwrap();
        let seq2 = Sequence::from_str("ACGT", test_alphabet()).unwrap();
        let seq3 = Sequence::from_str("TGCA", test_alphabet()).unwrap();
        
        assert_eq!(seq1, seq2);
        assert_ne!(seq1, seq3);
    }

    #[test]
    fn test_sequence_to_string_long() {
        let seq = Sequence::from_str("ACGTACGTACGT", test_alphabet()).unwrap();
        assert_eq!(seq.to_string(), "ACGTACGTACGT");
    }

    // ===== SharedSequence Tests =====

    #[test]
    fn test_sequence_to_shared() {
        let seq = Sequence::from_str("ACGT", test_alphabet()).unwrap();
        let shared = seq.to_shared();
        
        assert_eq!(shared.len(), 4);
        assert_eq!(shared.to_string(), "ACGT");
    }

    #[test]
    fn test_sequence_into_shared() {
        let seq = Sequence::from_str("ACGT", test_alphabet()).unwrap();
        let shared = seq.into_shared();
        
        assert_eq!(shared.len(), 4);
        assert_eq!(shared.to_string(), "ACGT");
    }

    #[test]
    fn test_shared_sequence_new() {
        let data: Arc<[u8]> = vec![0, 1, 2, 3].into();
        let shared = SharedSequence::new(data, test_alphabet());
        
        assert_eq!(shared.len(), 4);
        assert_eq!(shared.to_string(), "ACGT");
    }

    #[test]
    fn test_shared_sequence_get() {
        let seq = Sequence::from_str("ACGT", test_alphabet()).unwrap();
        let shared = seq.to_shared();
        
        assert_eq!(shared.get(0), Some(Nucleotide::A));
        assert_eq!(shared.get(1), Some(Nucleotide::C));
        assert_eq!(shared.get(2), Some(Nucleotide::G));
        assert_eq!(shared.get(3), Some(Nucleotide::T));
        assert_eq!(shared.get(4), None);
    }

    #[test]
    fn test_shared_sequence_indices() {
        let seq = Sequence::from_str("ACGT", test_alphabet()).unwrap();
        let shared = seq.to_shared();
        
        assert_eq!(shared.indices(), &[0, 1, 2, 3]);
    }

    #[test]
    fn test_shared_sequence_alphabet() {
        let alphabet = test_alphabet();
        let seq = Sequence::from_str("ACGT", alphabet.clone()).unwrap();
        let shared = seq.to_shared();
        
        assert_eq!(shared.alphabet(), &alphabet);
    }

    #[test]
    fn test_shared_sequence_clone_is_cheap() {
        let seq = Sequence::from_str("ACGT", test_alphabet()).unwrap();
        let shared1 = seq.to_shared();
        let shared2 = shared1.clone();
        
        // Verify they share the same Arc
        assert!(Arc::ptr_eq(&shared1.data, &shared2.data));
        assert_eq!(shared1.strong_count(), 2);
    }

    #[test]
    fn test_shared_sequence_strong_count() {
        let seq = Sequence::from_str("ACGT", test_alphabet()).unwrap();
        let shared1 = seq.to_shared();
        assert_eq!(shared1.strong_count(), 1);
        
        let shared2 = shared1.clone();
        assert_eq!(shared1.strong_count(), 2);
        assert_eq!(shared2.strong_count(), 2);
        
        drop(shared2);
        assert_eq!(shared1.strong_count(), 1);
    }

    #[test]
    fn test_shared_sequence_to_mutable() {
        let seq = Sequence::from_str("ACGT", test_alphabet()).unwrap();
        let shared = seq.to_shared();
        let mut mutable = shared.to_mutable();
        
        assert_eq!(mutable.to_string(), "ACGT");
        
        // Verify we can mutate
        mutable.set(1, Nucleotide::T).unwrap();
        assert_eq!(mutable.to_string(), "ATGT");
        
        // Original shared sequence unchanged
        assert_eq!(shared.to_string(), "ACGT");
    }

    #[test]
    fn test_shared_sequence_equality_same_arc() {
        let seq = Sequence::from_str("ACGT", test_alphabet()).unwrap();
        let shared1 = seq.to_shared();
        let shared2 = shared1.clone();
        
        assert_eq!(shared1, shared2);
    }

    #[test]
    fn test_shared_sequence_equality_different_arc() {
        let seq1 = Sequence::from_str("ACGT", test_alphabet()).unwrap();
        let seq2 = Sequence::from_str("ACGT", test_alphabet()).unwrap();
        let shared1 = seq1.to_shared();
        let shared2 = seq2.to_shared();
        
        assert_eq!(shared1, shared2);
    }

    #[test]
    fn test_shared_sequence_inequality() {
        let seq1 = Sequence::from_str("ACGT", test_alphabet()).unwrap();
        let seq2 = Sequence::from_str("TGCA", test_alphabet()).unwrap();
        let shared1 = seq1.to_shared();
        let shared2 = seq2.to_shared();
        
        assert_ne!(shared1, shared2);
    }

    #[test]
    fn test_shared_sequence_empty() {
        let seq = Sequence::new(test_alphabet());
        let shared = seq.to_shared();
        
        assert!(shared.is_empty());
        assert_eq!(shared.len(), 0);
    }

    // ===== Conversion Tests =====

    #[test]
    fn test_roundtrip_mutable_to_shared_to_mutable() {
        let seq1 = Sequence::from_str("ACGT", test_alphabet()).unwrap();
        let shared = seq1.to_shared();
        let seq2 = shared.to_mutable();
        
        assert_eq!(seq1, seq2);
    }

    #[test]
    fn test_into_shared_consumes_original() {
        let seq = Sequence::from_str("ACGT", test_alphabet()).unwrap();
        let _shared = seq.into_shared();
        // seq is now consumed, cannot use it anymore
    }

    // ===== Error Tests =====

    #[test]
    fn test_invalid_sequence_display() {
        let err = InvalidSequence::InvalidChar('N');
        let msg = format!("{}", err);
        assert!(msg.contains("Invalid character"));
        assert!(msg.contains("N"));
    }

    #[test]
    fn test_invalid_sequence_empty_display() {
        let err = InvalidSequence::EmptySequence;
        let msg = format!("{}", err);
        assert!(msg.contains("Empty sequence"));
    }

    #[test]
    fn test_out_of_bounds_display() {
        let err = OutOfBounds { index: 10, len: 5 };
        let msg = format!("{}", err);
        assert!(msg.contains("10"));
        assert!(msg.contains("5"));
        assert!(msg.contains("out of bounds"));
    }

    // ===== Performance-oriented Tests =====

    #[test]
    fn test_large_sequence() {
        let bases = "ACGT".repeat(1000);
        let seq = Sequence::from_str(&bases, test_alphabet()).unwrap();
        assert_eq!(seq.len(), 4000);
        assert_eq!(seq.to_string(), bases);
    }

    #[test]
    fn test_shared_sequence_no_clone_overhead() {
        let seq = Sequence::from_str(&"ACGT".repeat(1000), test_alphabet()).unwrap();
        let shared = seq.to_shared();
        
        // Multiple clones should all point to same data
        let clones: Vec<_> = (0..100).map(|_| shared.clone()).collect();
        
        assert_eq!(shared.strong_count(), 101); // Original + 100 clones
        
        // All should be equal
        for clone in &clones {
            assert_eq!(&shared, clone);
        }
    }

    #[test]
    fn test_alphabet_sharing() {
        let alphabet = test_alphabet();
        let seq1 = Sequence::new(alphabet.clone());
        let seq2 = Sequence::new(alphabet.clone());
        
        // Both sequences should use the same alphabet (verified by equality check)
        assert_eq!(seq1.alphabet(), seq2.alphabet());
    }

    #[test]
    fn test_sequence_mutation_operations() {
        let mut seq = Sequence::new(test_alphabet());
        
        // Build sequence
        for base in [Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T] {
            seq.push(base);
        }
        assert_eq!(seq.to_string(), "ACGT");
        
        // Mutate
        seq.set(0, Nucleotide::T).unwrap();
        assert_eq!(seq.to_string(), "TCGT");
        
        // Insert
        seq.insert(2, Nucleotide::A);
        assert_eq!(seq.to_string(), "TCAGT");
        
        // Remove
        seq.remove(2);
        assert_eq!(seq.to_string(), "TCGT");
        
        // Pop
        seq.pop();
        assert_eq!(seq.to_string(), "TCG");
    }
}
