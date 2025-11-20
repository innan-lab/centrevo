use super::Nucleotide;
use std::fmt;
use std::sync::Arc;
use std::str::FromStr;
use crate::errors::{InvalidSequence, OutOfBounds};

/// Mutable biological sequence backed by a vector of Nucleotides.
///
/// `Sequence` stores a vector of `Nucleotide`s. It is intended for active,
/// in-place operations such as mutation, insertion, and deletion. For read-only,
/// shareable views convert to `SharedSequence` using `to_shared` or `into_shared`.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Sequence(Vec<Nucleotide>);

impl Sequence {
    /// Create a new, empty `Sequence`.
    ///
    /// Example:
    ///
    /// ```rust
    /// # use centrevo::base::{Sequence, Nucleotide};
    /// let seq = Sequence::new();
    /// assert_eq!(seq.len(), 0);
    /// ```
    pub fn new() -> Self {
        Self(Vec::new())
    }

    /// Create a `Sequence` with reserved capacity for `capacity` bases.
    pub fn with_capacity(capacity: usize) -> Self {
        Self(Vec::with_capacity(capacity))
    }

    /// Create a `Sequence` from a vector of `Nucleotide`s.
    pub fn from_nucleotides(nucleotides: Vec<Nucleotide>) -> Self {
        Self(nucleotides)
    }

    /// Create a `Sequence` from a vector of indices (0-3).
    /// Indices outside 0-3 are treated as A (0).
    pub fn from_indices(indices: Vec<u8>) -> Self {
        let nucleotides = indices
            .into_iter()
            .map(|i| Nucleotide::from_index(i).unwrap_or(Nucleotide::A))
            .collect();
        Self(nucleotides)
    }

    /// Return the length of the sequence in bases.
    #[inline(always)]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Return `true` if the sequence contains no bases.
    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Return the `Nucleotide` at `index`, or `None` if out of range.
    #[inline]
    pub fn get(&self, index: usize) -> Option<Nucleotide> {
        self.0.get(index).copied()
    }

    /// Set the base at `index` to `base`.
    ///
    /// Returns `OutOfBounds` if `index` is greater than or equal to the
    /// sequence length.
    #[inline]
    pub fn set(&mut self, index: usize, base: Nucleotide) -> Result<(), OutOfBounds> {
        self.0
            .get_mut(index)
            .map(|slot| *slot = base)
            .ok_or(OutOfBounds {
                index,
                len: self.len(),
            })
    }

    /// Borrow the underlying `Nucleotide` slice.
    #[inline]
    pub fn as_slice(&self) -> &[Nucleotide] {
        &self.0
    }

    /// Borrow the mutable underlying `Nucleotide` slice.
    #[inline]
    pub fn as_mut_slice(&mut self) -> &mut [Nucleotide] {
        &mut self.0
    }

    /// Append `base` to the end of the sequence.
    #[inline]
    pub fn push(&mut self, base: Nucleotide) {
        self.0.push(base);
    }

    /// Remove and return the last base, or `None` if the sequence is empty.
    #[inline]
    pub fn pop(&mut self) -> Option<Nucleotide> {
        self.0.pop()
    }

    /// Insert `base` at position `index`, shifting subsequent elements.
    #[inline]
    pub fn insert(&mut self, index: usize, base: Nucleotide) {
        self.0.insert(index, base);
    }

    /// Remove and return the base at `index`.
    ///
    /// Panics if `index` is out of bounds (matching the behavior of
    /// `Vec::remove`).
    #[inline]
    pub fn remove(&mut self, index: usize) -> Nucleotide {
        self.0.remove(index)
    }

    /// Consume this `Sequence` and produce an immutable `SharedSequence`.
    ///
    /// This avoids cloning the internal buffer by reusing the owned `Vec`
    /// storage as an `Arc<[Nucleotide]>` where possible.
    pub fn into_shared(self) -> SharedSequence {
        SharedSequence(self.0.into())
    }

    /// Create an immutable `SharedSequence` by cloning the internal data.
    ///
    /// Use this when you need a shared, read-only view but still keep the
    /// original mutable `Sequence`.
    pub fn to_shared(&self) -> SharedSequence {
        SharedSequence(self.0.clone().into())
    }
}

impl Default for Sequence {
    fn default() -> Self {
        Self::new()
    }
}

impl fmt::Display for Sequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for &nuc in &self.0 {
            write!(f, "{}", nuc.to_char())?;
        }
        Ok(())
    }
}

impl FromStr for Sequence {
    type Err = InvalidSequence;

    /// Parse a textual representation (e.g. "ACGT") into a `Sequence`.
    ///
    /// Characters not present in the standard DNA alphabet produce an
    /// `InvalidSequence` error. This function is case-insensitive for ASCII letters.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let data: Result<Vec<Nucleotide>, _> = s
            .chars()
            .map(|c| Nucleotide::from_ascii(c as u8).ok_or(InvalidSequence::InvalidChar(c)))
            .collect();

        Ok(Self(data?))
    }
}

/// Immutable, shareable sequence view.
///
/// `SharedSequence` holds its data in a reference-counted `Arc<[Nucleotide]>`.
/// Cloning a `SharedSequence` is cheap and the structure is safe to share
/// across threads for read-only operations.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct SharedSequence(Arc<[Nucleotide]>);

impl SharedSequence {
    /// Create a `SharedSequence` from an `Arc<[Nucleotide]>`.
    pub fn new(data: Arc<[Nucleotide]>) -> Self {
        Self(data)
    }

    /// Return the number of bases in the shared sequence.
    #[inline(always)]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Return true if the shared sequence has no bases.
    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Return the `Nucleotide` at `index`, or `None` if out of bounds.
    #[inline]
    pub fn get(&self, index: usize) -> Option<Nucleotide> {
        self.0.get(index).copied()
    }

    /// Borrow the `Nucleotide` slice.
    #[inline]
    pub fn as_slice(&self) -> &[Nucleotide] {
        &self.0
    }

    /// Clone the shared data into a new mutable `Sequence`.
    ///
    /// This performs a copy of the underlying data into owned `Vec<Nucleotide>` so
    /// the returned `Sequence` can be mutated independently of the shared
    /// view.
    pub fn to_mutable(&self) -> Sequence {
        Sequence(self.0.to_vec())
    }

    /// Return the current strong reference count to the shared data (useful
    /// for debugging or assertions about copying behavior).
    pub fn strong_count(&self) -> usize {
        Arc::strong_count(&self.0)
    }
}

impl fmt::Display for SharedSequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for &nuc in self.0.iter() {
            write!(f, "{}", nuc.to_char())?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ===== Sequence Tests =====

    #[test]
    fn test_sequence_new() {
        let seq = Sequence::new();
        assert_eq!(seq.len(), 0);
        assert!(seq.is_empty());
    }

    #[test]
    fn test_sequence_with_capacity() {
        let seq = Sequence::with_capacity(100);
        assert_eq!(seq.len(), 0);
        assert!(seq.is_empty());
        // Capacity is at least 100 (actual capacity might be larger)
        assert!(seq.0.capacity() >= 100);
    }

    #[test]
    fn test_sequence_from_nucleotides() {
        let nucs = vec![Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T];
        let seq = Sequence::from_nucleotides(nucs);
        assert_eq!(seq.len(), 4);
        assert_eq!(seq.get(0), Some(Nucleotide::A));
        assert_eq!(seq.get(1), Some(Nucleotide::C));
        assert_eq!(seq.get(2), Some(Nucleotide::G));
        assert_eq!(seq.get(3), Some(Nucleotide::T));
    }

    #[test]
    fn test_sequence_from_indices() {
        let indices = vec![0, 1, 2, 3]; // A, C, G, T
        let seq = Sequence::from_indices(indices);
        assert_eq!(seq.len(), 4);
        assert_eq!(seq.get(0), Some(Nucleotide::A));
        assert_eq!(seq.get(1), Some(Nucleotide::C));
        assert_eq!(seq.get(2), Some(Nucleotide::G));
        assert_eq!(seq.get(3), Some(Nucleotide::T));
    }

    #[test]
    fn test_sequence_from_indices_with_invalid() {
        let indices = vec![0, 1, 4, 3]; // A, C, A(default), T
        let seq = Sequence::from_indices(indices);
        assert_eq!(seq.len(), 4);
        assert_eq!(seq.get(0), Some(Nucleotide::A));
        assert_eq!(seq.get(1), Some(Nucleotide::C));
        assert_eq!(seq.get(2), Some(Nucleotide::A));
        assert_eq!(seq.get(3), Some(Nucleotide::T));
    }

    #[test]
    fn test_sequence_from_str_valid() {
        let seq = Sequence::from_str("ACGT").unwrap();
        assert_eq!(seq.len(), 4);
        assert_eq!(seq.to_string(), "ACGT");
    }

    #[test]
    fn test_sequence_from_str_lowercase() {
        let seq = Sequence::from_str("acgt").unwrap();
        assert_eq!(seq.len(), 4);
        assert_eq!(seq.to_string(), "ACGT");
    }

    #[test]
    fn test_sequence_from_str_mixed_case() {
        let seq = Sequence::from_str("AcGt").unwrap();
        assert_eq!(seq.len(), 4);
        assert_eq!(seq.to_string(), "ACGT");
    }

    #[test]
    fn test_sequence_from_str_invalid() {
        let result = Sequence::from_str("ACGN");
        assert!(result.is_err());

        match result.unwrap_err() {
            InvalidSequence::InvalidChar(c) => assert_eq!(c, 'N'),
            _ => panic!("Expected InvalidChar error"),
        }
    }

    #[test]
    fn test_sequence_from_str_empty() {
        let seq = Sequence::from_str("").unwrap();
        assert_eq!(seq.len(), 0);
        assert!(seq.is_empty());
    }

    #[test]
    fn test_sequence_get() {
        let seq = Sequence::from_str("ACGT").unwrap();
        assert_eq!(seq.get(0), Some(Nucleotide::A));
        assert_eq!(seq.get(1), Some(Nucleotide::C));
        assert_eq!(seq.get(2), Some(Nucleotide::G));
        assert_eq!(seq.get(3), Some(Nucleotide::T));
        assert_eq!(seq.get(4), None);
    }

    #[test]
    fn test_sequence_set() {
        let mut seq = Sequence::from_str("ACGT").unwrap();
        seq.set(1, Nucleotide::T).unwrap();
        assert_eq!(seq.to_string(), "ATGT");
    }

    #[test]
    fn test_sequence_set_out_of_bounds() {
        let mut seq = Sequence::from_str("ACGT").unwrap();
        let result = seq.set(10, Nucleotide::A);
        assert!(result.is_err());

        let err = result.unwrap_err();
        assert_eq!(err.index, 10);
        assert_eq!(err.len, 4);
    }

    #[test]
    fn test_sequence_push() {
        let mut seq = Sequence::new();
        seq.push(Nucleotide::A);
        seq.push(Nucleotide::C);
        seq.push(Nucleotide::G);
        seq.push(Nucleotide::T);

        assert_eq!(seq.len(), 4);
        assert_eq!(seq.to_string(), "ACGT");
    }

    #[test]
    fn test_sequence_pop() {
        let mut seq = Sequence::from_str("ACGT").unwrap();

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
        let mut seq = Sequence::from_str("ACT").unwrap();
        seq.insert(2, Nucleotide::G);
        assert_eq!(seq.to_string(), "ACGT");
    }

    #[test]
    fn test_sequence_insert_at_beginning() {
        let mut seq = Sequence::from_str("CGT").unwrap();
        seq.insert(0, Nucleotide::A);
        assert_eq!(seq.to_string(), "ACGT");
    }

    #[test]
    fn test_sequence_insert_at_end() {
        let mut seq = Sequence::from_str("ACG").unwrap();
        seq.insert(3, Nucleotide::T);
        assert_eq!(seq.to_string(), "ACGT");
    }

    #[test]
    fn test_sequence_remove() {
        let mut seq = Sequence::from_str("ACGT").unwrap();
        let removed = seq.remove(1);
        assert_eq!(removed, Nucleotide::C);
        assert_eq!(seq.to_string(), "AGT");
    }

    #[test]
    #[should_panic]
    fn test_sequence_remove_out_of_bounds() {
        let mut seq = Sequence::from_str("ACGT").unwrap();
        seq.remove(10);
    }

    #[test]
    fn test_sequence_as_slice() {
        let seq = Sequence::from_str("ACGT").unwrap();
        let slice = seq.as_slice();
        assert_eq!(
            slice,
            &[Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T]
        );
    }

    #[test]
    fn test_sequence_as_mut_slice() {
        let mut seq = Sequence::from_str("ACGT").unwrap();
        let slice = seq.as_mut_slice();
        slice[1] = Nucleotide::T; // Change C to T
        assert_eq!(seq.to_string(), "ATGT");
    }

    #[test]
    fn test_sequence_clone() {
        let seq1 = Sequence::from_str("ACGT").unwrap();
        let seq2 = seq1.clone();

        assert_eq!(seq1, seq2);
        assert_eq!(seq1.to_string(), seq2.to_string());
    }

    #[test]
    fn test_sequence_equality() {
        let seq1 = Sequence::from_str("ACGT").unwrap();
        let seq2 = Sequence::from_str("ACGT").unwrap();
        let seq3 = Sequence::from_str("TGCA").unwrap();

        assert_eq!(seq1, seq2);
        assert_ne!(seq1, seq3);
    }

    #[test]
    fn test_sequence_to_string_long() {
        let seq = Sequence::from_str("ACGTACGTACGT").unwrap();
        assert_eq!(seq.to_string(), "ACGTACGTACGT");
    }

    // ===== SharedSequence Tests =====

    #[test]
    fn test_sequence_to_shared() {
        let seq = Sequence::from_str("ACGT").unwrap();
        let shared = seq.to_shared();

        assert_eq!(shared.len(), 4);
        assert_eq!(shared.to_string(), "ACGT");
    }

    #[test]
    fn test_sequence_into_shared() {
        let seq = Sequence::from_str("ACGT").unwrap();
        let shared = seq.into_shared();

        assert_eq!(shared.len(), 4);
        assert_eq!(shared.to_string(), "ACGT");
    }

    #[test]
    fn test_shared_sequence_new() {
        let data: Arc<[Nucleotide]> =
            vec![Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T].into();
        let shared = SharedSequence::new(data);

        assert_eq!(shared.len(), 4);
        assert_eq!(shared.to_string(), "ACGT");
    }

    #[test]
    fn test_shared_sequence_get() {
        let seq = Sequence::from_str("ACGT").unwrap();
        let shared = seq.to_shared();

        assert_eq!(shared.get(0), Some(Nucleotide::A));
        assert_eq!(shared.get(1), Some(Nucleotide::C));
        assert_eq!(shared.get(2), Some(Nucleotide::G));
        assert_eq!(shared.get(3), Some(Nucleotide::T));
        assert_eq!(shared.get(4), None);
    }

    #[test]
    fn test_shared_sequence_as_slice() {
        let seq = Sequence::from_str("ACGT").unwrap();
        let shared = seq.to_shared();

        assert_eq!(
            shared.as_slice(),
            &[Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T]
        );
    }

    #[test]
    fn test_shared_sequence_clone_is_cheap() {
        let seq = Sequence::from_str("ACGT").unwrap();
        let shared1 = seq.to_shared();
        let shared2 = shared1.clone();

        // Verify they share the same Arc
        assert!(Arc::ptr_eq(&shared1.0, &shared2.0));
        assert_eq!(shared1.strong_count(), 2);
    }

    #[test]
    fn test_shared_sequence_strong_count() {
        let seq = Sequence::from_str("ACGT").unwrap();
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
        let seq = Sequence::from_str("ACGT").unwrap();
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
        let seq = Sequence::from_str("ACGT").unwrap();
        let shared1 = seq.to_shared();
        let shared2 = shared1.clone();

        assert_eq!(shared1, shared2);
    }

    #[test]
    fn test_shared_sequence_equality_different_arc() {
        let seq1 = Sequence::from_str("ACGT").unwrap();
        let seq2 = Sequence::from_str("ACGT").unwrap();
        let shared1 = seq1.to_shared();
        let shared2 = seq2.to_shared();

        assert_eq!(shared1, shared2);
    }

    #[test]
    fn test_shared_sequence_inequality() {
        let seq1 = Sequence::from_str("ACGT").unwrap();
        let seq2 = Sequence::from_str("TGCA").unwrap();
        let shared1 = seq1.to_shared();
        let shared2 = seq2.to_shared();

        assert_ne!(shared1, shared2);
    }

    #[test]
    fn test_shared_sequence_empty() {
        let seq = Sequence::new();
        let shared = seq.to_shared();

        assert!(shared.is_empty());
        assert_eq!(shared.len(), 0);
    }

    // ===== Conversion Tests =====

    #[test]
    fn test_roundtrip_mutable_to_shared_to_mutable() {
        let seq1 = Sequence::from_str("ACGT").unwrap();
        let shared = seq1.to_shared();
        let seq2 = shared.to_mutable();

        assert_eq!(seq1, seq2);
    }

    #[test]
    fn test_into_shared_consumes_original() {
        let seq = Sequence::from_str("ACGT").unwrap();
        let _shared = seq.into_shared();
        // seq is now consumed, cannot use it anymore
    }

    // ===== Error Tests =====

    #[test]
    fn test_invalid_sequence_display() {
        let err = InvalidSequence::InvalidChar('N');
        let msg = format!("{err}");
        assert!(msg.contains("Invalid character"));
        assert!(msg.contains("N"));
    }

    #[test]
    fn test_invalid_sequence_empty_display() {
        let err = InvalidSequence::EmptySequence;
        let msg = format!("{err}");
        assert!(msg.contains("Empty sequence"));
    }

    #[test]
    fn test_out_of_bounds_display() {
        let err = OutOfBounds { index: 10, len: 5 };
        let msg = format!("{err}");
        assert!(msg.contains("10"));
        assert!(msg.contains("5"));
        assert!(msg.contains("out of bounds"));
    }

    // ===== Performance-oriented Tests =====

    #[test]
    fn test_large_sequence() {
        let bases = "ACGT".repeat(1000);
        let seq = Sequence::from_str(&bases).unwrap();
        assert_eq!(seq.len(), 4000);
        assert_eq!(seq.to_string(), bases);
    }

    #[test]
    fn test_shared_sequence_no_clone_overhead() {
        let seq = Sequence::from_str(&"ACGT".repeat(1000)).unwrap();
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
    fn test_sequence_mutation_operations() {
        let mut seq = Sequence::new();

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
