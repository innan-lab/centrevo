use std::sync::Arc;
use std::collections::HashMap;
use serde::{Serialize, Deserialize, Serializer, Deserializer};

/// Shared, immutable alphabet.
///
/// `Alphabet` defines the mapping between characters and compact indices used
/// to store sequences efficiently. It is intentionally immutable and
/// shareable via `Arc` so a single `Alphabet` instance can be reused across
/// many sequences or chromosomes without cloning the underlying data.
#[derive(Debug, Clone)]
pub struct Alphabet {
    /// Character representation of bases
    chars: Arc<[char]>,
    /// Mapping from char to index for fast lookup
    char_to_index: Arc<HashMap<char, u8>>,
}

impl Alphabet {
    /// Create a new, case-sensitive `Alphabet` from the provided characters.
    ///
    /// The order of `chars` determines the integer index assigned to each
    /// character (starting at 0). Duplicate characters are allowed but the
    /// last occurrence wins for the index mapping.
    ///
    /// Example:
    ///
    /// ```rust
    /// # use centrevo::Alphabet;
    /// let alph = Alphabet::new(vec!['A', 'C', 'G', 'T']);
    /// assert_eq!(alph.len(), 4);
    /// ```
    pub fn new(chars: impl Into<Vec<char>>) -> Self {
        let chars: Vec<char> = chars.into();
        let char_to_index = chars
            .iter()
            .enumerate()
            .map(|(i, &c)| (c, i as u8))
            .collect();

        Self {
            chars: chars.into(),
            char_to_index: Arc::new(char_to_index),
        }
    }

    /// Standard DNA alphabet (A, C, G, T).
    ///
    /// This is the crate's default alphabet and is returned by `Default`.
    pub fn dna() -> Self {
        Self::new(vec!['A', 'C', 'G', 'T'])
    }

    /// Standard RNA alphabet (A, C, G, U).
    pub fn rna() -> Self {
        Self::new(vec!['A', 'C', 'G', 'U'])
    }

    /// Standard amino acid alphabet (20 amino acids).
    pub fn amino_acids() -> Self {
        Self::new(vec![
            'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
            'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V',
        ])
    }

    /// Return the number of distinct characters in this alphabet.
    #[inline(always)]
    pub fn len(&self) -> usize {
        self.chars.len()
    }

    /// Return true if the alphabet contains no characters.
    ///
    /// Most code expects a non-empty alphabet; creating an empty alphabet is
    /// allowed but unusual and many higher-level operations will be invalid.
    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.chars.is_empty()
    }

    /// Return the character at `index`, or `None` if out of range.
    ///
    /// Indices are the compact integer representation used by `Sequence` and
    /// `SharedSequence` storage.
    #[inline]
    pub fn get_char(&self, index: u8) -> Option<char> {
        self.chars.get(index as usize).copied()
    }

    /// Return the compact index for character `c`, or `None` if `c` is not in
    /// the alphabet.
    #[inline]
    pub fn get_index(&self, c: char) -> Option<u8> {
        self.char_to_index.get(&c).copied()
    }

    /// Return the backing character slice for this alphabet.
    ///
    /// This slice is reference-counted internally; callers should not mutate
    /// it (it is owned by the `Alphabet`).
    #[inline]
    pub fn chars(&self) -> &[char] {
        &self.chars
    }

    /// Return `true` if the alphabet contains `c`.
    ///
    /// This is a convenience wrapper around `get_index` used when quick
    /// membership checks are needed.
    #[inline]
    pub fn contains(&self, c: char) -> bool {
        self.char_to_index.contains_key(&c)
    }
}

impl Default for Alphabet {
    fn default() -> Self {
        Self::dna()
    }
}

impl PartialEq for Alphabet {
    fn eq(&self, other: &Self) -> bool {
        // Fast path: check if they point to the same Arc
        Arc::ptr_eq(&self.chars, &other.chars)
            || self.chars == other.chars
    }
}

impl Eq for Alphabet {}

// Custom serialization for Alphabet
impl Serialize for Alphabet {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        // Serialize as the character list
        self.chars.as_ref().serialize(serializer)
    }
}

impl<'de> Deserialize<'de> for Alphabet {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        // Deserialize as a Vec and create new Alphabet
        let chars = Vec::<char>::deserialize(deserializer)?;
        Ok(Alphabet::new(chars))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_alphabet_new() {
        let alphabet = Alphabet::new(vec!['A', 'C', 'G', 'T']);
        assert_eq!(alphabet.len(), 4);
        assert!(!alphabet.is_empty());
    }

    #[test]
    fn test_alphabet_dna() {
        let alphabet = Alphabet::dna();
        assert_eq!(alphabet.len(), 4);
        assert_eq!(alphabet.chars(), &['A', 'C', 'G', 'T']);
    }

    #[test]
    fn test_alphabet_rna() {
        let alphabet = Alphabet::rna();
        assert_eq!(alphabet.len(), 4);
        assert_eq!(alphabet.chars(), &['A', 'C', 'G', 'U']);
    }

    #[test]
    fn test_alphabet_default() {
        let alphabet = Alphabet::default();
        assert_eq!(alphabet, Alphabet::dna());
    }

    #[test]
    fn test_alphabet_get_char() {
        let alphabet = Alphabet::dna();
        assert_eq!(alphabet.get_char(0), Some('A'));
        assert_eq!(alphabet.get_char(1), Some('C'));
        assert_eq!(alphabet.get_char(2), Some('G'));
        assert_eq!(alphabet.get_char(3), Some('T'));
        assert_eq!(alphabet.get_char(4), None);
        assert_eq!(alphabet.get_char(255), None);
    }

    #[test]
    fn test_alphabet_get_index() {
        let alphabet = Alphabet::dna();
        assert_eq!(alphabet.get_index('A'), Some(0));
        assert_eq!(alphabet.get_index('C'), Some(1));
        assert_eq!(alphabet.get_index('G'), Some(2));
        assert_eq!(alphabet.get_index('T'), Some(3));
        assert_eq!(alphabet.get_index('N'), None);
        assert_eq!(alphabet.get_index('X'), None);
    }

    #[test]
    fn test_alphabet_contains() {
        let alphabet = Alphabet::dna();
        assert!(alphabet.contains('A'));
        assert!(alphabet.contains('C'));
        assert!(alphabet.contains('G'));
        assert!(alphabet.contains('T'));
        assert!(!alphabet.contains('N'));
        assert!(!alphabet.contains('a')); // Case sensitive
    }

    #[test]
    fn test_alphabet_custom() {
        let alphabet = Alphabet::new(vec!['0', '1']);
        assert_eq!(alphabet.len(), 2);
        assert_eq!(alphabet.get_index('0'), Some(0));
        assert_eq!(alphabet.get_index('1'), Some(1));
        assert_eq!(alphabet.get_char(0), Some('0'));
        assert_eq!(alphabet.get_char(1), Some('1'));
    }

    #[test]
    fn test_alphabet_equality_same_arc() {
        let alphabet1 = Alphabet::dna();
        let alphabet2 = alphabet1.clone();
        
        // Should be equal (same Arc pointer)
        assert_eq!(alphabet1, alphabet2);
    }

    #[test]
    fn test_alphabet_equality_different_arc() {
        let alphabet1 = Alphabet::dna();
        let alphabet2 = Alphabet::dna();
        
        // Should be equal (same content, different Arc)
        assert_eq!(alphabet1, alphabet2);
    }

    #[test]
    fn test_alphabet_inequality() {
        let dna = Alphabet::dna();
        let binary = Alphabet::new(vec!['0', '1']);
        
        assert_ne!(dna, binary);
    }

    #[test]
    fn test_alphabet_chars_slice() {
        let alphabet = Alphabet::dna();
        let chars = alphabet.chars();
        
        assert_eq!(chars.len(), 4);
        assert_eq!(chars[0], 'A');
        assert_eq!(chars[3], 'T');
    }

    #[test]
    fn test_alphabet_clone_is_cheap() {
        let alphabet1 = Alphabet::dna();
        let alphabet2 = alphabet1.clone();
        
        // Verify they share the same Arc
        assert!(Arc::ptr_eq(&alphabet1.chars, &alphabet2.chars));
        assert!(Arc::ptr_eq(&alphabet1.char_to_index, &alphabet2.char_to_index));
    }

    #[test]
    fn test_alphabet_empty() {
        let alphabet = Alphabet::new(vec![] as Vec<char>);
        assert!(alphabet.is_empty());
        assert_eq!(alphabet.len(), 0);
    }

    #[test]
    fn test_alphabet_large() {
        // Test with larger alphabet
        let chars: Vec<char> = (b'A'..=b'Z').map(|b| b as char).collect();
        let alphabet = Alphabet::new(chars.clone());
        
        assert_eq!(alphabet.len(), 26);
        assert_eq!(alphabet.get_index('A'), Some(0));
        assert_eq!(alphabet.get_index('Z'), Some(25));
        assert_eq!(alphabet.get_char(0), Some('A'));
        assert_eq!(alphabet.get_char(25), Some('Z'));
    }

    #[test]
    fn test_alphabet_ordering_matters() {
        let alphabet1 = Alphabet::new(vec!['A', 'C', 'G', 'T']);
        let alphabet2 = Alphabet::new(vec!['T', 'G', 'C', 'A']);
        
        // Different ordering means different indices
        assert_eq!(alphabet1.get_index('A'), Some(0));
        assert_eq!(alphabet2.get_index('A'), Some(3));
        
        // But alphabets themselves are still different
        assert_ne!(alphabet1, alphabet2);
    }

    #[test]
    fn test_alphabet_duplicate_chars() {
        // Test behavior with duplicate characters (later one wins)
        let alphabet = Alphabet::new(vec!['A', 'B', 'A']);
        assert_eq!(alphabet.len(), 3);
        
        // The last 'A' at index 2 should be in the map
        assert_eq!(alphabet.get_index('A'), Some(2));
    }

    #[test]
    fn test_alphabet_unicode() {
        let alphabet = Alphabet::new(vec!['α', 'β', 'γ', 'δ']);
        assert_eq!(alphabet.len(), 4);
        assert!(alphabet.contains('α'));
        assert!(alphabet.contains('δ'));
        assert_eq!(alphabet.get_index('β'), Some(1));
    }
}
