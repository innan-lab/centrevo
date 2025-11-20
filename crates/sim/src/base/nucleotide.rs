use core::fmt;

use serde::{Serialize, Deserialize};
use crate::errors::InvalidNucleotide;

/// A DNA nucleotide base.
///
/// `Nucleotide` is a compact, Copyable representation of DNA bases backed by
/// a single byte (u8). The mapping of variants to integers is stable and used
/// throughout the crate (A=0, C=1, G=2, T=3). Use the convenience conversion
/// functions to go between bytes/chars and `Nucleotide`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[repr(u8)]
pub enum Nucleotide {
    A = 0,
    C = 1,
    G = 2,
    T = 3,
}

impl Nucleotide {
    /// Convert from u8 index (0-3)
    #[inline(always)]
    pub const fn from_index(idx: u8) -> Option<Self> {
        match idx {
            0 => Some(Self::A),
            1 => Some(Self::C),
            2 => Some(Self::G),
            3 => Some(Self::T),
            _ => None,
        }
    }

    /// Convert to the compact u8 index (0-3).
    #[inline(always)]
    pub const fn to_index(self) -> u8 {
        self as u8
    }

    /// Convert from an ASCII byte (`b'A'`, `b'C'`, `b'G'`, `b'T'`) and also
    /// accepts lowercase bytes. Returns `None` for non-standard characters.
    #[inline]
    pub const fn from_ascii(byte: u8) -> Option<Self> {
        match byte {
            b'A' | b'a' => Some(Self::A),
            b'C' | b'c' => Some(Self::C),
            b'G' | b'g' => Some(Self::G),
            b'T' | b't' => Some(Self::T),
            _ => None,
        }
    }

    /// Convert to an uppercase ASCII byte representing this nucleotide.
    #[inline(always)]
    pub const fn to_ascii(self) -> u8 {
        match self {
            Self::A => b'A',
            Self::C => b'C',
            Self::G => b'G',
            Self::T => b'T',
        }
    }

    /// Convert to an uppercase `char` representing this nucleotide.
    #[inline(always)]
    pub const fn to_char(self) -> char {
        self.to_ascii() as char
    }

    /// Return the complementary base (A <-> T, C <-> G).
    #[inline(always)]
    pub const fn complement(self) -> Self {
        match self {
            Self::A => Self::T,
            Self::T => Self::A,
            Self::C => Self::G,
            Self::G => Self::C,
        }
    }

    /// Return true if the nucleotide is a purine (A or G).
    #[inline(always)]
    pub const fn is_purine(self) -> bool {
        matches!(self, Self::A | Self::G)
    }

    /// Return true if the nucleotide is a pyrimidine (C or T).
    #[inline(always)]
    pub const fn is_pyrimidine(self) -> bool {
        matches!(self, Self::C | Self::T)
    }
}

impl TryFrom<u8> for Nucleotide {
    type Error = InvalidNucleotide;

    fn try_from(byte: u8) -> Result<Self, Self::Error> {
        Self::from_ascii(byte).ok_or(InvalidNucleotide(byte))
    }
}

impl From<Nucleotide> for u8 {
    #[inline(always)]
    fn from(nuc: Nucleotide) -> u8 {
        nuc.to_index()
    }
}

impl From<Nucleotide> for char {
    #[inline(always)]
    fn from(nuc: Nucleotide) -> char {
        nuc.to_char()
    }
}

impl fmt::Display for Nucleotide {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_char())
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nucleotide_from_index() {
        assert_eq!(Nucleotide::from_index(0), Some(Nucleotide::A));
        assert_eq!(Nucleotide::from_index(1), Some(Nucleotide::C));
        assert_eq!(Nucleotide::from_index(2), Some(Nucleotide::G));
        assert_eq!(Nucleotide::from_index(3), Some(Nucleotide::T));
        assert_eq!(Nucleotide::from_index(4), None);
        assert_eq!(Nucleotide::from_index(255), None);
    }

    #[test]
    fn test_nucleotide_to_index() {
        assert_eq!(Nucleotide::A.to_index(), 0);
        assert_eq!(Nucleotide::C.to_index(), 1);
        assert_eq!(Nucleotide::G.to_index(), 2);
        assert_eq!(Nucleotide::T.to_index(), 3);
    }

    #[test]
    fn test_nucleotide_from_ascii() {
        // Uppercase
        assert_eq!(Nucleotide::from_ascii(b'A'), Some(Nucleotide::A));
        assert_eq!(Nucleotide::from_ascii(b'C'), Some(Nucleotide::C));
        assert_eq!(Nucleotide::from_ascii(b'G'), Some(Nucleotide::G));
        assert_eq!(Nucleotide::from_ascii(b'T'), Some(Nucleotide::T));

        // Lowercase
        assert_eq!(Nucleotide::from_ascii(b'a'), Some(Nucleotide::A));
        assert_eq!(Nucleotide::from_ascii(b'c'), Some(Nucleotide::C));
        assert_eq!(Nucleotide::from_ascii(b'g'), Some(Nucleotide::G));
        assert_eq!(Nucleotide::from_ascii(b't'), Some(Nucleotide::T));

        // Invalid
        assert_eq!(Nucleotide::from_ascii(b'N'), None);
        assert_eq!(Nucleotide::from_ascii(b'X'), None);
        assert_eq!(Nucleotide::from_ascii(b'5'), None);
        assert_eq!(Nucleotide::from_ascii(b' '), None);
    }

    #[test]
    fn test_nucleotide_to_ascii() {
        assert_eq!(Nucleotide::A.to_ascii(), b'A');
        assert_eq!(Nucleotide::C.to_ascii(), b'C');
        assert_eq!(Nucleotide::G.to_ascii(), b'G');
        assert_eq!(Nucleotide::T.to_ascii(), b'T');
    }

    #[test]
    fn test_nucleotide_to_char() {
        assert_eq!(Nucleotide::A.to_char(), 'A');
        assert_eq!(Nucleotide::C.to_char(), 'C');
        assert_eq!(Nucleotide::G.to_char(), 'G');
        assert_eq!(Nucleotide::T.to_char(), 'T');
    }

    #[test]
    fn test_nucleotide_complement() {
        assert_eq!(Nucleotide::A.complement(), Nucleotide::T);
        assert_eq!(Nucleotide::T.complement(), Nucleotide::A);
        assert_eq!(Nucleotide::C.complement(), Nucleotide::G);
        assert_eq!(Nucleotide::G.complement(), Nucleotide::C);

        // Double complement returns original
        assert_eq!(Nucleotide::A.complement().complement(), Nucleotide::A);
        assert_eq!(Nucleotide::C.complement().complement(), Nucleotide::C);
    }

    #[test]
    fn test_nucleotide_is_purine() {
        assert!(Nucleotide::A.is_purine());
        assert!(!Nucleotide::C.is_purine());
        assert!(Nucleotide::G.is_purine());
        assert!(!Nucleotide::T.is_purine());
    }

    #[test]
    fn test_nucleotide_is_pyrimidine() {
        assert!(!Nucleotide::A.is_pyrimidine());
        assert!(Nucleotide::C.is_pyrimidine());
        assert!(!Nucleotide::G.is_pyrimidine());
        assert!(Nucleotide::T.is_pyrimidine());
    }

    #[test]
    fn test_nucleotide_try_from_u8() {
        assert_eq!(Nucleotide::try_from(b'A'), Ok(Nucleotide::A));
        assert_eq!(Nucleotide::try_from(b'c'), Ok(Nucleotide::C));
        assert!(Nucleotide::try_from(b'N').is_err());
        
        let err = Nucleotide::try_from(b'X').unwrap_err();
        assert_eq!(err.0, b'X');
    }

    #[test]
    fn test_nucleotide_into_u8() {
        let idx: u8 = Nucleotide::A.into();
        assert_eq!(idx, 0);
        
        let idx: u8 = Nucleotide::T.into();
        assert_eq!(idx, 3);
    }

    #[test]
    fn test_nucleotide_into_char() {
        let c: char = Nucleotide::A.into();
        assert_eq!(c, 'A');
        
        let c: char = Nucleotide::G.into();
        assert_eq!(c, 'G');
    }

    #[test]
    fn test_invalid_nucleotide_display() {
        let err = InvalidNucleotide(b'X');
        let msg = format!("{err}");
        assert!(msg.contains("Invalid"));
        assert!(msg.contains("88")); // ASCII value of 'X'
        assert!(msg.contains("X"));
    }

    #[test]
    fn test_nucleotide_equality() {
        assert_eq!(Nucleotide::A, Nucleotide::A);
        assert_ne!(Nucleotide::A, Nucleotide::C);
        
        // Test copy semantics
        let n1 = Nucleotide::G;
        let n2 = n1;
        assert_eq!(n1, n2);
    }

    #[test]
    fn test_nucleotide_hash() {
        use std::collections::HashSet;
        
        let mut set = HashSet::new();
        set.insert(Nucleotide::A);
        set.insert(Nucleotide::C);
        set.insert(Nucleotide::A); // Duplicate
        
        assert_eq!(set.len(), 2);
        assert!(set.contains(&Nucleotide::A));
        assert!(set.contains(&Nucleotide::C));
        assert!(!set.contains(&Nucleotide::T));
    }

    #[test]
    fn test_nucleotide_size() {
        // Ensure Nucleotide is exactly 1 byte
        assert_eq!(std::mem::size_of::<Nucleotide>(), 1);
    }
}
