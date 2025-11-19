use std::fmt;
use std::error;

/// Error returned when attempting to convert an invalid byte/character into
/// a `Nucleotide`.
///
/// The inner `u8` is the original byte that failed to parse. This type
/// implements `error::Error` and `Display` to provide helpful messages
/// when surfaced to callers or upstream libraries.
///
/// Example:
///
/// ```rust
/// # use centrevo::Nucleotide;
/// let err = Nucleotide::try_from(b'X').unwrap_err(); 
/// println!("{err}");
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct InvalidNucleotide(pub u8);

impl fmt::Display for InvalidNucleotide {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Invalid nucleotide byte: {} ('{}')", self.0, self.0 as char)
    }
}

impl error::Error for InvalidNucleotide {}

/// Error type for failures when constructing or manipulating a `Sequence`.
#[derive(Debug, Clone)]
pub enum InvalidSequence {
    /// A character was not recognized as a valid nucleotide.
    InvalidChar(char),

    /// The sequence was empty when a non-empty sequence was required.
    EmptySequence,
}

impl fmt::Display for InvalidSequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidChar(c) => write!(f, "Invalid character in sequence: '{c}'"),
            Self::EmptySequence => write!(f, "Empty sequence not allowed"),
        }
    }
}

impl error::Error for InvalidSequence {}

/// Error returned when an index is outside the valid range for a sequence.
#[derive(Debug, Clone, Copy)]
pub struct OutOfBounds {
    /// The index that was requested
    pub index: usize,

    /// The current length of the sequence (upper bound)
    pub len: usize,
}

impl fmt::Display for OutOfBounds {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Index {} out of bounds (len = {})", self.index, self.len)
    }
}

impl error::Error for OutOfBounds {}

