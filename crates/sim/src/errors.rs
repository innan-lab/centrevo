use std::fmt;
use std::error;

/// Error returned when attempting to convert an invalid byte/character into
/// a `Nucleotide`.
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

/// Errors that can occur in fitness calculations.
#[derive(Debug, Clone, PartialEq)]
pub enum FitnessError {
    /// Invalid parameter value
    InvalidParameter(String),
}

impl fmt::Display for FitnessError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FitnessError::InvalidParameter(msg) => {
                write!(f, "Invalid fitness parameter: {msg}")
            }
        }
    }
}

impl error::Error for FitnessError {}

/// Errors that can occur during recombination operations.
#[derive(Debug, Clone, PartialEq)]
pub enum RecombinationError {
    /// Invalid probability value
    InvalidProbability(&'static str, f64),
    /// Sequences have different lengths
    LengthMismatch { len1: usize, len2: usize },
    /// Invalid position for operation
    InvalidPosition { position: usize, length: usize },
    /// Invalid range for gene conversion
    InvalidRange { start: usize, end: usize },
}

impl fmt::Display for RecombinationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            RecombinationError::InvalidProbability(name, val) => {
                write!(
                    f,
                    "Invalid probability for {name}: {val} (must be between 0.0 and 1.0)"
                )
            }
            RecombinationError::LengthMismatch { len1, len2 } => {
                write!(f, "Sequence length mismatch: {len1} vs {len2}")
            }
            RecombinationError::InvalidPosition { position, length } => {
                write!(
                    f,
                    "Invalid position {position} for sequence of length {length}"
                )
            }
            RecombinationError::InvalidRange { start, end } => {
                write!(f, "Invalid range [{start}, {end}) for gene conversion")
            }
        }
    }
}

impl error::Error for RecombinationError {}

/// Database error types.
#[derive(Debug, Clone)]
pub enum DatabaseError {
    Connection(String),
    Initialization(String),
    Transaction(String),
    Query(String),
    Insert(String),
    Close(String),
    Vacuum(String),
}

impl fmt::Display for DatabaseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Connection(e) => write!(f, "Database connection error: {e}"),
            Self::Initialization(e) => write!(f, "Database initialization error: {e}"),
            Self::Transaction(e) => write!(f, "Transaction error: {e}"),
            Self::Query(e) => write!(f, "Query error: {e}"),
            Self::Insert(e) => write!(f, "Insert error: {e}"),
            Self::Close(e) => write!(f, "Close error: {e}"),
            Self::Vacuum(e) => write!(f, "Vacuum error: {e}"),
        }
    }
}

impl error::Error for DatabaseError {}

/// Errors that can occur during simulation building.
#[derive(Debug)]
pub enum BuilderError {
    /// A required parameter is missing
    MissingRequired(&'static str),
    /// An invalid parameter value was provided
    InvalidParameter(String),
    /// Error importing sequences
    SequenceImport(String),
}

impl fmt::Display for BuilderError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingRequired(param) => {
                write!(f, "Missing required parameter: {param}")
            }
            Self::InvalidParameter(msg) => {
                write!(f, "Invalid parameter: {msg}")
            }
            Self::SequenceImport(msg) => {
                write!(f, "Failed to import sequences: {msg}")
            }
        }
    }
}

impl error::Error for BuilderError {}

/// Errors that can occur during mutation operations.
#[derive(Debug, Clone, PartialEq)]
pub enum MutationError {
    /// Invalid mutation rate (must be between 0.0 and 1.0)
    InvalidMutationRate(f64),
}

impl fmt::Display for MutationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MutationError::InvalidMutationRate(mu) => {
                write!(
                    f,
                    "Invalid mutation rate: {mu} (must be between 0.0 and 1.0)"
                )
            }
        }
    }
}

impl error::Error for MutationError {}

/// Error types for sequence initialization.
#[derive(Debug)]
pub enum InitializationError {
    /// IO error
    Io(std::io::Error),
    /// Parse error
    Parse(String),
    /// Validation error
    Validation(String),
    /// Database error
    Database(String),
}

impl fmt::Display for InitializationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(e) => write!(f, "IO error: {e}"),
            Self::Parse(msg) => write!(f, "Parse error: {msg}"),
            Self::Validation(msg) => write!(f, "Validation error: {msg}"),
            Self::Database(msg) => write!(f, "Database error: {msg}"),
        }
    }
}

impl error::Error for InitializationError {}

impl From<std::io::Error> for InitializationError {
    fn from(e: std::io::Error) -> Self {
        Self::Io(e)
    }
}

impl From<serde_json::Error> for InitializationError {
    fn from(e: serde_json::Error) -> Self {
        Self::Parse(format!("JSON error: {e}"))
    }
}

/// Errors related to RepeatMap operations
#[derive(Debug, Clone, PartialEq)]
pub enum RepeatMapError {
    InvalidOffsets,
    InvalidHorOffsets,
    IndexOutOfBounds,
    SplitError,
}

impl fmt::Display for RepeatMapError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            RepeatMapError::InvalidOffsets => write!(f, "Invalid RU offsets (must be sorted and start with 0)"),
            RepeatMapError::InvalidHorOffsets => write!(f, "Invalid HOR offsets (must be sorted and within RU range)"),
            RepeatMapError::IndexOutOfBounds => write!(f, "Index out of bounds"),
            RepeatMapError::SplitError => write!(f, "Error splitting map"),
        }
    }
}

impl error::Error for RepeatMapError {}

/// Errors that can occur during Chromosome construction or manipulation.
#[derive(Debug, Clone, PartialEq)]
pub enum ChromosomeError {
    /// Invalid nucleotide byte encountered
    InvalidNucleotide(InvalidNucleotide),
    /// Error creating the repeat map
    RepeatMapError(RepeatMapError),
    /// Generic construction error
    Construction(String),
}

impl fmt::Display for ChromosomeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidNucleotide(e) => write!(f, "Invalid nucleotide: {e}"),
            Self::RepeatMapError(e) => write!(f, "Repeat map error: {e}"),
            Self::Construction(msg) => write!(f, "Chromosome construction error: {msg}"),
        }
    }
}

impl error::Error for ChromosomeError {}

impl From<InvalidNucleotide> for ChromosomeError {
    fn from(e: InvalidNucleotide) -> Self {
        Self::InvalidNucleotide(e)
    }
}

impl From<RepeatMapError> for ChromosomeError {
    fn from(e: RepeatMapError) -> Self {
        Self::RepeatMapError(e)
    }
}
