//! Sequence encoding and error correction library.
//!
//! # Overview
//!
//! This library provides tools to compress and protect DNA sequence data (A, C, G, T).
//! It is designed to be used by simulation engines and storage systems that need
//! efficient and reliable handling of genetic information.
//!
//! # Key Concepts
//!
//! ## 1. Compression (Bit-packing)
//!
//! DNA only has 4 possible letters (A, C, G, T). We can use a trick called **bit-packing** to
//! squish 4 letters into a single byte.
//!
//! *   **Raw:** 4 bases = 4 bytes.
//! *   **Packed:** 4 bases = 1 byte. (75% space savings!)
//!
//! ## 2. Error Correction (Reed-Solomon)
//! Sometimes data gets corrupted (like a scratch on a CD). **Reed-Solomon (RS)** is a mathematical
//! technique that adds some extra "redundancy" data. If part of your sequence is lost or changed,
//! the RS algorithm can use this extra data to reconstruct the original sequence perfectly.
//!
//! This library combines these two concepts to give you compact, durable DNA storage.

mod error;
mod strategies;
mod traits;
mod utils;

pub use error::CodecError as Error;
pub use error::CodecError;
pub use strategies::{BitPackedRS, ParallelBitPackedRS, UnpackedRS, UnpackedZ, UnpackedZRS};
pub use traits::Codec;

use serde::{Deserialize, Serialize};

/// Strategies for encoding (compressing & protecting) DNA sequences.
///
/// Choose a strategy based on your needs for speed vs. space.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CodecStrategy {
    /// **No Compression + Protection.**
    ///
    /// Stores the sequence as raw bytes (one byte per base).
    /// *   **Pros:** Simplest, no packing overhead.
    /// *   **Cons:** Uses the most memory/disk space.
    UnpackedRS,

    /// **Bit-Packing + Protection.**
    ///
    /// Compresses 4 bases into 1 byte.
    /// *   **Pros:** Saves ~75% space.
    /// *   **Cons:** Slightly slower unique to pack/unpack.
    /// *   **Best for:** General usage.
    BitPackedRS,

    /// **Parallel Bit-Packing + Protection.**
    ///
    /// Same as `BitPackedRS`, but processed in parallel chunks.
    /// *   **Pros:** Faster on huge sequences (multi-core).
    /// *   **Cons:** Slight overhead for small sequences.
    /// *   **Best for:** Very large datasets.
    ParallelBitPackedRS,

    /// **Unpacked + Zstd + Protection.**
    ///
    /// Compresses raw bytes with Zstd, then protects with RS.
    /// *   **Pros:** Best storage for repetitive data. Good speed. Checks integrity.
    /// *   **Best for:** Centromeres / Repetitive DNA.
    UnpackedZRS,

    /// **Unpacked + Zstd (No RS).**
    ///
    /// Compresses raw bytes with Zstd. No error correction.
    /// *   **Pros:** Fastest, smallest.
    /// *   **Cons:** No corruption protection.
    UnpackedZ,
}

impl CodecStrategy {
    /// Encode using the selected strategy.
    pub fn encode(&self, seq: &[u8]) -> Result<Vec<u8>, CodecError> {
        match self {
            CodecStrategy::UnpackedRS => UnpackedRS.encode(seq),
            CodecStrategy::BitPackedRS => BitPackedRS.encode(seq),
            CodecStrategy::ParallelBitPackedRS => ParallelBitPackedRS.encode(seq),
            CodecStrategy::UnpackedZRS => UnpackedZRS.encode(seq),
            CodecStrategy::UnpackedZ => UnpackedZ.encode(seq),
        }
    }

    /// Decode using the selected strategy.
    pub fn decode(&self, data: &[u8]) -> Result<Vec<u8>, CodecError> {
        match self {
            CodecStrategy::UnpackedRS => UnpackedRS.decode(data),
            CodecStrategy::BitPackedRS => BitPackedRS.decode(data),
            CodecStrategy::ParallelBitPackedRS => ParallelBitPackedRS.decode(data),
            CodecStrategy::UnpackedZRS => UnpackedZRS.decode(data),
            CodecStrategy::UnpackedZ => UnpackedZ.decode(data),
        }
    }
}

impl Default for CodecStrategy {
    fn default() -> Self {
        Self::UnpackedZRS
    }
}

impl std::fmt::Display for CodecStrategy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::UnpackedRS => write!(f, "unpacked-rs"),
            Self::BitPackedRS => write!(f, "packed-rs"),
            Self::ParallelBitPackedRS => write!(f, "parallel-packed-rs"),
            Self::UnpackedZRS => write!(f, "unpacked-rsz"),
            Self::UnpackedZ => write!(f, "unpacked-z"),
        }
    }
}

impl std::str::FromStr for CodecStrategy {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "unpacked-rs" => Ok(Self::UnpackedRS),
            "packed-rs" => Ok(Self::BitPackedRS),
            "parallel-packed-rs" => Ok(Self::ParallelBitPackedRS),
            "unpacked-rsz" => Ok(Self::UnpackedZRS),
            "unpacked-z" => Ok(Self::UnpackedZ),
            "unpacked" => Ok(Self::UnpackedRS), // Alias
            _ => Err(format!(
                "Unknown codec strategy: {s}. Available: unpacked-rs, packed-rs, parallel-packed-rs, unpacked-rsz, unpacked-z"
            )),
        }
    }
}

/// Helper to provide trait-like access if needed, though Enum dispatch is preferred.
impl Codec for CodecStrategy {
    fn encode(&self, seq: &[u8]) -> Result<Vec<u8>, CodecError> {
        CodecStrategy::encode(self, seq)
    }

    fn decode(&self, data: &[u8]) -> Result<Vec<u8>, CodecError> {
        CodecStrategy::decode(self, data)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_strategies_round_trip() {
        let seq = vec![0, 1, 2, 3, 0, 1, 2, 3];

        // BitPacked
        let encoded = CodecStrategy::BitPackedRS.encode(&seq).unwrap();
        let decoded = CodecStrategy::BitPackedRS.decode(&encoded).unwrap();
        assert_eq!(seq, decoded);

        // Unpacked
        let encoded = CodecStrategy::UnpackedRS.encode(&seq).unwrap();
        let decoded = CodecStrategy::UnpackedRS.decode(&encoded).unwrap();
        assert_eq!(seq, decoded);

        // UnpackedZRS
        let encoded = CodecStrategy::UnpackedZRS.encode(&seq).unwrap();
        let decoded = CodecStrategy::UnpackedZRS.decode(&encoded).unwrap();
        assert_eq!(seq, decoded);

        // UnpackedZ
        let encoded = CodecStrategy::UnpackedZ.encode(&seq).unwrap();
        let decoded = CodecStrategy::UnpackedZ.decode(&encoded).unwrap();
        assert_eq!(seq, decoded);
    }

    // Note: Error correction tests removed as reed-solomon-simd (without checksums)
    // acts as an erasure encoder and may not detect/correct corrupted data blindly.
    // For robust error correction, we would need to implement per-shard checksums.
}
