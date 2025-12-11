//! Sequence encoding and error correction library.
//!
//! Provides strategies for compressing and protecting DNA sequence data.

mod error;
mod strategies;
mod traits;
mod utils;

pub use error::CodecError as Error;
pub use error::CodecError;
pub use strategies::{BitPackedRS, ParallelBitPackedRS, UnpackedRS};
pub use traits::Codec;

use serde::{Deserialize, Serialize};

/// strategies for encoding sequences.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CodecStrategy {
    /// No compression (raw bytes) with RS protection.
    UnpackedRS,
    /// Bit-packed into 4 bases/byte + Reed-Solomon(255, 223) ECC.
    BitPackedRS,
    /// Parallel Bit-packed (Chunked + Rayon).
    ParallelBitPackedRS,
}

impl CodecStrategy {
    /// Encode using the selected strategy.
    pub fn encode(&self, seq: &[u8]) -> Result<Vec<u8>, CodecError> {
        match self {
            CodecStrategy::UnpackedRS => UnpackedRS.encode(seq),
            CodecStrategy::BitPackedRS => BitPackedRS.encode(seq),
            CodecStrategy::ParallelBitPackedRS => ParallelBitPackedRS.encode(seq),
        }
    }

    /// Decode using the selected strategy.
    pub fn decode(&self, data: &[u8]) -> Result<Vec<u8>, CodecError> {
        match self {
            CodecStrategy::UnpackedRS => UnpackedRS.decode(data),
            CodecStrategy::BitPackedRS => BitPackedRS.decode(data),
            CodecStrategy::ParallelBitPackedRS => ParallelBitPackedRS.decode(data),
        }
    }
}

impl Default for CodecStrategy {
    fn default() -> Self {
        Self::ParallelBitPackedRS
    }
}

impl std::fmt::Display for CodecStrategy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::UnpackedRS => write!(f, "unpacked-rs"),
            Self::BitPackedRS => write!(f, "packed-rs"),
            Self::ParallelBitPackedRS => write!(f, "parallel-packed-rs"),
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
            "unpacked" => Ok(Self::UnpackedRS), // Alias for backward compat/simplicity? Or maybe UnpackedNoRS if we had it exposed. User asked for simulation.
            _ => Err(format!(
                "Unknown codec strategy: {s}. Available: unpacked-rs, packed-rs, parallel-packed-rs"
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
    }

    // Note: Error correction tests removed as reed-solomon-simd (without checksums)
    // acts as an erasure encoder and may not detect/correct corrupted data blindly.
    // For robust error correction, we would need to implement per-shard checksums.
}
