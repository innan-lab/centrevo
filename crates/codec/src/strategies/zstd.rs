use crate::error::CodecError;
use crate::traits::Codec;
use crate::utils::{protect_data, unprotect_data};

/// Strategy: Unpacked + Zstd Compression (No RS).
///
/// # How it works
///
/// The raw sequence bytes are compressed using Zstd.
///
/// # Why use this?
///
/// *   **Efficiency:** Extremely efficient storage for repetitive sequences (like centromeres).
/// *   **Speed:** Fast compression and decompression.
/// *   **Constraint:** No error correction. If data is corrupted, decompression may fail.
pub struct UnpackedZ;

impl Codec for UnpackedZ {
    fn encode(&self, seq: &[u8]) -> Result<Vec<u8>, CodecError> {
        // 1. Prepare header (Original length)
        let mut data_with_header = Vec::with_capacity(8 + seq.len());
        data_with_header.extend_from_slice(&(seq.len() as u64).to_le_bytes());
        data_with_header.extend_from_slice(seq);

        // 2. Compress
        // Using level 3 as a good default balance
        zstd::bulk::compress(&data_with_header, 3)
            .map_err(|e| CodecError::Encode(format!("Zstd compression failed: {e}")))
    }

    fn decode(&self, data: &[u8]) -> Result<Vec<u8>, CodecError> {
        // 1. Decompress
        // We need to guess the size or stream decode. But zstd::bulk::decompress requires size?
        // Actually zstd::stream::decode_all works without size.
        // Or we can use a large buffer if we knew the max size.
        // Wait, for bulk::decompress we normally need the size.
        // Do we store decompressed size? Zstd frames usually contain it if enabled?
        // Let's use `decode_all` which handles dynamic resizing.
        let decompressed = zstd::stream::decode_all(std::io::Cursor::new(data))
            .map_err(|e| CodecError::Decode(format!("Zstd decompression failed: {e}")))?;

        // 2. Extract Length
        if decompressed.len() < 8 {
            return Err(CodecError::Decode("Payload too short for header".into()));
        }
        let mut len_bytes = [0u8; 8];
        len_bytes.copy_from_slice(&decompressed[0..8]);
        let original_len = u64::from_le_bytes(len_bytes) as usize;

        // 3. Extract Sequence
        if decompressed.len() < 8 + original_len {
            return Err(CodecError::Decode(
                "Payload shorter than expected sequence".into(),
            ));
        }
        let seq = decompressed[8..8 + original_len].to_vec();
        Ok(seq)
    }
}

/// Strategy: Unpacked + Zstd Compression + Reed-Solomon.
///
/// # How it works
///
/// 1.  **Compress**: Raw sequence is compressed with Zstd.
/// 2.  **Protect**: Compressed data is protected with Reed-Solomon.
///
/// # Why use this?
///
/// *   **Best of both worlds**: High compression for repetitive data AND validation/recovery.
/// *   **Order matters**: Compressing *before* RS ensures Zstd works on low-entropy data.
pub struct UnpackedZRS;

impl Codec for UnpackedZRS {
    fn encode(&self, seq: &[u8]) -> Result<Vec<u8>, CodecError> {
        // 1. Compress (via UnpackedZ logic)
        // We reuse the exact logic: Header + Seq -> Zstd
        let compressed = UnpackedZ.encode(seq)?;

        // 2. Add header for COMPRESSED length
        // We need this because protect_data adds padding (zeros) and Zstd decode_all tries to read them.
        let mut data_with_len = Vec::with_capacity(8 + compressed.len());
        data_with_len.extend_from_slice(&(compressed.len() as u64).to_le_bytes());
        data_with_len.extend_from_slice(&compressed);

        // 3. Protect
        protect_data(&data_with_len)
    }

    fn decode(&self, data: &[u8]) -> Result<Vec<u8>, CodecError> {
        // 1. Unprotect
        let unprotected = unprotect_data(data)?;

        // 2. Extract COMPRESSED Length
        if unprotected.len() < 8 {
            return Err(CodecError::Decode(
                "Payload too short for compressed header".into(),
            ));
        }
        let mut len_bytes = [0u8; 8];
        len_bytes.copy_from_slice(&unprotected[0..8]);
        let compressed_len = u64::from_le_bytes(len_bytes) as usize;

        if unprotected.len() < 8 + compressed_len {
            return Err(CodecError::Decode(
                "Payload shorter than expected compressed data".into(),
            ));
        }

        // 3. Decompress (via UnpackedZ logic)
        let compressed = &unprotected[8..8 + compressed_len];
        UnpackedZ.decode(compressed)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

    #[test]
    fn test_unpacked_z_round_trip() {
        let codec = UnpackedZ;
        let mut rng = rand::thread_rng();
        for _ in 0..5 {
            let len = rng.gen_range(100..1000);
            let input: Vec<u8> = (0..len).map(|_| rng.gen_range(0..4)).collect();
            let encoded = codec.encode(&input).expect("Encode failed");
            let decoded = codec.decode(&encoded).expect("Decode failed");
            assert_eq!(decoded, input);
        }
    }

    #[test]
    fn test_unpacked_rsz_round_trip() {
        let codec = UnpackedZRS;
        let mut rng = rand::thread_rng();
        for _ in 0..5 {
            let len = rng.gen_range(100..1000);
            let input: Vec<u8> = (0..len).map(|_| rng.gen_range(0..4)).collect();
            let encoded = codec.encode(&input).expect("Encode failed");
            let decoded = codec.decode(&encoded).expect("Decode failed");
            assert_eq!(decoded, input);
        }
    }

    #[test]
    fn test_compression_efficiency_repetitive() {
        let codec = UnpackedZ;
        // Verify that it actually compresses repetitive data driven by expectation
        let len = 10000;
        let input: Vec<u8> = vec![0; len]; // All A's
        let encoded = codec.encode(&input).expect("Encode failed");

        // Zstd should crush this. Header + overhead should be small.
        assert!(
            encoded.len() < len / 10,
            "Should achieve >10x compression on all-zeros"
        );
    }
}
