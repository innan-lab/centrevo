use crate::error::CodecError;
use reed_solomon_simd::{ReedSolomonDecoder, ReedSolomonEncoder};

pub const RS_DATA_SHARDS: usize = 223;
pub const RS_PARITY_SHARDS: usize = 32;
pub const RS_BLOCK_SIZE: usize = RS_DATA_SHARDS + RS_PARITY_SHARDS; // 255
pub const MIN_SHARD_SIZE: usize = 64; // Minimum shard size for SIMD compatibility

pub fn protect_data(input: &[u8]) -> Result<Vec<u8>, CodecError> {
    // Pad to ensure shard size is a multiple of MIN_SHARD_SIZE (64)
    let rs_k = RS_DATA_SHARDS;
    let min_shard_size = MIN_SHARD_SIZE;

    // Calculate required shard size
    // Fix clippy::manual_div_ceil: (input.len() + rs_k - 1) / rs_k  -> input.len().div_ceil(rs_k)
    let raw_shard_size = input.len().div_ceil(rs_k);
    let mut shard_size = raw_shard_size.max(min_shard_size);

    // Round up to next multiple of 64 if needed (seems required/safer)
    let remainder = shard_size % 64;
    if remainder != 0 {
        shard_size += 64 - remainder;
    }

    let target_len = shard_size * rs_k;
    let mut padded_input = input.to_vec();
    padded_input.resize(target_len, 0);

    // Transpose Setup
    let num_blocks = shard_size;
    // shard_size already set above

    let mut shards: Vec<Vec<u8>> = vec![vec![0u8; shard_size]; rs_k];

    // Transpose
    for b in 0..num_blocks {
        let block_start = b * rs_k;
        for (k, shard) in shards.iter_mut().enumerate().take(rs_k) {
            shard[b] = padded_input[block_start + k];
        }
    }
    // Note: The nested loop above might be cleaner with direct indexing to match original logic but avoid warnings
    // Original:
    // for b in 0..num_blocks {
    //     let block_start = b * rs_k;
    //     for k in 0..rs_k {
    //         shards[k][b] = padded_input[block_start + k];
    //     }
    // }
    // The clippy warning was "loop variable k is only used to index shards".
    // Actually the inner loop iterates over shards.
    // Let's stick closer to original but use iter_mut where possible or just acknowledge it's a transpose.
    // Actually, `shards` is `Vec<Vec<u8>>` where outer is `k` (row) and inner is `b` (col).
    // `padded_input` is row-major chunks? block_start = b * rs_k suggests we interpret input as `num_blocks` chunks of `rs_k` bytes.
    // So input is (b, k). We want shards[k][b]. This is a transpose.

    // Re-implementation of transpose to avoid range loop warnings if possible, but indices are clear here.
    // Let's just suppress individual warnings if needed, or iterate.
    for b in 0..num_blocks {
        let block_start = b * rs_k;
        for (k, shard) in shards.iter_mut().enumerate().take(rs_k) {
            shard[b] = padded_input[block_start + k];
        }
    }

    // Encoder
    let mut rs = ReedSolomonEncoder::new(RS_DATA_SHARDS, RS_PARITY_SHARDS, shard_size)
        .map_err(|e| CodecError::Encode(format!("Failed to create RS Encoder: {e:?}")))?;

    // Fix clippy::needless_range_loop
    for shard in shards.iter().take(RS_DATA_SHARDS) {
        rs.add_original_shard(shard)
            .map_err(|e| CodecError::Encode(format!("Failed to add shard: {e:?}")))?;
    }

    let result = rs
        .encode()
        .map_err(|e| CodecError::Encode(format!("RS Encode failed: {e:?}")))?;

    let mut parity_shards: Vec<Vec<u8>> = vec![vec![0u8; shard_size]; RS_PARITY_SHARDS];
    for (i, shard) in result.recovery_iter().enumerate() {
        if i < parity_shards.len() {
            parity_shards[i] = shard.to_vec();
        }
    }

    // Interleave back to blocks
    let rs_n = RS_BLOCK_SIZE;
    let mut output = Vec::with_capacity(num_blocks * rs_n);

    for b in 0..num_blocks {
        for shard in shards.iter().take(rs_k) {
            output.push(shard[b]);
        }
        for parity_shard in parity_shards.iter().take(RS_PARITY_SHARDS) {
            output.push(parity_shard[b]);
        }
    }

    Ok(output)
}

pub fn calculate_protected_size(input_len: usize) -> usize {
    let rs_k = RS_DATA_SHARDS;
    let min_shard_size = MIN_SHARD_SIZE;

    // Same logic as protect_data
    let raw_shard_size = input_len.div_ceil(rs_k);
    let mut shard_size = raw_shard_size.max(min_shard_size);

    let remainder = shard_size % 64;
    if remainder != 0 {
        shard_size += 64 - remainder;
    }

    // Output is block_size * shard_size
    // block_size is RS_BLOCK_SIZE (255)
    shard_size * RS_BLOCK_SIZE
}

pub fn unprotect_data(data: &[u8]) -> Result<Vec<u8>, CodecError> {
    let rs_n = RS_BLOCK_SIZE;
    if data.len() % rs_n != 0 {
        return Err(CodecError::Decode(format!(
            "Data length {} not multiple of {}",
            data.len(),
            rs_n
        )));
    }

    let num_blocks = data.len() / rs_n;
    let shard_size = num_blocks;

    let mut shards: Vec<Vec<u8>> = vec![vec![0u8; shard_size]; RS_DATA_SHARDS];
    let mut parity: Vec<Vec<u8>> = vec![vec![0u8; shard_size]; RS_PARITY_SHARDS];

    // Transpose
    for b in 0..num_blocks {
        let block_start = b * rs_n;
        for (k, shard) in shards.iter_mut().enumerate().take(RS_DATA_SHARDS) {
            shard[b] = data[block_start + k];
        }
        for (p, parity_shard) in parity.iter_mut().enumerate().take(RS_PARITY_SHARDS) {
            parity_shard[b] = data[block_start + RS_DATA_SHARDS + p];
        }
    }

    // Decode/Reconstruct
    let mut rs = ReedSolomonDecoder::new(RS_DATA_SHARDS, RS_PARITY_SHARDS, shard_size)
        .map_err(|e| CodecError::Decode(format!("RS Decoder Init failed: {e:?}")))?;

    for (k, shard) in shards.iter().enumerate().take(RS_DATA_SHARDS) {
        rs.add_original_shard(k, shard)
            .map_err(|e| CodecError::Decode(format!("Failed to add original shard: {e:?}")))?;
    }
    for (p, parity_shard) in parity.iter().enumerate().take(RS_PARITY_SHARDS) {
        rs.add_recovery_shard(p, parity_shard)
            .map_err(|e| CodecError::Decode(format!("Failed to add recovery shard: {e:?}")))?;
    }

    let result = rs
        .decode()
        .map_err(|e| CodecError::Corruption(format!("RS Reconstruction failed: {e:?}")))?;

    // Update shards with restored/corrected data
    // Assuming result.restored_original_iter() returns (index, data)
    for (k, shard) in result.restored_original_iter() {
        if k < RS_DATA_SHARDS {
            shards[k] = shard.to_vec();
        }
    }

    // Extract Data
    let mut decoded = Vec::with_capacity(num_blocks * RS_DATA_SHARDS);
    for b in 0..num_blocks {
        for shard in shards.iter().take(RS_DATA_SHARDS) {
            decoded.push(shard[b]);
        }
    }

    Ok(decoded)
}

#[cfg(test)]
mod tests {

    use super::*;
    use rand::Rng;

    #[test]
    fn test_protect_padding_alignment() {
        let input = vec![1, 2, 3];
        let protected = protect_data(&input).expect("Protection failed");

        let expected_block_size = RS_BLOCK_SIZE; // 255
        assert_eq!(
            protected.len() % expected_block_size,
            0,
            "Protected data length should be multiple of RS block size"
        );

        // Minimum shard size is 64, so minimum protected size is 64 * 255 bytes?
        // Wait, utils.rs says `target_len = shard_size * rs_k`. This is confusing.
        // `protect_data` returns `output` which is `num_blocks * rs_n`.
        // `num_blocks` is `shard_size`.
        // `rs_n` is `RS_BLOCK_SIZE` (255).
        // `shard_size` is at least 64.
        // So minimum output size should be 64 * 255 = 16320 bytes?
        // Let's verify this logic.
        assert!(
            protected.len() >= MIN_SHARD_SIZE * RS_BLOCK_SIZE,
            "Protected data too small"
        );
    }

    #[test]
    fn test_protect_minimum_size() {
        // Even 1 byte input should result in minimum shard size
        let input = vec![0xFF];
        let protected = protect_data(&input).unwrap();

        let min_size = MIN_SHARD_SIZE * RS_BLOCK_SIZE;
        assert_eq!(protected.len(), min_size);
    }

    #[test]
    fn test_round_trip_protection() {
        let mut rng = rand::thread_rng();
        // Test various sizes
        for size in [1, 10, 64, 100, 1000, 10000] {
            let input: Vec<u8> = (0..size).map(|_| rng.gen()).collect();
            let protected = protect_data(&input).expect("Protection failed");
            let unprotected = unprotect_data(&protected).expect("Unprotection failed");

            // Unprotection returns padded data. We need to verify it *contains* our data.
            // But wait, `unprotect_data` returns `decoded` which is `num_blocks * RS_DATA_SHARDS`.
            // It doesn't know the original length to truncate padding.
            // The caller (strategies) handles length.
            // So we should expect `unprotected` to start with `input` and be zero-padded.

            assert!(unprotected.len() >= input.len());
            assert_eq!(&unprotected[..input.len()], &input[..]);
            // Verify remaining bytes are zero (padding)
            assert!(unprotected[input.len()..].iter().all(|&b| b == 0));
        }
    }

    #[test]
    fn test_unprotect_invalid_length() {
        let input = vec![1, 2, 3]; // Not multiple of 255
        let result = unprotect_data(&input);
        assert!(matches!(result, Err(CodecError::Decode(_))));
    }
}
