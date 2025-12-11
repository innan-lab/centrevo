use crate::error::CodecError;
use crate::strategies::BitPackedRS;
use crate::traits::Codec;
use crate::utils::calculate_protected_size;
use rayon::prelude::*;

/// Strategy: Parallel Bit-packed (2 bits/base) + Reed-Solomon.
///
/// # How it works
///
/// This strategy uses the same **2-bit packing** as `BitPackedRS`, but it processes
/// the data in "chunks" (1MB each) simultaneously using multiple CPU cores.
///
/// # Analogy
///
/// Imagine a supermarket checkout:
/// *   `BitPackedRS` is a single cashier scanning items one by one.
/// *   `ParallelBitPackedRS` opens 10 lanes so 10 parts of the order are scanned at once.
///
/// # Why use this?
///
/// *   **Speed:** Much faster for very large files.
/// *   **Overhead:** Slightly slower for tiny files (due to setting up the "lanes").
pub struct ParallelBitPackedRS;

// Size of input chunks (bases) to process in parallel
const CHUNK_SIZE: usize = 1024 * 1024; // 1MB

impl Codec for ParallelBitPackedRS {
    fn encode(&self, seq: &[u8]) -> Result<Vec<u8>, CodecError> {
        if seq.is_empty() {
            return BitPackedRS.encode(seq);
        }

        let chunks: Vec<&[u8]> = seq.chunks(CHUNK_SIZE).collect();

        // Parallel Encoded
        let encoded_portions: Result<Vec<Vec<u8>>, CodecError> = chunks
            .par_iter()
            .map(|chunk| BitPackedRS.encode(chunk))
            .collect();

        let encoded_portions = encoded_portions?;

        // Calculate total size and pre-allocate
        let total_size: usize = encoded_portions.iter().map(|p| p.len()).sum();
        let mut result = Vec::with_capacity(total_size);

        for portion in encoded_portions {
            result.extend_from_slice(&portion);
        }

        Ok(result)
    }

    fn decode(&self, data: &[u8]) -> Result<Vec<u8>, CodecError> {
        if data.is_empty() {
            return BitPackedRS.decode(data);
        }

        // Determine chunk boundaries
        // We know that for all full chunks (size CHUNK_SIZE), the encoded size is constant.
        // Step 1: Calculate packed size for a full chunk
        // BitPacked logic: headers (8 bytes) + packed (ceil(len/4)) + RS padding

        let chunk_packed_len = CHUNK_SIZE.div_ceil(4);
        let chunk_payload_len = 8 + chunk_packed_len; // Header + Data
                                                      // Protect data padding logic:
        let chunk_encoded_len = calculate_protected_size(chunk_payload_len);

        // Slice up the data
        // We have N full chunks and maybe 1 remainder.
        // We can't rely just on divisibility because remainder might accidentally be same size (unlikely but possible).
        // But since we enforced encoding with CHUNK_SIZE, we know the stream structure.

        let mut data_slices = Vec::new();
        let mut cursor = 0;

        while cursor < data.len() {
            let remaining = data.len() - cursor;
            if remaining >= chunk_encoded_len {
                // Is this a full chunk or the last chunk?
                // If there's more data after this block, it MUST be a full chunk.
                // If this is exactly the end, it COULD be a full chunk.
                // BUT wait, what if the last chunk was *also* 1MB?
                // Then we have 2 full chunks.
                // So we can greedily take `chunk_encoded_len` as long as there is more data or exactly that amount left?
                // Yes, because `encode` produces exactly `chunk_encoded_len` for 1MB input.
                // Any input < 1MB will produce < `chunk_encoded_len`?
                // Check `calculate_protected_size` monotonicity.
                // Yes, smaller input -> smaller protected size (mostly, simpler step function).
                // Is it possible for `calculate_protected_size` to return `chunk_encoded_len` for a smaller input?
                // `protect_data` aligns to `MAX(64, input/223)`.
                // If input is slightly smaller, `input/223` drops.
                // `chunk_payload_len` = 8 + 262144 = 262152.
                // `raw_shard` = 262152 / 223 = 1175.5 -> 1176.
                // `shard_size` = 1176. Remainder 1176%64 = 24. padded += 40 -> 1216.
                // `target` = 1216 * 255 = 310080.

                // If input was 1 byte less (262151).
                // `raw` = 1175.56 -> 1176. Same shard size.
                // So yes, a slightly smaller input could result in same encoded size.
                // HOWEVER, `encode` only produces smaller inputs for the *last* chunk.
                // So we can safely act greedily:
                // If `remaining > chunk_encoded_len`, then the current block IS a full chunk.
                // If `remaining == chunk_encoded_len`, it IS a full chunk (and the last one).
                // If `remaining < chunk_encoded_len`, it IS the last chunk (partial).

                // What if the last chunk produced exactly `chunk_encoded_len` but wasn't full?
                // Then it's the last chunk anyway.
                // The only ambiguity is: Could `LastChunk` be bigger than `FullChunk`? No.

                if remaining >= chunk_encoded_len {
                    data_slices.push(&data[cursor..cursor + chunk_encoded_len]);
                    cursor += chunk_encoded_len;
                } else {
                    // Last partial chunk
                    data_slices.push(&data[cursor..]);
                    cursor = data.len();
                }
            } else {
                // Last partial chunk (smaller than full)
                data_slices.push(&data[cursor..]);
                cursor = data.len();
            }
        }

        // Parallel decode
        let decoded_chunks: Result<Vec<Vec<u8>>, CodecError> = data_slices
            .par_iter()
            .map(|slice| BitPackedRS.decode(slice))
            .collect();

        let decoded_chunks = decoded_chunks?;

        // Concatenate
        let total_decoded_len: usize = decoded_chunks.iter().map(|c| c.len()).sum();
        let mut result = Vec::with_capacity(total_decoded_len);
        for chunk in decoded_chunks {
            result.extend_from_slice(&chunk);
        }

        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

    #[test]
    fn test_parallel_round_trip() {
        let codec = ParallelBitPackedRS;
        let mut rng = rand::thread_rng();
        // Test with data larger than CHUNK_SIZE
        let size = CHUNK_SIZE + 500; // 1MB + 500 bases
        let input: Vec<u8> = (0..size).map(|_| rng.gen_range(0..4)).collect();

        let encoded = codec.encode(&input).expect("Encoding failed");
        let decoded = codec.decode(&encoded).expect("Decoding failed");

        assert_eq!(decoded, input);
    }

    #[test]
    fn test_parallel_exact_chunk() {
        let codec = ParallelBitPackedRS;
        let mut rng = rand::thread_rng();
        let size = CHUNK_SIZE * 2; // Exactly 2 chunks
        let input: Vec<u8> = (0..size).map(|_| rng.gen_range(0..4)).collect();

        let encoded = codec.encode(&input).expect("Encoding failed");
        let decoded = codec.decode(&encoded).expect("Decoding failed");

        assert_eq!(decoded, input);
    }
}
