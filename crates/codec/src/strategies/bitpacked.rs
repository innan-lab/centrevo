use crate::error::CodecError;
use crate::traits::Codec;
use crate::utils::{protect_data, unprotect_data};

/// Strategy: Bit-packed (2 bits/base) + Reed-Solomon.
pub struct BitPackedRS;

impl Codec for BitPackedRS {
    fn encode(&self, seq: &[u8]) -> Result<Vec<u8>, CodecError> {
        // 1. Pack bases
        // Fix clippy::manual_div_ceil: (seq.len() + 3) / 4 -> seq.len().div_ceil(4)
        let packed_len = seq.len().div_ceil(4);
        // Header: 8 bytes for original length
        let mut packed_data = Vec::with_capacity(8 + packed_len);

        packed_data.extend_from_slice(&(seq.len() as u64).to_le_bytes());

        let mut byte: u8 = 0;
        for (i, &base) in seq.iter().enumerate() {
            let shift = 2 * (3 - (i % 4));
            byte |= (base & 0x03) << shift;
            if (i + 1) % 4 == 0 {
                packed_data.push(byte);
                byte = 0;
            }
        }
        if seq.len() % 4 != 0 {
            packed_data.push(byte);
        }

        // 2. Protect with RS
        protect_data(&packed_data)
    }

    fn decode(&self, data: &[u8]) -> Result<Vec<u8>, CodecError> {
        // 1. Unprotect
        let packed_data = unprotect_data(data)?;

        // 2. Extract Length
        if packed_data.len() < 8 {
            return Err(CodecError::Decode("Payload too short for header".into()));
        }
        let mut len_bytes = [0u8; 8];
        len_bytes.copy_from_slice(&packed_data[0..8]);
        let original_len = u64::from_le_bytes(len_bytes) as usize;

        // 3. Unpack
        let mut seq = Vec::with_capacity(original_len);
        let start_offset = 8;
        let mut bases_read = 0;

        for &byte in &packed_data[start_offset..] {
            if bases_read >= original_len {
                break;
            }
            for i in 0..4 {
                if bases_read >= original_len {
                    break;
                }
                seq.push((byte >> (2 * (3 - i))) & 0x03);
                bases_read += 1;
            }
        }

        if seq.len() != original_len {
            return Err(CodecError::Decode(format!(
                "Length mismatch: expected {original_len}, got {}",
                seq.len()
            )));
        }

        Ok(seq)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

    #[test]
    fn test_encode_empty() {
        let rs = BitPackedRS;
        let empty: Vec<u8> = vec![];
        let encoded = rs.encode(&empty).expect("Encoding failed");
        let decoded = rs.decode(&encoded).expect("Decoding failed");
        assert_eq!(decoded, empty);
    }

    #[test]
    fn test_encode_single_base() {
        let rs = BitPackedRS;
        let input = vec![2]; // C
        let encoded = rs.encode(&input).expect("Encoding failed");
        let decoded = rs.decode(&encoded).expect("Decoding failed");
        assert_eq!(decoded, input);
    }

    #[test]
    fn test_encode_full_byte() {
        let rs = BitPackedRS;
        let input = vec![0, 1, 2, 3]; // A, C, G, T
        let encoded = rs.encode(&input).expect("Encoding failed");
        let decoded = rs.decode(&encoded).expect("Decoding failed");
        assert_eq!(decoded, input);
    }

    #[test]
    fn test_encode_multi_byte() {
        let rs = BitPackedRS;
        let input = vec![0, 1, 2, 3, 0]; // 5 bases
        let encoded = rs.encode(&input).expect("Encoding failed");
        let decoded = rs.decode(&encoded).expect("Decoding failed");
        assert_eq!(decoded, input);
    }

    #[test]
    fn test_round_trip_random() {
        let rs = BitPackedRS;
        let mut rng = rand::thread_rng();
        for _ in 0..10 {
            let len = rng.gen_range(1..1000);
            let input: Vec<u8> = (0..len).map(|_| rng.gen_range(0..4)).collect();
            let encoded = rs.encode(&input).expect("Encoding failed");
            let decoded = rs.decode(&encoded).expect("Decoding failed");
            assert_eq!(decoded, input);
        }
    }

    #[test]
    fn test_decode_too_short() {
        let rs = BitPackedRS;
        let invalid = vec![0u8; 4]; // Too short for header
        let result = rs.decode(&invalid);
        assert!(result.is_err());
    }

    #[test]
    fn test_rs_corruption_persistence() {
        let rs = BitPackedRS;
        let input = vec![0; 100];
        let mut encoded = rs.encode(&input).expect("Encoding failed");

        // Corrupt data (single byte)
        if let Some(byte) = encoded.get_mut(10) {
            *byte ^= 0xFF;
        }

        // The decoder treats provided original shards as authoritative if not marked missing.
        // So it won't correct the error; it will return the corrupted data.
        let decoded = rs
            .decode(&encoded)
            .expect("Decoding should succeed even with corruption");
        assert_ne!(
            decoded, input,
            "Corruption should persist (no magic correction without checksums)"
        );
    }
}
