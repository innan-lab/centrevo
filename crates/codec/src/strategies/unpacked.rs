use crate::error::CodecError;
use crate::traits::Codec;
use crate::utils::{protect_data, unprotect_data};

/// Strategy: Unpacked (1 byte/base) + Reed-Solomon.
pub struct UnpackedRS;

impl Codec for UnpackedRS {
    fn encode(&self, seq: &[u8]) -> Result<Vec<u8>, CodecError> {
        // 1. Prepare data (Add header for safety/consistency)
        let mut data_with_header = Vec::with_capacity(8 + seq.len());
        data_with_header.extend_from_slice(&(seq.len() as u64).to_le_bytes());
        // For Unpacked, we just verify inputs are valid bases?
        // Or just store raw bytes. Let's store raw bytes to be generic.
        // But for consistency with seq being &[u8] indices, we just copy.
        data_with_header.extend_from_slice(seq);

        // 2. Protect
        protect_data(&data_with_header)
    }

    fn decode(&self, data: &[u8]) -> Result<Vec<u8>, CodecError> {
        // 1. Unprotect
        let decoded_data = unprotect_data(data)?;

        // 2. Extract Length
        if decoded_data.len() < 8 {
            return Err(CodecError::Decode("Payload too short for header".into()));
        }
        let mut len_bytes = [0u8; 8];
        len_bytes.copy_from_slice(&decoded_data[0..8]);
        let original_len = u64::from_le_bytes(len_bytes) as usize;

        // 3. Extract Sequence
        if decoded_data.len() < 8 + original_len {
            return Err(CodecError::Decode(
                "Payload shorter than expected sequence".into(),
            ));
        }
        let seq = decoded_data[8..8 + original_len].to_vec();
        Ok(seq)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

    #[test]
    fn test_encode_empty() {
        let rs = UnpackedRS;
        let empty: Vec<u8> = vec![];
        let encoded = rs.encode(&empty).expect("Encoding failed");
        let decoded = rs.decode(&encoded).expect("Decoding failed");
        assert_eq!(decoded, empty);
    }

    #[test]
    fn test_round_trip_random() {
        let rs = UnpackedRS;
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
        let rs = UnpackedRS;
        let invalid = vec![0u8; 4]; // Too short for header
        let result = rs.decode(&invalid);
        assert!(result.is_err());
    }

    #[test]
    fn test_rs_corruption_persistence() {
        let rs = UnpackedRS;
        let input = vec![0; 100];
        let mut encoded = rs.encode(&input).expect("Encoding failed");

        // Corrupt data (single byte)
        if let Some(byte) = encoded.get_mut(50) {
            *byte ^= 0xFF;
        }

        // Expect persistence
        let decoded = rs
            .decode(&encoded)
            .expect("Decoding should succeed even with corruption");
        assert_ne!(
            decoded, input,
            "Corruption should persist (no magic correction without checksums)"
        );
    }
}
