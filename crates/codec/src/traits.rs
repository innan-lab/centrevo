use crate::error::CodecError;

/// Core trait for sequence codecs.
///
/// This trait defines the "contract" that all compression strategies must follow.
/// Any strategy (like `BitPackedRS` or `UnpackedRS`) must be able to:
/// 1.  `encode`: Take a raw sequence and turn it into protected binary data.
/// 2.  `decode`: Take protected binary data and turn it back into the original sequence.
pub trait Codec {
    fn encode(&self, seq: &[u8]) -> Result<Vec<u8>, CodecError>;
    fn decode(&self, data: &[u8]) -> Result<Vec<u8>, CodecError>;
}
