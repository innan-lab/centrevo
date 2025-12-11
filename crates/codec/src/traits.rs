use crate::error::CodecError;

/// Core trait for sequence codecs.
pub trait Codec {
    fn encode(&self, seq: &[u8]) -> Result<Vec<u8>, CodecError>;
    fn decode(&self, data: &[u8]) -> Result<Vec<u8>, CodecError>;
}
