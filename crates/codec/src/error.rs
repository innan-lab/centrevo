use thiserror::Error;

/// Error type for codec operations.
#[derive(Debug, Error)]
pub enum CodecError {
    #[error("Encoding error: {0}")]
    Encode(String),
    #[error("Decoding error: {0}")]
    Decode(String),
    #[error("Data corruption detected: {0}")]
    Corruption(String),
}
