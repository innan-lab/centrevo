use thiserror::Error;

/// Error type for codec operations.
///
/// This lists everything that can go wrong when compressing or decompressing.
#[derive(Debug, Error)]
pub enum CodecError {
    /// **Encoding Failed**: Something went wrong while trying to compress the data.
    #[error("Encoding error: {0}")]
    Encode(String),
    /// **Decoding Failed**: Something went wrong while trying to retrieve the original data.
    #[error("Decoding error: {0}")]
    Decode(String),
    /// **Data Corruption**: The data was damaged, and even the Reed-Solomon protection
    /// couldn't save it (or checksums failed).
    #[error("Data corruption detected: {0}")]
    Corruption(String),
}
