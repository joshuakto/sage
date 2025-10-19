use thiserror::Error;

#[derive(Error, Debug)]
pub enum MaxLfqError {
    #[error("Insufficient peptides for quantification: {0}")]
    InsufficientPeptides(String),

    #[error("Disconnected protein graph")]
    DisconnectedProtein,

    #[error("Numerical optimization failed: {0}")]
    OptimizationError(String),

    #[error("Invalid matrix dimensions")]
    InvalidDimensions,

    #[error("Failed to compute global normalization offsets")]
    NormalizationFailed,

    #[error("Invalid reference sample index: {0}")]
    InvalidReferenceSample(usize),
}
