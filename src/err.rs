use thiserror::Error;

use rmp_serde;

#[derive(Error, Debug)]
pub enum PileupError {
    #[error("Chromosome not found: {0}")]
    ChromosomeNotFound(String),
    #[error("Position not found: {0}")]
    PositionNotFound(u32),
}

#[derive(Error, Debug)]
pub enum UmiVarCalError {
    #[error(transparent)]
    Io(#[from] std::io::Error),

    #[error(transparent)]
    Serde(#[from] rmp_serde::encode::Error),

    #[error(transparent)]
    Utf8Decode(#[from] std::string::FromUtf8Error),

    #[error(transparent)]
    PileupError(#[from] PileupError),
}
