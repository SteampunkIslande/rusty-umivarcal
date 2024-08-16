use thiserror::Error;

use rmp_serde;

#[derive(Error, Debug)]
pub enum UmiVarCalError {
    #[error(transparent)]
    Io(#[from] std::io::Error),

    #[error(transparent)]
    Serde(#[from] rmp_serde::encode::Error),
}
