use std::{error::Error, fmt::{self, Display}};

/// Used for internal error management
#[derive(Debug)]
pub enum PhyloError {
    FileReadError(String),
    FileTooSmall(String),
    FileOpenError(String),
    KTooBig(u32),
    ConversionError,
    SearchGenomeError,
    SearchNodeError,
    GenomeInsertError,
}
impl Error for PhyloError {}
impl Display for PhyloError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self {
            Self::FileReadError(s) => {
                write!(f, "{}", s)
            },
            Self::FileTooSmall(s) => {
                write!(f, "{}", s)
            },
            Self::FileOpenError(s) => {
                write!(f, "{}", s)
            },
            Self::KTooBig(s) => {
                write!(f, "{}", s)
            },
            Self::ConversionError => {
                write!(f, "Conversion error!")
            },
            Self::SearchGenomeError => {
                write!(f, "Error when searching for a genome!")
            },
            Self::SearchNodeError => {
                write!(f, "Error when searching for a node!")
            },
            Self::GenomeInsertError => {
                write!(f, "Error when inserting a genome!")
            },
        }
    }
}


