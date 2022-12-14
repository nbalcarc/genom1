use std::{error::Error, fmt::{self, Display}};

/// Used for internal error management
#[derive(Debug)]
pub enum PhyloError {
    FileReadError(String),
    FileTooSmall(String),
    FileOpenError(String),
    KTooBig(u32),
    //ConversionError,
    SearchGenomeError(String),
    SearchNodeError(String),
    GenomeInsertError(String),
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
            //Self::ConversionError => {
            //    write!(f, "ConversionError")
            //},
            Self::SearchGenomeError(s) => {
                write!(f, "SearchGenomeError ({})", s)
            },
            Self::SearchNodeError(s) => {
                write!(f, "SearchNodeError ({})", s)
            },
            Self::GenomeInsertError(s) => {
                write!(f, "GenomeInsertError ({})", s)
            },
        }
    }
}


