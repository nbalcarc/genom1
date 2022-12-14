use std::{error::Error, fmt::{self, Display}};

/// Used for internal error management
#[derive(Debug)]
pub enum PhyloError {
    FileReadError(String),
    FileTooSmall(String),
    FileOpenError(String),
    FileWriteError,
    KTooBig(u32),
    //ConversionError,
    SearchGenomeError(String),
    SearchNodeError(String),
    GenomeInsertError(String),
    FileDeleteError,
    PathError(String),
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
            Self::FileWriteError => {
                write!(f, "Error when trying to write to the output file")
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
            Self::FileDeleteError => {
                write!(f, "Error when attempting to delete phylo_tree.txt")
            },
            Self::PathError(s) => {
                write!(f, "PathError ({})", s)
            }
        }
    }
}


