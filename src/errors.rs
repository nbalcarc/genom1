use std::{error::Error, fmt::{self, Display}};

/// Used for internal error management
#[derive(Debug)]
enum PhyloError {
    FileNotFound(String),
}
impl Error for PhyloError {}
impl Display for PhyloError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self {
            Self::FileNotFound(s) => {
                write!(f, "{}", s)
            }
        }
    }
}