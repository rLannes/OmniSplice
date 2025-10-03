//use crate::errors::{Error, Result};
use CigarParser::cigar::CigarError;

#[derive(Debug, thiserror::Error)]
pub enum OmniError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Parse error: {0}")]
    Parse(String),


    #[error("Value absent when expected: {0}")]
    Expect(String),

    #[error("output file {0} should not exist: {1}")]
    OutputExists(String, #[source] std::io::Error),

    #[error("GTF Parsing error: {0}")]
    GTFParse(String),

    #[error("Invalid CIGAR: {0}")]
    Cigar(#[from] CigarError),

    #[error("rust HTS-LIB error: {0}")]
    HtsLib(#[from] rust_htslib::errors::Error),

    #[error("ParseIntError error: {0}")]
    ParseIntError(#[from] std::num::ParseIntError),

    //#[error("rust bio error: {0}")]
    //RustBio(#[from] rust_htslib::errors::Error),

}

