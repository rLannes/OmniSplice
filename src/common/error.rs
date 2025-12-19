//use crate::errors::{Error, Result};
use CigarParser::cigar::CigarError;

#[derive(Debug, thiserror::Error)]
pub enum OmniError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Parse error: {0}")]
    WrongFile(String),

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

    #[error("Error converting: {0}")]
    fromutf8(#[from] std::string::FromUtf8Error),
    //#[error("rust bio error: {0}")]
    //RustBio(#[from] rust_htslib::errors::Error),
}



#[cfg(test)]
mod error_tests {
    use crate::common::error::OmniError;
    
    #[test]
    fn test_io_error_conversion() {
        let io_err = std::io::Error::new(std::io::ErrorKind::NotFound, "file not found");
        let omni_err: OmniError = io_err.into();
        assert!(matches!(omni_err, OmniError::Io(_)));
    }

    #[test]
    fn test_parse_error_display() {
        let err = OmniError::Parse("test parse error".to_string());
        let display = format!("{}", err);
        assert!(display.contains("Parse error"));
        assert!(display.contains("test parse error"));
    }

    #[test]
    fn test_expect_error() {
        let err = OmniError::Expect("expected value not found".to_string());
        assert!(format!("{}", err).contains("Value absent when expected"));
    }

    #[test]
    fn test_gtf_parse_error() {
        let err = OmniError::GTFParse("invalid GTF format".to_string());
        assert!(format!("{}", err).contains("GTF Parsing error"));
    }

    #[test]
    fn test_output_exists_error() {
        let io_err = std::io::Error::new(std::io::ErrorKind::AlreadyExists, "exists");
        let err = OmniError::OutputExists("output.txt".to_string(), io_err);
        assert!(format!("{}", err).contains("output.txt"));
        assert!(format!("{}", err).contains("should not exist"));
    }
}