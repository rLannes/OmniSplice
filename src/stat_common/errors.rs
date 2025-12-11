use thiserror::Error;



#[derive(Error, Debug, Clone)]
pub enum LogisticRegressionError {
    /// Perfect or quasi-complete separation detected
    #[error("SeparationError: Perfect or quasi-complete separation detected: {message}")]
    PerfectSeparation { message: String },
    
    /// Probability value outside valid range [0, 1]
    #[error("InvalidProbabilityError: Invalid probability value: {0} (must be in range [0, 1])")]
    InvalidProbability(f64),
    
    /// IRLS failed to converge within maximum iterations
    #[error("FailToConvergeError: IRLS failed to converge after {iterations} iterations (final norm: {final_norm:.6})")]
    ConvergenceFailure { iterations: usize, final_norm: f64 },
    
    /// Matrix is singular and cannot be inverted
    #[error("SingluarMatrixError: Singular matrix encountered: {0}")]
    SingularMatrix(String),
    
    /// Invalid input dimensions
    #[error("DimensionError: Dimension mismatch: expected {expected}, got {got}")]
    DimensionMismatch { expected: String, got: String },
    
    /// Invalid input data
    #[error("InvalidInputError: Invalid input data: {0}")]
    InvalidData(String),
    
    /// Numerical instability detected
    #[error("NumericalInsatbilityError: Numerical instability detected: {0}")]
    NumericalInstability(String),

    #[error("EmptyDataError: No observation (all zeroes)")]
    EmptyData,
    
    #[error("ControlIsNullDataError: No observation for control (all zeroes)")]
    ControlIsNull,

    #[error("TreatmentIsNullDataError: No observation for treatment (all zeroes)")]
    TreatmentIsNull
}



