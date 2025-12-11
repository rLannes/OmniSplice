use super::super::common::error::OmniError;
use super::errors::LogisticRegressionError;
use super::common::{Genotype, Tester, CountsStats, SplicingCategory, JunctionStats, TestResults, TestStatus};
//use fishers_exact::fishers_exact;
use statrs::stats_tests::fishers_exact;
use statrs::stats_tests::Alternative;
use statrs::stats_tests::fishers_exact_with_odds_ratio;

use core::f64;
use std::collections::{HashMap, HashSet};
use std::fmt::format;
use std::iter::Sum;
use std::mem::min_align_of;
use std::path::Path;
use std::{fs, result};
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::prelude::*;
use std::io::{BufReader, BufWriter};
use std::convert::TryInto;

use nalgebra::{DMatrix, DVector};
use statrs::distribution::{ChiSquared, ContinuousCDF};
use flexi_logger::{FileSpec, Logger, WriteMode};
use log::{debug, error, info, trace, warn};

pub struct GLM{

    success: Vec<u32>,
    failures: Vec<u32>,
    groups: Vec<Genotype>,
    identifier: String

}

impl Tester for GLM {

    fn success(&self) -> &Vec<u32>{&self.success}
    fn success_mut(&mut self) -> &mut Vec<u32>{&mut self.success}
    fn failures(&self) -> &Vec<u32>{&self.failures}
    fn failures_mut(&mut self) -> &mut Vec<u32>{&mut self.failures}
    fn groups(&self) -> &Vec<Genotype>{&self.groups}
    fn groups_mut(&mut self) -> &mut Vec<Genotype>{&mut self.groups}

    fn test(&self, donotrun: bool) -> TestResults {

        let mut test_res = TestResults::get_empty();
        let (ctrl_suc, ctrl_fail, treat_suc, treat_fail) = self.get_proportion();
        test_res.control_failure = ctrl_fail;
        test_res.control_sucess = ctrl_suc;
        test_res.treatment_failure = treat_fail;
        test_res.treatment_sucess = treat_suc;

        test_res.string_count = self.get_proportion_string();

        //info!("Counts prop: {} {} {} {}", test_res.control_failure, test_res.control_sucess, test_res.treatment_failure, test_res.treatment_sucess);
        

        let treat_trial = treat_suc + treat_fail;
        let ctrl_trial = ctrl_fail + ctrl_suc;
        if (treat_trial == ctrl_trial) && (ctrl_trial == 0){
            test_res.status = Some(TestStatus::EmptyData);
            return test_res;
        }
        else if           ctrl_trial == 0{
            test_res.status = Some(TestStatus::ControlIsNull);
            return test_res;
        }
        else if treat_trial == 0{
            test_res.status = Some(TestStatus::TreatmentIsNull);
            return test_res;
        }
        else{
            test_res.control_prop = Some(ctrl_suc as f32 / ctrl_trial as f32);
            test_res.treatment_prop = Some(treat_suc as f32 / treat_trial as f32);
            
            if donotrun == true{
                test_res.status = Some(TestStatus::Ambigious);
                return test_res;
            }

            match self.actual_test() {
                Ok((status, pval, odr, ord_lower, odr_high)) => {
                    test_res.status = Some(status);
                    test_res.p_value = Some(pval);
                    test_res.odd_ratio = Some(odr);
                    test_res.or_ci_lower = Some(ord_lower);
                    test_res.or_ci_upper = Some(odr_high);

                },
                Err(e) => {
                    warn!("error: {} {}", self.identifier, e);
                    test_res.status = Some(e.into())
                }
            }

        }

        test_res

    }
    
}
fn vec_to_array<T, const N: usize>(v: Vec<T>) -> [T; N] {
    v.try_into()
        .unwrap_or_else(|v: Vec<T>| panic!("Expected a Vec of length {} but it was {}", N, v.len()))
}
impl GLM{

   

    pub fn actual_test(&self) -> Result<(TestStatus, f64, f64, f64, f64), LogisticRegressionError> {

        let mut current_testStatus = TestStatus::Ok;

        let mut successes: Vec<f64> = self.success().into_iter().map(|x|  *x as f64).collect();
        let mut failures: Vec<f64> = self.failures().into_iter().map(|x|  *x as f64).collect();

        // if low count we fallback to Fischer Test instead:
        let contegency: Vec<u64> = self.to_contengency();
        if (*contegency.iter().min().unwrap_or(&0) < 5) || 
            (contegency.iter().sum::<u64>() < 30 ) || 
            (contegency[0] + contegency[1] < 10) || (contegency[2] + contegency[3] < 10) || 
            (contegency[0] + contegency[2] < 10) || (contegency[1] + contegency[3] < 10)
            {
                warn!("low count, fall back to Fischer {:?}", contegency);
                //println!("low count, fall back to Fischer {:?}", contegency);
                let p = fishers_exact(&vec_to_array(contegency.clone()), Alternative::TwoSided).unwrap();
                let mut odd_ratio = 0. ;
                if  (contegency[1] != 0) && (contegency[3] != 0){
                    odd_ratio = (contegency[0] as f64 / contegency[1] as f64 ) /(contegency[2] as f64 / contegency[3] as f64 );
                }
            return Ok((TestStatus::FISCHERFallBack, p, odd_ratio,  0. as f64 , 0. as f64)); 
        }


        let n_obs = self.groups().len();
        let y = DVector::from_vec(successes.clone());
        let n_trials: Vec<f64> = successes.iter()
        .zip(failures.iter())
        .map(|(s, f)| s + f)
        .collect();
        let n = DVector::from_vec(n_trials.clone());
    
        // Fit null model: intercept only
    
        let x_null = DMatrix::from_element(n_obs, 1, 1.0);
        let (beta_null, status) = GLM::irls_binomial(&x_null, &y, &n, 150, 1e-8)?;
    
        if status > current_testStatus{
            current_testStatus = status;
        }
        // Compute null model predictions and log-likelihood
        let eta_null = &x_null * &beta_null;
        let mu_null: Vec<f64> = eta_null.iter().map(|&e| GLM::inv_logit(e)).collect();

        let ll_null = GLM::binomial_log_likelihood(&successes, &n_trials, &mu_null);
    

    
    // Fit full model: intercept + treatment indicator
    let treatment: Vec<f64> = self.groups().iter()
        .map(|g| match g {
            Genotype::TREATMENT  => { 1.0 },
            Genotype::CONTROL => {0.0}
        })
        .collect();
    
    // Build design matrix [1, treatment]
    let mut x_full_data = vec![1.0; n_obs];
    x_full_data.extend(treatment);
    let x_full = DMatrix::from_vec(n_obs, 2, x_full_data);
    
    let (beta_full, status) = GLM::irls_binomial(&x_full, &y, &n, 150, 1e-8)?;
    if status > current_testStatus{
        current_testStatus = status
    }
    
    // Compute full model predictions and log-likelihood
    let eta_full = &x_full * &beta_full;
    let mu_full: Vec<f64> = eta_full.iter().map(|&e| GLM::inv_logit(e)).collect();
    let ll_full = GLM::binomial_log_likelihood(&successes, &n_trials, &mu_full);
    // Perform likelihood ratio test
    let lr_stat = 2.0 * (ll_full - ll_null);
    let df = 1.0;
    
    // Calculate p-value using chi-squared distribution
    let chi_sq = ChiSquared::new(df).expect("Failed to create chi-squared distribution");
    let p_value = chi_sq.sf(lr_stat);// 1.0 - chi_sq.cdf(lr_stat);
    
    let se = GLM::standard_errors(&x_full, &beta_full, &n);
    let (or, or_lower, or_upper) = GLM::odds_ratio_with_ci(beta_full[1], se[1]);

    // if or is really big > 100 most likely model collapsed fall back to Fischer
    if or > 100. || or < 0.01 {
                warn!("odd ratio extrem value, fall back to Fischer {:?}", contegency);
                //println!("low count, fall back to Fischer {:?}", contegency);
                let p = fishers_exact(&vec_to_array(contegency.clone()), Alternative::TwoSided).unwrap();
                let mut odd_ratio = 0. ;
                if  (contegency[1] != 0) && (contegency[3] != 0){
                    odd_ratio = (contegency[0] as f64 / contegency[1] as f64 ) /(contegency[2] as f64 / contegency[3] as f64 );
                }
                else{
                    odd_ratio = f64::NAN;
                }
            return Ok((TestStatus::FISCHERFallBack, p
                , odd_ratio,  0. as f64 , 0. as f64)); 
    }
    //let elapsed_time = now.elapsed();

    Ok((TestStatus::Ok, p_value, or, or_lower, or_upper))

    }

    pub fn new(treatment: &Vec<CountsStats>,
            control: &Vec<CountsStats>,
            successes_cat: &Vec<SplicingCategory>,
            failures_cat: &Vec<SplicingCategory>,
            identifier: String) -> Self {
        let mut glm = GLM{
            success: Vec::new(),
            failures: Vec::new(),
            groups: Vec::new(),
            identifier: identifier
        };

    glm.format_data(treatment, control, successes_cat, failures_cat);
    glm

    }





    /// Fits a binomial logistic regression model using Iteratively Reweighted Least Squares (IRLS).
    ///
    /// This implements the Fisher scoring algorithm for maximum likelihood estimation
    /// of generalized linear models with a binomial response and logit link.
    ///
    /// # Arguments
    /// * `x` - Design matrix (n_obs × n_vars) containing predictor variables
    /// * `y` - Response vector of success counts
    /// * `n` - Vector of trial counts for each observation
    /// * `max_iter` - Maximum number of IRLS iterations
    /// * `tol` - Convergence tolerance for coefficient estimates
    ///
    /// # Returns
    /// Vector of estimated regression coefficients (β)
    ///
    /// # Algorithm
    /// 1. Initialize β = 0
    /// 2. Repeat until convergence:
    ///    - Compute linear predictor: η = Xβ
    ///    - Compute fitted probabilities: μ = inv_logit(η)
    ///    - Compute weights: W = diag(n * μ * (1-μ))
    ///    - Compute working response: z = η + (y - nμ) / (nμ(1-μ))
    ///    - Update: β_new = (X'WX)^(-1) X'Wz
    ///
    /// # Panics
    /// Panics if the weighted least squares system cannot be solved
    fn irls_binomial(x: &DMatrix<f64>, y: &DVector<f64>, n: &DVector<f64>, max_iter: usize, tol: f64) -> Result<(DVector<f64>, TestStatus), LogisticRegressionError>  {
        
        let mut large_sep_warn = false;
        let n_obs = x.nrows();
        let n_vars = x.ncols();
        let mut beta = DVector::zeros(n_vars);

        let mut current_state = TestStatus::Ok;

        GLM::validate_irls_inputs(x, y, n)?;
        
        for iter in 0..max_iter {
            // Compute linear predictor
            let eta = x * &beta;
            
            for (i, &e) in eta.iter().enumerate() {
                if !e.is_finite() {
                    return Err(LogisticRegressionError::NumericalInstability(
                        format!("Non-finite linear predictor at iteration {} observation {}: eta={}",
                            iter + 1, i, e)
                    ));
                }
                
                // Warn about extreme values that might cause problems
                if (!large_sep_warn) && (e.abs() > 20.0) {
                    large_sep_warn = true;
                    if current_state < TestStatus::QuasiPerfectSeparation{
                        current_state = TestStatus::QuasiPerfectSeparation
                    }
                    //eprintln!("Warning: Large linear predictor at iteration {} observation {}: eta={:.2}. Possible separation.", 
                    //    iter + 1, i, e);
                }
            }

            // Compute fitted probabilities
            let mu: Vec<f64> = eta.iter().map(|&e| GLM::inv_logit(e)).collect();
            
            // Compute weights: W = n * μ * (1 - μ)
            let w: Vec<f64> = mu.iter()
                .zip(n.iter())
                .map(|(&m, &ni)| ni * m * (1.0 - m).max(1e-10))
                .collect();
            
            let w_mat = DMatrix::from_diagonal(&DVector::from_vec(w.clone()));
            
            // Compute working response: z = η + (y - nμ) / (nμ(1-μ))
            let z: Vec<f64> = eta.iter()
                .zip(y.iter())
                .zip(n.iter())
                .zip(mu.iter())
                .map(|(((e, yi), ni), mi)| {
                    e + (yi - ni * mi) / (ni * mi * (1.0 - mi)).max(1e-10)
                })
                .collect();


            
            // Check working response for numerical issues
            for (i, zi) in z.iter().enumerate() {
                if !zi.is_finite() {
                    return Err(LogisticRegressionError::NumericalInstability(
                        format!("Non-finite working response at iteration {} observation {}: z={}",
                            iter + 1, i, zi)
                    ));
                }
            }
            let z_vec = DVector::from_vec(z);
            
            // Solve weighted least squares: β_new = (X'WX)^(-1) X'Wz
            let xtw = x.transpose() * &w_mat;
            let xtwx = &xtw * x;
            let xtwz = &xtw * z_vec;
            
            let beta_new = match xtwx.lu().solve(&xtwz)  {
                Some(solution) => solution,
                None => {
                    // Try to provide more informative error message
                    return Err(LogisticRegressionError::SingularMatrix(
                        format!(
                            "Cannot solve X'WX at iteration {}. Matrix is singular or near-singular. \
                            This often indicates: (1) perfect multicollinearity among predictors, \
                            (2) a predictor with zero variance, or (3) perfect separation. \
                            Check for: redundant variables, constant columns, or perfect prediction.",
                            iter + 1
                        )
                    ));
                }
            };

            for (i, &b) in beta_new.iter().enumerate() {
                if !b.is_finite() {
                    return Err(LogisticRegressionError::NumericalInstability(
                        format!("Non-finite coefficient at iteration {} parameter {}: β={}",
                            iter + 1, i, b)
                    ));
                }
            }
            
            let delta = &beta_new - &beta;
            let delta_norm = delta.norm();
            
            if !delta_norm.is_finite() {
                return Err(LogisticRegressionError::NumericalInstability(
                    format!("Non-finite convergence criterion at iteration {}: ||Δβ||={}",
                        iter + 1, delta_norm)
                ));
            }


            // Check convergence
            if (&beta_new - &beta).norm() < tol {
                return Ok((beta_new, current_state));
            }
            beta = beta_new;
        }
        /// change return to include if it converged
        if TestStatus::ConvergenceFailed > current_state{
            current_state = TestStatus::ConvergenceFailed
        }
        Ok((beta, current_state))
        // return a value if not converged
    }


    /// Computes the logit (log-odds) transformation of a probability.
    ///
    /// # Arguments
    /// * `p` - A probability value between 0 and 1
    ///
    /// # Returns
    /// The log-odds: ln(p / (1-p))
    ///
    /// # Note
    /// Input values should be strictly between 0 and 1 to avoid infinity/NaN
    fn logit(p: f64) -> Result<f64, LogisticRegressionError> {
        if p <= 0. || p >= 1. {return Err(LogisticRegressionError::InvalidProbability(p))};
        Ok((p / (1.0 - p)).ln())

    }

    /// Computes the inverse logit (logistic) function.
    ///
    /// Transforms a real number to a probability in [0, 1].
    ///
    /// # Arguments
    /// * `x` - Any real number
    ///
    /// # Returns
    /// A probability: 1 / (1 + exp(-x))
    fn inv_logit(x: f64) -> f64 {
        1.0 / (1.0 + (-x).exp())
    }

    /// Calculates the binomial log-likelihood for observed data.
    ///
    /// # Arguments
    /// * `y` - Vector of success counts
    /// * `n` - Vector of trial counts (total observations per group)
    /// * `mu` - Vector of predicted probabilities
    ///
    /// # Returns
    /// The sum of log-likelihoods: Σ[y*ln(μ) + (n-y)*ln(1-μ)]
    ///
    /// # Note
    /// Probabilities are clamped to [1e-10, 1-1e-10] to prevent log(0)
    fn binomial_log_likelihood(y: &[f64], n: &[f64], mu: &[f64]) -> f64 {

        y.iter()
            .zip(n.iter())
            .zip(mu.iter())
            .map(|((yi, ni), mui)| {
                // Clamp probability to avoid numerical issues
                let mui = mui.max(1e-10).min(1.0 - 1e-10);
                yi * mui.ln() + (ni - yi) * (1.0 - mui).ln()
            })
            .sum()
    }

    /// Computes the odds ratio and its 95% confidence interval.
    ///
    /// # Arguments
    /// * `beta` - The log odds ratio (coefficient from logistic regression)
    /// * `std_error` - Standard error of the coefficient
    ///
    /// # Returns
    /// Tuple of (odds_ratio, lower_95_ci, upper_95_ci)
    ///
    /// # Note
    /// The odds ratio is exp(β), and the CI is computed as exp(β ± 1.96·SE)
    fn odds_ratio_with_ci(beta: f64, std_error: f64) -> (f64, f64, f64) {
        let or = beta.exp();
        let lower = (beta - 1.96 * std_error).exp();
        let upper = (beta + 1.96 * std_error).exp();
        (or, lower, upper)
    }


    /// Computes standard errors for regression coefficients.
    ///
    /// Uses the inverse of the Fisher information matrix: Var(β) = (X'WX)^(-1)
    ///
    /// # Arguments
    /// * `x` - Design matrix
    /// * `beta` - Fitted coefficients
    /// * `n` - Vector of trial counts
    ///
    /// # Returns
    /// Vector of standard errors for each coefficient
    fn standard_errors(x: &DMatrix<f64>, beta: &DVector<f64>, n: &DVector<f64>) -> DVector<f64> {
        let eta = x * beta;
        let mu: Vec<f64> = eta.iter().map(|&e| GLM::inv_logit(e)).collect();
        
        // Compute weights
        let w: Vec<f64> = mu.iter()
            .zip(n.iter())
            .map(|(&m, &ni)| ni * m * (1.0 - m).max(1e-10))
            .collect();
        
        let w_mat = DMatrix::from_diagonal(&DVector::from_vec(w));
        
        // Variance-covariance matrix: (X'WX)^(-1)
        let xtw = x.transpose() * &w_mat;
        let xtwx = &xtw * x;
        let var_covar = xtwx.try_inverse().expect("Failed to invert matrix");
        
        // Standard errors are square roots of diagonal elements
        DVector::from_iterator(
            var_covar.nrows(),
            var_covar.diagonal().iter().map(|&v| v.sqrt())
        )
    }



    /// Validate the input to shield against the most common error
    fn validate_irls_inputs(x: &DMatrix<f64>, y: &DVector<f64>, n: &DVector<f64>) -> Result<(), LogisticRegressionError>{


        let n_obs = x.nrows();
        let n_vars = x.ncols();
        
        // Check dimension compatibility
        if y.len() != n_obs {
            return Err(LogisticRegressionError::DimensionMismatch {
                expected: format!("y.len() = {}", n_obs),
                got: format!("y.len() = {}", y.len()),
            });
        }
        
        if n.len() != n_obs {
            return Err(LogisticRegressionError::DimensionMismatch {
                expected: format!("n.len() = {}", n_obs),
                got: format!("n.len() = {}", n.len()),
            });
        }
        
        // Need more observations than parameters
        if n_obs <= n_vars {
            return Err(LogisticRegressionError::InvalidData(
                format!("Insufficient observations: {} obs for {} parameters. Need at least {} observations.",
                    n_obs, n_vars, n_vars + 1)
            ));
        }

            // Validate data values
        for i in 0..n_obs {
            // Check for non-finite values
            if !y[i].is_finite() || !n[i].is_finite() {
                return Err(LogisticRegressionError::InvalidData(
                    format!("Non-finite value at observation {}: y={}, n={}", i, y[i], n[i])
                ));
            }
            
            // Check for negative values
            if y[i] < 0.0 {
                return Err(LogisticRegressionError::InvalidData(
                    format!("Negative success count at observation {}: y={}", i, y[i])
                ));
            }
            
            if n[i] <= 0.0 {
                return Err(LogisticRegressionError::InvalidData(
                    format!("Non-positive trial count at observation {}: n={}", i, n[i])
                ));
            }
            
            // Check that successes don't exceed trials
            if y[i] > n[i] {
                return Err(LogisticRegressionError::InvalidData(
                    format!("Success count exceeds trials at observation {}: y={} > n={}", i, y[i], n[i])
                ));
            }
        }

        // Check design matrix for non-finite values
        for i in 0..x.nrows() {
            for j in 0..x.ncols() {
                if !x[(i, j)].is_finite() {
                    return Err(LogisticRegressionError::InvalidData(
                        format!("Non-finite value in design matrix at ({}, {}): {}", i, j, x[(i, j)])
                    ));
                }
            }
        }
        Ok(())
    }

}




#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn glm_test() {
        let mut glm= GLM{
            groups: vec![Genotype::CONTROL, Genotype::CONTROL, Genotype::CONTROL,
                     Genotype::TREATMENT, Genotype::TREATMENT, Genotype::TREATMENT], 
            success: vec![1,0,0,22,1,2],
            failures: vec![4,1,2,0,0,0],
            identifier: "Non".to_string()
        };
        // should be non significant but I found it at 10-7
        let x = glm.test(false);
        println!("x: {:?}", x);
    }

}