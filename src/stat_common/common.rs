
use std::hash::{Hash, Hasher};
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::{fs, result};
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::io::prelude::*;
use std::convert::From;
use std::fmt;
use std::convert::TryFrom;

use flexi_logger::{FileSpec, Logger, WriteMode};
use log::{debug, error, info, trace, warn};


use crate::common::utils::ReadAssign;

use super::super::common::error::OmniError;
use super::errors::LogisticRegressionError;
 #[derive(Debug)]
pub struct CountsStats{
    spliced: u32,
    unspliced: u32,
    clipped: u32,
    exon_other: u32,
    skipped: u32,
    skipped_unrelated: u32,
    wrong_strand: u32,
    e_isoform: u32,
}

impl CountsStats{

    pub fn new(string: &[&str]) -> Self {
        let c = string
            .iter()
            .map(|x| x.parse::<u32>().unwrap())
            .collect::<Vec<u32>>();
        CountsStats {
            spliced: c[0],
            unspliced: c[1],
            clipped: c[2],
            exon_other: c[3],
            skipped: c[4],
            skipped_unrelated: c[5],
            wrong_strand: c[6],
            e_isoform: c[7],
        }
    }

    pub fn extract_(&self, category: &Vec<SplicingCategory>) -> u32{
        let mut count = 0;

        for cat in category{
            match cat{
                SplicingCategory::Spliced => {count += self.spliced}, 
                SplicingCategory::Unspliced => {count += self.unspliced},
                SplicingCategory::Clipped => {count += self.clipped},
                SplicingCategory::ExonOther => {count += self.exon_other},
                SplicingCategory::Skipped => {count += self.skipped}, 
                SplicingCategory::SkippedUnrelated => {count += self.skipped_unrelated}, 
                SplicingCategory::WrongStrand => {count += self.wrong_strand},
                SplicingCategory::EIsoform => {count += self.e_isoform}
            }
        }
        count
    }
}

#[derive(Debug)]
pub enum SplicingCategory{
    Spliced, 
    Unspliced,
    Clipped,
    ExonOther,
    Skipped, 
    SkippedUnrelated, 
    WrongStrand,
    EIsoform
}


impl TryFrom<&str> for SplicingCategory {
    type Error = &'static str;
    fn try_from(item: &str) -> Result<Self, Self::Error> {
        match item {
            "Spliced" | "SPLICED" => Ok(SplicingCategory::Spliced),
            "Unspliced" | "UNSPLICED" => Ok(SplicingCategory::Unspliced),
            "Clipped" | "CLIPPED" => Ok(SplicingCategory::Clipped),
            "Exon_other" | "EXONOTHER" => Ok(SplicingCategory::ExonOther),
            "Skipped" | "SKIPPED" => Ok(SplicingCategory::Skipped),
            "SkippedUnrelated"| "SKIPPEDUNRELATED" => Ok(SplicingCategory::SkippedUnrelated),
            "Wrong_strand" | "WRONGSTRAND" => Ok(SplicingCategory::WrongStrand),
            "E_isoform" | "EISOFORM" => Ok(SplicingCategory::EIsoform),
            _ => Err("input does not match one of the accepted values: Spliced Unspliced Clipped Exon_other Skipped SkippedUnrelated Wrong_strand  E_isoform")
        }
    }
}

impl fmt::Display for SplicingCategory {
    // This trait requires `fmt` with this exact signature.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {

        match self {
            SplicingCategory::Spliced => { write!(f, "SPLICED") }, 
            SplicingCategory::Unspliced => { write!(f, "UNSPLICED") },
            SplicingCategory::Clipped => { write!(f, "CLIPPED") },
            SplicingCategory::ExonOther => { write!(f, "EXONOTHER") },
            SplicingCategory::Skipped => { write!(f, "SKIPPED") }, 
            SplicingCategory::SkippedUnrelated => { write!(f, "SKIPPEDUNRELATED") }, 
            SplicingCategory::WrongStrand => { write!(f, "WRONGSTRAND") },
            SplicingCategory::EIsoform => { write!(f, "EISOFORM") }
        }
    }
}

#[derive(Debug)]
pub struct JunctionStats{
    pub contig: String, 
    pub start: String, 
    pub end: String, 
    pub strand: String,
    pub ambiguous: bool,
    pub control_count: Vec<CountsStats>,
    pub treat_count: Vec<CountsStats>,
    pub gene_tr: HashSet<String>,
    pub sample_done: HashSet<String>
}

impl JunctionStats{
    pub fn get_pos_string(&self) -> Vec<String>{
        vec![self.contig.clone(), self.strand.clone(), self.start.clone(), self.end.clone()]
    }

   
}

impl Hash for JunctionStats{
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.contig.hash(state); 
        self.strand.hash(state);
        self.start.hash(state); 
        self.end.hash(state);  
    }
}


fn header_to_map(header:  &str) -> Result<HashMap<String, usize>, OmniError>{
    Ok(header
        .trim()
        .split_whitespace()
        .enumerate()
        .map(|(i, v)| (v.to_string(), i))
        .collect())
    }

#[derive(Debug, PartialEq, Eq)]
pub enum Genotype{
    CONTROL,
    TREATMENT
}

pub fn parse_js_file(file_path: &str, result: &mut HashMap<String, JunctionStats>, genotype: Genotype) -> Result<(), OmniError>{



        let path = Path::new(file_path);
        let file_stem = format!("{:?}", path.file_stem().unwrap());
        info!("opening file");
        let f = File::open(file_path)?;
        let mut reader = BufReader::new(f);
        let mut header_line = String::new();

        info!("reading header");
        let len = reader.read_line(&mut header_line)?;
        let header = header_to_map(&header_line)?;

        
        let mut line: String = "".to_string();
        let mut spt: Vec<&str> = Vec::new();
        let mut trimed : &str = "";
        let mut gene: String = "".to_string();
        let mut contig: String = "".to_string();
        let mut strand : String = "".to_string();
        let mut start: String = "".to_string();
        let mut end: String = "".to_string();
        let mut key: String = "".to_string();

        for iterline in reader.lines(){
            line = iterline?;
            trimed = line.trim();

            if trimed.is_empty(){continue}
            spt = trimed.split("\t").collect::<Vec<&str>>();

            let ambi = if spt[*header.get("Ambiguous").unwrap()] == "true" {true} else {false}; 
            //if (ambiguous == false) && (ambi == true){
            //    continue
            //}
            gene = format!("{}_{}_{}", spt[*header.get("Gene").unwrap()], spt[*header.get("Transcript").unwrap()], spt[*header.get("Intron").unwrap()]);
            contig = spt[*header.get("Contig").unwrap()].to_string();
            strand = spt[*header.get("Strand").unwrap()].to_string();
            start = if strand == "+" {spt[*header.get("Donnor").unwrap()].to_string()} else {spt[*header.get("Acceptor").unwrap()].to_string()};
            end = if strand == "+" {spt[*header.get("Acceptor").unwrap()].to_string()} else {spt[*header.get("Donnor").unwrap()].to_string()};
            
            key = format!("{} {} {} {}", contig, strand, start, end);


            // 3. case
            // never seen 
            // made in previous file
            // already seen in this file

            let mut never = false;
            let mut made = false; 
            let mut visited = false;

            if !result.contains_key(&key){
                never = true
            }
            else if result.get(&key).unwrap().sample_done.contains(&file_stem) {
                visited=true
            }
            else {
                made = true
            }

            if never{
                let counts = CountsStats::new(&spt[8..]);
                let mut junction = JunctionStats{
                        contig: contig, 
                        start: start, 
                        end: end, 
                        strand: strand,
                        ambiguous: ambi,
                        control_count: Vec::new(),
                        treat_count: Vec::new(),
                        gene_tr: HashSet::new(),
                        sample_done: HashSet::new()
                };
                junction.sample_done.insert(file_stem.clone());
                junction.gene_tr.insert(gene);
                match genotype {
                    Genotype::CONTROL => {junction.control_count.push(counts);}
                    Genotype::TREATMENT => {junction.treat_count.push(counts);}
                }
                result.insert(key, junction);
            }
            else if made{
                let counts = CountsStats::new(&spt[8..]);
                match result.get_mut(&key){
                    Some(j) => {
                        j.sample_done.insert(file_stem.clone());
                        j.gene_tr.insert(gene);
                        match genotype {
                        Genotype::CONTROL => {j.control_count.push(counts);}
                        Genotype::TREATMENT => {j.treat_count.push(counts);}
                        }
                    },
                     None => {warn!("unreachable 'made' case in parse js"); ()}
                }
                
            }
            else if visited{
                match result.get_mut(&key){
                    Some(j) => {
                        j.sample_done.insert(file_stem.clone());
                        j.gene_tr.insert(gene);
                    },
                    None => {warn!("unreachable 'visited' case in parse js"); ()}
            }
        }
    }
    Ok(())
}


pub trait Tester{

    fn to_contengency(&self) -> Vec<u64>{
        let mut res = vec![0 as u64; 4];
        for (i, e) in self.groups().iter().enumerate(){
            match e{
                Genotype::CONTROL => {
                    res[0] += self.success()[i]  as u64 ;
                    res[1] += self.failures()[i] as u64 ;
                },
                Genotype::TREATMENT => {    
                    res[2] += self.success()[i]  as u64 ;
                    res[3] +=  self.failures()[i]  as u64 ;

                }
            }
        }
        res
    }


    fn extract_counts(&self, samples: &Vec<CountsStats>, cat: &Vec<SplicingCategory>, result: &mut Vec<u32>) -> (){
        for count in samples{
            result.push(count.extract_(cat));
        }
    }


    fn format_data(&mut self, treatment: &Vec<CountsStats>, control: &Vec<CountsStats>,
                     successes_cat: &Vec<SplicingCategory>,
                      failures_cat: &Vec<SplicingCategory>) -> (){
        
        let mut success: Vec<u32> = Vec::new();
        let mut failures: Vec<u32>  = Vec::new();
        let mut groups: Vec<Genotype> = Vec::new();
        for _ in 0..control.len(){
            groups.push(Genotype::CONTROL);
        }
        for _ in 0..treatment.len(){
            groups.push(Genotype::TREATMENT);
        }

        self.extract_counts(control, successes_cat, &mut success);
        self.extract_counts(treatment, successes_cat, &mut success);
        self.extract_counts(control, failures_cat, &mut failures);
        self.extract_counts(treatment, failures_cat, &mut failures);

        *self.success_mut() = success;
        *self.failures_mut() = failures;
        *self.groups_mut() = groups;
    }

    fn get_proportion(&self) -> (u32, u32, u32, u32){

        
        let mut ctrl_suc = 0;
        let mut ctrl_fail = 0;
        let mut treat_suc = 0;
        let mut treat_fail = 0;

        for i in 0..self.groups().len(){
            match  self.groups()[i]{
                Genotype::TREATMENT => {          
                        treat_suc += self.success()[i];
                        treat_fail += self.failures()[i];
                    },
                Genotype::CONTROL => {
                        ctrl_suc += self.success()[i];
                        ctrl_fail += self.failures()[i];
                }
            } 

        }
        (ctrl_suc, ctrl_fail, treat_suc, treat_fail)
    }

    fn get_proportion_string(&self) -> (String, String, String, String){

        
        let mut ctrl_suc: Vec<String> = Vec::new();
        let mut ctrl_fail: Vec<String> = Vec::new();
        let mut treat_suc: Vec<String> = Vec::new();
        let mut treat_fail: Vec<String> = Vec::new();

        for i in 0..self.groups().len(){
            match self.groups()[i] {
                Genotype::TREATMENT => {
                    treat_suc.push(self.success()[i].to_string());
                    treat_fail.push(self.failures()[i].to_string());
                }
                Genotype::CONTROL => {
                    ctrl_suc.push(self.success()[i].to_string());
                    ctrl_fail.push(self.failures()[i].to_string());
                }
            }
        }
        (ctrl_suc.join(","), ctrl_fail.join(","), treat_suc.join(","), treat_fail.join(","))
    }

    fn test(&self, donotrun: bool) -> TestResults; // donotrun in case we just need to recover prop data but not run test i.e. ambiguous sample
    fn success(&self) -> &Vec<u32>;
    fn success_mut(&mut self) -> &mut Vec<u32>;
    fn failures(&self) -> &Vec<u32>;
    fn failures_mut(&mut self) -> &mut Vec<u32>;
    fn groups(&self) -> &Vec<Genotype>;
    fn groups_mut(&mut self) -> &mut Vec<Genotype>;
}



#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum TestStatus {
    Ok,
    QuasiPerfectSeparation,
    DimensionMistmatch,
    InvalidData,
    InsufficientObservation,
    ambiguous,
    EmptyData,
    ControlIsNull,
    TreatmentIsNull,
    ConvergenceFailed,
    NumericalInstability,
    SingularMatrix,
    PerfectSeparation,
    FisherFallBack
}


impl fmt::Display for TestStatus {
    // This trait requires `fmt` with this exact signature.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {

        match self {
            TestStatus::Ok => { write!(f, "Ok") },
            TestStatus::QuasiPerfectSeparation => { write!(f, "QuasiPerfectSeparation") },
            TestStatus::DimensionMistmatch => { write!(f, "DimensionMistmatch") },
            TestStatus::InvalidData => { write!(f, "InvalidData") },
            TestStatus::InsufficientObservation => { write!(f, "InsufficientObservation") },
            TestStatus::ambiguous => { write!(f, "ambiguous") },
            TestStatus::EmptyData => { write!(f, "EmptyData") },
            TestStatus::ControlIsNull => { write!(f, "ControlIsNull") },

            TestStatus::TreatmentIsNull => { write!(f, "TreatmentIsNull") },

            TestStatus::ConvergenceFailed => { write!(f, "ConvergenceFailed") },
            TestStatus::NumericalInstability => { write!(f, "NumericalInstability") },
            TestStatus::SingularMatrix => { write!(f, "SingularMatrix") },
            TestStatus::PerfectSeparation => { write!(f, "PerfectSeparation") },
            TestStatus::FisherFallBack => { write!(f, "FisherTest")}
        }
    }
}


impl From<LogisticRegressionError> for TestStatus {
    fn from(item: LogisticRegressionError) -> Self {
        match item{

            LogisticRegressionError::PerfectSeparation{message: _} => TestStatus::PerfectSeparation,
            
            LogisticRegressionError::InvalidProbability(f64) => TestStatus::NumericalInstability ,
            
            LogisticRegressionError::ConvergenceFailure  { iterations: _, final_norm: _ } => TestStatus::ConvergenceFailed,
        
            LogisticRegressionError::SingularMatrix(String) => TestStatus::SingularMatrix,

            LogisticRegressionError::DimensionMismatch { expected: _, got: _ } => TestStatus::DimensionMistmatch,
            
            LogisticRegressionError::InvalidData(String) => TestStatus::InvalidData,
            
            LogisticRegressionError::NumericalInstability(String) => TestStatus::NumericalInstability,

            LogisticRegressionError::EmptyData => TestStatus::EmptyData,

            LogisticRegressionError::ControlIsNull => TestStatus::ControlIsNull,

            LogisticRegressionError::TreatmentIsNull => TestStatus::TreatmentIsNull

        }
        
    }
}
#[derive(Debug)]
pub struct TestResults{

    pub control_success: u32,
    pub control_failure: u32,
    pub control_prop: Option<f32>,
    pub treatment_success: u32,
    pub treatment_failure: u32,
    pub treatment_prop: Option<f32>,
    pub status: Option<TestStatus>,
    // Only available if model fit succeeded
    pub p_value: Option<f64>,
    pub odd_ratio: Option<f64>,
    pub or_ci_lower: Option<f64>,
    pub or_ci_upper: Option<f64>,

    pub string_count: (String, String, String, String)
}

impl TestResults{
    pub fn get_empty() -> Self{
        TestResults { control_success: 0, control_failure: 0, control_prop: None,
                      treatment_success: 0, treatment_failure: 0, treatment_prop: None,
                      status: None, p_value: None,
                      odd_ratio: None, or_ci_lower: None, or_ci_upper: None, string_count: ("".to_string(), "".to_string(), "".to_string(), "".to_string())}
    }

    fn helper_(value: Option<f32>) -> String{
        match value{
            Some(v) => format!("{:e}", v),
            None => "nan".to_string()
        }
    }
    fn helper_6(value: Option<f64>) -> String{
        match value{
            Some(v) => format!("{:e}", v),
            None => "nan".to_string()
        }
    }
    fn helper_t(value: &Option<TestStatus>) -> String{
        match value{
            Some(v) => v.to_string(),
            None => "nan".to_string()
        }
    }
    
    pub fn dump_stats(&self, q_val: Option<f32>) -> Vec<String>{

        //format!("{}\t{}\t{}\t{}", 
         vec![           Self::helper_6(self.odd_ratio),
                    Self::helper_6(self.p_value),
                    Self::helper_(q_val),
                    Self::helper_t(&self.status)]
        //        )
    }
}