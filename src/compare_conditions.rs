#![allow(warnings)]

use std::collections::{HashMap, HashSet};
use std::fmt::format;
use std::path::Path;
use std::sync::Arc;
use std::{fs, result};
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::prelude::*;
use std::io::{BufReader, BufWriter};
use std::time::Instant;

use clap::CommandFactory;
use clap::{Parser, Subcommand};
use std::path::PathBuf;


use rayon::ThreadPoolBuilder;
use rayon::prelude::*;

mod stat_common;
use stat_common::errors::LogisticRegressionError;
use stat_common::common::{Tester, parse_js_file, Genotype, JunctionStats, SplicingCategory};
use stat_common::glm_logistic::{GLM};

mod common;

use common::error::OmniError;

use nalgebra::{DMatrix, DVector};
use statrs::distribution::{ChiSquared, ContinuousCDF};
use adjustp::{adjust, Procedure};

use crate::common::junction_file;
use crate::stat_common::common::{TestResults, TestStatus};

use flexi_logger::{FileSpec, Logger, WriteMode};
use log::{debug, error, info, trace, warn};



/// keep Fisher in python?
pub struct Fisher{
    
}
//}
// If this ahppend to be to expansive I could order the indice of qvalue and then retorder both list accordingly.
fn sort_by_f32_copy<T>(data: &mut Vec<T>, scores: &mut Vec<f32>) {
    // Create pairs, sort them, then unzip
    let mut pairs: Vec<(T, f32)> = data.drain(..).zip(scores.drain(..)).collect();
    
    pairs.sort_by(|a, b| {
        a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal)
    });
    
    let (sorted_data, sorted_scores): (Vec<T>, Vec<f32>) = pairs.into_iter().unzip();
    
    *data = sorted_data;
    *scores = sorted_scores;
}

fn sort_by_f32_permutation<T: std::fmt::Debug>(data: &mut Vec<T>, scores: &mut Vec<f32>) {

    let mut indices: Vec<usize> =  (0..scores.len()).collect();
    let mut rank = vec![0usize; scores.len()];
    {
        indices.sort_by(|&a, &b| scores[a].partial_cmp(&scores[b]).unwrap_or_else(|| std::cmp::Ordering::Equal));
        
        
        for (sorted_pos, &orig_idx) in indices.iter().enumerate() {
            rank[orig_idx] = sorted_pos;
        }
    }
        
    let n = rank.len();
    //argsort
    
    for i in 0..n{
        let mut j = i;
        while rank[j] != usize::MAX && rank[j] != j{
            data.swap(i, rank[j]);
            scores.swap(i, rank[j]);
            let next = rank[j];
            rank[j] = usize::MAX;
            j = next;

        }
    }
}


fn parse_results_update_vec(vec_r: &Vec<(&JunctionStats, TestResults)>, result: &mut Vec<Vec<String>>, q_value: Option<Vec<f32>>){
    let mut value: (String, String, String, String);
    for (i, (j, t)) in vec_r.into_iter().enumerate(){
        
        let mut f = j.get_pos_string();
        match q_value{
            Some(ref q) => f.extend_from_slice(&t.dump_stats(Some(q[i]))),
            None =>  f.extend_from_slice(&t.dump_stats(None))

        }

        value = t.string_count.clone();

        f.push(value.2);
        f.push(value.3);
        f.push(  match t.treatment_prop {
            Some(c) => c.to_string(),
            None => "nan".to_string(),}
        );

        f.push(value.0);
        f.push(value.1);
        f.push(   match t.control_prop {
                        Some(c) => c.to_string(),
                        None => "nan".to_string()}
        );



        f.push(j.gene_tr.iter().map(|x| x.to_owned()).collect::<Vec<String>>().join(";"));
        result.push(f)
    }
}


fn run_one_test(junction: &HashMap<String, JunctionStats>, successes_cat: Vec<SplicingCategory>,
                 failures_cat: Vec<SplicingCategory>, out_file_path: &str, ambi: bool) -> Result<(), Box<dyn std::error::Error + Send + Sync>>{

    
    info!("Starting test!");
    let mut vec_ok: Vec<(&JunctionStats, TestResults)> = Vec::new();
    let mut vec_err: Vec<(&JunctionStats, TestResults)> = Vec::new();
    let mut x: TestResults; 
    for (k, j ) in junction.iter(){
        let mut glm = GLM::new(  &j.control_count,
                                        &j.treat_count,
                                                &successes_cat,
                                                &failures_cat,
                                            k.to_owned());
        if ambi == false && j.ambiguous == true{
            x = glm.test(true);
        }
        else {
            x = glm.test(false);
        }
        if x.p_value.is_some(){
            vec_ok.push((&j, x));
        }
        else{
            match x.status {
                Some(TestStatus::TreatmentIsNull) | Some(TestStatus::ControlIsNull) => {
                    vec_err.push((&j, x));
                } ,
                _ => (),
                None => ()
            }

        }
    }
    
    let pval = vec_ok.iter().map(|x| x.1.p_value.unwrap() as f32).collect::<Vec<f32>>();
    let mut qvalues = adjust(&pval, Procedure::BenjaminiHochberg);

    let mut final_: Vec<Vec<String>> = Vec::new();
    let mut value: (String, String, String, String);
    
    sort_by_f32_permutation(&mut vec_ok, &mut qvalues);

    parse_results_update_vec(&vec_ok, &mut final_, Some(qvalues));
    parse_results_update_vec(&vec_err, &mut final_, None);


    let mut out_file_open =
        File::create_new(out_file_path.clone()) //presorted out_file.clone()
            .unwrap_or_else(|_| panic!("output file {} should not exist.", &out_file_path)); //expect(&format!("output file {} should not exist.", &table));
    let mut out_stream = BufWriter::new(out_file_open);

    out_stream.write("#success: spliced\n".to_string().as_bytes());
    out_stream.write("#failures: unspliced\n".to_string().as_bytes());
    out_stream.write("#test: GLM Binomial: (successes + failures) ~ group\n".to_string().as_bytes());
    let header = vec!["chr", "strand", "start", "end", "statistic", "p_value", "q_value", "status", "control_success", "control_failures", "control_ratio", "treatment_success", "treatment_failures", "treatment_ratio",  "gene_transcript_intron"];
    out_stream.write(format!("{}\n", header.join("\t")).as_bytes());
    for e in final_{
        out_stream.write(format!("{}\n", e.join("\t")).as_bytes())?;
    }

    let _ = out_stream.flush();
    Ok(())
}




#[derive(Parser)]
#[command(name = "compare")]
#[command(about = "Allows to compare condition from Omnisplice junction file", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Run a single comparison
    Run {
        /// Output file prefix path
        #[arg(short, long,  required = true)]
        outfile: PathBuf,

        /// Control Condition
        #[arg(short, long, num_args = 1.., required = true)]
        control_files: Vec<PathBuf>,

        /// Control Condition
        #[arg(short, long, num_args = 1.., required = true)]
        treatment_files: Vec<PathBuf>,

        /// Splicing type considered as good
        /// must be one or more of the follwing: Spliced Unspliced  Clipped Exon_other Skipped SkippedUnrelated Wrong_strand  E_isoform
        #[arg(short, long, num_args = 1.., required = true)]
        splicing_ok: Vec<String>,

        /// Splicing type considered as failures
        /// must be one or more of the follwing: Spliced Unspliced       Clipped Exon_other      Skipped SkippedUnrelated        Wrong_strand    E_isoform
        #[arg(short, long, num_args = 1.., required = true)]
        splicing_fail: Vec<String>,

       /// Do you want to consider ambigious junction (overlaping exon)
       #[arg( long,)]
       ambigious: bool,


    },
    /// Run all comparisons against splices
    RunAll {
        /// Output file path
        #[arg(short, long, required = true)]
        outfile_prefix: String,

        /// Control Condition
        #[arg(short, long, num_args = 1.., required = true)]
        control_files: Vec<PathBuf>,

        /// Control Condition
        #[arg(short, long, num_args = 1.., required = true)]
        treatment_files: Vec<PathBuf>,
    },
}

fn parse_cat(input: Vec<String>) -> Result<Vec<SplicingCategory>, &'static str >{
    let mut res = Vec::new();
    for e in input{
        if e.trim().is_empty(){
            continue
        }
        res.push(SplicingCategory::try_from(e.trim())?)
    }
    Ok(res)
}

fn main() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {


    Logger::try_with_str("error").unwrap().start().unwrap();


    let cli = Cli::parse();

    match cli.command {
        Commands::Run { outfile, control_files, treatment_files, splicing_ok, splicing_fail, ambigious } => {
        
            println!("Running comparison, output: {:?}", outfile);
            let control = parse_cat(splicing_ok)?;
            let treatment = parse_cat(splicing_fail)?;
    
            let mut res:  HashMap<String, JunctionStats> = HashMap::with_capacity(1_000_000);

            for file in control_files{
                    info!("parsing {:?}", file);
                    parse_js_file(file.to_str().unwrap(), &mut res, Genotype::CONTROL).unwrap();
                    info!("done reading");
                }

            for file in treatment_files{
                    info!("parsing {:?}", file);
                    parse_js_file(file.to_str().unwrap(), &mut res, Genotype::TREATMENT).unwrap();
                    info!("done reading");
                }


        run_one_test( &res, control,
                 treatment,  outfile.to_str().unwrap(),  ambigious)?;
                //run_one_test
            // Your run logic here
        }
        Commands::RunAll { control_files, treatment_files, outfile_prefix } => {
            println!("Running all single comparisons, output: {:?}", outfile_prefix);
                
            let mut res:  HashMap<String, JunctionStats> = HashMap::with_capacity(1_000_000);

            for file in control_files{
                    info!("parsing {:?}", file);
                    parse_js_file(file.to_str().unwrap(), &mut res, Genotype::CONTROL).unwrap();
                    info!("done reading");
                }

            for file in treatment_files{
                    info!("parsing {:?}", file);
                    parse_js_file(file.to_str().unwrap(), &mut res, Genotype::TREATMENT).unwrap();
                    info!("done reading");
                }

            ThreadPoolBuilder::new()
        .num_threads(6)          // ← set the limit here
        .build_global()
        .expect("Failed to initialise Rayon thread‑pool");

    info!("All junction file parsed");
    let shared = Arc::new(res);

    let mut p = Path::new(&outfile_prefix).to_path_buf();


    let mut jobs: Vec<(Vec<SplicingCategory>, Vec<SplicingCategory>, &str, bool)> = Vec::new();
    let mut p = Path::new(&outfile_prefix).to_path_buf();

    let _ = p.set_extension("Unspliced.tsv");
    jobs.push((vec![SplicingCategory::Spliced], vec![SplicingCategory::Unspliced],
       p.to_str().unwrap() , false));

    let mut p = Path::new(&outfile_prefix).to_path_buf();
    let _ = p.set_extension("WrongStrand.tsv");
    jobs.push((vec![SplicingCategory::Spliced], vec![SplicingCategory::WrongStrand],
       p.to_str().unwrap() , false));

    let mut p = Path::new(&outfile_prefix).to_path_buf();
    let _ = p.set_extension("Skipped.tsv");
    jobs.push((vec![SplicingCategory::Spliced], vec![SplicingCategory::Skipped],
       p.to_str().unwrap() , false));
    
    let mut p = Path::new(&outfile_prefix).to_path_buf();
    let _ = p.set_extension("SkippedUnrelated.tsv");
    jobs.push((vec![SplicingCategory::Spliced], vec![SplicingCategory::SkippedUnrelated],
       p.to_str().unwrap() , false));
    
    let mut p = Path::new(&outfile_prefix).to_path_buf();
    let _ = p.set_extension("Clipped.tsv");
    jobs.push((vec![SplicingCategory::Spliced], vec![SplicingCategory::Clipped],
       p.to_str().unwrap() , false));

    let mut p = Path::new(&outfile_prefix).to_path_buf();
    let _ = p.set_extension("ExonOther.tsv");
    jobs.push((vec![SplicingCategory::Spliced], vec![SplicingCategory::ExonOther],
       p.to_str().unwrap() , false));
       
    let mut p = Path::new(&outfile_prefix).to_path_buf();
    let _ = p.set_extension("Isoform.tsv");
    jobs.push((vec![SplicingCategory::Spliced], vec![SplicingCategory::EIsoform],
       p.to_str().unwrap() , false));


    //println!("{:?}", );
    let now = Instant::now();
    let ok = jobs.into_par_iter().
    try_for_each(|(a, b, c, d)| {
        let shared_ref = Arc::clone(&shared);
        run_one_test(&shared_ref, a, b, c, d)
    });
    match ok {
        Ok(()) => println!("All four tests finished successfully."),
        Err(e) => eprintln!("A test failed: {}", e),
    }

    let elapsed_time = now.elapsed();
    println!("Running slow_function() took {} seconds.", elapsed_time.as_secs());   

        }
    }

    Ok(())
}
    /* let junction_file_ctrl = vec!["/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/TDP43/omniSplice/SRR22002170_R1_001.out.sorted.bam.junctions",
                                         "/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/TDP43/omniSplice/SRR22002171_R1_001.out.sorted.bam.junctions", 
                                        "/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/TDP43/omniSplice/SRR22002172_R1_001.out.sorted.bam.junctions"];
    let junction_file_treat = vec!["/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/TDP43/omniSplice/SRR22002167_R1_001.out.sorted.bam.junctions",
                                        "/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/TDP43/omniSplice/SRR22002168_R1_001.out.sorted.bam.junctions",
                                        "/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/TDP43/omniSplice/SRR22002169_R1_001.out.sorted.bam.junctions"];

    let mut res:  HashMap<String, JunctionStats> = HashMap::with_capacity(1_000_000);




    for file in junction_file_ctrl{
        info!("parsing {}", file);
        parse_js_file(file, &mut res, Genotype::CONTROL).unwrap();
        info!("done reading");
    }

    for file in junction_file_treat{
        info!("parsing {}", file);
        parse_js_file(file, &mut res, Genotype::TREATMENT).unwrap();
        info!("done reading");
    }


    ThreadPoolBuilder::new()
        .num_threads(6)          // ← set the limit here
        .build_global()
        .expect("Failed to initialise Rayon thread‑pool");

    info!("All junction file parsed");
    let shared = Arc::new(res);

    let jobs: Vec<(Vec<SplicingCategory>, Vec<SplicingCategory>, &str, bool)> = vec![
    (vec![SplicingCategory::Spliced], vec![SplicingCategory::Unspliced], "/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/TDP43/OS_run_oh3/test/tdp43_omni_female_Cortex.Unspliced.tsv", false),
    (vec![SplicingCategory::Spliced], vec![SplicingCategory::WrongStrand], "/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/TDP43/OS_run_oh3/test/tdp43_omni_female_Cortex.WrongStrand.tsv", true),
    (vec![SplicingCategory::Spliced], vec![SplicingCategory::Skipped], "/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/TDP43/OS_run_oh3/test/tdp43_omni_female_Cortex.Skipped.tsv", true),
    (vec![SplicingCategory::Spliced], vec![SplicingCategory::SkippedUnrelated], "/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/TDP43/OS_run_oh3/test/tdp43_omni_female_Cortex.SkippedUnrelated.tsv", true),
    (vec![SplicingCategory::Spliced], vec![SplicingCategory::Clipped], "/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/TDP43/OS_run_oh3/test/tdp43_omni_female_Cortex.Clipped.tsv", true),
    (vec![SplicingCategory::Spliced], vec![SplicingCategory::ExonOther], "/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/TDP43/OS_run_oh3/test/tdp43_omni_female_Cortex.ExonOther.tsv", true),
    (vec![SplicingCategory::Spliced], vec![SplicingCategory::EIsoform], "/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/TDP43/OS_run_oh3/test/tdp43_omni_female_Cortex.Isoform.tsv", true),
    ];
    let now = Instant::now();
    let ok = jobs.into_par_iter().
    try_for_each(|(a, b, c, d)| {
        let shared_ref = Arc::clone(&shared);
        run_one_test(&shared_ref, a, b, c, d)
    });
    match ok {
        Ok(()) => println!("All four tests finished successfully."),
        Err(e) => eprintln!("A test failed: {}", e),
    }

    let elapsed_time = now.elapsed();
    println!("Running slow_function() took {} seconds.", elapsed_time.as_secs());   */

       
