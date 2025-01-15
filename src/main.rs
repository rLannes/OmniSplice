#![allow(unused)]
extern crate CigarParser;
use clap::Parser;

use rust_htslib::bam::{IndexedReader, Read};
use strand_specifier_lib::{check_flag, LibType};
use CigarParser::cigar::Cigar;
//use rust_htslib::errors::Error;
use rust_htslib::bam::record::Record;
use std::collections::HashMap;
use std::fmt::format;
use std::fs::File;
use std::hash::Hash;
use std::io::prelude::*;
use std::io::BufWriter;

use std::path::{Path, PathBuf};
use std::ffi::OsStr;
use std::str::FromStr;
mod common;
use std::fs;
use std::str;
use std::convert::From;

//use crate::common::utils::ReadAssign;
use crate::common::it_approches::{
    dump_tree_to_cat_results, gtf_to_tree, update_tree_with_bamfile,
};
use crate::common::point::{read_gtf, InsideCounter, PointContainer};
use crate::common::read_record::file_to_table;
use crate::common::utils;
use crate::common::utils::{ExonType, ReadAssign, update_ReadToWriteHandle, ReadToWriteHandle, ReadsToWrite};
mod splicing_efficiency;


// TODO add LibType



// impl From<ReadAssign> for ReadsToWrite{
//     fn from(value: ReadAssign) -> Self {
//         match value {
//             ReadAssign::Empty => ReadsToWrite::Empty,
//             ReadAssign::FailPosFilter => ReadsToWrite::FailPosFilter,
//             ReadAssign::FailQc => ReadsToWrite::FailQc,
//             ReadAssign::ReadThrough => ReadsToWrite::ReadThrough,
//             ReadAssign::ReadJunction(n, m) => ReadsToWrite::ReadJunction,
//             ReadAssign::Unexpected => ReadsToWrite::Unexpected,
//             ReadAssign::EmptyPileup => ReadsToWrite::EmptyPileup,
//             ReadAssign::OverhangFail => ReadsToWrite::OverhangFail,
//             ReadAssign::SoftClipped => ReadsToWrite::SoftClipped,
//             ReadAssign::Skipped(n, m) => ReadsToWrite::Skipped,
//             ReadAssign::WrongStrand => ReadsToWrite::WrongStrand            
//         }
//     }
// }

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Name of Input file
    #[arg(short, long, required = true)]
    input: String,

    /// Prefix name  to be used for Output file
    #[arg(short, long, required = true)]
    outputFilePrefix: String,
    /// Name of GTF Input file define the feature to look at
    /// (v1) only consider feature annotated as exon
    /// if you use output_write_read with the whole genome the output can be very large,
    /// you may want to subset genes / features you are interested in.
    #[arg(short, long, required = true)]
    gtf: String,
    /// size of overhang
    #[arg(long, default_value_t = 1)]
    overhang: i64,
    /// path to a file (must not exist) Uses if you want to output the reads with their category.
    /// by default output all reads, this behaviour can be change using flags: clipped...
    #[arg(long)]
    output_write_read: Option<String>,

    #[arg(long, default_value_t = 0)]
    flag_in: u16,
    #[arg(long, default_value_t = 3840)]
    flag_out: u16,
    #[arg(long, default_value_t = 13)]
    mapq: u8,
    /// space separated list of the annotated read you xant to extract; i.e. all clipped read or all spliced read ...
    #[clap(long, value_parser, value_delimiter = ' ', num_args = 1..)]
    readToWrite: Vec<ReadsToWrite>,
    /// space separated list the column to use for "unspliced" for the splicing defect table.
    /// you can regenrate this using the splicing_efficiency exe
    /// What to consider as unspliced? unspliced: 9, clipped: 10, exon_intron: 11, exon_other: 12, skipped: 13, wrong_strand:14\n
    /// by default only use "-u 9" ->  unspliced (readthrough) reads \n
    /// to use unspliced and clipped : "-u 9 10" 
    #[clap(long, value_parser, default_value = "9", value_delimiter = ' ', num_args = 1..)]
    unspliced_def: Vec<usize>,
    /// splicing defect only: only consider donnor exons and not the acceptor exon
    #[arg(long, default_value_t = false, action)]
    one_side : bool,

}

/// This run the core of the program, will parse a gtf and a bam file and write a category file and if requested a read file.
fn main_loop(
    output: String,
    gtf: String,
    bam_input: String,
    overhang: i64,
    flag_in: u16,
    flag_out: u16,
    mapq: u8,
    output_write_read_handle: &mut ReadToWriteHandle// Option<String>,
    //clipped: bool,
) -> () {

    let bam_file = bam_input;

    let gtf_file = gtf;

    /*let mut output_read_stream: Option<BufWriter<File>> = None; 
    if let Some(file_path) = output_write_read {
        let file_ =
            File::create_new(file_path.clone()).expect("read output file should not exist.");
        fs::write(file_path, "read_name\tcig\tflag\taln_start\tread_assign\tfeature.pos\tnext_exon\tfeature.exon_type\tfeature.strand\tsequence\n".as_bytes()).expect("Unable to write file");
        output_read_stream = Some(BufWriter::new(file_));
    }*/

    // parse the gtf and return a hashmap<chromosome> -> intervalTree(exon(start, end), associated_data(gene_name...))
    let mut hash_tree = gtf_to_tree(gtf_file.as_str()).unwrap();
    
    update_tree_with_bamfile(
        &mut hash_tree,
        &bam_file,
        LibType::frFirstStrand,
        overhang,
        flag_in,
        flag_out,
        mapq,
        output_write_read_handle,
    );

    dump_tree_to_cat_results(&hash_tree, &output);
    output_write_read_handle.flush();

}

fn main() {
    let args = Args::parse();
    println!("{:#?}", args);
    let mut outputFilePrefix = args.outputFilePrefix;
    let table  = format!("{}{}", outputFilePrefix, ".table");
    let output = format!("{}{}", outputFilePrefix, ".cat");
    let splicing_defect = format!("{}{}", outputFilePrefix, ".sd");
    let mut clipped = false;


    // very verbose but this allow us to contains all those potential files handle in one struct to moove around
    let headerReadsHandle = "read_name\tcig\tflag\taln_start\tread_assign\tfeature.pos\tnext_exon\tfeature.exon_type\tfeature.strand\tsequence\n".as_bytes();//.expect("Unable to write file");
    let mut readouthandle = ReadToWriteHandle::new();
    update_ReadToWriteHandle(&mut readouthandle, args.readToWrite, headerReadsHandle, &outputFilePrefix);

    main_loop(
        output.clone(),
        args.gtf.clone(),
        args.input,
        args.overhang,
        args.flag_in,
        args.flag_out,
        args.mapq,
        &mut readouthandle,
    );


    let file = File::create_new(table.clone())
        .unwrap_or_else(|_| panic!("output file {} should not exist.", &table));//expect(&format!("output file {} should not exist.", &table));
    let mut stream = BufWriter::new(file);
    println!("table");
    let _ = stream.write("contig\tgene_name\ttranscript_name\texon_number\tambiguous\tstrand\tpos\tnext\texon_type\tspliced\tunspliced\tclipped\texon_intron\texon_other\tskipped\twrong_strand\te_isoform\n".as_bytes());
    file_to_table(output.clone(), &mut stream, args.gtf.as_str());

    splicing_efficiency::to_se_from_table(&table, &splicing_defect, args.unspliced_def, ! args.one_side);



}


