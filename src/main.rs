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
use crate::common::gtf_::{gtf_to_hashmap, get_junction_from_gtf};
use crate::common::it_intron::{interval_tree_from_gtfmap, update_tree_from_bam, dump_tree_to_cat_results};
//use crate::common::it_approches::{
//    dump_tree_to_cat_results, gtf_to_tree, update_tree_with_bamfile,
//};
//use crate::common::point::{read_gtf, InsideCounter, PointContainer};
use crate::common::read_record::file_to_table;
use crate::common::utils;
use crate::common::junction_file::junction_file_from_table;
use crate::common::utils::{ExonType, ReadAssign, update_read_to_write_handle, ReadToWriteHandle, ReadsToWrite};
mod splicing_efficiency;


// TODO add LibType

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Name of Input file
    #[arg(short, long, required = true)]
    input: String,
    /// Prefix name  to be used for Output file
    #[arg(short, long, required = true)]
    output_file_prefix: String,
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
    /// space separated list of the annotated read you want to extract; i.e. all clipped read or all spliced read ...
    #[clap(long, value_parser, value_delimiter = ' ', num_args = 1..)]
    readToWrite: Vec<ReadsToWrite>,
    /// space separated list the column to use for "unspliced" for the splicing defect table.
    /// you can regenrate this using the splicing_efficiency exe
    /// What to consider as unspliced? spliced, unspliced, clipped, exon_other, skipped,
    /// wrong_strand, isoform\n
    /// by default only use "-u unspliced" ->  unspliced (readthrough) reads \n
    /// to use unspliced and clipped : "-u unspliced clipped" 
    #[clap(long, value_parser, default_value = "unspliced", value_delimiter = ' ', num_args = 1..)]
    unspliced_def: Vec<String>,

    /// What to consider as spliced? spliced, unspliced, clipped, exon_other, skipped,
    /// wrong_strand, isoform\n
    /// by default only use "-u spliced" -> spliced (readthrough) reads \n
    /// to use spliced and isoform : "-u spliced isoform"
    #[clap(long, value_parser, default_value = "spliced", value_delimiter = ' ', num_args = 1..)]
    spliced_def: Vec<String>,

    /// Librairy types used for the RNAseq most modern stranded RNAseq are frFirstStrand which is the default value.
    /// acceptable value: frFirstStrand, frSecondStrand, fFirstStrand, fSecondStrand, ffFirstStrand, ffSecondStrand, rfFirstStrand,
    ///  rfSecondStrand, rFirstStrand, rSecondStrand, Unstranded, PairedUnstranded
    #[clap(long, value_parser, default_value = "frFirstStrand")]
    libtype: LibType,
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
    output_write_read_handle: &mut ReadToWriteHandle, 
    librairy_type: LibType
) -> () {

    let bam_file = bam_input;

    let gtf_file = gtf;
    println!("parse gtf file");
    let gtf_hashmap = gtf_to_hashmap(&gtf_file).expect("failed to parse gtf");
    println!("gtf file parsed");
    println!("building the interval tree");
    let mut hash_tree = interval_tree_from_gtfmap(&gtf_hashmap).expect("failed to generate the hash tree from gtf");
    println!("tree built");
    
    // parse the gtf and return a hashmap<chromosome> -> intervalTree(intron(start, end), associated_data(gene_name...))
    //let mut hash_tree = gtf_to_tree(gtf_file.as_str()).unwrap();
    let junction_ = get_junction_from_gtf(&gtf_file, &librairy_type);
    println!("valid junctions loaded");

    update_tree_from_bam(
        &mut hash_tree,
        &bam_file,
        librairy_type, //LibType::frFirstStrand,
        overhang,
        flag_in,
        flag_out,
        mapq,
        output_write_read_handle,
        &junction_
    );

    dump_tree_to_cat_results(&hash_tree, &output);
    output_write_read_handle.flush();

}

fn main() {
    let args = Args::parse();

    match args.libtype{
        LibType::Invalid => panic!("invalid librairy type"),
        _ => ()
    }
    let mut output_file_prefix = args.output_file_prefix;
    let table  = format!("{}{}", output_file_prefix, ".table");
    let output = format!("{}{}", output_file_prefix, ".cat");
    let splicing_defect = format!("{}{}", output_file_prefix, ".sd");
    let junction_file = format!("{}{}", output_file_prefix, ".junction");
    let mut clipped = false;


    let header_reads_handle = "read_name\tcig\tflag\taln_start\tread_assign\tfeature.pos\tnext_exon\tfeature.exon_type\tfeature.strand\tsequence\n".as_bytes();//.expect("Unable to write file");
    let mut readouthandle = ReadToWriteHandle::new();
    update_read_to_write_handle(&mut readouthandle, args.readToWrite, header_reads_handle, &output_file_prefix);



    main_loop(
        output.clone(),
        args.gtf.clone(),
        args.input,
        args.overhang,
        args.flag_in,
        args.flag_out,
        args.mapq,
        &mut readouthandle,
        args.libtype
    );


    let file = File::create_new(table.clone())
        .unwrap_or_else(|_| panic!("output file {} should not exist.", &table));//expect(&format!("output file {} should not exist.", &table));
    let mut stream = BufWriter::new(file);
    let _ = stream.write("contig\tgene_name\ttranscript_name\texon_number\tambiguous\tstrand\tpos\tnext\texon_type\tspliced\tunspliced\tclipped\texon_other\tskipped\twrong_strand\te_isoform\n".as_bytes());
    file_to_table(output.clone(), &mut stream, args.gtf.as_str(), args.libtype);
    junction_file_from_table(&table, &junction_file);
    splicing_efficiency::to_se_from_table(&table, &splicing_defect, args.spliced_def, args.unspliced_def);
    
}


