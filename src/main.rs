#![allow(unused)]
extern crate CigarParser;
use clap::Parser;
use rust_htslib::bam::{IndexedReader, Read};
use strand_specifier_lib::{check_flag, LibType};
use CigarParser::cigar::Cigar;
//use rust_htslib::errors::Error;
use rust_htslib::bam::record::Record;
use std::collections::HashMap;
use std::fs::File;
use std::hash::Hash;
use std::io::prelude::*;
use std::io::BufWriter;

use std::str::FromStr;
mod common;
use std::fs;
use std::str;
//use crate::common::utils::ReadAssign;
use crate::common::it_approches::{
    dump_tree_to_cat_results, gtf_to_tree, update_tree_with_bamfile,
};
use crate::common::point::{read_gtf, InsideCounter, PointContainer};
use crate::common::read_record::file_to_table;
use crate::common::utils;
use crate::common::utils::{ExonType, ReadAssign};

// TODO add specific subcommand to retrieve only specific reads;
// TODO add LibType
// TODO refactor to use a outFilePrefix and always output both .cat and .table file, and output read in the correspoding files.

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Name of Input file
    #[arg(short, long, required = true)]
    input: String,
    /// Name of Output file
    #[arg(short, long, required = true)]
    output: String,
    /// Name of GTF Input file define the feature to look at
    /// (v1) only consider feature annotated as exon
    /// if you use output_write_read with the whole genome the output can be very large,
    /// you may want to subset genes / features you are interested in.
    #[arg(short, long, required = true)]
    gtf: String,
    /// size of overhang
    #[arg(short, long, default_value_t = 1)]
    overhang: i64,
    /// path to a file (must not exist) Uses if you want to output the reads with their category.
    /// by default output all reads, this behaviour can be change using flags: clipped...
    #[arg(long)]
    output_write_read: Option<String>,
    /// path to a file (must not exist) Uses if you want to output the table of category
    #[arg(long)]
    table: Option<String>,

    #[arg(long, default_value_t = 0)]
    flag_in: u16,
    #[arg(long, default_value_t = 3840)]
    flag_out: u16,
    #[arg(long, default_value_t = 13)]
    mapq: u8,
    /// only do something if --output-write-read is set
    #[clap(long, short, action)]
    clipped: bool,
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
    output_write_read: Option<String>,
    clipped: bool,
) -> () {

    let bam_file = bam_input;

    let gtf_file = gtf;

    let mut output_read_stream: Option<BufWriter<File>> = None; 
    if let Some(file_path) = output_write_read {
        let file_ =
            File::create_new(file_path.clone()).expect("read output file should not exist.");
        fs::write(file_path, "read_name\tcig\tflag\taln_start\tread_assign\tfeature.pos\tnext_exon\tfeature.exon_type\tfeature.strand\tsequence\n".as_bytes()).expect("Unable to write file");
        output_read_stream = Some(BufWriter::new(file_));
    }

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
        &mut output_read_stream,
        clipped,
    );

    dump_tree_to_cat_results(&hash_tree, &output);
    //stream.flush().unwrap();
    if let Some(ref mut out_stream)= output_read_stream {
        out_stream
            .flush()
            .unwrap();

      
    }
}

fn main() {
    let args = Args::parse();

    let output = args.output;
    let mut clipped = false;
    if args.clipped {
        clipped = true;
    }

    main_loop(
        output.clone(),
        args.gtf.clone(),
        args.input,
        args.overhang,
        args.flag_in,
        args.flag_out,
        args.mapq,
        args.output_write_read,
        clipped,
    );

    if let Some(table) = args.table {
        let file = File::create_new(table.clone())
            .unwrap_or_else(|_| panic!("output file {} should not exist.", &table));//expect(&format!("output file {} should not exist.", &table));
        let mut stream = BufWriter::new(file);
        println!("table");
        let _ = stream.write("contig\tgene_name\ttranscript_name\texon_number\tambiguous\tstrand\tpos\tnext\texon_type\tspliced\tunspliced\tclipped\texon_intron\texon_other\tskipped\twrong_strand\te_isoform\n".as_bytes());
        file_to_table(output.clone(), &mut stream, args.gtf.as_str());
    }
}
//            skipped,
//wrong_strand,
//e_isoform
