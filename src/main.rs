
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
use std::io::{BufWriter};

use std::str::FromStr;
mod common;
use std::fs;
use std::str;
//use crate::common::utils::ReadAssign;
use crate::common::point::{read_gtf, PointContainer, InsideCounter};
use crate::common::read_record::file_to_table;
use crate::common::utils::{ExonType, ReadAssign};
use crate::common::utils;


fn parse_bam<T>  (
    // required parameter
    bam_file: &str,
    library_type: LibType,
    final_hash: &mut HashMap<String, PointContainer<utils::ReadAssign>>,
    overhang: i64,

    // TODO not yet implemented QC optional option (None for no filter) // TODO
    // default 0
    flag_in: u16,
    // default 1024(DUPLICATE) + 256 (NOT_PRIMARY_ALN) + 2048(SUPPLEMENTARY) + 512(FAIL_QC)
    flag_out: u16,
    // mapq (must be) >=   default 13
    mapq: u8,
    out_file_read_buffer: &mut Option<BufWriter<File>>,
    clipped: bool
) -> () {
    //let mut p_to_check: i64;
    let mut counter: i64;
    let mut record: Record;
    let mut pos_s: i64;
    let mut pos_e: i64;
    let mut cig: Cigar;
    let mut flag: u16;
    //let mut strand: Strand = Strand::Plus;
    //let mut readname: String = "".to_string();
    //let mut readnameset = HashSet::new();
    //let mut flag_test_strand = false;
    //let mut read_strand: Strand = Strand::Plus;

    let mut bam = IndexedReader::from_path(bam_file).unwrap();
    //let header = bam::Header::from_template(bam.header());
    // we gonna iterate thtough all contigs!
    //let mut contig: String = "".to_string();

    for (contig, vec_points) in final_hash.iter_mut() {
        counter = 0;
        println!("Contig: {}", contig);

        match bam.fetch(&contig) {
            Ok(_) => (),
            Err(_) => {
                println!("WARNING {} not found in the bam file check", contig);
                continue;
            }
        }
        for r in bam.records() {
            counter += 1;
            if counter % 1_000_000 == 0 {
                println!("Contig: {}; {} reads done", contig, counter);
            }
            record = r.unwrap();
            pos_s = record.pos();
            cig = Cigar::from_str(&record.cigar().to_string()).unwrap();
            pos_e = cig.get_end_of_aln(&pos_s);
            flag = record.flags();

            // QC
            if (!check_flag(flag, flag_in, flag_out)) || (record.mapq() < mapq) {
                continue;
            }

            if let Some(read_strand) = library_type.get_strand(flag) {
                let _ = vec_points.parse_reads(
                    &contig,
                    pos_s,
                    pos_e,
                    &cig,
                    &flag,
                    &read_strand,
                    overhang,
                    out_file_read_buffer,
                    &record,
                    clipped
                );
            }
            //read_strand = library_type.get_strand(flag).expect(&format!("LibType: {:?} {}", library_type,  flag) );
            //vec_points.parse_reads(pos_s, pos_e, &cig, &flag, &read_strand, overhang);
        }
    }
}


// TODO add specific subcommand to retrieve only specific reads;
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
    /// you may want to subset gene / feature you are interested in.
    #[arg(short, long, required = true)]
    gtf: String,
    /// size of overhang
    #[arg(short, long, default_value_t = 1)]
    overhang: i64,
    /// path to a file (must not exist) Uses if you want to output the reads with their category.
    /// by default output all reads, this behaviour can be change using flags: junction, clipped...
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
    clipped: bool
}

fn main_loop(
    output: String,
    gtf: String,
    bam_input: String,
    overhang: i64,
    flag_in:u16,
    flag_out:u16,
    mapq:u8,
    output_write_read: Option<String>,
    clipped: bool
) -> () {
    let file = File::create_new(output).expect("output file should not exist.");
    let mut stream = BufWriter::new(file);

    let bam_file = bam_input;

    let gtf_file = gtf;
    let overhang = overhang;


    let mut results = read_gtf(&gtf_file).unwrap();

    let mut output_read_stream: Option<BufWriter<File>> = None; // args.output_write_read;
    if let Some(file_path) = output_write_read {
        let file_ =
            File::create_new(file_path.clone()).expect("read output file should not exist.");
        fs::write(file_path, "read_name\tcig\tflag\taln_start\tread_assign\tfeature.pos\tfeature.exon_type\tfeature.strand\tsequence\n".as_bytes()).expect("Unable to write file");
        output_read_stream = Some(BufWriter::new(file_));
    }


    parse_bam::<utils::ReadAssign>(
        &bam_file,
        LibType::frFirstStrand,
        &mut results,
        overhang as i64,
        flag_in,
        flag_out,
        mapq,
        &mut output_read_stream,
        clipped
    );

    // Improovment better sort
    for (contig, vec_point) in  results.iter_mut() {
        let mut sorted_point = &mut vec_point.points;
        sorted_point.sort_unstable_by_key(|item| (item.transcript_id.clone(), item.pos));
        for point in sorted_point.iter(){
        //for point in vec_point.iter() {
            if point.counter.is_empty() {
                let _ = stream.write(
                    format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\tEmpty\t0\n",
                        contig,
                        point.pos,
                        point.gene_name,
                        point.transcript_id,
                        point.strand,
                        point.exon_type,
                    )
                    .as_bytes(),
                );
            } else {
                for (k, v) in point.counter.iter() {
                    let _ = stream.write(
                        format!(
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                            contig,
                            point.pos,
                            point.gene_name,
                            point.transcript_id,
                            point.strand,
                            point.exon_type,
                            k,
                            v
                        )
                        .as_bytes(),
                    );
                }
            }
        }
    }
    stream.flush().unwrap();
    if output_read_stream.is_some() {
        output_read_stream
            .expect("bufer does not exist")
            .flush()
            .unwrap();
    }
}

fn main() {
    let args = Args::parse();

    let output = args.output;
    let mut clipped = false;
    if args.clipped{
        clipped = true;
    }


    main_loop(
        output.clone(),
        args.gtf,
        args.input,
        args.overhang,
        args.flag_in,
        args.flag_out,
        args.mapq,
        args.output_write_read,
        clipped
    );

    if let Some(table) = args.table {
        let file = File::create_new(table.clone()).expect(&format!("output file {} should not exist.", &table));
        let mut stream = BufWriter::new(file);
        println!("table");
        let _ = stream.write("contig\tgene_name\ttranscript_name\texon_number\tstrand\tpos\texon_type\tspliced\tunspliced\tclipped\texon_intron\texon_other\tskipped\twrong_strand\te_isoform\n".as_bytes());
        file_to_table(output.clone(), &mut stream);
    }
}
//            skipped,
//wrong_strand, 
//e_isoform