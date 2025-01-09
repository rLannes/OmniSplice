#![allow(unused)]

use crate::common::point::InsideCounter;
use clap::builder::Str;
use clap::Parser;
use rust_htslib::bam::{Header, HeaderView, IndexedReader, Read, Reader, Record};
use std::collections::{hash_map, HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::io::BufRead;
use std::io::Write;
use std::io::{BufReader, BufWriter};
use std::process::{Command, Stdio};
use std::str;
use std::str::{from_utf8, FromStr};
use std::{thread, time};
use strand_specifier_lib::Strand;
use CigarParser::cigar::{Cigar, CigarOperation};
mod common;
use crate::common::point::{get_attr_id, Point, PointContainer};
use crate::common::utils::{read_toassign, ExonType, ReadAssign};
use itertools::Itertools;
use std::cmp::Reverse;
use std::fmt::{self, format};
use std::hash::Hash;
use std::path::{Path, PathBuf};
use strand_specifier_lib::{check_flag, LibType};

fn aln_bw(fa: &str, reference: &str, out_bam: &str) {
    let bowt_child = Command::new("bowtie2")
        .args(["-p", "2", "-x", reference, "-f", fa])
        .stdout(Stdio::piped())
        .spawn()
        .expect("bowtie2 command failed to start");

    let view_child = Command::new("samtools")
        .args(["view", "-bS", "-"])
        .stdin(Stdio::from(bowt_child.stdout.unwrap()))
        .stdout(Stdio::piped())
        .spawn()
        .expect("samtools view command failed to start");

    let mut sort_child = Command::new("samtools")
        .args(["sort", "-o", out_bam, "-"])
        .stdin(Stdio::from(view_child.stdout.unwrap()))
        .spawn()
        .expect("samtools sort command failed to start");
    let _result = sort_child.wait().unwrap();
    let index_child = Command::new("samtools")
        .args(["index", out_bam])
        .spawn()
        .expect("samtools index command failed to start");
}

fn get_sofclipped_seq(
    cig: &Cigar,
    exon_type: &ExonType,
    strand: &Strand,
    seq: &str,
) -> Option<String> {
    let l = seq.len();
    let mut r: Option<String> = Some("".to_string());
    match (exon_type, strand) {
        (ExonType::Donnor, Strand::Plus) | (ExonType::Acceptor, Strand::Minus) => {
            if let Some(n) = cig.get_soft_clipped_n(&Strand::Plus) {
                //println!("P {:?}", n);
                r = Some(seq[l - n as usize..l].to_string())
            }
        }
        (ExonType::Acceptor, Strand::Plus) | (ExonType::Donnor, Strand::Minus) => {
            if let Some(n) = cig.get_soft_clipped_n(&Strand::Minus) {
                //println!("M {:?} {:?}", n, cig);
                r = Some(seq[0..n as usize].to_string())
            }
        }
        (_, _) => unreachable!(),
    };
    r
}

#[derive(Debug, Clone)]
pub struct ReadInfo {
    gene: String,
    transcript: String,
    pos: i64,
    strand: String,
    exon_type: ExonType,
    chr_: String,
}

fn clipped_to_fasta(
    clipped_file: &str,
    clipped_fasta: &str,
    min_size: usize,
) -> HashMap<String, ReadInfo> {
    let mut seen: HashSet<String> = HashSet::new();
    let mut read_map: HashMap<String, ReadInfo> = HashMap::new();
    let f = File::open(clipped_file).expect("cannot open clipped file");
    let f_in = BufReader::new(f);
    let mut this_line: String;

    let file =
        File::create_new(clipped_fasta).expect("output clipped fasta file should not exist.");
    let mut f_out = BufWriter::new(file);
    let mut spt: Vec<&str>;

    let mut cig: Cigar;
    let mut exon_type: ExonType;
    let mut strand: Strand;
    let mut clip_seq: String = "".to_string();

    let mut id_: String = "".to_string();
    let mut l: usize = 0;

    for line in f_in.lines() {
        this_line = line.expect("cannot read line");
        spt = this_line.trim().split('\t').collect::<Vec<&str>>();
        id_ = spt[6].to_string();
        if seen.contains(&id_) {
            continue;
        } else {
            seen.insert(id_.clone());
        }

        cig = Cigar::from(spt[7]); //.expect(&format!("{}, invalid cigar", spt[1]));
        exon_type = ExonType::from(spt[5]);

        strand = Strand::from(spt[4]);

        clip_seq = get_sofclipped_seq(&cig, &exon_type, &strand, spt[8]).unwrap();

        l = clip_seq.len();
        if l < min_size {
            continue;
        }

        // FIXME ? overwritting?
        read_map.insert(
            id_.clone(),
            ReadInfo {
                gene: spt[2].to_string(),
                transcript: spt[3].to_string(),
                pos: spt[1].parse::<i64>().unwrap(),
                strand: spt[4].to_string(),
                exon_type: exon_type,
                chr_: spt[0].to_string(),
            },
        );

        let _ = f_out.write(format!(">{}\n{}\n", id_, clip_seq).as_bytes());
    }
    read_map
}

fn get_gtf_clipped(gtf: &str) -> Result<HashMap<String, PointContainer>, Box<dyn Error>> {
    let f = File::open(gtf).expect("cannot open gtf file");
    let reader = BufReader::new(f);

    let mut this_line: String; //::new();

    let mut chr_: String; // = "".to_string();
    let mut start: i64; //; = 0;
    let mut end: i64; // = 0;
    let mut strand: Strand;

    let mut gene_name: String; // = "".to_string();
    let mut transcript_id: String; // = "".to_string();

    // (gene_id, start, end)
    let mut been_seen: HashSet<(String, i64, ExonType)> = HashSet::new();
    let mut results: HashMap<String, PointContainer> = HashMap::new();

    for line in reader.lines() {
        this_line = line.expect("failed to read line");
        let spt = this_line.trim().split('\t').collect::<Vec<&str>>();
        if spt.len() < 8 {
            continue;
        }
        if spt[2] != "exon" {
            continue;
        }

        chr_ = spt[0].to_string();
        start = spt[3].parse::<i64>().unwrap();
        end = spt[4].parse::<i64>().unwrap();
        strand = Strand::from(spt[6]);

        if let Some((gene_tmp, tr_tmp)) =
            get_attr_id(spt[8], "gene_id").zip(get_attr_id(spt[8], "transcript_id"))
        {
            gene_name = gene_tmp;
            transcript_id = tr_tmp;
        } else if let Some((gene_tmp, tr_tmp)) =
            get_attr_id(spt[8], "gene_name").zip(get_attr_id(spt[8], "transcript_id"))
        {
            gene_name = gene_tmp;
            transcript_id = tr_tmp;
        } else {
            println!("warning gene id or transcript id not found: {:?}", spt);
            continue; // skip this gene
        }

        let start_type = match strand {
            Strand::Plus | Strand::NA => {
                start -= 1;
                ExonType::Acceptor
            }
            Strand::Minus => {
                start -= 1;
                ExonType::Donnor
            }
        };

        let end_type = match strand {
            Strand::Plus | Strand::NA => ExonType::Donnor,
            Strand::Minus => ExonType::Acceptor,
        };

        results
            .entry(gene_name.to_string().clone())
            .or_insert_with(|| PointContainer::new());
        if let Some(val) = results.get_mut(&gene_name.to_string()) {
            if !been_seen.contains(&(gene_name.to_string(), start, start_type)) {
                val.push(Point {
                    pos: start,
                    gene_name: gene_name.clone(), //get_gene_id(spt[8]),
                    transcript_id: transcript_id.clone(),
                    strand: strand,
                    counter: HashMap::new(),
                    exon_type: start_type,
                });
                been_seen.insert((gene_name.to_string(), start, start_type));
            }
            if !been_seen.contains(&(gene_name.to_string(), end, end_type)) {
                val.push(Point {
                    pos: end,
                    gene_name: gene_name.clone(), //get_gene_id(spt[8]),
                    transcript_id: transcript_id.clone(),
                    strand: strand,
                    counter: HashMap::new(),
                    exon_type: end_type,
                });
                been_seen.insert((gene_name.to_string(), end, end_type));
            }
        }
    }
    for (_k, v) in results.iter_mut() {
        v.sort()
    }
    Ok(results)
}

/* fn get_contig_map(header: &mut Header) -> HashMap<String, String>{

    let mut map: HashMap<String, String> = HashMap::new();
    for (key, records) in header.to_hashmap() {
        for record in records {
            if record.contains_key("SN"){
                map.insert(key, record["SN"].to_string());
                //seq.push(record["SN"].to_string());
            }
        }
    }
    map
} */

#[derive(Clone, Debug)]
pub struct BackSplicingCounter {
    gene: String,
    support_donnor: i64, // should I use hashset read name?
    support_acceptor: i64,
    tot_donnor: i64,
    tot_acceptor: i64,
    strand: Strand, //
    prime3: i64,
    prime5: i64,
    contig: String,
}

impl BackSplicingCounter {
    fn new(strand: Strand, prime5: i64, prime3: i64, contig: String, gene: String) -> Self {
        BackSplicingCounter {
            support_donnor: 0,
            support_acceptor: 0,
            tot_donnor: 0,
            tot_acceptor: 0,
            strand: strand,
            prime3: prime3,
            prime5: prime5,
            contig: contig,
            gene: gene,
        }
    }
    fn add_support_donnor(&mut self) -> () {
        self.support_donnor += 1;
    }
    fn add_support_acceptor(&mut self) -> () {
        self.support_acceptor += 1;
    }
    fn add_tot_acceptor(&mut self) -> () {
        self.tot_acceptor += 1;
    }
    fn add_tot_donnor(&mut self) -> () {
        self.tot_donnor += 1;
    }
    fn get_support(&self) -> i64 {
        self.support_donnor + self.support_acceptor
    }
    fn prettydump(&self) -> () {
        println!("{}", &format!("Position: {}:{}-{}({}), Total Support: {}, Donnor Support: {}, Acceptor Support: {}",
         self.contig, self.prime5, self.prime3, self.strand, self.support_acceptor + self.support_donnor, self.support_donnor, self.support_acceptor));
    }
}

impl fmt::Display for BackSplicingCounter {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}",
            &format!(
                "{}\t{}:{}-{}({})\t{}\t{}\t{}",
                self.gene,
                self.contig,
                self.prime5,
                self.prime3,
                self.strand,
                self.get_support(),
                self.support_donnor,
                self.support_acceptor
            )
        )
    }
}

// What do we want?
// circRNA
// When do we want it?
// now!
// TODO refactor this mess!
// Defined a clipped read struct?
pub fn parse_bam(
    bam_file: &str,
    transcript_map: HashMap<String, PointContainer>,
    readmap: &HashMap<String, ReadInfo>,
) -> HashMap<(String, i64, i64), BackSplicingCounter> {
    let mut counter: i64;
    let mut record: Record;
    let mut pos_s: i64;
    let mut pos_e: i64;
    let mut cig: Cigar;
    let mut flag: u16;
    let mut contig: String;
    let mut read_name: String;
    let mut original_mapping: ReadInfo;

    let mut results: HashMap<(String, i64, i64), BackSplicingCounter> = HashMap::new();

    let mut bam = Reader::from_path(bam_file).unwrap();
    let header = Header::from_template(bam.header());
    let header_view: HeaderView = HeaderView::from_header(&header);
    let mut cpt = 0;
    let mut back_flag = false;

    for r in bam.records() {
        back_flag = false;
        record = r.unwrap();
        pos_s = record.pos();
        cig = Cigar::from_str(&record.cigar().to_string()).unwrap();
        pos_e = cig.get_end_of_aln(&pos_s);
        flag = record.flags();

        if !check_flag(flag, 0, 256) {
            continue;
        }
        if record.tid() < 0 {
            continue;
        }
        if record.mapq() < 13 {
            continue;
        }
        contig = from_utf8(header_view.tid2name(record.tid() as u32))
            .expect("failed to parse contig name")
            .to_string();
        read_name = String::from_utf8(record.qname().to_vec()).expect("cannot parse read name");

        //read_info = readmap.get(&read_name).unwrap().clone();
        if let Some(original_map_info) = readmap.get(&read_name) {
            println!("{:?} {:?}", original_map_info, read_name);
            let transcript_container = transcript_map.get(&original_map_info.gene).unwrap().clone();
            for (indice, tr_junction) in transcript_container.iter().enumerate() {
                if (tr_junction.exon_type == original_map_info.exon_type) | (contig != original_map_info.chr_){
                    continue;
                }
                 else if (tr_junction.pos == pos_s) & (pos_s < original_map_info.pos) {
                    if (tr_junction.strand == Strand::Plus)
                        & (tr_junction.exon_type == ExonType::Acceptor)
                    {
                        results
                            .entry((contig.clone(), pos_s, original_map_info.pos))
                            .or_insert_with(|| {
                                BackSplicingCounter::new(
                                    tr_junction.strand,
                                    pos_s,
                                    original_map_info.pos,
                                    contig.clone(),
                                    original_map_info.gene.clone(),
                                )
                            })
                            .add_support_acceptor();
                    } else if (tr_junction.strand == Strand::Minus)
                        & (tr_junction.exon_type == ExonType::Donnor)
                    {
                        results
                            .entry((contig.clone(), pos_s, original_map_info.pos))
                            .or_insert_with(|| {
                                BackSplicingCounter::new(
                                    tr_junction.strand,
                                    pos_s,
                                    original_map_info.pos,
                                    contig.clone(),
                                    original_map_info.gene.clone(),
                                )
                            })
                            .add_support_donnor();
                    } else {
                        println!("{} {}", tr_junction.strand, tr_junction.exon_type);
                        continue;
                    }
                } else if (tr_junction.pos == pos_e) & (pos_s > original_map_info.pos) {
                    if (tr_junction.strand == Strand::Plus)
                        & (tr_junction.exon_type == ExonType::Donnor)
                    {
                        results
                            .entry((contig.clone(), original_map_info.pos, pos_e))
                            .or_insert_with(|| {
                                BackSplicingCounter::new(
                                    tr_junction.strand,
                                    original_map_info.pos,
                                    pos_e,
                                    contig.clone(),
                                    original_map_info.gene.clone(),
                                )
                            })
                            .add_support_donnor();
                    } else if (tr_junction.strand == Strand::Minus)
                        & (tr_junction.exon_type == ExonType::Acceptor)
                    {
                        results
                            .entry((contig.clone(), original_map_info.pos, pos_e))
                            .or_insert_with(|| {
                                BackSplicingCounter::new(
                                    tr_junction.strand,
                                    original_map_info.pos,
                                    pos_e,
                                    contig.clone(),
                                    original_map_info.gene.clone(),
                                )
                            })
                            .add_support_acceptor();
                    } else {
                        println!("{} {}", tr_junction.strand, tr_junction.exon_type);
                        continue;
                    }
                }
                /*                 else if (original_map_info.strand == "-") & ((tr_junction.pos == pos_e) | (tr_junction.pos == pos_s)){
                    println!("{:?}", transcript_container);
                    println!("{:?} {:?}", tr_junction, original_map_info);
                    println!("{:?} {:?} {:?}\n", pos_s, pos_e, read_name);

                 //   continue;
                } */
            }
        }
    }
    results
}

// TODO add specific subcommand to retrieve only specific reads;
/// OmniSplice Backsplicing
/// a Software to identify backslicing read.
/// you need bowtie2 in the path to use this software.
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Input file must be clipped read results from omnisplice
    #[arg(short, long, required = true)]
    input_clipped_read: String,
    /// Name of the base for output file will generate a .fna .bam .bai .backspliced.tsv
    #[arg(short, long, required = true)]
    output_base: String,
    // path to a bowtie2 reference for yor genome
    #[arg(short, long, required = true)]
    bowtie2_ref: String,
    /*     // path to a bowtie2 reference for yor genome
    #[arg(short, long, required = true)]
    bowtie2_ref:String, */
    /// Name of GTF Input file define the feature to look at
    /// (v1) only consider feature annotated as exon
    /// you may want to subset gene / feature you are interested in.
    #[arg(short, long, required = true)]
    gtf: String,
    /// clipped part of the read are only consider if thei length is longet than this value. defualt 20
    #[arg(short, long, default_value_t = 20)]
    min_clipped_size: usize,
}

fn main() {
    let args = Args::parse();

    let output_base = args.output_base;

    let mut clipped_fasta = PathBuf::from(&output_base);
    clipped_fasta.set_extension("fna");
    let mut bw_bam = PathBuf::from(&output_base);
    bw_bam.set_extension("bam");
    let mut output_file = PathBuf::from(&output_base);
    output_file.set_extension("backspliced.tsv");
    let gtf = args.gtf;
    let clipped_size_min = args.min_clipped_size;
    let bw2_ref = args.bowtie2_ref;
    let clipped_file = args.input_clipped_read;

    //let clipped_file = "/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/Jackie_MD6_OmniSplice/omniRun/trimmed_L459_1501_S2_L001_Aligned.read_through.read.tsv";
    //let clipped_fasta = "/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/Jackie_MD6_OmniSplice/trimmed_L459_1501_S2_L001_Aligned.clipped.fna";
    //let bw2_ref = "/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/Jackie_MD6_OmniSplice/bwRef/dm6_bw2";
    //let bw_bam = "/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/Jackie_MD6_OmniSplice/trimmed_L459_1501_S2_L001_Aligned.clipped.bam";
    //let gtf = "/lab/solexa_yamashita/people/Romain/References/MD6/gtf/dmel-all-r6.48.gtf";
    //let clipped_size_min = 20;
    //let output_file = "/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/Jackie_MD6_OmniSplice/Backsplicing/trimmed_L459_1501_S2_L001.backspliced.tsv";

    let map_read = clipped_to_fasta(
        &clipped_file,
        clipped_fasta.as_path().to_str().unwrap(),
        clipped_size_min,
    ); // read map -> HashMap<String, ReadInfo>
    aln_bw(
        clipped_fasta.as_path().to_str().unwrap(),
        &bw2_ref,
        bw_bam.as_path().to_str().unwrap(),
    );
    let mut point_cont = get_gtf_clipped(&gtf).expect("Failed to parse GTF"); // transcript is -> point container not optimal but I reuse what exist!

    let file = File::create_new(output_file).expect("output clipped fasta file should not exist.");
    let mut f_out = BufWriter::new(file);

    let result = parse_bam(bw_bam.as_path().to_str().unwrap(), point_cont, &map_read);
    for (key, val) in result.iter().sorted_by_key(|x| Reverse(x.1.get_support())) {
        f_out.write_all(format!("{val}\n").as_bytes());
    }

    f_out.flush();
}
