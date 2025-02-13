#![allow(irrefutable_let_patterns)]
use crate::common::point::{get_attr_id, InsideCounter};
use crate::common::utils::{read_toassign, ExonType, ReadAssign, ReadToWriteHandle};
use bio::bio_types::annot::contig;
use bio::data_structures::interval_tree;
use bio::utils::Interval;
use rust_htslib::bam::record::Record;
use rust_htslib::bam::{IndexedReader, Read};
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::{File, OpenOptions};
use std::io::BufRead;
use std::io::Write;
use std::io::{BufReader, BufWriter};
use std::process::{Command, Stdio};
use std::str::FromStr;
use strand_specifier_lib::Strand;
use strand_specifier_lib::{check_flag, LibType};
use CigarParser::cigar::Cigar;

#[derive(Debug)]
pub struct TreeData {
    gene_name: String,
    transcript: Vec<String>,
    start: i64,
    end: i64,
    strand: Strand,
    contig: String,
    start_type: ExonType,
    end_type: ExonType,
    counter_start: HashMap<ReadAssign, i32>,
    counter_end: HashMap<ReadAssign, i32>,
}

impl TreeData {
    fn new(
        gene_name: String,
        transcript: String,
        start: i64,
        end: i64,
        contig: String,
        strand: Strand,
    ) -> Self {
        let mut tr = Vec::new();
        tr.push(transcript);
        let start_type = match strand {
            Strand::Minus => ExonType::Donnor,
            Strand::Plus | Strand::NA => ExonType::Acceptor,
        };
        let end_type = match strand {
            Strand::Minus => ExonType::Acceptor,
            Strand::Plus | Strand::NA => ExonType::Donnor,
        };
        TreeData {
            gene_name: gene_name,
            transcript: tr,
            start: start,
            end: end,
            contig: contig,
            strand: strand,
            start_type: start_type,
            end_type: end_type,
            counter_start: HashMap::new(),
            counter_end: HashMap::new(),
        }
    }

    fn update_start_counter(&mut self, item: ReadAssign) -> () {
        if let Some(val) = self.counter_start.get_mut(&item) {
            *val += 1;
        } else {
            self.counter_start.insert(item, 1);
        }
    }

    fn update_end_counter(&mut self, item: ReadAssign) -> () {
        if let Some(val) = self.counter_end.get_mut(&item) {
            *val += 1;
        } else {
            self.counter_end.insert(item, 1);
        }
    }

    fn parse_read(
        &mut self,
        aln_start: i64,
        aln_end: i64,
        cigar: &Cigar,
        read_strand: &Strand,
        overhang: i64,
        out_file_read_buffer: &mut ReadToWriteHandle,
        read_name: Option<String>,
        sequence: Option<String>,
    ) -> () {
        let seq = match sequence {
            Some(seq) => seq,
            _ => "NoSeq".to_string()
        };
        let r_name = match read_name {
            Some(seq) => seq,
            _ => "NoSeq".to_string()
        };

        if let Some(start_map) = read_toassign(
            self.strand,
            self.start,
            self.start_type,
            aln_start,
            aln_end,
            cigar,
            read_strand,
            overhang,
        ) {
            self.update_start_counter(start_map);


            // I am not a fna of this it is very verbose but it seems to do the job.
            // also I am going to test that match case millions of time... but it is setup at start..
            match start_map {
                ReadAssign::Empty =>  match  &mut out_file_read_buffer.empty {
                        Some(handle) => handle.write(self.dump_reads_seq(&seq, &r_name, false, cigar).as_bytes())
                        .unwrap_or_else(|_| panic!("cannot write reads")),
                        _ => 0
                },
                ReadAssign::ReadThrough => match &mut out_file_read_buffer.read_through {
                    Some(handle) => handle.write(self.dump_reads_seq(&seq, &r_name, false, cigar).as_bytes())
                    .unwrap_or_else(|_| panic!("cannot write reads")),
                    _ => 0
                },
                ReadAssign::ReadJunction(_, _) => match &mut out_file_read_buffer.read_junction {
                    Some(handle) => handle.write(self.dump_reads_seq(&seq, &r_name, false, cigar).as_bytes())
                    .unwrap_or_else(|_| panic!("cannot write reads")),
                    _ => 0
                },
                ReadAssign::Unexpected => match &mut out_file_read_buffer.unexpected {
                    Some(handle) => handle.write(self.dump_reads_seq(&seq, &r_name, false, cigar).as_bytes())
                    .unwrap_or_else(|_| panic!("cannot write reads")),
                    _ => 0
                },
                ReadAssign::FailPosFilter => match &mut out_file_read_buffer.fail_pos_filter {
                    Some(handle) => handle.write(self.dump_reads_seq(&seq, &r_name, false, cigar).as_bytes())
                    .unwrap_or_else(|_| panic!("cannot write reads")),
                    _ => 0
                },
                ReadAssign::WrongStrand => match &mut out_file_read_buffer.wrong_strand {
                    Some(handle) => handle.write(self.dump_reads_seq(&seq, &r_name, false, cigar).as_bytes())
                    .unwrap_or_else(|_| panic!("cannot write reads")),
                    _ => 0
                },
                ReadAssign::FailQc => match &mut out_file_read_buffer.fail_qc {
                    Some(handle) => handle.write(self.dump_reads_seq(&seq, &r_name, false, cigar).as_bytes())
                    .unwrap_or_else(|_| panic!("cannot write reads")),
                    _ => 0
                },
                ReadAssign::EmptyPileup => match &mut out_file_read_buffer.empty_pileup {
                    Some(handle) => handle.write(self.dump_reads_seq(&seq, &r_name, false, cigar).as_bytes())
                    .unwrap_or_else(|_| panic!("cannot write reads")),
                    _ => 0
                },
                ReadAssign::Skipped(_, _) => match &mut out_file_read_buffer.skipped {
                    Some(handle) => handle.write(self.dump_reads_seq(&seq, &r_name, false, cigar).as_bytes())
                    .unwrap_or_else(|_| panic!("cannot write reads")),
                    _ => 0
                },
                ReadAssign::SoftClipped => match &mut out_file_read_buffer.soft_clipped {
                    Some(handle) => handle.write(self.dump_reads_seq(&seq, &r_name, false, cigar).as_bytes())
                    .unwrap_or_else(|_| panic!("cannot write reads")),
                    _ => 0
                },
                ReadAssign::OverhangFail => match &mut out_file_read_buffer.overhang_fail {
                    Some(handle) => handle.write(self.dump_reads_seq(&seq, &r_name, false, cigar).as_bytes())
                    .unwrap_or_else(|_| panic!("cannot write reads")),
                    _ => 0
                },
                
            };
            /*if let Some(handle) = out_file_read_buffer {
                match (clipped, start_map, &read_name, &sequence) {
                    (true, ReadAssign::SoftClipped, Some(r_name), Some(seq)) => {
                        handle.write(self.dump_reads_seq(seq, r_name, false, cigar).as_bytes());
                    }
                    (false, _, Some(r_name), Some(seq)) => {
                        handle.write(self.dump_reads_seq(seq, r_name, false, cigar).as_bytes());
                    }
                    (_, _, _, _) => (),
                }
            }*/
        }

        if let Some(end_map) = read_toassign(
            self.strand,
            self.end,
            self.end_type,
            aln_start,
            aln_end,
            cigar,
            read_strand,
            overhang,
        ) {
            self.update_end_counter(end_map);
            
            match end_map {
                ReadAssign::Empty =>  match &mut out_file_read_buffer.empty {
                        Some(handle) => handle.write(self.dump_reads_seq(&seq, &r_name,
                             true, cigar).as_bytes())
                        .unwrap_or_else(|_| panic!("cannot write reads")),
                        _ => 0
                },
                ReadAssign::ReadThrough => match &mut out_file_read_buffer.read_through {
                    Some(handle) => handle.write(self.dump_reads_seq(&seq, &r_name, true, cigar).as_bytes())
                    .unwrap_or_else(|_| panic!("cannot write reads")),
                    _ => 0
                },
                ReadAssign::ReadJunction(_, _) => match &mut out_file_read_buffer.read_junction {
                    Some(handle) => handle.write(self.dump_reads_seq(&seq, &r_name, true, cigar).as_bytes())
                    .unwrap_or_else(|_| panic!("cannot write reads")),
                    _ => 0
                },
                ReadAssign::Unexpected => match &mut out_file_read_buffer.unexpected {
                    Some(handle) => handle.write(self.dump_reads_seq(&seq, &r_name, true, cigar).as_bytes())
                    .unwrap_or_else(|_| panic!("cannot write reads")),
                    _ => 0
                },
                ReadAssign::FailPosFilter => match &mut out_file_read_buffer.fail_pos_filter {
                    Some(handle) => handle.write(self.dump_reads_seq(&seq, &r_name, true, cigar).as_bytes())
                    .unwrap_or_else(|_| panic!("cannot write reads")),
                    _ => 0
                },
                ReadAssign::WrongStrand => match &mut  out_file_read_buffer.wrong_strand {
                    Some(handle) => handle.write(self.dump_reads_seq(&seq, &r_name, true, cigar).as_bytes())
                    .unwrap_or_else(|_| panic!("cannot write reads")),
                    _ => 0
                },
                ReadAssign::FailQc => match &mut out_file_read_buffer.fail_qc {
                    Some(handle) => handle.write(self.dump_reads_seq(&seq, &r_name, true, cigar).as_bytes())
                    .unwrap_or_else(|_| panic!("cannot write reads")),
                    _ => 0
                },
                ReadAssign::EmptyPileup => match &mut out_file_read_buffer.empty_pileup {
                    Some(handle) => handle.write(self.dump_reads_seq(&seq, &r_name, true, cigar).as_bytes())
                    .unwrap_or_else(|_| panic!("cannot write reads")),
                    _ => 0
                },
                ReadAssign::Skipped(_, _) => match &mut out_file_read_buffer.skipped {
                    Some(handle) => handle.write(self.dump_reads_seq(&seq, &r_name, true, cigar).as_bytes())
                    .unwrap_or_else(|_| panic!("cannot write reads")),
                    _ => 0
                },
                ReadAssign::SoftClipped => match &mut out_file_read_buffer.soft_clipped {
                    Some(handle) => handle.write(self.dump_reads_seq(&seq, &r_name, true, cigar).as_bytes())
                    .unwrap_or_else(|_| panic!("cannot write reads")),
                    _ => 0
                },
                ReadAssign::OverhangFail => match &mut out_file_read_buffer.overhang_fail {
                    Some(handle) => handle.write(self.dump_reads_seq(&seq, &r_name, true, cigar).as_bytes())
                    .unwrap_or_else(|_| panic!("cannot write reads")),
                    _ => 0
                },
                
            };
            /*if let Some(handle) = out_file_read_buffer {
                match (clipped, end_map, &read_name, &sequence) {
                    (true, ReadAssign::SoftClipped, Some(r_name), Some(seq)) => {
                        handle.write(self.dump_reads_seq(seq, r_name, true, cigar).as_bytes());
                    }
                    (false, _, Some(r_name), Some(seq)) => {
                        handle.write(self.dump_reads_seq(seq, r_name, true, cigar).as_bytes());
                    }
                    (_, _, _, _) => (),
                }
            }
        }*/

        ()
        }
    }
    fn dump_base(&self, end: bool) -> Vec<String> {
        let mut res = Vec::new();
        let mut exontype = ExonType::Acceptor;
        let mut pos = 0;
        if end {
            exontype = self.end_type;
            pos = self.end;
        } else {
            exontype = self.start_type;
            pos = self.start;
        }
        for tr in &self.transcript {
            res.push(format!(
                "{}\t{}\t{}\t{}\t{}\t{}",
                self.contig, pos, self.gene_name, tr, self.strand, exontype
            ));
        }
        res
    }

    fn dump_counter(&self, end: bool) -> String {
        let base_vec = self.dump_base(end);
        let mut results = Vec::new();

        if end {
            if self.counter_end.is_empty() {
                for base in &base_vec {
                    results.push(format!("{}\tempty\t0", base))
                }
            } else {
                for (assign, value) in &self.counter_end {
                    for base in &base_vec {
                        results.push(format!("{}\t{}\t{}", base, assign, value))
                        //base.push_str(&format!("\t{}\t{}", assign, value))
                    }
                }
            }
        } else {
            if self.counter_start.is_empty() {
                for base in &base_vec {
                    results.push(format!("{}\tempty\t0", base))
                }
            } else {
                for (assign, value) in &self.counter_start {
                    for base in &base_vec {
                        results.push(format!("{}\t{}\t{}", base, assign, value))
                        //base.push_str(&format!("\t{}\t{}", assign, value))
                    }
                }
            }
        }

        results.join("\n")
    }

    fn dump_reads_seq(&self, sequence: &String, seqname: &String, end: bool, cigar: &Cigar) -> String {
        ///
        /// :param end: if set to true look at the end of the sequence (end > start)
        let mut base_vec = self.dump_base(end);

        for e in &mut base_vec {
            e.push_str(&format!("\t{}\t{}\t{}", seqname, cigar, sequence))
        }

        return format!("{}\n", base_vec.join("\n"));
    }
}

pub fn dump_tree_to_cat_results(
    hash_tree: &HashMap<String, interval_tree::IntervalTree<i64, TreeData>>,
    out_file: &str,
) -> () {
    let presorted = format!("{}.presorted", out_file);
    let header = "Contig\tPosition\tGeneID\tTranscriptID\tstrand\tExonEndType\tCategory\tReadCount\n".as_bytes();
    {
        let file = File::create_new(presorted.clone()).expect("output file should not exist.");
        let mut stream = BufWriter::new(file);

        //stream.write(header);

        for (contig, subtree) in hash_tree {
            // get ALL entry
            for node in subtree.find(i64::MIN..i64::MAX) {
                stream.write(node.data().dump_counter(false).as_bytes());
                stream.write("\n".as_bytes());
                stream.write(node.data().dump_counter(true).as_bytes());
                stream.write("\n".as_bytes());
            }
        }
        stream.flush().unwrap();
    }

    let file = File::create_new(out_file).expect("output file should not exist.");

    let mut file = OpenOptions::new()
        .write(true)
        .append(true)
        .open(out_file).expect("could not open file");
    file.write(header);
    //let mut stream = BufWriter::new(file);

    let sort = Command::new("sort")
        .args([
            "-k1,1",
            "-k3,3",
            "-k4,4",
            "-k2,2n",
            "-k7,7",
            presorted.as_str(),
        ])
        .stdout(file)
        .spawn()
        .expect("sort command failed ")
        .wait()
        .expect("sort command failed ");
    //stream.flush().unwrap();

    let rm_out = Command::new("rm")
        .args(["-f", presorted.as_str()])
        .output()
        .expect("failed to remove presorted");
}


/// parse the gtf and return a hashmap<chromosome> -> intervalTree(exon(start, end), associated_data(gene_name...))
pub fn gtf_to_tree(
    file: &str,
) -> Result<HashMap<String, interval_tree::IntervalTree<i64, TreeData>>, Box<dyn Error>> {
    let f = File::open(file)?;
    let reader = BufReader::new(f);
    let mut this_line: String; //::new();

    let mut chr_: String; // = "".to_string();
    let mut start: i64; //; = 0;
    let mut end: i64; // = 0;
    let mut strand: Strand;

    let mut gene_name: String; // = "".to_string();
    let mut transcript_id: String; // = "".to_string();
    let mut flag_exon_already_found = false;

    // (gene_id, start, end)
    let mut been_seen: HashSet<(String, i64, ExonType)> = HashSet::new();
    let mut results: HashMap<String, interval_tree::IntervalTree<i64, TreeData>> = HashMap::new();
    let mut spt: Vec<&str> = Vec::new();

    for line in reader.lines() {
        this_line = line?;
        spt = this_line.trim().split('\t').collect::<Vec<&str>>();
        if spt.len() < 8 {
            continue;
        }
        if spt[2] != "exon" {
            continue;
        }

        chr_ = spt[0].to_string();
        start = spt[3].parse::<i64>()? - 1;
        end = spt[4].parse::<i64>()?;
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
            continue;
        }

        results
            .entry(chr_.clone())
            .or_insert_with(|| interval_tree::IntervalTree::new());
        flag_exon_already_found = false;

        if let Some(this_tree) = results.get_mut(&chr_) {
            for (ref mut node) in this_tree.find_mut(start..end) {
                if (node.interval().start == start) && (node.interval().end == end) {
                    if let n = node.data() {
                        if n.gene_name == gene_name {
                            n.transcript.push(transcript_id.clone());
                            flag_exon_already_found = true;
                        }
                    }
                }
            }
            if !flag_exon_already_found {
                this_tree.insert(
                    start..end,
                    TreeData::new(gene_name, transcript_id, start, end, chr_, strand),
                );
            }
        }
    }
    Ok(results)
}

pub fn update_tree_with_bamfile(
    hash_tree: &mut HashMap<String, interval_tree::IntervalTree<i64, TreeData>>,
    bam_file: &str,
    // required parameter
    library_type: LibType,
    overhang: i64,

    // TODO not yet implemented QC optional option (None for no filter) // TODO
    // default 0
    flag_in: u16,
    // default 1024(DUPLICATE) + 256 (NOT_PRIMARY_ALN) + 2048(SUPPLEMENTARY) + 512(FAIL_QC)
    flag_out: u16,
    // mapq (must be) >=   default 13
    mapq: u8,
    out_file_read_buffer: &mut ReadToWriteHandle
) -> Result<(), Box<dyn Error>> {
    let mut read_name: Option<String> = None;
    let mut seq: Option<Vec<u8>> = None;
    let mut sequence: Option<String> = None;

    let mut counter: i64;
    let mut record: Record;
    let mut pos_s: i64;
    let mut pos_e: i64;
    let mut cig: Cigar;
    let mut flag: u16;
    let mut bam = IndexedReader::from_path(bam_file).unwrap();

    for (contig, subtree) in hash_tree.iter_mut() {
        //counter = 0;
        //println!("Contig: {}", contig);

        match bam.fetch(&contig) {
            Ok(_) => (),
            Err(_) => {
                println!("WARNING {} not found in the bam file check", contig);
                continue;
            }
        }
        for r in bam.records() {
            //counter += 1;
            //if counter % 1_000_000 == 0 {
            //    println!("Contig: {}; {} reads done", contig, counter);
            //}
            record = r.unwrap();
            pos_s = record.pos();
            cig = Cigar::from_str(&record.cigar().to_string()).unwrap();
            pos_e = cig.get_end_of_aln(&pos_s);
            flag = record.flags();

            // QC
            if (!check_flag(flag, flag_in, flag_out)) || (record.mapq() < mapq) {
                continue;
            }

            read_name = Some(
                String::from_utf8(record.qname().to_vec()).expect("cannot parse read name"),
            );
            let seq = record.seq().as_bytes();
            sequence = Some(String::from_utf8(seq).expect("cannot parse sequence"));


            if let Some(read_strand) = library_type.get_strand(flag) {
                for (ref mut exon) in subtree.find_mut(pos_s..pos_e + 1) {
                    exon.data().parse_read(
                        pos_s,
                        pos_e,
                        &cig,
                        &read_strand,
                        overhang,
                        out_file_read_buffer,
                        read_name.clone(),
                        sequence.clone(),
                    );
                }
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests_it {
    use super::*;

    #[test]
    fn test_1() {
        let mut tree = interval_tree::IntervalTree::new();
        tree.insert(11..20, "Range_1");
        tree.insert(11..20, "Range_3");
        tree.insert(25..30, "Range_2");
        for r in tree.find(15..25) {
            println!("{}", r.data());
        }
    }

    #[test]
    fn test_2() {
        let mut tree = interval_tree::IntervalTree::new();
        tree.insert(11..20, "Range_1");
        tree.insert(11..20, "Range_3");
        tree.insert(25..30, "Range_2");
        for r in tree.find(15..25 + 1) {
            println!("{}", r.data());
        }
    }
    /*#[test]
    fn test_2() {
    let f = "/lab/solexa_yamashita/people/Romain/Projets/Adrienne/OmniSplice/X_auto/test_gtf_chr.gtf".to_string();
    let mut tree = gtf_to_tree(f).unwrap();

    for (chr_, subtree) in tree.iter(){
        for inter in subtree.find(0..100_000_000){
            println!("{:?}", inter);
        }
        }
    }*/
}
