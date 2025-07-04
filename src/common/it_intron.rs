#![allow(irrefutable_let_patterns)]
//use crate::common::point::{get_attr_id, InsideCounter};
use crate::common::utils::{
    read_toassign, Exon, ExonType, ReadAssign, ReadToWriteHandle, SplicingEvent,
};
use bio::alignment::sparse::HashMapFx;
use bio::bio_types::annot::contig;
use bio::data_structures::interval_tree;
use bio::seq_analysis;
use bio::utils::Interval;
use itertools::Itertools;
use rust_htslib::bam::record::Record;
use rust_htslib::bam::{header, IndexedReader, Read};
use std::collections::{hash_map, HashMap, HashSet};
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
/// This structure is used in the intervall tree HASHMAP(contig) -> IT(intron(start, end)) -> TreeDataIntron
/// Represent an intron or the start of end of the transcript.
/// start < end (genome position + strand).
/// for this reason start and end are optional but at least one of those must be set
/// start_type and end_type  are deduced from the strand and positions, similarly from start and end they are also optional with one needed to be setup
/// the three counters are the things actually counting the reads.
pub struct TreeDataIntron {
    gene_name: String,
    transcript_intron: Vec<(String, i32)>,
    start: Option<i64>,
    end: Option<i64>,
    strand: Strand,
    contig: String,
    start_type: Option<ExonType>,
    end_type: Option<ExonType>,
    counter_start: HashMap<ReadAssign, i32>,
    counter_end: HashMap<ReadAssign, i32>,
    counter_intron: HashMap<SplicingEvent, i32>,
}

impl TreeDataIntron {
    /// if a record overlap this TreeDataIntron,
    /// parse the record and update the TreeDataIntron accordingly
    fn update_from_read(&mut self, record: &Record) {}

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
        valid_junction: &HashMap<(String, i64, i64), Strand>,
        valid_j_gene: &HashMap<String, HashSet<(i64, i64)>>,
    ) -> () {
        //println!("{:?}", self);

        let seq = match sequence {
            Some(seq) => seq,
            _ => "NoSeq".to_string(),
        };
        let r_name = match read_name {
            Some(seq) => seq,
            _ => "NoSeq".to_string(),
        };

        let start_map = read_toassign(
            self.strand,
            self.start,
            self.start_type,
            aln_start,
            aln_end,
            cigar,
            read_strand,
            overhang,
        );
        self.write_to_read_file(start_map, out_file_read_buffer, false, &seq, &r_name, cigar);
        TreeDataIntron::update_counter(&mut self.counter_start, start_map);

        let end_map = read_toassign(
            self.strand,
            self.end,
            self.end_type,
            aln_start,
            aln_end,
            cigar,
            read_strand,
            overhang,
        );
        //panic!("");
        TreeDataIntron::update_counter(&mut self.counter_end, end_map);
        self.write_to_read_file(end_map, out_file_read_buffer, true, &seq, &r_name, cigar);

        // valid_j_gene.get(&self.gene_id).unwrap_or(HashSet::new()),

        match (
            SplicingEvent::from_read_assign(
                start_map,
                self.contig.clone(),
                //valid_junction,
                valid_j_gene.get(&self.gene_name).unwrap_or(&HashSet::new()),
                self.start,
                self.end,
                &self.strand,
                read_strand,
            ),
            SplicingEvent::from_read_assign(
                end_map,
                self.contig.clone(),
                //valid_junction,
                valid_j_gene.get(&self.gene_name).unwrap_or(&HashSet::new()),
                self.start,
                self.end,
                &self.strand,
                read_strand,
            ),
        ) {
            (Some(x), None) => {
                TreeDataIntron::update_splicing_counter(&mut self.counter_intron, Some(x))
            }
            (None, Some(x)) => {
                TreeDataIntron::update_splicing_counter(&mut self.counter_intron, Some(x))
            }
            (Some(SplicingEvent::Spliced), Some(SplicingEvent::Spliced)) => {
                TreeDataIntron::update_splicing_counter(
                    &mut self.counter_intron,
                    Some(SplicingEvent::Spliced),
                )
            }
            (Some(left), Some(right)) => TreeDataIntron::update_splicing_counter(
                &mut self.counter_intron,
                SplicingEvent::merge(Some(left), Some(right)),
            ),
            (None, None) => (),
            (_, _) => (),
        }
    }
    fn update_counter(counter: &mut HashMap<ReadAssign, i32>, item: Option<ReadAssign>) -> () {
        match item {
            Some(read_assign) => {
                if let Some(val) = counter.get_mut(&read_assign) {
                    *val += 1;
                } else {
                    counter.insert(read_assign, 1);
                }
            }
            _ => (),
        }
    }

    fn update_splicing_counter(
        counter: &mut HashMap<SplicingEvent, i32>,
        item: Option<SplicingEvent>,
    ) -> () {
        match item {
            Some(read_assign) => {
                if let Some(val) = counter.get_mut(&read_assign) {
                    *val += 1;
                } else {
                    counter.insert(read_assign, 1);
                }
            }
            _ => (),
        }
    }

    fn dump_base(&self, end: bool) -> Vec<String> {
        let mut res = Vec::new();
        let mut exontype = ExonType::Acceptor;
        let mut pos = 0;
        if end {
            exontype = self.end_type.unwrap();
            pos = self.end.unwrap();
        } else {
            exontype = self.start_type.unwrap();
            pos = self.start.unwrap();
        }
        for (tr, intron) in &self.transcript_intron {
            res.push(format!(
                "{}\t{}\t{}\t{}\t{}\t{}",
                self.contig, pos, self.gene_name, tr, self.strand, exontype
            ));
        }
        res
    }

    fn dump_reads_seq(&self, sequence: &str, seqname: &str, end: bool, cigar: &Cigar) -> String {
        ///
        /// :param end: if set to true look at the end of the sequence (end > start)
        let mut base_vec = self.dump_base(end);

        for e in &mut base_vec {
            e.push_str(&format!("\t{}\t{}\t{}", seqname, cigar, sequence))
        }

        return format!("{}\n", base_vec.join("\n"));
    }

    fn write_to_read_file(
        &self,
        read_assign: Option<ReadAssign>,
        out_file_read_buffer: &mut ReadToWriteHandle,
        end: bool,
        seq: &str,
        r_name: &str,
        cigar: &Cigar,
    ) {
        match read_assign {
            None => 0,
            Some(ReadAssign::Empty) => match &mut out_file_read_buffer.empty {
                Some(handle) => handle
                    .write(self.dump_reads_seq(seq, r_name, end, cigar).as_bytes())
                    .unwrap_or_else(|_| panic!("cannot write reads")),
                _ => 0,
            },
            Some(ReadAssign::ReadThrough) => match &mut out_file_read_buffer.read_through {
                Some(handle) => handle
                    .write(self.dump_reads_seq(seq, r_name, end, cigar).as_bytes())
                    .unwrap_or_else(|_| panic!("cannot write reads")),
                _ => 0,
            },
            Some(ReadAssign::ReadJunction(_, _)) => match &mut out_file_read_buffer.read_junction {
                Some(handle) => handle
                    .write(self.dump_reads_seq(seq, r_name, end, cigar).as_bytes())
                    .unwrap_or_else(|_| panic!("cannot write reads")),
                _ => 0,
            },
            Some(ReadAssign::Unexpected) => match &mut out_file_read_buffer.unexpected {
                Some(handle) => handle
                    .write(self.dump_reads_seq(seq, r_name, end, cigar).as_bytes())
                    .unwrap_or_else(|_| panic!("cannot write reads")),
                _ => 0,
            },
            Some(ReadAssign::FailPosFilter) => match &mut out_file_read_buffer.fail_pos_filter {
                Some(handle) => handle
                    .write(self.dump_reads_seq(seq, r_name, end, cigar).as_bytes())
                    .unwrap_or_else(|_| panic!("cannot write reads")),
                _ => 0,
            },
            Some(ReadAssign::WrongStrand) => match &mut out_file_read_buffer.wrong_strand {
                Some(handle) => handle
                    .write(self.dump_reads_seq(seq, r_name, end, cigar).as_bytes())
                    .unwrap_or_else(|_| panic!("cannot write reads")),
                _ => 0,
            },
            Some(ReadAssign::FailQc) => match &mut out_file_read_buffer.fail_qc {
                Some(handle) => handle
                    .write(self.dump_reads_seq(seq, r_name, end, cigar).as_bytes())
                    .unwrap_or_else(|_| panic!("cannot write reads")),
                _ => 0,
            },
            Some(ReadAssign::EmptyPileup) => match &mut out_file_read_buffer.empty_pileup {
                Some(handle) => handle
                    .write(self.dump_reads_seq(seq, r_name, end, cigar).as_bytes())
                    .unwrap_or_else(|_| panic!("cannot write reads")),
                _ => 0,
            },
            Some(ReadAssign::Skipped(_, _)) => match &mut out_file_read_buffer.skipped {
                Some(handle) => handle
                    .write(self.dump_reads_seq(seq, r_name, end, cigar).as_bytes())
                    .unwrap_or_else(|_| panic!("cannot write reads")),
                _ => 0,
            },
            Some(ReadAssign::SoftClipped) => match &mut out_file_read_buffer.soft_clipped {
                Some(handle) => handle
                    .write(self.dump_reads_seq(seq, r_name, end, cigar).as_bytes())
                    .unwrap_or_else(|_| panic!("cannot write reads")),
                _ => 0,
            },
            Some(ReadAssign::OverhangFail) => match &mut out_file_read_buffer.overhang_fail {
                Some(handle) => handle
                    .write(self.dump_reads_seq(seq, r_name, end, cigar).as_bytes())
                    .unwrap_or_else(|_| panic!("cannot write reads")),
                _ => 0,
            },
        };
    }

    fn get_acceptor_donor(&self) -> (String, String) {
        let (mut left, mut right) = match (self.start, self.end) {
            (Some(left), Some(right)) => (left.to_string(), right.to_string()),
            (Some(left), None) => (left.to_string(), ".".to_string()),
            (None, Some(right)) => (".".to_string(), right.to_string()),
            (None, None) => (".".to_string(), ".".to_string()),
        };
        if self.strand == Strand::Plus {
            return (left, right);
        }
        return (right, left);
    }

    fn dump_junction_base(&self, ambigious_set: &HashSet<(String, i64)>) -> Vec<String> {
        let mut res = Vec::new();
        let mut exontype = ExonType::Acceptor;
        // get acceptor donnor!

        let mut pos = 0;
        let (donnor, acceptor) = self.get_acceptor_donor();

        let ambigious = if ambigious_set.contains(&(self.contig.clone(), self.start.unwrap_or(0)))
            | ambigious_set.contains(&(self.contig.clone(), self.end.unwrap_or(0)))
        {
            true
        } else {
            false
        };

        for (tr, intron) in &self.transcript_intron {
            res.push(format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                self.contig,
                self.gene_name,
                tr,
                intron,
                donnor,
                acceptor,
                self.strand,
                ambigious.to_string()
            ));
        }
        res
    }

    fn dump_junction_counter(
        &self,
        ambigious_set: &HashSet<(String, i64)>,
        junction_order: &Vec<SplicingEvent>,
    ) -> String {
        let base_vec = self.dump_junction_base(ambigious_set);
        let mut results = Vec::new();
        // TODO  order moove the allocation out

        let mut sub = Vec::new();
        for base in &base_vec {
            sub.clear();
            sub.push(base.to_owned());
            for even in junction_order {
                sub.push(self.counter_intron.get(&even).unwrap_or(&0).to_string());
            }
            results.push(sub.join("\t"))
        }
        results.join("\n")
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
}

/// we reads a gtf to generate an interval tree
/// HASHMAP(contig) -> IT(intron(start, end)) -> TreeDataIntron
pub fn interval_tree_from_gtfmap(
    gtf_map: &HashMap<String, HashMap<String, Vec<Exon>>>,
) -> Result<HashMap<String, interval_tree::IntervalTree<i64, TreeDataIntron>>, Box<dyn Error>> {
    let mut results: HashMap<String, interval_tree::IntervalTree<i64, TreeDataIntron>> =
        HashMap::new();
    let mut flag_exon_already_found = false;

    // by Iterating through the exon I can gather the Intron
    for (gene_id, tr_dico) in gtf_map {
        for (tr_id, exon_vec) in tr_dico.iter() {
            if exon_vec.len() == 1 {
                continue;
            }
            let strand = exon_vec[0].strand;
            let contig = &exon_vec[0].contig;
            flag_exon_already_found = false;

            results
                .entry(contig.to_string())
                .or_insert_with(|| interval_tree::IntervalTree::new());

            let exons = exon_vec
                .iter()
                .sorted_by(|x, y| x.start.cmp(&y.start))
                .collect::<Vec<&Exon>>();

            let mut cpt = 0;

            while cpt < exons.len() {
                let mut end_type = Some(ExonType::Acceptor);
                if strand == Strand::Minus {
                    end_type = Some(ExonType::Donnor);
                }

                let mut start_type = Some(ExonType::Donnor);
                if strand == Strand::Minus {
                    start_type = Some(ExonType::Acceptor);
                }

                let mut start_ = Some(0);
                let mut end_ = Some(0);
                let mut start_i = 0;
                let mut end_i = 0;
                let s = exons.len() - 1;

                if cpt == 0 {
                    start_i = exons[cpt].start - 1;
                    end_i = exons[cpt].start;
                    end_ = Some(exons[cpt].start);
                    start_type = None;
                    start_ = None;
                    update_tree(
                        &mut results,
                        &contig,
                        start_i,
                        end_i,
                        start_,
                        end_,
                        &strand,
                        gene_id.clone(),
                        tr_id.clone(),
                        start_type,
                        end_type,
                        Some(0),
                    );
                } else if cpt == exons.len() - 1 {
                    start_i = exons[cpt].end;
                    end_i = exons[cpt].end + 1;
                    end_type = None;
                    end_ = None;
                    start_ = Some(exons[cpt].end);

                    update_tree(
                        &mut results,
                        &contig,
                        start_i,
                        end_i,
                        start_,
                        end_,
                        &strand,
                        gene_id.clone(),
                        tr_id.clone(),
                        start_type,
                        end_type,
                        Some(1000),
                    );
                    cpt += 1;
                    continue;
                }

                let mut _i = cpt + 1;
                if strand == Strand::Minus {
                    _i = exons.len() - _i;
                }

                start_i = exons[cpt].end;
                end_i = exons[cpt + 1].start;
                start_ = Some(exons[cpt].end);
                end_ = Some(exons[cpt + 1].start);
                let mut end_type = Some(ExonType::Acceptor);
                if strand == Strand::Minus {
                    end_type = Some(ExonType::Donnor);
                }

                let mut start_type = Some(ExonType::Donnor);
                if strand == Strand::Minus {
                    start_type = Some(ExonType::Acceptor);
                }

                update_tree(
                    &mut results,
                    &contig,
                    start_i,
                    end_i,
                    start_,
                    end_,
                    &strand,
                    gene_id.clone(),
                    tr_id.clone(),
                    start_type,
                    end_type,
                    Some(_i as i32),
                );
                cpt += 1;
            }
        }
    }
    Ok(results)
}

pub fn update_tree(
    tree_map: &mut HashMap<String, interval_tree::IntervalTree<i64, TreeDataIntron>>,
    chr_: &str,
    interval_start: i64,
    interval_end: i64,
    start: Option<i64>,
    end: Option<i64>,
    strand: &Strand,
    gene_name: String,
    transcript_id: String,
    left_exon_type: Option<ExonType>,
    right_exon_type: Option<ExonType>,
    intron_n: Option<i32>,
) {
    let mut flag_exon_already_found = false;
    if let Some(this_tree) = tree_map.get_mut(chr_) {
        for (ref mut node) in this_tree.find_mut(interval_start..interval_end) {
            if (node.interval().start == interval_start) && (node.interval().end == interval_end) {
                if let n = node.data() {
                    if n.gene_name == gene_name {
                        n.transcript_intron
                            .push((transcript_id.clone(), intron_n.unwrap()));
                        flag_exon_already_found = true;
                    }
                }
            }
        }
        if !flag_exon_already_found {
            this_tree.insert(
                interval_start..interval_end,
                TreeDataIntron {
                    gene_name: gene_name,
                    transcript_intron: vec![(transcript_id, intron_n.unwrap())],
                    start: start,
                    end: end,
                    strand: *strand,
                    contig: chr_.to_string(),
                    start_type: left_exon_type,
                    end_type: right_exon_type,
                    counter_start: HashMap::new(),
                    counter_end: HashMap::new(),
                    counter_intron: HashMap::new(),
                },
            );
        }
    }
}

// HashMap<String, interval_tree::IntervalTree<i64, TreeDataIntron>>
pub fn update_tree_from_bam(
    hash_tree: &mut HashMap<String, interval_tree::IntervalTree<i64, TreeDataIntron>>,
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
    out_file_read_buffer: &mut ReadToWriteHandle,
    junction_valid: &HashMap<(String, i64, i64), Strand>,
    valid_j_gene: &HashMap<String, HashSet<(i64, i64)>>,
) -> Result<(), Box<dyn Error>> {
    let mut read_name: Option<String> = None;
    let mut seq: Option<Vec<u8>> = None;
    let mut sequence: Option<String> = None;

    let mut counter: i64;
    let mut record: Record = Record::new();
    let mut pos_s: i64;
    let mut pos_e: i64;
    let mut cig: Cigar;
    let mut flag: u16;
    let mut bam = IndexedReader::from_path(bam_file).unwrap();

    for (contig, subtree) in hash_tree.iter_mut() {
        match bam.fetch(&contig) {
            Ok(_) => (),
            Err(_) => {
                println!("WARNING {} not found in the bam file check", contig);
                continue;
            }
        }
        while let Some(r) = bam.read(&mut record) {
            //for r in bam.records() {
            //record = r.unwrap();
            pos_s = record.pos();
            cig = Cigar::from_str(&record.cigar().to_string()).unwrap();
            pos_e = cig.get_end_of_aln(&pos_s);
            flag = record.flags();

            // QC
            if (!check_flag(flag, flag_in, flag_out)) || (record.mapq() < mapq) {
                continue;
            }

            read_name =
                Some(String::from_utf8(record.qname().to_vec()).expect("cannot parse read name"));
            let seq = record.seq().as_bytes();
            sequence = Some(String::from_utf8(seq).expect("cannot parse sequence"));

            if let Some(read_strand) = library_type.get_strand(flag) {
                for (ref mut exon) in subtree.find_mut((pos_s - 1)..(pos_e + 1)) {
                    exon.data().parse_read(
                        pos_s,
                        pos_e,
                        &cig,
                        &read_strand,
                        overhang,
                        out_file_read_buffer,
                        read_name.clone(),
                        sequence.clone(),
                        junction_valid,
                        valid_j_gene,
                    );
                }
            }
        }
    }
    Ok(())
}

/*
pub fn dump_tree_junction (hash_tree: &HashMap<String, interval_tree::IntervalTree<i64, TreeDataIntron>>,
        out_file: &str,
        ) -> () {
            let presorted = format!("{}.presorted", out_file);
            let header = "Contig\tGene\tTranscript\tIntron\tDonnor\tAcceptor\tStrand\tAmbiguous\tspliced\tunspliced\tclipped\texon_other\tskipped\twrong_strand\te_isoform\n".as_bytes();

            for (contig, subtree) in hash_tree {
                for node in subtree.find(i64::MIN..i64::MAX) {

            }
        }
    }
*/

pub fn dump_tree_to_cat_results(
    hash_tree: &HashMap<String, interval_tree::IntervalTree<i64, TreeDataIntron>>,
    out_cat: &str,
    out_junction: &str,
    junction_ambigious: &HashSet<(String, i64)>,
    junction_order: &Vec<SplicingEvent>,
) -> () {
    let presorted_cat = format!("{}.presorted", out_cat);
    let presorted_j = format!("{}.presorted", out_junction);
    let header_cat =
        "Contig\tPosition\tGeneID\tTranscriptID\tstrand\tExonEndType\tCategory\tReadCount\n"
            .as_bytes();
    let header_j = "Contig\tGene\tTranscript\tIntron\tDonnor\tAcceptor\tStrand\tAmbiguous\tspliced\tunspliced\tclipped\texon_other\tskipped\twrong_strand\te_isoform\n".as_bytes();

    {
        let file_cat =
            File::create_new(presorted_cat.clone()).expect("output file should not exist.");
        let mut stream_cat = BufWriter::new(file_cat);

        let file_j = File::create_new(presorted_j.clone()).expect("output file should not exist.");
        let mut stream_j = BufWriter::new(file_j);
        //stream.write(header);

        for (contig, subtree) in hash_tree {
            // get ALL entry
            for node in subtree.find(i64::MIN..i64::MAX) {
                if node.data().start.is_some() {
                    stream_cat.write_all(node.data().dump_counter(false).as_bytes());
                    stream_cat.write_all("\n".as_bytes());
                }
                if node.data().end.is_some() {
                    stream_cat.write_all(node.data().dump_counter(true).as_bytes());
                    stream_cat.write_all("\n".as_bytes());
                }

                stream_j.write_all(
                    node.data()
                        .dump_junction_counter(junction_ambigious, junction_order)
                        .as_bytes(),
                );
                stream_j.write_all("\n".as_bytes());
            }
        }
        stream_cat.flush().unwrap();
        stream_j.flush().unwrap();
    }

    sort_and_clear_cat(&presorted_cat, out_cat, header_cat);
    sort_and_clear_junction(&presorted_j, out_junction, header_j)

    /*
    let file = File::create_new(out_cat).expect("output file should not exist.");

    let mut file = OpenOptions::new()
        .write(true)
        .append(true)
        .open(out_cat).expect("could not open file");
    file.write_all(header_cat);
    //let mut stream = BufWriter::new(file);

    let sort = Command::new("sort")
        .args([
            "-k1,1",
            "-k3,3",
            "-k4,4",
            "-k2,2n",
            "-k7,7",
            presorted_cat.as_str(),
        ])
        .stdout(file)
        .spawn()
        .expect("sort command failed ")
        .wait()
        .expect("sort command failed ");
    //stream.flush().unwrap();

    let rm_out = Command::new("rm")
        .args(["-f", presorted_cat.as_str()])
        .output()
        .expect("failed to remove presorted");*/
}

fn sort_and_clear_junction(presorted: &str, out: &str, header: &[u8]) -> () {
    let file = File::create_new(out).expect("output file should not exist.");
    let mut file = OpenOptions::new()
        .write(true)
        .append(true)
        .open(out)
        .expect("could not open file");
    file.write_all(header);
    //let mut stream = BufWriter::new(file);

    let sort = Command::new("sort")
        .args(["-k1,1", "-k2,2", "-k3,3", "-k4,4n", presorted])
        .stdout(file)
        .spawn()
        .expect("sort command failed ")
        .wait()
        .expect("sort command failed ");
    //stream.flush().unwrap();

    let rm_out = Command::new("rm")
        .args(["-f", presorted])
        .output()
        .expect("failed to remove presorted");
}

fn sort_and_clear_cat(presorted: &str, out: &str, header: &[u8]) -> () {
    let file = File::create_new(out).expect("output file should not exist.");

    let mut file = OpenOptions::new()
        .write(true)
        .append(true)
        .open(out)
        .expect("could not open file");
    file.write_all(header);
    //let mut stream = BufWriter::new(file);

    let sort = Command::new("sort")
        .args(["-k1,1", "-k3,3", "-k4,4", "-k2,2n", "-k7,7", presorted])
        .stdout(file)
        .spawn()
        .expect("sort command failed ")
        .wait()
        .expect("sort command failed ");
    //stream.flush().unwrap();

    let rm_out = Command::new("rm")
        .args(["-f", presorted])
        .output()
        .expect("failed to remove presorted");
}

#[cfg(test)]
mod tests_it {
    use super::*;
    use crate::common::gtf_::{get_junction_from_gtf, gtf_to_hashmap};

    #[test]
    fn test_1() {
        let file = "/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/debug/fbgn0001313.gtf"
            .to_string();
        let gtf_hashmap = gtf_to_hashmap(&file).expect("failed to parse gtf");
        let mut hash_tree = interval_tree_from_gtfmap(&gtf_hashmap)
            .expect("failed to generate the hash tree from gtf");

        assert_eq!(true, true);
    }
    #[test]
    fn test_2() {
        let file = "/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/debug/FBgn0000052.gtf"
            .to_string();
        let gtf_hashmap = gtf_to_hashmap(&file).expect("failed to parse gtf");
        let mut hash_tree = interval_tree_from_gtfmap(&gtf_hashmap)
            .expect("failed to generate the hash tree from gtf");

        assert_eq!(true, true);
    }

    #[test]
    fn test_3() {
        let x = read_toassign(
            Strand::Minus,
            Some(21681343),
            Some(ExonType::Donnor),
            21681343,
            21681343 + 188,
            &Cigar::from("63S188M"),
            &Strand::Minus,
            1,
        );
    }
} //21589349, 21589613
  /*        let start_map= read_toassign(
      self.strand,
      self.start,
      self.start_type,
      aln_start,
      aln_end,
      cigar,
      read_strand,
      overhang,
  ); */

/*
&mut self,
aln_start: i64,
aln_end: i64,
cigar: &Cigar,
read_strand: &Strand,
overhang: i64,
out_file_read_buffer: &mut ReadToWriteHandle,
read_name: Option<String>,
sequence: Option<String>,
valid_junction: &HashMap<(String, i64, i64),  Strand>
 */
