#![allow(irrefutable_let_patterns)]
use crate::common::error::OmniError;
//use crate::common::point::{get_attr_id, InsideCounter};
use crate::common::utils::{
    Exon, ExonType, ReadAssign, ReadToWriteHandle, SplicingEvent, read_toassign,
};

use CigarParser::cigar::Cigar;
use bio::alignment::sparse::HashMapFx;
use bio::bio_types::annot::contig;
use bio::data_structures::interval_tree;
use bio::seq_analysis;
use bio::utils::Interval;
use itertools::Itertools;
use rust_htslib::bam::record::Record;
use rust_htslib::bam::{IndexedReader, Read, header};
use std::collections::{HashMap, HashSet, hash_map};
use std::error::Error;
use std::fs::{File, OpenOptions};
use std::io::BufRead;
use std::io::Write;
use std::io::{BufReader, BufWriter};
use std::process::{Command, Stdio};
use std::str::FromStr;
use strand_specifier_lib::Strand;
use strand_specifier_lib::{LibType, check_flag};
use log::{info, debug, error, trace, warn};

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

    start_bis: Option<i64>,
    start: Option<i64>,
    end: Option<i64>,
    end_bis: Option<i64>,

    strand: Strand,
    contig: String,
    start_type: Option<ExonType>,
    end_type: Option<ExonType>,
    /// count exon end level Raw event with position included
    counter_start: HashMap<ReadAssign, i32>,
    /// count  exon end level Raw event with position included
    counter_end: HashMap<ReadAssign, i32>,
    /// count  exon end level SplicingEvent event with position included
    counter_splicingevent_start: HashMap<SplicingEvent, i32>,
    /// count  exon end level SplicingEvent event with position included
    counter_splicingevent_end: HashMap<SplicingEvent, i32>,
    /// count Intron level splicingevent event with position included
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
    ) -> Result<(), OmniError> {

        let seq = match sequence {
            Some(seq) => seq,
            _ => "NoSeq".to_string(),
        };
        let r_name = match read_name {
            Some(seq) => seq,
            _ => "NoSeq".to_string(),
        };
        
        // TODO thought could delegate that to treeDATA even more. HAve to think about it.
        let start_map = read_toassign(
            self.strand,
            self.start,
            self.start_bis,
            self.start_type,
            aln_start,
            aln_end,
            cigar,
            read_strand,
            overhang,
        )?;
        self.write_to_read_file(start_map, out_file_read_buffer, false, &seq, &r_name, cigar);
        TreeDataIntron::update_counter(&mut self.counter_start, start_map);

        let end_map = read_toassign(
            self.strand,
            self.end,
            self.end_bis,
            self.end_type,
            aln_start,
            aln_end,
            cigar,
            read_strand,
            overhang,
        )?;

        TreeDataIntron::update_counter(&mut self.counter_end, end_map);
        self.write_to_read_file(end_map, out_file_read_buffer, true, &seq, &r_name, cigar);



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
                //counter_splicingevent_start: HashMap::new(),
                TreeDataIntron::update_splicing_counter(
                    &mut self.counter_splicingevent_start,
                    Some(x),
                );
                TreeDataIntron::update_splicing_counter(&mut self.counter_intron, Some(x));
            }
            (None, Some(x)) => {
                TreeDataIntron::update_splicing_counter(
                    &mut self.counter_splicingevent_end,
                    Some(x),
                );
                TreeDataIntron::update_splicing_counter(&mut self.counter_intron, Some(x))
            }
            (Some(SplicingEvent::Spliced), Some(SplicingEvent::Spliced)) => {
                TreeDataIntron::update_splicing_counter(
                    &mut self.counter_splicingevent_start,
                    Some(SplicingEvent::Spliced),
                );
                TreeDataIntron::update_splicing_counter(
                    &mut self.counter_splicingevent_end,
                    Some(SplicingEvent::Spliced),
                );

                TreeDataIntron::update_splicing_counter(
                    &mut self.counter_intron,
                    Some(SplicingEvent::Spliced),
                )
            }
            (Some(SplicingEvent::Unspliced), Some(SplicingEvent::Unspliced)) => {
                TreeDataIntron::update_splicing_counter(
                    &mut self.counter_splicingevent_start,
                    Some(SplicingEvent::Unspliced),
                );
                TreeDataIntron::update_splicing_counter(
                    &mut self.counter_splicingevent_end,
                    Some(SplicingEvent::Unspliced),
                );

                TreeDataIntron::update_splicing_counter(
                    &mut self.counter_intron,
                    Some(SplicingEvent::Unspliced),
                )
            }
            (Some(left), Some(right)) => {
                TreeDataIntron::update_splicing_counter(
                    &mut self.counter_splicingevent_start,
                    Some(left),
                );
                TreeDataIntron::update_splicing_counter(
                    &mut self.counter_splicingevent_end,
                    Some(right),
                );

                TreeDataIntron::update_splicing_counter(
                    &mut self.counter_intron,
                    SplicingEvent::merge(Some(left), Some(right)),
                )
            }
            (None, None) => (),
            (_, _) => (),
        }
        Ok(())
    }

    /// Update the counter if item is None do Nothing.
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

    fn dump_base(&self, end: bool) -> Result<Vec<String>, OmniError> {
        let mut res = Vec::new();
        let mut exontype = ExonType::Acceptor;
        let mut pos = 0;
        if end {
            exontype = self.end_type.ok_or(OmniError::Expect("interger expect find None in dump base".to_string()))?;
            pos = self.end.ok_or(OmniError::Expect("interger expect find None in dump base".to_string()))?;
        } else {
            exontype = self.start_type.ok_or(OmniError::Expect("interger expect find None in dump base".to_string()))?;
            pos = self.start.ok_or(OmniError::Expect("interger expect find None in dump base".to_string()))?;
        }
        for (tr, intron) in &self.transcript_intron {
            res.push(format!(
                "{}\t{}\t{}\t{}\t{}\t{}",
                self.contig, pos, self.gene_name, tr, self.strand, exontype
            ));
        }
        Ok(res)
    }

    fn dump_reads_seq(&self, sequence: &str, seqname: &str, end: bool, cigar: &Cigar) -> Result<String, OmniError> {
        ///
        /// :param end: if set to true look at the end of the sequence (end > start)
        let mut base_vec = self.dump_base(end)?;

        for e in &mut base_vec {
            e.push_str(&format!("\t{}\t{}\t{}", seqname, cigar, sequence))
        }

        return Ok(format!("{}\n", base_vec.join("\n")));
    }

    fn write_to_read_file(
        &self,
        read_assign: Option<ReadAssign>,
        out_file_read_buffer: &mut ReadToWriteHandle,
        end: bool,
        seq: &str,
        r_name: &str,
        cigar: &Cigar,
    ) -> Result<(), OmniError>{
        match read_assign {
            None => 0,
            Some(ReadAssign::SkippedUnrelated(_, _)) => match &mut out_file_read_buffer.empty {
                Some(handle) => handle
                    .write(self.dump_reads_seq(seq, r_name, end, cigar)?.as_bytes())?,
                _ => 0,
            },
            Some(ReadAssign::Empty) => match &mut out_file_read_buffer.empty {
                Some(handle) => handle
                    .write(self.dump_reads_seq(seq, r_name, end, cigar)?.as_bytes())?,
                _ => 0,
            },
            Some(ReadAssign::ReadThrough) => match &mut out_file_read_buffer.read_through {
                Some(handle) => handle
                    .write(self.dump_reads_seq(seq, r_name, end, cigar)?.as_bytes())?,
                _ => 0,
            },
            Some(ReadAssign::ReadJunction(_, _)) => match &mut out_file_read_buffer.read_junction {
                Some(handle) => handle
                    .write(self.dump_reads_seq(seq, r_name, end, cigar)?.as_bytes())?,
                _ => 0,
            },
            Some(ReadAssign::Unexpected) => match &mut out_file_read_buffer.unexpected {
                Some(handle) => handle
                    .write(self.dump_reads_seq(seq, r_name, end, cigar)?.as_bytes())?,
                _ => 0,
            },
            Some(ReadAssign::FailPosFilter) => match &mut out_file_read_buffer.fail_pos_filter {
                Some(handle) => handle
                    .write(self.dump_reads_seq(seq, r_name, end, cigar)?.as_bytes())?,
                _ => 0,
            },
            Some(ReadAssign::WrongStrand) => match &mut out_file_read_buffer.wrong_strand {
                Some(handle) => handle
                    .write(self.dump_reads_seq(seq, r_name, end, cigar)?.as_bytes())?,
                _ => 0,
            },
            Some(ReadAssign::FailQc) => match &mut out_file_read_buffer.fail_qc {
                Some(handle) => handle
                    .write(self.dump_reads_seq(seq, r_name, end, cigar)?.as_bytes())?,
                _ => 0,
            },
            Some(ReadAssign::EmptyPileup) => match &mut out_file_read_buffer.empty_pileup {
                Some(handle) => handle
                    .write(self.dump_reads_seq(seq, r_name, end, cigar)?.as_bytes())?,
                _ => 0,
            },
            Some(ReadAssign::Skipped(_, _)) => match &mut out_file_read_buffer.skipped {
                Some(handle) => handle
                    .write(self.dump_reads_seq(seq, r_name, end, cigar)?.as_bytes())?,
                _ => 0,
            },
            Some(ReadAssign::SoftClipped) => match &mut out_file_read_buffer.soft_clipped {
                Some(handle) => handle
                    .write(self.dump_reads_seq(seq, r_name, end, cigar)?.as_bytes())?,
                _ => 0,
            },
            Some(ReadAssign::OverhangFail) => match &mut out_file_read_buffer.overhang_fail {
                Some(handle) => handle
                    .write(self.dump_reads_seq(seq, r_name, end, cigar)?.as_bytes())?,
                _ => 0,
            },
        };
        Ok(())
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
/*    counter_splicingevent_start: HashMap<SplicingEvent, i32>,
    /// count  exon end level SplicingEvent event with position included
    counter_splicingevent_end: HashMap<SplicingEvent, i32>, */



//dump_exon_counter_loop(counter: &HashMap<SplicingEvent, i32>, base_vec: &Vec<String>, results: &mut Vec<String>, sub: &mut Vec<String>, junction_order: &Vec<SplicingEvent>)

    fn dump_exon_counter(&self,
                        ambigious_set: &HashSet<(String, i64)>,
                        junction_order: &Vec<SplicingEvent>
                    ) -> Result<String, OmniError> {

        let base_vec_end = self.dump_exon_base(ambigious_set, true)?;
        let mut results = Vec::new();
        let mut sub = Vec::new();
        dump_exon_counter_loop(&self.counter_splicingevent_end, &self.dump_exon_base(ambigious_set, true)?, &mut results, &mut sub, junction_order);
        dump_exon_counter_loop(&self.counter_splicingevent_start, &self.dump_exon_base(ambigious_set, false)?, &mut results, &mut sub, junction_order);

        Ok(results.join("\n"))

    }

    fn dump_exon_base(&self, ambigious_set: &HashSet<(String, i64)>, end: bool) -> Result<Vec<String>, OmniError> {


        let mut res = Vec::new();
        let mut exontype = ExonType::Acceptor;
        // get acceptor donnor!

        let mut pos = 0;
        //let (donnor, acceptor) = self.get_acceptor_donor();

        let ambigious = if ambigious_set.contains(&(self.contig.clone(), self.start.unwrap_or(0)))
            | ambigious_set.contains(&(self.contig.clone(), self.end.unwrap_or(0)))
        {
            true
        } else {
            false
        };


/*    gene_name: String,
    transcript_intron: Vec<(String, i32)>,

    start_bis: Option<i64>,
    start: Option<i64>,
    end: Option<i64>,
    end_bis: Option<i64>,

    strand: Strand,
    contig: String,
    start_type: Option<ExonType>,
    end_type: Option<ExonType>,
*/

        // Contig Gene Transcript Exon_number Exon_type Position Next Strand Ambiguous
      
    if end == true{

        if self.end.is_none() {
            return Ok(Vec::new())
        };

        for (tr, intron) in &self.transcript_intron {
            res.push(format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                self.contig,
                self.gene_name,
                tr,
                match self.end_type{
                    Some(ExonType::Acceptor) => *intron + 1,
                    Some(ExonType::Donnor) => *intron,
                    None => {error!("expected value for end type found None"); return Err(OmniError::Expect("value expected None found".to_string()))}
                },
                self.end_type.ok_or(OmniError::Expect("value expected None found".to_string()))?,
                self.end.ok_or(OmniError::Expect("value expected None found".to_string()))?,
                self.start.ok_or(OmniError::Expect("value expected None found".to_string()))?,
                self.strand,
                ambigious.to_string()
            ));
        }
    }
    else{
        if self.start.is_none() {
            return Ok(Vec::new())
        };

                for (tr, intron) in &self.transcript_intron {
            res.push(format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                self.contig,
                self.gene_name,
                tr,
                match self.start_type{
                    Some(ExonType::Acceptor) => *intron + 1,
                    Some(ExonType::Donnor) => *intron,
                    None => {error!("expected value for start type found None"); return Err(OmniError::Expect("value expected None found".to_string()))}
                },
                self.start_type.ok_or(OmniError::Expect("value expected None found".to_string()))?,
                self.start.ok_or(OmniError::Expect("value expected None found".to_string()))?,
                self.end.ok_or(OmniError::Expect("value expected None found".to_string()))?,
                
                self.strand,
                ambigious.to_string()
            ));
        }

    }
        Ok(res)
    }

    fn dump_counter(&self, end: bool) -> Result<String, OmniError> {
        let base_vec = self.dump_base(end)?;
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

        Ok(results.join("\n"))
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
                let mut start_ = Some(0);
                let mut end_ = Some(0);
                let mut start_b = Some(0);
                let mut end_b_ = Some(0);
                let mut start_i = 0;
                let mut end_i = 0;
                let mut s = 0;

            while cpt < exons.len() {
                let mut end_type = Some(ExonType::Acceptor);
                if strand == Strand::Minus {
                    end_type = Some(ExonType::Donnor);
                }

                let mut start_type = Some(ExonType::Donnor);
                if strand == Strand::Minus {
                    start_type = Some(ExonType::Acceptor);
                }
                
                s = exons.len() - 1;


                if cpt == 0 {
                    start_i = exons[cpt].start - 1;
                    end_i = exons[cpt].start + 1;
                    end_ = Some(exons[cpt].start);
                    end_b_ = Some(exons[cpt].end);
                    start_type = None;
                    start_ = None;
                    start_b = None;

                    update_tree(
                        &mut results,
                        &contig,
                        start_i,
                        end_i,
                        start_,
                        start_b, 
                        end_,
                        end_b_,
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
                    end_b_ = None;
                    start_ = Some(exons[cpt].end);
                    start_b = Some(exons[cpt].start);

                    update_tree(
                        &mut results,
                        &contig,
                        start_i,
                        end_i,
                        start_,
                        start_b,
                        end_,
                        end_b_,
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
                start_b = Some(exons[cpt].start);
                end_ = Some(exons[cpt + 1].start);
                end_b_ = Some(exons[cpt + 1].end);

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
                    start_b,
                    end_,
                    end_b_,
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
    start_bis: Option<i64>,
    end: Option<i64>,
    end_bis: Option<i64>,
    strand: &Strand,
    gene_name: String,
    transcript_id: String,
    left_exon_type: Option<ExonType>,
    right_exon_type: Option<ExonType>,
    intron_n: Option<i32>,
) -> Result<(), OmniError> {
    let mut flag_exon_already_found = false;
    if let Some(this_tree) = tree_map.get_mut(chr_) {
        for (ref mut node) in this_tree.find_mut(interval_start..interval_end) {
            if (node.interval().start == interval_start) && (node.interval().end == interval_end) {
                // TODO FIX this
                // THIS cause a bug in case of ambigious junction that are very close to each others
                if let n = node.data() {
                    if n.gene_name == gene_name && n.start == start && n.end == end {
                        n.transcript_intron
                            .push((transcript_id.clone(),
                                 intron_n.ok_or(OmniError::Expect("expect integer None found".to_string()))?)
                                );
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
                    transcript_intron: vec![(transcript_id,
                                            intron_n
                                            .ok_or(OmniError::Expect("expect integer None found".to_string()))?)],
                    start: start,
                    end: end,
                    strand: *strand,
                    start_bis: start_bis,
                    end_bis: end_bis,
                    contig: chr_.to_string(),
                    start_type: left_exon_type,
                    end_type: right_exon_type,

                    counter_start: HashMap::new(),
                    counter_end: HashMap::new(),

                    counter_intron: HashMap::new(),

                    counter_splicingevent_end: HashMap::new(),
                    counter_splicingevent_start: HashMap::new(),
                },
            );
        }
    }
    Ok(())
}


fn dump_exon_counter_loop(counter: &HashMap<SplicingEvent, i32>, base_vec: &Vec<String>, results: &mut Vec<String>, sub: &mut Vec<String>, junction_order: &Vec<SplicingEvent>) -> (){
    for base in base_vec {
        sub.clear();
        sub.push(base.to_owned());
        for even in junction_order {
            sub.push(counter.get(&even).unwrap_or(&0).to_string());
        }
        results.push(sub.join("\t"))
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
) -> Result<(), OmniError> {
    let mut read_name: Option<String> = None;
    let mut seq: Vec<u8>;
    let mut sequence: Option<String> = None;

    let mut counter: i64;
    let mut record: Record = Record::new();
    let mut pos_s: i64;
    let mut pos_e: i64;
    let mut cig: Cigar;
    let mut flag: u16;
    let mut bam = IndexedReader::from_path(bam_file)?;

    for (contig, subtree) in hash_tree.iter_mut() {
        
        match bam.fetch(&contig) {
            Ok(_) => (),
            Err(_) => {
                println!("WARNING {} not found in the bam file check", contig);
                continue;
            }
        }
        
        info!("parsing reads mapping to contig: {}", contig);
        while let Some(r) = bam.read(&mut record) {

            pos_s = record.pos();
            cig = match Cigar::from_str(&record.cigar().to_string()) {
                Ok(c) => c,
                Err(_) => {
                    error!("unproper CigarString skipping this reads");
                    continue;
                }
            };

            pos_e = cig.get_end_of_aln(pos_s);
            flag = record.flags();

            // QC
            if (!check_flag(flag, flag_in, flag_out)) || (record.mapq() < mapq) {
                continue;
            }

            read_name =
                Some(String::from_utf8(record.qname().to_vec())?);

            seq = record.seq().as_bytes();
            sequence = Some(String::from_utf8(seq)?);

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
            } else {error!("failed to determined strand : {:?} flag:{:?}", read_name, flag);}
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


/*pub fn dump_raw(    hash_tree: &HashMap<String, interval_tree::IntervalTree<i64, TreeDataIntron>>,
    out_cat: &str,
    out_junction: &str,
    junction_ambigious: &HashSet<(String, i64)>,
    junction_order: &Vec<SplicingEvent>) -> Result<(), OmniError>{Ok(())}


pub fn dump_exons() -> Result<(), OmniError>{Ok(())}
pub fn dump_junctions() -> Result<(), OmniError>{Ok(())}*/


pub fn dump_tree_into_raw_exon_junction(    
    hash_tree: &HashMap<String, interval_tree::IntervalTree<i64, TreeDataIntron>>,
    out_raw: &str,
    out_exons: &str,
    out_junction: &str,
    junction_ambigious: &HashSet<(String, i64)>,
    junction_order: &Vec<SplicingEvent>,) -> Result<(), OmniError> {
    
    let presorted_raw = format!("{}.presorted", out_raw);
    let presorted_junction = format!("{}.presorted", out_junction);
    let presorted_exons = format!("{}.presorted", out_exons);

    /*        SplicingEvent::Spliced,
        SplicingEvent::Unspliced,
        SplicingEvent::Clipped,
        SplicingEvent::ExonOther,
        SplicingEvent::Skipped,
        SplicingEvent::SkippedUnrelated,
        SplicingEvent::WrongStrand,
        SplicingEvent::Isoform, */

    let header_raw =
        "Contig\tPosition\tGeneID\tTranscriptID\tstrand\tExonEndType\tCategory\tReadCount\n"
            .as_bytes();
    let header_junction = "Contig\tGene\tTranscript\tIntron\tDonnor\tAcceptor\tStrand\tAmbiguous\tSpliced\tUnspliced\tClipped\tExon_other\tSkipped\tSkippedUnrelated\tWrong_strand\tE_isoform\n".as_bytes();
    let header_exons = "Contig\tGene\tTranscript\tExon_number\tExon_type\tPosition\tNext\tStrand\tAmbiguous\tSpliced\tUnspliced\tClipped\tExon_other\tSkipped\tSkippedUnrelated\tWrong_strand\tE_isoform\n".as_bytes();

    {

        let file_raw =
                File::create_new(presorted_raw.clone()).map_err(|e| OmniError::Expect("output file should not exist".to_string()))?;
        let mut stream_raw = BufWriter::new(file_raw);

        let file_junction = 
                File::create_new(presorted_junction.clone()).map_err(|e| OmniError::Expect("output file should not exist".to_string()))?;
        let mut stream_junction = BufWriter::new(file_junction);

        let file_exons = 
                File::create_new(presorted_exons.clone()).map_err(|e| OmniError::Expect("output file should not exist".to_string()))?;
        let mut stream_exons = BufWriter::new(file_exons);

        for (contig, subtree) in hash_tree {
            // get ALL entry
            for node in subtree.find(i64::MIN..i64::MAX) {
                if node.data().start.is_some() {
                    stream_raw.write_all(node.data().dump_counter(false)?.as_bytes());
                    stream_raw.write_all("\n".as_bytes());
                }
                if node.data().end.is_some() {
                    stream_raw.write_all(node.data().dump_counter(true)?.as_bytes());
                    stream_raw.write_all("\n".as_bytes());
                }

                stream_junction.write_all(
                    node.data()
                        .dump_junction_counter(junction_ambigious, junction_order)
                        .as_bytes(),
                );
                stream_junction.write_all("\n".as_bytes());

                stream_exons.write_all(
                    node.data().dump_exon_counter(junction_ambigious, junction_order)?
                    .as_bytes(),

                );
                stream_exons.write_all("\n".as_bytes());


            }
        }
        stream_raw.flush()?;
        stream_junction.flush()?;
        stream_exons.flush()?;

    }

    info!("sorting the results files");
    sort_and_clear_cat(&presorted_raw, out_raw, header_raw);
    sort_and_clear_junction(&presorted_junction, out_junction, header_junction);
    sort_and_clear_junction(&presorted_exons, out_exons, header_exons);
    info!("done");

    Ok(())
}

pub fn dump_tree_to_cat_results(
    hash_tree: &HashMap<String, interval_tree::IntervalTree<i64, TreeDataIntron>>,
    out_cat: &str,
    out_junction: &str,
    junction_ambigious: &HashSet<(String, i64)>,
    junction_order: &Vec<SplicingEvent>,
) -> Result<(), OmniError> {
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
                    stream_cat.write_all(node.data().dump_counter(false)?.as_bytes());
                    stream_cat.write_all("\n".as_bytes());
                }
                if node.data().end.is_some() {
                    stream_cat.write_all(node.data().dump_counter(true)?.as_bytes());
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
        stream_cat.flush()?;
        stream_j.flush()?;
    }

    sort_and_clear_cat(&presorted_cat, out_cat, header_cat);
    sort_and_clear_junction(&presorted_j, out_junction, header_j);
    Ok(())

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
            Some(21681345),
            Some(ExonType::Donnor),
            21681343,
            21681343 + 188,
            &Cigar::from("63S188M"),
            &Strand::Minus,
            1,
        );
    }
} 


