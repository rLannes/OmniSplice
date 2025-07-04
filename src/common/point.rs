use std::cmp::Eq;
use std::cmp::Ordering;
use std::fs::File;
use std::hash::Hash;
//mod common;
use crate::common::utils::{read_toassign, ExonType, ReadAssign};
use rust_htslib::bam::record::Record;
use std::collections::HashMap;
use std::collections::HashSet;
use std::error::Error;
use std::io::BufRead;
use std::io::Write;
use std::io::{BufReader, BufWriter};
use strand_specifier_lib::Strand;
use CigarParser::cigar::Cigar;

#[derive(Debug, Clone)]
pub struct PointContainer {
    pub points: Vec<Point>,
}

pub struct PointContainerIterator<'a> {
    my_struct: &'a PointContainer,
    index: usize,
}

impl<'a> Iterator for PointContainerIterator<'a> {
    // We can refer to this type using Self::Item
    type Item = &'a Point;
    fn next(&mut self) -> Option<Self::Item> {
        if self.index < self.my_struct.points.len() {
            let item = &self.my_struct.points[self.index];
            self.index += 1;
            Some(item)
        } else {
            None
        }
    }
}
/*
impl PointContainer {
    pub fn parse_reads(
        self: &mut Self,
        contig: &str,
        aln_start: i64,
        aln_end: i64,
        cigar: &Cigar,
        flag: &u16,
        read_strand: &Strand,
        overhang: i64,
        output_file_read: &mut Option<BufWriter<File>>,
        record: &Record,
        clipped: bool,
    ) -> Result<(), ()> {
        //get the index of Point matching plus 1 delta to capture Softclipped
        let lower_bound = self.lower_bound(aln_start - 1);
        let upper_bound = self.upper_bound(aln_end + 1);
        //println!("{} {}" , lower_bound, upper_bound);

        //if (lower_bound == upper_bound) | (upper_bound <= 0) | (upper_bound >= self.points.len()) {

        //    return Err(());
        // }

        for p_indices in lower_bound..upper_bound {
            //feature_strand: Strand, feature_pos: i64, feature_exontype
            //self.points[p_indices].parse_read(aln_start, aln_end, cigar, flag, read_strand, overhang);
            if let Some(read_assign) = read_toassign(
                self.points[p_indices].strand,
                self.points[p_indices].pos,
                self.points[p_indices].exon_type,
                aln_start,
                aln_end,
                cigar,
                read_strand,
                overhang,
            ) {
                self.points[p_indices].insert_in_counter(read_assign);
                if let Some(buf_writer) = output_file_read {
                    let read_name =
                        String::from_utf8(record.qname().to_vec()).expect("cannot parse read name");
                    let seq = record.seq().as_bytes();
                    let sequence = String::from_utf8(seq).expect("cannot parse sequence");
                    if clipped {
                        if read_assign == ReadAssign::SoftClipped {
                            let _ = buf_writer.write(
                                format!(
                                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                                    read_name,
                                    cigar,
                                    flag,
                                    contig,
                                    aln_start,
                                    read_assign,
                                    self.points[p_indices].pos,
                                    self.points[p_indices].exon_type,
                                    self.points[p_indices].strand,
                                    self.points[p_indices].gene_name,
                                    self.points[p_indices].transcript_id,
                                    sequence
                                )
                                .as_bytes(),
                            );
                        }
                    } else {
                        let _ = buf_writer.write(
                            format!(
                                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                                read_name,
                                cigar,
                                flag,
                                contig,
                                aln_start,
                                read_assign,
                                self.points[p_indices].pos,
                                self.points[p_indices].exon_type,
                                self.points[p_indices].strand,
                                self.points[p_indices].gene_name,
                                self.points[p_indices].transcript_id,
                                sequence
                            )
                            .as_bytes(),
                        );
                    }
                }
            }
        }
        Ok(())
    }
}
*/
impl PointContainer {
    pub fn iter(&self) -> PointContainerIterator {
        PointContainerIterator {
            my_struct: self,
            index: 0,
        }
    }

    pub fn new() -> Self {
        PointContainer {
            points: Vec::with_capacity(1_000),
        }
    }

    pub fn lower_bound(self: &Self, thr: i64) -> usize {
        //let length = self.points.len();
        let i = match my_bs(&self.points, &thr) {
            Ok(ref mut i) => {
                while (*i > 0) && (self.points[*i - 1] >= thr) {
                    *i -= 1;
                }
                *i
            }
            Err(i) => i,
        };
        i
    }

    pub fn upper_bound(self: &Self, thr: i64) -> usize {
        let length = self.points.len();
        let j = match my_bs(&self.points, &thr) {
            Ok(ref mut i) => {
                while (*i < length - 2) && (self.points[*i] <= thr) {
                    *i += 1;
                }
                *i
            }
            Err(ref mut i) => {
                while (*i < length - 2) && (self.points[*i] <= thr) {
                    *i += 1;
                }
                *i
            }
        };

        j.min(length - 1)
    }

    pub fn push(self: &mut Self, point: Point) {
        self.points.push(point)
    }

    pub fn sort(self: &mut Self) {
        self.points.sort();
    }
}

#[derive(Debug, Clone)]
pub struct Point {
    pub pos: i64,
    pub gene_name: String,
    pub transcript_id: String,
    pub strand: Strand,
    pub counter: HashMap<ReadAssign, i32>,
    pub exon_type: ExonType,
}

impl Point {
    fn new(
        pos: i64,
        gene_name: String,
        transcript_id: String,
        strand: Strand,
        exon_type: ExonType,
    ) -> Self {
        Point {
            pos,
            gene_name,
            transcript_id,
            strand,
            counter: HashMap::new(),
            exon_type,
        }
    }
}

pub trait InsideCounter {
    fn insert_in_counter(&mut self, item: ReadAssign) -> ();
}

impl InsideCounter for Point {
    fn insert_in_counter(self: &mut Self, item: ReadAssign) -> () {
        if let Some(val) = self.counter.get_mut(&item) {
            *val += 1;
        } else {
            self.counter.insert(item, 1);
        }
    }
}

impl PartialEq<i64> for Point {
    fn eq(&self, other: &i64) -> bool {
        self.pos == *other
    }
}

impl PartialEq<Point> for Point {
    fn eq(&self, other: &Point) -> bool {
        self.pos == other.pos
    }
}

impl PartialOrd<i64> for Point {
    fn partial_cmp(&self, other: &i64) -> Option<Ordering> {
        Some(self.pos.cmp(other))
    }
}

impl PartialOrd<Point> for Point {
    fn partial_cmp(&self, other: &Point) -> Option<Ordering> {
        Some(self.pos.cmp(&other.pos))
    }
}
impl Eq for Point {}

impl Ord for Point {
    fn cmp(&self, other: &Point) -> Ordering {
        self.pos.cmp(&other.pos)
    }
}

pub fn get_attr_id(attr: &str, toget: &str) -> Option<String> {
    let mut result: String; // = "".to_string();
    let x = attr.split(';').collect::<Vec<&str>>();
    for e in x {
        let spl = e.trim().split(' ').collect::<Vec<&str>>();
        if spl[0].trim() == toget {
            result = spl[1].trim().trim_matches('\"').to_string();
            return Some(result);
        }
    }
    None
}

pub fn read_gtf(file: &str) -> Result<HashMap<String, PointContainer>, Box<dyn Error>> {
    let f = File::open(file)?;
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
        this_line = line?;
        let spt = this_line.trim().split('\t').collect::<Vec<&str>>();
        if spt.len() < 8 {
            continue;
        }
        if spt[2] != "exon" {
            continue;
        }

        chr_ = spt[0].to_string();
        start = spt[3].parse::<i64>()?;
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
            .entry(chr_.clone())
            .or_insert_with(|| PointContainer::new());
        if let Some(val) = results.get_mut(&chr_) {
            if !been_seen.contains(&(transcript_id.to_string(), start, start_type)) {
                val.push(Point {
                    pos: start,
                    gene_name: gene_name.clone(), //get_gene_id(spt[8]),
                    transcript_id: transcript_id.clone(),
                    strand: strand,
                    counter: HashMap::new(),
                    exon_type: start_type,
                });
                been_seen.insert((transcript_id.to_string(), start, start_type));
            }
            if !been_seen.contains(&(transcript_id.to_string(), end, end_type)) {
                val.push(Point {
                    pos: end,
                    gene_name: gene_name.clone(), //get_gene_id(spt[8]),
                    transcript_id: transcript_id.clone(),
                    strand: strand,
                    counter: HashMap::new(),
                    exon_type: end_type,
                });
                been_seen.insert((transcript_id.to_string(), end, end_type));
            }
        }
    }
    for (_k, v) in results.iter_mut() {
        v.sort()
    }
    Ok(results)
}

pub fn my_bs(vec: &Vec<Point>, thr: &i64) -> Result<usize, usize> {
    let length = vec.len();
    let mut i = length / 2;
    let mut left = 0;
    let mut right = length - 1;

    // first check lower and upper bound
    if vec[length - 1] < *thr {
        return Err(length);
    }
    if vec[0] > *thr {
        return Err(0);
    }

    // so know we are free!
    //let mut it = 0;
    while left <= right {
        i = (left + right) / 2;
        if vec[i] == *thr {
            return Ok(i);
        }
        if vec[i] < *thr {
            left = i + 1;
        } else {
            right = i - 1;
        }
    }
    // while (i >= 1) && (i < length) && (vec[i] < *thr) {  ///
    //    i += 1
    // }
    // because we want the insert  ///
    return Err(i.max(0));
}

#[cfg(test)]
mod test_points {
    use clap::builder::Str;

    use super::*;

    #[test]
    fn lower_bound() {
        let mut cont = PointContainer::new();
        cont.push(Point::new(
            35,
            "ee".to_string(),
            "ee".to_string(),
            Strand::Plus,
            ExonType::Acceptor,
        ));
        cont.push(Point::new(
            35,
            "ee".to_string(),
            "ee".to_string(),
            Strand::Plus,
            ExonType::Acceptor,
        ));
        cont.push(Point::new(
            35,
            "ee".to_string(),
            "ee".to_string(),
            Strand::Plus,
            ExonType::Acceptor,
        ));
        cont.push(Point::new(
            45,
            "ee".to_string(),
            "ee".to_string(),
            Strand::Plus,
            ExonType::Acceptor,
        ));
        cont.push(Point::new(
            45,
            "ee".to_string(),
            "ee".to_string(),
            Strand::Plus,
            ExonType::Acceptor,
        ));
        cont.push(Point::new(
            45,
            "ee".to_string(),
            "ee".to_string(),
            Strand::Plus,
            ExonType::Acceptor,
        ));
        cont.push(Point::new(
            55,
            "ee".to_string(),
            "ee".to_string(),
            Strand::Plus,
            ExonType::Acceptor,
        ));
        cont.push(Point::new(
            65,
            "ee".to_string(),
            "ee".to_string(),
            Strand::Plus,
            ExonType::Acceptor,
        ));
        cont.push(Point::new(
            85,
            "ee".to_string(),
            "ee".to_string(),
            Strand::Plus,
            ExonType::Acceptor,
        ));
        cont.push(Point::new(
            105,
            "ee".to_string(),
            "ee".to_string(),
            Strand::Plus,
            ExonType::Acceptor,
        ));

        assert_eq!(cont.lower_bound(34), 0);
        assert_eq!(cont.lower_bound(40), 3);
        assert_eq!(cont.lower_bound(45), 3);
        assert_eq!(cont.lower_bound(35), 0);
        assert_eq!(cont.lower_bound(105), 9);

        //println!("{}", cont.lower_bound(110));
    }
    #[test]
    fn upper_bound() {
        let mut cont = PointContainer::new();
        cont.push(Point::new(
            35,
            "ee".to_string(),
            "ee".to_string(),
            Strand::Plus,
            ExonType::Acceptor,
        ));
        cont.push(Point::new(
            45,
            "ee".to_string(),
            "ee".to_string(),
            Strand::Plus,
            ExonType::Acceptor,
        ));
        cont.push(Point::new(
            45,
            "ee".to_string(),
            "ee".to_string(),
            Strand::Plus,
            ExonType::Acceptor,
        ));
        cont.push(Point::new(
            45,
            "ee".to_string(),
            "ee".to_string(),
            Strand::Plus,
            ExonType::Acceptor,
        ));
        cont.push(Point::new(
            55,
            "ee".to_string(),
            "ee".to_string(),
            Strand::Plus,
            ExonType::Acceptor,
        ));
        cont.push(Point::new(
            65,
            "ee".to_string(),
            "ee".to_string(),
            Strand::Plus,
            ExonType::Acceptor,
        ));
        cont.push(Point::new(
            85,
            "ee".to_string(),
            "ee".to_string(),
            Strand::Plus,
            ExonType::Acceptor,
        ));
        cont.push(Point::new(
            105,
            "ee".to_string(),
            "ee".to_string(),
            Strand::Plus,
            ExonType::Acceptor,
        ));

        println!("{:?}", my_bs(&cont.points, &115));
        assert_eq!(cont.upper_bound(34), 0);
        assert_eq!(cont.upper_bound(40), 1);
        assert_eq!(cont.upper_bound(45), 4);
        assert_eq!(cont.upper_bound(35), 1);
        assert_eq!(cont.upper_bound(105), 7);
        assert_eq!(cont.upper_bound(115), 7);
    }
}
