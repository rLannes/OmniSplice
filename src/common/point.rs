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

/* pub struct PointContainerMutableIterator<'a, T> where T: Clone + Eq + Hash + PartialEq {
    data: &'a mut [Point<T>],
    index: usize,
} */

/* impl<'a, T> IntoIterator for &'a mut PointContainer<T> where T: Clone + Eq + Hash + PartialEq  {
    type Item = &'a mut Point<T>;
    type IntoIter = PointContainerMutableIterator<'a, T>;
    fn into_iter(self) -> Self::IntoIter {
        PointContainerMutableIterator {
            data: &mut self.points,
            index: 0 }
    }
} */
/* 
impl<'a, T> Iterator for PointContainerMutableIterator<'a, T> where T: Clone + Eq + Hash + PartialEq  {
    // We can refer to this type using Self::Item
    type Item = &'a mut Point<T>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.index < self.data.len() {

            let item: *mut Point<T> = &mut self.data[self.index] as *mut Point<T>;
            self.index += 1;
            unsafe { Some(&mut *item) }

        } else {
            None
        }
    }
} */
/* 
impl<T> PointContainer<T> where T: Clone + Eq + Hash + PartialEq {
    pub fn iter_mut(&mut self) -> PointContainerMutableIterator<'_, T> {
        self.into_iter()
    }
} */

impl<'a> Iterator for PointContainerIterator<'a>   {
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

/*impl <'a, F> for PointContainerIterator<'a> {
    pub fn sort_by_key<K, F>(&mut self, mut f: F)
    where
        F: FnMut(&T) -> K,
        K: Ord,
        {
            stable_sort(self, |a, b| f(a).lt(&f(b)));
        }
}*/
/*
where
    F: FnMut(&T, &T) -> Ordering{

    pub fn sort_unstable_by_key(&mut self, compare: F){
        self.my_struct.points.cmp(F);
    }
}*/



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

        if (lower_bound == upper_bound) | (upper_bound <= 0) | (upper_bound >= self.points.len()) {
            return Err(());
        }

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
                //flag,
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
                        if read_assign == ReadAssign::SoftClipped{
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
                    else{
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
                while (*i < length - 1) && (self.points[*i + 1] <= thr) {
                    *i += 1;
                }
                *i
            }
            Err(i) => i,
        };
        j
    }

    pub fn push(self: &mut Self, point: Point) {
        self.points.push(point)
    }

    pub fn sort(self: &mut Self) {
        self.points.sort();
        //self.points.sort_by_key(|x| x.pos);
    }

    /* pub fn parse_reads(
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

        if (lower_bound == upper_bound) | (upper_bound <= 0) | (upper_bound >= self.points.len()) {
            return Err(());
        }

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
                //flag,
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
                        if read_assign == ReadAssign::SoftClipped{
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
                    else{
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
    } */
}

#[derive(Debug, Clone)]
pub struct Point{
    pub pos: i64,
    pub gene_name: String,
    pub transcript_id: String,
    pub strand: Strand,
    pub counter: HashMap<ReadAssign, i32>,
    pub exon_type: ExonType,
}


pub trait InsideCounter{
    fn insert_in_counter(&mut self, item: ReadAssign) -> ();
}

impl InsideCounter for Point{
    fn insert_in_counter(self: &mut Self, item: ReadAssign) -> () {
        if let Some(val) = self.counter.get_mut(&item) {
            *val += 1;
        } else {
            self.counter.insert(item, 1);
        }
    }
}

/* impl<T> Point<T> where T: Clone + Eq + Hash + PartialEq{

} */

impl PartialEq<i64> for Point  {
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
impl Eq for Point  {}

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
    let mut this_line: String;//::new();

    let mut chr_: String; // = "".to_string();
    let mut start: i64;//; = 0;
    let mut end: i64;// = 0;
    let mut strand: Strand;

    let mut gene_name: String; // = "".to_string();
    let mut transcript_id: String;// = "".to_string();

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
        } else if  let Some((gene_tmp, tr_tmp)) = get_attr_id(spt[8], "gene_name")
        .zip(get_attr_id(spt[8], "transcript_id")) {
            gene_name = gene_tmp;
            transcript_id = tr_tmp;
        }
        else{
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

    while (i >= 1) && (i < length) && (vec[i] <= *thr) {
        i += 1
    }
    // because we want the insert

    return Err(i.max(0));
}


/*     fn insert_in_counter(self: &mut Self, key: ReadAssign) -> () {
        if let Some(val) = self.counter.get_mut(&key) {
            *val += 1;
        } else {
            self.counter.insert(key, 1);
        }
    }
 */
   /*  pub fn parse_read(
        self: &mut Self,
        aln_start: i64,
        aln_end: i64,
        cigar: &Cigar,
        flag: &u16,
        read_strand: &Strand,
        overhang: i64,
    ) -> () {
        if *read_strand != self.strand {
            self.insert_in_counter(ReadAssign::WrongStrand);
            return ();
        }

        let junction = cigar.get_skipped_pos_on_ref(&aln_start);
        if let Some(y) = junction {
            if (y
                .iter()
                .enumerate()
                .step_by(2)
                .any(|(i, &x)| (x < self.pos) & (y[i + 1] > self.pos)))
            {
                self.insert_in_counter(ReadAssign::Skipped);
                return ();
            }
        }

        match (self.strand, self.exon_type) {
            (Strand::Plus, ExonType::Donnor) | (Strand::Minus, ExonType::Acceptor) => {
                if !cigar.does_it_match_an_intervall(&aln_start, self.pos - overhang, self.pos) {
                    self.insert_in_counter(ReadAssign::OverhangFail);
                    return ();
                }

                if cigar.does_it_match_an_intervall(&aln_start, self.pos - overhang, self.pos + 1) {
                    //overhang) {
                    self.insert_in_counter(ReadAssign::ReadThrough);
                    return ();
                }

                if cigar.soft_clipped_end(&Strand::Plus, 10) {
                    self.insert_in_counter(ReadAssign::SoftClipped);
                    return ();
                }
            }
            (Strand::Plus, ExonType::Acceptor) | (Strand::Minus, ExonType::Donnor) => {
                if !cigar.does_it_match_an_intervall(&aln_start, self.pos, self.pos + overhang) {
                    //println!("{:?} {} {} {:?}", self, aln_start, aln_end, cigar);
                    self.insert_in_counter(ReadAssign::OverhangFail);
                    return ();
                }

                if cigar.does_it_match_an_intervall(&aln_start, self.pos - 1, self.pos + overhang) {
                    //overhang) {
                    self.insert_in_counter(ReadAssign::ReadThrough);
                    return ();
                }

                if cigar.soft_clipped_end(&Strand::Minus, 10) {
                    self.insert_in_counter(ReadAssign::SoftClipped);
                    return ();
                }
            }
            (Strand::NA, _) => {
                return ();
            }
        };

        if !((aln_start < self.pos) & (aln_end > self.pos)) {
            self.insert_in_counter(ReadAssign::FailPosFilter);
            return ();
        }

        match cigar.get_skipped_pos_on_ref(&aln_start) {
            Some(j) => {
                match (
                    j.iter().position(|&x| x == self.pos),
                    self.strand,
                    self.exon_type,
                ) {
                    (Some(p), Strand::Plus, ExonType::Donnor) => {
                        self.insert_in_counter(ReadAssign::ReadJunction(j[p], j[p + 1]));
                        return ();
                    }
                    (Some(p), Strand::Plus, ExonType::Acceptor) => {
                        self.insert_in_counter(ReadAssign::ReadJunction(j[p - 1], j[p]));
                        return ();
                    }
                    (Some(p), Strand::Minus, ExonType::Donnor) => {
                        self.insert_in_counter(ReadAssign::ReadJunction(j[p - 1], j[p]));
                        return ();
                    }
                    (Some(p), Strand::Minus, ExonType::Acceptor) => {
                        self.insert_in_counter(ReadAssign::ReadJunction(j[p], j[p + 1]));
                        return ();
                    }

                    (_, _, _) => {
                        //TODO add skipped

                        if (j
                            .iter()
                            .enumerate()
                            .step_by(2)
                            .any(|(i, &x)| (x < self.pos) & (j[i + 1] > self.pos)))
                        {
                            self.insert_in_counter(ReadAssign::Skipped);
                            return ();
                        } else {
                            //println!("Unexp: {:?} {:?} {:?} {:?} {:?}", j, aln_start ,aln_end, cigar, self);
                            self.insert_in_counter(ReadAssign::Unexpected);
                            return ();
                        }
                    }
                }
            }
            None => {
                self.insert_in_counter(ReadAssign::Unexpected);
                return ();
            }
        }
    }

    */

    /*fn get_gene_id(attr: &str) -> String {
    let mut result = "".to_string();
    let x = attr.split(';').collect::<Vec<&str>>();
    for e in x {
        let spl = e.trim().split(' ').collect::<Vec<&str>>();
        if spl[0].trim() == "gene_id" {
            result = spl[1].trim().trim_matches('\"').to_string()
        }
    }
    result*/