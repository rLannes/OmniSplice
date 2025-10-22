use CigarParser::cigar::Cigar;
use bio::bio_types::annot::spliced::Spliced;
use clap::builder::Str;
use lazy_static::lazy_static;
use log::{debug, error, info, trace, warn};
use regex::Regex;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::collections::HashSet;
use std::convert::TryFrom;
use std::convert::TryInto;
use std::fmt;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use std::str::FromStr;
use strand_specifier_lib::Strand;

#[derive(Debug)]
pub struct Intervall {
    pub start: i64,
    pub end: i64,
}

impl Intervall {
    pub fn from_zero_based(start: i64, end: i64) -> Self {
        Intervall { start, end }
    }

    pub fn from_one_based(start: i64, end: i64) -> Self {
        Intervall {
            start: start - 1,
            end: end,
        }
    }

    pub fn intervall(&self) -> (i64, i64) {
        (self.start, self.end)
    }

    pub fn intervall_zero(&self) -> (i64, i64) {
        (self.start + 1, self.end)
    }
}

#[derive(Debug)]
pub struct Exon {
    pub start: i64,
    pub end: i64,
    pub strand: Strand,
    pub contig: String,
}

pub fn parse_header(header: &str) -> HashMap<String, usize> {
    let mut map: HashMap<String, usize> = HashMap::new();
    header
        .split('\t')
        .collect::<Vec<&str>>()
        .iter()
        .enumerate()
        .for_each(|(i, v)| {
            map.insert(v.to_string(), i);
        });
    map
}

#[macro_export]
macro_rules! get_header {
    ($a:expr, $b:expr) => {
        *$b.get($a).expect(&format!("cannot find {} in header", $a))
    };
}

pub(crate) use get_header;

use crate::common::error::OmniError;

/*macro_rules! create_readfile {
    ($a:expr, $b:expr) => {
        *$b.get($a).expect(&format!("cannot find {} in header", $a ))
    }
}*/

/// Enum representing all the differents reads at a splicing junction
#[derive(clap::ValueEnum, Clone, Debug)]
pub enum ReadsToWrite {
    ReadThrough,
    ReadJunction,
    Unexpected,
    FailPosFilter,
    WrongStrand,
    FailQc,
    EmptyPileup,
    Skipped,
    SoftClipped,
    OverhangFail,
    Empty,
    All,
}

pub fn update_read_to_write_handle(
    read_out_handle: &mut ReadToWriteHandle,
    read_to_write: Vec<ReadsToWrite>,
    header_reads_handle: &[u8],
    output_file_prefix: &str,
) {
    for e in read_to_write {
        match e {
            ReadsToWrite::All => {
                read_out_handle.all = Some(BufWriter::new(
                    File::create_new(format!("{}{}", output_file_prefix, ".readAll"))
                        .unwrap_or_else(|_| panic!("readAll file should not exist.")),
                ));
                read_out_handle
                    .all
                    .as_mut()
                    .unwrap()
                    .write_all(header_reads_handle)
                    .expect("Unable to write file");
            }
            ReadsToWrite::ReadThrough => {
                read_out_handle.read_through = Some(BufWriter::new(
                    File::create_new(format!("{}{}", output_file_prefix, ".readThrough"))
                        .unwrap_or_else(|_| panic!("readThrough file should not exist.")),
                ));
                read_out_handle
                    .read_through
                    .as_mut()
                    .unwrap()
                    .write_all(header_reads_handle)
                    .expect("Unable to write file");
            }
            ReadsToWrite::ReadJunction => {
                read_out_handle.read_junction = Some(BufWriter::new(
                    File::create_new(format!("{}{}", output_file_prefix, ".readJunction"))
                        .unwrap_or_else(|_| panic!("readJunction file should not exist.")),
                ));
                read_out_handle
                    .read_junction
                    .as_mut()
                    .unwrap()
                    .write_all(header_reads_handle)
                    .expect("Unable to write file");
            }
            ReadsToWrite::Unexpected => {
                read_out_handle.unexpected = Some(BufWriter::new(
                    File::create_new(format!("{}{}", output_file_prefix, ".unexpected"))
                        .unwrap_or_else(|_| panic!("unexpected file should not exist.")),
                ));
                read_out_handle
                    .unexpected
                    .as_mut()
                    .unwrap()
                    .write_all(header_reads_handle)
                    .expect("Unable to write file");
            }
            ReadsToWrite::FailPosFilter => {
                read_out_handle.fail_pos_filter = Some(BufWriter::new(
                    File::create_new(format!("{}{}", output_file_prefix, ".FailPosFilter"))
                        .unwrap_or_else(|_| panic!("FailPosFilter file should not exist.")),
                ));
                read_out_handle
                    .fail_pos_filter
                    .as_mut()
                    .unwrap()
                    .write_all(header_reads_handle)
                    .expect("Unable to write file");
            }
            ReadsToWrite::WrongStrand => {
                read_out_handle.wrong_strand = Some(BufWriter::new(
                    File::create_new(format!("{}{}", output_file_prefix, ".WrongStrand"))
                        .unwrap_or_else(|_| panic!("WrongStrand file should not exist.")),
                ));
                read_out_handle
                    .wrong_strand
                    .as_mut()
                    .unwrap()
                    .write_all(header_reads_handle)
                    .expect("Unable to write file");
            }
            ReadsToWrite::FailQc => {
                read_out_handle.fail_qc = Some(BufWriter::new(
                    File::create_new(format!("{}{}", output_file_prefix, ".FailQC"))
                        .unwrap_or_else(|_| panic!("FailQC file should not exist.")),
                ));
                read_out_handle
                    .fail_qc
                    .as_mut()
                    .unwrap()
                    .write_all(header_reads_handle)
                    .expect("Unable to write file");
            }
            ReadsToWrite::EmptyPileup => {
                read_out_handle.empty_pileup = Some(BufWriter::new(
                    File::create_new(format!("{}{}", output_file_prefix, ".EmptyPileup"))
                        .unwrap_or_else(|_| panic!("EmptyPileup file should not exist.")),
                ));
                read_out_handle
                    .empty_pileup
                    .as_mut()
                    .unwrap()
                    .write_all(header_reads_handle)
                    .expect("Unable to write file");
            }
            ReadsToWrite::Skipped => {
                read_out_handle.skipped = Some(BufWriter::new(
                    File::create_new(format!("{}{}", output_file_prefix, ".Skipped"))
                        .unwrap_or_else(|_| panic!("Skipped file should not exist.")),
                ));
                read_out_handle
                    .skipped
                    .as_mut()
                    .unwrap()
                    .write_all(header_reads_handle)
                    .expect("Unable to write file");
            }
            ReadsToWrite::SoftClipped => {
                read_out_handle.soft_clipped = Some(BufWriter::new(
                    File::create_new(format!("{}{}", output_file_prefix, ".SoftClipped"))
                        .unwrap_or_else(|_| panic!("SoftClipped file should not exist.")),
                ));
                read_out_handle
                    .soft_clipped
                    .as_mut()
                    .unwrap()
                    .write_all(header_reads_handle)
                    .expect("Unable to write file");
            }
            ReadsToWrite::OverhangFail => {
                read_out_handle.overhang_fail = Some(BufWriter::new(
                    File::create_new(format!("{}{}", output_file_prefix, ".OverhangFail"))
                        .unwrap_or_else(|_| panic!("OverhangFail file should not exist.")),
                ));
                read_out_handle
                    .overhang_fail
                    .as_mut()
                    .unwrap()
                    .write_all(header_reads_handle)
                    .expect("Unable to write file");
            }
            ReadsToWrite::Empty => {
                read_out_handle.empty = Some(BufWriter::new(
                    File::create_new(format!("{}{}", output_file_prefix, ".Empty"))
                        .unwrap_or_else(|_| panic!("Empty file should not exist.")),
                ));
                read_out_handle
                    .empty
                    .as_mut()
                    .unwrap()
                    .write_all(header_reads_handle)
                    .expect("Unable to write file");
            }
        }
    }
}

#[derive(Clone, Debug, Copy, Eq, Hash, PartialEq)]
pub enum ExonType {
    Donnor,
    Acceptor,
}

impl fmt::Display for ExonType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            ExonType::Donnor => {
                write!(f, "Donnor")
            }
            ExonType::Acceptor => {
                write!(f, "Acceptor")
            }
        }
    }
}

impl From<&str> for ExonType {
    fn from(item: &str) -> Self {
        match item {
            "Donnor" => ExonType::Donnor,
            "Acceptor" => ExonType::Acceptor,
            &_ => {
                println!("{}", item);
                unreachable!()
            }
        }
    }
}

pub struct ReadToWriteHandle {
    pub read_through: Option<BufWriter<File>>,
    pub read_junction: Option<BufWriter<File>>,
    pub unexpected: Option<BufWriter<File>>,
    pub fail_pos_filter: Option<BufWriter<File>>,
    pub wrong_strand: Option<BufWriter<File>>,
    pub fail_qc: Option<BufWriter<File>>,
    pub empty_pileup: Option<BufWriter<File>>,
    pub skipped: Option<BufWriter<File>>,
    pub soft_clipped: Option<BufWriter<File>>,
    pub overhang_fail: Option<BufWriter<File>>,
    pub empty: Option<BufWriter<File>>,
    pub all: Option<BufWriter<File>>,
}
impl ReadToWriteHandle {
    pub fn new() -> Self {
        ReadToWriteHandle {
            all: None,
            read_through: None,
            read_junction: None,
            unexpected: None,
            fail_pos_filter: None,
            wrong_strand: None,
            fail_qc: None,
            empty_pileup: None,
            skipped: None,
            soft_clipped: None,
            overhang_fail: None,
            empty: None,
        }
    }
    pub fn flush(&mut self) -> Result<(), std::io::Error> {
        match self.read_through {
            Some(ref mut handle) => handle.flush(),
            _ => Ok(()),
        };
        match self.read_junction {
            Some(ref mut handle) => handle.flush(),
            _ => Ok(()),
        };
        match self.unexpected {
            Some(ref mut handle) => handle.flush(),
            _ => Ok(()),
        };
        match self.fail_pos_filter {
            Some(ref mut handle) => handle.flush(),
            _ => Ok(()),
        };
        match self.wrong_strand {
            Some(ref mut handle) => handle.flush(),
            _ => Ok(()),
        };
        match self.fail_qc {
            Some(ref mut handle) => handle.flush(),
            _ => Ok(()),
        };
        match self.empty_pileup {
            Some(ref mut handle) => handle.flush(),
            _ => Ok(()),
        };
        match self.skipped {
            Some(ref mut handle) => handle.flush(),
            _ => Ok(()),
        };
        match self.soft_clipped {
            Some(ref mut handle) => handle.flush(),
            _ => Ok(()),
        };
        match self.overhang_fail {
            Some(ref mut handle) => handle.flush(),
            _ => Ok(()),
        };
        match self.empty {
            Some(ref mut handle) => handle.flush(),
            _ => Ok(()),
        };
        match self.all {
            Some(ref mut handle) => handle.flush(),
            _ => Ok(()),
        };
        Ok(())
    }
}

#[derive(Clone, Debug, Copy, Eq, Hash, PartialEq)]
pub enum SplicingEvent {
    Spliced,
    Unspliced,
    Clipped,
    ExonOther,
    WrongStrand,
    Skipped,
    SkippedUnrelated,
    Isoform,
}

impl SplicingEvent {
    pub fn merge(
        left: Option<SplicingEvent>,
        right: Option<SplicingEvent>,
    ) -> Option<SplicingEvent> {
        match (left, right) {
            (None, Some(x)) | (None, Some(x)) => Some(x),
            (Some(SplicingEvent::Spliced), Some(SplicingEvent::Spliced)) => {
                Some(SplicingEvent::Spliced)
            }
            (Some(SplicingEvent::Isoform), _) | (_, Some(SplicingEvent::Isoform)) => {
                Some(SplicingEvent::Isoform)
            }
            (Some(SplicingEvent::ExonOther), _) | (_, Some(SplicingEvent::ExonOther)) => {
                Some(SplicingEvent::ExonOther)
            }
            (Some(SplicingEvent::Unspliced), _) | (_, Some(SplicingEvent::Unspliced)) => {
                Some(SplicingEvent::Unspliced)
            }
            (Some(SplicingEvent::Skipped), _) | (_, Some(SplicingEvent::Skipped)) => {
                Some(SplicingEvent::Skipped)
            }
            (Some(SplicingEvent::Clipped), _) | (_, Some(SplicingEvent::Clipped)) => {
                Some(SplicingEvent::Clipped)
            }
            (Some(SplicingEvent::SkippedUnrelated), _)
            | (_, Some(SplicingEvent::SkippedUnrelated)) => Some(SplicingEvent::SkippedUnrelated),
            (Some(SplicingEvent::WrongStrand), Some(SplicingEvent::WrongStrand)) => {
                Some(SplicingEvent::WrongStrand)
            }
            (None, None) => None,
            (None, Some(SplicingEvent::WrongStrand) | Some(SplicingEvent::Spliced))
            | (Some(SplicingEvent::WrongStrand) | Some(SplicingEvent::Spliced), None) => {
                error!("Unreachable code!");
                unreachable!()
            }
            (_, _) => {
                error!("Unreachable code!");
                unreachable!()
            }
        }
    }

    pub fn from_read_assign(
        item: Option<ReadAssign>,
        contig: String,
        valid_junction: &HashSet<(i64, i64)>,
        feature_start: Option<i64>,
        feature_end: Option<i64>,
        feature_strand: &Strand,
        read_strand: &Strand,
    ) -> Option<Self> {
        match item {
            None => None,
            Some(ReadAssign::WrongStrand) => Some(SplicingEvent::WrongStrand),
            Some(ReadAssign::ReadThrough) => Some(SplicingEvent::Unspliced),
            Some(ReadAssign::SoftClipped) => Some(SplicingEvent::Clipped),
            Some(ReadAssign::Skipped(n, m)) => {
                if valid_junction.contains(&(n, m)) || valid_junction.contains(&(m, n)) {
                    Some(SplicingEvent::Isoform)
                } else {
                    // if fe
                    Some(SplicingEvent::Skipped)
                }
            }
            Some(ReadAssign::SkippedUnrelated(n, m)) => {
                if valid_junction.contains(&(n, m)) || valid_junction.contains(&(m, n)) {
                    Some(SplicingEvent::Isoform)
                } else {
                    // if fe
                    Some(SplicingEvent::SkippedUnrelated)
                }
            }
            Some(ReadAssign::ReadJunction(n, m)) => {
                if (feature_start.is_some() && feature_end.is_some())
                    && (((n == feature_end.unwrap()) & (m == feature_start.unwrap()))
                        || ((n == feature_start.unwrap()) & (m == feature_end.unwrap())))
                {
                    Some(SplicingEvent::Spliced)
                } else if valid_junction.contains(&(n, m)) || valid_junction.contains(&(m, n)) {
                    Some(SplicingEvent::Isoform)
                } else {
                    Some(SplicingEvent::ExonOther)
                }
            }

            _ => None,
        }
    }
}

#[derive(Clone, Debug, Copy, Eq, Hash, PartialEq)]
pub enum ReadAssign {
    ReadThrough,
    ReadJunction(i64, i64),
    Unexpected,
    FailPosFilter,
    WrongStrand,
    FailQc,
    EmptyPileup,
    Skipped(i64, i64),
    SkippedUnrelated(i64, i64),
    SoftClipped,
    OverhangFail,
    Empty,
}

impl fmt::Display for ReadAssign {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            ReadAssign::ReadThrough => {
                write!(f, "Unspliced")
            }
            ReadAssign::ReadJunction(n, m) => {
                write!(f, "{}", format!("ReadJunction({}, {})", n, m))
            }
            ReadAssign::Skipped(n, m) => {
                write!(f, "{}", format!("Skipped({}, {})", n, m))
            }
            ReadAssign::SkippedUnrelated(n, m) => {
                write!(f, "{}", format!("SkippedUnrelated({}, {})", n, m))
            }
            ReadAssign::Unexpected => {
                write!(f, "Unexpected")
            }
            ReadAssign::Empty => {
                write!(f, "Empty")
            }
            ReadAssign::FailPosFilter => {
                write!(f, "FailPosFilter")
            }
            /*             ReadAssign::DoesNotMatchP1P => {
                write!(f, "DoesNotMatchP1P")
            } */
            ReadAssign::WrongStrand => {
                write!(f, "WrongStrand")
            }
            ReadAssign::FailQc => {
                write!(f, "FailQc")
            }
            ReadAssign::EmptyPileup => {
                write!(f, "EmptyPileup")
            }

            ReadAssign::SoftClipped => {
                write!(f, "SoftClipped")
            }
            ReadAssign::OverhangFail => {
                write!(f, "OverhangFail")
            }
        }
    }
}

lazy_static! {
    static ref reg_junction: Regex =
        Regex::new(r"(\d+),\s+(\d+)").expect("Failed to compile junction regexp");
}
//L/azyLock<Regex> =
//LazyLock::new(||

impl From<&str> for ReadAssign {
    fn from(item: &str) -> Self {
        match item {
            "Empty" => ReadAssign::Empty,
            "Unspliced" => ReadAssign::ReadThrough,
            "Unexpected" => ReadAssign::Unexpected,
            "FailPosFilter" => ReadAssign::FailPosFilter,
            //"DoesNotMatchP1P" => ReadAssign::DoesNotMatchP1P,
            "WrongStrand" => ReadAssign::WrongStrand,
            "FailQc" => ReadAssign::FailQc,
            "EmptyPileup" => ReadAssign::EmptyPileup,
            //"Skipped" => ReadAssign::Skipped,
            //"SkippedUnrelated" => ReadAssign::SkippedUnrelated,
            "SoftClipped" => ReadAssign::SoftClipped,
            "OverhangFail" => ReadAssign::OverhangFail,
            "empty" => ReadAssign::Empty,
            s if s.starts_with("SkippedUnrelated") => {
                let v = reg_junction.captures(&s).unwrap();
                ReadAssign::SkippedUnrelated(
                    v.get(1)
                        .unwrap()
                        .as_str()
                        .parse::<i64>()
                        .expect("failed to parse junction"), // Should never failed
                    v.get(2)
                        .unwrap()
                        .as_str()
                        .parse::<i64>()
                        .expect("failed to parse junction"),
                )
            }
            s if s.starts_with("Skipped") => {
                let v = reg_junction.captures(&s).unwrap();
                ReadAssign::Skipped(
                    v.get(1)
                        .unwrap()
                        .as_str()
                        .parse::<i64>()
                        .expect("failed to parse junction"),
                    v.get(2)
                        .unwrap()
                        .as_str()
                        .parse::<i64>()
                        .expect("failed to parse junction"),
                )
            }
            s if s.starts_with("ReadJunction") => {
                let v = reg_junction.captures(&s).unwrap();
                ReadAssign::ReadJunction(
                    v.get(1)
                        .unwrap()
                        .as_str()
                        .parse::<i64>()
                        .expect("failed to parse junction"),
                    v.get(2)
                        .unwrap()
                        .as_str()
                        .parse::<i64>()
                        .expect("failed to parse junction"),
                )
            }
            _ => {
                println!("{}", item);
                unreachable!();
            }
        }
    }
}

pub fn out_of_range(feature: i64, aln_start: i64, aln_end: i64, overhang: i64) -> bool {
    if !((aln_start <= feature) & (aln_end >= feature)) {
        return true;
    }
    false
}

/// Return Wrong Strand or None if Strand is NA return None.
pub fn test_wrong_strand(read_strand: &Strand, feature_strand: &Strand) -> Option<ReadAssign> {
    match read_strand {
        Strand::NA => (None),
        read => {
            if (*read_strand != *feature_strand) {
                return Some(ReadAssign::WrongStrand);
            } else {
                return None;
            }
        }
    }
}

fn test_skipped(
    cigar: &Cigar,
    aln_start: i64,
    feature_pos: i64,
    feature_bis: Option<i64>,
) -> Option<ReadAssign> {
    let junction = cigar.get_skipped_pos_on_ref(aln_start);
    if let Some(y) = junction {
        if let Some(i) = y
            .iter()
            .enumerate()
            .step_by(2)
            .find(|&(ref i, &x)| (x < feature_pos) & (y[i + 1] > feature_pos))
        {
            if feature_bis.is_some() {
                let (start, end) = match feature_pos.cmp(&feature_bis.unwrap()) {
                    // cannot failed tested just before
                    Ordering::Less | Ordering::Equal => (feature_pos, feature_bis.unwrap()), // cannot failed tested just before
                    Ordering::Greater => (feature_bis.unwrap(), feature_pos), // cannot failed tested just before
                };
                if cigar.does_it_overlap_an_intervall(aln_start, start, end) {
                    return Some(ReadAssign::Skipped(*i.1, y[i.0 + 1]));
                }
            }
            return Some(ReadAssign::SkippedUnrelated(*i.1, y[i.0 + 1]));
        }
    }
    None
}

pub fn read_toassign(
    feature_strand: Strand,
    feature_pos: Option<i64>,
    feature_pos_bis: Option<i64>,
    feature_exontype: Option<ExonType>,
    aln_start: i64,
    aln_end: i64,
    cigar: &Cigar,
    //flag: &u16,
    read_strand: &Strand,
    overhang: i64,
) -> Result<Option<ReadAssign>, OmniError> {
    if feature_pos.is_none() {
        return Ok(None);
    }

    let feature_pos = feature_pos.ok_or(OmniError::Expect(
        "Expected a value for featurePos got None".to_string(),
    ))?;
    let feature_exontype = feature_exontype.ok_or(OmniError::Expect(
        "Expected a value for feature_exontype got None".to_string(),
    ))?;
    /// sanity check
    if out_of_range(feature_pos, aln_start, aln_end, overhang) {
        return Ok(None);
    }

    match test_wrong_strand(&read_strand, &feature_strand) {
        Some(x) => return Ok(Some(x)),
        _ => (),
    };

    match test_skipped(&cigar, aln_start, feature_pos, feature_pos_bis) {
        Some(x) => return Ok(Some(x)),
        _ => (),
    }

    match (feature_strand, feature_exontype) {
        (Strand::Plus, ExonType::Donnor) | (Strand::Minus, ExonType::Acceptor) => {
            if !cigar.does_it_match_an_intervall(aln_start, feature_pos - overhang, feature_pos) {
                return Ok(Some(ReadAssign::OverhangFail));
            }

            if cigar.does_it_match_an_intervall(
                aln_start,
                feature_pos - overhang,
                feature_pos + overhang,
            ) {
                return Ok(Some(ReadAssign::ReadThrough));
            }

            if cigar.soft_clipped_end(&Strand::Plus, 10) && aln_end == feature_pos {
                return Ok(Some(ReadAssign::SoftClipped));
            }
        }
        (Strand::Plus, ExonType::Acceptor) | (Strand::Minus, ExonType::Donnor) => {
            if !cigar.does_it_match_an_intervall(aln_start, feature_pos, feature_pos + overhang) {
                return Ok(Some(ReadAssign::OverhangFail));
            }

            //if cigar.does_it_match_an_intervall(&aln_start, feature_pos - 1, feature_pos + overhang)
            if cigar.does_it_match_an_intervall(
                aln_start,
                feature_pos - overhang,
                feature_pos + overhang,
            ) {
                return Ok(Some(ReadAssign::ReadThrough));
            }

            if cigar.soft_clipped_end(&Strand::Minus, 10) && aln_start == feature_pos {
                return Ok(Some(ReadAssign::SoftClipped));
            }
        }
        //feature strand!s
        (Strand::NA, _) => {
            return Ok(None);
        }
    };

    match cigar.get_skipped_pos_on_ref(aln_start) {
        Some(j) => {
            match (
                j.iter().position(|&x| x == feature_pos),
                feature_strand,
                feature_exontype,
            ) {
                // TODO add overhang check!
                (Some(p), Strand::Plus, ExonType::Donnor)
                | (Some(p), Strand::Minus, ExonType::Acceptor) => {
                    if !((cigar.does_it_match_an_intervall(aln_start, j[p] - overhang, j[p]))
                        & (cigar.does_it_match_an_intervall(
                            aln_start,
                            j[p + 1],
                            j[p + 1] + overhang,
                        )))
                    {
                        return Ok(Some(ReadAssign::OverhangFail));
                    }
                    return Ok(Some(ReadAssign::ReadJunction(j[p], j[p + 1])));
                }

                (Some(p), Strand::Plus, ExonType::Acceptor)
                | (Some(p), Strand::Minus, ExonType::Donnor) => {
                    if !((cigar.does_it_match_an_intervall(
                        aln_start,
                        j[p - 1] - overhang,
                        j[p - 1],
                    )) & (cigar.does_it_match_an_intervall(aln_start, j[p], j[p] + overhang)))
                    {
                        return Ok(Some(ReadAssign::OverhangFail));
                    }
                    return Ok(Some(ReadAssign::ReadJunction(j[p - 1], j[p])));
                }
                (_, _, _) => {
                    if let Some(i) = j
                        .iter()
                        .enumerate()
                        .step_by(2)
                        .find(|&(ref i, &x)| (x < feature_pos) & (j[i + 1] > feature_pos))
                    {
                        //println!("Skipped: {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?}",\
                        // j, aln_start ,aln_end, cigar, read_strand, feature_strand, feature_pos, feature_exontype, overhang);
                        return Ok(Some(ReadAssign::Skipped(*i.1, j[i.0 + 1])));
                    } else {
                        // TODO this should be in the overhang?
                        // right now this is overhang 0 edge case.
                        if (aln_start == feature_pos) | (aln_end == feature_pos) {
                            return Ok(None);
                        }
                        //println!("Unexp: {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?}",\
                        // j, aln_start ,aln_end, cigar, read_strand, feature_strand, feature_pos, feature_exontype, overhang);
                        warn!(
                            "read unexpected alignment -> cigar:{:?} feature_strand:{:?}, read_strand: {:?}, feature_pos:{:?}, aln_start:{:?}, aln_end:{:?}",
                            cigar, feature_strand, read_strand, feature_pos, aln_start, aln_end
                        );
                        return Ok(Some(ReadAssign::Unexpected));
                    }
                }
            }
        }
        None => {
            warn!(
                "read unexpected alignment -> cigar:{:?} feature_strand:{:?}, read_strand: {:?}, feature_pos:{:?}, aln_start:{:?}, aln_end:{:?}",
                cigar, feature_strand, read_strand, feature_pos, aln_start, aln_end
            );
            return Ok(Some(ReadAssign::Unexpected));
        }
    }
}

/*
feature_strand: Strand,
feature_pos: i64,
feature_exontype: ExonType,
aln_start: i64,
aln_end: i64,
cigar: &Cigar,
//flag: &u16,
read_strand: &Strand,
overhang: i64 */

#[cfg(test)]
mod tests_it {
    use super::*;
    use crate::common::gtf_::{get_junction_from_gtf, gtf_to_hashmap};
    use crate::common::it_intron::TreeDataIntron;

    #[test]
    fn parse_strand_1() {
        // test_wrong_strand(read_strand: &Strand, feature_strand: &Strand)  -> Option<ReadAssign>
        let read_strand = Strand::Plus;
        let feature_strand = Strand::Plus;
        assert_eq!(test_wrong_strand(&read_strand, &feature_strand), None);
        let read_strand = Strand::Minus;
        let feature_strand = Strand::Minus;
        assert_eq!(test_wrong_strand(&read_strand, &feature_strand), None);

        let read_strand = Strand::Plus;
        let feature_strand = Strand::Minus;
        assert_eq!(
            test_wrong_strand(&read_strand, &feature_strand),
            Some(ReadAssign::WrongStrand)
        );
        let read_strand = Strand::Minus;
        let feature_strand = Strand::Plus;
        assert_eq!(
            test_wrong_strand(&read_strand, &feature_strand),
            Some(ReadAssign::WrongStrand)
        );
    }

    /*

               if cigar.soft_clipped_end(&Strand::Plus, 10) && aln_end == feature_pos {

               return Some(ReadAssign::SoftClipped);
           }

    */

    #[test]
    fn parse_clipped_1() {
        let cigar = Cigar::from_str("50M30S").unwrap();
        let aln_start = 1;
        let feature_pos = 51;
        let aln_end = cigar.get_end_of_aln(aln_start);
        assert_eq!(
            (cigar.soft_clipped_end(&Strand::Plus, 10) && aln_end == feature_pos),
            true
        );
        let cigar = Cigar::from_str("50S30M").unwrap();
        let aln_start = 51;
        assert_eq!(
            (cigar.soft_clipped_end(&Strand::Minus, 10) && aln_start == feature_pos),
            true
        );
    }

    #[test]
    fn test_1() {
        let aln_start = 21589327;
        let feature_pos = 21589347;
        let overhang = 5;
        let cigar = Cigar::from_str("22M264N78M").unwrap();

        assert_eq!(
            cigar.does_it_match_an_intervall(aln_start, feature_pos - overhang, feature_pos),
            true
        );
    }
    #[test]
    fn test_2() {
        let aln_start = 21589345;
        let feature_pos = 21589347;
        let overhang = 5;
        let cigar = Cigar::from_str("22M264N78M").unwrap();

        assert_eq!(
            cigar.does_it_match_an_intervall(aln_start, feature_pos - overhang, feature_pos),
            false
        );
    }

    #[test]
    fn test_3() {
        let aln_start = 21597482;
        let feature_pos = 21611098;
        let overhang = 5;
        let cigar = Cigar::from_str("89M13524N11M").unwrap();
        println!("{:?}", cigar.get_skipped_pos_on_ref(aln_start));
        assert_eq!(true, true);
    }

    #[test]
    fn test_read_toassign() {
        let aln_start = 21597482;
        let feature_pos = 21611098;
        let overhang = 5;

        let x = read_toassign(
            Strand::Minus,
            Some(5334258),
            Some(5334402),
            Some(ExonType::Donnor),
            5333978,
            5334283,
            &Cigar::from("225M55N26M"),
            &Strand::Minus,
            1,
        );

        info!("read assign {:?}", x);
        println!("read assign {:?}", x);

        let y = Cigar::from("225M55N26M").get_skipped_pos_on_ref(5333978);
        println!(" {:?}", y);

        assert_eq!(true, true);
    }
} //21589349, 21589613
