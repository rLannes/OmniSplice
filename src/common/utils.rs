use lazy_static::lazy_static;
use regex::Regex;
//use std::collections::HashMap;
use std::fmt;
//use std::fs::File;
//use std::io::BufRead;
//use std::io::Write;
//use std::io::{BufReader, BufWriter};
//use std::str::FromStr;
use strand_specifier_lib::Strand;
use CigarParser::cigar::Cigar;

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
            &_ => unreachable!(),
        }
    }
}

#[derive(Clone, Debug, Copy, Eq, Hash, PartialEq)]
pub enum ReadAssign {
    ReadThrough,
    ReadJunction(i64, i64),
    Unexpected,
    FailPosFilter,
    //DoesNotMatchP1P,
    WrongStrand,
    FailQc,
    EmptyPileup,
    Skipped(i64, i64),
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
            "SoftClipped" => ReadAssign::SoftClipped,
            "OverhangFail" => ReadAssign::OverhangFail,
            "empty" => ReadAssign::Empty,
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

pub fn read_toassign(
    feature_strand: Strand,
    feature_pos: i64,
    feature_exontype: ExonType,
    aln_start: i64,
    aln_end: i64,
    cigar: &Cigar,
    //flag: &u16,
    read_strand: &Strand,
    overhang: i64,
) -> Option<ReadAssign> {
    // TODO Should I report it? no...
    if !((aln_start <= feature_pos) & (aln_end >= feature_pos)) {
        return None; //return Some(ReadAssign::FailPosFilter);
    }

    if *read_strand != feature_strand {
        return Some(ReadAssign::WrongStrand);
    }

    let junction = cigar.get_skipped_pos_on_ref(&aln_start);
    if let Some(y) = junction {
        if let Some(i) = y
            .iter()
            .enumerate()
            .step_by(2)
            .find(|(i, &x)| (x < feature_pos) & (y[i + 1] > feature_pos))
        {
            return Some(ReadAssign::Skipped(*i.1, y[i.0 + 1]));
        }
    }

    match (feature_strand, feature_exontype) {
        (Strand::Plus, ExonType::Donnor) | (Strand::Minus, ExonType::Acceptor) => {
            if !cigar.does_it_match_an_intervall(&aln_start, feature_pos - overhang, feature_pos) {
                return Some(ReadAssign::OverhangFail);
            }

            if cigar.does_it_match_an_intervall(&aln_start, feature_pos - overhang, feature_pos + 1)
            {
                //overhang) {
                return Some(ReadAssign::ReadThrough);
            }

            if cigar.soft_clipped_end(&Strand::Plus, 10) && aln_end == feature_pos {
                return Some(ReadAssign::SoftClipped);
            }
        }
        (Strand::Plus, ExonType::Acceptor) | (Strand::Minus, ExonType::Donnor) => {
            if !cigar.does_it_match_an_intervall(&aln_start, feature_pos, feature_pos + overhang) {
                //println!("{:?} {} {} {:?}", self, aln_start, aln_end, cigar);
                return Some(ReadAssign::OverhangFail);
            }

            if cigar.does_it_match_an_intervall(&aln_start, feature_pos - 1, feature_pos + overhang)
            {
                //overhang) {
                return Some(ReadAssign::ReadThrough);
            }

            if cigar.soft_clipped_end(&Strand::Minus, 10) && aln_start == feature_pos {
                return Some(ReadAssign::SoftClipped);
            }
        }
        (Strand::NA, _) => {
            return None;
        }
    };

    match cigar.get_skipped_pos_on_ref(&aln_start) {
        Some(j) => {
            match (
                j.iter().position(|&x| x == feature_pos),
                feature_strand,
                feature_exontype,
            ) {
                (Some(p), Strand::Plus, ExonType::Donnor) => {
                    return Some(ReadAssign::ReadJunction(j[p], j[p + 1]));
                }
                (Some(p), Strand::Plus, ExonType::Acceptor) => {
                    return Some(ReadAssign::ReadJunction(j[p - 1], j[p]));
                }
                (Some(p), Strand::Minus, ExonType::Donnor) => {
                    return Some(ReadAssign::ReadJunction(j[p - 1], j[p]));
                }
                (Some(p), Strand::Minus, ExonType::Acceptor) => {
                    return Some(ReadAssign::ReadJunction(j[p], j[p + 1]));
                }

                (_, _, _) => {
                    //TODO add skipped

                    if let Some(i) = j
                        .iter()
                        .enumerate()
                        .step_by(2)
                        .find(|(i, &x)| (x < feature_pos) & (j[i + 1] > feature_pos))
                    {
                        return Some(ReadAssign::Skipped(*i.1, j[i.0 + 1]));
                    } else {
                        //println!("Unexp: {:?} {:?} {:?} {:?} {:?}", j, aln_start ,aln_end, cigar, self);
                        return Some(ReadAssign::Unexpected);
                    }
                }
            }
        }
        None => {
            return Some(ReadAssign::Unexpected);
        }
    }
}
