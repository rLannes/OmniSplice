/// This is a self contain file that parse the table file and output an splicing effiency table.
/// i am unsure if I keep it separate or integrate it to the main software or both
/// it is to replicate and more the exact formual of spilcing e

use clap::Parser;

use std::collections::HashMap;
use std::ffi::OsStr;
use std::fmt::format;
use std::fs::File;
use std::hash::Hash;
use std::io::prelude::*;
use std::io::BufWriter;
use std::path::{Path, PathBuf};

//
// we want to support:
// - either Donnor only Acceptor only or both
// - we want to support the combination of different fields. "unspliced" + "exon_other" ...
//

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
////
// either counted once or sum
// spliced is only counted once
// unspliced is sum
////
// 0 header contig gene_name transcript_name exon_number ambiguous
// 5 strand pos next exon_type spliced
// 10 unspliced clipped exon_intron exon_other
// 15 skipped wrong_strand e_isoform
////




