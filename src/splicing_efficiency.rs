use clap::Parser;
use std::collections::HashMap;
use std::ffi::OsStr;
use std::fmt::format;
use std::fs::File;
use std::hash::Hash;
use std::io::prelude::*;
use std::io::{BufReader, BufWriter};
use std::path::{Path, PathBuf};

/// This is a self contain file that parse the table file and output an splicing effiency table.
/// i am unsure if I keep it separate or integrate it to the main software or both
/// it is to replicate and more the exact formual of spilcing e
//
// we want to support:
// - either Donnor only Acceptor only or both
// - we want to support the combination of different fields. "unspliced" + "exon_other" ...
//

//hard coded field but very fast!
#[non_exhaustive]
pub struct Field_;
impl Field_{
    pub const CONTIG: usize = 0;
    pub const GENE: usize = 1;
    pub const TRAN: usize = 2;
    pub const EXON_N: usize = 3;
    pub const AMBIG: usize = 4;
    pub const STRAND: usize = 5;
    pub const POS: usize = 6;
    pub const NEXT: usize = 7;
    pub const EXON_T: usize = 8;
    pub const SPLICED: usize = 9;
    pub const UNSPLICED: usize = 10;
    pub const CLIPPED: usize = 11;
    pub const E_INT: usize = 12;
    pub const E_OTH: usize = 13;
    pub const SKIPPED: usize = 14;
    pub const WRONG_STRAND: usize = 15;
    pub const E_ISO: usize = 16;
}


fn main(){

    let g2 = [Field_::UNSPLICED];
    let test_amb= true;
    let mut line: String = "".to_string();
    let file = "/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/alzeihmer_madeleine/omnisplice/SRR22002266.table";
    let f = File::open(file).unwrap();
    let mut reader = BufReader::new(f);
    
    let _ = reader.read_line(&mut line);
    let header = line;
    let mut spt : Vec<String> = Vec::new();
    let mut previous_line: Option<Vec<String>> = None;

    let transcript_id = "".to_string();
    let gene_id = "".to_string();
    let mut spliced : u32;
    let mut unspliced: u32;



    for l in reader.lines(){
        line = l.unwrap();
        let spt = line.trim().split('\t').map(|s| s.to_owned()).collect::<Vec<String>>();

        if previous_line.is_some(){
                if (previous_line.as_ref().unwrap()[Field_::NEXT] == ".") | (spt[Field_::NEXT] == "."){
                    previous_line = None;
                    continue;
                }
                if test_amb & ( (previous_line.as_ref().unwrap()[Field_::AMBIG] == "true") | (spt[Field_::AMBIG] == "true") ){
                    previous_line = None;
                    continue;
                }
                spliced = spt[Field_::SPLICED].parse::<u32>().unwrap();
                // will be helpfull in the future
                unspliced = g2.iter()
                                    .map(|x| spt[*x].parse::<u32>().unwrap())
                                    .fold(0, |acc, x| acc + x);
                // need to update the hashmap
            }
            else{
                    previous_line = Some(spt.clone());
                }   

            }

}
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




