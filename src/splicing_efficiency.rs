use clap::{command, Parser, arg, value_parser, ArgAction, Command};
use std::collections::{HashMap, HashSet};
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

pub fn to_se_from_table(table_file: &str, out_file: &str, unspliced_def: Vec<usize>, two_sided: bool) -> (){

    let g2 = [Field_::UNSPLICED];
    let test_amb= true;
    let two_sided = true;

    let mut line: String = "".to_string();
    let file = "/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/alzeihmer_madeleine/omnisplice/SRR22002266.table";
    let out_file = "/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/alzeihmer_madeleine/RR22002266.test.dev.se";
    let out_file_open = File::create_new(out_file.clone())
        .unwrap_or_else(|_| panic!("output file {} should not exist.", &out_file));//expect(&format!("output file {} should not exist.", &table));
    let mut out_stream = BufWriter::new(out_file_open);

    let f = File::open(file).unwrap();
    let mut reader = BufReader::new(f);
    

    let mut spt : Vec<String> = Vec::new();
    let mut previous_line: Option<Vec<String>> = None;

    let transcript_id = "".to_string();
    let gene_id = "".to_string();
    let mut spliced : u32;
    let mut unspliced: u32;


    let mut results: HashMap<(u32, u32), (String, u32, u32, String, Vec<String>, f32, u32, u32, String)> = HashMap::new();

    let mut current_genes : Option<String> = None;
    let mut start :u32 = 0;
    let mut end: u32 = 0;

    let _ = reader.read_line(&mut line);
    let header = line;
    // may be an iterator over gene and transcript would make the code easier to follows.
    for l in reader.lines().skip(1){
        line = l.unwrap();
        let spt = line.trim().split('\t').map(|s| s.to_owned()).collect::<Vec<String>>();
        
        // dumping the hashmap at each new genes;
        match current_genes{
            Some(ref gene_id)=> { if *gene_id != spt[Field_::GENE]{
                
                for value in results.values(){
                    let _ = out_stream.write(format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", value.0, value.1, value.2, value.3, value.5, value.6, value.7, value.4.join(";")).as_bytes());
                }

                results.clear();
                current_genes = Some(spt[Field_::GENE].to_string());
            }},
            _ => {current_genes = Some(spt[Field_::GENE].to_string());}
        }


        if previous_line.is_some(){
                if (previous_line.as_ref().unwrap()[Field_::NEXT] == ".") | (spt[Field_::NEXT] == "."){
                    previous_line = None;
                    continue;
                }
                if test_amb & ( (previous_line.as_ref().unwrap()[Field_::AMBIG] == "true") | (spt[Field_::AMBIG] == "true") ){
                    previous_line = None;
                    continue;
                }
                // spliced value is the same in spt or previous line. we do not sum them.
                spliced = spt[Field_::SPLICED].parse::<u32>().unwrap();
                // this one liner will be helpfull in the future, as it can take any a list of fields.
                unspliced = g2.iter()
                                    .map(|x| previous_line.as_ref().unwrap()[*x].parse::<u32>().unwrap())       
                                    .fold(0, |acc, x| acc + x);
                if two_sided == true{
                    unspliced += g2.iter()
                        .map(|x|  spt[*x].parse::<u32>().unwrap())       
                        .fold(0, |acc, x| acc + x);
                }
                start = spt[Field_::NEXT].parse::<u32>().unwrap();
                end = spt[Field_::POS].parse::<u32>().unwrap();
                results.entry((start, end))
                       .or_insert((spt[Field_::CONTIG].to_owned(), start, end, spt[Field_::GENE].to_owned(),
                         Vec::new(), 
                        spliced as f32/ (unspliced as f32 + spliced as f32), spliced, unspliced, spt[Field_::STRAND].to_owned()))
                        .4
                        .push(format!("{}_{}", spt[Field_::TRAN], previous_line.as_deref().unwrap()[Field_::EXON_N]));

                previous_line = None;
            }
            else{
                    previous_line = Some(spt.clone());
                }   

            }

    for value in results.values(){
                let _ = out_stream.write(format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", value.0, value.1, value.2, value.3, value.4.join(";"), value.5, value.6, value.7 ).as_bytes());
        }
    let _ = out_stream.flush();

}


#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {

    /// Name of the Input file
    #[arg(short, long, required = true)]
    input: String,

    /// Name of the Output file
    #[arg(short, long, required = true)]
    output: String,
    

    /// only consider donnor exons. and not the acceptor exon
    #[arg(long, default_value_t = false, action)]
    one_side : bool,

    /// What to consider as unspliced? unspliced: 9, clipped: 10, exon_intron: 11, exon_other: 12, skipped: 13, wrong_strand:14\n
    /// by default only use "-u 9" ->  unspliced (readthrough) reads \n
    /// to use unspliced and clipped : "-u 9 10"
    #[clap(long, value_parser, default_value = "9", value_delimiter = ' ', num_args = 1..)]
    unspliced_def: Vec<usize>,

}


fn main(){

    let args = Args::parse();
    let two_sided = !args.one_side;
    
    to_se_from_table(&args.input, &args.output, args.unspliced_def, two_sided);

}

