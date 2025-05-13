use clap::{command, Parser, arg, value_parser, ArgAction, Command};
use std::collections::{HashMap, HashSet};
use std::ffi::OsStr;
use std::fmt::format;
use std::fs::{self, File};
use std::hash::Hash;
use std::io::prelude::*;
use std::io::{BufReader, BufWriter};
use std::path::{Path, PathBuf};
use std::process::{Command as Std_Command, Stdio};
use std::io::Write;
use std::fs::OpenOptions;
use strand_specifier_lib::Strand;
use std::str::FromStr;


pub enum SplicingChoice {
    Spliced(usize),
    Unspliced(usize),
    Clipped(usize),
    ExonOther(usize),
    Skipped(usize),
    WrongStrand(usize),
    Isoform(usize)
}

    #[derive(Debug, PartialEq, Eq)]
    pub struct SplicingChoiceError;

    impl FromStr for SplicingChoice {
    type Err = SplicingChoiceError;


    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s{
            "spliced" => Ok(SplicingChoice::Spliced(9)),
            "unspliced" => Ok(SplicingChoice::Unspliced(10)),
            "clipped" => Ok(SplicingChoice::Clipped(11)),
            "exon_other" => Ok(SplicingChoice::ExonOther(12)),
            "skipped" => Ok(SplicingChoice::Skipped(13)),
            "wrong_strand" => Ok(SplicingChoice::WrongStrand(14)),
            "isoform" => Ok(SplicingChoice::Isoform(15)),
            _ => Err(SplicingChoiceError)
        }
    }
}

impl SplicingChoice{
    pub fn index(&self) -> usize{
        match self{
            SplicingChoice::Spliced(n) => *n,
            SplicingChoice::Unspliced(n) => *n,
            SplicingChoice::Clipped(n) => *n,
            SplicingChoice::ExonOther(n) => *n,
            SplicingChoice::Skipped(n) => *n,
            SplicingChoice::WrongStrand(n) => *n,
            SplicingChoice::Isoform(n) => *n,
        }
    }
}

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
    pub const E_OTH: usize = 12;
    pub const SKIPPED: usize = 13;
    pub const WRONG_STRAND: usize = 14;
    pub const E_ISO: usize = 15;
}

struct Intron{
    contig: String,
    gene: String,
    transcript: String,
    intron_number: u16,
    donnor: String,
    acceptor: String,
    strand: Strand,
    spliced: u32,
    unspliced_donnor: u32,
    unspliced_acceptor: u32,
}

impl Intron{
    fn new(spt: &Vec<String>, counter: &Counter) -> Self{
        let mut donnor: u32;// = 0;
        let mut acceptor: u32;// = 0;
        let mut intron_n: u16;// = 0;


        let (spliced, unspliced)  = counter.count(spt);

        let mut donnor_uns = 0;
        let mut acceptor_uns = 0;
        

        intron_n = spt[3].split("_").collect::<Vec<&str>>()[1].parse::<u16>().unwrap();
        if spt[8] == "Acceptor"{
            acceptor_uns = unspliced;
            donnor = spt[7].parse::<u32>().unwrap();
            acceptor = spt[6].parse::<u32>().unwrap();
            intron_n -= 1;
        }
        else{
            donnor_uns = unspliced;
            donnor = spt[6].parse::<u32>().unwrap();
            acceptor = spt[7].parse::<u32>().unwrap();
        }

        Intron { contig: spt[0].to_string(),
            gene: spt[1].to_string(),
            transcript: spt[2].to_string(),
            intron_number: intron_n,
            donnor: donnor.to_string(),
            acceptor: acceptor.to_string(),
            strand: Strand::from(spt[5].as_str()),
            spliced: spliced,
            unspliced_donnor: donnor_uns,
            unspliced_acceptor: acceptor_uns }
    }

    fn update(&mut self, spt: &Vec<String>, counter: &Counter){
        let (_, unspliced)  = counter.count(spt);

        if spt[8] == "Acceptor"{
            self.unspliced_acceptor;
        }
        else{
            self.unspliced_donnor = unspliced;
        }
    }

    fn dump(&self) -> String{

        if self.spliced + self.unspliced_acceptor + self.unspliced_donnor == 0{
            format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", 
            self.contig, self.donnor, self.acceptor,
            self.strand, self.gene, "NA", self.spliced,
            self.unspliced_donnor, self.unspliced_acceptor, self.transcript,
            self.intron_number)
        }
        else{
        let ratio = (self.spliced) as f32 / (self.spliced + self.unspliced_acceptor + self.unspliced_donnor) as f32 ;
        format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", 
                self.contig, self.donnor, self.acceptor,
                self.strand, self.gene, ratio, self.spliced,
                self.unspliced_donnor, self.unspliced_acceptor, self.transcript,
                self.intron_number)
        }
    }
 }


 struct Counter{
    spliced: Vec<usize>,
    unspliced: Vec<usize>
 }

 impl Counter{
    fn new(spliced: Vec<usize>, unspliced: Vec<usize>) -> Self{
        Counter{
            spliced: spliced,
            unspliced: unspliced,
        }
    }

    fn count(&self, spt: &Vec<String>) -> (u32, u32){
        let spliced = self.spliced.iter()
        .map(|x| spt[*x].parse::<u32>().unwrap())       
        .fold(0, |acc, x| acc + x);
        let unspliced = self.unspliced.iter()
        .map(|x| spt[*x].parse::<u32>().unwrap())       
        .fold(0, |acc, x| acc + x);
        (spliced, unspliced)
    }
}


pub fn to_se_from_table(table_file: &str,
                        out_file: &str,
                        spliced_def: Vec<String>,
                        unspliced_def: Vec<String> ) -> (){

    println!("{}", table_file);
    let mut hashid: HashMap<usize, &str> = HashMap::new();
    hashid.insert(9, "spliced");
    hashid.insert(10, "unspliced");
    hashid.insert(11, "clipped");
    hashid.insert(12, "E_other");
    hashid.insert(13, "skipped");
    hashid.insert(14, "wrongStrand");   
    hashid.insert(15, "E_Isoform");    




    let g1: Vec<usize> = spliced_def.iter().map(|x| SplicingChoice::from_str(x).unwrap().index()).collect();
    let g2: Vec<usize> = unspliced_def.iter().map(|x| SplicingChoice::from_str(x).unwrap().index()).collect();
  

    
    let counter = Counter::new(g1.clone(), g2.clone());
    //let test_amb= true;

    let mut line: String;// = "".to_string();
    let file = table_file;
    let presorted = format!("{}.presorted", &out_file);
    
    let mut out_file_open = File::create_new(presorted.clone()) 
        .unwrap_or_else(|_| panic!("output file {} should not exist.", &presorted)); 

    let mut dict : HashMap<String, Intron> = HashMap::new();

    let mut donnor: u32;// = 0;
    let mut acceptor: u32; // = 0;

    let f = File::open(file).unwrap();
    let reader =  BufReader::with_capacity(64 * 1024, f); 


    let mut key: String; //= "".to_string();
    for (_i,l) in reader.lines().enumerate().skip(1){

        line = l.unwrap();
       // println!("{}: {}", i,  line);
        let spt = line.trim().split('\t').map(|s| s.to_owned()).collect::<Vec<String>>();
        if spt.len() <=1{
            continue
        }
        if spt[4] == "true"{
            continue
        }
        if spt[6] == "." {
            continue
        }
        if spt[7] == "."{
            continue
        }

        if spt[8] == "Acceptor"{
            
            donnor = spt[7].parse::<u32>().unwrap();
            acceptor = spt[6].parse::<u32>().unwrap();
        }
        else{
            donnor = spt[6].parse::<u32>().unwrap();
            acceptor = spt[7].parse::<u32>().unwrap();
        }

        key = format!("{}_{}_{}_{}_{}", spt[0], spt[1], spt[2], donnor, acceptor);
        dict.entry(key)
        .and_modify(|x| x.update(&spt, &counter))
        .or_insert(Intron::new(&spt, &counter));

    }

    for value in dict.values(){
        let _ = out_file_open.write_all(value.dump().as_bytes());

    }
    let _ = out_file_open.flush();


    {
        let mut out_final = File::create_new(out_file)
        .unwrap_or_else(|_| panic!("output file {} should not exist.", &out_file));
        let header=format!("# spliced definition: {};\n# unspliced definition: {};\nContig\tstart\tend\tstrand\tgeneID\tRatio\tspliced\tUnspliceDonnor\tUnsplicedAcceptor\ttranscriptID\tintronN\n",
          
          g1.iter().map(|x| hashid.get(x).unwrap().to_string()).collect::<Vec<String>>().join(" "),
          g2.iter().map(|x| hashid.get(x).unwrap().to_string()).collect::<Vec<String>>().join(" "));

        out_final.write_all(header.as_bytes()).unwrap();

    }

    {

    let output_cmd = Std_Command::new("sort") 
        .args(["-k", "1g", "-k", "5,5", "-k", "2n", &presorted])
        .output()
        .expect(" sort command failed to start");

    
    let mut file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(out_file).unwrap();

        // Write the output to file
        file.write_all(&output_cmd.stdout).unwrap();

    }

    

    let _ = Std_Command::new("rm")
        .args([&presorted.clone()])
        .spawn()
      .expect(" rm command failed to start").wait();



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
    /// space separated list the column to use for "unspliced" for the splicing defect table.
    /// you can regenrate this using the splicing_efficiency exe
    /// What to consider as unspliced? spliced, unspliced, clipped, exon_other, skipped,
    /// wrong_strand, isoform\n
    /// by default only use "-u unspliced" ->  unspliced (readthrough) reads \n
    /// to use unspliced and clipped : "-u unspliced clipped" 
    #[clap(long, value_parser, default_value = "unspliced", value_delimiter = ' ', num_args = 1..)]
    unspliced_def: Vec<String>,
    
    /// What to consider as spliced? spliced, unspliced, clipped, exon_other, skipped,
    /// wrong_strand, isoform\n
    /// by default only use "-u spliced" -> spliced (readthrough) reads \n
    /// to use spliced and isoform : "-u spliced isoform"
    #[clap(long, value_parser, default_value = "spliced", value_delimiter = ' ', num_args = 1..)]
    spliced_def: Vec<String>,
}


fn main(){
    let args = Args::parse();
    to_se_from_table(&args.input, &args.output, args.spliced_def, args.unspliced_def);
}

