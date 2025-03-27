use clap::{command, Parser, arg, value_parser, ArgAction, Command};
use std::collections::{HashMap, HashSet};
use std::ffi::OsStr;
use std::fmt::format;
use std::fs::File;
use std::hash::Hash;
use std::io::prelude::*;
use std::io::{BufReader, BufWriter};
use std::path::{Path, PathBuf};
use std::process::{Command as Std_Command, Stdio};
use std::io::Write;
use std::fs::OpenOptions;
use strand_specifier_lib::Strand;

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
    fn new(spt: &Vec<&str>, counter: Counter) -> Self{
        let mut donnor = 0;
        let mut acceptor = 0;
        let mut intron_n = 0;


        let (spliced, unspliced)  = counter.count(spt);

        let mut donnor_uns = 0;
        let mut acceptor_uns = 0;
        

        intron_n = spt[3].split("_").collect::<Vec<&str>>()[1].parse::<u16>().unwrap();
        if spt[8] == "Acceptor"{
            acceptor_uns = unspliced;
            donnor = spt[6].parse::<u32>().unwrap();
            acceptor = spt[7].parse::<u32>().unwrap();
            intron_n -= 1;
        }
        else{
            donnor_uns = unspliced;
            donnor = spt[7].parse::<u32>().unwrap();
            acceptor = spt[6].parse::<u32>().unwrap();
        }

        Intron { contig: spt[0].to_string(),
            gene: spt[1].to_string(),
            transcript: spt[2].to_string(),
            intron_number: intron_n,
            donnor: donnor.to_string(),
            acceptor: acceptor.to_string(),
            strand: Strand::from(spt[5]),
            spliced: spliced,
            unspliced_donnor: donnor_uns,
            unspliced_acceptor: acceptor_uns }
    }

    fn update(&mut self, spt: &Vec<&str>, counter: Counter){
        let (_, unspliced)  = counter.count(spt);

        if spt[8] == "Acceptor"{
            self.unspliced_acceptor;
        }
        else{
            self.unspliced_donnor = unspliced;
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

    fn count(self, spt: &Vec<&str>) -> (u32, u32){
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
                        spliced_def: Vec<usize>,
                        unspliced_def: Vec<usize> ) -> (){


    let mut hashid: HashMap<usize, &str> = HashMap::new();
    hashid.insert(9, "spliced");
    hashid.insert(10, "unspliced");
    hashid.insert(11, "clipped");
    hashid.insert(12, "E_intron");
    hashid.insert(13, "E_other");
    hashid.insert(14, "skipped");
    hashid.insert(15, "wrongStrand");   
    hashid.insert(16, "E_Isoform");    



    //let g1 = spliced_def;
    //let g2 = unspliced_def;

    let counter = Counter::new(spliced_def, unspliced_def);
    let test_amb= true;

    let mut line: String = "".to_string();
    let file = table_file;
    let presorted = format!("{}.presorted", &out_file);
    
    let mut out_file_open = File::create_new(presorted.clone()) 
        .unwrap_or_else(|_| panic!("output file {} should not exist.", &presorted)); 
    //let mut out_stream = BufWriter::new(out_file_open);

    
    let f = File::open(file).unwrap();
    let mut reader = BufReader::new(f);
    

    //let mut spt : Vec<String> = Vec::new();
    let mut previous_line: Option<Vec<String>> = None;

    let mut dict : HashMap<String, Intron> = HashMap::new();
    //let transcript_id = "".to_string();
    //let gene_id = "".to_string();
    let mut spliced : u32;
    let mut unspliced: u32;
    let mut unspliced_acceptor: u32;


    let mut results: HashMap<(u32, u32), (String, u32, u32, String, Vec<String>, f32, u32, u32, u32, String)> = HashMap::new();

    let mut current_genes : Option<String> = None;
    let mut start :u32 = 0;
    let mut end: u32 = 0;

    //let _ = reader.read_line(&mut line);
    //let header = line;
    // may be an iterator over gene and transcript would make the code easier to follows.

    for l in reader.lines().skip(1){
        line = l.unwrap();
        let spt = line.trim().split('\t').map(|s| s.to_owned()).collect::<Vec<String>>();
        
        // dumping the hashmap at each new genes;
        match current_genes{
            Some(ref gene_id)=> { if *gene_id != spt[Field_::GENE]{
                for value in results.values(){
                    let _ = out_file_open.write_all(format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", value.0, value.1, value.2, value.9, value.3, value.5, value.6, value.7, value.8, value.4.join(";")).as_bytes());
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
                spliced = g1.iter()
                        .map(|x| spt[*x].parse::<u32>().unwrap())       
                        .fold(0, |acc, x| acc + x);
                // this one liner will be helpfull in the future, as it can take any a list of fields.
                unspliced = g2.iter()
                                    .map(|x| previous_line.as_ref().unwrap()[*x].parse::<u32>().unwrap())       
                                    .fold(0, |acc, x| acc + x);
                unspliced_acceptor = g2.iter()
                        .map(|x|  spt[*x].parse::<u32>().unwrap())       
                        .fold(0, |acc, x| acc + x);

                start = spt[Field_::NEXT].parse::<u32>().unwrap();
                end = spt[Field_::POS].parse::<u32>().unwrap();
                results.entry((start, end))
                       .or_insert((spt[Field_::CONTIG].to_owned(), start, end, spt[Field_::GENE].to_owned(),
                         Vec::new(), 
                        spliced as f32/ (unspliced as f32 + spliced as f32 + unspliced_acceptor as f32), spliced, unspliced, unspliced_acceptor, spt[Field_::STRAND].to_owned()))
                        .4
                        .push(format!("{}_{}", spt[Field_::TRAN], previous_line.as_deref().unwrap()[Field_::EXON_N]));

                previous_line = None;
            }
            else{
                    previous_line = Some(spt.clone());
                }   
            }

    for value in results.values(){
        let _ = out_file_open.write_all(format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", value.0, value.1, value.2, value.9, value.3, value.5, value.6, value.7, value.8, value.4.join(";")).as_bytes());
        //let _ = out_stream.flush().unwrap();
        }
    //out_stream.flush().unwrap();
    //out_stream.into_inner().unwrap().sync_all().unwrap();

    //drop(out_stream);

    {
        let mut out_final = File::create_new(out_file)
        .unwrap_or_else(|_| panic!("output file {} should not exist.", &out_file));
        let header=format!("# spliced definition: {};\n# unspliced definition: {};\nContig\tstart\tend\tstrand\tgeneID\tRatio\tspliced\tUnspliceDonnor\tUnsplicedAcceptor\ttranscriptID_exonN\n",
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
    /// What to consider as unspliced? unspliced: 10, clipped: 11, exon_intron: 12,
    /// exon_other: 13, skipped: 14, wrong_strand:15, isoform:16\n
    /// by default only use "-u 10" ->  unspliced (readthrough) reads \n
    /// to use unspliced and clipped : "-u 10 11"
    #[clap(long, value_parser, default_value = "10", value_delimiter = ' ', num_args = 1..)]
    unspliced_def: Vec<usize>,
    
    /// What to consider as unspliced? unspliced: 10, clipped: 11, exon_intron: 12,
    /// exon_other: 13, skipped: 14, wrong_strand:15, isoform: 16\n
    /// by default only use "-u 9" -> spliced (readthrough) reads \n
    /// to use spliced and isform : "-u 9 16"
    #[clap(long, value_parser, default_value = "9", value_delimiter = ' ', num_args = 1..)]
    spliced_def: Vec<usize>,
}


fn main(){
    let args = Args::parse();
    to_se_from_table(&args.input, &args.output, args.spliced_def, args.unspliced_def);
}

