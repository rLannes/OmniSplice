use crate::common::utils;
use crate::get_header;
use bio::stats::bayesian::bayes_factors::evidence::KassRaftery;
use rust_htslib::bam::header;
use rust_htslib::bam::record::Record;
use std::collections::HashMap;
use std::fmt::format;
use std::fs;
use std::fs::File;
use std::hash::Hash;
use std::io::prelude::*;
use std::io::{BufReader, BufWriter};
use strand_specifier_lib::Strand;

struct Counts {
    spliced: u32,
    unspliced: u32,
    clipped: u32,
    exon_other: u32,
    skipped: u32,
    wrong_strand: u32,
    e_isoform: u32,
}

impl Counts {
    fn new(string: &[&str]) -> Self {
        let c = string
            .iter()
            .map(|x| x.parse::<u32>().unwrap())
            .collect::<Vec<u32>>();
        Counts {
            spliced: c[0],
            unspliced: c[1],
            clipped: c[2],
            exon_other: c[3],
            skipped: c[4],
            wrong_strand: c[5],
            e_isoform: c[6],
        }
    }

    fn update(&mut self, string: &[&str]) -> () {
        let c = string
            .iter()
            .map(|x| x.parse::<u32>().unwrap())
            .collect::<Vec<u32>>();

        // no c[0] because spliced represents the same event on both side of the junctions
        self.unspliced += c[1];
        self.clipped += c[2];
        self.exon_other += c[3];
        self.skipped += c[4];
        self.wrong_strand += c[5];
        self.e_isoform += c[6];
    }

    fn dump(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.spliced,
            self.unspliced,
            self.clipped,
            self.exon_other,
            self.skipped,
            self.wrong_strand,
            self.e_isoform
        )
    }
}

struct Junction {
    contig: String,
    gene: String,
    transcript: String,
    ambiguous: bool,
    intron_number: u16,
    donnor: String,
    acceptor: String,
    strand: Strand,
    counts: Counts,
}

impl Junction {
    fn dump(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.contig,
            self.gene,
            self.transcript,
            self.intron_number,
            self.donnor,
            self.acceptor,
            self.strand,
            self.ambiguous,
            self.counts.dump()
        )
    }

    fn new(string: &Vec<&str>, header: &HashMap<String, usize>) -> Self {
        let exon_type = string[utils::get_header!("exon_type", header)];
        let exon_n: u16 = string[utils::get_header!("exon_number", header)]
            .split("_")
            .collect::<Vec<&str>>()[1]
            .parse::<u16>()
            .expect("unable to parse exon number junction");

        let pos = if exon_type == "Donnor" {
            string[utils::get_header!("pos", header)]
        } else {
            string[utils::get_header!("next", header)]
        };
        let next = if exon_type == "Donnor" {
            string[utils::get_header!("next", header)]
        } else {
            string[utils::get_header!("pos", header)]
        };

        let exon_n: u16 = string[utils::get_header!("exon_number", header)]
            .split("_")
            .collect::<Vec<&str>>()[1]
            .parse::<u16>()
            .expect("unable to parse exon number junction");

        Junction {
            contig: string[utils::get_header!("contig", header)].to_string(),
            gene: string[utils::get_header!("gene_name", header)].to_string(),
            transcript: string[utils::get_header!("transcript_name", header)].to_string(),
            donnor: pos.to_string(),
            acceptor: next.to_string(),
            ambiguous: if string[utils::get_header!("ambiguous", header)] == "true" {
                true
            } else {
                false
            },
            intron_number: if exon_type == "Donnor" {
                exon_n
            } else {
                exon_n - 1
            },
            strand: Strand::from(string[utils::get_header!("strand", header)]),
            counts: Counts::new(&string[utils::get_header!("spliced", &header)..]),
        }
    }

    fn update(&mut self, spt: &Vec<&str>, header: &HashMap<String, usize>) -> () {
        if spt[utils::get_header!("ambiguous", &header)] == "true" {
            self.ambiguous = true;
        }
        self.counts.update(&spt[get_header!("spliced", &header)..]);
    }
}

pub fn junction_file_from_table(table_file: &str, junction_file: &str) {
    let mut out_file_open =
        File::create_new(junction_file.clone()) //presorted out_file.clone()
            .unwrap_or_else(|_| panic!("output file {} should not exist.", &junction_file)); //expect(&format!("output file {} should not exist.", &table));
                                                                                             //let mut out_stream = BufWriter::new(out_file_open);

    let f = File::open(table_file).unwrap();
    let mut reader = BufReader::new(f);

    let mut header_line: String = "".to_string();
    let _ = reader
        .read_line(&mut header_line)
        .expect("fail to read table file for junction");
    let header = utils::parse_header(&header_line);

    let mut line: String = "".to_string();
    let mut exon_type = "";
    let mut intron_n: u16 = 0;
    let mut exon_n: u16 = 0;
    let mut key = "".to_string();

    let mut map: HashMap<String, Junction> = HashMap::new();

    for l in reader.lines() {
        line = l.expect("cannot parse line table to junction");

        let spt = line.trim().split('\t').collect::<Vec<&str>>();
        if spt.len() <= 1 {
            break;
        }
        if (spt[get_header!("pos", header)] == ".") || (spt[get_header!("next", header)] == ".") {
            continue;
        }

        let exon_type = spt[utils::get_header!("exon_type", header)];
        exon_n = spt[utils::get_header!("exon_number", header)]
            .split("_")
            .collect::<Vec<&str>>()[1]
            .parse::<u16>()
            .expect("unable to parse exon number junction");
        intron_n = if exon_type == "Donnor" {
            exon_n
        } else {
            exon_n - 1
        };

        key = format!(
            "{}_{}_{}",
            spt[utils::get_header!("gene_name", header)],
            spt[utils::get_header!("transcript_name", header)],
            intron_n
        );
        //println!("{}", key);
        map.entry(key)
            .and_modify(|e| e.update(&spt, &header))
            .or_insert(Junction::new(&spt, &header));
    }

    let mut junctions = map.values().collect::<Vec<&Junction>>();

    junctions.sort_by_key(|x| (&x.gene, &x.transcript, x.intron_number));

    out_file_open.write_all("Contig\tGene\tTranscript\tIntron\tDonnor\tAcceptor\tStrand\tAmbiguous\tspliced\tunspliced\tclipped\texon_other\tskipped\twrong_strand\te_isoform\n".as_bytes());
    for j in junctions {
        out_file_open.write_all(format!("{}\n", j.dump()).as_bytes());
    }

    //out_stream.flush().unwrap();
    //out_stream.into_inner().unwrap().sync_all().unwrap();
}

mod tests_it {
    use super::*;

    #[test]
    fn test_1() {
        let table_file = "/lab/solexa_yamashita/people/Romain/Projets/Adrienne/PAPER/checking/omni_Y_mau_onlyY/trimmed_sim_mau_1_S13_L002_VS_Mau.table";
        let out_file = "test.test";
        junction_file_from_table(table_file, out_file);
    }
}
