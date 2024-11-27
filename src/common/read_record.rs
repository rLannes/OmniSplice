use crate::common::utils::{ExonType, ReadAssign};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use strand_specifier_lib::Strand;

struct ReadtRecord {
    gene_id: String,
    transcript_id: String,
    contig: String,
    strand: Strand,
    exon_type: ExonType,
    pos: i64,
    assigned: Vec<(ReadAssign, u32)>,
}

impl ReadtRecord {
    fn new(spt: &Vec<&str>) -> Self {
        ReadtRecord {
            gene_id: spt[2].to_string(),
            transcript_id: spt[3].to_string(),
            contig: spt[0].to_string(),
            strand: Strand::from(spt[4]),
            exon_type: ExonType::from(spt[5]),
            pos: spt[1].parse::<i64>().expect("failed to parse number"),
            assigned: vec![(
                ReadAssign::from(spt[6]),
                spt[7].parse::<u32>().expect("failed to parse number"),
            )],
        }
    }
    fn update(&mut self, ass: ReadAssign, value: u32) {
        self.assigned.push((ass, value));
    }

    fn dump_junction(
        &self,
        next_acceptor: Option<i64>,
        junction_set: &HashSet<(i64, i64)>,
    ) -> String {
        let mut unspliced = 0;
        let mut spliced = 0;
        let mut clipped = 0;
        let mut exon_intron = 0;
        let mut exon_other = 0;
        let mut wrong_strand = 0;
        let mut skipped = 0;
        let mut e_isoform = 0;

        // adding Wrong Strand and Skipped

        for elem in &self.assigned {
            match elem.0 {
                ReadAssign::ReadThrough => unspliced = elem.1,
                ReadAssign::SoftClipped => clipped = elem.1,
                ReadAssign::WrongStrand => wrong_strand = elem.1,
                ReadAssign::Skipped(start, end) => skipped = elem.1,
                ReadAssign::ReadJunction(start, end) => {
                    // it is super verbose... but it works!
                    match (self.strand, self.exon_type, next_acceptor) {
                        (_, _, None) => (),
                        (Strand::Plus | Strand::NA, ExonType::Donnor, Some(next))
                        | (Strand::Minus, ExonType::Acceptor, Some(next)) => {
                            if next == end {
                                spliced = spliced + elem.1
                            } else if junction_set.contains(&(start, end)) {
                                e_isoform = e_isoform + elem.1
                            } else if end < next {
                                exon_intron = exon_intron + elem.1;
                            } else {
                                exon_other = exon_other + elem.1
                            }
                        }
                        (Strand::Minus, ExonType::Donnor, Some(next))
                        | (Strand::Plus | Strand::NA, ExonType::Acceptor, Some(next)) => {
                            if next == start {
                                spliced = spliced + elem.1;
                                //println!("junction: {} {}; pos {}; next: {}; spliced", start, end, self.pos, next);
                            } else if junction_set.contains(&(start, end)) {
                                e_isoform = e_isoform + elem.1
                            } else if start > next {
                                exon_intron = exon_intron + elem.1;
                            }
                            //println!("junction: {} {}; pos {};; next: {}; exon_intron", start, end, self.pos, next);
                            else {
                                exon_other = exon_other + elem.1;
                                //println!("junction: {} {}; pos {};; next: {}; exon_other", start, end, self.pos, next);
                            }
                        } //(_, _, _) => (),
                    }
                }
                _ => (),
            }
        }

        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.strand,
            self.pos,
            self.exon_type,
            spliced,
            unspliced,
            clipped,
            exon_intron,
            exon_other,
            skipped,
            wrong_strand,
            e_isoform
        )
    }
}

struct ReadtRecordContainer {
    container: Vec<ReadtRecord>,
    strand: Strand,
}

impl ReadtRecordContainer {
    fn new(strand: Strand) -> Self {
        ReadtRecordContainer {
            container: Vec::new(),
            strand: strand,
        }
    }

    fn sort(&mut self) -> () {
        match self.strand {
            Strand::Plus | Strand::NA => self.container.sort_by(|a, b| a.pos.cmp(&b.pos)),
            Strand::Minus => self.container.sort_by(|a, b| b.pos.cmp(&a.pos)),
        }
    }

    fn get_junction(&self, gene_junction_set: &mut HashSet<(i64, i64)>) -> () {
        //let mut result: HashSet<(i64, i64)> = HashSet::new();
        let ll = self.container.len();
        //let mut i = 1;
        for i in (1..(ll - 1)).step_by(2) {
            gene_junction_set.insert((self.container[i].pos, self.container[i + 1].pos));
        }
    }

    fn find(&mut self, exon_type: ExonType, pos: i64) -> Option<&mut ReadtRecord> {
        for record in &mut self.container {
            if record.pos == pos && record.exon_type == exon_type {
                return Some(record);
            }
        }
        None
    }

    fn insert(&mut self, spt: Vec<&str>) -> () {
        if let Some(record) = self.find(
            ExonType::from(spt[5]),
            spt[1].parse::<i64>().expect("failed to parse number"),
        ) {
            record.update(
                ReadAssign::from(spt[6]),
                spt[7].parse::<u32>().expect("failed to parse number"),
            );
        } else {
            self.container.push(ReadtRecord::new(&spt));
        }
    }

    fn dump(&mut self, junction_set: &HashSet<(i64, i64)>) -> String {
        self.sort();
        let mut result: Vec<String> = Vec::new();
        let mut exon_number: usize; // = 0;
        let size = self.container.len();
        let mut next: Option<i64>; // = None;
        for (indice, junction) in self.container.iter().enumerate() {
            next = None;
            exon_number = indice / 2;

            if (indice + 1 < size) && (indice > 0) {
                match (self.strand, &junction.exon_type) {
                    (Strand::Plus | Strand::NA, ExonType::Donnor)
                    | (Strand::Minus, ExonType::Donnor) => {
                        next = Some(self.container[indice + 1].pos);
                    }
                    (Strand::Plus | Strand::NA, ExonType::Acceptor)
                    | (Strand::Minus, ExonType::Acceptor) => {
                        next = Some(self.container[indice - 1].pos);
                    }
                }
            }
            result.push(format!(
                "{}\t{}\t{}\texon_{}\t{}",
                junction.contig,
                junction.gene_id,
                junction.transcript_id,
                exon_number + 1,
                junction.dump_junction(next, &junction_set)
            ));
        }

        format!("{}\n", result.join("\n"))
        //format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", self.strand, self.pos, self.exon_type,
        //unspliced, spliced, clipped, exon_intron, exon_other)
    }
}

pub fn file_to_table(file: String, out_file: &mut BufWriter<File>) -> () {
    let mut mymap = parse_file(file);

    for (_gene_name, container) in &mut mymap {
        let mut gene_junction_set: HashSet<(i64, i64)> = HashSet::new();
        for (_transcript_name, cont) in &mut *container {
            cont.sort();
            cont.get_junction(&mut gene_junction_set);
        }
        for (_transcript_name, cont) in container {
            let _ = out_file.write(cont.dump(&gene_junction_set).as_bytes());
        }
    }
}

fn parse_file(file: String) -> HashMap<String, HashMap<String, ReadtRecordContainer>> {
    // Refactor String -> String -> ReadtRecordContainer
    let mut resultamap: HashMap<String, HashMap<String, ReadtRecordContainer>> = HashMap::new();
    // read the file
    let f = File::open(file).expect("cannot open file");
    let reader = BufReader::new(f);
    let mut this_line: String; // = String::new();
    let mut spt: Vec<&str>;
    let mut gene_name: String;
    let mut tr_name: String;
    let mut gene_strand: Strand;

    for line in reader.lines() {
        this_line = line.expect("failed to read line in make_table");
        spt = this_line.trim().split('\t').collect::<Vec<&str>>();
        if &spt.len() < &5 {
            continue;
        }
        gene_name = spt[2].to_string();
        tr_name = spt[3].to_string();
        gene_strand = Strand::from(spt[4]);
        resultamap
            .entry(gene_name.clone())
            .or_insert_with(|| {
                let mut t = HashMap::new();
                t.insert(tr_name.clone(), ReadtRecordContainer::new(gene_strand));
                t
            })
            .entry(tr_name.clone())
            .or_insert_with(|| ReadtRecordContainer::new(gene_strand))
            .insert(spt);
    }
    resultamap
}
