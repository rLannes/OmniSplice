use crate::common::error::OmniError;
//use crate::common::point::get_attr_id;/
use crate::common::gtf_::{get_attr_id, get_junction_from_gtf};
use crate::common::utils::{ExonType, ReadAssign};
use bio::data_structures::interval_tree::IntervalTree;
use clap::builder::Str;
use log::{debug, error, info, trace, warn};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::hash::Hash;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::mem;
use std::sync::Arc;
use strand_specifier_lib::{LibType, Strand};
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
        junction_set: &HashMap<(String, i64, i64), Strand>,
        valid_j_gene: &HashSet<(i64, i64)>,
    ) -> String {
        let mut current_strand = Strand::Plus;
        let mut unspliced = 0;
        let mut spliced = 0;
        let mut clipped = 0;
        let mut exon_intron = 0;
        let mut exon_other = 0;
        let mut wrong_strand = 0;
        let mut skipped = 0;
        let mut e_isoform = 0;
        let mut next_id: String;

        // adding Wrong Strand and Skipped
        let mut start_order: i64 = 0;
        let mut end_order: i64 = 0;

        for elem in &self.assigned {
            match elem.0 {
                ReadAssign::ReadThrough => unspliced = elem.1,
                ReadAssign::SoftClipped => clipped = elem.1,
                ReadAssign::WrongStrand => wrong_strand = elem.1,
                ReadAssign::Skipped(start, end) => {
                    start_order = start;
                    end_order = end;

                    //if valid_j_gene.contains(&(start_order, end_order))  | valid_j_gene.contains(&(end_order, start_order))
                    if valid_j_gene.contains(&(start_order, end_order))
                        | valid_j_gene.contains(&(end_order, start_order))
                    {
                        e_isoform = e_isoform + elem.1;
                        continue;
                    } else {
                        skipped = elem.1
                    }
                }
                ReadAssign::ReadJunction(start, end) => {
                    start_order = start;
                    end_order = end;
                    current_strand = self.strand.clone();

                    match (self.strand, self.exon_type, next_acceptor) {
                        (_, _, None) => (),
                        (Strand::Plus | Strand::NA, ExonType::Donnor, Some(next))
                        | (Strand::Minus, ExonType::Acceptor, Some(next)) => {
                            if next == end {
                                spliced = spliced + elem.1;
                            } else if valid_j_gene.contains(&(start_order, end_order))
                                | valid_j_gene.contains(&(end_order, start_order))
                            {
                                e_isoform = e_isoform + elem.1;
                            } else {
                                exon_other = exon_other + elem.1;
                            }
                        }
                        (Strand::Minus, ExonType::Donnor, Some(next))
                        | (Strand::Plus | Strand::NA, ExonType::Acceptor, Some(next)) => {
                            if next == start {
                                spliced = spliced + elem.1;
                                continue;
                            } else if valid_j_gene.contains(&(start_order, end_order))
                                | valid_j_gene.contains(&(end_order, start_order))
                            {
                                e_isoform = e_isoform + elem.1;
                            } else {
                                exon_other = exon_other + elem.1;
                            }
                        }
                    }
                }
                _ => (),
            }
        }
        next_id = match next_acceptor {
            Some(i) => i.to_string(),
            _ => ".".to_string(),
        };
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.strand,
            self.pos,
            next_id,
            self.exon_type,
            spliced,
            unspliced,
            clipped,
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

    fn dump(
        &mut self,
        junction_set: &HashMap<(String, i64, i64), Strand>,
        invalid: &HashSet<(String, i64)>,
        valid_j_gene: &HashSet<(i64, i64)>,
    ) -> String {
        self.sort();
        let mut ambigous = false;
        let mut result: Vec<String> = Vec::new();
        let mut exon_number: usize; // = 0;
        let size = self.container.len();
        let mut next: Option<i64>; // = None;

        for (indice, junction) in self.container.iter().enumerate() {
            ambigous = false;
            next = None;
            ambigous = false;
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

            ambigous = if invalid.contains(&(junction.contig.clone(), junction.pos)) {
                true
            } else {
                false
            };

            result.push(format!(
                "{}\t{}\t{}\texon_{}\t{}\t{}",
                junction.contig,
                junction.gene_id,
                junction.transcript_id,
                exon_number + 1,
                ambigous,
                junction.dump_junction(next, junction_set, valid_j_gene),
            ));
        }

        format!("{}\n", result.join("\n"))
    }
}

pub fn file_to_table(
    file: String,
    out_file: &mut BufWriter<File>,
    gtf: &str,
    libtype: LibType,
    valid_j_gene: &HashMap<String, HashSet<(i64, i64)>>,
) -> Result<(), OmniError> {
    let mut mymap = parse_file(file.as_str());

    // ambigious position
    let invalid_pos = get_invalid_pos(gtf)?;
    //isoform
    let gene_junction_set = get_junction_from_gtf(gtf, &libtype)?;

    for (_gene_name, container) in &mut mymap {
        let all_gene_junction = valid_j_gene.get(_gene_name).ok_or(OmniError::GTFParse(
            "failed to recover junction".to_string(),
        ))?;

        for (_transcript_name, cont) in &mut *container {
            cont.sort();
        }

        for (_transcript_name, cont) in container {
            //let sub_set: Option<&HashSet<(i64, i64, Strand)>> = gene_junction_set.get(&cont.container[0].contig);

            let _ = out_file.write(
                cont.dump(&gene_junction_set, &invalid_pos, all_gene_junction)
                    .as_bytes(),
            );
        }
    }
    out_file.write("\n".as_bytes());
    out_file.flush();
    Ok(())
}

fn parse_file(file: &str) -> HashMap<String, HashMap<String, ReadtRecordContainer>> {
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

    for line in reader.lines().skip(1) {
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

#[derive(Eq, Hash, PartialEq, Clone, Copy, Debug)]
struct Intervall<T: Ord + Copy + Eq + Hash + PartialEq> {
    start: T,
    end: T,
}

impl<T> Intervall<T>
where
    T: Ord + Copy + Eq + Hash + PartialEq,
{
    fn new(start: T, end: T) -> Self {
        Intervall { start, end }
    }

    fn overlap(&self, other: &Intervall<T>) -> bool {
        if (self.start > other.end) | (self.end < other.start) {
            false
        } else {
            true
        }
    }
}

fn gtf_to_it(file: &str) -> HashMap<String, IntervalTree<i64, String>> {
    //let file = "genomic.gtf";
    let f = File::open(file).unwrap();
    let reader = BufReader::new(f);
    let mut this_line: String; //::new();

    let mut start: i64; //; = 0;
    let mut chr_: String; // = "".to_string();
    let mut end: i64; // = 0;
    let mut gene_name: String = "".to_string();
    let mut transcript_id: String; // = "".to_string();

    let mut result: HashMap<String, IntervalTree<i64, String>> = HashMap::new();

    for line in reader.lines() {
        this_line = line.unwrap();
        let spt = this_line.trim().split('\t').collect::<Vec<&str>>();
        if spt.len() < 8 {
            continue;
        }
        if spt[2] != "exon" {
            continue;
        }

        chr_ = spt[0].to_string();
        start = spt[3].parse::<i64>().unwrap() - 1;
        end = spt[4].parse::<i64>().unwrap();

        if let Some(gene_tmp) = get_attr_id(spt[8], "gene_id") {
            gene_name = gene_tmp;
        } else if let Some(gene_tmp) = get_attr_id(spt[8], "gene_name") {
            gene_name = gene_tmp;
        } else {
            print!("WARNING Cannot find {:?}", this_line);
            continue;
        }
        result
            .entry(chr_)
            .or_insert_with(|| IntervalTree::new())
            .insert(start..end, gene_name);
    }
    result
}

fn graph_from_gtf(
    file: &str,
) -> Result<HashMap<String, HashMap<Intervall<i64>, HashSet<Intervall<i64>>>>, OmniError> {
    //  TO remove clutter ca use a hashset to remove identical  exon
    // from same genes

    //let file = "genomic.gtf";
    let mut it_dico: HashMap<String, IntervalTree<i64, String>> = gtf_to_it(file);

    let f = File::open(file)?;
    let reader = BufReader::new(f);
    let mut this_line: String; //::new();

    let mut start: i64; //; = 0;
    let mut chr_: String; // = "".to_string();
    let mut end: i64; // = 0;
    let mut gene_name: String = "".to_string();
    let mut transcript_id: String; // = "".to_string();

    let mut seen: HashSet<(String, i64, i64, String)> = HashSet::new();

    let mut g: HashMap<String, HashMap<Intervall<i64>, HashSet<Intervall<i64>>>> = HashMap::new();

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

        start = spt[3].parse::<i64>()? - 1;
        end = spt[4].parse::<i64>()?; //- 1;

        if let Some(gene_tmp) = get_attr_id(spt[8], "gene_id") {
            gene_name = gene_tmp;
        } else if let Some(gene_tmp) = get_attr_id(spt[8], "gene_name") {
            gene_name = gene_tmp;
        } else {
            print!("WARNING Cannot find {:?}", this_line);
            continue;
        }

        if let Some(subtree) = it_dico.get(&chr_) {
            for inter in subtree.find(start..end) {
                // exactly same intervall.
                if (inter.interval().start == start) & (inter.interval().end == end) {
                    continue;
                } else {
                    g.entry(chr_.clone())
                        .or_insert_with(|| HashMap::new())
                        .entry(Intervall::new(start, end))
                        .or_insert_with(|| HashSet::new())
                        .insert(Intervall::new(inter.interval().start, inter.interval().end));
                }
            }
        }
    }
    Ok(g)
}

fn get_invalid_pos(file: &str) -> Result<HashSet<(String, i64)>, OmniError> {
    let g = graph_from_gtf(file)?;
    let mut seen = HashSet::new();
    let mut inter_vec = Vec::new();
    let mut e1 = 0;
    let mut borne: &Intervall<i64>;
    let mut results: HashSet<(String, i64)> = HashSet::new();
    let mut file_: Vec<Intervall<i64>> = Vec::new();

    let mut current_neihbors: Vec<&Intervall<i64>>;
    for (chr_, subdict) in g.iter() {
        for current_node in subdict.keys() {
            inter_vec.clear();
            if seen.contains(current_node) {
                continue;
            }

            file_.push(current_node.clone());
            while let Some(bfs_node) = file_.pop() {
                if seen.contains(&bfs_node) {
                    continue;
                }
                inter_vec.push(bfs_node);
                seen.insert(bfs_node);
                if let Some(current_neihbors_hash) = subdict.get(&bfs_node) {
                    current_neihbors = current_neihbors_hash
                        .iter()
                        .collect::<Vec<&Intervall<i64>>>();
                    for n in current_neihbors {
                        if seen.contains(n) {
                            continue;
                        }
                        file_.push(n.clone())
                    }
                }
            }

            e1 = inter_vec[0].start;
            if !(inter_vec.iter().all(|x| x.start == e1)) {
                //println!("{} {:?}", chr_, inter_vec );
                borne = inter_vec
                    .iter()
                    .min_by_key(|x| x.start)
                    .ok_or(OmniError::Expect(
                        "value expected in get_invalid_pos".to_string(),
                    ))?;
                for i in &inter_vec {
                    if i == borne {
                        continue;
                    } else {
                        results.insert((chr_.clone(), i.start));
                    }
                }
            }

            e1 = inter_vec[0].end;
            if !(inter_vec.iter().all(|x| x.end == e1)) {
                //println!("{} {:?}", chr_, inter_vec );
                borne = inter_vec
                    .iter()
                    .max_by_key(|x| x.end)
                    .ok_or(OmniError::Expect(
                        "value expected in get_invalid_pos".to_string(),
                    ))?;
                for i in &inter_vec {
                    if i == borne {
                        continue;
                    } else {
                        results.insert((chr_.clone(), i.end));
                    }
                }
            }
            // limit
            //inter_vec.clear();
            //}
        }
    }

    Ok(results)
}
