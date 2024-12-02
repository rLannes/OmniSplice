use crate::common::point::get_attr_id;
use crate::common::utils::{ExonType, ReadAssign};
use bio::data_structures::interval_tree::IntervalTree;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::hash::Hash;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::sync::Arc;
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
        junction_set: Option<&HashSet<(i64, i64)>>,
    ) -> String {
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
                            //if self.gene_id == "gene-LOC27208819"{
                            //    println!("js {:?}",  junction_set);
                            //    println!("this: {}, next {:?}",  self.pos, next_acceptor);
                            //};
                            //println!("Acc junction: {} {}; pos {}; next: {}; n: {} all", start, end, self.pos, next, elem.1);
                            if next == end {
                                spliced = spliced + elem.1;
                                //println!("spliced");
                            } else if junction_set.is_some() && junction_set.unwrap().contains(&(start, end)) {
                                e_isoform = e_isoform + elem.1;
                                //println!("iso");
                            } else if end < next {
                                //println!("intron");
                                exon_intron = exon_intron + elem.1;
                                //println!("intron");
                            } else {
                                exon_other = exon_other + elem.1;
                                //println!("other");
                            }
                        }
                        (Strand::Minus, ExonType::Donnor, Some(next))
                        | (Strand::Plus | Strand::NA, ExonType::Acceptor, Some(next)) => {
                            //if self.gene_id == "gene-LOC27208819"{
                            //    println!("js {:?}",  junction_set);
                            //    println!("this: {}, next {:?}",  self.pos, next_acceptor);
                            ///};
                            //println!(" Doo junction: {} {}; pos {}; next: {}; n: {} all", start, end, self.pos, next, elem.1);
                            if next == start {
                                //println!("spliced");
                                spliced = spliced + elem.1;
                                //println!("spliced");
                                //println!("junction: {} {}; pos {}; next: {}; spliced", start, end, self.pos, next);
                            } else if junction_set.is_some() && junction_set.unwrap().contains(&(start, end)) {
                                e_isoform = e_isoform + elem.1;
                                //println!("iso");
                            } else if start > next {
                                //println!("intron");
                                exon_intron = exon_intron + elem.1;
                                //println!("intron");
                            }
                            //println!("junction: {} {}; pos {};; next: {}; exon_intron", start, end, self.pos, next);
                            else {
                                //println!("other");
                                exon_other = exon_other + elem.1;
                                //println!("other");
                                //println!("junction: {} {}; pos {};; next: {}; exon_other", start, end, self.pos, next);
                            }
                        } //(_, _, _) => (),
                    }
                }
                _ => (),
            }
        }
        next_id = match next_acceptor{
            Some(i) => i.to_string(),
            _ => ".".to_string()
        };
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.strand,
            self.pos,
            next_id,
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

    fn dump(&mut self, junction_set: Option<&HashSet<(i64, i64)>>, invalid: &HashSet<(String, i64)>) -> String {
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
                junction.dump_junction(next, junction_set)
            ));
        }

        format!("{}\n", result.join("\n"))
        //format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", self.strand, self.pos, self.exon_type,
        //unspliced, spliced, clipped, exon_intron, exon_other)
    }
}

pub fn file_to_table(file: String, out_file: &mut BufWriter<File>, gtf: &str) -> () {
    let mut mymap = parse_file(file.as_str());

    let invalid_pos = get_invalid_pos(gtf);

    let gene_junction_set = get_junction_from_gtf(gtf);

    let gene_junction_set = get_junction_from_gtf(gtf.clone());

    for (_gene_name, container) in &mut mymap {
        //if !gene_junction_set.contains_key(_gene_name){
        //    println!("{}", _gene_name);
        //}
        let sub_set = gene_junction_set.get(_gene_name);
        for (_transcript_name, cont) in &mut *container {
            cont.sort();
        }

        for (_transcript_name, cont) in container {
            let _ = out_file.write(cont.dump(sub_set, &invalid_pos).as_bytes());
        }
    }
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



#[derive(Eq, Hash, PartialEq, Clone, Copy)]
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

fn get_junction_from_gtf(file: &str) -> HashMap<String, HashSet<(i64, i64)>> {
    //let file = "genomic.gtf";
    let f = File::open(file).unwrap();
    let reader = BufReader::new(f);
    let mut this_line: String; //::new();

    let mut start: i64; //; = 0;
    let mut chr_: String; // = "".to_string();
    let mut end: i64; // = 0;
    let mut gene_name: String = "".to_string();
    let mut transcript_id: String; // = "".to_string();

    let mut result: HashMap<String, HashSet<(i64, i64)>> = HashMap::new();
    let mut transcript_vec: Vec<(i64, i64)> = Vec::new();
    let mut current_gene_id = "".to_string();
    let mut current_transcript_id = "".to_string();
    let mut myset: HashSet<(i64, i64)> = HashSet::new();

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

        if let Some((gene_tmp, tr_tmp)) =
            get_attr_id(spt[8], "gene_id").zip(get_attr_id(spt[8], "transcript_id"))
        {
            gene_name = gene_tmp;
            transcript_id = tr_tmp;
        } else if let Some((gene_tmp, tr_tmp)) =
            get_attr_id(spt[8], "gene_name").zip(get_attr_id(spt[8], "transcript_id"))
        {
            gene_name = gene_tmp;
            transcript_id = tr_tmp;
        } else {
            println!("warning gene id or transcript id not found: {:?}", spt);
            continue; // skip this gene
        }



        if gene_name != current_gene_id {

            let ll = transcript_vec.len();
            if ll > 1 {
                for i in 0..(ll - 1) {
                    myset.insert((transcript_vec[i].1, transcript_vec[i + 1].0));
                }
            }
            //if gene_name == "gene-LOC27208819"{
            //    println!("G {:?}", transcript_vec);
            //    println!("G {:?}", myset);
            //};


            transcript_vec.clear();

            if !myset.is_empty(){
                result.insert(current_gene_id.clone(), myset.clone());
                //if gene_name == "gene-LOC27208819"{
                //    println!("G b {:?}", result.get("gene-LOC27208819"));
                //};
                myset.clear();
                //if gene_name == "gene-LOC27208819"{
                //    println!("G a {:?}", result.get("gene-LOC27208819"));
                //};
            }
            current_gene_id = gene_name.clone();
            current_transcript_id = transcript_id;



        } else if transcript_id != current_transcript_id {


            let ll = transcript_vec.len();
            if ll > 1 {
                for i in 0..(ll - 1) {
                    myset.insert((transcript_vec[i].1, transcript_vec[i + 1].0));
                }
            }

            transcript_vec.clear();
            current_transcript_id = transcript_id;
        }

        transcript_vec.push((start, end));
    }


    let ll = transcript_vec.len();
    if ll > 1 {
        for i in 0..(ll - 1) {
            myset.insert((transcript_vec[i].1, transcript_vec[i + 1].0));
        }
    }
    //transcript_vec.clear();

    result.insert(gene_name.clone(), myset.clone());
    //myset.clear();
    result
}


fn gtf_toIT(file: &str) -> HashMap<String, IntervalTree<i64, String>> {
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
            print!("WARNING Cannort find {:?}", this_line);
            continue;
        }
        result
            .entry(chr_)
            .or_insert_with(|| IntervalTree::new())
            .insert(start..end, gene_name);
    }
    result
}

fn graph_from_gtf(file: &str) -> HashMap<String, HashMap<Intervall<i64>, Vec<Intervall<i64>>>> {
    //  TO remove clutter ca use a hashset to remove identical  exon
    // from same genes

    //let file = "genomic.gtf";
    let mut it_dico: HashMap<String, IntervalTree<i64, String>> = gtf_toIT(file.clone());

    let f = File::open(file).unwrap();
    let reader = BufReader::new(f);
    let mut this_line: String; //::new();

    let mut start: i64; //; = 0;
    let mut chr_: String; // = "".to_string();
    let mut end: i64; // = 0;
    let mut gene_name: String = "".to_string();
    let mut transcript_id: String; // = "".to_string();

    let mut seen: HashSet<(String, i64, i64, String)> = HashSet::new();

    let mut g: HashMap<String, HashMap<Intervall<i64>, Vec<Intervall<i64>>>> = HashMap::new();

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
            print!("WARNING Cannort find {:?}", this_line);
            continue;
        }
        if seen.contains(&(chr_.clone(), start, end, gene_name.clone())) {
            continue;
        } else {
            seen.insert((chr_.clone(), start, end, gene_name.clone()));
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
                        .or_insert_with(|| Vec::new())
                        .push(Intervall::new(inter.interval().start, inter.interval().end));
                }
            }
        }
    }
    g
}

fn get_invalid_pos(file: &str) -> HashSet<(String, i64)> {
    let g = graph_from_gtf(file.clone());
    let mut seen = HashSet::new();

    let mut inter_vec = Vec::new();
    let mut e1 = 0;
    let mut borne: &Intervall<i64>;
    let mut results: HashSet<(String, i64)> = HashSet::new();
    let mut file_: Vec<Intervall<i64>> = Vec::new();

    for (chr_, subdict) in g.iter() {
        for (current_node, neigbors) in subdict.iter() {
            inter_vec.clear();
            if seen.contains(current_node) {
                continue;
            }
            // get the CC
            file_.push(current_node.clone());
            while let Some(bfs_node) = file_.pop() {
                if seen.contains(&bfs_node) {
                    continue;
                }
                seen.insert(bfs_node);
                inter_vec.push(bfs_node);
                for n in neigbors {
                    if seen.contains(n) {
                        continue;
                    }
                    file_.push(n.clone())
                }
            }

            // here find "bad pos"
            // here find "bad pos"
            e1 = inter_vec[0].start;
            if !(inter_vec.iter().all(|x| x.start == e1)) {
                //println!("{} {:?}", chr_, inter_vec );
                borne = inter_vec.iter().min_by_key(|x| x.start).unwrap();
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
                borne = inter_vec.iter().max_by_key(|x| x.end).unwrap();
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
        }
    }
    results
}


#[cfg(test)]
mod tests_it {
    use super::*;

    #[test]
    fn test_1() {
        let junction = 
        get_junction_from_gtf("/lab/solexa_yamashita/people/Romain/Projets/Adrienne/OmniSplice/X_auto/gtr_test_1.gtf");
        println!("{:?}", junction);
        println!("{:?}", junction.get("gene-LOC27208819"));
    }

    #[test]
    fn test_2() {
        let junction = 
        get_junction_from_gtf("/lab/solexa_yamashita/people/Romain/Projets/Adrienne/Reference/refseq/Simulans_NCBI/data/GCF_016746395.2/genomic.gtf");
        println!("{:?}", junction.get("gene-LOC27208819"));
        println!("{:?}", junction.get("LOC27208819"));
    }
}
