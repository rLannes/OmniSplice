use crate::common::utils::Exon;
use crate::common::error::OmniError;
use CigarParser::cigar::{Cigar, CigarError};
use bio::data_structures::interval_tree::IntervalTree;
use itertools::Itertools;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::{File, OpenOptions};
use std::hash::Hash;
use std::io::BufRead;
use std::io::BufReader;
use strand_specifier_lib::Strand;
use strand_specifier_lib::{LibType, check_flag};


pub fn get_attr_id(attr: &str, toget: &str) -> Option<String> {
    let mut result: String; 
    let x = attr.split(';').collect::<Vec<&str>>();
    for e in x {
        let spl = e.trim().split(' ').collect::<Vec<&str>>();
        if spl[0].trim() == toget {
            result = spl[1].trim().trim_matches('\"').to_string();
            return Some(result);
        }
    }
    None
}

pub fn gtf_to_hashmap(
    gtf_file: &str,
) -> Result<HashMap<String, HashMap<String, Vec<Exon>>>, OmniError> {
    // gene -> transcript -> Vec<Exon>
    let mut results: HashMap<String, HashMap<String, Vec<Exon>>> = HashMap::new();

    let f = File::open(gtf_file)?;
    let reader = BufReader::new(f);
    let mut this_line: String; 

    let mut chr_: String; 
    let mut start: i64; 
    let mut end: i64; 
    let mut strand: Strand;

    let mut gene_name: String; 
    let mut transcript_id: String; 
    let mut spt: Vec<&str> = Vec::new();

    for line in reader.lines() {
        this_line = line?;
        spt = this_line.trim().split('\t').collect::<Vec<&str>>();

        if spt.len() < 8 {
            continue;
        }
        if spt[2] != "exon" {
            continue;
        }

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
            continue;
        }

        chr_ = spt[0].to_string();
        start = spt[3].parse::<i64>()? - 1; // GTF are one based but BED are 0 based. eg. [3,7] -> [2,7[
        end = spt[4].parse::<i64>()?;
        strand = Strand::from(spt[6]);

        results
            .entry(gene_name)
            .or_insert_with(|| HashMap::new())
            .entry(transcript_id)
            .or_insert_with(|| Vec::new())
            .push(Exon {
                start: start,
                end: end,
                strand: strand,
                contig: chr_.clone(),
            });
    }
    Ok(results)
}

//TODO likely will need refactor
pub fn get_junction_from_gtf(file: &str, libtype: &LibType) -> Result<HashMap<(String, i64, i64), Strand>, OmniError> {
    let f = File::open(file)?;
    let reader = BufReader::new(f);
    let mut this_line: String;

    let mut start: i64;
    let mut chr_: String = "".to_string();
    let mut end: i64;
    let mut strand: Strand = Strand::Plus;
    let mut gene_name: String = "".to_string();
    let mut transcript_id: String;

    let mut transcript_vec: Vec<(i64, i64)> = Vec::new();
    let mut current_gene_id = "".to_string();
    let mut current_transcript_id = "".to_string();
    let mut myset: HashSet<(i64, i64, Strand)> = HashSet::new();

    //let mut result: HashMap<String, HashMap<(String, i64, i64), Strand>> = HashMap::new();
    let mut my_map: HashMap<(String, i64, i64), Strand> = HashMap::new();

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
        end = spt[4].parse::<i64>()?;
        strand = Strand::from(spt[6]);

        strand = match libtype {
            LibType::PairedUnstranded | LibType::Unstranded => Strand::NA,
            _ => Strand::from(spt[6]),
        };

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

        // another approach!
        if gene_name != current_gene_id {
            let ll = transcript_vec.len();
            transcript_vec.sort_by_key(|x| x.0);
            if ll > 1 {
                for i in 0..(ll - 1) {
                    my_map.insert(
                        (chr_.clone(), transcript_vec[i].1, transcript_vec[i + 1].0),
                        strand,
                    );
                    //myset.insert((transcript_vec[i].1, transcript_vec[i + 1].0, strand));
                }
            }

            transcript_vec.clear();

            current_gene_id = gene_name.clone();
            current_transcript_id = transcript_id;
        } else if transcript_id != current_transcript_id {
            let ll = transcript_vec.len();
            transcript_vec.sort_by_key(|x| x.0);
            if ll > 1 {
                for i in 0..(ll - 1) {
                    my_map.insert(
                        (chr_.clone(), transcript_vec[i].1, transcript_vec[i + 1].0),
                        strand,
                    );
                    //myset.insert((transcript_vec[i].1, transcript_vec[i + 1].0, strand));
                }
            }

            transcript_vec.clear();
            current_transcript_id = transcript_id;
        }

        transcript_vec.push((start, end));
    }

    let ll = transcript_vec.len();
    transcript_vec.sort_by_key(|x| x.0);
    if ll > 1 {
        for i in 0..(ll - 1) {
            my_map.insert(
                (chr_.clone(), transcript_vec[i].1, transcript_vec[i + 1].0),
                strand,
            );
            //myset.insert((transcript_vec[i].1, transcript_vec[i + 1].0, strand));
        }
    }

    Ok(my_map)
}

pub fn get_all_junction_for_a_gene(
    gtf_map: &HashMap<String, HashMap<String, Vec<Exon>>>,
) -> Result<HashMap<String, HashSet<(i64, i64)>>, OmniError> {
    let mut res: HashMap<String, HashSet<(i64, i64)>> = HashMap::new();
    let mut tmp = Vec::new();
    for (gene, tr_dict) in gtf_map {
        res.insert(gene.to_string(), HashSet::new());
        for (tr_id, exon) in tr_dict {
            for e in exon {
                tmp.push((e.start, e.end));
            }
            if tmp.len() > 1 {
                tmp.sort_by(|x, y| x.cmp(&y));
                for i in 0..=(tmp.len() - 2) {



                    res.get_mut(gene)
                        .ok_or_else(|| OmniError::Expect("gene unfound in hashmap, cannot recover from this".to_string()))?
                        .insert((tmp[i].1, tmp[i + 1].0));

                }
            }
            tmp.clear();
        }
    }
    Ok(res)
}

// may need to refactor this was part of an older moduler keep it because it work but could crefactor to improve code redibility and performance
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

fn gtf_to_it(file: &str) -> Result<HashMap<String, IntervalTree<i64, String>>, OmniError> {
    //let file = "genomic.gtf";
    let f = File::open(file)?;
    let reader = BufReader::new(f);
    let mut this_line: String; //::new();

    let mut start: i64; //; = 0;
    let mut chr_: String; // = "".to_string();
    let mut end: i64; // = 0;
    let mut gene_name: String = "".to_string();
    let mut transcript_id: String; // = "".to_string();

    let mut result: HashMap<String, IntervalTree<i64, String>> = HashMap::new();

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
        end = spt[4].parse::<i64>()?;

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
    Ok(result)
}

fn graph_from_gtf(file: &str) -> Result<HashMap<String, HashMap<Intervall<i64>, HashSet<Intervall<i64>>>>, OmniError> {
    //  TO remove clutter ca use a hashset to remove identical  exon
    // from same genes

    //let file = "genomic.gtf";
    let mut it_dico: HashMap<String, IntervalTree<i64, String>> = gtf_to_it(file)?;

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

        start = spt[3].parse::<i64>()? - 1; // 1 to 0 based 
        end = spt[4].parse::<i64>()?;

        if start == end {
            continue;
        }

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

pub fn get_invalid_pos(file: &str) -> Result<HashSet<(String, i64)>, OmniError> {
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
                borne = inter_vec.iter().min_by_key(|x| x.start)
                    .ok_or( OmniError::Expect("get invalid pos expecting value".to_string()))?;

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
                borne = inter_vec.iter().max_by_key(|x| x.end)
                .ok_or( OmniError::Expect("get invalid pos expecting value".to_string()))?;
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
