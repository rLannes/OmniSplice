Welcome to the omnisplice read me.

Omnisplice is a tool to categorize all reads present at exon's extremity.
get more info reading the paper: https://www.biorxiv.org/content/10.1101/2025.04.06.647416v1

This file follow the following organisation:

    - Quick uses
    - Installation
    - Input Files
    - Command details
    - Output Files
    - Context and discussion.
    - Algorithm overview


## To keep in mind
The same read will be counted as many time as a read overlap different feature.
even it is the same position (different transcript same gene).


## Quick uses
```
omni_splice --gtf <gtf file> --input <indexBamFile> -o <outFilePrefix>

# specify the lib type for unstranded data
omni_splice --gtf <gtf file> --input <indexBamFile> -o <outFilePrefix> --Libtype Unstranded

# extending the overhang to 5 pb
omni_splice --gtf <gtf file> --input <indexBamFile> -o <outFilePrefix> --overhang 5 

# for backsplicing you need to first extract the clipped reads
omni_splice --gtf <gtf file> --input <indexBamFile> -o <outFilePrefix> --read-to-write soft-clipped

# then you need bowtie 2 in your OS path, plus a bowtie 2 reference for your genome
backsplicing -i <omnisplice_out.clipped> -o <outputPrefix>> -b <bowtie2 ref> -g <gtf> -m <min clipped size, default 20> 

# for comparison (Fischer test, GLM is available but we are still testing it)
# select defect as appropriate for your study
python omnisplice/compare_conditions.py --control control1.junction control2.junction --treatment treatment1.junction treatment2.junction --stat FISCHER --defect  SPLICED UNSPLICED CLIPPED EXON_OTHER SKIPPED  --out <out file>


# quick plot of genes (same as in manuscript)
coming out soon!
```




## Installation

to use omni splice you need access to a linux machine (mac never tested).
you need rust, to install rust go to (https://www.rust-lang.org/tools/install). 
you also need cmake. (brew install cmake on mac.)

```
# Clone this repository using git
> git clone <####>
> cd omnisplice
> cargo build --release
# this is a compilation command it will output a bunch of text including some wanring
# and should finish with: "Finished `release` profile [optimized] target(s) "
```

omnisplice executables will then be at <path>/omnisplice/target/release/omni_splice
if you want to add it to your Path, you can add the release directory to your ~/.bashrc.

to do so 
```
> cd target/release/
> pwd
# copy the results of the pwd command.
> vim  ~/.bashrc # on mac vim ~/.zshrc

# 
```

if you encouter any error, have any questions, want to propose improvment
or let us now by opening a new issue on this github.


### Rust installation
    To see if you have rust installed type "cargo" on the command line, if the terminal return a line with "command not found" you need to install it.
    It is very easy just follow this link: 
        https://www.rust-lang.org/tools/install
    
    if you have rust installed you may need to update it:   
        rustup update


## Input Files
    - gtf: The gtf must have a valid "gene_id" and "transcript_id" for every feature annotated as "exon".
    omnisplice will look at all the exon in this file. if you want to limit your search to subset of gene or exon.
    only include the one your are interested in. (this is particularly helpfull when using the "readToWrite" option as to limit the size of the output).
    - input: a valid position sorted indexed **bam file**


## Command details
    required:
        --input -i  <sortedIndexedBamFile>
        --output -o <outFilePrefix>
        --gtf -g <gtfFile>
    optional: 
        --mapq <Minimum mapq default 13>
        --LibType <Rna seq libtype, default frFirstStrand>
        --overhang <default 1>
        --flag_in <Bam flag a read must have>
        --flag_out <Bam flag a read must not have>
        --readToWrite <used to output read of a specific category>
        --unspliced_def <used for splicing efficiency default unspliced>
        --spliced_def <used for splicing efficiency default spliced>

```
❯ omni_splice -h
Usage: 
omni_splice [OPTIONS] --input <INPUT> --output-file-prefix <OUTPUT_FILE_PREFIX> --gtf <GTF>

Options:
  -i, --input <INPUT>
          Name of Input file
  -o, --output-file-prefix <OUTPUT_FILE_PREFIX>
          Prefix name  to be used for Output file
  -g, --gtf <GTF>
          Name of GTF Input file define the feature to look at (v1) only consider 
          feature annotated as exon if you use output_write_read with the whole genome
          the output can be very large, you may want to subset genes / features you are
          interested in
      --overhang <OVERHANG>
          size of overhang [default: 1]
      --output-write-read <OUTPUT_WRITE_READ>
          path to a file (must not exist) Uses if you want to output the reads with their category.
          by default output all reads, this behaviour can be change using flags: clipped...
      --flag-in <FLAG_IN>
          [default: 0]
      --flag-out <FLAG_OUT>
          [default: 3840]
      --mapq <MAPQ>
          [default: 13]
      --read-to-write <READ_TO_WRITE>...
          space separated list of the annotated read you want to extract; i.e.
          all clipped read or all spliced read ... [possible values: read-through,
          read-junction, unexpected, fail-pos-filter, wrong-strand, fail-qc,
          empty-pileup, skipped, soft-clipped, overhang-fail, empty, all]
      --unspliced-def <UNSPLICED_DEF>...
          space separated list the column to use for "unspliced" for the splicing defect table.
          you can regenrate this using the splicing_efficiency exe What to consider as unspliced?
          unspliced: 10, clipped: 11, exon_intron: 12, exon_other: 13, skipped: 14, wrong_strand:15,
          isoform:16
          by default only use "-u 10" ->  unspliced (readthrough) reads \n to use unspliced and clipped : "-u 10 11" [default: 10]
      --spliced-def <SPLICED_DEF>...
          What to consider as unspliced? unspliced: 10, clipped: 11, exon_intron: 12, exon_other: 13, skipped: 14,
          wrong_strand:15, isoform: 16
          by default only use "-u 9" -> spliced (readthrough) reads \n to use spliced and isform : "-u 9 16" [default: 9]
      --libtype <LIBTYPE>
          Librairy types used for the RNAseq most modern stranded RNAseq are frFirstStrand which is the default value.
          acceptable value: frFirstStrand, frSecondStrand, fFirstStrand, fSecondStrand, ffFirstStrand, ffSecondStrand,
           rfFirstStrand, rfSecondStrand, rFirstStrand, rSecondStrand, Unstranded, PairedUnstranded [default: frFirstStrand]
  -h, --help
          Print help
  -V, --version
          Print version
```

## Context and discussion.
    OmniSplice surprisingly does not consume much ressource. Execution time will depend on the size of the bam and memory consumption will depend on the number of "exon" Feature in the gtf .



## Output Files
    category file (.tsv file): 
        Chromosome | Position of the feature | gene ID | transcript ID | Strand | ExonType | ReadCategory | Read Support
        for each feature they are as many line as ReadCategory found for this specific feature

    with:
        ExonType:   [Donnor | Acceptor] indicate if the exon extremity is a junction donnor of acceptor.
            Note: the start of the gene is Acceptor, the end of the gene is Donnor.
        
        ReadCategory:
            * ReadJunction(n, m) -> junction read from n, m
            * Unspliced -> read is unspliced
            * Skipped(n, m) -> skipped , read has junction from  n to m that include the feature. 
            * SoftClipped -> Read is clipped at the junction.
            * Empty -> not read found at this junction
            * WrongStrand -> The read come from a fragment which strand is different feom the feature
            * FailPosFilter -> the way the algorithm works to reduce the number of comparison may give some false positive. 
                you can discard it, no meaningfull biologicaly.
            * OverhangFail -> The read fail the overhang test.
            * Unexpected -> if you see this please open a bug report.

    both table and junction file describe events. But the table file give the detail by exon end whereas the junction file give the results per junction.

    the table file (.tsv file):
        contig | gene_name | transcript_name | exon_number | ambiguous | strand | pos | next | exon_type | spliced | unspliced | clipped | exon_other | skipped  | wrong_strand | e_isoform
    
    the junction file (.tsv file):
        contig | gene_name | transcript_name | intron_number | strand | ambiguous | Donnor | Acceptor | spliced | unspliced | clipped | exon_other | skipped  | wrong_strand | e_isoform



## Algorithm overview

```
For each read find all exon extremity it does overlap.
for each read R:
    for each extremity E:
        Test if fail pos filter.
        Test if the strand match. -> if False Wrong Strand
        Test if the read skipped the junction.
        Tesf if the read fail the overhang.
        Test if the read is unspliced at the feature specifically.
        Test if the read is SoftClipped.
        Test if the read is a Junction read.
```

