Welcome to the omnisplice read me.

Omnisplice is a tool to categorize all reads present at exon's extremity.

This file follow the following organisation:

    - Quick uses
    - Installation
    - Input Files
    - Output Files
    - Command details
    - Algoritm overview
    - Context and discussion.



## To keep in mind
The same read will be counted as many time as a read overlap different feature.
even it is the same position (different transcript same gene).


## Quick uses
```
omni_splice --gtf <gtf file> --input <indexBamFile> --output <outFilePrefix>
```


## Installation

to use omni splice you need access to a linux machine (mac never tested).
you need rust, to install rust go to (https://www.rust-lang.org/tools/install). 
```
# Clone this repository using git
git clone <####>\
cd omnisplice\
# this is a compilation command it will output a bunch of text including some wanring
# and should finish with: "Finished `release` profile [optimized] target(s) "
cargo build --release
```

omnisplice executables will then be at <path>/omnisplice/target/release/omni_splice\
if you want to add it to your Path, you can add the realease directory to your ~/.bashrc.

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
    - input: a valid position sorted indexed bam file

## Output Files
    category file (tsv file): 
        Chromosome | Position of the feature | gene ID | transcript ID | Strand | ExonType | ReadCategory | Read Support
        for each feature they are as many line as ReadCategory

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

    the table file:


    the backsplicing file:
    


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

    

## Algoritm overview


For each read find all exon extremity it does overlap.

for each read R:
    for each extremity E:

        Test if the strand match. -> if False Wrong Strand
        Test if the read skipped the junction.
        Tesf if the read fail the overhang.
        Test if the read is unspliced at the feature specifically.
        Test if the read is SoftClipped.
        // changing put it first!
        Test if fail pos filter.
        Test if the read is a Junction read.


## Context and discussion.
    OmniSplice surprisingly does not consume much ressource. Execution time will depend on the size of the bam and memory consumption will depend on the number of "exon" Feature in the gtf .


