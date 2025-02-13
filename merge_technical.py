import argparse
from pathlib import Path
import os
import re 
import sys
import polars as pl

def get_file_by_match(dir_path, reg):
        dico = {}
        for file in dir_path.iterdir():
            
            if m := reg.search(str(Path(file).name)):
                key = m.group(1)
                if key not in dico:
                    dico[key] = []
                dico[key].append(file)
        return dico


## refactor using polars
def parse_rep(rep, out):
    dfs = [pl.read_csv(f, separator='\t') for f in rep]
    concat = pl.concat(dfs, how="diagonal")
    concat.group_by(["contig", "gene_name", "transcript_name", "exon_number", "ambiguous", "strand", "pos", "next", "exon_type"])\
        .agg(pl.col(["spliced", "unspliced", "clipped", "exon_intron", "exon_other", "skipped", "wrong_strand", "e_isoform"]).sum())\
        .sort([
                pl.col(["contig", "gene_name", "transcript_name"]),
                pl.col("exon_number").str.split('_').list.get(1).cast(pl.Int64),
                pl.col("exon_type")
            ])\
        .write_csv(out, separator="\t")
    
    #pl.col('name_value').map(lambda x: int(x.split('_')[1]))

#     dico = {}

#     for file in rep:
#         with open(file) as fi:
#             header = fi.readline().strip().split("\t")
#             spliced_i = header.index("spliced") 
#             for l in fi:
#                 spt = l.strip().split("\t")
#                 if not spt:
#                     continue
#                 key = tuple(spt[:-8])
#                 value = list(map(int, spt[-8:]))

#                 if key not in dico:
#                     dico[key] = [0]* len(value)
#                 pr = dico[key]
#                 dico[key] = [value[i] + v for i, v in enumerate(pr)]
#         return header, dico

# def write_output(out, header, dico):
#     try:
#         assert not os.path.isfile(out)
#     except:
#         raise AssertionError("outputfile exist {}".format(out))
#     with open(out, "w") as fo:
#         fo.write("\t".join(header) + "\n")

#         for k, v in sorted(dico.items(), key = lambda x : (x[0][1], x[0][2], int(x[0][3].split("_")[-1]), x[0][6] )):
#             fo.write("\t".join(k) + "\t")
#             fo.write("\t".join(list(map(str, v))) + "\n")






if __name__ == "__main__":
    parse = argparse.ArgumentParser()
    parse.add_argument("--rep", nargs="+", help="""usage, separate replicate you want to merge by coma,
                        you can give multiple job by separating them by a space exmaple:
                       the two control technical ctrl1_1 ctrl1_2
                       biologicql replicate ctrl2_1 ctr2_2
                       python3 merge_techincal.py --rep ctrl1_1,ctrl1_2 ctrl2_1,ctr2_2
                       will merge ctrl1_1 with ctrl1_2 and will merge ctrl2_1 with ctrl2_2 
                       """)
    parse.add_argument("--out", nargs="+", help="""space separated list of output, number must match with the --rep argument,
                       or be a directory if you use the match argument""")
    parse.add_argument("--match", help="use regexp, usefull when multiple groups")
    parse.add_argument("--dir", help="directory too look for match")
    args = parse.parse_args()

    if args.rep:
        print(args.rep)
        assert len(args.rep) == len(args.out)
        

        for i, reps in enumerate(args.rep):
            parse_rep(reps.split(','), args.out[i])
            #header, dico = parse_rep(reps)
            #write_output(out=args.out[i], header=header, dico=dico)

    elif args.match:

        assert len(args.out) == 1
        args.out = args.out[0]
        assert os.path.isdir(args.out) and os.path.isdir(args.dir)
        reg = re.compile(args.match)

        
        dir_path = Path(args.dir)  
        dico = get_file_by_match(dir_path, reg)
        print("match identified the following groups:")
        for k, v in dico.items():
            print("{} -> {}".format(str(k), '\t'.join(list(map(str, v)))))

        out_p = Path(args.out)
        
        for k, reps in dico.items():
            parse_rep(reps, out=out_p / "{}.merged.table.tsv".format(k))
            #write_output(out=out_p / "{}.merged.table.tsv".format(k), header=header, dico=dico_)
