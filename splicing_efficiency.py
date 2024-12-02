from pathlib import Path
import argparse
import os
import re


def iter_gene(open_file):
    current_gene = ""
    genes = []
    for l in open_file:
        spt = l.strip().split()
        this_gene = spt[1]
        
        if current_gene != this_gene:
            if genes:
                yield genes
            genes = []
            current_gene = this_gene
        genes.append(spt)
    yield genes

def filter_gene(gene):
    res = []
    for g in gene:
        if (g[8] != "Donnor") or (g[4] == "true") or (g[7] == "."):
            continue
        res.append(g)
    return res


def parse_file(file):
    dico_res = {}
    with open(file) as f_in:
        for genes in iter_gene(f_in):
            if not genes:
                continue
                
            junctions  = filter_gene(genes)
            if not junctions:
                continue

            dico_res[junctions[0][1]] = {}

            for g in junctions:
                this_pos = int(g[6])
                this_next = int(g[7])

                if this_pos not in dico_res[junctions[0][1]]:
                    dico_res[junctions[0][1]][this_pos] = {}
                dico_res[junctions[0][1]][this_pos][this_next] = [int(x) for x in g[9:]]

    return dico_res


def write_result(out, dico_res):
    with open(out, "w") as f_o:
        for g, v  in dico_res.items():

            for start, end in sorted(v.items(), key =lambda x: x[0]):
            
                spliced = sum([x[0] for x in end.values()])
                unspliced = sum(list(end.values())[0][1:4])
                
                kk = sorted(end.keys())
                #if len(end) > 1:
                ratio = '.'
                if unspliced + spliced > 0:
                    ratio =   spliced / (unspliced + spliced)
                    
                to_write = "{}\t{}\t{}\t{}\t{}\t{}\n".format(g, start, ','.join(list(map(str, kk))), str(spliced), str(unspliced), str(ratio))
                f_o.write(to_write)



if __name__ == "__main__":
    parse = argparse.ArgumentParser()
    parse.add_argument("--input", help="comma separated list of file to parse")
    parse.add_argument("--out", nargs="+", help="""space separated list of output, number rmust match with the --rep argument,
                       or be a directory if you use the match argument""")
    parse.add_argument("--match", help="use regexp")
    parse.add_argument("--dir", help="directory too look for match")
    args = parse.parse_args()


    if args.input:
        print(args.input)
        assert len(args.rep) == len(args.out)

        for i, reps in enumerate(args.input):
            dico = parse_file(reps)
            write_result(out=args.out[i], dico_res=dico)


    elif args.match:

        assert len(args.out) == 1
        args.out = args.out[0]
        assert os.path.isdir(args.out) and os.path.isdir(args.dir)
        reg = re.compile(args.match)
        dir_path = Path(args.dir)  
        dico = {}

        for file in dir_path.iterdir():
            
            if m := reg.search(str(Path(file).name)):
                key = m.group(1)
                if key not in dico:
                    dico[key] = []
                dico[key].append(file)
        out_p = Path(args.out)
        
        for k, reps in dico.items():
            dico = parse_file(reps[0])
            write_result(out=out_p / "{}.SE.table.tsv".format(k), dico_res=dico)



