import argparse
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from pathlib import Path
#import statsmodels.api as sm
#import statsmodels.formula.api as smf
import re


mpl.rcParams["axes.spines.right"] = False
mpl.rcParams["axes.spines.top"] = False

"""
def data_to_df_for_glm(control, treatment):
    df = pd.DataFrame(control + treatment)
    df["failures"] = df.apply(lambda x: sum(x[1:]), axis=1)
    df["successes"] = df[0]
    df["group"] = ["control"] * len(control) + ["treatment"] * len(treatment)
    return df


def gln_binomial_test(control, treatment):
    df = data_to_df_for_glm(control, treatment)
    mod = smf.glm('successes + failures ~ group', family=sm.families.Binomial(), data=df).fit()
    # convert intercept to probability
    odds_i = math.exp(mod.params[0])
    control_p = odds_i / (1 + odds_i)

    # convert stable="stable" to probability
    odds_p = math.exp( mod.params[0]) * math.exp(mod.params[1])
    treatment_p = (odds_p / (1 + odds_p))
    return (mod.pvalues[1], control_p , treatment_p)
"""

def get_attr(string):
    dico = {}
    spt = string.split(";")
    
    for x in spt:
        if x:
            dico[x.split()[0].strip()] =  x.split()[1].replace('"', "").strip()
    return dico


def gtf_to_dict(gtf_file, main_ = "gene_id"):
    dico = {}
    with open(gtf_file) as f_in:
        """atm reads only genes"""
        for line in f_in:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue
            spt = line.strip().split("\t")
            try:
                chr_ = spt[0]
                start = int(spt[3]) #- 1 # gtf are 1 based
                end = int(spt[4]) #-1 
                strand =  spt[6]
                attr = get_attr(spt[-1])
                """                 if main_ not in attr:
                                    print(attr)
                continue """
                try:
                    gene_id = attr[main_]
                except:
                    print("main not found")
                    print(attr)
                    raise

                gene_symbol = attr.get("gene_symbol", None)
                if not gene_symbol:
                    gene_symbol = attr.get("gene_name", None)
                if not gene_symbol:
                    gene_symbol = "None"

                type_ = spt[2]
            except:
                print(attr)
                print(line)
                print(spt)
                raise
            
            if gene_id not in dico:
                dico[gene_id] = {
                    "chr": chr_,
                    "symbol" :  gene_symbol,
                    "strand" : strand,
                    "transcript": {}
                }
                
            if type_ == "gene":
    
                dico[gene_id]["start"] = start
                dico[gene_id]["end"] = end
            
            else:
                try:
                    transcript_id  = attr.get("transcript_id", "None")
                    transcript_symbol = attr.get("transcript_symbol", "None")

                    # gene_id  = attr["gene_id"]
                    # gene_symbol = attr.get("gene_symbol", "None")
                    # strand =  spt[6]
                    
                    if transcript_id not in dico[gene_id]["transcript"]:
                        dico[gene_id]["transcript"][transcript_id] = {
                        "transcript_symbol" : transcript_symbol,
                        "transcript_id" : transcript_id
                        }
                    
                    if type_ not in dico[gene_id]["transcript"][transcript_id]: 
                        dico[gene_id]["transcript"][transcript_id][type_] = []
                        
                    dico[gene_id]["transcript"][transcript_id][type_].append({
                        "chr": chr_,
                        "start": start, 
                        "end" : end,
                        "strand" : strand,
                    })
                except: 
                    continue
    return dico


def collapse(values):
    res = []
    for i, v in enumerate(values):
        if i % 2 == 0:
            res.append([x+values[i+1][ii] for ii, x in enumerate(v)])
    return res


def parse_cat_file(genotype, files, dico_result, gene_list):

    gene_set = set(gene_list)
    for file in files:
        geno = genotype
        
        with open(file) as fi:
            fi.readline()
            for l in fi:
                l = l.strip()
                if not l:
                    continue
                spt = l.split('\t')
                
                value = list(map(int, spt[9: 9 +5]))
                
                gene = spt[1]
                transcript = spt[2]
                if gene not in gene_set:
                    continue
                if gene not in dico_result:
                    dico_result[gene] = {}
                if transcript not in dico_result[gene]:
                    dico_result[gene][transcript] = {}
                if spt[6] == "Acceptor":
                    continue
                exon_n = spt[3]
                if exon_n not in dico_result[gene][transcript] :
                    dico_result[gene][transcript][exon_n] = {}
                
                if geno not in dico_result[gene][transcript][exon_n]:
                    dico_result[gene][transcript][exon_n][geno] = []
                dico_result[gene][transcript][exon_n][geno].append(value)


def plot(dico_result, color, out_base, format):
    for gene, sub_gene in dico_result.items():
        #gene_name = dico_gtf.get(gene).get("symbol")
        gene_name = gene
        for transcript_id, tr_dico in sub_gene.items(): 

            len_exon = len(dico_result[gene][transcript_id])
            z = 0
            
            exon_order = list(sorted(dico_result[gene][transcript_id].keys(), key=lambda x: x.split('_')[0], reverse=True))[:-1]
            for i_exon, exon in enumerate(exon_order):

                group1 = dico_result[gene][transcript_id][exon]["group1"]
                control = group1
                group1 = np.sum(np.array(group1), axis=0)
                s_group1 = sum(group1) 
                group1 = group1 / s_group1

                group2 = dico_result[gene][transcript_id][exon]["group2"]
                treatment = group2
                group2 = np.sum(np.array(group2), axis=0)
                s_group2 = sum(group2) 
                group2 = group2 / s_group2

                this = exon_order[i_exon]
                this = this.split('_')
                this[1] = str(int(this[1]))
                exon_order[i_exon] = " ".join(this)
                """
                try:
                    pval, control_p, treat_p = gln_binomial_test(control, treatment)
                    #exon_order[i_exon] += "_ctrlp={}_treatp={}".format(round(control_p,2), round(treat_p, 2))
                    if pval < 0.001:
                        
                        exon_order[i_exon] += " *"#.format(round(-math.log(pval, 10), 2))
                    #print(i_exon, pval)    
                except:
                    #print(control, treatment)
                    exon_order[i_exon] += "_na"
                """
                bottom = 0
                #print(group1)
                #print(group2)
                for i,e in enumerate(group1):
                    plt.barh(y=z, width=e,left=bottom, color = color[i])
                    bottom += e
                bottom += 0.05
                for i,e in enumerate(group2[::-1]):
                    plt.barh(y=z, width=e,left=bottom, color = color[4-i])
                    bottom += e
        
                z += 1
            plt.yticks(ticks=range(0,z), labels=exon_order)
            handles, labels = plt.gca().get_legend_handles_labels()
            line1 = Line2D([], [], label="Spliced",   color=color[0], linewidth=6)
            line2 = Line2D([], [], label="Unspliced",  color=color[1], linewidth=6)
            line3 = Line2D([], [], label="Clipped",   color=color[2], linewidth=6)
            line4 = Line2D([], [], label="Exon_Intron",   color=color[3], linewidth=6)
            line5 = Line2D([], [], label="Exon_other",   color=color[4], linewidth=6)
            
            handles.extend([ line1, line2, line3, line4, line5])
            plt.gcf().legend(handles=handles, loc='outside right upper', bbox_to_anchor=(1.25, 0.8))
            plt.gca().set_xticks([0.5, 1.5], labels=['group1', "group2"])
            plt.title(gene_name,     fontsize=24)
            plt.xlim(-0.1, 2.2)

            plt.gca().set_facecolor('aliceblue')
            plt.tick_params(axis='x', labelsize=20)
            plt.tick_params(axis='y', labelsize=12)
            
            plt.savefig("{}_{}_{}.{}".format(out_base, gene_name, transcript_id, format), bbox_inches="tight")
            plt.close()


if __name__ == "__main__":
    parse = argparse.ArgumentParser()
    parse.add_argument("--gtf", required=True)
    parse.add_argument("--gtf_gene_id", default="gene_id")
    parse.add_argument("--group1", nargs="+", required=True, help="table file for group 1")
    parse.add_argument("--group2", nargs="+", required=True, help="table file for group 2")
    parse.add_argument("--color", nargs="+")
    
    parse.add_argument("--gene_list", nargs="+", required=True, help='gene list to plot')
    parse.add_argument("--outfile_prefix", required=True, help='basename for output')
    parse.add_argument("--format",  default='.pdf', help='output format default: ".pdf"')
    
    args = parse.parse_args()

    dico_gtf = gtf_to_dict(args.gtf, args.gtf_gene_id)
    dico_result = {}
    parse_cat_file("group1", args.group1, dico_result, args.gene_list)
    parse_cat_file("group2", args.group2, dico_result, args.gene_list)
    color = ["#D3D3D385", "#214d4e",  "#cc655b", '#41bbc5', "#c6dbae", '#069668',  '#b3e61c',  '#d55e00', '#cc78bc',
         '#ca9161', '#fbafe4',  '#029e73']
    if args.color:
        color = args.color
    plot(dico_result, color, args.outfile_prefix, args.format)
    

