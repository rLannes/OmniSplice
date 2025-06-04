import argparse
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from pathlib import Path
import re
from common_python.junction_class import Junction, parse_js_file
import logging

mpl.rcParams["axes.spines.right"] = False
mpl.rcParams["axes.spines.top"] = False

logging.basicConfig(format='%(asctime)s %(levelname)s => %(message)s', level=logging.INFO, datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger()


def get_attr(string):
    dico = {}
    spt = string.rstrip(";").split(";")
    
    for x in spt:
        if x:
            try:
                dico[x.split()[0].strip()] =  x.split()[1].replace('"', "").strip()
            except:
                print(x, x.split())
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



def ax_plot_bar(ax, counts_, genotype, defect_index, defect_to_plot, colors, reverse=False):


    counts = np.array(counts_)
    masks = [defect_index[defect] for defect in defect_to_plot]

    counts = counts[:, masks]
    sum_ = np.sum(np.array(counts), axis=1).reshape(counts.shape[0])
    prop = counts /  sum_[:, np.newaxis]  


    bottom = np.zeros(len(counts_))
    for i in range(len(masks)):

        ax.barh(y=range(len(counts_)), width=prop[:,i], left=bottom, color = colors[i])
        bottom += prop[:,i]
    
    ax.set_xlabel(genotype)

    if reverse:
        ax.xaxis.set_inverted(True)
    
    return ax


def plot_1(out, defect_index, colors, counts_intron, order, defect_to_plot, width, height, title):

    counts_intron.sort(key=lambda x: x[1])
    counts, intron = list(zip(*counts_intron))
    fig = plt.figure(figsize=(width, height))
    ax1 = fig.add_axes((0, 0, 1, 1))
    ax1 = ax_plot_bar(ax1, [x[order[0]] for x in counts], order[0], defect_index, defect_to_plot, colors, reverse=False)
    ax1.set(ylabel = "introns")
    ax1.set_xticks([])

    custom_lines = [Line2D([0], [0], color=colors[indice], lw=4) for indice, labels in enumerate(defect_to_plot)]
    ax1.legend(custom_lines, defect_to_plot, bbox_to_anchor=(-0.18 , 0.9))

    fig.suptitle(title, y = 1.02)
    plt.savefig(out, bbox_inches="tight")

def plot_2(out, defect_index, colors, counts_intron, order, defect_to_plot, width, height, title):

    counts_intron.sort(key=lambda x: x[1])
    counts, intron = list(zip(*counts_intron))
    

    fig = plt.figure(figsize=(width, height))
    ax1 = fig.add_axes((0, 0, 0.48, 1))
    ax2 = fig.add_axes((0.52, 0, 0.48, 1), sharey=ax1)

    ax1 = ax_plot_bar(ax1, [x[order[0]] for x in counts], order[0], defect_index, defect_to_plot, colors, reverse=False)
    ax2 = ax_plot_bar(ax2, [x[order[1]] for x in counts], order[1], defect_index, defect_to_plot, colors, reverse=True)

    ax1.set(ylabel = "introns")
    ax1.set_xticks([])
    ax2.set_xticks([])
    ax2.yaxis.set_tick_params(labelleft=False, width=0, length=0, left=False)
    ax2.spines['left'].set_visible(False)

    custom_lines = [Line2D([0], [0], color=colors[indice], lw=4) for indice, labels in enumerate(defect_to_plot)]
    ax1.legend(custom_lines, defect_to_plot, bbox_to_anchor=(-0.18 , 0.9))

    fig.suptitle(title, y = 1.02)
    plt.savefig(out, bbox_inches="tight")


def plot_3(out, defect_index, colors, counts_intron, order, defect_to_plot, width, height, title):
    
    fig = plt.figure(figsize=(width, height))

    counts_intron.sort(key=lambda x: x[1])
    counts, intron = list(zip(*counts_intron))
    #(left, bottom, width, height)
    space = 0.05 * 2 
    h = 1 - space
    axes_h = h / 3

    left = 0
    ax1 = fig.add_axes((left, 0, axes_h, 1))
    left += axes_h + 0.05
    ax2 = fig.add_axes((left, 0, axes_h, 1), sharey=ax1)
    left += axes_h + 0.05
    ax3 = fig.add_axes((left, 0, axes_h, 1), sharey=ax1)

    ax1 = ax_plot_bar(ax1, [x[order[0]] for x in counts], order[0], defect_index, defect_to_plot, colors, reverse=False)
    ax2 = ax_plot_bar(ax2, [x[order[1]] for x in counts], order[1], defect_index, defect_to_plot, colors, reverse=False)
    ax3 = ax_plot_bar(ax3, [x[order[2]] for x in counts], order[2], defect_index, defect_to_plot, colors, reverse=False)

    ax1.set(ylabel = "introns")
    ax1.set_xticks([])
    ax2.set_xticks([])
    ax3.set_xticks([])


    #ax2.spines['right'].set_visible(False)
    #ax2.spines['top'].set_visible(False)

    ax3.yaxis.set_tick_params(labelleft=False, width=0, length=0, left=False)
    ax2.yaxis.set_tick_params(labelleft=False, width=0, length=0, left=False)

    ax2.spines['left'].set_visible(False)
    ax3.spines['left'].set_visible(False)

    custom_lines = [Line2D([0], [0], color=colors[indice], lw=4) for indice, labels in enumerate(defect_to_plot)]
    ax1.legend(custom_lines, defect_to_plot, bbox_to_anchor=(-0.18 , 0.9))
    fig.suptitle(title, y = 1.02)
    plt.savefig(out, bbox_inches="tight")


def checkinput_args(args):
    try:
        assert args.gene_list or args.transcript_list and not (args.gene_list and args.transcript_list)
    except AssertionError:
        logging.error("either one or the other but not both arguments [--gene_list] [--transcript_list] should be set")
        raise

    try:
        assert not args.group3 or (args.group3 and (args.group1 and args.group2))
    except AssertionError:
        logging.error("group3 should only be set if group 1 and 2 are set")
        raise
    
    try: 
        assert not args.group2 or (args.group1 and args.group2)
    except AssertionError:
        logging.error("group2 should only be set if group 1 is set")

    try:
        assert len(args.color) >= len(args.splicing_defect)
    except:
        logging.error("you shoudl provides more color than splicing defect")
        raise


if __name__ == "__main__":
    parse = argparse.ArgumentParser(description="to plot junction defect from Omnisplice junction file accept from one up to three genotypes." \
    " it merge and sums file from the same genotype." \
    " you can either submit file using glob or space seprate them."
    "example usage:"
    "python omnisplice/gene_visualisation.py --gtf <gtf file> --group1 genotype1*junctions --group1_name <genotype1 name> --gene_list FBgn0001313 --outfile_prefix my_plot") # Can do more but will need a some works to make it able to handle X genotypes
    parse.add_argument("--gtf", required=True)
    parse.add_argument("--gtf_gene_id", default="gene_id")
    
    parse.add_argument("--group1", nargs="+", required=True, help="junctions file for group 1")
    parse.add_argument("--group2", nargs="+", required=False, help="junctions file for group 2")
    parse.add_argument("--group3", nargs="+", required=False, help="Optional junctions file for group 3")

    parse.add_argument("--group1_name", required=False, default='group1', help="name group 1")
    parse.add_argument("--group2_name", required=False, default='group2', help="name for group 2")
    parse.add_argument("--group3_name", required=False, default='group3', help="Optional name group 3")
    parse.add_argument("--color", nargs="+", default=["#D3D3D385", "#214d4e",  "#cc655b", '#41bbc5',
                                                      "#c6dbae", '#069668',  '#b3e61c',  '#d55e00', '#cc78bc',
                                                      '#ca9161', '#fbafe4',  '#029e73'])
    
    parse.add_argument("--gene_list", nargs="+", required=False, help='gene list to plot')
    parse.add_argument("--transcript_list", nargs="+", required=False, help='gene list to plot')
    parse.add_argument("--outfile_prefix", required=True, help='basename for output')

    parse.add_argument("--splicing_defect", nargs="+", help="""default = SPLICED, UNSPLICED, CLIPPED, EXON_OTHER""", 
                        default=["SPLICED" , "UNSPLICED",  "CLIPPED",  "EXON_OTHER"],
                        choices=["SPLICED", 
                                 "UNSPLICED", 
                                 "CLIPPED", 
                                 "EXON_OTHER",
                                 "SKIPPED",
                                 "WRONG_STRAND",
                                 "E_ISOFORM"] )

    parse.add_argument("--format",  default='.pdf', help='output format default: ".pdf"')
    parse.add_argument("--fig_width",  default=6.4, help="")
    parse.add_argument("--fig_height",  default=4.8, help='')
    parse.add_argument("--logging_level", "-v", default="INFO",choices=["DEBUG", "INFO", "ERROR"])

    args = parse.parse_args()


    defect_index = ["SPLICED",      
                    "UNSPLICED", 
                    "CLIPPED", 
                    "EXON_OTHER",
                    "SKIPPED",
                    "WRONG_STRAND",
                    "E_ISOFORM"]

    defect_index = dict((e, i) for (i, e) in enumerate(defect_index))

    if args.logging_level == "DEBUG":
        logger.setLevel(logging.DEBUG)
    elif args.logging_level == "INFO": 
        logger.setLevel(logging.INFO)
    elif args.logging_level == "ERROR":
        logger.setLevel(logging.ERROR)

    checkinput_args(args)


    dico_gtf = gtf_to_dict(args.gtf, args.gtf_gene_id)
    dico_result = {}

    for file in args.group1:
        parse_js_file(file=file, results=dico_result, genotype=args.group1_name, ambigious=True, gene_list=args.gene_list, transcript_list=args.transcript_list)
    for file in args.group2:
        parse_js_file(file=file, results=dico_result, genotype=args.group2_name, ambigious=True, gene_list=args.gene_list, transcript_list=args.transcript_list)
    if args.group3:
            logging.debug("group3")
            for file in args.group3:
                parse_js_file(file=file, results=dico_result, genotype=args.group3_name, ambigious=True, gene_list=args.gene_list, transcript_list=args.transcript_list)


    dico_j = {}

    for hashkey, junction in dico_result.items():
    
        transcript_list = junction.transcript
        gene_list = junction.gene
        intron_list = junction.intron_number

        zipped = zip(gene_list, transcript_list, intron_list)
        for (gene, transcript, intron) in zipped:
            if gene not in dico_j:
                dico_j[gene] = {}
            if transcript not in dico_j[gene]:
                dico_j[gene][transcript] = []

            dico_j[gene][transcript].append((junction.get_count_per_genotype_summed(), intron))
     
    # build legend



    for gene_id, tr_dico in dico_j.items():
        for tr_id, counts_dico in tr_dico.items():
            if not args.group2:
                logging.debug(" 1 samples ")
                order = [args.group1_name]
                plot_1(out="{}_{}_{}{}".format(args.outfile_prefix, gene_id, tr_id, args.format), defect_index=defect_index, colors=args.color, 
                       counts_intron=counts_dico, order=order, defect_to_plot=args.splicing_defect, width=args.fig_width, height=args.fig_height, title='{} {}'.format( gene_id, tr_id))
            elif  not args.group3:
                logging.debug(" 2 samples ")
                order = [args.group1_name, args.group2_name]
                plot_2(out="{}_{}_{}{}".format(args.outfile_prefix, gene_id, tr_id, args.format), defect_index=defect_index, colors=args.color, 
                       counts_intron=counts_dico, order=order, defect_to_plot=args.splicing_defect, width=args.fig_width, height=args.fig_height, title='{} {}'.format( gene_id, tr_id))
                # plot2
            elif args.group3 :
                logging.debug(" 3 samples ")
                order = [args.group1_name, args.group2_name,  args.group3_name]
                plot_3(out="{}_{}_{}{}".format(args.outfile_prefix, gene_id, tr_id, args.format), defect_index=defect_index, colors=args.color, 
                       counts_intron=counts_dico, order=order, defect_to_plot=args.splicing_defect, width=args.fig_width, height=args.fig_height, title='{} {}'.format( gene_id, tr_id))
            else:
                pass
            plt.close()


