import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path
import math
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
from  statsmodels.stats.multitest import fdrcorrection
import numpy as np
import pandas as pd
import re
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import chi2
from abc import abstractmethod
from enum import Enum
from abc import ABC
from scipy.stats import chisquare

from math import isclose
from scipy.stats import fisher_exact
from fast_fisher import fast_fisher_exact_compatibility
import time

import sys
from common_python.counter_junction import Counter, ReadJunction
from common_python.junction_class import Junction

def timer_decorator(func):
    def inner1(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)

        print('Function', func.__name__, 'time:', time.time() - start)

        return result
    return inner1


class Stat_test(ABC):

    @abstractmethod
    def test_sample(self, control, treatment):
        pass

    @abstractmethod
    def load_data(self, control, treatment):
        pass

    @abstractmethod
    def get_prop(self, val, group):
        pass

    @abstractmethod
    def repr(self):
        pass

    

class Fischer_test(Stat_test):

    def __init__(self): 
        super().__init__()
    
    def repr(self):
        return "FISCHER"
    
    def test_sample(self, control, treatment):
        data = self.load_data(control, treatment)
        (odds_f, pval_f )= fast_fisher_exact_compatibility(data, alternative='two-sided')
        return (pval_f, "Fischer_statistic:{};".format(odds_f) + "\t" + self.get_prop(data[0], "control") + "\t" + self.get_prop(data[1], "treatment"))

   
    def load_data(self, control, treatment):
        cc = [sum([x for x in control["successes"]]), sum([x for x in control["failures"]])]
        tt = [sum([x for x in treatment["successes"]]), sum([x for x in treatment["failures"]])]
        return [cc, tt]


    def get_prop(self, val, group):
        successes = val[0]
        failures = val[1]
        p = -1
        if successes + failures > 0:
            p = successes / (successes + failures)
        
        return "{}:{}/{}/{}".format(group, successes, failures, p)

class GLM_model(Stat_test):

    def __init__(self):
        super().__init__()
    
    def repr(self):
        return "GLM"
    
    
    def test_sample(self, control, treatment):
        
        data = self.load_data(control, treatment)
        try:
            mod = smf.glm('successes + failures ~ group', family=sm.families.Binomial(), data=data).fit()
            return (mod.pvalues[1], self.get_prop(data, "control")  + ";" + self.get_prop(data, "treatment"))
        
        except:
            return ("NA", self.get_prop(data, "control") + ";" + self.get_prop(data, "treatment"))


    def load_data(self, control, treatment):

        try:
            df = pd.DataFrame({"successes":[x for x in control["successes"]] + [x for x in treatment["successes"]],
                        "failures":[x for x in control["failures"]] + [x for x in treatment["failures"]],
                            "group":["control"]*len(control) + ["treatment"]*len(treatment)})
        except:
            print(control_flatten, treatment_flatten, control, treatment)
        return df
    

    def get_prop(self, df, group):

        successes = sum(df[df["group"] == group]["successes"])  
        failures =  sum(df[df["group"] == group]["failures"])
        p = -1

        if successes + failures > 0:
            p = successes / (successes + failures)        
        
        return "{}:{}/{}/{}".format(group, successes, failures, p)



class Chi2(Stat_test):

    def repr(self):
        return "Chi2"
    
    def test_sample(self, control, treatment):
        pass

   
    def load_data(self, control, treatment):
        cc = [sum([x[0] for x in control]), sum([x[1] for x in control])]
        tt = [sum([x[0] for x in treatment]), sum([x[1] for x in treatment])]
        return [cc, tt]

    def get_prop(self, df, group):
        pass
"""

class Junction():

    def __init__(self, spt, header):
        self.contig = spt[header["contig"]]
        self.gene = spt[header["gene_name"]]
        self.transcript = spt[header["transcript_name"]]
        self.ambigious = True if spt[header["ambiguous"]] == "true" else False
        exon_number = int(re.search("(\d+)", spt[header["exon_number"]]).group(1))
        self.intron_number = exon_number  if spt[header["exon_type"]] == "Donnor" else exon_number - 1
        self.strand = spt[header["strand"]]
        self.pos =  spt[header["pos"]] if spt[header["exon_type"]] == "Donnor" else spt[header["next"]] 
        self.next = spt[header["next"]] if spt[header["exon_type"]] == "Donnor" else spt[header["pos"]]
        # need to add the files names. for the both counts. 
        self.count = {"control": {"Donnor": {}, "Acceptor": {}}, 
                      "treatment": {"Donnor": {}, "Acceptor": {}}}
    

    def dump(self):
        l = [self.contig, self.gene,
            self.transcript, self.strand, "intron:{}".format(self.intron_number),
            self.pos, self.next]
        return l
        #return "{}\t{}\t{}\t{}\tintron:{}\tdonnor:{}\tacceptor:{}".format(self.contig, self.gene,
        #                    self.transcript, self.strand, self.intron_number,
        #                    self.pos, self.next)

    def __hash__(self):
        hash((self.contig, self.gene, self.transcript, self.intron_number, self.strand)) # , self.pos, self.next,

    def get_hash_key(self):
        return (self.contig, self.gene, self.transcript, self.intron_number, self.strand)

    def update_count(self, spt, genotype, header, sample):
        #(principal, secondary) = 
        self.count[genotype][spt[header["exon_type"]]][sample] = list(map(int, spt[header["spliced"]: ]))
        if spt[header["ambiguous"]] == "true":
            self.ambigious = True 

    
    def get_junction_count(self, counter):
        dico_count = {"control": [], "treatment": []}
        for genotype, geno_dict in self.count.items():

            for samples in geno_dict["Acceptor"]:
                try:
                    acceptor_count = geno_dict["Acceptor"][samples]
                    donnor_count = geno_dict["Donnor"][samples]
                except:
                    print(geno_dict)
                    print(self.dump(), self.ambigious)
                    raise
                (sucesses, failures) = counter.get_count(acceptor_count, donnor_count)
                dico_count[genotype].append((sucesses, failures))
            
            # dico_count[] = #
            # for exon_type, exon_dict in geno_dict.items():
            #    for file, counts in exon_dict.items():
    
        return dico_count
"""



def parse_js_file(file, results, genotype, ambigious=False):
    basename = Path(file).stem
    with open(file) as fi:
        header= fi.readline()
        header = dict((k, i) for i, k in enumerate(header.strip().split()))

        for l in fi:
            l = l.strip()
            if not l:
                continue
            spt = l.split('\t')
            if not ambigious and spt[header["Ambiguous"]] == "true":
                continue

                        
            contig = spt[header["Contig"]]
            strand = spt[header["Strand"]]
            pos =  spt[header["Donnor"]] if spt[header["Strand"]] == "+" else spt[header["Acceptor"]] 
            next = spt[header["Acceptor"]] if spt[header["Strand"]] == "+" else spt[header["Donnor"]]
            
            hash_key = (contig, strand, pos, next)

            if hash_key not in results:
                results[hash_key] = Junction(spt=spt, header=header, genotype=genotype, basename=basename)
            else:
                results[hash_key].update_count(spt=spt, genotype=genotype, header=header, basename=basename)
    pass


def main(condition_1, condition_2, successes, failures, tester, out_file,
          condition_1_Name="control", condition_2_Name="treatment", 
          ambigious=False):
    """
    This program test if a set of category(ies) (successes) has different proportion in the population
    successes + failures. Fischer Test, glm binomial, (chi2 not implemented)
    - Condition_1 and 2 are list of OmniSplice table files(s).

    In the results condition 1 is named Control and condition 2 is named treatment.
    # TODO add header with this info.
    - The comparision is condition_1 vs condition_2.

    condition_1: [Path] | [str]
    condition_2: [Path] | [str]
    principal: [str]
    secondary: [str]
    """

    counter = Counter(successes, failures)
    results = {}
    for file in condition_1:
        parse_js_file(file, results, genotype=condition_1_Name, ambigious=ambigious)
    
    for file in condition_2:
        parse_js_file(file, results, genotype=condition_2_Name, ambigious=ambigious)

    # now filter and test!
    
    dico_r = {}
    for e, v in sorted(results.items(), key=lambda x: (x[0][2], x[0][3])):
        if not v.ambigious:
            try:
                counts = v.get_junction_count(counter)
            except:
                print("failed")
                print(v.__dict__)
                raise
                continue
            if counts == (-1, -1):
                print(v.__dict__)
                continue
            #else:
            #
            #     print(counts)

            (p_value, data) = tester.test_sample(counts["control"], counts['treatment'])
            dump_list = v.dump()
            key = "{}_{}_{}_{}".format(dump_list[0], dump_list[1], dump_list[3], dump_list[4])
            dico_r[key] = [dump_list, p_value, data]
    
    # I can avoid this loop by putting it in the previsou one
    keys = []
    keys_na = []
    pvalues = []
    for k, v in dico_r.items():
        if v[1] == "NA":
            keys_na.append(k)
            continue
        keys.append(k)
        pvalues.append(v[1])


    q_values = fdrcorrection(pvalues)
    header = ["chr", "strand", "gene_transcrit_intron", "start", "end", "statistic",
               "control_success_failures_prop", "treatment_control_success_failures", "p_value", "q_value\n"]
    
    with open(out_file, "w") as fo:
        fo.write("#successes: {}\n".format(','.join(successes)))
        fo.write("#failures: {}\n".format(','.join(failures)))
        fo.write("#test: {}\n".format(tester.repr()))
        fo.write("\t".join(header))
        for i, k in enumerate(keys):
            v = dico_r[k]
            to_print = "\t".join(["\t".join(v[0]), v[2], str(v[1]), str(q_values[1][i])])
            fo.write(to_print + '\n')
        for i, k in enumerate(keys_na):
            v = dico_r[k]
            to_print = "\t".join(["\t".join(v[0]), v[2], str(v[1]), "NA"])
            fo.write(to_print + '\n')


           

#         l = [self.contig, self.gene,
#            self.transcript, self.strand, "intron:{}".format(self.intron_number),
#           self.pos, self.next]

# filter => start simple at least X principal in one condition.
# the 

# Running mode separate donnor accptor VS grouped
# Filter:
# Junction detected with at least X read per sample.
# Junction detected in at least X sample.
# Junction dected in total with X reads in a condition

# TODO Filter impl 
class Filter():
    def __init__(self):
    
        self.min_read_principal_conditions = 8

    def filter(self,junction):
        pass

"""
c1=["/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/February_2025/with_junction_file/Control_S1_L001.table",
"/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/February_2025/with_junction_file/Control_S2_L001.table",
"/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/February_2025/with_junction_file/Control_S3_L001.table",
"/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/February_2025/with_junction_file/Control_S1_L002.table",
"/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/February_2025/with_junction_file/Control_S2_L002.table",
"/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/February_2025/with_junction_file/Control_S3_L002.table"]

c2=["/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/February_2025/with_junction_file/U2af38_S1_L001.table",
"/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/February_2025/with_junction_file/U2af38_S2_L001.table",
"/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/February_2025/with_junction_file/U2af38_S3_L001.table",
"/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/February_2025/with_junction_file/U2af38_S1_L002.table",
"/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/February_2025/with_junction_file/U2af38_S2_L002.table",
"/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/February_2025/with_junction_file/U2af38_S3_L002.table"]

#main(c1, c2, ['spliced'], ['unspliced', 'exon_other', 'exon_intron', 'clipped'])
"""

if __name__ == "__main__":
    parse = argparse.ArgumentParser(description="""
    this program will compare two set of ".junction" files and compare the proportion of 2 groups of splicing groups.
    Allowing to identify splicing junction perturbed in one condition vs. the others.
    
    3 statistical test are available Fischer, [Chi2 not implement yet TODO] and GLM. 
    Fischer test and Chi2 pool samples together. Chi2 limitations will be enforced junction with low count will output NA.
    GLM is based on a binomial distribution and use this formula: 'successes + failures ~ group'
    GLM will take into account the varriance between sample from the same group.
                                    
    usage: > python3 compare_conditions.py -c  control_1.junction control_2.junction -t treatment_1.junction treatment_2.junction --spliced SPLICED 
            --defect UNSPLICED EXON_INTRON EXON_OTHER --stat FISCHER --out my_comparison_file.tsv
                                
    """)
    parse.add_argument("--control", "-c", required=True,
                        help="space separated list of control file. eg: --control control_1.junction control_2.junction",
                        nargs="+")
    parse.add_argument("--treatment", "-t", required=True,
                        help="space separated list of treatment file. eg: --treatment treatment_1.junction treatment_2.junction",
                        nargs="+")
    
    parse.add_argument("--spliced", help="what category to count as 'normal' default spliced", default="spliced", nargs="+", choices=["SPLICED", 
                                           "UNSPLICED", 
                                           "CLIPPED",  
                                           "SKIPPED",
                                           "EXON_OTHER",
                                           "WRONG_STRAND",
                                           "E_ISOFORM"])
    parse.add_argument("--defect", help="", required=True, nargs="+", choices=["SPLICED", 
                                           "UNSPLICED", 
                                           "CLIPPED", 
                                           "EXON_OTHER",
                                           "SKIPPED",
                                           "WRONG_STRAND",
                                           "E_ISOFORM"])
    
    parse.add_argument("--stat", required=True, choices=["GLM", "FISCHER", "CHI2"])
    parse.add_argument("--out", required=True)
    parse.add_argument("--control_name", default="control", help="replace control condition by the name of your choice")
    parse.add_argument("--treatment_name", default="treatment", help="replace treatment condition by the name of your choice")
    parse.add_argument("--ambigious", action='store_true', help="flag, do you want to take into account abigious junction, default is False. to set to true add '--amigious' to the command")
    args = parse.parse_args()
    #parse.add_argument()


    for e in args.control:
        try:
            assert Path(e).is_file()
        except AssertionError:
            print("file {e} not found please check you inputed the correct path".format(e))
            raise
    for e in args.treatment:
        try:
            assert Path(e).is_file()
        except AssertionError:
            print("file {e} not found please check you inputed the correct path".format(e))
            raise
    if type(args.spliced) == str:
        args.spliced = [args.spliced]
    if type(args.defect) == str:
        args.defect = [args.defect]

    cond1 = list(set([x.lower() for x in args.spliced]))
    cond2 = list(set([x.lower() for x in args.defect]))

    tester = None
    if args.stat == "GLM":
        try:
            assert(len(args.control) > 3 and len(args.treatment) > 3)
        except AssertionError:
            print("for glm you need at least tree replicates per conditions")
            raise AssertionError
        except:
            raise
        tester = GLM_model()
    elif args.stat == "FISCHER":
        tester = Fischer_test()
    elif args.stat == "CHI2":
        print("CHI2 test not yet avalaible")
        raise AssertionError
    if not tester:
        print("cannot find {} please checks your input".format(args.stat))
        raise AssertionError

    main(condition_1=args.control,condition_2=args.treatment, successes=cond1, failures=cond2,
          tester=tester, condition_1_Name=args.control_name, condition_2_Name=args.treatment_name,
           ambigious=args.ambigious, out_file=args.out)
    

