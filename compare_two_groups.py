import argparse
from pathlib import Path
import numpy as np
from  statsmodels.stats.multitest import fdrcorrection
import logging
import time
from common_python.counter_junction import Counter
from common_python.junction_class import Junction, parse_js_file
from common_python.tester import FisherTest, GLMModel


logging.basicConfig(format='%(asctime)s %(levelname)s => %(message)s', level=logging.INFO, datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("compare_conditions")


def timer_decorator(func):
    def inner1(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)

        print('Function', func.__name__, 'time:', time.time() - start)

        return result
    return inner1



def main(condition_1, condition_2, successes, failures, tester, out_file,
          ambiguous=False):
    """
    This program test if a set of category(ies) (successes) has different proportion in the population
    successes + failures. Fisher Test, glm binomial, (chi2 not implemented)
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
        logging.info("parsing control file {}; ".format(file))
        parse_js_file(file, results, genotype="control", ambiguous=ambiguous)
    
    for file in condition_2:
        logging.info("parsing treatment file {}; ".format(file))
        parse_js_file(file, results, genotype="treatment", ambiguous=ambiguous)

    # now filter and test!
    
    dico_r = {}
    for e, v in sorted(results.items(), key=lambda x: (x[0][2], x[0][3])):
        if ambiguous or not v.ambiguous:
            try:
                v.stat_test(counter, tester)
            except:
                print("failed")
                print(e, v.__dict__)
                raise

    junctions = list(results.values())
    junctions = [j for j in junctions if not np.isnan(j.p_value) ]
    #y  = fdrcorrection([x.p_value for x in junctions])
    
    q_values = fdrcorrection([x.p_value for x in junctions])[1] #[ y[0] for y in fdrcorrection([x.p_value for x in junctions]) ]
    #print(q_values)

    for i, v in enumerate(q_values):
        junctions[i].q_value = v

    if args.signif_level:
        signif_level = float(args.signif_level)
        junctions = [x for x in junctions if x.q_value < signif_level]

    junctions = sorted(junctions, key = lambda x: x.q_value if x.q_value else 2)
        
    header = ["chr", "strand", "start", "end", "statistic", "p_value", "q_value", "status", 
               "control_success", "control_failures", "control_ratio", "treatment_success", "treatment_failures", "treatment_ratio", "gene_transcript_intron\n"]
    
    with open(out_file, "w") as fo:
        fo.write("#successes: {}\n".format(','.join(successes)))
        fo.write("#failures: {}\n".format(','.join(failures)))
        fo.write("#test: {}\n".format(tester.repr()))
        fo.write("\t".join(header))
        for junction in junctions:
            #print([junction.dump(), junction.data_stats, str(junction.p_value), str(junction.q_value) ])
            d = junction.dump()
            try:
                fo.write("\t".join([d[0], d[1], d[3], d[4], str(junction.testResult.statistics), \
                                    str(junction.p_value), str(junction.q_value), junction.testResult.status.value, \
                                    junction.testResult.control_prop.dump(), junction.testResult.treatment_prop.dump(),  \
                                    d[2] ]) + "\n" )
            except:
                print("failed", junction.__dict__)
                continue

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


if __name__ == "__main__":
    parse = argparse.ArgumentParser(description="""
    this program will compare two set of ".junction" files and compare the proportion of 2 groups of splicing groups.
    Allowing to identify splicing junction perturbed in one condition vs. the others.
    
                                    
    usage: > python3 compare_conditions.py -c  control_1.junction control_2.junction -t treatment_1.junction treatment_2.junction --spliced SPLICED 
            --defect UNSPLICED EXON_INTRON EXON_OTHER --out my_comparison_file.tsv
                                
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
                                           "SKIPPEDUNRELATED",
                                           "EXON_OTHER",
                                           "WRONG_STRAND",
                                           "E_ISOFORM"])
    parse.add_argument("--defect", help="", required=True, nargs="+", choices=["SPLICED", 
                                           "UNSPLICED", 
                                           "CLIPPED", 
                                           "EXON_OTHER",
                                           "SKIPPED",
                                            "SKIPPEDUNRELATED",
                                           "WRONG_STRAND",
                                           "E_ISOFORM"])
    
    parse.add_argument("--stat", required=True, choices=["GLM", "FISHER"])
    parse.add_argument("--signif_level", "-s", help="if set filter junction with q_vlaue higher than this option" )
    parse.add_argument("--out", required=True)
    #parse.add_argument("--control_name", default="control", help="replace control condition by the name of your choice")
    #parse.add_argument("--treatment_name", default="treatment", help="replace treatment condition by the name of your choice")
    parse.add_argument("--ambiguous", action='store_true', help="flag, do you want to take into account abigious junction, default is False. to set to true add '--amigious' to the command")
    parse.add_argument("--logging_level", "-v", default="ERROR",choices=["DEBUG", "INFO", "ERROR"] )
    args = parse.parse_args()


    if args.logging_level == "DEBUG":
        logger.setLevel(logging.DEBUG)
    elif args.logging_level == "INFO": 
        logger.setLevel(logging.INFO)
    elif args.logging_level == "ERROR":
        logger.setLevel(logging.ERROR)


    for e in args.control:
        try:
            assert Path(e).is_file()
        except AssertionError:
            logger.error("file {e} not found please check you inputed the correct path".format(e))
            raise
    for e in args.treatment:
        try:
            assert Path(e).is_file()
        except AssertionError:
            logger.error("file {e} not found please check you inputed the correct path".format(e))
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
            assert(len(args.control) >= 3 and len(args.treatment) >= 3)
        except AssertionError:
            print("for glm you need at least tree replicates per conditions")
            raise AssertionError
        except:
            raise
        tester = GLMModel()

    elif args.stat == "FISHER":
        tester = FisherTest()

    #elif args.stat == "CHI2":
    #    print("CHI2 test not yet avalaible")
    #    raise AssertionError
    
    if not tester:
        print("cannot find {} please checks your input".format(args.stat))
        raise AssertionError

    #tester = FisherTest()
    
    logger.debug("control file: {}".format(args.control))
    logger.debug("treatment file: {}".format(args.treatment))
    main(condition_1=args.control,condition_2=args.treatment, successes=cond1, failures=cond2,
          tester=tester,
           ambiguous=args.ambiguous, out_file=args.out)
    

