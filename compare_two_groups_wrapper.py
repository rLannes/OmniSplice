import argparse
from pathlib import Path
import numpy as np
from  statsmodels.stats.multitest import fdrcorrection
import logging
import os
import time
from common_python.counter_junction import Counter
from common_python.junction_class import Junction, parse_js_file
from common_python.tester import Chi2, Fischer_test, GLM_model
import subprocess
from pathlib import Path


logging.basicConfig(format='%(asctime)s %(levelname)s => %(message)s', level=logging.INFO, datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("compare_conditions")



if __name__ == "__main__":
    parse = argparse.ArgumentParser(description="""
    this program will compare two set of ".junction" files/
    
                                    
    usage: > python3 compare_two_groups.py -c  control_1.junction control_2.junction -t treatment_1.junction treatment_2.junction --out my_comparison_file.tsv
                                
    """)
    parse.add_argument("--control", "-c", required=True,
                        help="space separated list of control file. eg: --control control_1.junction control_2.junction",
                        nargs="+")
    parse.add_argument("--treatment", "-t", required=True,
                        help="space separated list of treatment file. eg: --treatment treatment_1.junction treatment_2.junction",
                        nargs="+")
    
    #parse.add_argument("--stat", required=True, choices=["GLM", "FISCHER", "CHI2"])
    parse.add_argument("--signif_level", "-s", help="if set filter junction with q_vlaue higher than this option" )
    parse.add_argument("--out", required=True, help="basename")
    #parse.add_argument("--control_name", default="control", help="replace control condition by the name of your choice")
    #parse.add_argument("--treatment_name", default="treatment", help="replace treatment condition by the name of your choice")
    #parse.add_argument("--ambigious", action='store_true', help="flag, do you want to take into account abigious junction, default is False. to set to true add '--amigious' to the command")
    parse.add_argument("--logging_level", "-v", default="ERROR",choices=["DEBUG", "INFO", "ERROR"] )
    args = parse.parse_args()


    if args.logging_level == "DEBUG":
        logger.setLevel(logging.DEBUG)
    elif args.logging_level == "INFO": 
        logger.setLevel(logging.INFO)
    elif args.logging_level == "ERROR":
        logger.setLevel(logging.ERROR)

    script_path = Path(os.path.realpath(__file__)).absolute()
    script_dir = script_path.parent
    compare_ = script_dir / "compare_two_groups.py"
    ctr = " ".join(args.control)
    treat = " ".join(args.treatment)
    for defect in ["CLIPPED", "EXON_OTHER", "SKIPPED", "WRONG_STRAND", "E_ISOFORM"]:
        child = subprocess.Popen(f"python3 {compare_} -c {ctr} --spliced SPLICED \
                                 --defect {defect} -t {treat} --out {args.out}_{defect}.tsv --ambigious --logging_level {args.logging_level}", shell=True)
        child.wait()

    child = subprocess.Popen(f"python3 {compare_} -c {ctr} --spliced SPLICED \
                             --defect UNSPLICED -t {treat} --out {args.out}_UNSPLICED.tsv --logging_level {args.logging_level}", shell=True)
    child.wait()


