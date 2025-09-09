import argparse
from multiprocessing import Process, Pool
from common_python.junction_class import Junction, parse_js_file
from  statsmodels.stats.multitest import fdrcorrection
from functools import partial
import logging
from pathlib import Path
import numpy as np
try:
    from scipy.stats import dirichlet_multinomial
except:
    print("you should update to a scipy version that support dirichlet_multinomial")
    raise
from scipy.optimize import minimize
from scipy.stats import chi2



logging.basicConfig(format='%(asctime)s %(levelname)s => %(message)s', level=logging.INFO, datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger()


# Null model: one alpha for all replicates
def nll_null(params, groups):
    alpha = np.exp(params)  # ensure positive
    ll = 0
    for g in groups:
        for counts in g:
            ll += dirichlet_multinomial.logpmf(counts, n=counts.sum(), alpha=alpha)
    return -ll

# Full model: one alpha per condition
def nll_full(params, groups):
    k = int( len(params) / 2 )
    alpha0 = np.exp(params[:k])
    alpha1 = np.exp(params[k:])
    ll = 0
    for counts in groups[0]:
        ll += dirichlet_multinomial.logpmf(counts, n=counts.sum(), alpha=alpha0)
    for counts in groups[1]:
        ll += dirichlet_multinomial.logpmf(counts, n=counts.sum(), alpha=alpha1)
    return -ll



def test_dirichlet_multinomial(g1, g2):


    Y0 = np.array(g1)
    Y1 = np.array(g2)

    # Combine into list for easy looping
    groups = [Y0, Y1]

    # Fit models
    k = Y0.shape[1]

    nll_null_bound = partial(nll_null, groups=groups)
    nll_full_bound = partial(nll_full, groups=groups)

    res_null = minimize(nll_null_bound, np.zeros(k), method='BFGS')
    res_full = minimize(nll_full_bound, np.zeros(2 * k), method='BFGS')

    ll_null = -res_null.fun
    ll_full = -res_full.fun

    # Likelihood ratio test
    lr_stat = 2 * (ll_full - ll_null)
    df = k  # extra parameters in full model
    p_value = chi2.sf(lr_stat, df)

    return(p_value, ll_null, ll_full, lr_stat, df)
    #print(f"Null log-likelihood: {ll_null:.3f}")
    #print(f"Full log-likelihood: {ll_full:.3f}")
    #print(f"LR stat: {lr_stat:.3f}, df={df}, p={p_value:.4g}")


def get_chuncksize(thread, n_junction):
    if n_junction %  thread == 0:
        return n_junction / thread
    else:
        return int(n_junction / thread) + 1


def run_mp(args):
    junctions_pairs = args

    for junction_id, junction in junctions_pairs:
        count_treatment = []
        count_control = []

        for sample_name, counts in junction.count['control'].items():
            count_control.append(counts)

        for sample_name, counts in junction.count["treatment"].items():
            count_treatment.append(counts)
        try:
            (p_value, ll_null, ll_full, lr_stat, df) = test_dirichlet_multinomial(count_control, count_treatment)
        except ValueError:
            print("junction {} is a vector with null value, skipping it.".format(junction_id))
            continue  
        junction.pvalue = p_value
        junction.arbitrary["ll_null"] = ll_null
        junction.arbitrary["ll_full"] = ll_full
        junction.arbitrary["lr_stat"] = lr_stat
        junction.arbitrary["df"] = df



def main(condition_1, condition_2, ambigious, out_file, thread, min_read, min_spliced):
    results = {}
    for file in condition_1:
        logging.info("parsing control file {}; ".format(file))
        parse_js_file(file, results, genotype="control", ambigious=ambigious)
    
    for file in condition_2:
        logging.info("parsing treatment file {}; ".format(file))
        parse_js_file(file, results, genotype="treatment", ambigious=ambigious)

    to_delete = []
    if min_read > 0 or min_spliced > 0:
        print('filtering out junction based on read counts')
        for junction_id, junction in results.items():
            if min_read > 0 and not junction.pass_min_read_filter(min_read):
                to_delete.append(junction_id)
            if min_spliced > 0 and not junction.pass_min_spliced_filter(min_spliced):
                to_delete.append(junction_id)
        print("filtering done")
    
    to_delete = list(set(to_delete))
    for key in to_delete:
        del results[key]



    # rewrite for multiprocess
    if thread > 1 :
        junctions_pairs = list(results.items())
        chunk = get_chuncksize(n_junction=len(junctions_pairs), thread=thread)
        junctions_pairs = [junctions_pairs[i:i+chunk] for i in range(0, len(junctions_pairs), chunk)]
        print("starting pool")
        with Pool(processes=thread) as pool:
            child = pool.map(run_mp, junctions_pairs)#, chunksize=100)

            pool.join()
        print('pool done')
    
    else:
        for junction_id, junction in results.items():
            count_treatment = []
            count_control = []

            for sample_name, counts in junction.count['control'].items():
                count_control.append(counts)

            for sample_name, counts in junction.count["treatment"].items():
                count_treatment.append(counts)
            try:
                (p_value, ll_null, ll_full, lr_stat, df) = test_dirichlet_multinomial(count_control, count_treatment)
            except ValueError:
                print("junction {} at least one replicate is a vector with null value, skipping it.".format(junction_id))
                continue  
            #print(p_value, ll_null, ll_full, lr_stat, df, junction.__dict__, count_control, count_treatment)
            junction.pvalue = p_value
            junction.arbitrary["ll_null"] = ll_null
            junction.arbitrary["ll_full"] = ll_full
            junction.arbitrary["lr_stat"] = lr_stat
            junction.arbitrary["df"] = df

    junctions = list(results.values())
    junctions = [j for j in junctions if not np.isnan(j.p_value) ]
    
    q_values = fdrcorrection([x.p_value for x in junctions])[1]

    for i, v in enumerate(q_values):
        junctions[i].q_value = v

    if args.signif_level:
        signif_level = float(args.signif_level)
        junctions = [x for x in junctions if x.q_value < signif_level]

    junctions = sorted(junctions, key = lambda x: x.q_value if x.q_value else 2)
    

    header = ["chr", "strand", "start", "end", "statistic",
               "control_success_failures_prop", "treatment_control_success_failures", "p_value", "q_value", "gene_transcript_intron\n"]
    
    with open(out_file, "w") as fo:
        fo.write("#Dirichlet multonomial test")
        fo.write("\t".join(header))
        for junction in junctions:
            d = junction.dump()
            try:
                fo.write("\t".join([d[0], d[1], d[3], d[4],  str(junction.p_value), str(junction.q_value),
                                    junction.arbitrary["ll_null"], junction.arbitrary["ll_full"],
                                    junction.arbitrary["lr_stat"], junction.arbitrary["df"], d[2] ]) + "\n" )
            except:
                print("failed", junction.__dict__)
                continue



if __name__ == "__main__":

    parse = argparse.ArgumentParser()

    parse.add_argument("--control", "-c", required=True,
                        help="space separated list of control file. eg: --control control_1.junction control_2.junction",
                        nargs="+")
    parse.add_argument("--treatment", "-t", required=True,
                        help="space separated list of treatment file. eg: --treatment treatment_1.junction treatment_2.junction",
                        nargs="+")
    parse.add_argument("--signif_level", "-s", help="if set filter junction with q_vlaue higher than this option" )
    parse.add_argument("--out", required=True)
    parse.add_argument("--thread", default=8, type=int)
    #parse.add_argument("--control_name", default="control", help="replace control condition by the name of your choice")
    #parse.add_argument("--treatment_name", default="treatment", help="replace treatment condition by the name of your choice")
    parse.add_argument("--logging_level", "-v", default="ERROR",choices=["DEBUG", "INFO", "ERROR"] )

    group = parse.add_argument_group('junction filtering')
    group.add_argument("--min_spliced", type=int, default=-1, help="default -1 deactivated. set to any positive integer (>0)to filter out junction. Where at least one replicate from both condition has spliced read count under this value")
    group.add_argument("--min_reads", type=int, default=-1, help="default -1 deactivated. set to any positive integer (>0) to filter out junction where any replicate has total read count under this value")
    group.add_argument("--ambigious", action='store_true', help="flag, do you want to take into account abigious junction, default is False. to set to true add '--amigious' to the command")

    
    args = parse.parse_args()

    # argument validation

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


    logger.debug("control file: {}".format(args.control))
    logger.debug("treatment file: {}".format(args.treatment))
    main(condition_1=args.control,condition_2=args.treatment,
           ambigious=args.ambigious, out_file=args.out, thread=args.thread, min_read=args.min_reads, min_spliced=args.min_spliced)
    