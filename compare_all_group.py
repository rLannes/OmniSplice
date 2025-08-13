import numpy as np
from scipy.stats import dirichlet_multinomial
from scipy.optimize import minimize
from scipy.stats import chi2
from functools import partial


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