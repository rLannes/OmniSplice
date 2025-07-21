from math import isclose
from scipy.stats import fisher_exact
from fast_fisher import fast_fisher_exact_compatibility
from scipy.stats import chisquare
import pandas as pd
import statsmodels
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import chi2
from abc import abstractmethod
from enum import Enum
from abc import ABC
import logging
import numpy as np

import warnings
warnings.filterwarnings('ignore')

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
        (odds_f, pval_f ) = fast_fisher_exact_compatibility(data, alternative='two-sided')
        return (pval_f, "{}".format(odds_f) + "\t" + self.get_prop(data[0], "control") + "\t" + self.get_prop(data[1], "treatment"))

   
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
        
        return "{}/{}/{}".format(successes, failures, p)

class GLM_model(Stat_test):

    def __init__(self):
        super().__init__()
    
    def repr(self):
        return "GLM Binomial: (successes + failures) ~ group"
    
    
    def test_sample(self, control, treatment):
        
        data = self.load_data(control, treatment)
        if sum(data["successes"].values) == 0:
            return (np.nan, "NA" + '\t' + self.get_prop(data, "control") + '\t' + self.get_prop(data, "treatment"))
        
        data = data[~((data["successes"] == 0) & (data["failures"] == 0))]
        
        if len(list(set(data["group"].values))) < 2 or data.empty:
            return (np.nan, "NA" + '\t' + self.get_prop(data, "control") + '\t' + self.get_prop(data, "treatment"))

        try:
            mod = smf.glm('successes + failures ~ group', family=sm.families.Binomial(), data=data).fit()
            odr = float(np.exp(mod.params["group"]))
            
            return (float(mod.pvalues[1]),  "{}".format(round(odr, 4)) + '\t' + self.get_prop(data, "control")  + '\t' + self.get_prop(data, "treatment"))
        
        # could use Fischer in case of perfect or quasi separation
        except statsmodels.tools.sm_exceptions.PerfectSeparationError:
            return (np.nan, "PS" + '\t' + self.get_prop(data, "control") + '\t' + self.get_prop(data, "treatment"))
        except KeyError:

            #print(data)
            #print(mod.params)
            odr = float(np.exp(mod.params.values[1]))
            
            #print(data, mod.summary(), float(mod.pvalues[1]),  odr)
            
            return (float(mod.pvalues[1]),  "{}".format(round(odr, 4)) + '\t' + self.get_prop(data, "control")  + '\t' + self.get_prop(data, "treatment"))

            #return (np.nan, "QPS" + '\t' + self.get_prop(data, "control") + '\t' + self.get_prop(data, "treatment"))
        except:
            print(data)
            print(mod.params)
            raise
            return (np.nan, "NA" + '\t' + self.get_prop(data, "control") + '\t' + self.get_prop(data, "treatment"))


    def load_data(self, control, treatment):

        df =  pd.DataFrame({"successes": [0,0,0]})
        try:
            df = pd.DataFrame( {
                               "successes": [x for x in control["successes"]] + [x for x in treatment["successes"]],
                               "failures": [x for x in control["failures"]] + [x for x in treatment["failures"]],
                               "group": ["control"]*len(control["successes"]) + ["treatment"]*len(treatment["successes"])
                               } )
        except:
            print( control, treatment)
        return df
    

    def get_prop(self, df, group):

        successes, failures = 0, 0 
        if successes in df["group"]:
            successes = sum(df[df["group"] == group]["successes"])  
        if failures in df["group"]:
            failures =  sum(df[df["group"] == group]["failures"])
        p = -1

        if successes + failures > 0:
            p = successes / (successes + failures)        
        
        return "{}/{}/{}".format(successes, failures, p)



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