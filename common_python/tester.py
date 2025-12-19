"""
Statistical testing module for comparing binary outcomes between control and treatment groups.

This module provides implementations of Fisher's exact test and generalized linear models (GLM)
for testing differences in proportions between two groups. It handles edge cases such as empty
datasets, zero counts, and perfect separation, returning standardized results with p-values,
statistics, and proportion summaries.

Classes:
    StatTest: Abstract base class defining the interface for statistical tests.
    TestResult: Dataclass containing test results (p-value, statistics, proportions).
    FisherTest: Implementation of Fisher's exact test for 2x2 contingency tables.
    GLMModel: Implementation using binomial GLM to test group differences.
    TestTester: Unit test suite validating test implementations.

Example:
    >>> control = {"successes": [10, 12, 8], "failures": [5, 3, 7]}
    >>> treatment = {"successes": [15, 18, 12], "failures": [2, 1, 3]}
    >>> tester = FisherTest()
    >>> result = tester.test_sample(control, treatment)
    >>> print(result.pvalue)
"""

from enum import Enum
from scipy.stats import fisher_exact
from dataclasses import dataclass
import pandas as pd
import statsmodels
import statsmodels.api as sm
import statsmodels.formula.api as smf
from abc import abstractmethod
from abc import ABC
import logging
import numpy as np
import unittest
import warnings

class TestStatus(Enum):
    OK = "Ok"
    ZEROES = "Zeroes"
    EMPTY = "Empty"
    MISMATCHED_LENGTH = "MismatchedLength"
    PERFECT_SEPARATION = "PerfectSeparation"
    POSSIBLE_PERFECT_SEPARATION = "possiblePerfectSeparation"
    CONVERGENCE_FAILED = "ConvergenceFailed"

INVALID_COUNT_THRESHOLD = 0
MIN_OBSERVATIONS_PER_GROUP = 1

@dataclass
class ProportionResult:
    successes: int
    failures: int
    proportion: float

    @staticmethod
    def failed():
        return ProportionResult(0,0,-1)

    def as_list(self):
        return [self.successes, self.failures, self.proportion]

    def dump(self):
        return "\t".join(map(str, self.as_list()))


@dataclass
class TestResult:
    """
    Container for statistical test results.
    
    Attributes:
        pvalue (float): P-value from the statistical test.
        statistics (float): Test statistic (e.g., odds ratio).
        control_prop (str): Tab-separated string of control group statistics.
        treatment_prop (str): Tab-separated string of treatment group statistics.
        any_data (dict): Dictionary for storing additional test-specific data.
    """
    pvalue: float
    statistics:  float
    control_prop: ProportionResult
    treatment_prop: ProportionResult
    status: TestStatus
    any_data: dict

    def dump(self):
        return "{}\t{}\t{}".format( str(self.statistics), self.control_prop.dump(), self.treatment_prop.dump())


logger = logging.getLogger(__name__)
fisher = fisher_exact
try:
    from fast_fisher import fast_fisher_exact_compatibility
    fisher = fast_fisher_exact_compatibility
except ImportError:
    logger.info("fast_Fisher is not installed, using scipy Fisher test instead (much slower)")

"""
Abstract base class for statistical tests comparing control and treatment groups.

Defines the interface that all statistical test implementations must follow,
including methods for testing samples, loading data, computing proportions,
and quality checking inputs.
"""
class StatTest(ABC):

    @abstractmethod
    def test_sample(self, control, treatment):
        """
        Perform statistical test comparing control and treatment groups.
        
        Args:
            control (dict): Dictionary with keys "successes" and "failures", each containing
                           lists of counts for the control group.
            treatment (dict): Dictionary with keys "successes" and "failures", each containing
                             lists of counts for the treatment group.
        
        Returns:
            TestResult: Object containing p-value, test statistic, and proportion strings.
        """
        pass

    @abstractmethod
    def load_data(self, control, treatment):
        """
        Load and format data from control and treatment groups for testing.
        
        Args:
            control (dict): Control group data with "successes" and "failures" keys.
            treatment (dict): Treatment group data with "successes" and "failures" keys.
        
        Returns:
            Data structure appropriate for the specific test (list or DataFrame).
        """
        pass

    @abstractmethod
    def get_prop(self, val, group):
        """
        Calculate and format proportion statistics for a group.
        
        Args:
            val: Data values for the group (format depends on test implementation).
            group (str): Group identifier ("control" or "treatment").
        
        Returns:
            str: Tab-separated string with format "successes\tfailures\tproportion".
        """
        pass

    @abstractmethod
    def repr(self):
        """
        Return human-readable name of the statistical test.
        
        Returns:
            str: Name/description of the test method.
        """
        pass

    def quality_check(self, control, treatment):
        """
        Validate input data quality and structure before running tests.
        
        Checks for mismatched lengths between successes and failures lists,
        and for empty datasets, and missing keys.
        
        Args:
            control (dict): Control group data.
            treatment (dict): Treatment group data.
        
        Returns:
            TestResult or None: Returns TestResult with NaN values if validation fails,
                               None if validation passes.
        """
        if "successes" not in control or "successes" not in treatment or "failures" not in control or "failures" not in treatment:
            return TestResult(np.nan, np.nan, ProportionResult.failed(), ProportionResult.failed(), TestStatus.EMPTY, {})
        
        if len(control["successes"]) != len(control["failures"]) or len(treatment["successes"]) != len(treatment["failures"]):
            return TestResult(np.nan, np.nan, ProportionResult.failed(), ProportionResult.failed(), TestStatus.MISMATCHED_LENGTH, {})

        if len(control["successes"]) == 0 or len(treatment["successes"]) == 0:
            return TestResult(np.nan, np.nan, ProportionResult.failed(), ProportionResult.failed(), TestStatus.EMPTY, {})
        
        if (sum(control["successes"]) == 0 and sum(control["failures"]) == 0) or (sum(treatment["successes"]) == 0 and sum(treatment["failures"]) == 0):
            return TestResult(np.nan, np.nan, ProportionResult.failed(), ProportionResult.failed(), TestStatus.ZEROES, {})

class FisherTest(StatTest):
    """
    Implementation of Fisher's exact test for 2x2 contingency tables.
    
    Tests for association between binary outcomes and group membership using
    Fisher's exact test. Sums all successes and failures across replicates
    to create a single 2x2 table for analysis.
    """
    def __init__(self): 
        """Initialize FisherTest instance."""
        super().__init__()
    
    def repr(self):
        """
        Return test identifier.
        
        Returns:
            str: "FISHER test"
        """
        return "FISHER test"
    
    def test_sample(self, control, treatment):
        """
        Perform Fisher's exact test on control and treatment groups.
        
        Aggregates successes and failures across all replicates and performs
        a two-sided Fisher's exact test on the resulting 2x2 table.
        
        Args:
            control (dict): Control group with "successes" and "failures" lists.
            treatment (dict): Treatment group with "successes" and "failures" lists.
        
        Returns:
            TestResult: Contains p-value, odds ratio, and proportion summaries.
        """

        if  x := self.quality_check(control, treatment):
            return x
        
        data = self.load_data(control, treatment)
        if sum(data[0]) == sum(data[1]) == 0:
            return TestResult(np.nan, np.nan, ProportionResult.failed(), ProportionResult.failed(), TestStatus.ZEROES, {})

        (odds_f, pval_f ) = fisher(data, alternative='two-sided')

        return TestResult(pval_f, odds_f, self.get_prop(data[0], "control"), self.get_prop(data[1], "treatment"), TestStatus.OK, {})

   
    def load_data(self, control, treatment):
        """
        Aggregate counts into 2x2 contingency table.
        
        Sums positive counts from successes and failures lists for both groups,
        filtering out zero or negative values.

        Filtere negative and zeroes values because < 0 indicate a potential issue in the count
        
        Args:
            control (dict): Control group data.
            treatment (dict): Treatment group data.
        
        Returns:
            list: [[control_successes, control_failures], [treatment_successes, treatment_failures]]
        """
        cc = [sum([x for x in control["successes"] if x >= INVALID_COUNT_THRESHOLD ]), sum([x for x in control["failures"] if x >= INVALID_COUNT_THRESHOLD])]
        tt = [sum([x for x in treatment["successes"] if x >= INVALID_COUNT_THRESHOLD]), sum([x for x in treatment["failures"] if x >= INVALID_COUNT_THRESHOLD])]
        return [cc, tt]


    def get_prop(self, val, group):
        """
        Calculate proportion from aggregated counts.
        
        Args:
            val (list): [successes, failures] for the group.
            group (str): Group name (for logging/debugging).
        
        Returns:
            str: "successes\tfailures\tproportion" or "0\t0\t-1" if no data.
        """
        successes = val[0]
        failures = val[1]
        p = -1
        if successes + failures > 0:
            p = successes / (successes + failures)
        
        return ProportionResult(successes, failures, p) #"{}\t{}\t{}".format(successes, failures, p)

class GLMModel(StatTest):
    """
    Generalized linear model with binomial family for testing group differences.
    
    Fits a binomial GLM with formula "(successes + failures) ~ group" to test
    whether success rates differ between control and treatment groups while
    accounting for replicate-level variation.
    """
    def __init__(self):
        """Initialize GLMModel instance."""
        super().__init__()
    
    def repr(self):
        """
        Return test identifier.
        
        Returns:
            str: "GLM Binomial: (successes + failures) ~ group"
        """
        return "GLM Binomial: (successes + failures) ~ group"
    
    
    def test_sample(self, control, treatment):
        """
        Fit binomial GLM to test group effect on success rate.
        
        Constructs a DataFrame from control and treatment data, filters out
        rows with zero counts, and fits a binomial GLM. Handles perfect
        separation and requires at least 3 non-zero observations per group.
        
        Args:
            control (dict): Control group with "successes" and "failures" lists.
            treatment (dict): Treatment group with "successes" and "failures" lists.
        
        Returns:
            TestResult: Contains p-value for group effect, odds ratio (exp(coefficient)),
                       and proportion summaries. Returns NaN/special codes for edge cases.
        """

        if  x := self.quality_check(control, treatment):
            return x

        data = self.load_data(control, treatment)

        
        if sum(data["successes"].values) == 0:
            return TestResult(np.nan, np.nan, ProportionResult.failed(), ProportionResult.failed(), TestStatus.ZEROES , {})

        # does not affect model fit
        data = data[~((data["successes"] == 0) & (data["failures"] == 0))]
        # 
        g = data.groupby("group").count()

        if g.loc["control"]["successes"] < MIN_OBSERVATIONS_PER_GROUP or g.loc["treatment"]["successes"] < MIN_OBSERVATIONS_PER_GROUP:
            return TestResult(np.nan, np.nan, ProportionResult.failed(),ProportionResult.failed(), TestStatus.EMPTY, {})
        
        if len(list(set(data["group"].values))) < 2 or data.empty :
            return TestResult(np.nan, np.nan, ProportionResult.failed(), ProportionResult.failed(), TestStatus.EMPTY, {})
        
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', 
                                message='divide by zero encountered in scalar divide')
            warnings.filterwarnings('ignore', 
                                message='invalid value encountered in scalar divide')


            try:
                mod = smf.glm('successes + failures ~ group', family=sm.families.Binomial(), data=data).fit()
                odr = float(np.exp(mod.params["group"]))
                logger.info(mod.params)
                if not mod.converged:
                    return TestResult(np.nan, np.nan,  self.get_prop(data, "control"), self.get_prop(data, "treatment"), TestStatus.CONVERGENCE_FAILED, {})
                if np.any(np.abs(mod.params) > 20): 
                    return TestResult(np.nan, np.nan,  self.get_prop(data, "control"), self.get_prop(data, "treatment"), TestStatus.POSSIBLE_PERFECT_SEPARATION, {})
                
                return TestResult(float(mod.pvalues[1]), odr, self.get_prop(data, "control"), self.get_prop(data, "treatment"), TestStatus.OK, {})

            except statsmodels.tools.sm_exceptions.PerfectSeparationError:
                return TestResult(np.nan, np.nan,  self.get_prop(data, "control"), self.get_prop(data, "treatment"), TestStatus.PERFECT_SEPARATION, {})
            except statsmodels.tools.sm_exceptions.ConvergenceWarning:
                return TestResult(np.nan, np.nan,  self.get_prop(data, "control"), self.get_prop(data, "treatment"), TestStatus.CONVERGENCE_FAILED, {})
            except KeyError:
                odr = float(np.exp(mod.params.values[1]))
                return TestResult(float(mod.pvalues[1]), odr, self.get_prop(data, "control"), self.get_prop(data, "treatment"), TestStatus.OK, {})
                
            except Exception:
                logger.error(data)
                logger.error(mod.params)
                raise


    def load_data(self, control, treatment):        
        """
        Create DataFrame with replicate-level data for GLM fitting.
        
        Args:
            control (dict): Control group data.
            treatment (dict): Treatment group data.
        
        Returns:
            pd.DataFrame: DataFrame with columns "successes", "failures", "group".
        """

        df =  pd.DataFrame({"successes": [0,0,0]})
        try:
            df = pd.DataFrame( {
                               "successes": [x for x in control["successes"]] + [x for x in treatment["successes"]],
                               "failures": [x for x in control["failures"]] + [x for x in treatment["failures"]],
                               "group": ["control"]*len(control["successes"]) + ["treatment"]*len(treatment["successes"])
                               } )
        except:
            logger.error( control, treatment)
            raise
        return df
    

    def get_prop(self, df, group):
        """
        Calculate proportion from DataFrame for specified group.
        
        Args:
            df (pd.DataFrame): DataFrame containing the data.
            group (str): Group identifier ("control" or "treatment").
        
        Returns:
            str: "successes\tfailures\tproportion" or "0\t0\t-1" if no data.
        """
        successes, failures = 0, 0 

        if "successes" in df[df["group"]==group].columns:
            successes = sum(df[df["group"] == group]["successes"])  
        if "failures" in df[df["group"]==group].columns:
            failures =  sum(df[df["group"] == group]["failures"])
        p = -1

        if successes + failures > 0:
            p = successes / (successes + failures)        
        
        return ProportionResult(successes, failures, p) #"{}\t{}\t{}".format(successes, failures, p)



class TestTester(unittest.TestCase):


    def test_zeroes_Fisher(self):
        ctrl = {"successes": [0,0,0,0], "failures": [0,0,0,0]}
        treatment = {"successes": [0,0,0,0], "failures": [0,0,0,0]}

        tester = FisherTest()
        x = tester.test_sample(ctrl, treatment)
        self.assertEqual(x, TestResult(np.nan, np.nan, ProportionResult.failed(),ProportionResult.failed(), TestStatus.EMPTY, {}))


    def test_empty_Fisher(self):
        ctrl = {"successes": [], "failures": []}
        treatment = {"successes": [10,10,0,0], "failures": [10,10,0,0]}

        tester = FisherTest()
        x = tester.test_sample(ctrl, treatment)
        self.assertEqual(x, TestResult(np.nan, np.nan, ProportionResult.failed(),ProportionResult.failed(), TestStatus.EMPTY, {}))


    def test_zeroes_GLM(self):
        ctrl = {"successes": [0,0,0,0], "failures": [0,0,0,0]}
        treatment = {"successes": [0,0,0,0], "failures": [0,0,0,0]}

        tester = GLMModel()
        x = tester.test_sample(ctrl, treatment)
        self.assertEqual(x, TestResult(np.nan, np.nan, ProportionResult.failed(),ProportionResult.failed(), TestStatus.ZEROES, {}))


    def test_empty_GLM(self):
        ctrl = {"successes": [], "failures": []}
        treatment = {"successes": [10,10,0,0], "failures": [10,10,0,0]}

        tester = GLMModel()
        x = tester.test_sample(ctrl, treatment)
        self.assertEqual(x, TestResult(np.nan, np.nan, ProportionResult.failed(),ProportionResult.failed(), TestStatus.EMPTY, {}))


    def test_rejectMismatchedLengths_GLM(self):
        ctrl = {"successes": [10,20,30,40, 60], "failures": [10,20,30,40, 60]}
        treatment = {"successes": [10,20,30,40, 60], "failures": [10,20,30,40]}

        tester = GLMModel()
        x = tester.test_sample(ctrl, treatment)
        self.assertEqual(x, TestResult(np.nan, np.nan, ProportionResult.failed(),ProportionResult.failed(), TestStatus.MISMATCHED_LENGTH, {}))
    

    def test_rejectMismatchedLengths_Fisher(self):
        ctrl = {"successes": [10,20,30,40, 60], "failures": [10,20,30,40, 60]}
        treatment = {"successes": [10,20,30,40, 60], "failures": [10,20,30,40]}
        tester = GLMModel()
        x = tester.test_sample(ctrl, treatment)
        self.assertEqual(x, TestResult(np.nan, np.nan, ProportionResult.failed(),ProportionResult.failed(), TestStatus.MISMATCHED_LENGTH, {}))


    def test_get_prop_GM(self):
        ctrl = {"successes": [34, 31, 30, 28, 30], "failures": [4, 7, 2, 7, 6]}
        treatment = {"successes": [23, 29, 17, 31, 24], "failures": [15, 13, 23, 18, 14]}
        tester = GLMModel()

        if  x := tester.quality_check(ctrl, treatment):
            return x

        data = tester.load_data(ctrl, treatment)
        
        self.assertEqual(tester.get_prop(data, "control"), ProportionResult(153,26,0.8547486033519553))
        self.assertEqual(tester.get_prop(data, "treatment"), ProportionResult(124,83,0.5990338164251208))


    def test_shouldwork_GLM(self):
        # should work no zeroes
        ctrl = {"successes": [34, 31, 30, 28, 30], "failures": [4, 7, 2, 7, 6]}
        treatment = {"successes": [23, 29, 17, 31, 24], "failures": [15, 13, 23, 18, 14]}
        tester = GLMModel()
        x = tester.test_sample(ctrl, treatment)
        print(x)
        self.assertEqual(x, TestResult(pvalue=7.762584700721925e-08, statistics=0.2538782581305612,
                                        control_prop=ProportionResult(153,26,0.8547486033519553),
                                        treatment_prop=ProportionResult(124,83,0.5990338164251208),
                                        status=TestStatus.OK, any_data={}))
    
    def test_filterZeroRows_GLM_empty(self):
        # should return empty, to many zeroes
        ctrl = {"successes": [34, 31, 30, 28, 30], "failures": [4, 7, 2, 7, 6]}
        treatment = {"successes": [23, 29, 0, 0, 0], "failures": [15, 13, 0, 0, 0]}
        tester = GLMModel()
        x = tester.test_sample(ctrl, treatment)
        self.assertEqual(x,  TestResult(np.nan, np.nan, ProportionResult.failed(),ProportionResult.failed(), TestStatus.EMPTY, {}))

    def test_filterZeroRows_no_success(self):
        ctrl = { "failures": [4, 7, 2, 7, 6]}
        treatment = {"successes": [23, 29, 0, 0, 0], "failures": [15, 13, 0, 0, 0]}
        tester = GLMModel()
        x = tester.test_sample(ctrl, treatment)
        self.assertEqual(x,  TestResult(np.nan, np.nan, ProportionResult.failed(),ProportionResult.failed(), TestStatus.EMPTY, {}))

    



# UNSURE TO INCLUDE HOLD FOR NOW
"""class Chi2Test(StatTest):

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
if __name__ == '__main__':
    unittest.main()