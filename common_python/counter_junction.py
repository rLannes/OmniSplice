from enum import Enum
import unittest

# This enum is useless and could be replaced by just reading the header
class ReadJunction(Enum):
    SPLICED = 0
    UNSPLICED = 1
    CLIPPED = 2 
    EXON_OTHER = 3
    SKIPPED = 4
    WRONG_STRAND = 5
    E_ISOFORM = 6


class Counter():

    def __init__(self, cmd_successes, cmd_failures):
        self.successes = self.parse_args(cmd_successes)
        self.failures = self.parse_args(cmd_failures)


    def get_count(self, junction):
        
        successes = 0
        failures = 0

        for i in self.successes:
            successes += junction[i]

        for i in self.failures:
            failures += junction[i]

        return (successes, failures)


    def parse_args(self, args):
        res = []
        for e in args:
            if e == "spliced":
                res.append(ReadJunction.SPLICED.value)
            elif e == "unspliced":
                res.append(ReadJunction.UNSPLICED.value)
            elif e == "clipped":
                res.append(ReadJunction.CLIPPED.value)
            elif e == "exon_other":
                res.append(ReadJunction.EXON_OTHER.value)
            elif e == "skipped":
                res.append(ReadJunction.SKIPPED.value)
            elif e == "wrong_strand":
                res.append(ReadJunction.WRONG_STRAND.value)
            elif e == "e_isoform":
                res.append(ReadJunction.E_ISOFORM.value)
            else:
                raise AssertionError("unknown {}".format(e))
        return res                                
   

