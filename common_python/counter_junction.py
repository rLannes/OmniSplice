from enum import Enum
import unittest

# This enum is useless and could be replaced by just reading the header
class ReadJunction(Enum):
    SPLICED = 0
    UNSPLICED = 1
    CLIPPED = 2 
    EXON_INTRON = 3
    EXON_OTHER = 4
    SKIPPED = 5
    WRONG_STRAND = 6
    E_ISOFORM = 7


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
            elif e == "exon_intron":
                res.append(ReadJunction.EXON_INTRON.value)
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
   


class TestCounter(unittest.TestCase):

    def test_1(self):
        donnor = [0,88,0,0,0,1,1228,0]
        acceptor = [0, 9, 0, 0, 1, 1, 1305]
        counter = Counter(["SPLICED".lower()], ["UNSPLICED".lower(), "EXON_OTHER".lower(), "EXON_INTRON".lower()])
        v = counter.get_count(donnor, acceptor)
        self.assertEqual(v, (0, 98))

    def test_2(self):
        donnor = [15,88,0,0,0,1,1228,0]
        acceptor = [15, 9, 0, 0, 1, 1, 1305]
        counter = Counter(["SPLICED".lower()], ["UNSPLICED".lower(), "EXON_OTHER".lower(), "EXON_INTRON".lower()])
        v = counter.get_count(donnor, acceptor)
        self.assertEqual(v, (15, 98))

    def test_3(self):
        donnor = [15,5,3,2,1,1,1228,0]
        acceptor = [15, 9, 5, 7, 1, 1, 1305]
        counter = Counter(["SPLICED".lower()], ["UNSPLICED".lower(), "EXON_OTHER".lower(), "EXON_INTRON".lower()])
        v = counter.get_count(donnor, acceptor)
        self.assertEqual(v, (15, 25))


"""
2L      FBgn0250906     FBtr0330716     exon_1  false   -       2750078 2748792 Donnor  0       88      0       0       0       1       1228    0
2L      FBgn0250906     FBtr0330716     exon_2  false   -       2748792 2750078 Acceptor        0       9       0       0       1       1       1305    9
2L      FBgn0250906     FBtr0330716     exon_2  false   -       2748319 2748255 Donnor  565     15      0       0       0       1       1361    0
2L      FBgn0250906     FBtr0330716     exon_3  false   -       2748255 2748319 Acceptor        565     36      0       0       0       1       1371    0
2L      FBgn0250906     FBtr0330716     exon_3  false   -       2747310 2747252 Donnor  573     2       0       0       0       1       1348    0
2L      FBgn0250906     FBtr0330716     exon_4  false   -       2747252 2747310 Acceptor        573     12      0       0       0       1       1301    0
2L      FBgn0250906     FBtr0330716     exon_4  true    -       2746879 .       Donnor  0       213     0       0       0       1       1451    0
2L      FBgn0032117     FBtr0346623     exon_1  false   -       9430696 9430429 Donnor  352     222     0       0       0       4       398     0
2L      FBgn0032117     FBtr0346623     exon_2  false   -       9430429 9430696 Acceptor        352     59      0       611     17      1       186     0
2L      FBgn0032117     FBtr0346623     exon_2  false   -       9430022 9429589 Donnor  1562    111     5       3       15      1       134     0
2L      FBgn0032117     FBtr0346623     exon_3  true    -       9429589 9430022 Acceptor        1562    26438   0       516     0       15      135     0
2L      FBgn0032117     FBtr0346623     exon_3  false   -       9429411 9429345 Donnor  20071   552     1       0       99      1       129     0
2L      FBgn0032117     FBtr0346623     exon_4  true    -       9429345 9429411 Acceptor        20071   40      34      84      3       1       126     0
2L      FBgn0032117     FBtr0346623     exon_4  true    -       9428668 .       Donnor  0       54      0       0       0       14      7058    0

"""