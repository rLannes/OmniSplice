import re
from pathlib import Path
from .counter_junction import Counter
import unittest


class Junction():

    def __init__(self, spt, header, basename, genotype):


        self.contig = spt[header["Contig"]]
        self.strand = spt[header["Strand"]]

        self.pos =  spt[header["Donnor"]] if spt[header["Strand"]] == "+" else spt[header["Acceptor"]] 
        self.next = spt[header["Acceptor"]] if spt[header["Strand"]] == "+" else spt[header["Donnor"]]


        self.ambigious = True if spt[header["Ambiguous"]] == "true" else False

        self.gene = [spt[header["Gene"]]]
        self.transcript = [spt[header["Transcript"]]]
        self.intron_number = [spt[header["Intron"]]]
        
        self.count = {genotype:  {basename:{}}}
        self.count[genotype][basename] = list(map( int, spt[8:]))
    

    def update_count(self, spt, genotype, header, basename):

        self.ambigious = True if (spt[header["Ambiguous"]] == "true" or self.ambigious == True) else self.ambigious
        if genotype not in self.count:
            self.count[genotype] = {basename: list(map( int, spt[8:]))}
        elif basename not in self.count[genotype]:
            self.count[genotype][basename] =  list(map( int, spt[8:]))

        if spt[header["Transcript"]] not in self.transcript:
            self.gene.append(spt[header["Gene"]])
            self.transcript.append(spt[header["Transcript"]])
            self.intron_number.append(spt[header["Intron"]])

    def __hash__(self):
        hash((self.contig, self.strand, self.pos, self.next)) # , self.pos, self.next,

    def get_hash_key(self):
        return (self.contig, self.strand, self.pos, self.next)


    def dump(self):
        g = []
        for i, v in enumerate(self.gene):
            g.append("{}_{}_{}".format(v, self.transcript[i], str(self.intron_number[i])))
        l = [self.contig, self.strand, ";".join(g),
            self.pos, self.next]
        
        return l


    def get_junction_count(self, counter: Counter):

        dico_count = {}
        for genotype, file_dict in self.count.items():
            if genotype not in dico_count:
                dico_count[genotype] = {"successes": [], "failures": []} # = {"control": [], "treatment": [], "file": []}
            for file, data in file_dict.items():
                try:
                    (sucesses, failures) = counter.get_count(junction=data)
                except:
                    print(self.__dict__)
                    raise
                dico_count[genotype]["successes"].append(sucesses)
                dico_count[genotype]["failures"].append(failures)
    
        return dico_count

    def pass_threshold(self, ):
        """
            manage QC filtering!
        """
        pass


class TestCounter(unittest.TestCase):

    def test_1(self):
        header = "contig  gene_name       transcript_name exon_number     ambiguous       strand  pos     next    exon_type       spliced unspliced       clipped exon_intron     exon_other      skipped wrong_strand    e_isoform".split()
        header = dict((k, i) for i, k in enumerate(header))

        spt = "3R      FBgn0037773     FBtr0082215     exon_2  true    -       10000572        10000368        Donnor  222     19      1       0       0       1       33   548".split()      
        pos_1 =  spt[header["pos"]] if spt[header["exon_type"]] == "Donnor" else spt[header["next"]] 
        next_1 = spt[header["next"]] if spt[header["exon_type"]] == "Donnor" else spt[header["pos"]]
        spt = "3R      FBgn0037773     FBtr0082215     exon_3  false   -       10000368        10000572        Acceptor        222     8       0       0       0       1       31      0".split()
        pos_2 =  spt[header["pos"]] if spt[header["exon_type"]] == "Donnor" else spt[header["next"]] 
        next_2 = spt[header["next"]] if spt[header["exon_type"]] == "Donnor" else spt[header["pos"]]
        self.assertEqual(pos_1, pos_2)
        self.assertEqual(next_1, next_2)
        print(pos_1, pos_2, next_1, next_2)

    def test2(self):
        header = "contig  gene_name       transcript_name exon_number     ambiguous       strand  pos     next    exon_type       spliced unspliced       clipped exon_intron     exon_other      skipped wrong_strand    e_isoform".split()
        header = dict((k, i) for i, k in enumerate(header))

        X = """3R      FBgn0037773     FBtr0082215     exon_2  true    -       10000572        10000368        Donnor  222     19      1       0       0       1       33      548""".split()
        X2 = """3R      FBgn0037773     FBtr0082215     exon_3  false   -       10000368        10000572        Acceptor        222     8       0       0       0       1       31      0""".split()
        j = Junction(spt=X, header=header, basename="f1", genotype="control")
        j.update_count(X2, genotype="control", header=header, basename="f1")
        print(j.__dict__)

    def test3(self):
        header = "contig  gene_name       transcript_name exon_number     ambiguous       strand  pos     next    exon_type       spliced unspliced       clipped exon_intron     exon_other      skipped wrong_strand    e_isoform".split()
        header = dict((k, i) for i, k in enumerate(header))

        X = """3R      FBgn0037773     FBtr0082215     exon_2  true    -       10000572        10000368        Donnor  222     19      1       0       0       1       33      548""".split()
        X2 = """3R      FBgn0037773     FBtr0082215     exon_3  false   -       10000368        10000572        Acceptor        222     8       0       0       0       1       31      0""".split()
        j = Junction(spt=X, header=header, basename="f2", genotype="treatment")
        j.update_count(X2, genotype="treatment", header=header, basename="f2")
        j.update_count(X, genotype="control", header=header, basename="f1")
        j.update_count(X2, genotype="control", header=header, basename="f1")
        print(j.__dict__)

# spt, genotype, header, basename


        # def __init__(self, spt, header, basename, genotype):

"""     
    def get_count(self, type):
    
        type = [Acceptor, Donnor, Both]
        
        control = []
        treatment = []

        if type != "Both":

            for sub_count in self.count["control"][type]:
                main_ = 0
                other = 0
                for i, v in enumerate(sub_count):
                    if i == tobecompared:
                        main_ += v
                    else: 
                        other += v
                control.append((main_, other))
            for sub_count in self.count["control"][type]:
                main_ = 0
                other = 0
                for i, v in enumerate(sub_count):
                    if i == tobecompared:
                        main_ += v
                    else: 
                        other += v
                treatment.append((main_, other))
        
        return (control, treatment)
            

  """

