from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

__author__ = 'ahbbollen'
from .parser import Exon, Record
from .generics import *

class ExonFactory(object):
    """
    This class defines a factory for producing Exons
    Takes the a record
    """
    def __init__(self, record):
        self.record = record

    def make(self):
        exons = []
        for i, (s, e) in enumerate(zip(self.record.exonStarts,
                                       self.record.exonEnds)):
            exons.append(Exon(self.record.gene, self.record.transcript,
                              self.record.chromosome, s, e, i))
        return exons


class RecordFactory(object):
    def __init__(self, *args):
        if len(args) == 1:
            self.items = args[0]._items
            self.line = args[0]._line
        elif len(args) > 1:
            raise TypeError("Too many arguments")
        else:
            self.line, self.items = empty_line()

    def check(self):
        for column in COLUMNS:
            try:
                _ = self.items[column]
            except KeyError:
                return False
        return True

    def setattribute(self, key, value):
        self.items[key] = value

    def set_gene(self, gene):
        assert isinstance(gene, basestring)

        self.items["geneName"] = gene

    def set_transcript(self, name):
        assert isinstance(name, basestring)
        self.items["name"] = name

    def set_chromosome(self, chrom):
        assert isinstance(chrom, basestring)

        self.items["chrom"] = chrom

    def set_transcription_start(self, txStart):
        assert isinstance(txStart, int)

        self.items["txStart"] = txStart

    def set_transcription_end(self, txEnd):
        assert isinstance(txEnd, int)

        self.items["txEnd"] = txEnd

    def set_cds_start(self, cdsStart):
        assert isinstance(cdsStart, int)

        self.items["cdsStart"] = cdsStart

    def set_cds_end(self, cdsEnd):
        assert isinstance(cdsEnd, int)

        self.items["cdsEnd"] = cdsEnd

    def set_exon_count(self, exonCount):
        assert isinstance(exonCount, int)

        self.items["exonCount"] = exonCount

    def set_exon_starts(self, exonStarts):
        assert isinstance(exonStarts, list)
        assert isinstance(exonStarts[0], int)
        assert len(exonStarts) == self.items["exonCount"]

        self.items["exonStarts"] = exonStarts

    def set_exon_ends(self, exonEnds):
        assert isinstance(exonEnds, list)
        assert isinstance(exonEnds[0], int)
        assert len(exonEnds) == self.items["exonCount"]

        self.items["exonEnds"] = exonEnds

    def make(self):
        normal_columns = set(COLUMNS) - set(NUMERIC_LIST_COLUMNS) # <-- remember, this is UNORDERED!
        line = [str(self.items[x]) for x in COLUMNS if x in normal_columns]
        line = "\t".join(line)
        for c in NUMERIC_LIST_COLUMNS:
            if isinstance(self.items[c], list) and all([isinstance(x, int) for x in self.items[c]]):
                line += "\t" + ",".join(map(str, self.items[c])) + ","
            elif isinstance(self.items[c], basestring) and self.items[c].endswith(","):
                line += "\t" + self.items[c]
            elif isinstance(self.items[c], basestring) and not self.items[c].endswith(","):
                line += "\t" + self.items[c] + ","
            else:
                raise ValueError
        r = Record(line, "internal")
        return r

