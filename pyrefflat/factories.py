__author__ = 'ahbbollen'

from .parser import Exon, Record
from .generics import COLUMNS, NO_EXON_COLUMNS

class ExonFactory(object):
    """
    This class defines a factory for producing Exons
    Takes the items dictionary of a record
    """
    def __init__(self, items):
        self.items = items

    def make(self):
        exons = []
        starts = self.items['exonStarts'].split(',')
        if starts[-1] == '':
            starts.pop()
        ends = self.items['exonEnds'].split(',')
        if ends[-1] == '':
            ends.pop()
        assert len(starts) == len(ends), "Supplied items don't match"
        for i, (s, e) in enumerate(zip(starts, ends)):
            exons.append(Exon(self.items['geneName'], self.items['name'], self.items['chrom'], s, e, i))
        return exons


class RecordFactory(object):
    def __init__(self, items, filename):
        self.items = items
        self._filename = filename
        self.check()

    def check(self):
        for column in COLUMNS:
            try:
                _ = self.items[column]
            except KeyError:
                raise KeyError("Malformed dictionary! Couldn't find required element {0}".format(column))

    def make(self):
        exons = ExonFactory(self.items).make()
        line = [self.items[x] for x in NO_EXON_COLUMNS]
        exonstarts = [x.start for x in exons]
        exonends = [x.stop for x in exons]
        exonsline = ",".join(exonstarts) + ","
        exonstline = ",".join(exonends) + ","
        line += [exonsline]
        line += [exonstline]
        r = Record("\t".join(line), self._filename)
        return r


