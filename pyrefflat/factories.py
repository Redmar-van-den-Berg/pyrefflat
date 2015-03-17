__author__ = 'ahbbollen'

from .parser import Exon, Record
from .generics import *

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
    def __init__(self, *args):
        if len(args) == 1:
            self.items = args[0]
            self.line = None
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

    def make(self):
        normal_columns = set(COLUMNS) - set(NUMERIC_LIST_COLUMNS) # <-- remember, this is UNORDERED!
        line = [str(self.items[x]) for x in COLUMNS if x in normal_columns]
        line = "\t".join(line)
        for c in NUMERIC_LIST_COLUMNS:
            line += "\t" + ",".join(map(str, self.items[c])) + ","
        r = Record(line, "internal")
        return r

