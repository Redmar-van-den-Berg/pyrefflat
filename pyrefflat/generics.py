from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *


__author__ = 'Sander Bollen'


COLUMNS = ["geneName", "name", "chrom", "strand", "txStart", "txEnd",
           "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds"]

NO_EXON_COLUMNS = ["geneName", "name", "chrom", "strand", "txStart", "txEnd",
           "cdsStart", "cdsEnd", "exonCount"]

NUMERIC_COLUMNS = ["txStart", "txEnd", "cdsStart", "cdsEnd"]

NUMERIC_LIST_COLUMNS = ["exonStarts", "exonEnds"]

STRING_COLUMNS = set(COLUMNS) - set(NUMERIC_COLUMNS) - set(NUMERIC_LIST_COLUMNS)

def empty_line():
    items = {}
    for n in NUMERIC_COLUMNS:
        items[n] = 0
    for s in STRING_COLUMNS:
        items[s] = "undefined"
    for l in NUMERIC_LIST_COLUMNS:
        items[l] = [0]

    line = "\t".join([str(items[x]) for x in set(COLUMNS) - set(NUMERIC_LIST_COLUMNS)])
    for c in NUMERIC_LIST_COLUMNS:
        line += "\t" + ",".join(map(str, items[c])) + ","

    return line, items

