__author__ = 'ahbbollen'


COLUMNS = ["geneName", "name", "chrom", "strand", "txStart", "txEnd",
           "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds"]

NO_EXON_COLUMNS = ["geneName", "name", "chrom", "strand", "txStart", "txEnd",
           "cdsStart", "cdsEnd", "exonCount"]

NUMERIC_COLUMNS = ["txStart", "txEnd", "cdsStart", "cdsEnd"]
STRING_COLUMNS = set(COLUMNS) - set(NUMERIC_COLUMNS)

def empty_line():
    items = {}
    for n in NUMERIC_COLUMNS:
        items[n] = 0
    for s in STRING_COLUMNS:
        items[s] = "undefined"

    line = "\t".join([items[x] for x in COLUMNS])

    return line

