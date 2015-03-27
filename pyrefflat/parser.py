from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

__author__ = 'ahbbollen'

import locale

COLUMNS = ["geneName", "name", "chrom", "strand", "txStart", "txEnd",
           "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds"]

ENCODING = locale.getdefaultlocale()[1]


class Reader(object):
    def __init__(self, filename):
        self._filename = filename
        self._handler = open(filename, 'rb')

    def __iter__(self):
        return self

    def next(self):
        line = next(self._handler)
        return Record(line, self._filename)

    # python 3 compatibility
    def __next__(self):
        return self.next()

    def close(self):
        self._handler.close()


class Writer(object):
    def __init__(self, filename):
        self._filename = filename
        self._handler = open(filename, 'wb')

    def write(self, record):
        self._handler.write((record._line + "\n").encode(ENCODING))

    def close(self):
        self._handler.close()


class Record(object):
    def __init__(self, line, filename):
        self._line = line.decode(ENCODING)
        self.filename = filename
        self._parse_line()

    def _parse_line(self):
        self._raw_items = self._line.strip().split('\t')
        assert len(self._raw_items) >= 11, "Contains less than 11 columns!"
        self._items = dict()
        for i in range(11):
            self._items[COLUMNS[i]] = self._raw_items[i]

    @property
    def gene(self):
        return str(self._items["geneName"])

    @property
    def transcript(self):
        return str(self._items["name"])

    @property
    def chromosome(self):
        return str(self._items["chrom"])

    @property
    def strand(self):
        return str(self._items["strand"])

    @property
    def txStart(self):
        return int(self._items["txStart"])

    @property
    def txEnd(self):
        return int(self._items["txEnd"])

    @property
    def cdsStart(self):
        return int(self._items["cdsStart"])

    @property
    def cdsEnd(self):
        return int(self._items["cdsEnd"])

    @property
    def n_exons(self):
        return int(self._items["exonCount"])

    @property
    def exonStarts(self):
        assert str(self._items["exonStarts"]).endswith(","), "Malformed refFlat line!"
        starts = str(self._items["exonStarts"]).split(",")
        # remove final unneccesary comma
        if starts[-1] == '':
            _ = starts.pop()
        return list(map(int, starts))

    @property
    def exonEnds(self):
        assert str(self._items["exonEnds"]).endswith(","), "Malformed refFlat line!"
        ends = str(self._items["exonEnds"]).split(",")
        # remove final unneccessary comma
        if ends[-1] == '':
            _ = ends.pop()
        return list(map(int, ends))

    @property
    def exons(self):
        return ExonFactory(self._items).make()


class Exon(object):
    """
    This class defines an exon inside a record
    """
    def __init__(self, gene, transcript, chr, start, stop, n):
        self._gene = gene
        self._transcript = transcript
        self._chr = chr
        self._start = start
        self._end = stop
        self._number = n

    @property
    def gene(self):
        return self._gene

    @property
    def transcript(self):
        return self._transcript

    @property
    def chr(self):
        return self._chr

    @property
    def start(self):
        return self._start

    @property
    def stop(self):
        return self._end

    @property
    def number(self):
        return self._number


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