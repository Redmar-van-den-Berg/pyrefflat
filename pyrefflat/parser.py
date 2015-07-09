from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

__author__ = 'ahbbollen'

import locale
import os.path

import logging

from .generics import *
from .models import *

ENCODING = locale.getdefaultlocale()[1]


class Reader(object):
    def __init__(self, filename):
        self._filename = os.path.basename(filename)
        self._handler = open(filename, 'rb')

    def __iter__(self):
        return self

    def next(self):
        try:
            line = next(self._handler)
        except ValueError:
            raise StopIteration
        return Record.fromline(line)

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
        nline = record.line + "\n"
        self._handler.write(nline.encode(ENCODING))

    def close(self):
        self._handler.close()


class Record(object):
    def __init__(self, geneName, name, chrom, strand, txStart, txEnd,
                 cdsStart, cdsEnd, exonCount, exonStarts, exonEnds):
        self._gene = geneName
        self._tr_name = name
        self._chr = chrom
        self._strand = strand
        self._tx_start = txStart
        self._tx_end = txEnd
        self._cds_start = cdsStart
        self._cds_end = cdsEnd
        self._exon_count = exonCount
        self._exon_start = exonStarts
        self._exon_ends = exonEnds

    @property
    def gene(self):
        return str(self._gene)

    @property
    def transcript(self):
        return str(self._tr_name)

    @property
    def chromosome(self):
        return str(self._chr)

    @property
    def strand(self):
        return str(self._strand)

    @property
    def txStart(self):
        return int(self._tx_start)

    @property
    def txEnd(self):
        return int(self._tx_end)

    @property
    def cdsStart(self):
        return int(self._cds_start)

    @property
    def cdsEnd(self):
        return int(self._cds_end)

    @property
    def n_exons(self):
        return int(self._exon_count)

    @property
    def exonStarts(self):
        return [int(x) for x in self._exon_start]

    @property
    def exonEnds(self):
        return [int(x) for x in self._exon_ends]

    @property
    def exons(self):
        return Exon.fromrecord(self)

    def to_dict(self):
        d = {}
        d["geneName"] = self.gene
        d["name"] = self.transcript
        d["chrom"] = self.chromosome
        d["strand"] = self.strand
        d["txStart"] = self.txStart
        d["txEnd"] = self.txEnd
        d["cdsStart"] = self.cdsStart
        d["cdsEnd"] = self.cdsEnd
        d["exonStarts"] = self.exonStarts
        d["exonEnds"] = self.exonEnds
        d["exonCount"] = self.n_exons

        return d

    @property
    def line(self):
        line = []
        d = self.to_dict()
        for nlc in NUMERIC_LIST_COLUMNS:
            d[nlc] = ",".join(map(str, d[nlc])) + ","
        for col in COLUMNS:
            line += [d[col]]
        return "\t".join(map(str, line))



    @classmethod
    def fromdict(cls, items):
        """
        Builds a record from a dictionary.
        This dictionary must contain all fields specified in generics.COLUMNS
        """
        normal_columns = set(COLUMNS) - set(NUMERIC_LIST_COLUMNS) # <-- remember, this is UNORDERED!

        # first check whether all columns are there and properly formatted
        for c in COLUMNS:
            if c not in items:
                raise ValueError("Item {c} must be given".format(c=c))
        for nlc in NUMERIC_LIST_COLUMNS:
            if not isinstance(items[nlc], list):
                raise ValueError("Item {nlc} must be a list of integers".format(nlc=nlc))
            elif not all([isinstance(x, int) for x in items[nlc]]):
                raise ValueError("Item {nlc} must be a list of integers".format(nlc=nlc))
        for nc in NUMERIC_COLUMNS:
            if not isinstance(items[nc], int):
                raise ValueError("Item {nc} must be an integer".format(nc=nc))

        #
        #
        # line = [str(items[x]) for x in COLUMNS if x in normal_columns]
        # line = "\t".join(line)
        # for c in NUMERIC_LIST_COLUMNS:
        #     if isinstance(items[c], list) and all([isinstance(x, int) for x in items[c]]):
        #         line += "\t" + ",".join(map(str, items[c])) + ","
        #     elif isinstance(items[c], basestring) and items[c].endswith(","):
        #         line += "\t" + items[c]
        #     elif isinstance(items[c], basestring) and not items[c].endswith(","):
        #         line += "\t" + items[c] + ","
        #     else:
        #         raise ValueError
        r = Record(**items)
        return r

    @classmethod
    def fromline(cls, line):
        """
        Builds a record from a line
        """
        raw_items = line.strip().split('\t')
        assert len(raw_items) >= 11, "Contains less than 11 columns!"
        items = dict()
        for i in xrange(11):
            items[COLUMNS[i]] = raw_items[i]
        for nc in NUMERIC_COLUMNS:
            items[nc] = int(items[nc])
        for lnc in NUMERIC_LIST_COLUMNS:
            if not items[lnc].endswith(','):
                raise ValueError("Malformed refFlat file! Value {lnc} must end in a comma".format(lnc))

            it = items[lnc].split(',')
            it.pop()
            items[lnc] = [int(x) for x in it]


        r = Record(**items)
        return r


class RefFlatProcessor(object):
    def __init__(self, filename, log=True, log_level="INFO"):
        self.filename = filename
        self._already_processed = False
        self.genes = {}
        self.transcripts = {}
        self.n_duplicates = 0
        self.duplicates = []
        self.log = log
        if log:
            self.logger = logging.getLogger("RefFlatProcessor")
            numeric_level = getattr(logging, log_level.upper(), None)
            self.logger.setLevel(numeric_level)
            ch = logging.StreamHandler()
            ch.setLevel(numeric_level)
            formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            ch.setFormatter(formatter)
            self.logger.addHandler(ch)
            self.logger.info("Initializing...")

    def process(self, remove_duplicates=True, flush=True):
        """
        Create gene and transcript tables for the input refflat file
        :param remove_duplicates: Boolean. True by default.
        If true, will remove duplicates from the tables.
        Duplicate is defined when the transcript name occurs more than one (aka transcript should be unique)
        Since we have no way of knowing what transcript is the correct one,
        the retained transcript is simply the first we encounter.
        Conflicting records will be placed in the `duplicates` list of this object
        Please note that if this variable is set to False, derived gene coordinates may no longer make sense
        (Yes, there are transcripts that are annotated on multiple chromosomes even!)
        :param flush: flush contents of object first
        :return: genes as dict to self.genes, transcript as dict to self.transcripts
        """

        if flush:
            self._already_processed = False
            self.genes = {}
            self.transcripts = {}
            self.duplicates = []
            self.n_duplicates = 0

        for i, record in enumerate(Reader(self.filename)):
            if i % 1000 == 0 and self.log:
                self.logger.info("Processed {0} records".format(i))
            if not record.gene in self.genes and not record.transcript in self.transcripts:
                tmptranscript = Transcript(record.transcript, record.chromosome, record.txStart, record.txEnd,
                                           record.cdsStart, record.cdsEnd, exons=record.exons, strand=record.strand)
                tmpgene = Gene(record.gene)
                tmpgene.update_transcripts(tmptranscript)
                tmptranscript.gene = tmpgene
                self.genes[record.gene] = tmpgene
                self.transcripts[tmptranscript.name] = tmptranscript

            elif record.gene in self.genes and not record.transcript in self.transcripts:
                tmptranscript = Transcript(record.transcript, record.chromosome, record.txStart, record.txEnd,
                                           record.cdsStart, record.cdsEnd, exons=record.exons, strand=record.strand)
                self.genes[record.gene].update_transcripts(tmptranscript)
                tmptranscript.gene = self.genes[record.gene]
                self.transcripts[tmptranscript.name] = tmptranscript

            elif record.gene in self.genes and record.transcript in self.transcripts:
                self.n_duplicates += 1
                self.duplicates += [record]
                if not remove_duplicates:
                    tmptranscript = Transcript(record.transcript, record.chromosome, record.txStart, record.txEnd,
                                               record.cdsStart, record.cdsEnd, exons=record.exons, strand=record.strand)
                    self.genes[record.gene].update_transcripts(tmptranscript)
                    tmptranscript.gene = self.genes[record.gene]
                    self.transcripts[tmptranscript.name] = tmptranscript

            else:
                # would happen in record.transcript in transcripts and not record.gene in genes
                # which should not be possible
                raise ValueError("Something odd happened!")

        self._already_processed = True
        if self.log:
            self.logger.info("Done processing")