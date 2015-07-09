from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

from generics import *

class Exon(object):
    """
    This class defines an exon inside a record
    """

    __slots__ = ["_gene", "_transcript", "_chr", "_start", "_end", "_number"]

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

    @classmethod
    def fromrecord(cls, record):
        exons = []
        for i, (s, e) in enumerate(zip(record.exonStarts,
                                       record.exonEnds)):
            exons.append(Exon(record.gene, record.transcript,
                              record.chromosome, s, e, i))
        return exons


class Transcript(object):
    __slots__ = ["name", "gene", "chr", "start", "end", "cds_start", "cds_end", "exons", "strand"]

    def __init__(self, name, chr, start, end, cds_start, cds_end, exons=None, gene=None, strand="+"):
        self.name = name
        self.gene = gene
        self.chr = chr
        self.start = start
        self.end = end
        self.cds_start = cds_start
        self.cds_end = cds_end
        self.exons = exons
        self.strand = strand

    def update_exons(self, exon):
        if exon.start < self.start:
            raise ValueError("Start of exon cannot be in front of start of transcript")
        if exon.end < self.end:
            raise ValueError("End of exon cannot be behind end of transcript")

        if self.exons:
            self.exons.append(exon)
        else:
            self.exons = [exon]

    @property
    def line(self):
        line = []
        d = self.to_dict()
        for nlc in NUMERIC_LIST_COLUMNS:
            d[nlc] = ",".join(map(str, d[nlc])) + ","
        for col in COLUMNS:
            line += [d[col]]
        return "\t".join(map(str, line))

    def to_dict(self):
        d = {}
        d["geneName"] = self.gene.name
        d["name"] = self.name
        d["chrom"] = self.chr
        d["strand"] = self.strand
        d["txStart"] = self.start
        d["txEnd"] = self.end
        d["cdsStart"] = self.cds_start
        d["cdsEnd"] = self.cds_end
        d["exonStarts"] = [int(x.start) for x in self.exons]
        d["exonEnds"] = [int(x.stop) for x in self.exons]
        d["exonCount"] = len(self.exons)

        return d


class Gene(object):
    __slots__ = ["name", "min_coord", "max_coord", "transcripts", "chr"]

    def __init__(self, name, chr=None, min_coord=None, max_coord=None, transcripts=None):
        self.name = name
        self.min_coord = min_coord
        self.max_coord = max_coord
        self.transcripts = transcripts
        self.chr = chr

    def update_transcripts(self, transcript):
        if self.min_coord:
            if transcript.start < self.min_coord:
                self.min_coord = transcript.start
        else:
            self.min_coord = transcript.start

        if self.max_coord:
            if transcript.end > self.max_coord:
                self.max_coord = transcript.end
        else:
            self.max_coord = transcript.end

        if self.transcripts:
            self.transcripts += [transcript]
            self.chr += [transcript.chr]
        else:
            self.transcripts = [transcript]
            self.chr = [transcript.chr]

