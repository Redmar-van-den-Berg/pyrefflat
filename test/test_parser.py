from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *
__author__ = 'ahbbollen'

import pytest
from pyrefflat.parser import *

class TestReader:
    def test_return(self):
        r = Reader(filename="data/mini.refFlat")
        re = next(r)
        assert isinstance(re, Record)

@pytest.fixture(scope="module")
def record():
    r = Reader(filename="data/mini.refFlat")
    return next(r)

@pytest.fixture(scope="module")
def exons(record):
    return record.exons

class TestRecord():

    def test_filename(self, record):
        assert record.filename == "mini.refFlat"

    def test_raw_items(self, record):
        assert len(record._raw_items) == 11

    def test_gene_type(self, record):
        assert isinstance(record.gene, basestring)

    def test_name_type(self, record):
        assert isinstance(record.transcript, basestring)

    def test_chrom_type(self, record):
        assert isinstance(record.chromosome, basestring)

    def test_strand_type(self, record):
        assert isinstance(record.strand, basestring)

    def test_txstart_type(self, record):
        assert isinstance(record.txStart, int)

    def test_txend_type(self, record):
        assert isinstance(record.txEnd, int)

    def test_cdsstart_type(self, record):
        assert isinstance(record.cdsStart, int)

    def test_cdsend_type(self, record):
        assert isinstance(record.cdsEnd, int)

    def test_exons_type(self, record):
        assert isinstance(record.exons, list)

    def test_nexons_type(self, record):
        assert isinstance(record.n_exons, int)

    def test_individual_exon_type(self, record):
        for ex in record.exons:
            assert isinstance(ex, Exon)

    def test_exonsstarts_type(self, record):
        assert isinstance(record.exonStarts, list)

    def test_ind_exonsstarts_type(self, record):
        for s in record.exonStarts:
            assert isinstance(s, int)

    def test_exonends_type(self, record):
        assert isinstance(record.exonEnds, list)

    def test_ind_exonends_type(self, record):
        for e in record.exonEnds:
            assert isinstance(e, int)

    def test_gene_contents(self, record):
        assert (record.gene == "MLH1")

    def test_name_contents(self, record):
        assert (record.transcript == "NM_000249")

    def test_chromosome_contents(self, record):
        assert (record.chromosome == "chr3")

    def test_strand_contents(self, record):
        assert (record.strand == "+")

    def test_txstart_contents(self, record):
        assert (record.txStart == 37034840)

    def test_txend_contents(self, record):
        assert (record.txEnd == 37092337)

    def test_cdsstart_contents(self, record):
        assert (record.cdsStart == 37035038)

    def test_cdsend_contents(self, record):
        assert (record.cdsEnd == 37092144)

    def test_nexons_contents(self, record):
        assert (record.n_exons == 19)

    def test_exonstarts_contents(self, record):
        assert (record.exonStarts == [37034840, 37038109, 37042445, 37045891,
                                      37048481, 37050304, 37053310, 37053501,
                                      37055922, 37058996, 37061800, 37067127,
                                      37070274, 37081676, 37083758, 37089009,
                                      37090007, 37090394, 37091976])

    def test_exonends_contents(self, record):
        assert (record.exonEnds == [37035154, 37038200, 37042544, 37045965,
                                    37048554, 37050396, 37053353, 37053590,
                                    37056035, 37059090, 37061954, 37067498,
                                    37070423, 37081785, 37083822, 37089174,
                                    37090100, 37090508, 37092337])

class TestExon():

    def test_gene_type(self, exons):
        for ex in exons:
            assert isinstance(ex.gene, basestring)

    def test_name_type(self, exons):
        for ex in exons:
            assert isinstance(ex.transcript, basestring)

    def test_chromosome_type(self, exons):
        for ex in exons:
            assert isinstance(ex.chr, basestring)

    def test_start_type(self, exons):
        for ex in exons:
            assert isinstance(ex.start, int)

    def test_end_type(self, exons):
        for ex in exons:
            assert isinstance(ex.stop, int)

    def test_n_type(self, exons):
        for ex in exons:
            assert isinstance(ex.number, int)

    def test_gene_contents(self, exons):
        for ex in exons:
            assert (ex.gene == "MLH1")

    def test_name_contents(self, exons):
        for ex in exons:
            assert (ex.transcript == "NM_000249")

    def test_chromosome_contents(self, exons):
        for ex in exons:
            assert (ex.chr == "chr3")

    def test_start_contents(self, exons):
        starts = [37034840, 37038109, 37042445, 37045891,
                  37048481, 37050304, 37053310, 37053501,
                  37055922, 37058996, 37061800, 37067127,
                  37070274, 37081676, 37083758, 37089009,
                  37090007, 37090394, 37091976]
        for i, ex in enumerate(exons):
            assert (ex.start == starts[i])

    def test_stop_contents(self, exons):
        ends = [37035154, 37038200, 37042544, 37045965,
                37048554, 37050396, 37053353, 37053590,
                37056035, 37059090, 37061954, 37067498,
                37070423, 37081785, 37083822, 37089174,
                37090100, 37090508, 37092337]

        for i, ex in enumerate(exons):
            assert (ex.stop == ends[i])

    def test_number_contents(self, exons):
        for i, ex in enumerate(exons):
            assert (ex.number == i)