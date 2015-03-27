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

