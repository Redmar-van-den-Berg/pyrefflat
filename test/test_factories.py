from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

__author__ = 'ahbbollen'

import pytest

from pyrefflat import Reader
from pyrefflat.factories import *

@pytest.fixture(scope="module")
def factory():
    factory = RecordFactory()
    return factory

@pytest.fixture(scope="module")
def record():
    r = Reader(filename="test/data/mini.refFlat")
    return next(r)

@pytest.fixture(scope="module")
def from_record(record):
    factory = RecordFactory(record)
    return factory


class Test_RecordFactory():
    def test_creation_empty(self):
        new_factory = RecordFactory()
        for k, v in new_factory.items.iteritems():
            assert v is not None

    def test_too_many_args_creation(self):
        with pytest.raises(TypeError):
            _ = RecordFactory("one", "two")

    def test_creation_from_record_basic(self, record):
        new_factory = RecordFactory(record)
        for k, v in new_factory.items.iteritems():
            assert v is not None

    def test_from_record_gene(self, from_record):
        nrecord = from_record.make()
        assert nrecord.gene == "MLH1"

    def test_from_record_name(self, from_record):
        nrecord = from_record.make()
        assert nrecord.transcript == "NM_000249"

    def test_from_record_chr(self, from_record):
        nrecord = from_record.make()
        assert nrecord.chromosome == "chr3"

    def test_from_record_strand(self, from_record):
        assert from_record.make().strand == "+"

    def test_from_record_txstart(self, from_record):
        assert from_record.make().txStart == 37034840

    def test_from_record_texend(self, from_record):
        assert from_record.make().txEnd == 37092337

    def test_from_record_cdsstart(self, from_record):
        assert from_record.make().cdsStart == 37035038

    def test_from_record_cdsend(self, from_record):
        assert from_record.make().cdsEnd == 37092144

    def test_from_record_nexons(self, from_record):
        assert from_record.make().n_exons == 19

    def test_from_record_exonstarts(self, from_record):
        assert from_record.make().exonStarts == [37034840, 37038109, 37042445, 37045891,
                                      37048481, 37050304, 37053310, 37053501,
                                      37055922, 37058996, 37061800, 37067127,
                                      37070274, 37081676, 37083758, 37089009,
                                      37090007, 37090394, 37091976]

    def test_from_record_exonends(self, from_record):
        assert from_record.make().exonEnds == [37035154, 37038200, 37042544, 37045965,
                                    37048554, 37050396, 37053353, 37053590,
                                    37056035, 37059090, 37061954, 37067498,
                                    37070423, 37081785, 37083822, 37089174,
                                    37090100, 37090508, 37092337]

    def test_empty(self, factory):
        for k, v in factory.items.iteritems():
            assert v is not None

    def test_set_attribute(self, factory):
        factory.setattribute("geneName", "MLH1")
        assert factory.items["geneName"] == "MLH1"

    def test_set_attribute_make_record(self, factory):
        factory.setattribute("geneName", "MLH1")
        record = factory.make()
        assert record.gene == "MLH1"

    def test_set_gene(self, factory):
        factory.set_gene("MLH1")
        assert factory.make().gene == "MLH1"

    def test_set_transcript(self, factory):
        factory.set_transcript("XM_test")
        assert factory.make().transcript == "XM_test"

    def test_set_chromosome(self, factory):
        factory.set_chromosome("chr3")
        assert factory.make().chromosome == "chr3"

    def test_set_transcription_start(self, factory):
        factory.set_transcription_start(100)
        assert factory.make().txStart == 100

    def test_set_transcription_end(self, factory):
        factory.set_transcription_end(1000)
        assert factory.make().txEnd == 1000

    def test_set_cds_start(self, factory):
        factory.set_cds_start(100)
        assert factory.make().cdsStart == 100

    def test_set_cds_end(self, factory):
        factory.set_cds_end(1000)
        assert factory.make().cdsEnd == 1000

    def test_set_exoncount(self, factory):
        factory.set_exon_count(3)
        assert factory.make().n_exons == 3

    def test_set_exonstarts(self, factory):
        factory.set_exon_starts([1, 2, 3])
        assert factory.make().exonStarts == [1, 2, 3]

    def test_exonends(self, factory):
        factory.set_exon_ends([2, 3, 4])
        assert factory.make().exonEnds == [2, 3, 4]

    def test_strand(self, factory):
        factory.set_strand("+")
        assert factory.make().strand == "+"

    def test_except_set_gene(self, factory):
        with pytest.raises(AssertionError):
            factory.set_gene(100)

    def test_except_set_transcript(self, factory):
        with pytest.raises(AssertionError):
            factory.set_transcript(100)

    def test_except_set_chromosome(self, factory):
        with pytest.raises(AssertionError):
            factory.set_chromosome(100)

    def test_except_set_txstart(self, factory):
        with pytest.raises(AssertionError):
            factory.set_transcription_start("100")

    def test_except_set_txend(self, factory):
        with pytest.raises(AssertionError):
            factory.set_transcription_end("100")

    def test_except_set_cdstart(self, factory):
        with pytest.raises(AssertionError):
            factory.set_cds_start("100")

    def test_except_set_cdsend(self, factory):
        with pytest.raises(AssertionError):
            factory.set_cds_end("100")

    def test_except_set_exoncount(self, factory):
        with pytest.raises(AssertionError):
            factory.set_exon_count("100")

    def test_except_set_exonstart(self, factory):
        with pytest.raises(AssertionError):
            factory.set_exon_starts("1, 2 ,3")

    def test_except_set_exonstarts_list(self, factory):
        with pytest.raises(AssertionError):
            factory.set_exon_starts(["1", "2", "3"])

    def test_except_set_exonends(self, factory):
        with pytest.raises(AssertionError):
            factory.set_exon_ends("2,3,4")

    def test_except_set_exonends_list(self, factory):
        with pytest.raises(AssertionError):
            factory.set_exon_ends(["2", "3", "4"])

    def test_except_set_strand(self, factory):
        with pytest.raises(AssertionError):
            factory.set_strand(1)


