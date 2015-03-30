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

