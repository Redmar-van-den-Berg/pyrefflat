from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *
__author__ = 'Sander Bollen'

import pytest
from pyrefflat.parser import *
from pyrefflat.models import *


class TestReader:

    def test_return(self):
        r = Reader(filename="test/data/mini.refFlat")
        re = next(r)
        assert isinstance(re, Record)

    def test_close(self):
        with pytest.raises(StopIteration):
            r = Reader(filename="test/data/mini.refFlat")
            r.close()
            _ = next(r)

class TestWriter:
    def test_write(self, writer, record):
        writer.write(record)
        writer.close()
        assert open("test/data/write.refFlat", "rb").\
                   readline().decode(ENCODING) == "MLH1\tNM_000249\t" \
                                                  "chr3\t+\t37034840\t" \
                                                  "37092337\t37035038\t" \
                                                  "37092144\t19\t" \
                                                  "37034840,37038109,37042445," \
                                                  "37045891,37048481,37050304," \
                                                  "37053310,37053501,37055922," \
                                                  "37058996,37061800,37067127," \
                                                  "37070274,37081676,37083758," \
                                                  "37089009,37090007,37090394," \
                                                  "37091976," \
                                                  "\t37035154,37038200,37042544" \
                                                  ",37045965,37048554,37050396" \
                                                  ",37053353,37053590,37056035" \
                                                  ",37059090,37061954,37067498," \
                                                  "37070423,37081785,37083822," \
                                                  "37089174,37090100,37090508," \
                                                  "37092337,\n"

    def test_close(self, writer, record):
        writer.close()
        with pytest.raises(ValueError):
            writer.write(record)


@pytest.fixture(scope="module")
def writer():
    w = Writer(filename="test/data/write.refFlat")
    return w


@pytest.fixture(scope="module")
def record():
    r = Reader(filename="test/data/mini.refFlat")
    return next(r)


@pytest.fixture(scope="module")
def reverse_record():
    r = Reader(filename="test/data/reverse.refFlat")
    return next(r)


@pytest.fixture(scope="module")
def reverse_exons(reverse_record):
    return reverse_record.exons


@pytest.fixture(scope="module")
def exons(record):
    return record.exons


@pytest.fixture(scope="module")
def proc():
    p = RefFlatProcessor("test/data/midi.refFlat")
    return p

@pytest.fixture(scope="module")
def mrecord():
    m = Reader(filename="test/data/myot.refFlat")
    return next(m)


class TestRecord():

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

    def test_fromdict(self, record):
        nrecrd = Record.fromdict(record.to_dict())
        assert isinstance(nrecrd, Record)

    def test_fromdict_notcolumn(self, record):
        d = record.to_dict()
        for c in COLUMNS:
            d.pop(c)
            with pytest.raises(ValueError):
                _ = Record.fromdict(d)

    def test_fromdict_not_nlccolumn(self, record):
        d = record.to_dict()
        for c in NUMERIC_LIST_COLUMNS:
            d[c] = "blah"
            with pytest.raises(ValueError):
                _ = Record.fromdict(d)

    def test_fromdict_not_nc_column(self, record):
        d = record.to_dict()
        for c in NUMERIC_COLUMNS:
            d[c] = "blah"
            with pytest.raises(ValueError):
                _ = Record.fromdict(d)

    def test_cds_exons(self, mrecord):
        cds_ex = mrecord.cds_exons
        assert len(cds_ex) == 9
        assert cds_ex[0].number == 2


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
            assert (ex.number == i+1)

    def test_reverse_numbers(self, reverse_exons):
        i = len(reverse_exons)
        for x in reverse_exons:
            assert(x.number == i)
            i -= 1


class TestProcessor():
    def test_process_init(self):
        p = RefFlatProcessor("test/data/midi.refFlat")
        assert isinstance(p, RefFlatProcessor)

    def test_processor_process(self, proc):
        proc.process()
        assert proc._already_processed

    def test_processor_transcripts(self, proc):
        assert len(proc.transcripts) == 30

    def test_processor_duplicates(self, proc):
        assert proc.n_duplicates == 30

    def test_process_do_not_remove(self, proc):
        proc.process(remove_duplicates=False)
        tr = []
        for x in proc.genes.values():
            tr += x.transcripts
        assert len(tr) == 60

