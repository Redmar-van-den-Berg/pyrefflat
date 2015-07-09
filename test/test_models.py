from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *
__author__ = 'Sander Bollen'

import pytest

from pyrefflat.generics import *
from pyrefflat.parser import RefFlatProcessor
from pyrefflat.models import *

@pytest.fixture(scope="module")
def proc():
    p = RefFlatProcessor("test/data/midi.refFlat")
    p.process()
    return p

class TestGene():

    def test_gene_type(self, proc):
        for gene in proc.genes.values():
            assert isinstance(gene, Gene)

    def test_gene_transcript_type(self, proc):
        for gene in proc.genes.values():
            for tr in gene.transcripts:
                assert isinstance(tr, Transcript)

    def test_gene_amount(self, proc):
        assert len(proc.genes) == 15

    def test_gene_add(self, proc):
        g = proc.genes.values()[0]
        t = g.transcripts[0]
        t.start = t.start - 1
        early_mincoord = g.min_coord
        g.update_transcripts(t)
        late_mincoord = g.min_coord
        assert late_mincoord < early_mincoord


class TestTranscript():

    def test_transcript_add(self, proc):
        proc.process()
        g = proc.genes.values()[0]
        t = g.transcripts[0]
        e = t.exons[0]
        early_len = len(t.exons)
        t.update_exons(e)
        assert early_len < len(t.exons)

    def test_transcript_add_failure(self, proc):
        proc.process()
        g = proc.genes.values()[0]
        t = g.transcripts[0]
        e = t.exons[0]
        e._start = e.start - 1
        with pytest.raises(ValueError):
            t.update_exons(e)

    def test_transcript_add_failure2(self, proc):
        proc.process()
        g = proc.genes.values()[0]
        t = g.transcripts[0]
        e = t.exons[0]
        e._end = t.end + 1
        with pytest.raises(ValueError):
            t.update_exons(e)

    def test_transcript_line(self, proc):
        proc.process()
        g = proc.genes.values()[0]
        t = g.transcripts[0]
        assert isinstance(t.line, basestring)

    def test_transcript_dict(self, proc):
        proc.process()
        g = proc.genes.values()[0]
        t = g.transcripts[0]
        assert isinstance(t.to_dict(), dict)