from random import random

import pandas as pd
import pytest
from genomicranges.genomic_ranges import GenomicRanges

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def test_ucsc():
    gr = GenomicRanges.fromUCSC("hg19")
    assert gr is not None
    assert gr.len() > 0
    assert len(gr) == gr.len()
    assert gr.mcols() is not None
    assert gr.mcols().shape[0] > 0
    assert gr.granges() is not None
