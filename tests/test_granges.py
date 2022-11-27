"""Test the `GenomicRanges` class."""

from random import random
from pytest import raises

from pandas import DataFrame

from genomicranges import GenomicRanges

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"

row_data = DataFrame(
    {
        "seqnames": [
            "chr1",
            "chr2",
            "chr2",
            "chr2",
            "chr1",
            "chr1",
            "chr3",
            "chr3",
            "chr3",
            "chr3",
        ],
        "starts": range(100, 110),
        "ends": range(110, 120),
        "strand": ["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"],
        "score": range(0, 10),
        "GC": [random() for _ in range(10)],
    }
)

row_ranges = GenomicRanges(row_data)


def test_bad_row_data():
    """Test that we can't create a GenomicRanges object with bad row data."""
    with raises(ValueError):
        GenomicRanges(
            DataFrame(
                {
                    "starts": range(100, 110),
                    "ends": range(110, 120),
                    "strand": [
                        "-",
                        "+",
                        "+",
                        "*",
                        "*",
                        "+",
                        "+",
                        "+",
                        "-",
                        "-",
                    ],
                    "score": range(0, 10),
                    "GC": [random() for _ in range(10)],
                }
            )
        )


def test_properties():
    """Test that we can get the properties of a GenomicRanges object."""
    assert len(row_ranges) == row_ranges.len()
    assert len(row_ranges) == len(row_data)
    mcols = row_ranges.mcols()
    assert mcols is not None
    assert len(mcols) == len(row_data)


def test_methods():
    """Test the methods of a `GenomicRanges` object."""
    assert isinstance(row_ranges.granges(), GenomicRanges)


def test_nearest():
    """Test that we can get the nearest ranges."""
    nearest_ranges = GenomicRanges(
        DataFrame(
            {
                "seqnames": ["chr1", "chr2", "chr3"],
                "starts": [100, 115, 119],
                "ends": [103, 116, 120],
            }
        )
    )

    hits = row_ranges.nearest(nearest_ranges)
    assert hits is not None


def test_slices():
    """Test that we can get slices of a `GenomicRanges` object."""
    subset_genomic_ranges = row_ranges[5:8]

    assert subset_genomic_ranges is not None
    assert len(subset_genomic_ranges) == 3


def test_export():
    """Test that we can export `GenomicRanges` data to a `DataFrame`."""
    new_row_data = row_ranges.to_df()

    assert new_row_data is not None
    assert new_row_data.shape[0] == len(row_ranges)


def test_ucsc():
    ranges_data = GenomicRanges.from_ucsc("hg19")
    assert ranges_data is not None
    assert ranges_data.len() > 0
    assert len(ranges_data) == ranges_data.len()
    mcols = row_ranges.mcols()
    assert mcols is not None
    assert mcols > 0
    assert ranges_data.granges() is not None
