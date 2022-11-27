"""Test the `GenomicRanges` class."""

from random import random
from pytest import raises

from pandas import DataFrame

from genomicranges import GenomicRanges

__author__ = "jkanche, whargrea"
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


def test_export():
    """Test that we can export `GenomicRanges` data to a `DataFrame`."""
    new_row_data = row_ranges.to_df()
    assert new_row_data is not None
    assert len(new_row_data) == len(row_ranges)


def test_granges():
    """Test the `grab_ranges` method."""
    granges = row_ranges.granges()
    assert len(granges) == len(row_ranges)
    assert granges.metadata is None


def test_minor_methods():
    """Test the methods of a `GenomicRanges` object."""
    assert row_ranges.len() == len(row_ranges)


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
    assert len(hits) > 0


def test_properties():
    """Test properties of a `GenomicRanges` object."""
    assert row_ranges.mcols is not None
    assert row_ranges.metadata is not None
    assert row_ranges.ranges is not None
    assert row_ranges.seqnames is not None  # type: ignore
    assert row_ranges.strand is not None  # type: ignore


def test_slicing():
    """Test that we can get slices of a `GenomicRanges` object."""
    row_ranges_subset = row_ranges[5:8]

    assert row_ranges_subset is not None
    assert len(row_ranges_subset) == 3


def test_ucsc():
    ucsc_ranges = GenomicRanges.from_ucsc("hg19")
    assert ucsc_ranges.mcols is not None
    assert len(ucsc_ranges.mcols) >= 0
