"""The `GenomicRanges` class."""

from pathlib import Path
from typing import List, MutableMapping, Optional, Tuple, Union, Iterator, Any

from ncls import NCLS, NCLS32, NCLS64  # type: ignore
from pandas import DataFrame, Series, concat  # type: ignore

from .utils import construct_ucsc_gtf_file_url, parse_gtf_file

__author__ = "jkanche, whargrea"
__copyright__ = "jkanche"
__license__ = "MIT"

# TODO: I made strand optional because one of your tests had no strand column.
required_ranges_columns = ["seqnames", "starts", "ends"]
strand = "strand"
all_ranges_columns = required_ranges_columns + [strand]
IndexType = Union[List[bool], List[int], slice]


class GenomicRangesIterator:
    """An iterator for a `GenomicRanges` object."""

    def __init__(self, genomic_ranges: "GenomicRanges") -> None:
        """Initialize a `GenomicRangesIterator` object."""
        self._genomic_ranges = genomic_ranges
        self._index = 0

    def __iter__(self) -> "GenomicRangesIterator":
        """Iterate over the intervals."""
        return self

    def __next__(self) -> Tuple[Series, Optional[Series]]:  # type: ignore
        """Return the next interval.

        Returns
        -------
        row : Series
            A `Series` object containing the next row of data.
        metadata : Series | None
            A `Series` object containing the next row of metadata, if any.
        """
        try:
            row = self._genomic_ranges.ranges.iloc[self._index]  # type: ignore
            metadata = (  # type: ignore
                None
                if self._genomic_ranges.metadata is None
                else self._genomic_ranges.metadata.iloc[self._index]
            )
        except IndexError as exc:
            raise StopIteration from exc

        self._index += 1
        return row, metadata  # type: ignore


# TODO: If re-constructing the NCLS indices is prohibitivelly expensive, we
# could make them an optional parameter to the constructor. It would mostly be
# used by methods that return a new GenomicRanges object. (This would be
# similar to how you had it before). Would require additional validation, and
# removal of the ranges setter.


class GenomicRanges:
    """The `GenomicRanges` class for representing genomic intervals."""

    def __init__(
        self,
        ranges: DataFrame,
        metadata: Optional[DataFrame] = None,
    ) -> None:
        """Initialize an instance of `GenomicRanges`.

        Parameters
        ----------
        ranges : DataFrame
            Intervals as a `DataFrame`. Required columns are `seqnames`,
            `starts`, `ends`, and `strand`. Extra columns will be considered
            metadata. Additional metadata can be provided separately using the
            `metadata` parameter.
        metadata : DataFrame | None
            Optional additional metadata in a `DataFrame`. Default = `None`.

        Raises
        ------
        ValueError
            If the `ranges` or `metadata` are not valid.
        """
        self._ranges, self._indices = self._validate_ranges(ranges)  # type: ignore # noqa: E501
        self._metadata = self._validate_metadata(
            self._combine_metadata(ranges, metadata)
        )

    @staticmethod
    def _create_indices(
        ranges: DataFrame,
    ) -> MutableMapping[str, Union[NCLS32, NCLS64]]:  # type: ignore
        """Create the indices from the `ranges` data.

        Parameters
        ----------
        ranges : DataFrame
            The `ranges` data.

        Returns
        -------
        indices : MutableMapping[str, NCLS32 | NCLS64]
            The indices.
        """
        return {
            group: NCLS(
                rows["starts"].astype(int).to_list(),  # type: ignore
                rows["ends"].astype(int).to_list(),  # type: ignore
                rows["starts"].astype(str).to_list(),  # type: ignore
            )
            for group, rows in ranges.groupby("seqnames")  # type: ignore
        }

    def _validate_ranges(
        self, ranges: DataFrame
    ) -> Tuple[DataFrame, MutableMapping[str, Union[NCLS32, NCLS64]]]:  # type: ignore # noqa: E501
        """Validate the `ranges` and create the `indices` data.

        Parameters
        ----------
        ranges : DataFrame
            Intervals as a `DataFrame`. Required columns are `seqnames`,
            `starts`, and `ends`, `strand` is optional. Extra columns will
            dropped.

        Returns
        -------
        ranges : DataFrame
            The `ranges` data.
        indices : MutableMapping[str, NCLS32 | NCLS64]
            The indices.

        Raises
        ------
        ValueError
            If the `ranges` data is not valid.
        """
        for name in required_ranges_columns:
            if name not in ranges.columns:
                raise ValueError(f"Missing column '{name}' in `ranges`.")

        ranges = (
            ranges[all_ranges_columns]
            if strand in ranges.columns
            else ranges[required_ranges_columns]
        )

        return ranges, self._create_indices(ranges)  # type: ignore

    def _validate_metadata(
        self, metadata: Optional[DataFrame]
    ) -> Optional[DataFrame]:
        """Validate the `metadata`.

        Parameters
        ----------
        metadata : DataFrame | None
            The metadata to validate.

        Returns
        -------
        metadata : DataFrame | None
            The validated `metadata` data.

        Raises
        ------
        ValueError
            If the `metadata` data is not valid.
        """
        if metadata is not None and len(metadata) != len(self._ranges):
            raise ValueError(
                "Metadata must have the same number of rows as `ranges`."
            )

        return metadata

    @staticmethod
    def _combine_metadata(
        ranges: DataFrame, metadata: Optional[DataFrame]
    ) -> Optional[DataFrame]:
        """Combine `ranges` and `metadata` metadata into a single `DataFrame`.

        Parameters
        ----------
        ranges : DataFrame
            The `ranges` data.
        metadata : DataFrame | None
            The metadata to set.

        Returns
        -------
        metadata : DataFrame | None
            The `metadata` data.

        Raises
        ------
        ValueError
            If the metadata is not valid.
        """
        ranges_metadata = ranges.drop(
            columns=all_ranges_columns
            if strand in ranges.columns
            else required_ranges_columns,
            axis=1,
        )
        if len(ranges_metadata.columns) > 0:
            if metadata is not None:
                metadata = concat([ranges_metadata, metadata], axis=1)
            else:
                metadata = ranges_metadata

        return metadata

    @property
    def ranges(self) -> DataFrame:
        """Get or set the `ranges` of the `GenomicRanges` object.

        Parameters
        ----------
        ranges : DataFrame
            A replacement ranges, any metadata columns will be dropped. To
            change the metadata, use the `metadata` property.

        Returns
        -------
        index : DataFrame
            The current ranges.

        Raises
        ------
        ValueError
            If the `ranges` data is not valid.
        """
        return self._ranges

    # TODO: accept a different number of ranges?
    @ranges.setter
    def ranges(self, ranges: DataFrame) -> None:
        if len(ranges) != len(self._ranges):
            raise ValueError(
                "New `ranges` must have the same number of rows as old "
                "`ranges`."
            )

        self._ranges, self._indices = self._validate_ranges(ranges)  # type: ignore # noqa: E501

    @property
    def metadata(self) -> Optional[DataFrame]:
        """Get or set the `metadata` of the `GenomicRanges` object.

        Parameters
        ----------
        metadata : Optional[DataFrame]
            A replacement metadata.

        Returns
        -------
        metadata : Optional[DataFrame]
            The current metadata.

        Raises
        ------
        ValueError
            If the `metadata` data is not valid.
        """
        return self._metadata

    @metadata.setter
    def metadata(self, metadata: Optional[DataFrame]) -> None:
        self._metadata = self._validate_metadata(metadata)

    # TODO: rename metadata property to mcols?
    @property
    def mcols(self) -> Optional[DataFrame]:
        """Get the metadata.

        Returns
        -------
        metadata : DataFrame
            Metadata across all intervals.
        """
        return self.metadata

    @property
    def seqnames(self) -> Series:  # type: ignore
        """Get the chromosome or sequence names.

        Returns
        -------
        seqnames : Series[str]
            All chromosome or sequence names.
        """
        return self._ranges["seqnames"]  # type: ignore

    @property
    def strand(self) -> Optional[Series]:  # type: ignore
        """Get the strand data.

        Returns
        -------
        seqnames : Series[str] | None
            Strand data by chromosome or sequence name, if any.
        """
        if strand in self._ranges.columns:
            return self._ranges[strand]  # type: ignore

        return None

    def granges(self) -> "GenomicRanges":
        """Get a new `GenomicRanges` object without the `metadata`.

        Returns
        -------
        genomic_ranges : GenomicRanges
            `GenomicRanges` without metadata.
        """
        return GenomicRanges(self._ranges)

    # TODO: is needed?
    def len(self) -> int:
        """Get the number of intervals.

        Returns
        -------
        int_intervals : int
            The number of intervals.
        """
        return len(self._ranges)

    def __len__(self) -> int:
        """Get the number of intervals.

        Returns
        -------
        int_intervals : int
            The number of intervals.
        """
        return self.len()

    # TODO: is needed?
    def length(self) -> int:
        """Alias for `len()`.

        Returns
        -------
        int_intervals : int
            The number of intervals.
        """
        return self.len()

    def __iter__(self) -> GenomicRangesIterator:
        """Iterate over the intervals.

        Returns
        -------
        genomic_ranges_iterator : GenomicRangesIterator
            An iterator over the intervals.
        """
        return GenomicRangesIterator(self)

    # TODO: revise this
    def __str__(self) -> str:
        """Represent the object as a `str`.

        Returns
        -------
        __str__ : str
            A description of the object.
        """
        return f"""\
Class: 'GenomicRanges'
    Num intervals: {self.ranges.shape[0]}
"""

    def __repr__(self) -> str:
        """Represent the object as a `str`.

        Returns
        -------
        __repr__ : str
            A description of the object.
        """
        return self.__str__()

    def __getitem__(self, indices: IndexType) -> "GenomicRanges":
        """Slice the `GenomicRanges` object.

        Parameters
        ----------
        indices : IndexType
            The indices to slice by.

        Returns
        -------
        genomic_ranges : GenomicRanges
            A new `GenomicRanges` object containg only the specified indices.

        Raises
        ------
        IndexError
            Indices exceeds the dataset dimensions
        """
        return GenomicRanges(
            self._ranges.iloc[indices, :],
            None
            if self._metadata is None
            else self._metadata.iloc[indices, :],
        )

    def nearest(
        self, x: "GenomicRanges", k: int = 1
    ) -> List[Optional[List[int]]]:
        """Find the nearest overlapping interval(s) to each interval in `x`.

        Parameters
        ----------
        x : GenomicRanges
            The `GenomicRanges` object to find the nearest intervals for.
        k : int
            The number of nearest neighbors to return. Default = `1`.

        Returns
        -------
        nearest : List[List[int] | None]
            The nearest interval(s) for each interval in `x`.
        """
        hits: List[Optional[List[int]]] = []
        for row, _ in x:  # type: ignore
            # TODO: figure out actual iterator return type
            row_hits: Iterator[Tuple[Any, Any, int]] = self._indices[  # type: ignore # noqa: E501
                row["seqnames"]
            ].find_overlap(  # type: ignore
                row["starts"], row["ends"]
            )

            sub_hits: List[int] = []
            for hit in row_hits:
                sub_hits.append(hit[2])
                if len(sub_hits) >= k:
                    break

            hits.append(sub_hits if len(sub_hits) > 0 else None)

        return hits

    def to_df(self) -> DataFrame:
        """Convert `GenomicRanges` to a `DataFrame`.

        Returns
        -------
        df : DataFrame
            A `DataFrame` representation of the `GenomicRanges` object.
        """
        df = self.ranges
        if self.metadata is not None:
            df = concat([df, self.metadata], ignore_index=False, axis=1)

        return df

    @staticmethod
    def from_gtf(gtf: Union[Path, DataFrame, str]) -> "GenomicRanges":
        """Read a `GTF` file into a `GenomicRanges` object.

        Parameters
        ----------
        gtf_file : Path | DataFrame | str
            A `GTF` file or a `DataFrame` representation of a `GTF` file made
            by the `parse_gtf_file` function.

        Returns
        -------
        genomic_ranges : GenomicRanges
            A `GenomicRanges` object representing the `GTF` file.
        """
        if isinstance(gtf, str):
            gtf = Path(gtf)

        if isinstance(gtf, Path):
            gtf = parse_gtf_file(gtf)

        return GenomicRanges(gtf)

    @staticmethod
    def from_ucsc(
        genome: str, genome_type: str = "refGene"
    ) -> "GenomicRanges":
        """Load a genome `GTF` file from UCSC as a `GenomicRanges` object.

        Parameters
        ----------
        genome : str
            The genome shortcode: `"hg19"`, `"hg38"`, `"mm10"`, etc.
        genome_type : str
            The type of genome to load, one of: `"refGene"`, `"ensGene"`,
            `"knownGene"` or `"ncbiRefSeq"`. Default = `"refGene"`.
        """
        return GenomicRanges(
            parse_gtf_file(
                construct_ucsc_gtf_file_url(genome, genome_type=genome_type)
            )
        )
