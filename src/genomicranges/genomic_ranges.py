"""The `GenomicRanges` class."""

import logging
from ast import Num, Slice
from typing import MutableMapping, List, Union, Optional, Tuple

from ncls import NCLS32, NCLS64  # type: ignore
from pandas import DataFrame, Series, concat

from .gtf import parse_gtf
from .ucsc import access_gtf_ucsc
from .utils import split_pandas_df
from .config import IndexType

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


class GenomicRanges:
    """The `GenomicRanges` class for representing genomic intervals."""

    def __init__(
        self,
        index: MutableMapping[str, Union[NCLS32, NCLS64]],  # type: ignore
        indices: List[Num],
        ranges: DataFrame,
        metadata: Optional[DataFrame] = None,
    ) -> None:
        """Initialize an instance of `GenomicRanges`.

        Parameters
        ----------
        index : MutableMapping[str, Union[NCLS32, NCLS64]]
            Nested Containment List (NCL) index.
        indices : List[Num]
            Indices that map between `index`, `ranges` and `metadata`.
        ranges : DataFrame
            Intervals as a `DataFrame`.
        metadata : DataFrame | None
            Optional metadata in a `DataFrame`. Default = `None`.
        """
        self._index = index  # type: ignore
        self._indices = indices
        self._ranges = ranges
        self._metadata = metadata
        self._iter_index = 0

        self._validate_data()

    def _validate_ranges(self, ranges: DataFrame) -> None:
        """Validate the `ranges` data."""
        # TODO

    def _validate_index(
        self,
        index: MutableMapping[str, Union[NCLS32, NCLS64]],  # type: ignore
    ) -> None:
        """Validate the `index` data."""
        # TODO

    def _validate_indices(self, indices: List[Num]) -> None:
        """Validate the `indices` data."""
        # TODO

    def _validate_metadata(self, metadata: Optional[DataFrame]) -> None:
        """Validate the `metadata` data."""
        # TODO

    def _validate_data(self) -> None:
        """Validate the coherency of data fields."""
        self._validate_ranges(self._ranges)
        self._validate_index(self._index)  # type: ignore
        self._validate_indices(self._indices)
        self._validate_metadata(self._metadata)

    @property
    def ranges(self) -> DataFrame:
        """Get or set the `ranges` of the `GenomicRanges` object.

        Parameters
        ----------
        ranges : DataFrame
            A replacement ranges.

        Returns
        -------
        index : DataFrame
            The current ranges.
        """
        return self._ranges

    @ranges.setter
    def ranges(self, ranges: DataFrame) -> None:
        self._validate_ranges(ranges)
        self._ranges = ranges

    @property
    def index(self) -> MutableMapping[str, Union[NCLS32, NCLS64]]:  # type: ignore # noqa: E501
        """Get or set the `index` of the `GenomicRanges` object.

        Parameters
        ----------
        index : MutableMapping[str, Union[NCLS32, NCLS64]]
            A replacement index.

        Returns
        -------
        index : MutableMapping[str, Union[NCLS32, NCLS64]]
            The current index.
        """
        return self._index  # type: ignore

    @index.setter
    def index(
        self,
        index: MutableMapping[str, Union[NCLS32, NCLS64]],  # type: ignore
    ) -> None:
        self._validate_index(index)  # type: ignore
        self._index = index  # type: ignore

    @property
    def indices(self) -> List[Num]:
        """Get or set the `indices` of the `GenomicRanges` object.

        Parameters
        ----------
        indices : List[Num]
            A replacement indices.

        Returns
        -------
        index : List[Num]
            The current indices.
        """
        return self._indices

    @indices.setter
    def indices(self, indices: List[Num]) -> None:
        self._validate_indices(indices)
        self._indices = indices

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
        """
        return self._metadata

    @metadata.setter
    def metadata(self, metadata: Optional[DataFrame]) -> None:
        self._validate_metadata(metadata)
        self._metadata = metadata

    def seqnames(self) -> Series:  # type: ignore
        """Get the chromosome or sequence names.

        Returns
        -------
        seqnames : Series[str]
            All chromosome or sequence names.
        """
        return self._ranges["seqnames"]  # type: ignore

    def strand(self) -> Series:  # type: ignore
        """Get the strand data.

        Returns
        -------
        seqnames : Series[str]
            Strand data by chromosome or sequence name.
        """
        return self._ranges["strand"]  # type: ignore

    def granges(self) -> "GenomicRanges":
        """Get a new `GenomicRanges` object without metadata.

        Returns
        -------
        genomic_ranges : GenomicRanges
            `GenomicRanges` without metadata.
        """
        return GenomicRanges(self._index, self._indices, self._ranges)  # type: ignore # noqa: E501

    def mcols(self) -> Optional[DataFrame]:
        """Get the metadata.

        Returns
        -------
        metadata : DataFrame
            Metadata across all intervals.
        """
        return self._metadata

    def len(self) -> int:
        """Get the number of intervals.

        Returns
        -------
        num_intervals : int
            The number of intervals.
        """
        return self.ranges.shape[0]

    def __len__(self) -> int:
        """Get the number of intervals.

        Returns
        -------
        num_intervals : int
            The number of intervals.
        """
        return self.len()

    def length(self) -> int:
        """Alias for `len()`.

        Returns
        -------
        num_intervals : int
            The number of intervals.
        """
        return self.len()

    # no __next__ because we want genomic ranges to be an iterable (repeatable)
    # and not an iterator (not repeatable), this way we can iterate multiple
    # times
    def __iter__(self) -> Tuple[Series, Optional[Series]]:  # type: ignore
        """Iterate over the intervals.

        Returns
        -------
        row : Series
            A `Series` object containing the next row of data.
        metadata : Series | None
            A `Series` object containing the paired row of metadata if any.
        """
        for name, row in self._ranges.iterrows():  # type: ignore
            yield row, None if self._metadata is None else self._metadata.loc[
                name, :
            ]

    def __str__(self) -> str:
        """Represent the object as a `str`.

        Returns
        -------
        __str__ : str
            A description of the object.
        """
        num_seqs = len(self._index.keys())  # type: ignore
        return f"""\
Class: 'GenomicRanges'
    Num sequences: {num_seqs}
    Num intervals: {self.ranges.shape[0]}
"""

    # HERE, also need to fill validators
    def __getitem__(self, indices: IndexType) -> "GenomicRanges":
        """Slice the object

        Args:
            args (Union[Slice, List[int]]): A slice object or a list of indices to subset

        Raises:
            Exception: args must be either a slice or list of indices less than the shape
            IndexError: indices exceeds the dataset dimensions

        Returns:
            GenomicRanges: a subset GenomicRanges object
        """
        indices_subset = (
            self._indices[indices]
            if isinstance(indices, slice)
            else [self._indices[index] for index in indices]
        )

        df = self.ranges.iloc[indices_subset, :]
        if self.metadata is not None:
            df = concat(
                [
                    df,
                    self.metadata.iloc[
                        indices_subset,
                    ],
                ],
                ignore_index=False,
                axis=1,
            )

        (indexes, indices, ranges, metadata) = split_pandas_df(df)
        return GenomicRanges(indexes, indices, ranges, metadata)

    def nearest(
        self, x: "GenomicRanges", k: int = 1
    ) -> List[Union[int, List[int]]]:
        """Find nearest neighbors that overlap with the positions in `x`

        Args:
            x (GenomicRanges): input genomic positions to find
            k (int, optional): find k nearest neighbors ?. Defaults to 1.

        Returns:
            List[Union[int, List[int]]]: Possible hits
        """
        hits = []
        for gr in x:
            grhits = self.index[gr[0]["seqnames"]].find_overlap(
                gr[0]["starts"], gr[0]["ends"]
            )
            counter = 0
            tmp_hits = []
            for i in grhits:
                if counter < k:
                    tmp_hits.append(i[2])
                else:
                    break
            hits.append(tmp_hits if len(tmp_hits) > 0 else None)
        return hits

    def toDF(self) -> DataFrame:
        """Export or Convert `GenomicRanges` to Pandas DataFrame

        Returns:
            pd.DataFrame: Pandas DataFrame
        """
        df = self.ranges
        if self.metadata is not None:
            df = concat([df, self.metadata], ignore_index=False, axis=1)

        return df

    @staticmethod
    def fromPandas(data: DataFrame) -> "GenomicRanges":
        """Convert a pandas Dataframe to GenomicRanges

        Args:
            data (pd.DataFrame): a Pandas DataFrame object containing genomic positions.
                Must contain `seqnames`, `starts` & `ends` columns.

        Raises:
            Exception: Validation Error

        Returns:
            GenomicRanges: An Object representing genomic positions
        """

        # validation:
        # if `seqnames`, `starts` and `ends` don't exist, abort!
        if not set(["seqnames", "starts", "ends"]).issubset(
            set(data.columns.tolist())
        ):
            logging.error(
                f"DataFrame does not contain columns: `seqnames`, `starts` and `ends`"
            )
            raise Exception(
                f"DataFrame does not contain columns: `seqnames`, `starts` and `ends`"
            )

        (indexes, indices, ranges, metadata) = split_pandas_df(data)
        return GenomicRanges(indexes, indices, ranges, metadata)

    @staticmethod
    def fromGTF(file: str) -> "GenomicRanges":
        """Load a GTF file as GenomicRanges

        Args:
            file (str): path to gtf file

        Returns:
            GenomicRanges:  An Object representing genomic positions
        """
        compressed = True if file.endswith("gz") else False
        data = parse_gtf(file, compressed=compressed)

        (indexes, indices, ranges, metadata) = split_pandas_df(data)

        return GenomicRanges(indexes, indices, ranges, metadata)

    @staticmethod
    def fromUCSC(genome: str, type: str = "refGene") -> "GenomicRanges":
        """Load a GTF file from UCSC as GenomicRanges

        Args:
            genome (str): genome shortcode; e.g. hg19, hg38, mm10 etc
            type (str): One of refGene, ensGene, knownGene or ncbiRefSeq

        Returns:
            GenomicRanges:  An Object representing genomic positions
        """
        path = access_gtf_ucsc(genome, type=type)
        compressed = True
        data = parse_gtf(path, compressed=compressed)

        (indexes, indices, ranges, metadata) = split_pandas_df(data)

        return GenomicRanges(indexes, indices, ranges, metadata)
