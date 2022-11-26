"""Configuration for summarized experiment types."""

from typing import List, Optional, Tuple, Union

IndexType = Union[List[int], slice]
OptionalIndexType = Optional[IndexType]
OptionalIndicesType = Union[List[IndexType], Tuple[IndexType, ...]]
