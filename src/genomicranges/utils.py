"""gtf_file Parser."""

from logging import info
from pathlib import Path
from typing import Any, Dict, Optional, Union

from joblib import Parallel, delayed  # type: ignore
from pandas import DataFrame, Series, read_csv  # type: ignore

# Variation of
# https://github.com/epiviz/epivizfileserver/src/epivizfileserver/cli.py

__author__ = "jkanche, whargrea"
__copyright__ = "jkanche"
__license__ = "MIT"

columns = [
    "seqnames",
    "source",
    "feature",
    "starts",
    "ends",
    "score",
    "strand",
    "frame",
    "group",
]


def _parse_all_attributes(row: Series) -> Dict[str, Any]:  # type: ignore
    """Parse all attributes in a `GTF` file row.

    Parameters
    ----------
    row : Series
        A `Series` object of a `GTF` file row.

    Returns
    -------
    attributes : Dict[str, Any]
        A `Dict` of all attributes of a `GTF` file row.
    """
    row: Dict[str, Any] = row.to_dict()  # type: ignore
    attr: str
    for attr in row["group"].split(";"):
        if len(attr) == 0:
            continue
        attr_name, attr_value = attr.strip().split(" ", 1)
        row[attr_name] = attr_value.strip().strip('"')

    del row["group"]

    return row


def parse_gtf_file(path: Union[Path, str]) -> DataFrame:
    """Parse a `GTF` file into a `DataFrame`.

    Parameters
    ----------
    path : Union[Path, str]
        The path to the gtf_file.
    compressed : bool, optional
        Whether the gtf_file is compressed,. Default = `True`.

    Returns
    -------
    data_frame : DataFrame
        A `DataFrame` of the gtf_file.
    """
    info(f"Reading File: '{path}'")

    df = read_csv(path, sep="\t", names=columns)

    rows: Optional[Dict[str, Any]] = Parallel(n_jobs=-2)(
        delayed(_parse_all_attributes)(row) for _, row in df.iterrows()  # type: ignore # noqa: E501
    )

    if rows is None:
        raise ValueError("Could not parse gtf_file.")

    return DataFrame.from_records(rows)  # type: ignore


def construct_ucsc_gtf_file_url(
    genome: str, genome_type: str = "refGene"
) -> str:
    """Construct a UCSC `GTF` file URL.

    Parameters
    ----------
    genome : str
        The genome shortcode: `"hg19"`, `"hg38"`, `"mm10"`, etc.
    genome_type : str
        The type of genome to load, one of: `"refGene"`, `"ensGene"`,
        `"knownGene"` or `"ncbiRefSeq"`. Default = `"refGene"`.

    Raises:
        Exception: TypeError, when `type` does not match with a valid input
    """
    if genome_type not in ["refGene", "ensGene", "knownGene", "ncbiRefSeq"]:
        raise TypeError(
            f"Provided genome_type: {genome_type} which is not one of "
            "'refGene', 'ensGene', 'knownGene' or 'ncbiRefSeq'."
        )

    base_path = (
        f"http://hgdownload.cse.ucsc.edu/goldenPath/{genome}/bigZips/genes/"
    )

    full_path = f"{base_path}/{genome}.{genome_type}.gtf.gz"

    return full_path
