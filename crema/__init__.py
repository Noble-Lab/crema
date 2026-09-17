"""Initialize the crema package."""

try:
    from importlib.metadata import version, PackageNotFoundError

    try:
        __version__ = version("crema-ms")
    except PackageNotFoundError:
        pass

except ImportError:
    from pkg_resources import get_distribution, DistributionNotFound

    try:
        __version__ = get_distribution("crema-ms").version
    except DistributionNotFound:
        pass

# Here is where we can export public functions and classes.
from .dataset import PsmDataset
from .parsers.tide import read_tide
from .parsers.msgf import read_msgf
from .parsers.msamanda import read_msamanda
from .parsers.msfragger import read_msfragger
from .parsers.comet import read_comet
from .parsers.txt import read_txt
from .parsers.mztab import read_mztab
from .parsers.pepxml import read_pepxml
from .confidence import TdcConfidence, DuckdbTdcConfidence, assign_confidence
from .writers.txt import to_txt
from .writers.parquet import to_parquet
from .parsers.parquet import read_parquet

__all__ = [
    "DuckdbTdcConfidence",
    "PsmDataset",
    "TdcConfidence",
    "assign_confidence",
    "read_comet",
    "read_msamanda",
    "read_msfragger",
    "read_msgf",
    "read_mztab",
    "read_parquet",
    "read_pepxml",
    "read_tide",
    "read_txt",
    "to_parquet",
    "to_txt",
]
