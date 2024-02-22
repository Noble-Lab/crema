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
from .confidence import TdcConfidence, assign_confidence
from .writers.txt import to_txt
