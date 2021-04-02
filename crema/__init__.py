"""Initialize the crema package."""
try:
    from importlib.metadata import version, PackageNotFoundError

    try:
        __version__ = version(__name__)
    except PackageNotFoundError:
        pass

except ImportError:
    from pkg_resources import get_distribution, DistributionNotFound

    try:
        __version__ = get_distribution(__name__).version
    except DistributionNotFound:
        pass

# Here is where we can export public functions and classes.
from crema.dataset import PsmDataset
from crema.methods import calculate_tdc
from crema.parsers import *
from crema.result import Result
from crema.params import Params
