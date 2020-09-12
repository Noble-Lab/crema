"""
Initialize the crema package.
"""
from pkg_resources import get_distribution, DistributionNotFound

try:
    __version__ = get_distribution("crema-ms").version
except DistributionNotFound:
    pass

# Here is where we can export public functions and classes.
