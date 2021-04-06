"""Utility functions that are used in multiple modules"""


def listify(obj):
    """Create list containing an object if it is not already a list."""
    try:
        assert not isinstance(obj, str)
        iter(obj)
    except (AssertionError, TypeError):
        obj = [obj]

    return list(obj)
