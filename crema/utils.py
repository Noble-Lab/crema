"""Utility functions that are used in multiple modules"""


def listify(obj):
    """Create list containing an object if it is not already a list."""
    try:
        assert not isinstance(obj, str)
        iter(obj)
    except (AssertionError, TypeError):
        obj = [obj]

    return list(obj)


def new_column(name, df):
    """Add a new column, ensuring a unique name"""
    new_name = name
    cols = set(df.columns)
    i = 0
    while new_name in cols:
        new_name = name + "_" + str(i)
        i += 1

    return new_name
