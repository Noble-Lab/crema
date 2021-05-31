"""A parser for generic delmited text files."""
import logging

import pandas as pd
from ..dataset import PsmDataset
from ..utils import listify

LOGGER = logging.getLogger(__name__)


def read_txt(
    txt_files,
    target_column,
    spectrum_columns,
    score_columns,
    peptide_column,
    sep="\t",
    copy_data=True,
):
    """Read peptide-spectrum matches (PSMs) from delimited text files.

    Parameters
    ----------
    txt_files : str, panda.DataFrame, or tuple of str
        One or more collection of PSMs in a tabular text format.
    target_column : str
        The column that indicates whether a PSM is a target or a decoy.
    spectrum_columns : str or tuple of str
        One or more columns that together define a unique mass spectrum.
    score_columns : str or tuple of str, optional
        One or more columns that indicate scores by which crema can rank PSMs.
    peptide_column : str
        The column that defines a unique peptide. Modifications should be
        indicated either in square brackets :code:`[]` or parentheses
        :code:`()`. The exact modification format within these entities does
        not matter, so long as it is consistent.
    sep : str, optional
        The delimiter to use.
    copy_data : bool, optional
        If true, a deep copy of the data is created. This uses more memory, but
        is safer because it prevents accidental modification of the underlying
        data. This argument only has an effect when `pin_files` is a
        :py:class:`pandas.DataFrame`

    Returns
    -------
    PsmDataset
        A :py:class:`~crema.dataset.PsmDataset` object containing the parsed
        PSMs.
    """
    # Store column names in a list to be used by read_csv function
    fields = [target_column, peptide_column]

    # Verify some arguments are lists:
    spectrum_columns = listify(spectrum_columns)
    score_columns = listify(score_columns)
    fields += spectrum_columns + score_columns

    # Parse the data
    if isinstance(txt_files, pd.DataFrame):
        data = txt_files.copy(deep=copy_data).loc[:, fields]
    else:
        data = pd.concat(
            [_parse_psms(f, sep, fields) for f in listify(txt_files)]
        )

    data[target_column] = _convert_target_col(data[target_column])
    return PsmDataset(
        psms=data,
        target_column=target_column,
        spectrum_columns=spectrum_columns,
        score_columns=score_columns,
        peptide_column=peptide_column,
        copy_data=False,
    )


def _parse_psms(txt_file, sep, cols):
    """Parse a single delimited txt file.

    Parameters
    ----------
    txt_file : str
        The crux tab-delimited file to read.
    sep : str
        The delimiter to use.
    cols : list of str
        The columns to parse.

    Returns
    -------
    pandas.DataFrame
        A :py:class:`pandas.DataFrame` containing the parsed PSMs
    """
    LOGGER.info("Reading PSMs from %s...", txt_file)
    return pd.read_csv(txt_file, sep=sep, usecols=cols)


def _convert_target_col(data):
    """Convert values in target column to boolean True/False.

    Parameters
    ----------
    data : pandas.Series
        The target column before the target/decoy column has been converted to
        boolean.

    Returns
    -------
    pandas.Series
        The target column after it has been converted to boolean.
    """
    if data.dtype == bool:
        return data
    elif data.dtype == "object":
        targets = {
            "target": True,
            "t": True,
            "decoy": False,
            "d": False,
            "f": False,
        }
        return data.map(targets)
    elif len(data.unique()) > 2:
        raise ValueError(
            "The specificed target column appears to contain more than 2 "
            "values."
        )

    return data > 0
