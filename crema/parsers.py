"""
This module contains the parsers for reading in PSMs
"""

import pandas as pd
import numpy as np
import pyteomics.mztab
from .dataset import PsmDataset


def read_file(
    input_files,
    spectrum_col=["scan"],
    score_col=["combined p-value"],
    target_col="target/decoy",
):
    """
    Read tab-delimited files.

    Parameters
    ----------
    input_files : str or tuple of str
        One or more tab-delimited file(s) to read
    spectrum_col : tuple of str, optional
        One or more column names that identify the psm. Defaults to ["scan"].
    score_col : str or tuple of str, optional
        One or more column names that identifies the scores (p-values) of the psms. Defaults to ["combined p-value"].
    target_col : str, optional
        Name of the column that indicates if a psm is a target/decoy. Defaults to "target/decoy".

    Returns
    -------
    PsmDataset
        A :py:class:`~crema.dataset.PsmDataset` object
        containing the PSM data from the given tab-delimited file.
    """
    # Store column names in a list to be used by read_csv method
    fields = []

    # Add spectrum column(s) to fields list
    if type(spectrum_col) == str:
        spectrum_col = [spectrum_col]
    for col in spectrum_col:
        fields.append(col)

    # Add score column(s) to fields list
    if type(score_col) == str:
        score_col = [score_col]
    for col in score_col:
        fields.append(col)

    # Add target column to fields list
    fields.append(target_col)

    # Create empty Pandas dataframe
    data = pd.DataFrame()

    # Loop through all given files
    if type(input_files) == str:
        input_files = [input_files]
    for file in input_files:
        data = data.append(
            pd.read_csv(file, sep=None, usecols=fields, engine="python"),
            ignore_index=True,
        )
    data = _convert_target_col(data, target_col)
    return PsmDataset(data, spectrum_col, score_col, target_col)


def read_mztab(input_file):
    """
    Read file in mzTab format.

    Parameters
    ----------
    input_file : str or tuple of str
        The mzTab file required to be read

    Returns
    -------
    PsmDataset
        A :py:class:`~crema.dataset.PsmDataset` object
        containing the PSM data from the given tab-delimited file.
    """

    # Read mzTab file using Pyteomics and extract the psm table
    psm_table = pyteomics.mztab.MzTab(input_file).spectrum_match_table

    # Initialize column names from mzTab standard specifications
    spectrum_col = ["spectra_ref"]
    score_col = []
    target_col = "opt_global_cv_MS:1002217_decoy_peptide"

    # Look through all psm_table columns for those relating to search engine score
    for col in psm_table.columns:
        if "search_engine_score" in col:
            score_col.append(col)

    # Check that all column headers are valid, otherwise, throw error
    for col in spectrum_col:
        if col not in psm_table.columns:
            raise KeyError(
                "Provided mzTab file does not contain the specified spectrum column"
            )
    for col in score_col:
        if col not in psm_table.columns:
            raise KeyError(
                "Provided mzTab file does not contain the specified score column"
            )
    if target_col not in psm_table.columns:
        raise KeyError(
            "Provided mzTab file does not contain the specified target/decoy column"
        )

    # Store column header names in a list for pandas dataframe to concat
    columns = []
    for col in spectrum_col:
        columns.append(psm_table[col])
    for col in score_col:
        columns.append(psm_table[col])
    columns.append(psm_table[target_col])

    # Create sub_table of psm containing only the necessary columns and reset index
    sub_table = pd.concat(columns, axis=1).reset_index(drop=True)
    sub_table = _convert_target_col(sub_table, target_col, decoy=True)

    return PsmDataset(
        sub_table,
        spectrum_col,
        score_col,
        target_col,
    )


def _convert_target_col(data, target_col, decoy=False):
    """
    Convert values in target column to boolean True/False.

    Parameters
    ----------
    data : pandas.DataFrame
        A pandas.DataFrame of the data before the target/decoy column has been converted to boolean
    target_col : str
        Name of the column that indicates if a psm is a target/decoy
    decoy : bool
        If true: 0 -> Target, 1 -> Decoy ; If false: 0 -> Decoy, 1 -> Target

    Returns
    -------
    data : pandas.DataFrame
        A pandas.DataFrame of the data after the target/decoy column has been converted to boolean
    """
    # Grab the first value in the target column
    identifier = data.iloc[0][target_col]
    # If the first value is already a boolean, return the data without manipulating anything
    if isinstance(identifier, bool):
        return data
    # If the first value is a numeric value, convert to boolean
    elif isinstance(identifier, (float, int, np.int64)):
        if decoy:
            targets = {
                # not (0 and 0.0 and -1 and -1.0): True,
                0 or 0.0: True,
                1 or 1.0: False,
            }
        else:
            targets = {
                not (0 and 0.0 and -1 and -1.0): True,
                0 or 0.0: False,
                -1 or -1.0: False,
            }
        data[target_col] = data[target_col].map(targets)
    # If the first value is a string, convert to boolean
    elif isinstance(identifier, str):
        targets = {
            "target": True,
            "t": True,
            "decoy": False,
            "d": False,
            "f": False,
        }
        data[target_col] = data[target_col].map(targets)
    return data
