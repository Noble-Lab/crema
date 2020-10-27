"""
This module contains the parsers for reading in PSMs
"""

import pandas as pd
from .dataset import PsmDataset


def read_file(input_files, spectrum_col, score_col, target_col, delimiter=","):
    """
    Read tab-delimited files.

    Parameters
    ----------
    input_files : tuple of str
        one or more tab-delimited file(s) to read
    spectrum_col : str
        name of the column that identifies the psm
    score_col : str
        name of the column that defines the scores (p-values) of the psms
    target_col : str
        name of the column that indicates if a psm is a target/decoy
    delimiter : str
        string character equal to what is used to separate columns
        within the tab-delimited file

    Returns
    -------
    PsmDataset
        A :py:class:`~creme_beta.dataset.PsmDataset` object
        containing the PSM data from the given tab-delimited file.
    """
    # Store column names in a list to be used by read_csv method
    fields = [spectrum_col, score_col, target_col]
    # Create empty Pandas dataframe
    data = pd.DataFrame()
    # Loop through all given files
    for file in input_files:
        data = data.append(
            pd.read_csv(file, sep=delimiter, usecols=fields), ignore_index=True
        )
    return PsmDataset(data, spectrum_col, score_col, target_col)


def read_crux(input_files):
    """
    Read crux formatted files.

    Parameters
    ----------
    input_files : tuple of str
        one or more crux formatted file(s) to read

    Returns
    -------
    PsmDataset
        A :py:class:`~creme_beta.dataset.PsmDataset` object
        containing the PSM data from the given tab-delimited file.
    """
    # Call generic read_file method with pre-determined column field names from crux format
    dataset = read_file(input_files, "scan", "combined p-value", "target/decoy", "\t")
    # Convert values in target_col to boolean True/False instead of String target/decoy
    targets = {'target': True, 'decoy': False}
    dataset.data['target/decoy'] = dataset.data['target/decoy'].map(targets)
    # Return the PsmDataset
    return dataset

