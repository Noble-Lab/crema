"""This module contains the parser for PSMs in mzTab format"""
import logging

import pandas as pd
from pyteomics.mztab import MzTab

from ..utils import listify
from ..dataset import PsmDataset

LOGGER = logging.getLogger(__name__)


def read_mztab(mztab_files):
    """Read peptide-spectrum matches (PSMs) from mzTab files.

    Parameters
    ----------
    mztab_files : str or tuple of str
        One or more collections of PSMs in the mzTab format.

    Returns
    -------
    PsmDataset
        A :py:class:`~crema.dataset.PsmDataset` object
        containing the PSMs from the mzTab file.
    """
    mztab_files = listify(mztab_files)

    # Create a dataframe from the PSMs in the mzTab files.
    psms = pd.concat([_parse_psms(f) for f in mztab_files])

    # Initialize column names from mzTab standard specifications
    spectrum_col = ["spectra_ref"]
    score_col = [c for c in psms.columns if "search_engine_score" in c]
    target_col = "opt_global_cv_MS:1002217_decoy_peptide"
    sequence_col = "sequence"
    mod_col = "modifications"

    # Check that all column headers are valid, otherwise, throw error
    if len(set(spectrum_col) & set(psms.columns)) < len(spectrum_col):
        raise KeyError(
            "The mzTab file does not contain the columns that define a "
            f"spectrum. These are: {', '.join(spectrum_col)}."
        )

    if sequence_col not in psms.columns or mod_col not in psms.columns:
        raise KeyError(
            "The mzTab file does not the columns to specify peptide sequence "
            "and modifications."
        )

    if target_col not in psms.columns:
        raise KeyError(
            "The mzTab file does not contain the column that specifies "
            f"whether a PSM is a target or decoy, {target_col}"
        )

    if not score_col:
        raise ValueError(
            "No columns containing search engine scores were detected. These "
            "start with 'search_engine_score*'."
        )

    # Create the necesssary columns
    psms["peptide"] = psms[sequence_col] + "[" + psms[mod_col] + "]"
    psms["target"] = ~psms[target_col].astype(bool)

    # Keep only the relevant columns
    columns = spectrum_col + score_col + ["peptide", "target"]
    psms = psms.loc[:, columns]

    return PsmDataset(
        psms=psms,
        target_column="target",
        spectrum_columns=spectrum_col,
        score_columns=score_col,
        peptide_column="peptide",
        copy_data=False,
    )


def _parse_psms(mztab_file):
    """Parse a single mzTab file using Pyteomics

    Parameters
    ----------
    mztab_file : str
        The mzTab file to read.

    Returns
    -------
    pandas.DataFrame
        A :py:class:`pandas.DataFrame` containing the parsed PSMs.
    """
    LOGGER.info("Reading PSMs from %s...", mztab_file)
    return MzTab(str(mztab_file)).spectrum_match_table
