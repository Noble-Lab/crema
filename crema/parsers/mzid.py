"""This module contains the parser for PSMs in mzID format"""
import logging

import pandas as pd
import pyteomics.mzid

from ..utils import listify
from ..dataset import PsmDataset

LOGGER = logging.getLogger(__name__)


def read_mzid(mzid_files):
    """Read peptide-spectrum matches (PSMs) from mzID files.

    Parameters
    ----------
    mziz_files : str or tuple of str
        One or more collections of PSMs in the mzID format.

    Returns
    -------
    PsmDataset
        A :py:class:`~crema.dataset.PsmDataset` object
        containing the PSMs from the mzTab file.
    """
    mzid_files = listify(mzid_files)

    # Create a dataframe from the PSMs in the mzID files.
    psms = pd.concat([_parse_psms(f) for f in mzid_files])
    print(psms.columns)

    # Determine which database search engine generated the mzID file.
    if "MS-GF:SpecEValue" in psms.columns:
        # Initialize column names from MSGF+ specifications
        spectrum_col = ["scan number(s)", "calculatedMassToCharge"]
        score_col = [
            c for c in psms.columns if "MS-GF:" in c and "QValue" not in c
        ]
    elif "Amanda:AmandaScore" in psms.columns:
        # Initialize column names from MSAmanda specifications
        spectrum_col = ["name", "calculatedMassToCharge"]
        score_col = [c for c in psms.columns if "Amanda:AmandaScore" in c]
    else:
        raise ValueError(
            "Unsupported database search engine "
            "generated mzid file. Please reach out to "
            "Crema team to add support."
        )

    target_col = "isDecoy"
    sequence_col = "PeptideSequence"
    mod_col = "Modification"

    # Check that all column headers are valid, otherwise, throw error
    if len(set(spectrum_col) & set(psms.columns)) < len(spectrum_col):
        raise KeyError(
            "The mzTab file does not contain the columns that define a "
            f"spectrum. These are: {', '.join(spectrum_col)}."
        )

    if sequence_col not in psms.columns or mod_col not in psms.columns:
        raise KeyError(
            "The mzID file does not the columns to specify peptide sequence "
            "and modifications."
        )

    if target_col not in psms.columns:
        raise KeyError(
            "The mzID file does not contain the column that specifies "
            f"whether a PSM is a target or decoy, {target_col}"
        )

    if not score_col:
        raise ValueError(
            "No columns containing search engine scores were detected."
        )

    # Create the necesssary columns
    psms["peptide"] = psms[sequence_col] + mod_col

    # Keep only the relevant columns
    columns = spectrum_col + score_col + ["peptide", target_col]
    psms = psms.loc[:, columns]

    return PsmDataset(
        psms=psms,
        target_column=target_col,
        spectrum_columns=spectrum_col,
        score_columns=score_col,
        peptide_column="peptide",
        copy_data=False,
    )


def _parse_psms(mzid_file):
    """Parse a single mzID file using Pyteomics

    Parameters
    ----------
    mzid_file : str
        The mzID file to read.

    Returns
    -------
    pandas.DataFrame
        A :py:class:`pandas.DataFrame` containing the parsed PSMs.
    """
    LOGGER.info("Reading PSMs from %s...", mzid_file)
    return pyteomics.mzid.DataFrame(str(mzid_file))
