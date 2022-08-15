"""A parser for the MSAmanda tab-delimited format"""
import re
import logging

import pandas as pd

from .txt import read_txt
from .. import utils

LOGGER = logging.getLogger(__name__)


def read_msamanda(txt_files, pairing_file_name=None, copy_data=True):
    """Read peptide-spectrum matches (PSMs) from MSAmanda tab-delimited files.

    Parameters
    ----------
    txt_files : str, pandas.DataFrame or tuple of str
        One or more collection of PSMs in the MSAmanda tab-delimited format.
    pairing_file_name : str, optional
        A tab-delimited file that explicity pairs target and decoy peptide
        sequences. Requires one column labled 'target' that contains target
        sequences and a second colun labeled 'decoy' that contains decoy
        sequences.
    copy_data : bool, optional
        If true, a deep copy of the data is created. This uses more memory, but
        is safer because it prevents accidental modification of the underlying
        data. This argument only has an effect when `txt_files` is a
        :py:class:`pandas.DataFrame`

    Returns
    -------
    PsmDataset
        A :py:class:`~crema.dataset.PsmDataset` object containing the parsed
        PSMs.
    """
    target = "target/decoy"
    peptide = "Sequence"
    spectrum = ["Filename", "Scan Number"]
    pairing = ""
    protein = "Protein Accessions"
    protein_delim = ";"

    # Possible score columns output by MSAmanda.
    scores = {
        "Amanda Score",
        "Weighted Probability",
    }

    # Keep only crux scores that exist in all of the files.
    if isinstance(txt_files, pd.DataFrame):
        scores = scores.intersection(set(txt_files.columns))
    else:
        txt_files = utils.listify(txt_files)
        for txt_file in txt_files:
            with open(txt_file) as txt_ref:
                # First line of MSAmanda output consists only of version line
                skipLine = txt_ref.readline()
                cols = txt_ref.readline().rstrip().split("\t")
                scores = scores.intersection(set(cols))

    if not scores:
        raise ValueError(
            "Could not find any of the MSAmanda score columns in all of the files."
            f"The columns crema looks for are {', '.join(list(scores))}"
        )

    scores = list(scores)

    # Read in the files:
    fields = spectrum + [peptide] + [target] + scores + [pairing] + [protein]
    if isinstance(txt_files, pd.DataFrame):
        data = txt_files.copy(deep=copy_data).loc[:, fields]
    else:
        data = pd.concat([_parse_psms(f, fields) for f in txt_files])

    data["target/decoy"] = ~data[protein].str.contains("REV_")

    psms = read_txt(
        data,
        target_column=target,
        spectrum_columns=spectrum,
        score_columns=scores,
        peptide_column=peptide,
        protein_column=protein,
        protein_delim=protein_delim,
        sep="\t",
        copy_data=False,
    )

    # pairing with MSAmanda not possible at this time

    # Remove decoy prefix from protein ID
    protein_column = psms.data[protein]
    new_protein_column = protein_column.str.replace("REV_", "", regex=True)
    psms.set_protein_column(new_protein_column)

    return psms


def _parse_psms(txt_file, cols, log=True):
    """Parse a single MSAmanda tab-delimited file

    Parameters
    ----------
    txt_file : str
        The MSAmanda tab-delimited file to read.
    cols : list of str
        The columns to parse.

    Returns
    -------
    pandas.DataFrame
        A :py:class:`pandas.DataFrame` containing the parsed PSMs
    """
    if log:
        LOGGER.info("Reading PSMs from %s...", txt_file)
    return pd.read_csv(
        txt_file, sep="\t", skiprows=1, usecols=lambda c: c in cols
    )