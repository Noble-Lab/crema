"""A parser for the MSGF+ tab-delimited format"""
import re
import logging

import pandas as pd

from .txt import read_txt
from .. import utils

LOGGER = logging.getLogger(__name__)


def read_msgf(txt_files, pairing_file_name=None, copy_data=True):
    """Read peptide-spectrum matches (PSMs) from MSGF+ tab-delimited files.

    Parameters
    ----------
    txt_files : str, pandas.DataFrame or tuple of str
        One or more collection of PSMs in the MSGF+ tab-delimited format.
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
    peptide = "Peptide"
    spectrum = ["#SpecFile", "ScanNum"]
    pairing = ""
    protein = "Protein"
    # TODO need to test case where protein are in diff row
    protein_delim = ";"

    # Possible score columns output by MSGF+.
    scores = {
        "DeNovoSCore",
        "MSGFScore",
        "SpecEValue",
        "Evalue",
    }

    # Keep only MSGF+ scores that exist in all of the files.
    if isinstance(txt_files, pd.DataFrame):
        scores = scores.intersection(set(txt_files.columns))
    else:
        txt_files = utils.listify(txt_files)
        for txt_file in txt_files:
            with open(txt_file) as txt_ref:
                cols = txt_ref.readline().rstrip().split("\t")
                scores = scores.intersection(set(cols))

    if not scores:
        raise ValueError(
            "Could not find any of the MSGF+ score columns in all of the files."
            f"The columns crema looks for are {', '.join(list(scores))}"
        )

    scores = list(scores)

    # Read in the files:
    fields = spectrum + [peptide] + [target] + scores + [pairing] + [protein]
    if isinstance(txt_files, pd.DataFrame):
        data = txt_files.copy(deep=copy_data).loc[:, fields]
    else:
        data = pd.concat([_parse_psms(f, fields) for f in txt_files])

    data["target/decoy"] = ~data[protein].str.contains("XXX_")

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

    if pairing_file_name != None:
        psms._peptide_pairing = utils.create_pairing_from_file(
            pairing_file_name
        )

    # Remove pre/post from protein ID
    # This looks like "sp|P0AC43|SDHA_ECO57(pre=R,post=G)"
    # Remove decoy prefix from protein ID
    protein_column = psms.data[protein]
    new_protein_column = protein_column.str.replace(
        "\\([^()]*\\)", "", regex=True
    )
    new_protein_column = new_protein_column.str.replace("XXX_", "", regex=True)
    psms.set_protein_column(new_protein_column)

    return psms


def _parse_psms(txt_file, cols, log=True):
    """Parse a single MSGF+ tab-delimited file

    Parameters
    ----------
    txt_file : str
        The MSGF+ tab-delimited file to read.
    cols : list of str
        The columns to parse.

    Returns
    -------
    pandas.DataFrame
        A :py:class:`pandas.DataFrame` containing the parsed PSMs
    """
    if log:
        LOGGER.info("Reading PSMs from %s...", txt_file)
    return pd.read_csv(txt_file, sep="\t", usecols=lambda c: c in cols)
