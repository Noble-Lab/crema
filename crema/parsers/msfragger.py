"""A parser for the MSFragger tab-delimited format"""
import re
import logging

import pandas as pd

from .txt import read_txt
from .pepxml import read_pepxml
from .. import utils

LOGGER = logging.getLogger(__name__)


def read_msfragger(txt_files, pairing_file_name=None, copy_data=True):
    """Read peptide-spectrum matches (PSMs) from MSFragger tab-delimited files.

    Parameters
    ----------
    txt_files : str, pandas.DataFrame or tuple of str
        One or more collection of PSMs in the MSFragger tab-delimited format.
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
    target = "label"
    peptide = "peptide"
    # TODO well annoying the column names for pepXML and tsv output are different
    spectrum = ["Filename", "start scan"]
    pairing = ""
    protein = "protein"
    protein_delim = ";"
    decoy_prefix = "rev_"

    # Possible score columns output by MSFragger.
    scores = {
        "hyperscore",
        "nextscore",
        "expect",
        "expectscore",
    }
    # scores_all = scores

    # TODO need to figure out how to for this
    # Keep only crux scores that exist in all of the files.
    # if isinstance(txt_files, pd.DataFrame):
    #    scores = scores.intersection(set(txt_files.columns))
    # else:
    #    txt_files = utils.listify(txt_files)
    #    for txt_file in txt_files:
    #        with open(txt_file) as txt_ref:
    #            cols = txt_ref.readline().rstrip().split("\t")
    #            scores = scores.intersection(set(cols))

    # if not scores:
    #    raise ValueError(
    #        "Could not find any of the MSFragger score columns in all of the files."
    #        f"The columns crema looks for are {', '.join(list(scores_all))}"
    #    )

    scores = list(scores)

    # Read in the files:
    if isinstance(txt_files, pd.DataFrame):
        data = txt_files.copy(deep=copy_data).loc[:, fields]
    else:
        data = _parse_psms(txt_files, decoy_prefix)

    if pairing_file_name != None:
        psms._peptide_pairing = utils.create_pairing_from_file(
            pairing_file_name
        )

    return data


def _parse_psms(txt_file, decoy_prefix):
    """Parse a single MSFragger pepXML file

    Parameters
    ----------
    txt_file : str
        The MSFragger pepXML file to read.
    decoy_prefix : str
        The prefix used to indicate a decoy protein in the
        description lines of the FASTA file.

    Returns
    -------
    pandas.DataFrame
        A :py:class:`pandas.DataFrame` containing the parsed PSMs
    """
    LOGGER.info("Reading PSMs from %s...", txt_file)
    return read_pepxml(txt_file, "rev_")
