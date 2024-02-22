"""A parser for the MSFragger tab-delimited format"""

import re
import logging

import pandas as pd

from .txt import read_txt
from .pepxml import read_pepxml
from .pepxml import _parse_pepxml
from ..dataset import PsmDataset
from .. import utils

LOGGER = logging.getLogger(__name__)


def read_msfragger(
    txt_files, pairing_file_name=None, decoy_prefix="rev_", copy_data=True
):
    """Read peptide-spectrum matches (PSMs) from MSFragger pepXML files.

    Parameters
    ----------
    txt_files : str, pandas.DataFrame or tuple of str
        One or more collection of PSMs in the MSFragger tab-delimited format.
    pairing_file_name : str, optional
        A tab-delimited file that explicity pairs target and decoy peptide
        sequences. Requires one column labeled 'target' that contains target
        sequences and a second column labeled 'decoy' that contains decoy
        sequences.
    decoy_prefix : str, optional
        The prefix used to indicate a decoy protein in the protein column.
        Default value is 'rev_'.
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
    # spectrum = ["Filename", "start scan"] # check this - is this for TSV?
    spectrum = ["ms_data_file", "scan"]
    pairing = ""
    # protein = "protein" # check for TSV
    protein = "proteins"
    protein_delim = ";"

    # The text below in any pepXML field identifies the field as a score field
    score_id = "search_engine_score:"

    # Possible score columns output by MSFragger.
    scores = {
        "hyperscore",
        "nextscore",
        "expect",
        "expectscore",
    }
    scores_all = scores

    # Read in the files:
    if isinstance(txt_files, pd.DataFrame):
        scores = scores.intersection(set(txt_files.columns))
    else:
        txt_files = utils.listify(txt_files)
        data_list = [_parse_pepxml(f, decoy_prefix) for f in txt_files]

        for data_file in data_list:
            score_col = [c for c in data_file.columns if score_id in c]
            score_col = [re.sub(score_id, "", c) for c in score_col]
            scores = scores.intersection(set(score_col))

        data = pd.concat(data_list)
        data.columns = data.columns.str.replace(score_id, "")

    if not scores:
        raise ValueError(
            "Could not find any of the MSFragger score columns in all of the "
            "files."
            f"The columns crema looks for are {', '.join(list(scores_all))}"
        )

    scores = list(scores)
    data[score_col] = data[score_col].astype(float)

    psms = PsmDataset(
        psms=data,
        target_column=target,
        spectrum_columns=spectrum,
        score_columns=scores,
        peptide_column=peptide,
        protein_column=protein,
        protein_delim=protein_delim,
        copy_data=False,
    )

    if pairing_file_name != None:
        psms._peptide_pairing = utils.create_pairing_from_file(
            pairing_file_name
        )

    return psms
