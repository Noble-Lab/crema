"""A parser for the MSAmanda tab-delimited format"""

import re
import logging

import pandas as pd
from pathlib import Path

from .txt import read_txt
from .. import utils

LOGGER = logging.getLogger(__name__)


def read_msamanda(
    txt_files, pairing_file_name=None, decoy_prefix="REV_", copy_data=True
):
    """Read peptide-spectrum matches (PSMs) from MSAmanda tab-delimited files.

    Parameters
    ----------
    txt_files : str, pandas.DataFrame or tuple of str
        One or more collection of PSMs in the MSAmanda tab-delimited format.
    pairing_file_name : str, optional
        A tab-delimited file that explicity pairs target and decoy peptide
        sequences. Requires one column labeled 'target' that contains target
        sequences and a second column labeled 'decoy' that contains decoy
        sequences.
    decoy_prefix : str, optional
        The prefix used to indicate a decoy protein in the protein column.
        Default value is 'REV_'.
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
    scores_all = scores

    # Keep only MSAmanda scores that exist in all of the files.
    skip_first_row = False
    if isinstance(txt_files, pd.DataFrame):
        scores = scores.intersection(set(txt_files.columns))
    else:
        txt_files = utils.listify(txt_files)
        for txt_file in txt_files:
            # check file type
            if Path(txt_file).suffix != ".csv":
                raise ValueError(f"{txt_file} must be in .csv format.")

            with open(txt_file) as txt_ref:
                # First line of MSAmanda output consists only of version line
                # If statement below in case first line is removed
                line = txt_ref.readline().rstrip()
                if line.startswith("#version"):
                    line = txt_ref.readline().rstrip()
                    skip_first_row = True
                cols = line.split("\t")
                scores = scores.intersection(set(cols))

    if not scores:
        raise ValueError(
            "Could not find any of the MSAmanda score columns in all of the files."
            f"The columns Crema looks for are {', '.join(list(scores_all))}"
        )

    scores = list(scores)

    # Read in the files:
    fields = [*spectrum, peptide, target, *scores, pairing, protein]
    if isinstance(txt_files, pd.DataFrame):
        data = txt_files.copy(deep=copy_data).loc[:, fields]
    else:
        data = pd.concat(
            [
                utils.parse_psms_txt(f, fields, skip_first_row)
                for f in txt_files
            ]
        )

    data["target/decoy"] = ~data[protein].str.contains(decoy_prefix)

    psms = read_txt(
        data,
        target_column=target,
        spectrum_columns=spectrum,
        score_columns=scores,
        peptide_column=peptide,
        protein_column=protein,
        protein_delim=protein_delim,
        sep="\t",
        pairing_file_name=pairing_file_name,
        copy_data=False,
    )

    # Remove decoy prefix from protein ID
    protein_column = psms.proteins
    new_protein_column = protein_column.str.replace(
        decoy_prefix, "", regex=True
    )
    psms.set_protein_column(new_protein_column)

    return psms
