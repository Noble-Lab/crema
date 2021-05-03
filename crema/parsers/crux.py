"""A parser for the crux tab-delimited format"""
import logging

import pandas as pd
from .txt import read_txt
from ..utils import listify

LOGGER = logging.getLogger(__name__)


def read_crux(txt_files, copy_data=True):
    """Read peptide-spectrum matches (PSMs) from Crux tab-delimited files.

    Parameters
    ----------
    txt_files : str, pandas.DataFrame or tuple of str
        One or more collection of PSMs in the Crux tab-delimited format.
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
    peptide = "sequence"
    spectrum = ["scan", "spectrum precursor m/z"]

    # Possible score columns output by Crux.
    scores = {
        "sp score",
        "delta_cn",
        "delta_lcn",
        "xcorr score",
        "exact p-value",
        "refactored xcorr",
        "res-ev p-value",
        "combined p-value",
    }

    # Keep only crux scores that exist in all of the files.
    if isinstance(txt_files, pd.DataFrame):
        scores = scores.intersection(set(txt_files.columns))
    else:
        txt_files = listify(txt_files)
        for txt_file in txt_files:
            with open(txt_file) as txt_ref:
                cols = txt_ref.readline().rstrip().split("\t")
                scores = scores.intersection(set(cols))

    if not scores:
        raise ValueError(
            "Could not find any of the Crux score columns in all of the files."
            f"The columns crema looks for are {', '.join(list(scores))}"
        )

    scores = list(scores)

    # Read in the files:
    fields = spectrum + [peptide] + [target] + scores
    if isinstance(txt_files, pd.DataFrame):
        data = txt_files.copy(deep=copy_data).loc[:, fields]
    else:
        data = pd.concat([_parse_psms(f, fields) for f in txt_files])

    psms = read_txt(
        data,
        target_column=target,
        spectrum_columns=spectrum,
        score_columns=scores,
        peptide_column=peptide,
        sep="\t",
        copy_data=False,
    )

    psms.add_peptide_pairing(_create_pairing(txt_files))

    return psms


def _parse_psms(txt_file, cols):
    """Parse a single Crux tab-delimited file

    Parameters
    ----------
    txt_file : str
        The crux tab-delimited file to read.
    cols : list of str
        The columns to parse.

    Returns
    -------
    pandas.DataFrame
        A :py:class:`pandas.DataFrame` containing the parsed PSMs
    """
    LOGGER.info("Reading PSMs from %s...", txt_file)
    return pd.read_csv(txt_file, sep="\t", usecols=lambda c: c in cols)


def _create_pairing(txt_files):
    """Parse a single Crux tab-delimited file

    Parameters
    ----------
    txt_files : str, pandas.DataFrame or tuple of str
        One or more collection of PSMs in the Crux tab-delimited format.

    Returns
    -------
    pairing : dict
        A map of target and decoy peptide sequence pairings
    """
    fields = ["peptide mass", "sequence", "target/decoy", "original target sequence"]
    targets = pd.DataFrame()
    decoys = pd.DataFrame()
    if isinstance(txt_files, pd.DataFrame):
        data = txt_files.copy(deep=True).loc[:, fields]
        targets = data[data["target/decoy"] == "target"]
        decoys = data[data["target/decoy"] == "decoy"]
    else:
        for f in txt_files:
            data = _parse_psms(f, fields)
            if "original target sequence" in data.columns:
                decoys = pd.concat([decoys, data])
            else:
                targets = pd.concat([targets, data])
    pairing = {}
    decoys = decoys.drop_duplicates(["peptide mass", "sequence"]).sample(frac=1).reset_index(drop=True)
    targets = targets.drop_duplicates(["peptide mass", "sequence"]).reset_index(drop=True)
    for mass, sequence in zip(targets["peptide mass"], targets["sequence"]):
        if sequence in pairing:
            continue
        raw_sequence = _remove_mod(sequence)
        sub_decoy = decoys[decoys["original target sequence"] == raw_sequence]
        for d_index, d_row in sub_decoy.iterrows():
            if _check_match(mass, sequence, d_row):
                pairing[sequence] = d_row["sequence"]
                pairing[d_row["sequence"]] = sequence
                decoys = decoys.drop(d_index).reset_index(drop=True)
                break
    return pairing


def _remove_mod(sequence):
    """Remove the modifications of a given peptide sequence

    Parameters
    ----------
    sequence : str
        A target peptide sequence

    Returns
    -------
    sequence : str
        A peptide sequence with its modifications removed
    """
    try:
        left = sequence.index("[")
        right = sequence.index("]")
    except ValueError:
        return sequence
    return sequence[:left]+sequence[right+1:]


def _check_match(target_mass, target_sequence, decoy_row):
    """Check if a target peptide and decoy peptide make a valid pair

    Parameters
    ----------
    target_row : pandas.Series
        A row from the pandas.Dataframe holding target psm data
    decoy_row : pandas.Series
        A row from the pandas.Dataframe holding decoy psm data

    Returns
    -------
    bool
        True if the peptides make a valid pair, false otherwise. Valid peptide pairs
        are shuffled versions of one another with identical mass and identical
        amino acid modification
    """
    # check that peptides have identical mass
    if target_mass != decoy_row["peptide mass"]:
        return False
    target_mod = target_sequence.find("[")
    decoy_mod = decoy_row["sequence"].find("[")
    # if neither peptide has modifications
    if target_mod == -1 and decoy_mod == -1:
        return True
    # if one peptide has modifications and the other doesn't
    if target_mod*decoy_mod < 0:
        return False
    # if both peptides have modifications
    if target_sequence[target_mod-1] != decoy_row["sequence"][decoy_mod-1]:
        return False
    return True

