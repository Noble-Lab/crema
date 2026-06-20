"""Parser for Parquet input files."""

import logging

import pandas as pd

from ..dataset import PsmDataset
from .. import utils

LOGGER = logging.getLogger(__name__)


def read_parquet(
    parquet_files,
    target_column,
    spectrum_columns,
    score_columns,
    peptide_column,
    protein_column,
    protein_delim,
    pairing_file_name=None,
    copy_data=False,
):
    """Read peptide-spectrum matches (PSMs) from Parquet files.

    Requires the ``pyarrow`` package (``pip install crema[fast]``).

    Parameters
    ----------
    parquet_files : str or list of str
        One or more Parquet files to read.
    target_column : str
        Column indicating target (True) or decoy (False) PSMs.
    spectrum_columns : str or list of str
        Column(s) that together uniquely identify a spectrum.
    score_columns : str or list of str
        Column(s) containing PSM scores.
    peptide_column : str
        Column containing peptide sequences.
    protein_column : str
        Column containing protein identifiers.
    protein_delim : str
        Delimiter separating multiple protein IDs in one cell.
    pairing_file_name : str, optional
        Tab-delimited file pairing target and decoy peptide sequences.
    copy_data : bool, optional
        If True, deep-copy the data on construction.

    Returns
    -------
    PsmDataset
    """
    try:
        import pyarrow.parquet as pq  # noqa: F401
    except ImportError as exc:
        raise ImportError(
            "The 'pyarrow' package is required to read Parquet files. "
            "Install it with: pip install crema[fast]"
        ) from exc

    spectrum_columns = utils.listify(spectrum_columns)
    score_columns = utils.listify(score_columns)
    fields = [
        target_column,
        peptide_column,
        protein_column,
        *spectrum_columns,
        *score_columns,
    ]

    frames = []
    for f in utils.listify(parquet_files):
        LOGGER.info("Reading PSMs from %s...", f)
        df = pd.read_parquet(f, columns=fields)
        frames.append(df)

    data = pd.concat(frames, ignore_index=True)

    # Reuse the same target-column conversion as the txt parser
    from ..parsers.txt import _convert_target_col

    data[target_column] = _convert_target_col(data[target_column])

    psms = PsmDataset(
        psms=data,
        target_column=target_column,
        spectrum_columns=spectrum_columns,
        score_columns=score_columns,
        peptide_column=peptide_column,
        protein_column=protein_column,
        protein_delim=protein_delim,
        copy_data=False,
    )

    if pairing_file_name is not None:
        psms._peptide_pairing = utils.create_pairing_from_file(
            pairing_file_name
        )

    return psms
