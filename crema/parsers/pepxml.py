"""
This module contains the parser for PSMs in pepXML format.
This code is heavily based on Will Fondrie's Mokapot pepxml parser code
"""

import logging
from lxml import etree
from functools import partial

import pandas as pd
import itertools
import re

from ..utils import listify
from ..dataset import PsmDataset

LOGGER = logging.getLogger(__name__)


def read_pepxml(pepxml_files, decoy_prefix):
    """Read peptide-spectrum matches (PSMs) from pepXML files.

    Parameters
    ----------
    pepxml_files : str or tuple of str
       One or more collections of PSMs in the pepXML format.
    decoy_prefix : str
       The prefix used to indicate a decoy protein in the
       description lines of the FASTA file.

    Returns
    -------
    PsmDataset
        A :py:class:`~crema.dataset.PsmDataset` object
        containing the PSMs from the pepxml file.
    """
    pepxml_files = listify(pepxml_files)

    # Create a dataframe from the PSMs in the pepXML files.
    psms = pd.concat([_parse_pepxml(f, decoy_prefix) for f in pepxml_files])

    # Initialize column names from pepXML standard specifications
    spectrum_col = ["ms_data_file", "scan"]
    score_col = [c for c in psms.columns if "search_engine_score" in c]
    target_col = "label"
    sequence_col = "peptide"
    protein_col = "proteins"
    protein_delim = ","

    # Check that all column headers are valid, otherwise, throw error
    if len(set(spectrum_col) & set(psms.columns)) < len(spectrum_col):
        raise KeyError(
            "The pepXML file does not contain the columns that define a "
            f"spectrum. These are: {', '.join(spectrum_col)}."
        )

    if sequence_col not in psms.columns:
        raise KeyError(
            "The pepXML file does not contain the columns to specify "
            "peptide sequence."
        )

    if target_col not in psms.columns:
        raise KeyError(
            "The pepXML file does not contain the column that specifies "
            f"whether a PSM is a target or decoy, {target_col}"
        )

    if not score_col:
        raise ValueError(
            "No columns containing search engine scores were detected. These "
            "start with 'search_engine_score*'."
        )

    if not protein_col:
        raise ValueError(
            "The pepXML file does not contain the columns to specify "
            "protein sequence."
        )

    # Remove "search_engine_score:" from column name and score_col
    # and convert scores to float
    psms.columns = psms.columns.str.replace("search_engine_score:", "")
    score_col = [re.sub("search_engine_score:", "", c) for c in score_col]
    psms[score_col] = psms[score_col].astype(float)

    return PsmDataset(
        psms=psms,
        target_column=target_col,
        spectrum_columns=spectrum_col,
        score_columns=score_col,
        peptide_column=sequence_col,
        protein_column=protein_col,
        protein_delim=protein_delim,
        copy_data=False,
    )


def _parse_pepxml(pepxml_file, decoy_prefix):
    """Parse a single pepXML file using lxml into a DataFrame

    Parameters
    ----------
    pepxml_file : str
        The pepXML file to read.
    decoy_prefix : str
        The prefix used to indicate a decoy protein in the
        description lines of the FASTA file.

    Returns
    -------
    pandas.DataFrame
        A :py:class:`pandas.DataFrame` containing the parsed PSMs.
    """
    LOGGER.info("Reading PSMs from %s...", pepxml_file)
    parser = etree.iterparse(str(pepxml_file), tag="{*}msms_run_summary")
    parse_fun = partial(_parse_msms_run, decoy_prefix=decoy_prefix)
    spectra = map(parse_fun, parser)
    try:
        psms = itertools.chain.from_iterable(spectra)
        df = pd.DataFrame.from_records(itertools.chain.from_iterable(psms))
        df["ms_data_file"] = df["ms_data_file"].astype("category")
    except etree.XMLSyntaxError:
        raise ValueError(
            f"{pepxml_file} is not a PepXML file or is malformed."
        )
    return df


def _parse_msms_run(msms_run, decoy_prefix):
    """Parse a single MS/MS run.

    Each of these corresponds to a raw MS data file.

    Parameters
    ----------
    msms_run: tuple of anything, lxml.etree.Element
        The second element of the tuple should be the XML element for a single
        msms_run. The first is not used, but is necessary for compatibility
        with using :code:`map()`.
    decoy_prefix : str
        The prefix used to indicate a decoy protein in the description lines of
        the FASTA file.
    Yields
    ------
    dict
        A dictionary describing all of the PSMs in a run.
    """
    msms_run = msms_run[1]
    ms_data_file = msms_run.get("base_name")
    run_ext = msms_run.get("raw_data")
    if not ms_data_file.endswith(run_ext):
        ms_data_file += run_ext

    run_info = {"ms_data_file": ms_data_file}
    for spectrum in msms_run.iter("{*}spectrum_query"):
        yield _parse_spectrum(spectrum, run_info, decoy_prefix)


def _parse_spectrum(spectrum, run_info, decoy_prefix):
    """Parse the PSMs for a single mass spectrum
    Parameters
    ----------
    spectrum : lxml.etree.Element
        The XML element for a single
    run_info : dict
        The parsed run data.
    decoy_prefix : str
        The prefix used to indicate a decoy protein in the description lines of
        the FASTA file.
    Yields
    ------
    dict
        A dictionary describing all of the PSMs for a spectrum.
    """
    spec_info = run_info.copy()
    spec_info["scan"] = int(spectrum.get("end_scan"))
    for psms in spectrum.iter("{*}search_result"):
        for psm in psms.iter("{*}search_hit"):
            yield _parse_psm(psm, spec_info, decoy_prefix=decoy_prefix)


def _parse_psm(psm_info, spec_info, decoy_prefix):
    """Parse a single PSM
    Parameters
    ----------
    psm_info : lxml.etree.Element
        The XML element containing information about the PSM.
    spec_info : dict
        The parsed spectrum data.
    decoy_prefix : str
        The prefix used to indicate a decoy protein in the description lines of
        the FASTA file.
    Returns
    -------
    dict
        A dictionary containing parsed data about the PSM.
    """
    psm = spec_info.copy()
    psm["peptide"] = psm_info.get("peptide")
    psm["proteins"] = [psm_info.get("protein").split(" ")[0]]
    psm["label"] = not psm["proteins"][0].startswith(decoy_prefix)

    queries = [
        "{*}modification_info",
        "{*}search_score",
        "{*}alternative_protein",
    ]
    for element in psm_info.iter(*queries):
        if "modification_info" in element.tag:
            offset = 0
            mod_pep = psm["peptide"]
            for mod in element.iter("{*}mod_aminoacid_mass"):
                idx = offset + int(mod.get("position"))
                mass = mod.get("mass")
                mod_pep = mod_pep[:idx] + "[" + mass + "]" + mod_pep[idx:]
                offset += 2 + len(mass)

            psm["peptide"] = mod_pep

        elif "alternative_protein" in element.tag:
            psm["proteins"].append(element.get("protein").split(" ")[0])
            if not psm["label"]:
                psm["label"] = not psm["proteins"][-1].startswith(decoy_prefix)

        else:
            psm["search_engine_score:" + element.get("name")] = element.get(
                "value"
            )

    psm["proteins"] = ",".join(psm["proteins"])
    return psm
