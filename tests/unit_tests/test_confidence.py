"""
These tests verify the confidence implementations function as expected
"""

import pytest
import numpy as np

from crema.confidence import TdcConfidence, MixmaxConfidence
from crema.dataset import PsmDataset

from .test_dataset import simple_df


@pytest.fixture
def simple_psms(simple_df):
    return PsmDataset(
        psms=simple_df,
        target_column="target",
        spectrum_columns=["scan", "spectrum precursor m/z"],
        score_columns=["combined p-value", "x"],
        peptide_column="sequence",
        protein_column="protein id",
        protein_delim=",",
    )


def test_tdc_confidence(simple_psms: PsmDataset):
    """Test that we can compute confidence with `TdcConfidence`"""
    conf = simple_psms.assign_confidence(
        score_column="x",
        method="tdc",
        pep_fdr_type="psm-only",
    )

    assert isinstance(conf, TdcConfidence), f"Unexpected result type: {conf}"
    # TODO: assertions


def test_mixmax_confidence(simple_psms: PsmDataset):
    """Test that we can compute confidence with `TdcConfidence`"""
    conf = simple_psms.assign_confidence(
        score_column="x",
        method="mixmax",
        pep_fdr_type="psm-only",
    )

    assert isinstance(
        conf, MixmaxConfidence
    ), f"Unexpected result type: {conf}"
    # TODO: assertions


def test_mixmax_confidence_desc(simple_psms: PsmDataset):
    """Test that we can compute confidence with `TdcConfidence`"""
    simple_psms.data["x"] = -1.0 * simple_psms.data["x"]

    conf = simple_psms.assign_confidence(
        score_column="x",
        method="mixmax",
        pep_fdr_type="psm-only",
        desc=True,
    )

    assert isinstance(
        conf, MixmaxConfidence
    ), f"Unexpected result type: {conf}"
    # TODO: assertions
