"""
These are unit tests for the PSM Dataset Class:
"""
import pytest
import numpy as np
import pandas as pd
from crema import PsmDataset

from ..fixtures import *


@pytest.fixture
def simple_df(basic_crux_df):
    """A simple dataframe of PSMs"""
    df = basic_crux_df
    df["target"] = df["target/decoy"].replace({"target": True, "decoy": False})
    return df.drop(columns="target/decoy")


def test_create_object(simple_df):
    """Ensures that a PsmDataset object can be initialized properly."""
    psms = PsmDataset(
        psms=simple_df,
        target_column="target",
        spectrum_columns=["scan", "spectrum precursor m/z"],
        score_columns=["combined p-value", "x"],
        peptide_column="sequence",
    )
    assert isinstance(psms, PsmDataset)


def test_properties(simple_df):
    """Check the public properties of the PsmDataset object."""
    psms = PsmDataset(
        psms=simple_df,
        target_column="target",
        spectrum_columns=["scan", "spectrum precursor m/z"],
        score_columns=["combined p-value", "x"],
        peptide_column="sequence",
    )

    pd.testing.assert_frame_equal(psms.data, simple_df, check_like=True)
    assert psms.spectrum_columns == ["scan", "spectrum precursor m/z"]
    assert psms.score_columns == ["combined p-value", "x"]
    assert psms.peptide_column == "sequence"
    pd.testing.assert_frame_equal(
        psms.scores, simple_df.loc[:, ["combined p-value", "x"]]
    )
    np.testing.assert_array_equal(
        psms.targets, simple_df.loc[:, "target"].values
    )


def test_getitem(simple_df):
    """Test the __getitem__ method"""
    psms = PsmDataset(
        psms=simple_df,
        target_column="target",
        spectrum_columns=["scan", "spectrum precursor m/z"],
        score_columns=["combined p-value", "x"],
        peptide_column="sequence",
    )

    pd.testing.assert_series_equal(psms["scan"], simple_df["scan"])
    pd.testing.assert_series_equal(
        psms["combined p-value"], simple_df["combined p-value"],
    )


def test_find_best_score(simple_df):
    """Test the find_best_score() method"""
    psms = PsmDataset(
        psms=simple_df,
        target_column="target",
        spectrum_columns=["scan", "spectrum precursor m/z"],
        score_columns=["combined p-value", "x"],
        peptide_column="sequence",
    )

    score, npass, desc = psms.find_best_score(eval_fdr=0.4)
    assert score == "x"
    assert npass == 4
    assert desc