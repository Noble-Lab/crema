"""
Tests for confidence estimation (TdcConfidence and MixmaxConfidence).

Includes type-checks carried over from the original file, plus numerical
correctness tests, invariant checks, and parameter-validation tests.
"""

import pytest
import numpy as np

import crema
from crema.confidence import TdcConfidence, MixmaxConfidence
from crema.dataset import PsmDataset

from .test_dataset import simple_df  # noqa: F401

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def simple_psms(simple_df):  # noqa: F811
    return PsmDataset(
        psms=simple_df,
        target_column="target",
        spectrum_columns=["scan", "spectrum precursor m/z"],
        score_columns=["combined p-value", "x"],
        peptide_column="sequence",
        protein_column="protein id",
        protein_delim=",",
    )


@pytest.fixture
def psms_with_pairing(mod_target_tide_txt, mod_decoy_tide_txt):
    """PsmDataset that includes a peptide pairing dict (from Tide files)."""
    return crema.read_tide([mod_decoy_tide_txt, mod_target_tide_txt])


# ---------------------------------------------------------------------------
# Existing type-check tests (preserved)
# ---------------------------------------------------------------------------


def test_tdc_confidence(simple_psms: PsmDataset):
    conf = simple_psms.assign_confidence(
        score_column="x",
        method="tdc",
        pep_fdr_type="psm-only",
    )
    assert isinstance(conf, TdcConfidence)


def test_mixmax_confidence(simple_psms: PsmDataset):
    conf = simple_psms.assign_confidence(
        score_column="x",
        method="mixmax",
        pep_fdr_type="psm-only",
    )
    assert isinstance(conf, MixmaxConfidence)


def test_mixmax_confidence_desc(simple_df):  # noqa: F811
    df = simple_df.copy()
    df["x"] = -1.0 * df["x"]
    psms = PsmDataset(
        psms=df,
        target_column="target",
        spectrum_columns=["scan", "spectrum precursor m/z"],
        score_columns=["combined p-value", "x"],
        peptide_column="sequence",
        protein_column="protein id",
        protein_delim=",",
    )
    conf = psms.assign_confidence(
        score_column="x",
        method="mixmax",
        pep_fdr_type="psm-only",
        desc=True,
    )
    assert isinstance(conf, MixmaxConfidence)


# ---------------------------------------------------------------------------
# Numerical correctness - PSM level (uses clean_psms from conftest.py)
# ---------------------------------------------------------------------------


def test_tdc_psm_qvalues_exact(clean_psms):
    """All four winning targets must have q-value == 0.25 (hand-computed)."""
    conf = clean_psms.assign_confidence(
        score_column="score",
        method="tdc",
        desc=True,
        pep_fdr_type="psm-only",
        threshold="q-value",
    )
    psm_df = conf.confidence_estimates["psms"]
    qvals = psm_df["crema q-value"].values
    np.testing.assert_allclose(qvals, 0.25, atol=1e-9)


def test_tdc_psm_target_count(clean_psms):
    """After competition 4 targets survive; all must appear in estimates."""
    conf = clean_psms.assign_confidence(
        score_column="score",
        method="tdc",
        desc=True,
        pep_fdr_type="psm-only",
        threshold="q-value",
    )
    assert len(conf.confidence_estimates["psms"]) == 4


def test_tdc_psm_decoy_count(clean_psms):
    """After competition 2 decoys survive and appear in decoy estimates."""
    conf = clean_psms.assign_confidence(
        score_column="score",
        method="tdc",
        desc=True,
        pep_fdr_type="psm-only",
        threshold="q-value",
    )
    assert len(conf.decoy_confidence_estimates["psms"]) == 2


def test_tdc_psm_count_at_fdr_threshold(clean_psms):
    """At FDR=0.25 all 4 targets pass; at FDR=0.01 none do."""
    conf_25 = clean_psms.assign_confidence(
        score_column="score",
        method="tdc",
        desc=True,
        pep_fdr_type="psm-only",
        threshold=0.25,
    )
    conf_01 = clean_psms.assign_confidence(
        score_column="score",
        method="tdc",
        desc=True,
        pep_fdr_type="psm-only",
        threshold=0.01,
    )
    assert conf_25.confidence_estimates["psms"]["accept"].sum() == 4
    assert conf_01.confidence_estimates["psms"]["accept"].sum() == 0


def test_tdc_psm_scores_descending(clean_psms):
    """PSM confidence estimates must be sorted best-score-first."""
    conf = clean_psms.assign_confidence(
        score_column="score",
        method="tdc",
        desc=True,
        pep_fdr_type="psm-only",
        threshold="q-value",
    )
    scores = conf.confidence_estimates["psms"]["score"].values
    assert list(scores) == sorted(scores, reverse=True)


# ---------------------------------------------------------------------------
# Q-value invariants across all levels
# ---------------------------------------------------------------------------


def _assert_qvals_valid(df, qval_col="crema q-value"):
    """Helper: q-values in [0,1] and non-decreasing as score decreases."""
    qvals = df[qval_col].values
    scores = df["score"].values if "score" in df.columns else None
    assert (qvals >= 0).all(), "Q-values must be >= 0"
    assert (qvals <= 1).all(), "Q-values must be <= 1"
    if scores is not None:
        # rows are sorted best-score-first; q-values must be non-decreasing
        assert all(
            qvals[i] <= qvals[i + 1] for i in range(len(qvals) - 1)
        ), "Q-values not monotone"


def test_psm_qvalues_valid(clean_psms):
    conf = clean_psms.assign_confidence(
        score_column="score",
        method="tdc",
        desc=True,
        pep_fdr_type="psm-only",
        threshold="q-value",
    )
    _assert_qvals_valid(conf.confidence_estimates["psms"])


def test_peptide_qvalues_valid(clean_psms):
    conf = clean_psms.assign_confidence(
        score_column="score",
        method="tdc",
        desc=True,
        pep_fdr_type="psm-only",
        threshold="q-value",
    )
    _assert_qvals_valid(conf.confidence_estimates["peptides"])


def test_protein_qvalues_valid(clean_psms):
    conf = clean_psms.assign_confidence(
        score_column="score",
        method="tdc",
        desc=True,
        pep_fdr_type="psm-only",
        threshold="q-value",
    )
    prot_df = conf.confidence_estimates["proteins"]
    qvals = prot_df["crema q-value"].values
    assert (qvals >= 0).all()
    assert (qvals <= 1).all()


# ---------------------------------------------------------------------------
# Protein level correctness
# ---------------------------------------------------------------------------


def test_protein_level_unique_proteins(clean_psms):
    """Each protein must appear at most once in protein-level estimates."""
    conf = clean_psms.assign_confidence(
        score_column="score",
        method="tdc",
        desc=True,
        pep_fdr_type="psm-only",
        threshold="q-value",
    )
    prot_df = conf.confidence_estimates["proteins"]
    assert prot_df["protein"].nunique() == len(prot_df)


def test_protein_level_has_score_column(clean_psms):
    """Protein-level results must include the score column used."""
    conf = clean_psms.assign_confidence(
        score_column="score",
        method="tdc",
        desc=True,
        pep_fdr_type="psm-only",
        threshold="q-value",
    )
    assert "score" in conf.confidence_estimates["proteins"].columns


@pytest.mark.parametrize("prot_fdr_type", ["best", "combine"])
def test_protein_score_aggregation(clean_psms, prot_fdr_type):
    """Both protein aggregation methods must run without error."""
    conf = clean_psms.assign_confidence(
        score_column="score",
        method="tdc",
        desc=True,
        pep_fdr_type="psm-only",
        prot_fdr_type=prot_fdr_type,
        threshold="q-value",
    )
    assert "proteins" in conf.confidence_estimates


# ---------------------------------------------------------------------------
# pep_fdr_type parameter handling
# ---------------------------------------------------------------------------


def test_pep_fdr_type_psm_only_no_pairing_required(clean_psms):
    """psm-only must work when no peptide_pairing is provided."""
    assert clean_psms.peptide_pairing is None
    conf = clean_psms.assign_confidence(
        score_column="score",
        method="tdc",
        desc=True,
        pep_fdr_type="psm-only",
    )
    assert "peptides" in conf.confidence_estimates


def test_pep_fdr_type_psm_peptide_requires_pairing(clean_psms):
    """psm-peptide without peptide_pairing must raise ValueError."""
    with pytest.raises(ValueError):
        clean_psms.assign_confidence(
            score_column="score",
            method="tdc",
            desc=True,
            pep_fdr_type="psm-peptide",
        )


def test_pep_fdr_type_peptide_only_requires_pairing(clean_psms):
    """peptide-only without peptide_pairing must raise ValueError."""
    with pytest.raises(ValueError):
        clean_psms.assign_confidence(
            score_column="score",
            method="tdc",
            desc=True,
            pep_fdr_type="peptide-only",
        )


@pytest.mark.parametrize(
    "pep_fdr_type", ["psm-only", "psm-peptide", "peptide-only"]
)
def test_pep_fdr_type_with_pairing(psms_with_pairing, pep_fdr_type):
    """All pep_fdr_type options must produce a peptides level result."""
    conf = psms_with_pairing.assign_confidence(
        score_column="combined p-value",
        method="tdc",
        desc=False,
        pep_fdr_type=pep_fdr_type,
        eval_fdr=0.5,
        threshold="q-value",
    )
    assert "peptides" in conf.confidence_estimates
    assert len(conf.confidence_estimates["peptides"]) > 0


# ---------------------------------------------------------------------------
# Invalid parameter validation
# ---------------------------------------------------------------------------


def test_invalid_pep_fdr_type_raises(clean_psms):
    with pytest.raises(ValueError):
        clean_psms.assign_confidence(
            score_column="score",
            method="tdc",
            pep_fdr_type="invalid_option",
        )


def test_invalid_prot_fdr_type_raises(clean_psms):
    with pytest.raises(ValueError):
        clean_psms.assign_confidence(
            score_column="score",
            method="tdc",
            pep_fdr_type="psm-only",
            prot_fdr_type="invalid_option",
        )


def test_eval_fdr_out_of_range_raises(clean_psms):
    with pytest.raises(ValueError):
        clean_psms.assign_confidence(
            score_column="score",
            method="tdc",
            pep_fdr_type="psm-only",
            eval_fdr=1.5,
        )


# ---------------------------------------------------------------------------
# Threshold parameter: "q-value" vs float
# ---------------------------------------------------------------------------


def test_threshold_float_gives_accept_column(clean_psms):
    """threshold=0.5 must produce an 'accept' boolean column."""
    conf = clean_psms.assign_confidence(
        score_column="score",
        method="tdc",
        desc=True,
        pep_fdr_type="psm-only",
        threshold=0.5,
    )
    assert "accept" in conf.confidence_estimates["psms"].columns
    assert conf.confidence_estimates["psms"]["accept"].dtype == bool


def test_threshold_qvalue_gives_qvalue_column(clean_psms):
    """threshold='q-value' must produce a 'crema q-value' float column."""
    conf = clean_psms.assign_confidence(
        score_column="score",
        method="tdc",
        desc=True,
        pep_fdr_type="psm-only",
        threshold="q-value",
    )
    assert "crema q-value" in conf.confidence_estimates["psms"].columns
    assert "accept" not in conf.confidence_estimates["psms"].columns


# ---------------------------------------------------------------------------
# Confidence levels present
# ---------------------------------------------------------------------------


def test_tdc_all_four_levels_present(clean_psms):
    """TDC must populate all four confidence levels."""
    conf = clean_psms.assign_confidence(
        score_column="score",
        method="tdc",
        desc=True,
        pep_fdr_type="psm-only",
    )
    for level in ("psms", "peptides", "proteins", "protein_groups"):
        assert level in conf.confidence_estimates


def test_mixmax_only_psm_level(simple_psms):
    """MixmaxConfidence must only produce PSM-level estimates."""
    conf = simple_psms.assign_confidence(
        score_column="x",
        method="mixmax",
        pep_fdr_type="psm-only",
    )
    assert "psms" in conf.confidence_estimates
    assert "peptides" not in conf.confidence_estimates
    assert "proteins" not in conf.confidence_estimates


# ---------------------------------------------------------------------------
# Targets only in confidence_estimates; decoys in decoy_confidence_estimates
# ---------------------------------------------------------------------------


def test_confidence_estimates_contains_only_targets(clean_psms):
    """confidence_estimates must contain only target PSMs.

    Note: _prettify_tables strips the target column from output DataFrames,
    so we use the peptide column to distinguish targets (PEP1-PEP6) from
    decoys (PEP1D-PEP6D) in the clean_psms fixture.
    """
    conf = clean_psms.assign_confidence(
        score_column="score",
        method="tdc",
        desc=True,
        pep_fdr_type="psm-only",
        threshold="q-value",
    )
    psm_df = conf.confidence_estimates["psms"]
    # In clean_psm_df, decoy peptides all end in "D" (e.g. PEP1D)
    assert not psm_df["peptide"].str.endswith("D").any()


def test_decoy_confidence_estimates_contains_only_decoys(clean_psms):
    """decoy_confidence_estimates must contain only decoy PSMs.

    Note: _prettify_tables strips the target column from output DataFrames,
    so we use the peptide column to distinguish targets (PEP1-PEP6) from
    decoys (PEP1D-PEP6D) in the clean_psms fixture.
    """
    conf = clean_psms.assign_confidence(
        score_column="score",
        method="tdc",
        desc=True,
        pep_fdr_type="psm-only",
        threshold="q-value",
    )
    dec_df = conf.decoy_confidence_estimates["psms"]
    # In clean_psm_df, decoy peptides all end in "D" (e.g. PEP5D, PEP6D)
    assert dec_df["peptide"].str.endswith("D").all()
