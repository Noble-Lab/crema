"""
Edge-case tests that cut across multiple modules.
"""

import pytest
import numpy as np
import pandas as pd

from crema import PsmDataset

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_psms(rows, spectrum_columns=("file", "scan")):
    """Build a minimal PsmDataset from a list of dicts."""
    df = pd.DataFrame(rows)
    return PsmDataset(
        psms=df,
        target_column="target",
        spectrum_columns=list(spectrum_columns),
        score_columns=["score"],
        peptide_column="peptide",
        protein_column="protein",
        protein_delim=",",
    )


# ---------------------------------------------------------------------------
# Competition: best PSM per spectrum survives
# ---------------------------------------------------------------------------


def test_competition_target_wins_spectrum():
    """When target outscores the decoy, the decoy must be eliminated."""
    psms = _make_psms(
        [
            {
                "file": "f1",
                "scan": 1,
                "score": 0.9,
                "target": True,
                "peptide": "AAA",
                "protein": "P1",
            },
            {
                "file": "f1",
                "scan": 1,
                "score": 0.1,
                "target": False,
                "peptide": "AXA",
                "protein": "P1",
            },
            {
                "file": "f1",
                "scan": 2,
                "score": 0.2,
                "target": True,
                "peptide": "BBB",
                "protein": "P2",
            },
            {
                "file": "f1",
                "scan": 2,
                "score": 0.8,
                "target": False,
                "peptide": "BXB",
                "protein": "P2",
            },
        ]
    )
    conf = psms.assign_confidence(
        score_column="score",
        method="tdc",
        desc=True,
        pep_fdr_type="psm-only",
        threshold="q-value",
    )
    psm_df = conf.confidence_estimates["psms"]
    # scan 1 target survives; scan 2 decoy survives (in decoy estimates)
    assert 1 in psm_df["scan"].values
    assert 2 not in psm_df["scan"].values


def test_competition_decoy_wins_spectrum():
    """When decoy outscores the target, the decoy must be the sole survivor."""
    psms = _make_psms(
        [
            {
                "file": "f1",
                "scan": 1,
                "score": 0.9,
                "target": False,
                "peptide": "AXA",
                "protein": "P1",
            },
            {
                "file": "f1",
                "scan": 1,
                "score": 0.1,
                "target": True,
                "peptide": "AAA",
                "protein": "P1",
            },
            {
                "file": "f1",
                "scan": 2,
                "score": 0.8,
                "target": True,
                "peptide": "BBB",
                "protein": "P2",
            },
            {
                "file": "f1",
                "scan": 2,
                "score": 0.2,
                "target": False,
                "peptide": "BXB",
                "protein": "P2",
            },
        ]
    )
    conf = psms.assign_confidence(
        score_column="score",
        method="tdc",
        desc=True,
        pep_fdr_type="psm-only",
        threshold="q-value",
    )
    dec_df = conf.decoy_confidence_estimates["psms"]
    # scan 1 decoy should appear in decoy estimates
    assert 1 in dec_df["scan"].values


def test_competition_multiple_candidates_per_spectrum():
    """With 3 candidates for one spectrum, only the highest scorer survives."""
    psms = _make_psms(
        [
            {
                "file": "f1",
                "scan": 1,
                "score": 0.5,
                "target": True,
                "peptide": "AA",
                "protein": "P1",
            },
            {
                "file": "f1",
                "scan": 1,
                "score": 0.9,
                "target": True,
                "peptide": "BB",
                "protein": "P1",
            },
            {
                "file": "f1",
                "scan": 1,
                "score": 0.3,
                "target": False,
                "peptide": "XY",
                "protein": "P1",
            },
            # Second spectrum so we have at least one decoy surviving
            {
                "file": "f1",
                "scan": 2,
                "score": 0.2,
                "target": False,
                "peptide": "CC",
                "protein": "P2",
            },
            {
                "file": "f1",
                "scan": 2,
                "score": 0.1,
                "target": True,
                "peptide": "DD",
                "protein": "P2",
            },
        ]
    )
    conf = psms.assign_confidence(
        score_column="score",
        method="tdc",
        desc=True,
        pep_fdr_type="psm-only",
        threshold="q-value",
    )
    psm_df = conf.confidence_estimates["psms"]
    # scan 1 winner should be BB (score 0.9)
    scan1_rows = psm_df[psm_df["scan"] == 1]
    assert len(scan1_rows) == 1
    assert scan1_rows["score"].iloc[0] == pytest.approx(0.9)
    assert scan1_rows["peptide"].iloc[0] == "BB"


# ---------------------------------------------------------------------------
# Shared peptides excluded from protein-level FDR
# ---------------------------------------------------------------------------


def test_shared_peptide_excluded_from_proteins():
    """A peptide mapping to two proteins (via delim) must be dropped at protein level."""
    psms = _make_psms(
        [
            # Peptide maps to two proteins → must be excluded at protein level
            {
                "file": "f1",
                "scan": 1,
                "score": 0.9,
                "target": True,
                "peptide": "AAA",
                "protein": "P1,P2",
            },
            {
                "file": "f1",
                "scan": 1,
                "score": 0.1,
                "target": False,
                "peptide": "AXA",
                "protein": "P1,P2",
            },
            # Peptide maps to one protein → included
            {
                "file": "f1",
                "scan": 2,
                "score": 0.8,
                "target": True,
                "peptide": "BBB",
                "protein": "P3",
            },
            {
                "file": "f1",
                "scan": 2,
                "score": 0.2,
                "target": False,
                "peptide": "BXB",
                "protein": "P3",
            },
        ]
    )
    conf = psms.assign_confidence(
        score_column="score",
        method="tdc",
        desc=True,
        pep_fdr_type="psm-only",
        threshold="q-value",
    )
    prot_df = conf.confidence_estimates["proteins"]
    # "P1,P2" should not appear as a protein entry
    protein_values = prot_df["protein"].tolist()
    assert "P1,P2" not in protein_values


# ---------------------------------------------------------------------------
# All-in-one: presence of all four levels
# ---------------------------------------------------------------------------


def test_tdc_all_levels_present(clean_psms):
    """TdcConfidence must populate all four confidence levels."""
    conf = clean_psms.assign_confidence(
        score_column="score",
        method="tdc",
        desc=True,
        pep_fdr_type="psm-only",
    )
    for level in ("psms", "peptides", "proteins", "protein_groups"):
        assert level in conf.confidence_estimates, f"Missing level: {level}"


# ---------------------------------------------------------------------------
# Empty and degenerate construction cases
# ---------------------------------------------------------------------------


def test_empty_dataframe_raises():
    """An empty DataFrame must raise ValueError."""
    empty_df = pd.DataFrame(
        columns=["file", "scan", "score", "target", "peptide", "protein"]
    )
    with pytest.raises(ValueError):
        PsmDataset(
            psms=empty_df,
            target_column="target",
            spectrum_columns=["file", "scan"],
            score_columns=["score"],
            peptide_column="peptide",
            protein_column="protein",
            protein_delim=",",
        )


def test_all_targets_raises():
    """A dataset with no decoys must raise ValueError."""
    df = pd.DataFrame(
        [
            {
                "file": "f1",
                "scan": 1,
                "score": 0.9,
                "target": True,
                "peptide": "AA",
                "protein": "P1",
            },
            {
                "file": "f1",
                "scan": 2,
                "score": 0.8,
                "target": True,
                "peptide": "BB",
                "protein": "P2",
            },
        ]
    )
    with pytest.raises(ValueError, match="[Nn]o decoy"):
        PsmDataset(
            psms=df,
            target_column="target",
            spectrum_columns=["file", "scan"],
            score_columns=["score"],
            peptide_column="peptide",
            protein_column="protein",
            protein_delim=",",
        )


def test_all_decoys_raises():
    """A dataset with no targets must raise ValueError."""
    df = pd.DataFrame(
        [
            {
                "file": "f1",
                "scan": 1,
                "score": 0.9,
                "target": False,
                "peptide": "AX",
                "protein": "P1",
            },
            {
                "file": "f1",
                "scan": 2,
                "score": 0.8,
                "target": False,
                "peptide": "BX",
                "protein": "P2",
            },
        ]
    )
    with pytest.raises(ValueError, match="[Nn]o target"):
        PsmDataset(
            psms=df,
            target_column="target",
            spectrum_columns=["file", "scan"],
            score_columns=["score"],
            peptide_column="peptide",
            protein_column="protein",
            protein_delim=",",
        )


# ---------------------------------------------------------------------------
# Monotonicity invariant on real-ish data
# ---------------------------------------------------------------------------


def test_qvalues_monotone_at_all_levels(clean_psms):
    """Q-values must be non-decreasing as score decreases at every level."""
    conf = clean_psms.assign_confidence(
        score_column="score",
        method="tdc",
        desc=True,
        pep_fdr_type="psm-only",
        threshold="q-value",
    )
    for level in ("psms", "peptides"):
        df = conf.confidence_estimates[level]
        scores = df["score"].values
        qvals = df["crema q-value"].values
        for i in range(len(qvals) - 1):
            assert qvals[i] <= qvals[i + 1] or scores[i] == scores[i + 1], (
                f"Non-monotone q-values at level={level}: "
                f"score[{i}]={scores[i]}, q[{i}]={qvals[i]} > q[{i+1}]={qvals[i+1]}"
            )
