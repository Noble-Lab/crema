"""
Tests for the DuckDB backend (crema/backends/duckdb.py) and the
DuckdbTdcConfidence class.

Tests are skipped automatically when the optional 'duckdb' package is not
installed.
"""

import pytest
import pandas as pd

pytest.importorskip(
    "duckdb", reason="duckdb not installed; skipping DuckDB tests"
)

from crema.confidence import DuckdbTdcConfidence
from crema.dataset import PsmDataset
from crema.confidence import TdcConfidence

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


# ---------------------------------------------------------------------------
# Type / API
# ---------------------------------------------------------------------------


def test_duckdb_backend_returns_duckdb_confidence(simple_psms):
    """backend='duckdb' must return a DuckdbTdcConfidence instance."""
    conf = simple_psms.assign_confidence(
        score_column="x",
        method="tdc",
        pep_fdr_type="psm-only",
        desc=True,
        backend="duckdb",
    )
    assert isinstance(conf, DuckdbTdcConfidence)


def test_duckdb_backend_invalid_method_raises(simple_psms):
    """backend='duckdb' with method='mixmax' must raise ValueError."""
    with pytest.raises(ValueError, match="backend='duckdb'"):
        simple_psms.assign_confidence(
            score_column="x",
            method="mixmax",
            pep_fdr_type="psm-only",
            desc=True,
            backend="duckdb",
        )


def test_duckdb_backend_default_is_pandas(simple_psms):
    """Omitting backend must still return the pandas TdcConfidence."""
    conf = simple_psms.assign_confidence(
        score_column="x",
        method="tdc",
        pep_fdr_type="psm-only",
        desc=True,
    )
    assert isinstance(conf, TdcConfidence)
    assert not isinstance(conf, DuckdbTdcConfidence)


# ---------------------------------------------------------------------------
# Output structure
# ---------------------------------------------------------------------------


def test_duckdb_psm_level_present(simple_psms):
    """confidence_estimates must contain a 'psms' key."""
    conf = simple_psms.assign_confidence(
        score_column="x",
        pep_fdr_type="psm-only",
        desc=True,
        backend="duckdb",
    )
    assert "psms" in conf.confidence_estimates
    assert conf.confidence_estimates["psms"] is not None


def test_duckdb_peptide_level_present(simple_psms):
    """confidence_estimates must contain a 'peptides' key."""
    conf = simple_psms.assign_confidence(
        score_column="x",
        pep_fdr_type="psm-only",
        desc=True,
        backend="duckdb",
    )
    assert "peptides" in conf.confidence_estimates


def test_duckdb_protein_level_present(simple_psms):
    """confidence_estimates must contain a 'proteins' key."""
    conf = simple_psms.assign_confidence(
        score_column="x",
        pep_fdr_type="psm-only",
        desc=True,
        backend="duckdb",
    )
    assert "proteins" in conf.confidence_estimates


def test_duckdb_no_protein_groups(simple_psms):
    """DuckDB backend must not produce a 'protein_groups' level."""
    conf = simple_psms.assign_confidence(
        score_column="x",
        pep_fdr_type="psm-only",
        desc=True,
        backend="duckdb",
    )
    assert "protein_groups" not in conf.confidence_estimates


def test_duckdb_psm_output_has_score_and_accept(simple_psms):
    """PSM output must contain the score column and an accept column."""
    conf = simple_psms.assign_confidence(
        score_column="x",
        pep_fdr_type="psm-only",
        desc=True,
        threshold=0.5,
        backend="duckdb",
    )
    df = conf.confidence_estimates["psms"]
    assert "x" in df.columns
    assert "accept" in df.columns


def test_duckdb_psm_output_has_spectrum_columns(simple_psms):
    """PSM output must contain the spectrum-identifying columns."""
    conf = simple_psms.assign_confidence(
        score_column="x",
        pep_fdr_type="psm-only",
        desc=True,
        backend="duckdb",
    )
    df = conf.confidence_estimates["psms"]
    assert "scan" in df.columns
    assert "spectrum precursor m/z" in df.columns


def test_duckdb_protein_output_has_protein_column(simple_psms):
    """Protein output must contain the protein column."""
    conf = simple_psms.assign_confidence(
        score_column="x",
        pep_fdr_type="psm-only",
        desc=True,
        backend="duckdb",
    )
    df = conf.confidence_estimates["proteins"]
    assert "protein id" in df.columns


# ---------------------------------------------------------------------------
# Numerical correctness — DuckDB must agree with pandas TdcConfidence
# ---------------------------------------------------------------------------


def _run_both(psms, **kwargs):
    """Run both backends with the same arguments and return (pandas_conf, duckdb_conf)."""
    pandas_conf = psms.assign_confidence(**kwargs, backend="pandas")
    duckdb_conf = psms.assign_confidence(**kwargs, backend="duckdb")
    return pandas_conf, duckdb_conf


def test_duckdb_psm_qvalues_match_pandas(simple_psms):
    """PSM-level q-values from DuckDB must match the pandas path."""
    p_conf, d_conf = _run_both(
        simple_psms,
        score_column="x",
        pep_fdr_type="psm-only",
        desc=True,
        threshold="q-value",
    )
    p_df = (
        p_conf.confidence_estimates["psms"]
        .sort_values("x")
        .reset_index(drop=True)
    )
    d_df = (
        d_conf.confidence_estimates["psms"]
        .sort_values("x")
        .reset_index(drop=True)
    )
    pd.testing.assert_series_equal(
        p_df["crema q-value"],
        d_df["crema q-value"],
        check_names=False,
        atol=1e-6,
    )


def test_duckdb_psm_count_matches_pandas(simple_psms):
    """DuckDB and pandas must produce the same number of target PSMs."""
    p_conf, d_conf = _run_both(
        simple_psms,
        score_column="x",
        pep_fdr_type="psm-only",
        desc=True,
    )
    assert len(p_conf.confidence_estimates["psms"]) == len(
        d_conf.confidence_estimates["psms"]
    )


def test_duckdb_peptide_count_matches_pandas(simple_psms):
    """DuckDB and pandas must produce the same number of target peptides."""
    p_conf, d_conf = _run_both(
        simple_psms,
        score_column="x",
        pep_fdr_type="psm-only",
        desc=True,
    )
    assert len(p_conf.confidence_estimates["peptides"]) == len(
        d_conf.confidence_estimates["peptides"]
    )


def test_duckdb_accept_matches_pandas_at_threshold(simple_psms):
    """DuckDB 'accept' counts must agree with the pandas path at the same threshold."""
    threshold = 0.5
    p_conf, d_conf = _run_both(
        simple_psms,
        score_column="x",
        pep_fdr_type="psm-only",
        desc=True,
        threshold=threshold,
    )
    p_accept = p_conf.confidence_estimates["psms"]["accept"].sum()
    d_accept = d_conf.confidence_estimates["psms"]["accept"].sum()
    assert p_accept == d_accept


def test_duckdb_desc_false_matches_pandas(simple_psms):
    """DuckDB results for desc=False must match the pandas path."""
    p_conf, d_conf = _run_both(
        simple_psms,
        score_column="combined p-value",
        pep_fdr_type="psm-only",
        desc=False,
    )
    assert len(p_conf.confidence_estimates["psms"]) == len(
        d_conf.confidence_estimates["psms"]
    )


# ---------------------------------------------------------------------------
# to_txt compatibility
# ---------------------------------------------------------------------------


def test_duckdb_to_txt_creates_files(simple_psms, tmp_path):
    """to_txt must work on a DuckdbTdcConfidence object."""
    conf = simple_psms.assign_confidence(
        score_column="x",
        pep_fdr_type="psm-only",
        desc=True,
        backend="duckdb",
    )
    paths = conf.to_txt(output_dir=tmp_path, file_root="duckdb")
    assert (tmp_path / "duckdb.crema.psms.txt").exists()
    assert (tmp_path / "duckdb.crema.peptides.txt").exists()
    assert (tmp_path / "duckdb.crema.proteins.txt").exists()
    assert isinstance(paths, list)
    assert len(paths) >= 3
