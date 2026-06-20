"""
Tests for Parquet I/O: crema/parsers/parquet.py and crema/writers/parquet.py.

Tests are skipped automatically when pyarrow is not installed.
"""

import pytest
import pandas as pd

pytest.importorskip(
    "pyarrow", reason="pyarrow not installed; skipping Parquet tests"
)

from crema.parsers.parquet import read_parquet
from crema.dataset import PsmDataset

from .test_dataset import simple_df  # noqa: F401

# ---------------------------------------------------------------------------
# Shared fixture: write a small PsmDataset to Parquet and read it back
# ---------------------------------------------------------------------------


@pytest.fixture
def parquet_file(simple_df, tmp_path):  # noqa: F811
    """Write simple_df to a Parquet file and return its path."""
    p = tmp_path / "psms.parquet"
    simple_df.to_parquet(p, index=False)
    return p


@pytest.fixture
def parquet_psms(parquet_file):
    """PsmDataset loaded from a Parquet file."""
    return read_parquet(
        str(parquet_file),
        target_column="target",
        spectrum_columns=["scan", "spectrum precursor m/z"],
        score_columns=["combined p-value", "x"],
        peptide_column="sequence",
        protein_column="protein id",
        protein_delim=",",
    )


# ---------------------------------------------------------------------------
# read_parquet — API and content
# ---------------------------------------------------------------------------


def test_read_parquet_returns_psm_dataset(parquet_psms):
    """read_parquet must return a PsmDataset."""
    assert isinstance(parquet_psms, PsmDataset)


def test_read_parquet_row_count(simple_df, parquet_psms):  # noqa: F811
    """Row count must match the original DataFrame."""
    assert len(parquet_psms.data) == len(simple_df)


def test_read_parquet_score_columns(parquet_psms):
    """Score columns must be preserved."""
    assert "combined p-value" in parquet_psms.data.columns
    assert "x" in parquet_psms.data.columns


def test_read_parquet_target_column_is_bool(parquet_psms):
    """Target column must be boolean after parsing."""
    assert parquet_psms.data["target"].dtype == bool


def test_read_parquet_has_targets_and_decoys(parquet_psms):
    """Dataset must contain both targets and decoys."""
    assert parquet_psms._num_targets > 0
    assert parquet_psms._num_decoys > 0


def test_read_parquet_missing_pyarrow(monkeypatch, parquet_file):
    """ImportError must be raised when pyarrow is unavailable."""
    import builtins

    real_import = builtins.__import__

    def mock_import(name, *args, **kwargs):
        if name == "pyarrow.parquet":
            raise ImportError("mocked")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", mock_import)
    with pytest.raises(ImportError, match="pyarrow"):
        read_parquet(
            str(parquet_file),
            target_column="target",
            spectrum_columns=["scan", "spectrum precursor m/z"],
            score_columns=["combined p-value", "x"],
            peptide_column="sequence",
            protein_column="protein id",
            protein_delim=",",
        )


def test_read_parquet_multiple_files(simple_df, tmp_path):  # noqa: F811
    """Reading two Parquet files must concatenate them."""
    p1 = tmp_path / "a.parquet"
    p2 = tmp_path / "b.parquet"
    simple_df.to_parquet(p1, index=False)
    simple_df.to_parquet(p2, index=False)
    psms = read_parquet(
        [str(p1), str(p2)],
        target_column="target",
        spectrum_columns=["scan", "spectrum precursor m/z"],
        score_columns=["combined p-value", "x"],
        peptide_column="sequence",
        protein_column="protein id",
        protein_delim=",",
    )
    assert len(psms.data) == 2 * len(simple_df)


# ---------------------------------------------------------------------------
# to_parquet (writer) — file creation and content
# ---------------------------------------------------------------------------


def test_to_parquet_creates_psm_peptide_protein_files(parquet_psms, tmp_path):
    """to_parquet must create psms, peptides, and proteins output files."""
    conf = parquet_psms.assign_confidence(
        score_column="x",
        method="tdc",
        desc=True,
        pep_fdr_type="psm-only",
    )
    conf.to_parquet(output_dir=tmp_path)
    assert (tmp_path / "crema.psms.parquet").exists()
    assert (tmp_path / "crema.peptides.parquet").exists()
    assert (tmp_path / "crema.proteins.parquet").exists()


def test_to_parquet_file_root_prefix(parquet_psms, tmp_path):
    """A file_root prefix must appear in every output file name."""
    conf = parquet_psms.assign_confidence(
        score_column="x",
        method="tdc",
        desc=True,
        pep_fdr_type="psm-only",
    )
    conf.to_parquet(output_dir=tmp_path, file_root="myrun")
    assert (tmp_path / "myrun.crema.psms.parquet").exists()
    assert (tmp_path / "myrun.crema.peptides.parquet").exists()
    assert (tmp_path / "myrun.crema.proteins.parquet").exists()


def test_to_parquet_returns_file_paths(parquet_psms, tmp_path):
    """to_parquet must return a list of created file paths."""
    conf = parquet_psms.assign_confidence(
        score_column="x",
        method="tdc",
        desc=True,
        pep_fdr_type="psm-only",
    )
    paths = conf.to_parquet(output_dir=tmp_path)
    assert isinstance(paths, list)
    assert len(paths) >= 3
    for p in paths:
        from pathlib import Path

        assert Path(p).exists()


def test_to_parquet_psm_row_count(parquet_psms, tmp_path):
    """PSM Parquet file must have the same row count as the in-memory result."""
    conf = parquet_psms.assign_confidence(
        score_column="x",
        method="tdc",
        desc=True,
        pep_fdr_type="psm-only",
    )
    conf.to_parquet(output_dir=tmp_path)
    result = pd.read_parquet(tmp_path / "crema.psms.parquet")
    assert len(result) == len(conf.confidence_estimates["psms"])


def test_to_parquet_psm_has_score_and_accept(parquet_psms, tmp_path):
    """PSM Parquet output must contain score and accept columns."""
    conf = parquet_psms.assign_confidence(
        score_column="x",
        method="tdc",
        desc=True,
        pep_fdr_type="psm-only",
        threshold=0.5,
    )
    conf.to_parquet(output_dir=tmp_path)
    result = pd.read_parquet(tmp_path / "crema.psms.parquet")
    assert "x" in result.columns
    assert "accept" in result.columns


def test_to_parquet_decoys_flag(parquet_psms, tmp_path):
    """decoys=True must also write decoy Parquet files."""
    conf = parquet_psms.assign_confidence(
        score_column="x",
        method="tdc",
        desc=True,
        pep_fdr_type="psm-only",
    )
    conf.to_parquet(output_dir=tmp_path, decoys=True)
    assert (tmp_path / "crema.decoy.psms.parquet").exists()


def test_to_parquet_missing_pyarrow(parquet_psms, tmp_path, monkeypatch):
    """ImportError must be raised when pyarrow is unavailable."""
    import builtins

    real_import = builtins.__import__

    def mock_import(name, *args, **kwargs):
        if name == "pyarrow":
            raise ImportError("mocked")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", mock_import)
    conf = parquet_psms.assign_confidence(
        score_column="x",
        method="tdc",
        desc=True,
        pep_fdr_type="psm-only",
    )
    with pytest.raises(ImportError, match="pyarrow"):
        conf.to_parquet(output_dir=tmp_path)


# ---------------------------------------------------------------------------
# Roundtrip: write Parquet → read_parquet → assign_confidence
# ---------------------------------------------------------------------------


def test_parquet_roundtrip_row_count(parquet_psms, tmp_path):
    """Confidence output written as Parquet must roundtrip with correct row count."""
    conf = parquet_psms.assign_confidence(
        score_column="x",
        method="tdc",
        desc=True,
        pep_fdr_type="psm-only",
    )
    conf.to_parquet(output_dir=tmp_path)
    result = pd.read_parquet(tmp_path / "crema.psms.parquet")
    assert len(result) == len(conf.confidence_estimates["psms"])


def test_parquet_roundtrip_scores_match_txt(parquet_psms, tmp_path):
    """Parquet and txt outputs must contain the same score values."""
    conf = parquet_psms.assign_confidence(
        score_column="x",
        method="tdc",
        desc=True,
        pep_fdr_type="psm-only",
    )
    conf.to_parquet(output_dir=tmp_path, file_root="pq")
    conf.to_txt(output_dir=tmp_path, file_root="tx")

    pq_df = (
        pd.read_parquet(tmp_path / "pq.crema.psms.parquet")
        .sort_values("x")
        .reset_index(drop=True)
    )
    tx_df = (
        pd.read_csv(tmp_path / "tx.crema.psms.txt", sep="\t")
        .sort_values("x")
        .reset_index(drop=True)
    )

    pd.testing.assert_series_equal(
        pq_df["x"], tx_df["x"], check_names=False, atol=1e-5
    )
