"""
Tests for writers/txt.py (the to_txt output writer).
"""

import pandas as pd
from pathlib import Path

# ---------------------------------------------------------------------------
# File creation
# ---------------------------------------------------------------------------


def test_to_txt_creates_psm_peptide_protein_files(clean_psms, tmp_path):
    """to_txt must create psms, peptides, and proteins output files."""
    conf = clean_psms.assign_confidence(
        score_column="score",
        desc=True,
        pep_fdr_type="psm-only",
    )
    conf.to_txt(output_dir=tmp_path)
    assert (tmp_path / "crema.psms.txt").exists()
    assert (tmp_path / "crema.peptides.txt").exists()
    assert (tmp_path / "crema.proteins.txt").exists()


def test_to_txt_file_root_prefix(clean_psms, tmp_path):
    """A file_root prefix must appear in every output file name."""
    conf = clean_psms.assign_confidence(
        score_column="score",
        desc=True,
        pep_fdr_type="psm-only",
    )
    conf.to_txt(output_dir=tmp_path, file_root="myrun")
    assert (tmp_path / "myrun.crema.psms.txt").exists()
    assert (tmp_path / "myrun.crema.peptides.txt").exists()
    assert (tmp_path / "myrun.crema.proteins.txt").exists()


def test_to_txt_returns_file_paths(clean_psms, tmp_path):
    """to_txt must return a list of the created file paths."""
    conf = clean_psms.assign_confidence(
        score_column="score",
        desc=True,
        pep_fdr_type="psm-only",
    )
    paths = conf.to_txt(output_dir=tmp_path)
    assert isinstance(paths, list)
    assert len(paths) >= 3
    for p in paths:
        assert Path(p).exists()


def test_to_txt_decoys_flag_creates_decoy_files(clean_psms, tmp_path):
    """decoys=True must also write decoy output files."""
    conf = clean_psms.assign_confidence(
        score_column="score",
        desc=True,
        pep_fdr_type="psm-only",
    )
    conf.to_txt(output_dir=tmp_path, decoys=True)
    assert (tmp_path / "crema.decoy.psms.txt").exists()


# ---------------------------------------------------------------------------
# Column names in output
# ---------------------------------------------------------------------------


def test_to_txt_psm_file_has_expected_columns(clean_psms, tmp_path):
    """The PSM output file must contain the score and accept columns."""
    conf = clean_psms.assign_confidence(
        score_column="score",
        desc=True,
        pep_fdr_type="psm-only",
        threshold=0.5,
    )
    conf.to_txt(output_dir=tmp_path)
    result = pd.read_csv(tmp_path / "crema.psms.txt", sep="\t")
    assert "score" in result.columns
    assert "accept" in result.columns


def test_to_txt_psm_file_has_spectrum_columns(clean_psms, tmp_path):
    """The PSM output file must contain spectrum-identifying columns."""
    conf = clean_psms.assign_confidence(
        score_column="score",
        desc=True,
        pep_fdr_type="psm-only",
    )
    conf.to_txt(output_dir=tmp_path)
    result = pd.read_csv(tmp_path / "crema.psms.txt", sep="\t")
    assert "file" in result.columns
    assert "scan" in result.columns


def test_to_txt_protein_file_has_protein_column(clean_psms, tmp_path):
    """The protein output file must contain the protein column."""
    conf = clean_psms.assign_confidence(
        score_column="score",
        desc=True,
        pep_fdr_type="psm-only",
    )
    conf.to_txt(output_dir=tmp_path)
    result = pd.read_csv(tmp_path / "crema.proteins.txt", sep="\t")
    assert "protein" in result.columns


# ---------------------------------------------------------------------------
# Content correctness
# ---------------------------------------------------------------------------


def test_to_txt_psm_row_count(clean_psms, tmp_path):
    """PSM output file must have exactly as many rows as the in-memory result."""
    conf = clean_psms.assign_confidence(
        score_column="score",
        desc=True,
        pep_fdr_type="psm-only",
    )
    conf.to_txt(output_dir=tmp_path)
    result = pd.read_csv(tmp_path / "crema.psms.txt", sep="\t")
    assert len(result) == len(conf.confidence_estimates["psms"])


def test_to_txt_roundtrip_scores(clean_psms, tmp_path):
    """Scores written to file and read back must match in-memory values."""
    conf = clean_psms.assign_confidence(
        score_column="score",
        desc=True,
        pep_fdr_type="psm-only",
        threshold=0.5,
    )
    conf.to_txt(output_dir=tmp_path)
    result = pd.read_csv(tmp_path / "crema.psms.txt", sep="\t")

    in_memory = conf.confidence_estimates["psms"].reset_index(drop=True)
    pd.testing.assert_series_equal(
        result["score"].reset_index(drop=True),
        in_memory["score"].reset_index(drop=True),
        check_names=False,
        atol=1e-5,
    )


def test_to_txt_threshold_filters_accept_column(clean_psms, tmp_path):
    """accept column must be True iff q-value <= threshold."""
    # At threshold=0.25, all 4 targets pass
    conf = clean_psms.assign_confidence(
        score_column="score",
        desc=True,
        pep_fdr_type="psm-only",
        threshold=0.25,
    )
    conf.to_txt(output_dir=tmp_path)
    result = pd.read_csv(tmp_path / "crema.psms.txt", sep="\t")
    assert result["accept"].sum() == 4

    # At threshold=0.01, no targets pass
    conf2 = clean_psms.assign_confidence(
        score_column="score",
        desc=True,
        pep_fdr_type="psm-only",
        threshold=0.01,
    )
    conf2.to_txt(output_dir=tmp_path, file_root="strict")
    result2 = pd.read_csv(tmp_path / "strict.crema.psms.txt", sep="\t")
    assert result2["accept"].sum() == 0


def test_to_txt_multi_confidence_row_count(clean_psms, tmp_path):
    """Writing a tuple of Confidence objects must produce the combined row count."""
    import crema

    conf = clean_psms.assign_confidence(
        score_column="score",
        desc=True,
        pep_fdr_type="psm-only",
        threshold=0.5,
    )
    conf.to_txt(output_dir=tmp_path, file_root="single")
    crema.to_txt((conf, conf), output_dir=tmp_path, file_root="combined")

    single = pd.read_csv(tmp_path / "single.crema.psms.txt", sep="\t")
    combined = pd.read_csv(tmp_path / "combined.crema.psms.txt", sep="\t")
    assert len(combined) == 2 * len(single)
