"""
These are unit tests for functions within parsers.py:
"""
import os
import pytest
import pandas as pd

import crema

from ..fixtures import *


def test_read_crux(real_crux_txt):
    """Test that we parse crux files correctly"""
    psms = crema.read_crux(real_crux_txt)
    assert isinstance(psms.data, pd.DataFrame)
    assert psms.data.shape == (21818, 11)
    assert psms.spectrum_columns == ["scan", "spectrum precursor m/z"]

    scores = {
        "refactored xcorr",
        "combined p-value",
        "exact p-value",
        "sp score",
        "delta_cn",
        "delta_lcn",
        "res-ev p-value",
    }
    assert set(psms.score_columns) == scores
    assert psms.scores.shape == (21818, len(scores))
    assert psms.targets.shape == (21818,)
    assert psms.targets.sum() == 10909
    assert (~psms.targets).sum() == 10909


def test_read_txt(basic_crux_csv):
    """Test that we can read generic delimited files"""
    psms = crema.read_txt(
        basic_crux_csv,
        target_column="target/decoy",
        spectrum_columns="scan",
        score_columns=["combined p-value", "x"],
        peptide_column="sequence",
        sep=",",
    )
    assert isinstance(psms.data, pd.DataFrame)
    assert psms.data.shape == (10, 5)
    assert psms.spectrum_columns == ["scan"]
    assert set(psms.score_columns) == {"combined p-value", "x"}
    assert psms.peptide_column == "sequence"
    assert psms.targets.shape == (10,)
    assert psms.targets.sum() == 6
    assert (~psms.targets).sum() == 4


def test_read_mztab(real_mztab):
    psms = crema.read_mztab(real_mztab)


def download_msv(msv_id, tmp_path):
    dat = ppx.MSVDataset(msv_id)
    mztab_file = dat.list_files("ccms_result")[0]
    file_path = "ccms_result/" + mztab_file
    file = dat.download(files=file_path, dest_dir=tmp_path)[0]
    return file


def verify_mztab(msv_id, num_scores, expected_df_head, tmp_path):
    file = download_msv(msv_id, tmp_path)
    psm_dataset = read_mztab(file)
    # result = calculate_tdc(psm_dataset, 2)
    # print(result.data)
    assert isinstance(psm_dataset, PsmDataset)
    assert psm_dataset.spectrum_col == ["spectra_ref"]
    score_col = []
    for i in range(1, num_scores + 1):
        score_col.append("search_engine_score[" + str(i) + "]")
    assert psm_dataset.score_col == score_col
    assert psm_dataset.target_col == "opt_global_cv_MS:1002217_decoy_peptide"
    pd.testing.assert_frame_equal(psm_dataset.data.head(), expected_df_head)
