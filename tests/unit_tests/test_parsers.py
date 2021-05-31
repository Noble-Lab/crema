"""
These are unit tests for functions within parsers.py:
"""
import pytest
import pandas as pd

import crema


def test_read_crux(real_crux_txt):
    """Test that we parse crux files correctly"""
    psms = crema.read_crux(real_crux_txt)
    assert isinstance(psms.data, pd.DataFrame)
    assert psms.data.shape == (21818, 11)
    assert list(psms.spectra.columns) == ["scan", "spectrum precursor m/z"]

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


def test_read_crux_peptide_pairing(mod_target_crux_txt, mod_decoy_crux_txt):
    """Test that peptide pairing is correctly created when parsing crux files"""
    expected_peptide_pairing = {
        "ALLSLR": "ALLLSR",
        "QTPPAR": "QTAPPR",
        "GEVPN[0.98]R": "GN[0.98]PEVR",
        "GGHMDR": "GDMGHR",
    }
    psms = crema.read_crux(
        [mod_decoy_crux_txt, mod_target_crux_txt],
        peptide_tdc=True,
    )
    assert expected_peptide_pairing == psms.peptide_pairing


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
    assert psms.spectra.columns == ["scan"]
    assert set(psms.score_columns) == {"combined p-value", "x"}
    assert psms.peptides.name == "sequence"
    assert psms.targets.shape == (10,)
    assert psms.targets.sum() == 6
    assert (~psms.targets).sum() == 4


def test_read_mztab(real_mztab):
    with pytest.raises(ValueError):
        psms = crema.read_mztab(real_mztab)
