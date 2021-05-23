"""Fixtures that are used in multiple tests"""
from pathlib import Path

import pytest
import pandas as pd


@pytest.fixture
def basic_crux_df():
    """A simple crux-like dataframe"""
    df = pd.DataFrame(
        [
            [1, 10, 0.7, "A", "target", 0.7],
            [2, 20, 0.4, "B", "decoy", 0.1],
            [3, 30, 0.1, "C", "target", 0.2],
            [4, 40, 0.55, "D", "target", 0.8],
            [5, 50, 0.25, "E", "target", 0.25],
            [1, 10, 0.6, "F", "target", 0.6],
            [2, 20, 0.2, "G", "decoy", 0.2],
            [3, 30, 0.7, "H", "decoy", 0.4],
            [4, 40, 0.56, "I", "target", 0.56],
            [5, 50, 0.3, "J", "decoy", 0.3],
        ],
        columns=[
            "scan",
            "spectrum precursor m/z",
            "combined p-value",
            "sequence",
            "target/decoy",
            "x",
        ],
    )
    return df


@pytest.fixture
def target_crux_df():
    """A crux-like dataframe of target psms"""
    df = pd.DataFrame(
        [
            [1, 10, 101, 0.5, "APPLE", "target"],
            [2, 20, 102, 0.4, "BANANA", "target"],
            [3, 30, 103, 0.1, "CHERRY", "target"],
            [4, 40, 104, 0.55, "DURIAN", "target"],
            [5, 50, 105, 0.25, "EGGPLANT", "target"],
            [6, 60, 106, 0.6, "APPLE", "target"],
            [7, 70, 107, 0.2, "BANANA", "target"],
            [8, 80, 108, 0.7, "CHERRY", "target"],
            [9, 90, 109, 0.56, "DURIAN", "target"],
            [10, 100, 110, 0.3, "EGGPLANT", "target"],
        ],
        columns=[
            "scan",
            "spectrum precursor m/z",
            "peptide mass",
            "combined p-value",
            "sequence",
            "target/decoy",
        ],
    )
    return df


@pytest.fixture
def decoy_crux_df():
    """A crux-like dataframe of decoy psms"""
    df = pd.DataFrame(
        [
            [1, 10, 101, 0.7, "APLPE", "decoy", "APPLE"],
            [2, 20, 102, 0.45, "ANANAB", "decoy", "BANANA"],
            [3, 30, 103, 0.3, "CHRREY", "decoy", "CHERRY"],
            [4, 40, 104, 0.6, "DRIUAN", "decoy", "DURIAN"],
            [5, 50, 105, 0.9, "EGPGLATN", "decoy", "EGGPLANT"],
            [6, 60, 106, 0.75, "APLPE", "decoy", "APPLE"],
            [7, 70, 107, 0.25, "ANANAB", "decoy", "BANANA"],
            [8, 80, 108, 0.9, "CHRREY", "decoy", "CHERRY"],
            [9, 90, 109, 0.5, "DRIUAN", "decoy", "DURIAN"],
            [10, 100, 110, 0.8, "EGPGLATN", "decoy", "EGGPLANT"],
        ],
        columns=[
            "scan",
            "spectrum precursor m/z",
            "peptide mass",
            "combined p-value",
            "sequence",
            "target/decoy",
            "original target sequence",
        ],
    )
    return df


@pytest.fixture
def mod_target_crux_df():
    """A crux-like dataframe of target psms"""
    df = pd.DataFrame(
        [
            [1, 10, 101, 0.1, "ALLSLR", "target"],
            [2, 20, 102, 0.2, "AILSIR", "target"],
            [3, 30, 103, 0.3, "LAVITR", "target"],
            [4, 40, 104, 0.4, "QTPPAR", "target"],
            [5, 50, 105, 0.5, "IPLVNL", "target"],
            [6, 60, 106, 0.6, "SRGPPR", "target"],
            [7, 70, 107, 0.7, "GEVPN[0.98]R", "target"],
            [8, 80, 108, 0.8, "GGHMDR", "target"],
            [9, 90, 109, 0.9, "NPANRT", "target"],
            [10, 100, 110, 0.99, "SADAGPR", "target"],
        ],
        columns=[
            "scan",
            "spectrum precursor m/z",
            "peptide mass",
            "combined p-value",
            "sequence",
            "target/decoy",
        ],
    )
    return df


@pytest.fixture
def mod_decoy_crux_df():
    """A crux-like dataframe of target psms with peptide modifications"""
    df = pd.DataFrame(
        [
            [1, 10, 101, 0.1, "ALLLSR", "decoy", "ALLSLR"],
            [2, 20, 102, 0.3, "ALLISR", "decoy", "ASLILR"],
            [3, 30, 103, 0.2, "LATVIR", "decoy", "LIVATR"],
            [4, 40, 104, 0.6, "QTAPPR", "decoy", "QTPPAR"],
            [5, 50, 105, 0.5, "IVLN[0.98]PL", "decoy", "IPLVNL"],
            [6, 60, 106, 0.9, "RGSPPR", "decoy", "RSGPPR"],
            [7, 70, 107, 0.4, "GN[0.98]PEVR", "decoy", "GEVPNR"],
            [8, 80, 108, 0.7, "GDMGHR", "decoy", "GGHMDR"],
            [9, 90, 109, 0.8, "N[0.98]NPTAR", "decoy", "NPANTR"],
            [10, 100, 110, 0.55, "SGDAPAR", "decoy", "SAADGPR"],
        ],
        columns=[
            "scan",
            "spectrum precursor m/z",
            "peptide mass",
            "combined p-value",
            "sequence",
            "target/decoy",
            "original target sequence",
        ],
    )
    return df


@pytest.fixture
def basic_crux_txt(basic_crux_df, tmp_path):
    """A simple crux-like txt file"""
    out_file = tmp_path / "crux.txt"
    basic_crux_df.to_csv(out_file, sep="\t", index=False)
    return out_file


@pytest.fixture
def target_crux_txt(target_crux_df, tmp_path):
    """A crux-like txt file of target psms"""
    out_file = tmp_path / "target_crux.txt"
    target_crux_df.to_csv(out_file, sep="\t", index=False)
    return out_file


@pytest.fixture
def decoy_crux_txt(decoy_crux_df, tmp_path):
    """A crux-like txt file of decoy psms"""
    out_file = tmp_path / "decoy_crux.txt"
    decoy_crux_df.to_csv(out_file, sep="\t", index=False)
    return out_file


@pytest.fixture
def mod_target_crux_txt(mod_target_crux_df, tmp_path):
    """A crux-like txt file of target psms"""
    out_file = tmp_path / "mod_target_crux.txt"
    mod_target_crux_df.to_csv(out_file, sep="\t", index=False)
    return out_file


@pytest.fixture
def mod_decoy_crux_txt(mod_decoy_crux_df, tmp_path):
    """A crux-like txt file of decoy psms"""
    out_file = tmp_path / "mod_decoy_crux.txt"
    mod_decoy_crux_df.to_csv(out_file, sep="\t", index=False)
    return out_file


@pytest.fixture
def basic_crux_csv(basic_crux_df, tmp_path):
    """A simple crux-like csv file"""
    out_file = tmp_path / "crux.csv"
    basic_crux_df.to_csv(out_file, sep=",", index=False)
    return out_file


@pytest.fixture
def real_crux_txt():
    """Return real crux txt files"""
    targets = Path("data/example_psms_target.txt")
    decoys = Path("data/example_psms_decoy.txt")
    return targets, decoys


@pytest.fixture
def real_mztab():
    """Return a real mzTab file"""
    return Path("data/MSV000085729.mzTab")
