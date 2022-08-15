"""Fixtures that are used in multiple tests"""
from pathlib import Path

import pytest
import pandas as pd


@pytest.fixture
def basic_crux_df():
    """A simple crux-like dataframe"""
    df = pd.DataFrame(
        [
            ["file1", 1, 10, 0.7, "APPLE", "target", 0.7, "prot1", "APPLE"],
            ["file1", 2, 20, 0.4, "ANANAB", "decoy", 0.1, "prot2", "BANANA"],
            ["file1", 3, 30, 0.1, "CHERRY", "target", 0.2, "prot3", "CHERRY"],
            ["file1", 4, 40, 0.55, "DURIAN", "target", 0.8, "prot4", "DURIAN"],
            [
                "file1",
                5,
                50,
                0.25,
                "EGGPLANT",
                "target",
                0.25,
                "prot5",
                "EGGPLANT",
            ],
            ["file1", 1, 10, 0.6, "FIG", "target", 0.6, "prot6", "FIG"],
            ["file1", 2, 20, 0.2, "GARPE", "decoy", 0.2, "prot7", "GRAPE"],
            [
                "file1",
                3,
                30,
                0.7,
                "HEONDEYW",
                "decoy",
                0.4,
                "prot8",
                "HONEYDEW",
            ],
            ["file1", 4, 40, 0.56, "ICE", "target", 0.56, "prot9", "ICE"],
            ["file1", 5, 50, 0.3, "JMA", "decoy", 0.3, "prot10", "JAM"],
        ],
        columns=[
            "file",
            "scan",
            "spectrum precursor m/z",
            "combined p-value",
            "sequence",
            "target/decoy",
            "x",
            "protein id",
            "original target sequence",
        ],
    )
    return df


@pytest.fixture
def target_crux_df():
    """A crux-like dataframe of target psms"""
    df = pd.DataFrame(
        [
            ["file1", 1, 10, 101, 0.5, "APPLE", "target", "prot1", "APPLE"],
            ["file1", 2, 20, 102, 0.4, "BANANA", "target", "prot2", "BANANA"],
            ["file1", 3, 30, 103, 0.1, "CHERRY", "target", "prot3", "CHERRY"],
            ["file1", 4, 40, 104, 0.55, "DURIAN", "target", "prot4", "DURIAN"],
            [
                "file1",
                5,
                50,
                105,
                0.25,
                "EGGPLANT",
                "target",
                "prot5",
                "EGGPLANT",
            ],
            ["file1", 6, 60, 106, 0.6, "APPLE", "target", "prot6", "APPLE"],
            ["file1", 7, 70, 107, 0.2, "BANANA", "target", "prot7", "BANANA"],
            ["file1", 8, 80, 108, 0.7, "CHERRY", "target", "prot8", "CHERRY"],
            ["file1", 9, 90, 109, 0.56, "DURIAN", "target", "prot9", "DURIAN"],
            [
                "file1",
                10,
                100,
                110,
                0.3,
                "EGGPLANT",
                "target",
                "prot10",
                "EGGPLANT",
            ],
        ],
        columns=[
            "file",
            "scan",
            "spectrum precursor m/z",
            "peptide mass",
            "combined p-value",
            "sequence",
            "target/decoy",
            "protein id",
            "original target sequence",
        ],
    )
    return df


@pytest.fixture
def decoy_crux_df():
    """A crux-like dataframe of decoy psms"""
    df = pd.DataFrame(
        [
            ["file1", 1, 10, 101, 0.7, "APLPE", "decoy", "APPLE", "prot1"],
            ["file1", 2, 20, 102, 0.45, "ANANAB", "decoy", "BANANA", "prot2"],
            ["file1", 3, 30, 103, 0.3, "CHRREY", "decoy", "CHERRY", "prot3"],
            ["file1", 4, 40, 104, 0.6, "DRIUAN", "decoy", "DURIAN", "prot4"],
            [
                "file1",
                5,
                50,
                105,
                0.9,
                "EGPGLATN",
                "decoy",
                "EGGPLANT",
                "prot5",
            ],
            ["file1", 6, 60, 106, 0.75, "APLPE", "decoy", "APPLE", "prot6"],
            ["file1", 7, 70, 107, 0.25, "ANANAB", "decoy", "BANANA", "prot7"],
            ["file1", 8, 80, 108, 0.9, "CHRREY", "decoy", "CHERRY", "prot8"],
            ["file1", 9, 90, 109, 0.5, "DRIUAN", "decoy", "DURIAN", "prot9"],
            [
                "file1",
                10,
                100,
                110,
                0.8,
                "EGPGLATN",
                "decoy",
                "EGGPLANT",
                "prot10",
            ],
        ],
        columns=[
            "file",
            "scan",
            "spectrum precursor m/z",
            "peptide mass",
            "combined p-value",
            "sequence",
            "target/decoy",
            "original target sequence",
            "protein id",
        ],
    )
    return df


@pytest.fixture
def mod_target_crux_df():
    """A crux-like dataframe of target psms"""
    df = pd.DataFrame(
        [
            ["file1", 1, 10, 101, 0.1, "ALLSLR", "target", "prot1", "ALLSLR"],
            ["file1", 2, 20, 102, 0.2, "AILSIR", "target", "prot2", "AILSIR"],
            ["file1", 3, 30, 103, 0.3, "LAVITR", "target", "prot3", "LAVITR"],
            ["file1", 4, 40, 104, 0.4, "QTPPAR", "target", "prot4", "QTPPAR"],
            ["file1", 5, 50, 105, 0.5, "IPLVNL", "target", "prot5", "IPLVNL"],
            ["file1", 6, 60, 106, 0.6, "SRGPPR", "target", "prot6", "SRGPPR"],
            [
                "file1",
                7,
                70,
                107,
                0.7,
                "GEVPN[0.98]R",
                "target",
                "prot7",
                "GEVPN[0.98]R",
            ],
            ["file1", 8, 80, 108, 0.8, "GGHMDR", "target", "prot8", "GGHMDR"],
            ["file1", 9, 90, 109, 0.9, "NPANRT", "target", "prot9", "NPANRT"],
            [
                "file1",
                10,
                100,
                110,
                0.99,
                "SADAGPR",
                "target",
                "prot10",
                "SADAGPR",
            ],
        ],
        columns=[
            "file",
            "scan",
            "spectrum precursor m/z",
            "peptide mass",
            "combined p-value",
            "sequence",
            "target/decoy",
            "protein id",
            "original target sequence",
        ],
    )
    return df


@pytest.fixture
def mod_decoy_crux_df():
    """A crux-like dataframe of target psms with peptide modifications"""
    df = pd.DataFrame(
        [
            ["file1", 1, 10, 101, 0.1, "ALLLSR", "decoy", "ALLSLR", "prot1"],
            ["file1", 2, 20, 102, 0.3, "ALLISR", "decoy", "ASLILR", "prot2"],
            ["file1", 3, 30, 103, 0.2, "LATVIR", "decoy", "LIVATR", "prot3"],
            ["file1", 4, 40, 104, 0.6, "QTAPPR", "decoy", "QTPPAR", "prot4"],
            [
                "file1",
                5,
                50,
                105,
                0.5,
                "IVLN[0.98]PL",
                "decoy",
                "IPLVNL",
                "prot5",
            ],
            ["file1", 6, 60, 106, 0.9, "RGSPPR", "decoy", "RSGPPR", "prot6"],
            [
                "file1",
                7,
                70,
                107,
                0.4,
                "GN[0.98]PEVR",
                "decoy",
                "GEVPNR",
                "prot7",
            ],
            ["file1", 8, 80, 108, 0.7, "GDMGHR", "decoy", "GGHMDR", "prot8"],
            [
                "file1",
                9,
                90,
                109,
                0.8,
                "N[0.98]NPTAR",
                "decoy",
                "NPANTR",
                "prot9",
            ],
            [
                "file1",
                10,
                100,
                110,
                0.55,
                "SGDAPAR",
                "decoy",
                "SAADGPR",
                "prot10",
            ],
        ],
        columns=[
            "file",
            "scan",
            "spectrum precursor m/z",
            "peptide mass",
            "combined p-value",
            "sequence",
            "target/decoy",
            "original target sequence",
            "protein id",
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


@pytest.fixture
def real_pepxml():
    """Return a real pepXML file"""
    return Path("data/tide-search.pep.xml")
