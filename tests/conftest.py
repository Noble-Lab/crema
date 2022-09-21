"""Fixtures that are used in multiple tests"""
from pathlib import Path

import pytest
import pandas as pd


@pytest.fixture
def basic_crux_df():
    """A simple crux-like dataframe"""
    df = pd.DataFrame(
        [
            ["f1", 1, 10, 0.7, "APPLE", "target", 0.7, "p1", "APPLE"],
            ["f1", 2, 20, 0.4, "ANANAB", "decoy", 0.1, "p2", "BANANA"],
            ["f1", 3, 30, 0.1, "CHERRY", "target", 0.2, "p3", "CHERRY"],
            ["f1", 4, 40, 0.55, "DURIAN", "target", 0.8, "p4", "DURIAN"],
            ["f1", 5, 50, 0.25, "EGGPLANT", "target", 0.25, "p5", "EGGPLANT"],
            ["f1", 1, 10, 0.6, "FIG", "target", 0.6, "p6", "FIG"],
            ["f1", 2, 20, 0.2, "GARPE", "decoy", 0.2, "p7", "GRAPE"],
            ["f1", 3, 30, 0.7, "HEONDEYW", "decoy", 0.4, "p8", "HONEYDEW"],
            ["f1", 4, 40, 0.56, "ICE", "target", 0.56, "p9", "ICE"],
            ["f1", 5, 50, 0.3, "JMA", "decoy", 0.3, "p10", "JAM"],
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
            ["f1", 1, 10, 101, 0.5, "APPLE", "target", "p1"],
            ["f1", 2, 20, 102, 0.4, "BANANA", "target", "p2"],
            ["f1", 3, 30, 103, 0.1, "CHERRY", "target", "p3"],
            ["f1", 4, 40, 104, 0.55, "DURIAN", "target", "p4"],
            ["f1", 5, 50, 105, 0.25, "EGGPLANT", "target", "p5"],
            ["f1", 6, 60, 106, 0.6, "APPLE", "target", "p6"],
            ["f1", 7, 70, 107, 0.2, "BANANA", "target", "p7"],
            ["f1", 8, 80, 108, 0.7, "CHERRY", "target", "p8"],
            ["f1", 9, 90, 109, 0.56, "DURIAN", "target", "p9"],
            ["f1", 10, 100, 110, 0.3, "EGGPLANT", "target", "p10"],
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
        ],
    )
    return df


@pytest.fixture
def decoy_crux_df():
    """A crux-like dataframe of decoy psms"""
    df = pd.DataFrame(
        [
            ["f1", 1, 10, 101, 0.7, "APLPE", "decoy", "APPLE", "p1"],
            ["f1", 2, 20, 102, 0.45, "ANANAB", "decoy", "BANANA", "p2"],
            ["f1", 3, 30, 103, 0.3, "CHRREY", "decoy", "CHERRY", "p3"],
            ["f1", 4, 40, 104, 0.6, "DRIUAN", "decoy", "DURIAN", "p4"],
            ["f1", 5, 50, 105, 0.9, "EGPGLATN", "decoy", "EGGPLANT", "p5"],
            ["f1", 6, 60, 106, 0.75, "APLPE", "decoy", "APPLE", "p6"],
            ["f1", 7, 70, 107, 0.25, "ANANAB", "decoy", "BANANA", "p7"],
            ["f1", 8, 80, 108, 0.9, "CHRREY", "decoy", "CHERRY", "p8"],
            ["f1", 9, 90, 109, 0.5, "DRIUAN", "decoy", "DURIAN", "p9"],
            ["f1", 10, 100, 110, 0.8, "EGPGLATN", "decoy", "EGGPLANT", "p10"],
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
            ["f1", 1, 10, 101, 0.1, "ALLSLR", "target", "p1"],
            ["f1", 2, 20, 102, 0.2, "AILSIR", "target", "p2"],
            ["f1", 3, 30, 103, 0.3, "LAVITR", "target", "p3"],
            ["f1", 4, 40, 104, 0.4, "QTPPAR", "target", "p4"],
            ["f1", 5, 50, 105, 0.5, "IPLVNL", "target", "p5"],
            ["f1", 6, 60, 106, 0.6, "SRGPPR", "target", "p6"],
            ["f1", 7, 70, 107, 0.7, "GEVPN[0.98]R", "target", "p7"],
            ["f1", 8, 80, 108, 0.8, "GGHMDR", "target", "p8"],
            ["f1", 9, 90, 109, 0.9, "NPANRT", "target", "p9"],
            ["f1", 10, 100, 110, 0.99, "SADAGPR", "target", "p10"],
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
        ],
    )
    return df


@pytest.fixture
def mod_decoy_crux_df():
    """A crux-like dataframe of target psms with peptide modifications"""
    df = pd.DataFrame(
        [
            ["f1", 1, 10, 101, 0.1, "ALLLSR", "decoy", "ALLSLR", "p1"],
            ["f1", 2, 20, 102, 0.3, "ALLISR", "decoy", "ASLILR", "p2"],
            ["f1", 3, 30, 103, 0.2, "LATVIR", "decoy", "LIVATR", "p3"],
            ["f1", 4, 40, 104, 0.6, "QTAPPR", "decoy", "QTPPAR", "p4"],
            ["f1", 5, 50, 105, 0.5, "IVLN[0.98]PL", "decoy", "IPLVNL", "p5"],
            ["f1", 6, 60, 106, 0.9, "RGSPPR", "decoy", "RSGPPR", "p6"],
            ["f1", 7, 70, 107, 0.4, "GN[0.98]PEVR", "decoy", "GEVPNR", "p7"],
            ["f1", 8, 80, 108, 0.7, "GDMGHR", "decoy", "GGHMDR", "p8"],
            ["f1", 9, 90, 109, 0.8, "N[0.98]NPTAR", "decoy", "NPANTR", "p9"],
            ["f1", 10, 100, 110, 0.55, "SGDAPAR", "decoy", "SAADGPR", "p10"],
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
def basic_comet_df():
    """A simple comet-like dataframe"""
    df = pd.DataFrame(
        [
            [1, 10, 0.7, "K.A[15.9949]PPLE.A", "prot1"],
            [2, 20, 0.4, "K.BNANAA.A", "DECOY_prot2"],
            [3, 30, 0.1, "K.CH[15.9949]E[15.9949]RRY.A", "prot3"],
            [4, 40, 0.55, "K.DURIAN.A", "prot4"],
            [5, 50, 0.25, "K.EGGPLANT.A", "prot5"],
            [1, 10, 0.6, "K.FIGS.A", "prot6"],
            [2, 20, 0.2, "K.GPARE.A", "DECOY_prot7"],
            [3, 30, 0.7, "K.H[15.9949]E[15.9949]DYENOW.A", "DECOY_prot8"],
            [4, 40, 0.56, "K.IC[15.9949]ES.A", "prot9"],
            [5, 50, 0.3, "K.JMAS.A", "DECOY_prot10"],
        ],
        columns=[
            "scan",
            "exp_neutral_mass",
            "e-value",
            "modified_peptide",
            "protein",
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


@pytest.fixture
def mod_comet_txt(basic_comet_df, tmp_path):
    """A crux-like txt file of target psms"""
    out_file = tmp_path / "mod_comet.txt"
    basic_comet_df.to_csv(out_file, sep="\t", index=False)
    return out_file
