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
def basic_crux_txt(basic_crux_df, tmp_path):
    """A simple crux-like txt file"""
    out_file = tmp_path / "crux.txt"
    basic_crux_df.to_csv(out_file, sep="\t", index=False)
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
