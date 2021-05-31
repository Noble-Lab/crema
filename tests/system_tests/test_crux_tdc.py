"""These tests verify that confidence estimates are correctly produced via tdc through various parsers"""
import unittest

import numpy as np
import pandas as pd
from crema.parsers.crux import read_crux


def test_crux_tdc(target_crux_txt, decoy_crux_txt, tmp_path):
    expected_peptide_pairing = {
        "APPLE": "APLPE",
        "BANANA": "ANANAB",
        "CHERRY": "CHRREY",
        "DURIAN": "DRIUAN",
        "EGGPLANT": "EGPGLATN",
    }
    expected_target_psms = pd.DataFrame(
        [
            [3, 30, "CHERRY", 0.1, 1 / 5],
            [7, 70, "BANANA", 0.2, 1 / 5],
            [5, 50, "EGGPLANT", 0.25, 1 / 5],
            [10, 100, "EGGPLANT", 0.3, 1 / 5],
            [2, 20, "BANANA", 0.4, 1 / 5],
            [1, 10, "APPLE", 0.5, 2 / 9],
            [4, 40, "DURIAN", 0.55, 2 / 9],
            [6, 60, "APPLE", 0.60, 2 / 9],
            [8, 80, "CHERRY", 0.70, 2 / 9],
        ],
        columns=[
            "scan",
            "spectrum precursor m/z",
            "sequence",
            "combined p-value",
            "crema q-value",
        ],
    )
    expected_decoy_psms = pd.DataFrame(
        [
            [9, 90, "DRIUAN", 0.5, 2 / 9],
        ],
        columns=[
            "scan",
            "spectrum precursor m/z",
            "sequence",
            "combined p-value",
            "crema q-value",
        ],
    )
    expected_target_peptides = pd.DataFrame(
        [
            [3, 30, "CHERRY", 0.1, 1 / 3],
            [7, 70, "BANANA", 0.2, 1 / 3],
            [5, 50, "EGGPLANT", 0.25, 1 / 3],
            [1, 10, "APPLE", 0.5, 1 / 2],
        ],
        columns=[
            "scan",
            "spectrum precursor m/z",
            "sequence",
            "combined p-value",
            "crema q-value",
        ],
    )
    expected_decoy_peptides = pd.DataFrame(
        [
            [9, 90, "DRIUAN", 0.5, 2 / 9],
        ],
        columns=[
            "scan",
            "spectrum precursor m/z",
            "sequence",
            "combined p-value",
            "crema q-value",
        ],
    )

    psms = read_crux([target_crux_txt, decoy_crux_txt], peptide_tdc=True)
    conf = psms.assign_confidence(score_column="combined p-value", desc=False)

    unittest.TestCase().assertDictEqual(
        expected_peptide_pairing, psms.peptide_pairing
    )
    np.array_equal(
        expected_target_psms.values,
        conf.confidence_estimates["psms"].values,
    )
    np.array_equal(
        expected_decoy_psms.values,
        conf.decoy_confidence_estimates["psms"].values,
    )
    np.array_equal(
        expected_target_peptides.values,
        conf.confidence_estimates["peptides"].values,
    )
    np.array_equal(
        expected_decoy_peptides.values,
        conf.decoy_confidence_estimates["peptides"].values,
    )
