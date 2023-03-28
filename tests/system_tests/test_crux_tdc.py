"""These tests verify that confidence estimates are correctly produced via tdc through various parsers"""
import unittest

import numpy as np
import pandas as pd
from crema.parsers.tide import read_tide


def test_tide_tdc(target_tide_txt, decoy_tide_txt, tmp_path):
    expected_peptide_pairing = {
        "APPLE": "APLPE",
        "BANANA": "ANANAB",
        "CHERRY": "CHRREY",
        "DURIAN": "DRIUAN",
        "EGGPLANT": "EGPGLATN",
    }
    expected_target_psms = pd.DataFrame(
        [
            ["f1", 3, "CHERRY", "p3", 0.1, 1 / 5],
            ["f1", 7, "BANANA", "p7", 0.2, 1 / 5],
            ["f1", 5, "EGGPLANT", "p5", 0.25, 1 / 5],
            ["f1", 10, "EGGPLANT", "p10", 0.3, 1 / 5],
            ["f1", 2, "BANANA", "p2", 0.4, 1 / 5],
            ["f1", 1, "APPLE", "p1", 0.5, 2 / 9],
            ["f1", 4, "DURIAN", "p4", 0.55, 2 / 9],
            ["f1", 6, "APPLE", "p6", 0.60, 2 / 9],
            ["f1", 8, "CHERRY", "p8", 0.70, 2 / 9],
        ],
        columns=[
            "file",
            "scan",
            "sequence",
            "protein id",
            "combined p-value",
            "crema q-value",
        ],
    )
    expected_decoy_psms = pd.DataFrame(
        [
            ["f1", 9, "DRIUAN", "p9", 0.5],
        ],
        columns=[
            "file",
            "scan",
            "sequence",
            "protein id",
            "combined p-value",
        ],
    )
    expected_target_peptides = pd.DataFrame(
        [
            ["f1", 3, "CHERRY", "p3", 0.1, 1 / 3],
            ["f1", 7, "BANANA", "p7", 0.2, 1 / 3],
            ["f1", 5, "EGGPLANT", "p5", 0.25, 1 / 3],
            ["f1", 1, "APPLE", "p1", 0.5, 1 / 2],
        ],
        columns=[
            "file",
            "scan",
            "sequence",
            "protein id",
            "combined p-value",
            "crema q-value",
        ],
    )
    expected_decoy_peptides = pd.DataFrame(
        [
            ["f1", 9, "DRIUAN", "p9", 0.5],
        ],
        columns=[
            "file",
            "scan",
            "sequence",
            "protein id",
            "combined p-value",
        ],
    )

    psms = read_tide([target_tide_txt, decoy_tide_txt])
    conf = psms.assign_confidence(score_column="combined p-value", desc=False)

    unittest.TestCase().assertDictEqual(
        expected_peptide_pairing, psms.peptide_pairing
    )
    np.testing.assert_array_equal(
        expected_target_psms.values,
        conf.confidence_estimates["psms"].values,
    )
    np.testing.assert_array_equal(
        expected_decoy_psms.values,
        conf.decoy_confidence_estimates["psms"].values,
    )
    np.testing.assert_array_equal(
        expected_target_peptides.values,
        conf.confidence_estimates["peptides"].values,
    )
    np.testing.assert_array_equal(
        expected_decoy_peptides.values,
        conf.decoy_confidence_estimates["peptides"].values,
    )
