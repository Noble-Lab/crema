"""
These tests verify that our q-value calculations are correct.
"""

import pytest
import numpy as np

from crema.qvalues import tdc


# TDC -------------------------------------------------------------------------
@pytest.fixture
def desc_scores():
    """Create a series of descending scores and their q-values"""
    scores = np.array([10, 10, 9, 8, 7, 7, 6, 5, 4, 3, 2, 2, 1, 1, 1, 1])
    target = np.array([1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0])
    qvals = np.array(
        [
            1 / 4,
            1 / 4,
            1 / 4,
            1 / 4,
            2 / 6,
            2 / 6,
            2 / 6,
            3 / 7,
            3 / 7,
            4 / 7,
            5 / 8,
            5 / 8,
            1,
            1,
            1,
            1,
        ]
    )
    return scores, target, qvals


def test_tdc_descending(desc_scores):
    """Test that q-values are correct for descending scores"""
    scores, target, true_qvals = desc_scores
    dtypes = [np.float64, np.uint8, np.int8, np.float32]
    for dtype in dtypes:
        qvals = tdc(scores.astype(dtype), target, desc=True)
        np.testing.assert_array_equal(qvals, true_qvals)

        qvals = tdc(scores, target.astype(dtype), desc=True)
        np.testing.assert_array_equal(qvals, true_qvals)


def test_tdc_ascending(desc_scores):
    """Test that q-values are correct for ascending scores"""
    scores, target, true_qvals = desc_scores
    scores = -scores
    dtypes = [np.float64, np.uint8, np.int8, np.float32]
    for dtype in dtypes:
        qvals = tdc(scores.astype(dtype), target, desc=False)
        np.testing.assert_array_equal(qvals, true_qvals)

        qvals = tdc(scores, target.astype(dtype), desc=False)
        np.testing.assert_array_equal(qvals, true_qvals)


def test_tdc_non_bool():
    """If targets is not boolean, should get a value error"""
    scores = np.array([1, 2, 3, 4, 5])
    targets = np.array(["1", "0", "1", "0", "blarg"])
    with pytest.raises(ValueError):
        tdc(scores, targets)


def test_tdc_diff_len():
    """If the arrays are different lengths, should get a ValueError"""
    scores = np.array([1, 2, 3, 4, 5])
    targets = np.array([True] * 3 + [False] * 3)
    with pytest.raises(ValueError):
        tdc(scores, targets)



# TDC edge cases ------------------------------------------------------------------
def test_tdc_target_wins_both():
    """Single spectrum: target beats decoy → q-value reflects (0+1)/1 = 1.0."""
    scores = np.array([2.0, 1.0])
    target = np.array([True, False])
    qvals = tdc(scores, target, desc=True)
    # Only one target, one decoy; FDR = (0+1)/1 = 1.0 → q = 1.0
    assert qvals[0] == pytest.approx(1.0)


def test_tdc_decoy_wins():
    """Single spectrum: decoy beats target → first entry FDR undefined (no targets yet)."""
    scores = np.array([2.0, 1.0])
    target = np.array([False, True])
    qvals = tdc(scores, target, desc=True)
    # Decoy is first: cum_targets=0, FDR = out=1.0 (where clause); target: FDR=(1+1)/1=2.0→q=1.0
    assert qvals[0] == pytest.approx(1.0)
    assert qvals[1] == pytest.approx(1.0)


def test_tdc_four_targets_two_decoys():
    """Hand-computed case: 4 targets win their spectra, 2 decoys win theirs.

    After competition (desc=True), sorted scores are:
      0.95 T, 0.90 T, 0.85 T, 0.80 T, 0.75 D, 0.70 D

    FDR:  1.0,  0.5, 0.333, 0.25,  0.5, 0.75
    Q:   0.25, 0.25, 0.25, 0.25,  0.50, 0.75
    """
    scores = np.array([0.95, 0.90, 0.85, 0.80, 0.75, 0.70])
    target = np.array([True, True, True, True, False, False])
    qvals = tdc(scores, target, desc=True)
    expected = np.array([0.25, 0.25, 0.25, 0.25, 0.5, 0.75])
    np.testing.assert_allclose(qvals, expected, atol=1e-6)


@pytest.mark.parametrize(
    "scores,target",
    [
        # All targets, then all decoys
        (
            np.array([5.0, 4.0, 3.0, 2.0, 1.0]),
            np.array([True, True, False, False, False]),
        ),
        # Alternating
        (
            np.array([5.0, 4.0, 3.0, 2.0, 1.0]),
            np.array([True, False, True, False, True]),
        ),
        # Decoys first
        (
            np.array([5.0, 4.0, 3.0, 2.0, 1.0]),
            np.array([False, False, True, True, True]),
        ),
    ],
)
def test_tdc_qvalues_non_decreasing(scores, target):
    """Q-values must be non-decreasing as score decreases (desc=True)."""
    qvals = tdc(scores, target, desc=True)
    sorted_qvals = qvals[np.argsort(-scores)]
    assert all(
        sorted_qvals[i] <= sorted_qvals[i + 1]
        for i in range(len(sorted_qvals) - 1)
    ), f"Q-values not monotone: {sorted_qvals}"


def test_tdc_qvalues_in_unit_interval():
    """All q-values must lie in [0, 1]."""
    scores = np.array([10, 9, 8, 7, 6, 5, 4, 3, 2, 1], dtype=float)
    target = np.array([1, 1, 1, 1, 0, 1, 0, 1, 0, 0])
    qvals = tdc(scores, target, desc=True)
    assert (qvals >= 0).all()
    assert (qvals <= 1).all()


def test_fdr2qvalue_monotone_output():
    """_fdr2qvalue must return a non-increasing sequence (worst→best ordering)."""
    from crema.qvalues import _fdr2qvalue

    scores = np.array([1.0, 2.0, 3.0, 4.0, 5.0])  # worst to best
    fdr = np.array([0.9, 0.6, 0.3, 0.4, 0.5])  # non-monotone FDR
    qvals = _fdr2qvalue(scores, fdr)
    # Output should be non-increasing (monotone minimum)
    assert all(qvals[i] >= qvals[i + 1] for i in range(len(qvals) - 1))
