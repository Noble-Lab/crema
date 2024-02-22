"""
These tests verify that our q-value calculations are correct.
"""

import logging

import pytest
import numpy as np

from crema.qvalues import tdc, mixmax


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


# MixMax -------------------------------------------------------------------------
@pytest.fixture
def mixmax_scores():
    """Increasing score/target arrays with enough depth for mixmax to estimate _q_-values"""
    N = 32
    tgt = 10 + 2 * np.random.randn(N)
    dec = 7 + 2 * np.random.randn(N)

    return (
        np.concatenate([tgt, dec]),
        np.array([True] * len(tgt) + [False] * len(dec)),
    )


@pytest.fixture
def desc_scores_mixmax():
    """Create a series of descending scores and their q-values"""
    scores = np.sort([10, 10, 9, 8, 7, 7, 6, 5, 4, 3, 2, 2, 1, 1, 1, 1])
    target = np.array([1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0])
    qvals = np.array(
        [
            0,
            0,
            0,
            0,
            np.nan,
            0.0833333,
            0.0833333,
            np.nan,
            0.1122449,
            np.nan,
            0.14671053,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ]
    )
    return scores, target, qvals


def do_mixmax(scores, target, desc):
    """Compute the necessary inputs for `mixmax()`"""
    # qvalues.py::calculate_mixmax_qval expects target scores
    # and decoy scores to be sorted from worst to best
    # qvalues.py::mixmax expected combined_sorted scores
    # to be sorted from best to worst
    tgt = np.array(sorted(scores[target.astype(bool)], reverse=desc))
    dec = np.array(sorted(scores[~target.astype(bool)], reverse=desc))
    all_scores = sorted(zip(scores, target), reverse=not desc)

    # Target/decoy score arrays must always be increasing, so invert if necessary
    if desc:
        tgt = -1 * tgt
        dec = -1 * dec

    res = mixmax(
        tgt,
        dec,
        np.array([s for s, _ in all_scores]),
        np.array([t for _, t in all_scores]),
        # **kwargs,
    )

    # Q-values are returned sorted worst to best; reverse them to match best-to-worst order
    # of the combined dataset
    return res[0], res[1][::-1]


def test_mixmax_descending(desc_scores_mixmax):
    """Test that q-values are correct for descending scores"""
    scores, target, true_qvals = desc_scores_mixmax

    tgt_qvals = [q for q, t in zip(true_qvals, target) if t]

    dtypes = [np.float64, np.uint8, np.int8, np.float32]
    for dtype in dtypes:
        qvals = do_mixmax(scores.astype(dtype), target, desc=True)[1]
        assert np.allclose(qvals, tgt_qvals, atol=1e-5)

        qvals = do_mixmax(scores, target.astype(dtype), desc=True)[1]
        assert np.allclose(qvals, tgt_qvals, atol=1e-5)


def test_mixmax_ascending(desc_scores_mixmax):
    """Test that q-values are correct for ascending scores"""
    scores, target, true_qvals = desc_scores_mixmax

    tgt_qvals = [q for q, t in zip(true_qvals, target) if t]

    scores = -scores
    dtypes = [np.float64, np.uint8, np.int8, np.float32]
    for dtype in dtypes:
        qvals = do_mixmax(scores.astype(dtype), target, desc=False)[1]
        assert np.allclose(qvals, tgt_qvals, atol=1e-5)

        qvals = do_mixmax(scores, target.astype(dtype), desc=False)[1]
        assert np.allclose(qvals, tgt_qvals, atol=1e-5)


def test_mixmax_singular(desc_scores):
    """Test that q-values are correct when pi0 == 1.0"""
    scores, target, _ = desc_scores

    dtypes = [np.float64, np.uint8, np.int8, np.float32]
    for dtype in dtypes:
        pi0, qvals = do_mixmax(scores.astype(dtype), target, desc=True)
        assert pi0 == 1.0
        assert all(q == 1.0 for q in qvals)

        pi0, qvals = do_mixmax(scores, target.astype(dtype), desc=True)
        assert pi0 == 1.0
        assert all(q == 1.0 for q in qvals)


def test_mixmax_nonsingular(mixmax_scores):
    """Test that q-values can be computed for pi0 != 1.0"""
    scores, target = mixmax_scores

    dtypes = [np.float64, np.uint8, np.int8, np.float32]
    for dtype in dtypes:
        pi0, qvals = do_mixmax(scores.astype(dtype), target, desc=False)
        assert pi0 != 1.0
        assert any(q != 1.0 for q in qvals)

        pi0, qvals = do_mixmax(scores, target.astype(dtype), desc=False)
        assert pi0 != 1.0
        assert any(q != 1.0 for q in qvals)


def test_mixmax_empty():
    """Test that q-values can be computed for pi0 != 1.0"""
    scores, target = np.array([], dtype=float), np.array([], dtype=float)

    dtypes = [np.float64, np.uint8, np.int8, np.float32]
    for dtype in dtypes:
        pi0, qvals = do_mixmax(scores.astype(dtype), target, desc=False)
        assert pi0 == 1.0
        assert all(q == 1.0 for q in qvals)

        pi0, qvals = do_mixmax(scores, target.astype(dtype), desc=False)
        assert pi0 == 1.0
        assert all(q == 1.0 for q in qvals)
