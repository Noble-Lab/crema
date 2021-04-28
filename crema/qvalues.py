"""
This module estimates q-values.
"""
import numpy as np
import numba as nb


def tdc(scores, target, desc=True):
    """
    Estimate q-values using target decoy competition.

    Estimates q-values using the simple target decoy competition method.
    For set of target and decoy PSMs meeting a specified score threshold,
    the false discovery rate (FDR) is estimated as:

    ...math:
        FDR = \frac{Decoys + 1}{Targets}

    More formally, let the scores of target and decoy PSMs be indicated as
    :math:`f_1, f_2, ..., f_{m_f}` and :math:`d_1, d_2, ..., d_{m_d}`,
    respectively. For a score threshold :math:`t`, the false discovery
    rate is estimated as:

    ...math:
        E\\{FDR(t)\\} = \frac{|\\{d_i > t; i=1, ..., m_d\\}| + 1}
        {\\{|f_i > t; i=1, ..., m_f|\\}}

    The reported q-value for each PSM is the minimum FDR at which that
    PSM would be accepted.

    Parameters
    ----------
    scores : numpy.ndarray of float
        A 1D array containing the score to rank by

    target : numpy.ndarray of bool
        A 1D array indicating if the entry is from a target or decoy
        hit. This should be boolean, where `True` indicates a target
        and `False` indicates a decoy. `target[i]` is the label for
        `metric[i]`; thus `target` and `metric` should be of
        equal length.

    desc : bool
        Are higher scores better? `True` indicates that they are,
        `False` indicates that they are not.

    Returns
    -------
    numpy.ndarray
        A 1D array with the estimated q-value for each entry. The
        array is the same length as the `scores` and `target` arrays.
    """
    scores = np.array(scores)

    try:
        target = np.array(target, dtype=bool)
    except ValueError:
        raise ValueError("'target' should be boolean.")

    if scores.shape[0] != target.shape[0]:
        raise ValueError("'scores' and 'target' must be the same length")

    # Unsigned integers can cause weird things to happen.
    # Convert all scores to floats to for safety.
    if np.issubdtype(scores.dtype, np.integer):
        scores = scores.astype(np.float_)

    # Sort and estimate FDR
    if desc:
        srt_idx = np.argsort(-scores)
    else:
        srt_idx = np.argsort(scores)

    scores = scores[srt_idx]
    target = target[srt_idx]
    cum_targets = target.cumsum()
    cum_decoys = (~target).cumsum()

    # Handles zeros in denominator
    fdr = np.divide(
        (cum_decoys + 1),
        cum_targets,
        out=np.ones_like(cum_targets, dtype=float),
        where=(cum_targets != 0),
    )

    # Some arrays need to be flipped so that we can loop through from
    # worse to best score.
    fdr = np.flip(fdr)
    scores = np.flip(scores)
    qvals = _fdr2qvalue(scores, fdr)
    qvals = np.flip(qvals)
    qvals = qvals[np.argsort(srt_idx)]

    return qvals


@nb.njit
def _fdr2qvalue(scores, fdr):
    """Quickly calculate q-values.

    Note that the if block logic in this function is to account for multiple
    tied scores. This function assumes that scores and fdr are sorted such
    that they are ordered from worst to best score. Additionally, for a series
    of tied scores, we assume that the first entry in these sorted arrays
    accounts for all of the PSMs that were tied.

    Parameters
    ----------
    scores : np.array
        The scores of the PSMs, sorted from worst to best.
    fdr : np.array
        The calculated FDR of the PSMs, sorted to match scores. Note that
        this is not the same thing as sorting the fdr array itself!

    Returns
    -------
    np.ndarray
        An array of q-values, in the same order as the input scores.
    """
    min_q = 1
    qvals = np.ones(len(fdr))
    start = 0
    for idx in range(len(qvals)):
        if idx < len(qvals) - 1 and scores[idx + 1] == scores[start]:
            continue

        if fdr[start] < min_q:
            min_q = fdr[start]

        qvals[start : idx + 1] = min_q
        start = idx + 1

    return qvals
