"""
This module estimates q-values.
"""

import logging

import numpy as np
import numba as nb

LOGGER = logging.getLogger(__name__)


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


def mixmax(target_scores, decoy_scores, combined_score, combined_score_target):
    """
    Estimate q-values using mix-max competition.

    TODO
    Estimates q-values using...

    Parameters
    ----------
    target_scores : numpy.ndarray of float
       An array of the best target PSM score per spectrum.
    decoy_scores : numpy.ndarray of float
       An array of the best decoy PSM score per spectrum.
    combined_score : pandas.DataFrame
       A :py:class:`pandas.DataFram` of PSMs. Dataframe includes
       best ranked target and decoy per spectrum.
    combined_score_target: numpy.ndarray of float
       The target/deoy column from the combined_score dataframe

    Returns
    ----------
    pi0 : float
        Estimated pi_zero.
    fdrmod : numpy.ndarray
        Array of calculated q-values.
    """
    # TODO try except and some error checking

    num_targets = target_scores.shape[0]
    num_decoys = decoy_scores.shape[0]

    # calculate p-values from scores
    n_decoys = 1
    pos_same = 0
    neg_same = 0
    pval_list = []
    cur_score = None
    for score, target in zip(combined_score, combined_score_target):
        if target:
            pos_same += 1
        else:
            neg_same += 1

        # TODO Unclear if need if statement. This appears in Percolator code
        # however, adding it makes Crema and Percolator pValList different length
        # if curScore != score:
        for ix in range(0, pos_same):
            pval_list.append(n_decoys + (neg_same * (ix + 1)) / (pos_same + 1))

        n_decoys += neg_same
        neg_same = 0
        pos_same = 0
        cur_score = score
    pval_list = np.array(pval_list) / n_decoys

    # calculate pi0
    if len(pval_list) > 0:
        pi0 = estimate_pi0(pval_list)
    else:
        # Corner case: if pval_list is empty there are no targets.
        # In this case pi0 is undefined, but we set it to 1.0, as all
        # provided targets are incorrect, and we avoid zeroes.
        pi0 = 1.0

    if pi0 == 1.0:
        # All targets are assumed to be incorrect! Algorithm 1, line 5-6
        LOGGER.debug("FALLBACK: pi0==1.0; all q-values will be 1.0")
        fdrmod = np.full(num_targets, 1.0)  # all q-values are 1
    elif pi0 < 0 or pi0 >= 1 or not np.isfinite(pi0):
        raise ValueError(
            f"Invalid pi0 estimate ({pi0}); unable to proceed FDR estimation!"
        )
    else:
        fdrmod = calculate_mixmax_qval(
            np.array(target_scores), np.array(decoy_scores), pi0
        )

    return (pi0, fdrmod)


@nb.njit
def estimate_pi0(pval_list):
    """
    Estimates pi0. Add description. TODO

    Parameters
    -------
    pval_list : numpy.ndarray of float
        A list of p-values sorted from smallest to largest.

    Returns
    -------
    pi0 : float
        Estimated pi_zero.
    """
    num_lambda = 100
    max_lambda = 0.5
    num_boot = 100

    # LOGGER.debug("pval_list=%s", pval_list)

    n_pval = pval_list.size
    lambda_list = []
    pi0s_list = []
    for idx in range(0, num_lambda):
        cur_lambda = ((idx + 1) / num_lambda) * max_lambda

        # Find the index of the first element in pval_list
        # that is > lambda.
        start = np.searchsorted(pval_list, cur_lambda)
        W1 = n_pval - start
        pi0 = W1 / n_pval / (1.0 - cur_lambda)

        if pi0 > 0.0:
            lambda_list.append(cur_lambda)
            pi0s_list.append(pi0)

    # LOGGER.debug("%s", lambda_list)
    # LOGGER.debug("%s", pi0s_list)

    assert len(pi0s_list) != 0, (
        "Error in the input data: "
        "too good separation between target and decoy PSMs."
    )

    minPi0 = min(pi0s_list)

    mse_list = np.zeros(len(pi0s_list))
    max_size = 1000
    # Examine which lambda level that is most stable under bootstrap
    for i in range(0, num_boot):
        # Create an array of bootstrapped p-values, and sort in ascending order.
        num_draw = min(n_pval, max_size)
        pBoot_list = np.random.choice(pval_list, size=num_draw, replace=True)
        pBoot_list.sort()

        for idx in range(0, len(lambda_list)):
            start = np.searchsorted(pBoot_list, lambda_list[idx])
            W1 = num_draw - start
            pi0Boot = W1 / num_draw / (1.0 - lambda_list[idx])

            # Estimated mean-squared error
            mse_list[idx] += (pi0Boot - minPi0) * (pi0Boot - minPi0)

    # Which index did the iterator go?
    minIdx = np.argmin(mse_list)

    # LOGGER.debug(f"Estimated pi0=%f at lambda=%f (MSE=%f)", pi0s_list[minIdx], lambda_list[minIdx], mse_list[minIdx])

    pi0 = max(min(pi0s_list[minIdx], 1.0), 0.0)
    return pi0


@nb.njit
def calculate_mixmax_qval(target_scores, decoy_scores, pi0):
    """
    Estimate q-values using mix-max.
    """
    # Note that the notation in this function follows the notation found in
    # Percolator, which itself follows the notation found in Supplementary Note
    # 3 (Keich et al., JPR. 2015.). This supplment can be found at
    # http://dx.doi.org/10.1021/acs.jproteome.5b00081.

    # assert pi0 >= 0 and pi0 < 1

    num_targets = target_scores.shape[0]
    num_decoys = decoy_scores.shape[0]

    h_w_le_z = np.zeros(num_decoys)  # histogram for N_{w<=z}
    h_z_le_z = np.zeros(num_decoys)  # histogram for N_{z<=z}

    for j in range(0, num_decoys):
        h_w_le_z[j] = np.searchsorted(
            target_scores, decoy_scores[j], side="right"
        )
        h_z_le_z[j] = np.searchsorted(
            decoy_scores, decoy_scores[j], side="right"
        )

    fdrmod = np.zeros(num_targets)
    estPx_lt_zj = 0.0
    E_f1_mod_run_tot = 0.0
    j = num_decoys - 1
    n_z_ge_w = 0
    n_w_ge_w = 0
    prev_fdr = -1

    for i in range(num_targets - 1, -1, -1):
        while j >= 0 and decoy_scores[j] >= target_scores[i]:
            cnt_w = h_w_le_z[j]
            cnt_z = h_z_le_z[j]

            estPx_lt_zj = (cnt_w - pi0 * cnt_z) / ((1.0 - pi0) * cnt_z)
            if estPx_lt_zj > 1:
                estPx_lt_zj = 1.0
            elif estPx_lt_zj < 0:
                estPx_lt_zj = 0.0

            E_f1_mod_run_tot += estPx_lt_zj * (1.0 - pi0)
            n_z_ge_w += 1
            j -= 1

        # not sure if this faster or if original while loop is faster
        # AssignConfidence.cpp:1225
        n_w_ge_w = (target_scores >= target_scores[i]).sum()
        fdr = (n_z_ge_w * pi0 + E_f1_mod_run_tot) / n_w_ge_w

        if fdr > 1:
            fdr = 1.0
        fdrmod[i] = fdr

    # convert qvalues to fdr
    return _fdr2qvalue(target_scores, fdrmod)
