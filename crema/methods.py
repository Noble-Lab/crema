"""
This module contains the various methods for calculating FDRs and q-values
"""

import random
from .result import Result


def calculate_tdc(psm):
    """
    Calculates FDR and Q-Value using Target-Decoy Competition

    Parameters
    ----------
    psm : A :py:class:`~crema.dataset.PsmDataset` object
        dataset to apply FDR calculation upon

    Returns
    -------
    data : pandas.DataFrame
        A pandas.DataFrame of the data from the original dataset with additional FDR and Q-Value columns
    """

    # copy psm dataframe from dataset object - don't want to manipulate the original
    data = psm.data.copy()

    # note down all column names from dataset object
    spectrum_col = psm.spectrum_col
    score_col = psm.score_col
    target_col = psm.target_col

    # determine items to sort by
    sort_order = []
    if type(spectrum_col) == str:
        sort_order = [spectrum_col, score_col, target_col]
    else:
        for col in spectrum_col:
            sort_order.append(col)
        sort_order.append(score_col)
        sort_order.append(target_col)

    # sort dataframe by spectrum and p-value ascending
    data = data.sort_values(by=sort_order)

    # look through and delete all duplicate psms with higher p-values
    data = _delete_duplicates(data, spectrum_col, score_col, target_col)

    # sort dataframe by spectrum and p-value ascending
    data = data.sort_values(by=[score_col])

    # calculate fdr at each psm
    _calculate_fdr(data, target_col)

    # calculate q value at each psm
    _calculate_q_value(data, "FDR")

    # return a result object containing the manipulated data and respective calculations
    return Result(data, spectrum_col, score_col, target_col)


def _compare_spectrum_col(spectrum_col, curr_row, next_row):
    if type(spectrum_col) == str:
        return curr_row[spectrum_col] != next_row[spectrum_col]
    else:
        for col in spectrum_col:
            if curr_row[col] != next_row[col]:
                return True
        return False


def _delete_duplicates(data, spectrum_col, score_col, target_col):
    """
    Deletes duplicate PSMs based on unique spectrum identifier

    Parameters
    ----------
    data : pandas.DataFrame
        dataframe of PSMs of data from original dataset object (MUST BE SORTED BY SPECTRUM, SCORE, TARGET)
    spectrum_col : str
        name of the column that identifies the psm
    score_col : str
        name of the column that defines the scores (p-values) of the psms
    target_col : str
        name of the column that indicates if a psm is a target/decoy

    Returns
    -------
    data : pandas.DataFrame
        A pandas.DataFrame of the data from the original dataset with additional FDR and Q-Value columns
    """

    prev_score = 1
    prev_target = True
    prev_index = data.iloc[-1].name
    for index, row in data.iterrows():
        # If spectrum col in current row is different from next row, increment curr variables and move on
        if _compare_spectrum_col(
            spectrum_col, data.loc[prev_index], data.loc[index]
        ):
            prev_score = row[score_col]
            prev_target = row[target_col]
            prev_index = index
        # If they are the same, check if the score/target values are the same and remove accordingly
        else:
            if row[score_col] == prev_score:
                if row[target_col] == prev_target:
                    data = data.drop(index)
                else:
                    prev_target = not prev_target
                    decision = random.randint(1, 2)
                    if decision == 1:
                        data.loc[prev_index, target_col] = True
                    else:
                        data.loc[index, target_col] = False
                    data = data.drop(index)
            else:
                data = data.drop(index)
    return data


def _calculate_fdr(data, target_col):
    """
    Calculates FDR at each PSM

    Parameters
    ----------
    data : pandas.DataFrame
        dataframe of PSMs of data from original dataset object (MUST BE SORTED BY SCORE ASCENDING)
    target_col : str
        name of the column that indicates if a psm is a target/decoy
    """
    # List to store fdr at each psm
    fdr = []
    # Keep track of # of targets/decoys along the way
    target = 0
    decoy = 0
    for value in data[target_col]:
        if value:
            target += 1
        else:
            decoy += 1
        fdr.append(min((decoy + 1) / target, 1))
    data["FDR"] = fdr


def _calculate_q_value(data, fdr_col):
    """
    Calculates Q_Value at each PSM

    Parameters
    ----------
    data : pandas.DataFrame
        dataframe of PSMs of data from original dataset object (MUST BE SORTED BY SCORE ASCENDING)
    fdr_col : str
        name of the column that indicates the FDR of each psm
    """
    # List to store q_value at each psm
    q_val = []
    lowest_q_value = data[fdr_col][data.index[-1]]
    for index, row in data[::-1].iterrows():
        curr_q_value = row[fdr_col]
        if curr_q_value < lowest_q_value:
            lowest_q_value = curr_q_value
        q_val.append(lowest_q_value)
    q_val.reverse()
    data["Q_Value"] = q_val
