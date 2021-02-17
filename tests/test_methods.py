"""
These are unit tests for functions within methods.py:
"""
import pytest
import pandas as pd

from crema.result import Result
from crema.dataset import PsmDataset
from crema.methods import (
    _select_columns,
    _delete_duplicates,
    _compare_spectrum_col,
    _calculate_fdr,
    _calculate_q_value,
    calculate_tdc,
)


@pytest.fixture
def dataframe_1():
    """
    Creates a pytest fixture of a pandas DataFrame
    to be converted into a csv file to test read_file().
    Used with test_read_file_add_score().

    Returns
    -------
    pandas.DataFrame
        A pandas.DataFrame containing psm data
    """
    return pd.DataFrame(
        {
            "scan": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],
            "combined p-value 0": [
                0.7,
                0.4,
                0.1,
                0.55,
                0.25,
                0.6,
                0.2,
                0.7,
                0.56,
                0.3,
            ],
            "combined p-value 1": [
                0.1,
                0.4,
                0.7,
                0.8,
                0.5,
                0.4,
                0.7,
                0.6,
                0.3,
                0.6,
            ],
            "combined p-value 2": [
                0.8,
                0.5,
                0.9,
                0.4,
                0.7,
                0.1,
                0.6,
                0.3,
                0.8,
                0.4,
            ],
            "target/decoy": [
                True,
                False,
                True,
                True,
                True,
                True,
                False,
                False,
                True,
                False,
            ],
        }
    )


@pytest.fixture
def dataframe_2():
    """
    Creates a pytest fixture of a pandas DataFrame
    to be used in test_calculate_fdr and test_calculate_q_value

    Returns
    -------
    pandas.DataFrame
        A pandas.DataFrame containing psm data
    """
    return pd.DataFrame(
        {
            "scan": [2, 4, 3, 6, 9, 1, 7, 8, 5],
            "combined p-value": [
                0.1,
                0.15,
                0.2,
                0.3,
                0.4,
                0.45,
                0.6,
                0.75,
                0.95,
            ],
            "target/decoy": [
                True,
                True,
                True,
                False,
                True,
                False,
                True,
                True,
                False,
            ],
        }
    )


@pytest.fixture
def dataframe_fdr(dataframe_2):
    """
    Creates a pytest fixture of a pandas DataFrame
    to be used in test_calculate_fdr and test_calculate_q_value.
    Identical to dataframe_2 with an appended fdr column

    Parameters
    ----------
    dataframe_2 : pytest fixture of a pandas dataframe
        Used as a base dataframe to append fdr column

    Returns
    -------
    pandas.DataFrame
        A pandas.DataFrame containing psm data
    """
    fdr = [
        1.0,
        1 / 2,
        1 / 3,
        2 / 3,
        1 / 2,
        3 / 4,
        3 / 5,
        1 / 2,
        2 / 3,
    ]
    dataframe_2["FDR"] = fdr
    return dataframe_2


@pytest.fixture
def dataframe_q_value(dataframe_fdr):
    """
    Creates a pytest fixture of a pandas DataFrame
    to be used in test_calculate_fdr and test_calculate_q_value.
    Identical to dataframe_fdr with an appended q_value column

    Parameters
    ----------
    dataframe_fdr : pytest fixture of a pandas dataframe
        Used as a base dataframe to append q_value column

    Returns
    -------
    pandas.DataFrame
        A pandas.DataFrame containing psm data
    """
    q_value = [
        1 / 3,
        1 / 3,
        1 / 3,
        1 / 2,
        1 / 2,
        1 / 2,
        1 / 2,
        1 / 2,
        2 / 3,
    ]
    dataframe_fdr["Q_Value"] = q_value
    return dataframe_fdr


@pytest.fixture
def psm_dataset(dataframe_1):
    """
    Creates a pytest fixture of a pandas DataFrame
    to be converted into a csv file to test read_file().
    Used with test_read_file_add_score().

    Parameters
    ----------
    dataframe_1 : pytest fixture of a pandas dataframe
        Used as the data parameter to create the PsmDataset object

    Returns
    -------
    pandas.DataFrame
        A pandas.DataFrame containing psm data
    """
    return PsmDataset(
        dataframe_1,
        ["scan"],
        ["combined p-value 0", "combined p-value 1", "combined p-value 2"],
        "target/decoy",
    )


@pytest.fixture
def psm_dataset_2(dataframe_2):
    """
    Creates a pytest fixture of a pandas DataFrame
    to be converted into a csv file to test read_file().
    Used with test_read_file_add_score().

    Parameters
    ----------
    dataframe_2 : pytest fixture of a pandas dataframe
        Used as the data parameter to create the PsmDataset object

    Returns
    -------
    pandas.DataFrame
        A pandas.DataFrame containing psm data
    """
    return PsmDataset(
        dataframe_2, ["scan"], ["combined p-value"], "target/decoy",
    )


def test_calculate_tdc_full(psm_dataset_2, dataframe_q_value):
    """
    Ensures that the _calculate_tdc() produces a Result object
    that holds all the expected fields; this is inclusive of
    object type, accuracy of dataframe, and correctness of col names

    Parameters
    ----------
    psm_dataset_2 : pytest fixture of a PsmDataset object
        The actual dataframe used for this test case
    dataframe_q_value : pytest fixture of a Pandas Dataframe
        The expected dataframe used to assert against the data
        held in the produced result object
    """
    actual = calculate_tdc(psm_dataset_2)
    # Asserts that the PsmDataset object is an instance of the PsmDataset class
    assert isinstance(actual, Result)
    # Asserts that the data saved in the Result object is as expected
    pd.testing.assert_frame_equal(actual.data, dataframe_q_value)
    # Asserts that the spectrum column name saved in the Result object is as expected
    assert actual.spectrum_col == ["scan"]
    # Asserts that the score column name saved in the Result object is as expected
    assert actual.score_col == "combined p-value"
    # Asserts that the target column name saved in the Result object is as expected
    assert actual.target_col == "target/decoy"
    # Asserts that the fdr column name saved in the Result object is as expected
    assert actual.fdr_col == "FDR"
    # Asserts that the q_value column name saved in the Result object is as expected
    assert actual.q_val_col == "Q_Value"


def test_calculate_tdc_data(psm_dataset):
    """
    Ensures that the _calculate_tdc() function properly
    executes all helper functions in order to generate a
    Result object (with FDR and q_value columns) with
    the correct dataframe values. Tests three different score
    columns of the provided psm_dataset.

    Parameters
    ----------
    psm_dataset : pytest fixture of a PsmDataset object
        The actual dataframe used for this test case
    """
    actual0 = calculate_tdc(psm_dataset, 0).data
    actual1 = calculate_tdc(psm_dataset, 1).data
    actual2 = calculate_tdc(psm_dataset, 2).data
    expected0 = pd.DataFrame(
        {
            "scan": [3, 2, 5, 4, 1],
            "combined p-value 0": [0.1, 0.2, 0.25, 0.55, 0.6],
            "target/decoy": [True, False, True, True, True],
            "FDR": [1, 1, 1, 2 / 3, 1 / 2],
            "Q_Value": [0.5, 0.5, 0.5, 0.5, 0.5],
        }
    )
    expected1 = pd.DataFrame(
        {
            "scan": [1, 4, 2, 5, 3],
            "combined p-value 1": [0.1, 0.3, 0.4, 0.5, 0.6],
            "target/decoy": [True, True, False, True, False],
            "FDR": [1, 0.5, 1, 2 / 3, 1],
            "Q_Value": [0.5, 0.5, 2 / 3, 2 / 3, 1],
        }
    )
    expected2 = pd.DataFrame(
        {
            "scan": [1, 3, 4, 5, 2],
            "combined p-value 2": [0.1, 0.3, 0.4, 0.4, 0.5],
            "target/decoy": [True, False, True, False, False],
            "FDR": [1.0, 1.0, 1.0, 1.0, 1.0],
            "Q_Value": [1.0, 1.0, 1.0, 1.0, 1.0],
        }
    )
    pd.testing.assert_frame_equal(actual0, expected0, check_index_type=False)
    pd.testing.assert_frame_equal(actual1, expected1, check_index_type=False)
    pd.testing.assert_frame_equal(actual2, expected2, check_index_type=False)


def test_select_columns(psm_dataset):
    """
    Ensures that the _select_columns() function properly
    extracts the specified columns from the psm dataset object

    Parameters
    ----------
    psm_dataset : pytest fixture of a PsmDataset object
        The actual dataframe used for this test case
    """
    expected = pd.DataFrame(
        {
            "scan": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],
            "combined p-value 0": [
                0.7,
                0.4,
                0.1,
                0.55,
                0.25,
                0.6,
                0.2,
                0.7,
                0.56,
                0.3,
            ],
            "target/decoy": [
                True,
                False,
                True,
                True,
                True,
                True,
                False,
                False,
                True,
                False,
            ],
        }
    )
    actual = _select_columns(
        psm_dataset, ["scan"], "combined p-value 0", "target/decoy"
    )
    # Index types are different due to construction of dataframes - ignore since we care only about correct values
    pd.testing.assert_frame_equal(actual, expected, check_index_type=False)


def test_compare_spectrum_columns(dataframe_1):
    """
    Ensures that the _compare_spectrum_columns() function properly
    compares the spectrum column value of two given dataframe rows.
    Should return True if they are different,
    and False if they are the same.

    Parameters
    ----------
    dataframe_1 : pytest fixture of a Pandas Dataframe
    """
    # scan column 0th index != 1st index -> return True
    assert (
        _compare_spectrum_col(
            dataframe_1[["scan"]], dataframe_1.loc[0], dataframe_1.loc[1]
        )
        is True
    )
    # scan column 1st index == 6th index -> return False
    assert (
        _compare_spectrum_col(
            dataframe_1[["scan"]], dataframe_1.loc[1], dataframe_1.loc[6]
        )
        is False
    )
    # scan column 3rd index == 8th index -> return False
    assert (
        _compare_spectrum_col(
            dataframe_1[["scan"]], dataframe_1.loc[3], dataframe_1.loc[8]
        )
        is False
    )


def test_delete_duplicates():
    """
    Ensures that the _delete_duplicates() function properly
    compares each row of a dataframe against the next row,
    and deletes the next row if it is a "duplicate" - that is,
    the spectrum column value is identical. Note that the dataframe
    passed to this function must be sorted by
    spectrum column, score (ascending), and target/decoy value.
    """
    data = pd.DataFrame(
        {
            "scan": [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5],
            "combined p-value": [
                0.1,
                0.1,
                0.3,
                0.4,
                0.5,
                0.6,
                0.7,
                0.8,
                0.9,
                0.2,
                0.5,
                0.7,
                0.3,
                0.6,
                0.9,
            ],
            "target/decoy": [
                True,
                False,
                False,
                True,
                True,
                False,
                True,
                True,
                True,
                False,
                False,
                False,
                True,
                True,
                False,
            ],
        }
    )
    expected = pd.DataFrame(
        {
            "scan": [1, 2, 3, 4, 5],
            "combined p-value": [
                0.1,
                0.4,
                0.7,
                0.2,
                0.3,
            ],
            "target/decoy": [
                True,
                True,
                True,
                False,
                True,
            ],
        }
    )
    actual = _delete_duplicates(
        data, ["scan"], "combined p-value", "target/decoy"
    ).reset_index(drop=True)
    pd.testing.assert_frame_equal(actual, expected, check_index_type=False)


def test_calculate_fdr(dataframe_2, dataframe_fdr):
    """
    Ensures that the _calculate_fdr() function properly
    calculates and appends an FDR column to the given dataframe

    Parameters
    ----------
    dataframe_2 : pytest fixture of a Pandas Dataframe
        The actual dataframe used for this test case
    dataframe_fdr : pytest fixture of a Pandas Dataframe
        The expected dataframe result to assert against
    """
    data = dataframe_2
    expected = dataframe_fdr
    actual = _calculate_fdr(data, "target/decoy")
    pd.testing.assert_frame_equal(actual, expected, check_index_type=False)


def test_calculate_q_value(dataframe_fdr, dataframe_q_value):
    """
    Ensures that the _calculate_q_value() function properly
    calculates and appends a q_value column to the given dataframe

    Parameters
    ----------
    dataframe_fdr : pytest fixture of a Pandas Dataframe
        The actual dataframe used for this test case
    dataframe_q_value : pytest fixture of a Pandas Dataframe
        The expected dataframe result to assert against
    """
    data = dataframe_fdr
    expected = dataframe_q_value
    actual = _calculate_q_value(data, "FDR")
    pd.testing.assert_frame_equal(actual, expected, check_index_type=False)
