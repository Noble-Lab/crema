"""
These are unit tests for the Result Class:
"""
import os
import pytest
import filecmp
import pandas as pd
from crema.result import Result


@pytest.fixture
def dataframe():
    """
    Creates a pytest fixture of a Pandas Dataframe to be used
    to test and create Result objects

    Returns
    -------
    Pandas Dataframe
        A pandas dataframe holding psm data with pre-calculated
        fdr and q_value columns
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
            "FDR": [
                1.0,
                1 / 2,
                1 / 3,
                2 / 3,
                1 / 2,
                3 / 4,
                3 / 5,
                1 / 2,
                2 / 3,
            ],
            "Q_Value": [
                1 / 3,
                1 / 3,
                1 / 3,
                1 / 2,
                1 / 2,
                1 / 2,
                1 / 2,
                1 / 2,
                2 / 3,
            ],
        }
    )


@pytest.fixture
def result(dataframe):
    """
    Creates a pytest fixture of a Result object to be used
    to check if internal methods and properties work properly

    Parameters
    ----------
    dataframe : pytest fixture of a pandas dataframe
        The dataframe used as the data parameter to create the Result Object

    Returns
    -------
    PsmDataset
        A Result object
    """
    return Result(
        dataframe,
        ["scan"],
        ["combined p-value0"],
        "target/decoy",
    )


@pytest.fixture
def expected_file_path(dataframe, tmp_path):
    """
    Creates a pytest fixture of a Result object to be used
    to check if internal methods and properties work properly

    Parameters
    ----------
    dataframe : pytest fixture of a pandas dataframe
        The dataframe used as the data parameter to create the Result Object
    tmp_path : pytest fixture of a temporary directory path object

    Returns
    -------
    PsmDataset
        A Result object
    : str
        The filepath used to check against the result object's write function
    """
    file_path = os.path.join(tmp_path, "expected.psm_results.txt")
    dataframe.to_csv(file_path)
    return file_path


def test_create_object(result):
    """
    Ensures that a Result object can be initialized
    properly.

    Parameters
    ----------
    result : pytest fixture of a Result object
        The Result object used for the test
    """
    assert isinstance(result, Result)


def test_properties(result, dataframe):
    """
    Ensures that the data property of the Result object
    is as expected.

    Parameters
    ----------
    result : pytest fixture of a Result object
        The Result object used for the test
    """
    pd.testing.assert_frame_equal(result.data, dataframe)


def test_get_col(result, dataframe):
    """
    Ensures that the get_col() method of the Result object
    properly retrieves the specified column.

    Parameters
    ----------
    result : pytest fixture of a Result object
        The Result object used for the test
    """
    pd.testing.assert_series_equal(result.get_col("scan"), dataframe["scan"])
    pd.testing.assert_series_equal(result.get_col("combined p-value"), dataframe["combined p-value"])
    pd.testing.assert_series_equal(result.get_col("target/decoy"), dataframe["target/decoy"])
    pd.testing.assert_series_equal(result.get_col("FDR"), dataframe["FDR"])
    pd.testing.assert_series_equal(result.get_col("Q_Value"), dataframe["Q_Value"])


def test_write_file(result, tmp_path, expected_file_path):
    """
    Ensures that the write_file() method of the Result object
    properly writes the data to the specified file location.

    Parameters
    ----------
    result : pytest fixture of a Result object
        The Result object used for the test
    tmp_path : pytest fixture of a temporary directory path object
    expected_file_path : str
        The pytest fixture of a filepath that contains the expected
        contents of the result object's write method
    """
    result.write_file(output_dir=tmp_path)
    result.write_file(output_dir=tmp_path, file_root="myFileRoot")
    assert os.path.isfile(os.path.join(tmp_path, "crema.psm_results.txt"))
    assert os.path.isfile(os.path.join(tmp_path, "myFileRootcrema.psm_results.txt"))
    assert filecmp.cmp(os.path.join(tmp_path, "crema.psm_results.txt"), expected_file_path)
