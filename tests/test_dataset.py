"""
These are unit tests for the PSM Dataset Class:
"""
import pytest
import pandas as pd
from crema.dataset import PsmDataset


@pytest.fixture
def dataframe():
    """
    Creates a pytest fixture of a Pandas Dataframe to be used
    to test and create PsmDataset objects

    Returns
    -------
    Pandas Dataframe
        A pandas dataframe holding psm data
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
def psm_dataset(dataframe):
    """
    Creates a pytest fixture of a Psm Dataset to be used
    to check if internal methods and properties work properly

    Parameters
    ----------
    dataframe : pytest fixture of a pandas dataframe
        The dataframe used as the data parameter to create the PsmDataset Object

    Returns
    -------
    PsmDataset
        A PsmDataset object
    """
    return PsmDataset(
        dataframe,
        ["scan"],
        ["combined p-value0", "combined p-value 1", "combined p-value 2"],
        "target/decoy",
    )


def test_create_object(psm_dataset):
    """
    Ensures that a PsmDataset object can be initialized
    properly.

    Parameters
    ----------
    psm_dataset : pytest fixture of a PsmDataset object
        The PsmDataset object used for the test
    """
    assert isinstance(psm_dataset, PsmDataset)


def test_properties(psm_dataset, dataframe):
    """
    Ensures that the data property of the PsmDataset object
    is as expected.

    Parameters
    ----------
    psm_dataset : pytest fixture of a PsmDataset object
        The PsmDataset object used for the test
    """
    pd.testing.assert_frame_equal(psm_dataset.data, dataframe)


def test_get_col(psm_dataset, dataframe):
    """
    Ensures that the get_col() method of the PsmDataset object
    properly retrieves the specified column.

    Parameters
    ----------
    psm_dataset : pytest fixture of a PsmDataset object
        The PsmDataset object used for the test
    """
    pd.testing.assert_series_equal(psm_dataset.get_col("scan"), dataframe["scan"])
    pd.testing.assert_series_equal(psm_dataset.get_col("combined p-value 0"), dataframe["combined p-value 0"])
    pd.testing.assert_series_equal(psm_dataset.get_col("combined p-value 1"), dataframe["combined p-value 1"])
    pd.testing.assert_series_equal(psm_dataset.get_col("combined p-value 2"), dataframe["combined p-value 2"])
    pd.testing.assert_series_equal(psm_dataset.get_col("target/decoy"), dataframe["target/decoy"])
