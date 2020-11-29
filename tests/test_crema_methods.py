"""
These tests verify that the everything works as intended:
1. dataset objects are properly created
2. fdr methods are properly executed
3. result objects are properly created
"""

import pytest
from crema.parsers import *
from crema.methods import *


# The following dataframes are what I use as "actual" in the following tests:
# i.e. this is what the pandas dataframe should look like after the dataset/result object is created

# Used in test_single_data
testframe_single = pd.DataFrame(
    {
        "scan": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],
        "p-value": [0.7, 0.4, 0.1, 0.55, 0.3, 0.6, 0.2, 0.7, 0.56, 0.3],
        "target": [
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

# Used in test_single_basic_data
testframe_single_basic = pd.DataFrame(
    {
        "scan": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],
        "combined p-value": [0.7, 0.4, 0.1, 0.55, 0.3, 0.6, 0.2, 0.7, 0.56, 0.3],
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

# Used in test_single_text_scan_data
testframe_single_text_scan = pd.DataFrame(
    {
        "scan": [
            "C:/test/string/fake/file/path/1",
            "C:/test/string/fake/file/path/2",
            "C:/test/string/fake/file/path/3",
            "C:/test/string/fake/file/path/4",
            "C:/test/string/fake/file/path/5",
            "C:/test/string/fake/file/path/1",
            "C:/test/string/fake/file/path/2",
            "C:/test/string/fake/file/path/3",
            "C:/test/string/fake/file/path/4",
            "C:/test/string/fake/file/path/5",
        ],
        "p-value": [0.7, 0.4, 0.1, 0.55, 0.3, 0.6, 0.2, 0.7, 0.56, 0.3],
        "target": [
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

# Used in test_single_tdc
# Result if TDC arbitrarily chooses True
testframe_single_tdc_true = pd.DataFrame(
    {
        "scan": [3, 2, 5, 4, 1],
        "p-value": [0.1, 0.2, 0.3, 0.55, 0.6],
        "target": [True, False, True, True, True],
        "FDR": [1, 1, 1, 2 / 3, 1 / 2],
        "Q_Value": [0.5, 0.5, 0.5, 0.5, 0.5],
    }
)

# Used in test_single_tdc
# Result if TDC arbitrarily chooses False
testframe_single_tdc_false = pd.DataFrame(
    {
        "scan": [3, 2, 5, 4, 1],
        "p-value": [0.1, 0.2, 0.3, 0.55, 0.6],
        "target": [True, False, False, True, True],
        "FDR": [1.0, 1.0, 1.0, 1.0, 1.0],
        "Q_Value": [1.0, 1.0, 1.0, 1.0, 1.0],
    }
)

# Used in test_single_text_scan_tdc
# Result if TDC arbitrarily chooses True
testframe_single_text_scan_tdc_true = pd.DataFrame(
    {
        "scan": [
            "C:/test/string/fake/file/path/3",
            "C:/test/string/fake/file/path/2",
            "C:/test/string/fake/file/path/5",
            "C:/test/string/fake/file/path/4",
            "C:/test/string/fake/file/path/1",
        ],
        "p-value": [0.1, 0.2, 0.3, 0.55, 0.6],
        "target": [True, False, True, True, True],
        "FDR": [1, 1, 1, 2 / 3, 1 / 2],
        "Q_Value": [0.5, 0.5, 0.5, 0.5, 0.5],
    }
)

# Used in test_single_text_scan_tdc
# Result if TDC arbitrarily chooses False
testframe_single_text_scan_tdc_false = pd.DataFrame(
    {
        "scan": [
            "C:/test/string/fake/file/path/3",
            "C:/test/string/fake/file/path/2",
            "C:/test/string/fake/file/path/5",
            "C:/test/string/fake/file/path/4",
            "C:/test/string/fake/file/path/1",
        ],
        "p-value": [0.1, 0.2, 0.3, 0.55, 0.6],
        "target": [True, False, False, True, True],
        "FDR": [1.0, 1.0, 1.0, 1.0, 1.0],
        "Q_Value": [1.0, 1.0, 1.0, 1.0, 1.0],
    }
)

# Used in test_multi_data
testframe_multi = pd.DataFrame(
    {
        "scan": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5],
        "p-value": [
            0.7,
            0.4,
            0.1,
            0.55,
            0.3,
            0.6,
            0.2,
            0.7,
            0.56,
            0.3,
            0.4,
            0.7,
            0.59,
            0.55,
            0.9,
            0.44,
            0.75,
            0.6,
            0.5,
            0.89,
        ],
        "target": [
            True,
            True,
            True,
            True,
            True,
            True,
            True,
            True,
            True,
            True,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
        ],
    }
)

# Used in test_multi_tdc
testframe_multi_tdc = pd.DataFrame(
    {
        "scan": [3, 2, 5, 1, 4],
        "p-value": [0.1, 0.2, 0.3, 0.4, 0.5],
        "target": [True, True, True, False, False],
        "FDR": [1, 1 / 2, 1 / 3, 2 / 3, 1],
        "Q_Value": [1 / 3, 1 / 3, 1 / 3, 2 / 3, 1],
    }
)


@pytest.fixture
def test_single_basic_dataset_class():
    """
    Creates a pytest fixture of a PsmDataset object
    by reading in "single_basic.csv" to use in
    subsequent test cases.

    Returns
    -------
    PsmDataset
        A :py:class:`~crema.dataset.PsmDataset` object
        containing the PSM data from the given tab-delimited file.
    """
    psm = read_file(["data/single_basic.csv"])
    return psm


@pytest.fixture
def test_single_dataset_class():
    """
    Creates a pytest fixture of a PsmDataset object
    by reading in "single.csv" to use in
    subsequent test cases.

    Returns
    -------
    PsmDataset
        A :py:class:`~crema.dataset.PsmDataset` object
        containing the PSM data from the given tab-delimited file.
    """
    psm = read_file(["data/single.csv"], "scan", "p-value", "target")
    return psm


@pytest.fixture
def test_single_int_targets_dataset_class():
    """
    Creates a pytest fixture of a PsmDataset object
    by reading in "single.csv" to use in
    subsequent test cases. Note that target column
    contain integer values (0 and 1), but also
    works with (-1 and 1).

    Returns
    -------
    PsmDataset
        A :py:class:`~crema.dataset.PsmDataset` object
        containing the PSM data from the given tab-delimited file.
    """
    psm = read_file(
        ["data/single_int_targets.csv"], "scan", "p-value", "target"
    )
    return psm


@pytest.fixture
def test_single_text_scan_dataset_class():
    """
    Creates a pytest fixture of a PsmDataset object
    by reading in "single.csv" to use in
    subsequent test cases. Note that target column
    contain integer values (0 and 1), but also
    works with (-1 and 1).

    Returns
    -------
    PsmDataset
        A :py:class:`~crema.dataset.PsmDataset` object
        containing the PSM data from the given tab-delimited file.
    """
    psm = read_file(["data/single_text_scan.csv"], "scan", "p-value", "target")
    return psm


@pytest.fixture
def test_multi_dataset_class():
    """
    Creates a pytest fixture of a PsmDataset object
    by reading in "multi_target.csv" and
    "multi_target.csv" to use in subsequent test cases.

    Returns
    -------
    PsmDataset
        A :py:class:`~crema.dataset.PsmDataset` object
        containing the PSM data from the given tab-delimited files.
    """
    files = ["data/multi_target.csv", "data/multi_decoy.csv"]
    psm = read_file(files, "scan", "p-value", "target")
    return psm


def test_multi_tdc(test_multi_dataset_class):
    """
    Checks whether or not a Result object is created properly
    after executing the calculate_tdc method on a PsmDataset object
    created with multiple files.

    Parameters
    ----------
    test_multi_dataset_class : pytest fixture of a PsmDataset object
        A psm object created by reading in the "multi_target.csv" file
        and the "multi_decoy.csv" file.

    Returns
    -------
    Pandas Assert Frame
        Asserts whether or not the data in Result object
        is equal to the dataframe named "testframe_multi_tdc"
    """
    actual = testframe_multi_tdc.copy()
    output = calculate_tdc(test_multi_dataset_class)
    # Following 2 lines tests if result can be properly written to file:
    # output.write_csv("C:/Users/donav/Github/TestDatasets/test_multi_tdc.csv")
    # output.write_excel("C:/Users/donav/Github/TestDatasets/test_multi_tdc.xlsx")
    compare = output.data
    pd.testing.assert_frame_equal(actual, compare)


def test_multi_data(test_multi_dataset_class):
    """
    Checks whether or not a PsmDataset object is created properly
    after reading in multiple tab delimited files.

    Parameters
    ----------
    test_multi_dataset_class : pytest fixture of a PsmDataset object
        A psm object created by reading in the "multi_target.csv" file
        and the "multi_decoy.csv" file.

    Returns
    -------
    Pandas Assert Frame
        Asserts whether or not the data in PsmDataset object
        is equal to the dataframe named "testframe_multi"
    """
    actual = testframe_multi.copy()
    compare = test_multi_dataset_class.data
    pd.testing.assert_frame_equal(actual, compare)


def test_single_tdc(test_single_dataset_class):
    """
    Checks whether or not a Result object is created properly
    after executing the calculate_tdc method on a PsmDataset object
    created with a single file.

    Parameters
    ----------
    test_single_dataset_class : pytest fixture of a PsmDataset object
        A psm object created by reading in the "single.csv" file

    Returns
    -------
    Pandas Assert Frame
        Asserts whether or not the data in Result object
        is equal to the dataframe named "testframe_single_tdc"
    """
    # Note that this tests a file where there are duplicate target/decoy PSMs with equal P-Value
    # The algorithm arbitrarily chooses, in which case there are two possible end states
    actual_true = testframe_single_tdc_true.copy()
    actual_false = testframe_single_tdc_false.copy()
    output = calculate_tdc(test_single_dataset_class)
    compare = output.data
    if compare.iloc[2, 2]:
        # If arbitrarily chooses True
        pd.testing.assert_frame_equal(actual_true, compare)
    else:
        # If arbitrarily chooses False
        pd.testing.assert_frame_equal(actual_false, compare)


def test_single_data(test_single_dataset_class):
    """
    Checks whether or not a PsmDataset object is created properly
    after reading in a single tab delimited file.

    Parameters
    ----------
    test_single_dataset_class : pytest fixture of a PsmDataset object
        A psm object created by reading in the "single.csv" file

    Returns
    -------
    Pandas Assert Frame
        Asserts whether or not the data in PsmDataset object
        is equal to the dataframe named "testframe_single"
    """
    actual = testframe_single.copy()
    compare = test_single_dataset_class.data
    pd.testing.assert_frame_equal(actual, compare)


def test_single_basic_data(test_single_basic_dataset_class):
    """
    Checks whether or not a PsmDataset object is created properly
    after reading in a single tab delimited file.

    Parameters
    ----------
    test_single_basic_dataset_class : pytest fixture of a PsmDataset object
        A psm object created by reading in the "single_basic.csv" file

    Returns
    -------
    Pandas Assert Frame
        Asserts whether or not the data in PsmDataset object
        is equal to the dataframe named "testframe_single"
    """
    actual = testframe_single_basic.copy()
    compare = test_single_basic_dataset_class.data
    pd.testing.assert_frame_equal(actual, compare)


def test_single_int_targets_data(test_single_int_targets_dataset_class):
    """
    Checks whether or not a PsmDataset object is created properly
    after reading in a single tab delimited file.

    Parameters
    ----------
    test_single_int_targets_dataset_class : pytest fixture of a PsmDataset object
        A psm object created by reading in the "single_int_targets.csv" file

    Returns
    -------
    Pandas Assert Frame
        Asserts whether or not the data in PsmDataset object
        is equal to the dataframe named "testframe_single"
    """
    actual = testframe_single.copy()
    compare = test_single_int_targets_dataset_class.data
    pd.testing.assert_frame_equal(actual, compare)


def test_single_text_scan_data(test_single_text_scan_dataset_class):
    """
    Checks whether or not a PsmDataset object is created properly
    after reading in a single tab delimited file.

    Parameters
    ----------
    test_single_text_scan_dataset_class : pytest fixture of a PsmDataset object
        A psm object created by reading in the "single_text_scan.csv" file

    Returns
    -------
    Pandas Assert Frame
        Asserts whether or not the data in PsmDataset object
        is equal to the dataframe named "testframe_single_text_scan"
    """
    actual = testframe_single_text_scan.copy()
    compare = test_single_text_scan_dataset_class.data
    pd.testing.assert_frame_equal(actual, compare)


def test_single_text_scan_tdc(test_single_text_scan_dataset_class):
    """
    Checks whether or not a Result object is created properly
    after executing the calculate_tdc method on a PsmDataset object
    created with a single file.

    Parameters
    ----------
    test_single_text_scan_dataset_class : pytest fixture of a PsmDataset object
        A psm object created by reading in the "single_text_scan.csv" file

    Returns
    -------
    Pandas Assert Frame
        Asserts whether or not the data in Result object
        is equal to the dataframe named "testframe_single_text_scan_tdc"
    """
    # Note that this tests a file where there are duplicate target/decoy PSMs with equal P-Value
    # The algorithm arbitrarily chooses, in which case there are two possible end states
    # Furthermore, this tests data that has a scan column of strings instead of int
    # It also tests a target column of numeric value instead of boolean
    actual_true = testframe_single_text_scan_tdc_true.copy()
    actual_false = testframe_single_text_scan_tdc_false.copy()
    output = calculate_tdc(test_single_text_scan_dataset_class)
    compare = output.data
    if compare.iloc[2, 2]:
        # If arbitrarily chooses True
        pd.testing.assert_frame_equal(actual_true, compare)
    else:
        # If arbitrarily chooses False
        pd.testing.assert_frame_equal(actual_false, compare)
