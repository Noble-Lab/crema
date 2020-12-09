"""
These tests verify that the crema api works as intended:
1. dataset objects are properly created
2. fdr methods are properly executed
3. result objects are properly created
"""

import pytest
from crema.parsers import *
from crema.methods import *


# The following dataframes are what I use as "actual" in the following tests:
# i.e. this is what the pandas dataframe should look like after the dataset/result object is created

# Used in test_single_basic_data, test_single_int_targets_data
testframe_single_basic_data = pd.DataFrame(
    {
        "scan": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],
        "combined p-value": [
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

# Used in test_single_basic_tdc, test_single_int_targets_tdc
testframe_single_basic_tdc = pd.DataFrame(
    {
        "scan": [3, 2, 5, 4, 1],
        "combined p-value": [0.1, 0.2, 0.25, 0.55, 0.6],
        "target/decoy": [True, False, True, True, True],
        "FDR": [1, 1, 1, 2 / 3, 1 / 2],
        "Q_Value": [0.5, 0.5, 0.5, 0.5, 0.5],
    }
)

# Used in test_single_text_scan_data
testframe_single_text_scan_data = pd.DataFrame(
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
        "combined p-value": [
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

# Used in test_single_text_scan_tdc
testframe_single_text_scan_tdc = pd.DataFrame(
    {
        "scan": [
            "C:/test/string/fake/file/path/3",
            "C:/test/string/fake/file/path/2",
            "C:/test/string/fake/file/path/5",
            "C:/test/string/fake/file/path/4",
            "C:/test/string/fake/file/path/1",
        ],
        "combined p-value": [0.1, 0.2, 0.25, 0.55, 0.6],
        "target/decoy": [True, False, True, True, True],
        "FDR": [1, 1, 1, 2 / 3, 1 / 2],
        "Q_Value": [0.5, 0.5, 0.5, 0.5, 0.5],
    }
)

# Used in test_single_add_spectrum
testframe_single_add_spectrum_data = pd.DataFrame(
    {
        "scan": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],
        "combined p-value": [
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
        "extras": [
            "a",
            "b",
            "c",
            "d",
            "e",
            "z",
            "b",
            "y",
            "d",
            "x",
        ],
    }
)

# Used in test_single_add_spectrum_tdc
testframe_single_add_spectrum_tdc = pd.DataFrame(
    {
        "scan": [3, 2, 5, 5, 4, 1, 1, 3],
        "combined p-value": [0.1, 0.2, 0.25, 0.3, 0.55, 0.6, 0.7, 0.7],
        "target/decoy": [True, False, True, False, True, True, True, False],
        "extras": ["c", "b", "e", "x", "d", "z", "a", "y"],
        "FDR": [1, 1, 1, 1, 1, 3 / 4, 3 / 5, 4 / 5],
        "Q_Value": [3 / 5, 3 / 5, 3 / 5, 3 / 5, 3 / 5, 3 / 5, 3 / 5, 4 / 5],
    }
)

# Used in test_single_arbitrary_data
testframe_single_arbitrary_data = pd.DataFrame(
    {
        "scan": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],
        "combined p-value": [
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

# Used in test_single_arbitrary_tdc
# Result if TDC arbitrarily chooses True
testframe_single_arbitrary_tdc_true = pd.DataFrame(
    {
        "scan": [3, 2, 5, 4, 1],
        "combined p-value": [0.1, 0.2, 0.3, 0.55, 0.6],
        "target/decoy": [True, False, True, True, True],
        "FDR": [1, 1, 1, 2 / 3, 1 / 2],
        "Q_Value": [0.5, 0.5, 0.5, 0.5, 0.5],
    }
)

# Used in test_single_arbitrary_tdc
# Result if TDC arbitrarily chooses False
testframe_single_arbitrary_tdc_false = pd.DataFrame(
    {
        "scan": [3, 2, 5, 4, 1],
        "combined p-value": [0.1, 0.2, 0.3, 0.55, 0.6],
        "target/decoy": [True, False, False, True, True],
        "FDR": [1.0, 1.0, 1.0, 1.0, 1.0],
        "Q_Value": [1.0, 1.0, 1.0, 1.0, 1.0],
    }
)

# Used in test_single_noncrux_data
testframe_single_noncrux_data = pd.DataFrame(
    {
        "scan": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],
        "p-value": [
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

# Used in test_single_noncrux_tdc
testframe_single_noncrux_tdc = pd.DataFrame(
    {
        "scan": [3, 2, 5, 4, 1],
        "p-value": [0.1, 0.2, 0.25, 0.55, 0.6],
        "target": [True, False, True, True, True],
        "FDR": [1, 1, 1, 2 / 3, 1 / 2],
        "Q_Value": [0.5, 0.5, 0.5, 0.5, 0.5],
    }
)

# Used in test_multi_data
testframe_multi_data = pd.DataFrame(
    {
        "scan": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5],
        "combined p-value": [
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
        "target/decoy": [
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
        "combined p-value": [0.1, 0.2, 0.3, 0.4, 0.5],
        "target/decoy": [True, True, True, False, False],
        "FDR": [1, 1 / 2, 1 / 3, 2 / 3, 1],
        "Q_Value": [1 / 3, 1 / 3, 1 / 3, 2 / 3, 1],
    }
)


@pytest.fixture
def test_single_basic_dataset_class():
    """
    Creates a pytest fixture of a PsmDataset object
    by reading in "single_basic.csv" to use in
    subsequent test cases. This file has crux
    default column names.

    Returns
    -------
    PsmDataset
        A :py:class:`~crema.dataset.PsmDataset` object
        containing the PSM data from the given comma separated value file.
    """
    psm = read_file(["data/single_basic.csv"])
    return psm


@pytest.fixture
def test_single_basic_tab_dataset_class():
    """
    Creates a pytest fixture of a PsmDataset object
    by reading in "single_basic_tab.csv" to use in
    subsequent test cases. This file has crux
    default column names.

    Returns
    -------
    PsmDataset
        A :py:class:`~crema.dataset.PsmDataset` object
        containing the PSM data from the given tab delimited file.
    """
    psm = read_file(["data/single_basic_tab.txt"])
    return psm


@pytest.fixture
def test_single_int_targets_dataset_class():
    """
    Creates a pytest fixture of a PsmDataset object
    by reading in "single_int_targets.csv" to use in
    subsequent test cases. This file has crux default
    column names. However, note that target column
    contains integer values (0 and 1), but also
    works with (-1 and 1).

    Returns
    -------
    PsmDataset
        A :py:class:`~crema.dataset.PsmDataset` object
        containing the PSM data from the given comma separated value file.
    """
    psm = read_file(["data/single_int_targets.csv"])
    return psm


@pytest.fixture
def test_single_text_scan_dataset_class():
    """
    Creates a pytest fixture of a PsmDataset object
    by reading in "single_text_scan.csv" to use in
    subsequent test cases. This file has crux default
    column names. However, note that scan column
    contains String identifiers.

    Returns
    -------
    PsmDataset
        A :py:class:`~crema.dataset.PsmDataset` object
        containing the PSM data from the given comma separated value file.
    """
    psm = read_file(["data/single_text_scan.csv"])
    return psm


@pytest.fixture
def test_single_add_spectrum_dataset_class():
    """
    Creates a pytest fixture of a PsmDataset object
    by reading in "single_add_spectrum.csv" to use in
    subsequent test cases. This file has crux default
    column names with an additional "extras" column
    to serve as an additional spectrum column input.

    Returns
    -------
    PsmDataset
        A :py:class:`~crema.dataset.PsmDataset` object
        containing the PSM data from the given comma separated value file.
    """
    psm = read_file(
        ["data/single_add_spectrum.csv"], spectrum_col=["scan", "extras"]
    )
    return psm


@pytest.fixture
def test_single_arbitrary_dataset_class():
    """
    Creates a pytest fixture of a PsmDataset object
    by reading in "single_arbitrary.csv" to use in
    subsequent test cases. This file has crux default
    column names. However, there are target/decoy psms
    with equal p-value to force an arbitrary choice.

    Returns
    -------
    PsmDataset
        A :py:class:`~crema.dataset.PsmDataset` object
        containing the PSM data from the given comma separated value file.
    """
    psm = read_file(["data/single_arbitrary.csv"])
    return psm


@pytest.fixture
def test_single_noncrux_dataset_class():
    """
    Creates a pytest fixture of a PsmDataset object
    by reading in "single_noncrux.csv" to use in
    subsequent test cases. This file has non-crux
    column names:
    "scan", "p-value", "target".

    Returns
    -------
    PsmDataset
        A :py:class:`~crema.dataset.PsmDataset` object
        containing the PSM data from the given comma separated value file.
    """
    psm = read_file(["data/single_noncrux.csv"], "scan", "p-value", "target")
    return psm


@pytest.fixture
def test_multi_dataset_class():
    """
    Creates a pytest fixture of a PsmDataset object
    by reading in "multi_target.csv" and
    "multi_target.csv" to use in subsequent test cases.
    Both of these files have crux default
    column names.

    Returns
    -------
    PsmDataset
        A :py:class:`~crema.dataset.PsmDataset` object
        containing the PSM data from the given comma separated value files.
    """
    files = ["data/multi_target.csv", "data/multi_decoy.csv"]
    psm = read_file(files)
    return psm


def test_single_basic_data(test_single_basic_dataset_class):
    """
    Checks whether or not a PsmDataset object is created properly
    after reading in a single tab delimited file in its most basic
    form with crux default column names. The point of this test is
    to ensure that there are no errors with the most basic input.

    Parameters
    ----------
    test_single_basic_dataset_class : pytest fixture of a PsmDataset object
        A psm object created by reading in the "single_basic.csv" file

    Returns
    -------
    Pandas Assert Frame
        Asserts whether or not the data in PsmDataset object
        is equal to the dataframe named "testframe_single_basic_data"
    """
    actual = testframe_single_basic_data.copy()
    compare = test_single_basic_dataset_class.data
    pd.testing.assert_frame_equal(actual, compare)


def test_single_basic_tdc(test_single_basic_dataset_class):
    """
    Checks whether or not a Result object is created properly
    after executing the calculate_tdc method on a PsmDataset object
    from test_single_basic_data.

    Parameters
    ----------
    test_single_basic_dataset_class : pytest fixture of a PsmDataset object
        A psm object created by reading in the "single_basic.csv" file

    Returns
    -------
    Pandas Assert Frame
        Asserts whether or not the data in Result object
        is equal to the dataframe named "testframe_single_basic_tdc"
    """
    actual = testframe_single_basic_tdc.copy()
    output = calculate_tdc(test_single_basic_dataset_class)
    compare = output.data
    pd.testing.assert_frame_equal(actual, compare)


def test_single_basic_tab_data(test_single_basic_tab_dataset_class):
    """
    Checks whether or not a PsmDataset object is created properly
    after reading in a single tab delimited file with crux default column names,
    but tab delimited instead of comma separated. The point of this test is
    to ensure that crema can read data separated by tabs as well as commas.

    Parameters
    ----------
    test_single_basic_tab_dataset_class : pytest fixture of a PsmDataset object
        A psm object created by reading in the "single_basic_tab.txt" file

    Returns
    -------
    Pandas Assert Frame
        Asserts whether or not the data in PsmDataset object
        is equal to the dataframe named "testframe_single_basic_data"
    """
    actual = testframe_single_basic_data.copy()
    compare = test_single_basic_tab_dataset_class.data
    pd.testing.assert_frame_equal(actual, compare)


def test_single_int_targets_data(test_single_int_targets_dataset_class):
    """
    Checks whether or not a PsmDataset object is created properly
    after reading in a single tab delimited file with crux default
    column names and target column values containing integers.
    The point of this test is to make sure calculations can be performed
    with integer identifiers (1/0 or 1/-1) for target/decoy psms.

    Parameters
    ----------
    test_single_int_targets_dataset_class : pytest fixture of a PsmDataset object
        A psm object created by reading in the "single_int_targets.csv" file

    Returns
    -------
    Pandas Assert Frame
        Asserts whether or not the data in PsmDataset object
        is equal to the dataframe named "testframe_single_basic_data"
    """
    actual = testframe_single_basic_data.copy()
    compare = test_single_int_targets_dataset_class.data
    pd.testing.assert_frame_equal(actual, compare)


def test_single_int_targets_tdc(test_single_int_targets_dataset_class):
    """
    Checks whether or not a Result object is created properly
    after executing the calculate_tdc method on a PsmDataset object
    from test_single_int_targets_data.

    Parameters
    ----------
    test_single_int_targets_dataset_class : pytest fixture of a PsmDataset object
        A psm object created by reading in the "single_int_targets.csv" file

    Returns
    -------
    Pandas Assert Frame
        Asserts whether or not the data in PsmDataset object
        is equal to the dataframe named "testframe_single_basic_tdc"
    """
    actual = testframe_single_basic_tdc.copy()
    output = calculate_tdc(test_single_int_targets_dataset_class)
    compare = output.data
    pd.testing.assert_frame_equal(actual, compare)


def test_single_text_scan_data(test_single_text_scan_dataset_class):
    """
    Checks whether or not a PsmDataset object is created properly
    after reading in a single tab delimited file with crux default
    column names and scan column values containing Strings.
    The point of this test is to make sure calculations can be performed
    with non integer scan identifiers, i.e. filename/#.

    Parameters
    ----------
    test_single_text_scan_dataset_class : pytest fixture of a PsmDataset object
        A psm object created by reading in the "single_text_scan.csv" file

    Returns
    -------
    Pandas Assert Frame
        Asserts whether or not the data in PsmDataset object
        is equal to the dataframe named "testframe_single_text_scan_data"
    """
    actual = testframe_single_text_scan_data.copy()
    compare = test_single_text_scan_dataset_class.data
    pd.testing.assert_frame_equal(actual, compare)


def test_single_text_scan_tdc(test_single_text_scan_dataset_class):
    """
    Checks whether or not a Result object is created properly
    after executing the calculate_tdc method on a PsmDataset object
    from test_single_text_scan_data.

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
    actual = testframe_single_text_scan_tdc.copy()
    output = calculate_tdc(test_single_text_scan_dataset_class)
    compare = output.data
    pd.testing.assert_frame_equal(actual, compare)


def test_single_add_spectrum_data(test_single_add_spectrum_dataset_class):
    """
    Checks whether or not a PsmDataset object is created properly
    after reading in a single tab delimited file with crux default
    column names and an additional spectrum column.
    The point of this test is to make sure the user can specify more
    than one spectrum columns if the need arises.

    Parameters
    ----------
    test_single_add_spectrum_dataset_class : pytest fixture of a PsmDataset object
        A psm object created by reading in the "single_add_spectrum.csv" file

    Returns
    -------
    Pandas Assert Frame
        Asserts whether or not the data in PsmDataset object
        is equal to the dataframe named "testframe_single_add_spectrum_data"
    """
    actual = testframe_single_add_spectrum_data.copy()
    compare = test_single_add_spectrum_dataset_class.data
    pd.testing.assert_frame_equal(actual, compare)


def test_single_add_spectrum_tdc(test_single_add_spectrum_dataset_class):
    """
    Checks whether or not a Result object is created properly
    after executing the calculate_tdc method on a PsmDataset object
    from test_single_add_spectrum_data.

    Parameters
    ----------
    test_single_add_spectrum_dataset_class : pytest fixture of a PsmDataset object
        A psm object created by reading in the "single_add_spectrum.csv" file

    Returns
    -------
    Pandas Assert Frame
        Asserts whether or not the data in Result object
        is equal to the dataframe named "testframe_single_add_spectrum_tdc"
    """
    actual = testframe_single_add_spectrum_tdc.copy()
    output = calculate_tdc(test_single_add_spectrum_dataset_class)
    compare = output.data
    pd.testing.assert_frame_equal(actual, compare)


def test_single_arbitrary_data(test_single_arbitrary_dataset_class):
    """
    Checks whether or not a PsmDataset object is created properly
    after reading in a single tab delimited file with crux default
    column names and duplicate target/decoy psms of equal p-value.
    The point of this test is to make sure that crema arbitrarily
    chooses which psm to keep.

    Parameters
    ----------
    test_single_arbitrary_dataset_class : pytest fixture of a PsmDataset object
        A psm object created by reading in the "single_arbitrary.csv" file

    Returns
    -------
    Pandas Assert Frame
        Asserts whether or not the data in PsmDataset object
        is equal to the dataframe named "testframe_single_arbitrary_datadata"
    """
    actual = testframe_single_arbitrary_data.copy()
    compare = test_single_arbitrary_dataset_class.data
    pd.testing.assert_frame_equal(actual, compare)


def test_single_arbitrary_tdc(test_single_arbitrary_dataset_class):
    """
    Checks whether or not a Result object is created properly
    after executing the calculate_tdc method on a PsmDataset object
    from test_single_arbitrary_data.

    Parameters
    ----------
    test_single_arbitrary_dataset_class : pytest fixture of a PsmDataset object
        A psm object created by reading in the "single_arbitrary.csv" file

    Returns
    -------
    Pandas Assert Frame
        Asserts whether or not the data in Result object
        is equal to the dataframe named "testframe_single_arbitrary_tdc"
    """
    actual_true = testframe_single_arbitrary_tdc_true.copy()
    actual_false = testframe_single_arbitrary_tdc_false.copy()
    output = calculate_tdc(test_single_arbitrary_dataset_class)
    compare = output.data
    if compare.iloc[2, 2]:
        # If arbitrarily chooses True
        pd.testing.assert_frame_equal(actual_true, compare)
    else:
        # If arbitrarily chooses False
        pd.testing.assert_frame_equal(actual_false, compare)


def test_single_noncrux_data(test_single_noncrux_dataset_class):
    """
    Checks whether or not a PsmDataset object is created properly
    after reading in a single tab delimited file with non crux
    column names. The point of this test is to make sure the user
    can supply an input file with differing column names by specifying
    the names of the new column names.

    Parameters
    ----------
    test_single_noncrux_dataset_class : pytest fixture of a PsmDataset object
        A psm object created by reading in the "single_noncrux.csv" file

    Returns
    -------
    Pandas Assert Frame
        Asserts whether or not the data in PsmDataset object
        is equal to the dataframe named "testframe_single_noncrux_data"
    """
    actual = testframe_single_noncrux_data.copy()
    compare = test_single_noncrux_dataset_class.data
    pd.testing.assert_frame_equal(actual, compare)


def test_single_noncrux_tdc(test_single_noncrux_dataset_class):
    """
    Checks whether or not a Result object is created properly
    after executing the calculate_tdc method on a PsmDataset object
    from test_single_noncrux_data.

    Parameters
    ----------
    test_single_noncrux_dataset_class : pytest fixture of a PsmDataset object
        A psm object created by reading in the "single_noncrux.csv" file

    Returns
    -------
    Pandas Assert Frame
        Asserts whether or not the data in Result object
        is equal to the dataframe named "testframe_single_noncrux_tdc"
    """
    actual = testframe_single_noncrux_tdc.copy()
    output = calculate_tdc(test_single_noncrux_dataset_class)
    compare = output.data
    pd.testing.assert_frame_equal(actual, compare)


def test_multi_data(test_multi_dataset_class):
    """
    Checks whether or not a PsmDataset object is created properly
    after reading in multiple (two) tab delimited files with crux default
    column names. The point of this test is to ensure that multiple files
    can be combined to perform analysis.

    Parameters
    ----------
    test_multi_dataset_class : pytest fixture of a PsmDataset object
        A psm object created by reading in the "multi_target.csv" file
        and the "multi_decoy.csv" file.

    Returns
    -------
    Pandas Assert Frame
        Asserts whether or not the data in PsmDataset object
        is equal to the dataframe named "testframe_multi_data"
    """
    actual = testframe_multi_data.copy()
    compare = test_multi_dataset_class.data
    pd.testing.assert_frame_equal(actual, compare)


def test_multi_tdc(test_multi_dataset_class):
    """
    Checks whether or not a Result object is created properly
    after executing the calculate_tdc method on a PsmDataset object
    from test_multi_data.

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
    compare = output.data
    pd.testing.assert_frame_equal(actual, compare)
