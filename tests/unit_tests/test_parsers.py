"""
These are unit tests for functions within parsers.py:
"""
import os
import pytest

# import ppx
from crema import *


@pytest.fixture
def tmp_file(tmp_path):
    """
    Creates a pytest fixture of temporary filepath
    used to write/read test files

    Parameters
    ----------
    tmp_path : pytest fixture of a temporary directory path object

    Returns
    -------
    : str
        The filepath to be used to create the testing file
    """
    return os.path.join(tmp_path, "test.txt")


@pytest.fixture
def dataframe_basic():
    """
    Creates a pytest fixture of a pandas DataFrame
    to be converted into a csv file to test read_file().
    Used with test_read_file_basic().

    Returns
    -------
    pandas.DataFrame
        A pandas.DataFrame containing psm data
    """
    return pd.DataFrame(
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


@pytest.fixture()
def dataframe_text_spectrum():
    """
    Creates a pytest fixture of a pandas DataFrame
    to be converted into a csv file to test read_file().
    Used with test_read_file_text_spectrum().

    Returns
    -------
    pandas.DataFrame
        A pandas.DataFrame containing psm data
    """
    return pd.DataFrame(
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


@pytest.fixture
def dataframe_add_spectrum():
    """
    Creates a pytest fixture of a pandas DataFrame
    to be converted into a csv file to test read_file().
    Used with test_read_file_add_spectrum().

    Returns
    -------
    pandas.DataFrame
        A pandas.DataFrame containing psm data
    """
    return pd.DataFrame(
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
            "extras": ["a", "b", "c", "d", "e", "z", "b", "y", "d", "x",],
        }
    )


@pytest.fixture
def dataframe_add_score():
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
def dataframe_non_crux():
    """
    Creates a pytest fixture of a pandas DataFrame
    to be converted into a csv file to test read_file().
    Used with test_read_file_non_crux().

    Returns
    -------
    pandas.DataFrame
        A pandas.DataFrame containing psm data
    """
    return pd.DataFrame(
        {
            "spectra_ref": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],
            "p-score": [0.7, 0.4, 0.1, 0.55, 0.25, 0.6, 0.2, 0.7, 0.56, 0.3,],
            "decoy_ind": [
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
def dataframe_multi_file_1():
    """
    Creates a pytest fixture of a pandas DataFrame
    to be converted into a csv file to test read_file().
    Used with test_read_file_multi_file().

    Returns
    -------
    pandas.DataFrame
        A pandas.DataFrame containing psm data
    """
    return pd.DataFrame(
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
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
            ],
        }
    )


@pytest.fixture
def dataframe_multi_file_2():
    """
    Creates a pytest fixture of a pandas DataFrame
    to be converted into a csv file to test read_file().
    Used with test_read_file_multi_file().

    Returns
    -------
    pandas.DataFrame
        A pandas.DataFrame containing psm data
    """
    return pd.DataFrame(
        {
            "scan": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],
            "combined p-value": [
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


def test_read_file_basic(dataframe_basic, tmp_file):
    """
    Ensures that a PsmDataset object is created properly
    after reading in a single file in its most basic
    form with crux default column names.

    Parameters
    ----------
    dataframe_basic : pytest fixture of a PsmDataset object
        The actual dataframe used for this test case
    tmp_file : str
        The filepath to be used to create the testing file
    """
    verify_dataset(
        act_dataframe=dataframe_basic,
        tmp_file=tmp_file,
        act_spectrum_col=["scan"],
        act_score_col=["combined p-value"],
        act_target_col="target/decoy",
    )


def test_read_file_text_spectrum(dataframe_text_spectrum, tmp_file):
    """
    Ensures that a PsmDataset object is created properly
    after reading in a single file with
    spectrum column values containing strings

    Parameters
    ----------
    dataframe_text_spectrum : pytest fixture of a PsmDataset object
        The actual dataframe used for this test case
    tmp_file : str
        The filepath to be used to create the testing file
    """
    verify_dataset(
        act_dataframe=dataframe_text_spectrum,
        tmp_file=tmp_file,
        act_spectrum_col=["scan"],
        act_score_col=["combined p-value"],
        act_target_col="target/decoy",
    )


def test_read_file_add_spectrum(dataframe_add_spectrum, tmp_file):
    """
    Ensures that a PsmDataset object is created properly
    after reading in a single file with an additional
    spectrum column

    Parameters
    ----------
    dataframe_add_spectrum : pytest fixture of a PsmDataset object
        The actual dataframe used for this test case
    tmp_file : str
        The filepath to be used to create the testing file
    """
    verify_dataset(
        act_dataframe=dataframe_add_spectrum,
        tmp_file=tmp_file,
        act_spectrum_col=["scan", "extras"],
        act_score_col=["combined p-value"],
        act_target_col="target/decoy",
    )


def test_read_file_add_score(dataframe_add_score, tmp_file):
    """
    Ensures that a PsmDataset object is created properly
    after reading in a single file with a few additional
    score columns

    Parameters
    ----------
    dataframe_add_score : pytest fixture of a PsmDataset object
        The actual dataframe used for this test case
    tmp_file : str
        The filepath to be used to create the testing file
    """
    verify_dataset(
        act_dataframe=dataframe_add_score,
        tmp_file=tmp_file,
        act_spectrum_col=["scan"],
        act_score_col=[
            "combined p-value 0",
            "combined p-value 1",
            "combined p-value 2",
        ],
        act_target_col="target/decoy",
    )


def test_read_file_non_crux(dataframe_non_crux, tmp_file):
    """
    Ensures that a PsmDataset object is created properly
    after reading in a single file with non crux default
    column names.

    Parameters
    ----------
    dataframe_non_crux : pytest fixture of a PsmDataset object
        The actual dataframe used for this test case
    tmp_file : str
        The filepath to be used to create the testing file
    """
    verify_dataset(
        act_dataframe=dataframe_non_crux,
        tmp_file=tmp_file,
        act_spectrum_col=["spectra_ref"],
        act_score_col=["p-score"],
        act_target_col="decoy_ind",
    )


def test_read_file_multi_file(
    dataframe_multi_file_1, dataframe_multi_file_2, tmp_path
):
    """
    Ensures that a PsmDataset object is created properly
    after reading in two separate files

    Parameters
    ----------
    dataframe_multi_file_1: pytest fixture of a PsmDataset object
        The actual dataframe used for this test case
    dataframe_multi_file_2: pytest fixture of a PsmDataset object
        The actual dataframe used for this test case
    tmp_path : str
        The filepath to be used to create the testing file
    """
    path1 = os.path.join(tmp_path, "test1.txt")
    path2 = os.path.join(tmp_path, "test2.txt")
    # Create a temporary file from the given dataframe
    dataframe_multi_file_1.to_csv(path1)
    dataframe_multi_file_2.to_csv(path2)
    # Create a PsmDataset object by calling the read_file() function
    psm_dataset = read_file(
        [path1, path2],
        spectrum_col=["scan"],
        score_col=["combined p-value"],
        target_col="target/decoy",
    )
    # Asserts that the PsmDataset object is an instance of the PsmDataset class
    assert isinstance(psm_dataset, PsmDataset)
    # Asserts that the data saved in the PsmDataset object is as expected
    expected_frame = dataframe_multi_file_1.append(
        dataframe_multi_file_2, ignore_index=True
    )
    pd.testing.assert_frame_equal(psm_dataset.data, expected_frame)
    # Asserts that the spectrum column name saved in the PsmDataset object is as expected
    assert psm_dataset.spectrum_col == ["scan"]
    # Asserts that the score column name saved in the PsmDataset object is as expected
    assert psm_dataset.score_col == ["combined p-value"]
    # Asserts that the target column name saved in the PsmDataset object is as expected
    assert psm_dataset.target_col == "target/decoy"


def verify_dataset(
    act_dataframe, tmp_file, act_spectrum_col, act_score_col, act_target_col,
):
    """
    Helper method to speed up test writing.
    Asserts that a PsmDataset object is properly created
    with the correct fields and attributes.

    Parameters
    ----------
    act_dataframe : pandas.Dataframe
        The actual dataframe that is being used to create the testing file
    tmp_file : str
        The filepath to be used to create the testing file
    act_spectrum_col : tuple of str
        The actual column name(s) of the spectrum column used in the testing file
    act_score_col : tuple of str
        The actual column name(s) of the score column used in the testing file
    act_target_col : str
        The actual column name of the target column used in the testing file
    """
    # Create a temporary file from the given dataframe
    act_dataframe.to_csv(tmp_file)
    # Create a PsmDataset object by calling the read_file() function
    psm_dataset = read_file(
        tmp_file,
        spectrum_col=act_spectrum_col,
        score_col=act_score_col,
        target_col=act_target_col,
    )
    # Asserts that the PsmDataset object is an instance of the PsmDataset class
    assert isinstance(psm_dataset, PsmDataset)
    # Asserts that the data saved in the PsmDataset object is as expected
    pd.testing.assert_frame_equal(psm_dataset.data, act_dataframe)
    # Asserts that the spectrum column name saved in the PsmDataset object is as expected
    assert psm_dataset.spectrum_col == act_spectrum_col
    # Asserts that the score column name saved in the PsmDataset object is as expected
    assert psm_dataset.score_col == act_score_col
    # Asserts that the target column name saved in the PsmDataset object is as expected
    assert psm_dataset.target_col == act_target_col


def download_msv(msv_id, tmp_path):
    dat = ppx.MSVDataset(msv_id)
    mztab_file = dat.list_files("ccms_result")[0]
    file_path = "ccms_result/" + mztab_file
    file = dat.download(files=file_path, dest_dir=tmp_path)[0]
    return file


def verify_mztab(msv_id, num_scores, expected_df_head, tmp_path):
    file = download_msv(msv_id, tmp_path)
    psm_dataset = read_mztab(file)
    # result = calculate_tdc(psm_dataset, 2)
    # print(result.data)
    assert isinstance(psm_dataset, PsmDataset)
    assert psm_dataset.spectrum_col == ["spectra_ref"]
    score_col = []
    for i in range(1, num_scores + 1):
        score_col.append("search_engine_score[" + str(i) + "]")
    assert psm_dataset.score_col == score_col
    assert psm_dataset.target_col == "opt_global_cv_MS:1002217_decoy_peptide"
    pd.testing.assert_frame_equal(psm_dataset.data.head(), expected_df_head)


# def test_read_mztab1(tmp_path):
#     expected_df_head = pd.DataFrame(
#         {
#             "spectra_ref": [
#                 "ms_run[1]:index=1492",
#                 "ms_run[1]:index=1493",
#                 "ms_run[2]:index=2260",
#                 "ms_run[2]:index=1740",
#                 "ms_run[1]:index=335",
#             ],
#             "search_engine_score[1]": [
#                 0.9999993,
#                 0.99999905,
#                 0.982506,
#                 0.99364066,
#                 0.9983276,
#             ],
#             "search_engine_score[2]": [
#                 5.20000089003336e-08,
#                 6.899999633797e-08,
#                 0.000249999949201165,
#                 0.0000590000104995094,
#                 0.000230000061930564,
#             ],
#             "opt_global_cv_MS:1002217_decoy_peptide": [
#                 True,
#                 True,
#                 True,
#                 True,
#                 True,
#             ],
#         }
#     )
#     verify_mztab("MSV000085943", 2, expected_df_head, tmp_path)
#     expected_df_head = pd.DataFrame(
#         {
#             "spectra_ref": [
#                 "ms_run[1]:index=111",
#                 "ms_run[1]:index=112",
#                 "ms_run[1]:index=448",
#                 "ms_run[1]:index=449",
#                 "ms_run[1]:index=790",
#             ],
#             "search_engine_score[1]": [
#                 28.5,
#                 29.02,
#                 60.55,
#                 60.04,
#                 103.69,
#             ],
#             "search_engine_score[2]": [
#                 25.0,
#                 25.0,
#                 25.0,
#                 25.0,
#                 25.0,
#             ],
#             "search_engine_score[3]": [
#                 np.nan,
#                 0.0109999999702295,
#                 0.000180000050470909,
#                 0.000659999852526374,
#                 0.0000430000170645866,
#             ],
#             "search_engine_score[4]": [
#                 0.730391,
#                 0.99606955,
#                 0.99999803,
#                 0.99999774,
#                 1,
#             ],
#             "opt_global_cv_MS:1002217_decoy_peptide": [
#                 True,
#                 True,
#                 True,
#                 True,
#                 True,
#             ],
#         }
#     )
#     verify_mztab("MSV000085729", 4, expected_df_head, tmp_path)
#     expected_df_head = pd.DataFrame(
#         {
#             "spectra_ref": [
#                 "ms_run[1]:index=6897",
#                 "ms_run[1]:index=6019",
#                 "ms_run[1]:index=1708",
#                 "ms_run[1]:index=4268",
#                 "ms_run[1]:index=2517",
#             ],
#             "search_engine_score[1]": [
#                 37.31,
#                 57.93,
#                 38.93,
#                 40.83,
#                 52.78,
#             ],
#             "search_engine_score[2]": [
#                 0.9998933,
#                 1,
#                 0.9999762,
#                 0.99999636,
#                 1,
#             ],
#             "opt_global_cv_MS:1002217_decoy_peptide": [
#                 True,
#                 True,
#                 True,
#                 True,
#                 True,
#             ],
#         }
#     )
#     verify_mztab("MSV000085211", 2, expected_df_head, tmp_path)
