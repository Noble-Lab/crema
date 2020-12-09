"""
These tests verify that the crema CLI works as expected.
"""

import os
import subprocess


def test_cli_basic(tmp_path):
    """
    Test that the basic cli works. Reads in a file
    with crux default column names. Does nothing else.

    Parameters
    ----------
    tmp_path : pytest fixture of a temporary directory
        A pytest temporary directory unique to the test invocation

    Returns
    -------
    Pandas Assert Frame
        Asserts whether or not the the results file and log file
        are created properly with the correct file path.
    """
    cmd = ["crema", "data/single_basic.csv", "--output_dir", tmp_path]
    subprocess.run(cmd, check=True)
    assert os.path.isfile(os.path.join(tmp_path, "crema.psm_results.txt"))
    assert os.path.isfile(os.path.join(tmp_path, "crema.logfile.log"))


def test_cli_basic_tab(tmp_path):
    """
    Test that the cli works with a tab delimited file. Reads in a file
    with crux default column names.

    Parameters
    ----------
    tmp_path : pytest fixture of a temporary directory
        A pytest temporary directory unique to the test invocation

    Returns
    -------
    Pandas Assert Frame
        Asserts whether or not the the results file and log file
        are created properly with the correct file path.
    """
    cmd = ["crema", "data/single_basic_tab.txt", "--output_dir", tmp_path]
    subprocess.run(cmd, check=True)
    assert os.path.isfile(os.path.join(tmp_path, "crema.psm_results.txt"))
    assert os.path.isfile(os.path.join(tmp_path, "crema.logfile.log"))


def test_cli_custom_root(tmp_path):
    """
    Test that the cli works with custom file root.
    Reads in a file with crux default column names.
    Specifies a custom file_root string to tag on
    as a prefix to the output files.

    Parameters
    ----------
    tmp_path : pytest fixture of a temporary directory
        A pytest temporary directory unique to the test invocation

    Returns
    -------
    Pandas Assert Frame
        Asserts whether or not the the results file and log file
        are created properly with the correct file path.
    """
    cmd = [
        "crema",
        "data/single_basic.csv",
        "--output_dir",
        tmp_path,
        "--file_root",
        "myFileRoot",
    ]
    subprocess.run(cmd, check=True)
    assert os.path.isfile(
        os.path.join(tmp_path, "myFileRootcrema.psm_results.txt")
    )
    assert os.path.isfile(
        os.path.join(tmp_path, "myFileRootcrema.logfile.log")
    )


def test_cli_custom_param(tmp_path):
    """
    Test that the cli works when specifying non-default column
    names. Reads in a file with non-crux default column names.
    Specifies different score and target columns to read.

    Parameters
    ----------
    tmp_path : pytest fixture of a temporary directory
        A pytest temporary directory unique to the test invocation

    Returns
    -------
    Pandas Assert Frame
        Asserts whether or not the the results file and log file
        are created properly with the correct file path.
    """
    cmd = [
        "crema",
        "data/single_noncrux.csv",
        "--output_dir",
        tmp_path,
        "--score",
        "p-value",
        "--target",
        "target",
    ]
    subprocess.run(cmd, check=True)
    assert os.path.isfile(os.path.join(tmp_path, "crema.psm_results.txt"))
    assert os.path.isfile(os.path.join(tmp_path, "crema.logfile.log"))


def test_cli_int_targets(tmp_path):
    """
    Test that the cli works with integer value target
    column values. Reads in a file
    with crux default column names. Does nothing else.

    Parameters
    ----------
    tmp_path : pytest fixture of a temporary directory
        A pytest temporary directory unique to the test invocation

    Returns
    -------
    Pandas Assert Frame
        Asserts whether or not the the results file and log file
        are created properly with the correct file path.
    """
    cmd = [
        "crema",
        "data/single_int_targets.csv",
        "--output_dir",
        tmp_path,
    ]
    subprocess.run(cmd, check=True)
    assert os.path.isfile(os.path.join(tmp_path, "crema.psm_results.txt"))
    assert os.path.isfile(os.path.join(tmp_path, "crema.logfile.log"))


def test_cli_text_scan(tmp_path):
    """
    Test that the cli works with String value scan
    column values. Reads in a file
    with crux default column names. Does nothing else.

    Parameters
    ----------
    tmp_path : pytest fixture of a temporary directory
        A pytest temporary directory unique to the test invocation

    Returns
    -------
    Pandas Assert Frame
        Asserts whether or not the the results file and log file
        are created properly with the correct file path.
    """
    cmd = [
        "crema",
        "data/single_text_scan.csv",
        "--output_dir",
        tmp_path,
    ]
    subprocess.run(cmd, check=True)
    assert os.path.isfile(os.path.join(tmp_path, "crema.psm_results.txt"))
    assert os.path.isfile(os.path.join(tmp_path, "crema.logfile.log"))


def test_cli_add_spectrum(tmp_path):
    """
    Test that the cli works with more than one spectrum
    column. Reads in a file with crux default column names.
    Specifies two scan columns to read in.

    Parameters
    ----------
    tmp_path : pytest fixture of a temporary directory
        A pytest temporary directory unique to the test invocation

    Returns
    -------
    Pandas Assert Frame
        Asserts whether or not the the results file and log file
        are created properly with the correct file path.
    """
    cmd = [
        "crema",
        "data/single_add_spectrum.csv",
        "--spectrum",
        "scan",
        "extras",
        "--output_dir",
        tmp_path,
    ]
    subprocess.run(cmd, check=True)
    assert os.path.isfile(os.path.join(tmp_path, "crema.psm_results.txt"))
    assert os.path.isfile(os.path.join(tmp_path, "crema.logfile.log"))


def test_cli_multi(tmp_path):
    """
    Test that the cli works with multiple
    file inputs. Reads in two files
    with crux default column names. Does nothing else.

    Parameters
    ----------
    tmp_path : pytest fixture of a temporary directory
        A pytest temporary directory unique to the test invocation

    Returns
    -------
    Pandas Assert Frame
        Asserts whether or not the the results file and log file
        are created properly with the correct file path.
    """
    cmd = [
        "crema",
        "data/multi_target.csv",
        "data/multi_decoy.csv",
        "--output_dir",
        tmp_path,
    ]
    subprocess.run(cmd, check=True)
    assert os.path.isfile(os.path.join(tmp_path, "crema.psm_results.txt"))
    assert os.path.isfile(os.path.join(tmp_path, "crema.logfile.log"))
