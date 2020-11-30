"""
These tests verify that the CLI works as expected.
"""

import os
import subprocess

FILES = [
    os.path.join("data", f)
    for f in os.listdir("data")
    if f.startswith("single_")
]


def test_basic_cli(tmp_path):
    """Test that basic cli works."""
    cmd = ["crema", "data/single_basic.csv", "--output_dir", tmp_path]
    subprocess.run(cmd, check=True)
    assert os.path.isfile(os.path.join(tmp_path, "crema.psm_results.txt"))
    assert os.path.isfile(os.path.join(tmp_path, "crema.logfile.log"))


def test_cli_custom_root(tmp_path):
    """Test that basic cli works."""
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
    """Test that cli works with non-standard crux parameters."""
    cmd = [
        "crema",
        "data/single.csv",
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
    """Test that cli works with non-standard crux parameters and targets identified via numbers rather than boolean."""
    cmd = [
        "crema",
        "data/single_int_targets.csv",
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


def test_cli_text_scan(tmp_path):
    """Test that cli works with non-standard crux parameters and scan column containing Strings instead of ints."""
    cmd = [
        "crema",
        "data/single_text_scan.csv",
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


def test_cli_multi(tmp_path):
    """Test that cli works with non-standard crux parameters and multiple file inputs."""
    cmd = [
        "crema",
        "data/multi_target.csv",
        "data/multi_decoy.csv",
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
