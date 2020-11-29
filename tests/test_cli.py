"""
These tests verify that the CLI works as expected.
"""

import os
import subprocess

FILES = [
    os.path.join("data", f)
    for f in os.listdir("data")
    if f.startswith("single")
]


def test_basic_cli(tmp_path):
    """Test that basic cli works."""
    cmd = ["crema", FILES[1], "--output_dir", tmp_path]
    subprocess.run(cmd, check=True)
    assert os.path.isfile(os.path.join(tmp_path, "crema.psm_results.txt"))
    assert os.path.isfile(os.path.join(tmp_path, "crema.logfile.log"))
