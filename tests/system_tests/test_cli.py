"""These tests verify that the crema CLI works as expected."""
from pathlib import Path
import subprocess


def test_cli_basic(basic_crux_txt, tmp_path):
    """Test that the basic cli works."""
    cmd = ["crema", "--output_dir", tmp_path, "-e", "0.5", basic_crux_txt]
    subprocess.run(cmd, check=True)
    assert Path(tmp_path, "crema.psms.txt").exists()
    assert Path(tmp_path, "crema.log.txt").exists()


def test_cli_custom_root(basic_crux_txt, tmp_path):
    """Test that the cli works with custom file root."""
    cmd = ["crema", "-o", tmp_path, "-f", "myFileRoot", "-e", "0.5"]
    cmd.append(basic_crux_txt)
    subprocess.run(cmd, check=True)
    assert Path(tmp_path, "myFileRoot.crema.psms.txt").exists()
    assert Path(tmp_path, "myFileRoot.crema.log.txt").exists()


def test_real(real_crux_txt, tmp_path):
    """Test that crema works on real Crux results"""
    cmd = ["crema", "--output_dir", tmp_path] + list(real_crux_txt)
    subprocess.run(cmd, check=True)
    assert Path(tmp_path, "crema.psms.txt").exists()
    assert Path(tmp_path, "crema.log.txt").exists()
