# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 09:25:17 2024

@author: veenstra
"""

import os
import pytest
from click.testing import CliRunner
from hatyan import cli
import shutil
import glob

dir_tests = os.path.dirname(__file__) #F9 doesnt work, only F5 (F5 also only method to reload external definition scripts)

@pytest.mark.unittest
def test_command_line_interface(tmp_path):
    """Test the CLI."""
    os.chdir(tmp_path)
    
    filename = os.path.join(dir_tests,'examples','minimal_example.py')
    assert os.path.exists(filename)
    
    runner = CliRunner()
    
    result = runner.invoke(cli.cli)
    assert result.exit_code == 2
    assert 'Missing argument' in result.output
    
    help_result = runner.invoke(cli.cli, ['--help'])
    assert help_result.exit_code == 0
    assert 'Show this message and exit.' in help_result.output
    
    dir_output = os.path.abspath(os.path.basename(filename).replace(".py",""))
    if os.path.exists(dir_output):
        shutil.rmtree(dir_output)
    filename_result = runner.invoke(cli.cli, [filename, '--redirect-stdout'])
    assert filename_result.exit_code == 0
    
    # assert file presence
    file_stdout = os.path.join(dir_output, "STDOUT.txt")
    file_script = os.path.join(dir_output, os.path.basename(filename))
    dia_list = glob.glob(os.path.join(dir_output, "*.dia"))
    png_list = glob.glob(os.path.join(dir_output, "*.png"))
    
    assert os.path.exists(dir_output)
    assert os.path.exists(file_stdout)
    assert os.path.exists(file_script)
    assert len(dia_list) > 0
    assert len(png_list) > 0
