# -*- coding: utf-8 -*-
"""
Created on Sat Apr  3 15:42:40 2021

@author: veenstra

"""

import pytest
import os
import glob
import subprocess

# ACCEPTANCE TESTS VIA EXAMPLE SCRIPTS, these are the ones who are only meant to generate output files

dir_tests = os.path.dirname(__file__) #F9 doesnt work, only F5 (F5 also only method to reload external definition scripts)
list_configfiles = glob.glob(os.path.join(dir_tests,'examples','*.py'))
dir_output_general = os.path.join(dir_tests,'examples_output')
os.makedirs(dir_output_general, exist_ok=True)

@pytest.mark.acceptance
@pytest.mark.parametrize("file_config", [pytest.param(file_config, id=os.path.basename(file_config).replace('.py','')) for file_config in list_configfiles])
def test_run_examples(file_config):
    # 1. Set up test data
    os.chdir(dir_output_general)

    p = subprocess.Popen(f"hatyan {file_config} --redirect-stdout --overwrite",
                          stderr=subprocess.STDOUT, # Merge stdout and stderr
                          stdout=subprocess.PIPE,
                          shell=True)
    # max 3 minutes per test, if it hangs longer the test is killed
    p.wait(180)
    
    if p.returncode:
        out, err = p.communicate()
        out_str = "\n".join([x.decode("utf-8") for x in out.split(b'\r\n')])
        raise RuntimeError(out_str)