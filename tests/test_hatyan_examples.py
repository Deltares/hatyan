# -*- coding: utf-8 -*-
"""
Created on Sat Apr  3 15:42:40 2021

@author: veenstra

"""


import pytest
import os
import glob

dir_tests = os.path.dirname(__file__) #F9 doesnt work, only F5 (F5 also only method to reload external definition scripts)
dir_output_general = os.path.join(dir_tests,'examples_output')
if not os.path.exists(dir_output_general):
    os.mkdir(dir_output_general)


# High level acceptance tests, these are the ones who are only meant to generate output files
# for the testers to verify (in Teamcity) whether the runs generate the expected files or not.
""" Run hatyan_main.py with test-configfiles as input """
list_examplescripts = glob.glob(os.path.join(dir_tests,'examples','*.py'))
list_examplescripts = [x for x in list_examplescripts if '_interactive' not in x]
#list_examplescripts = ['predictie_2019_b02ex2_19Ycomp4Ydia_CUXHVN_test.py']


@pytest.mark.acceptance
@pytest.mark.parametrize("file_config", [pytest.param(file_config, id=os.path.basename(file_config).replace('.py','')) for file_config in list_examplescripts])
def test_examplescripts(file_config): #FROM DFM_TOOLS
    """
    file_config = os.path.join(dir_tests,'configfiles','predictie_2019_b02ex2_19Ycomp4Ydia_CUXHVN_test.py')
    """
    # 1. Set up test data
    dir_output = os.path.join(dir_output_general,os.path.basename(file_config).replace('.py',''))
    if not os.path.exists(dir_output):
        os.mkdir(dir_output)
    os.chdir(dir_output)
    test = os.system(f"python {file_config} > {dir_output}/STDOUT.txt")
    #test = os.system(f"python {file_config}")
    
    if test:
        raise Exception('execution did not finish properly')



