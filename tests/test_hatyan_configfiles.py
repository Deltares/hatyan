# -*- coding: utf-8 -*-
"""
Created on Sat Apr  3 15:42:40 2021

@author: veenstra

"""


import pytest
import os
import glob
import hatyan

dir_scriptfile = os.path.realpath(__file__) #F9 doesnt work, only F5 (F5 also only method to reload external definition scripts)
dir_tests = os.path.abspath(os.path.join(dir_scriptfile,os.pardir))  #1 level up from dir_scripts

dir_output = os.path.join(dir_tests,'output_configfiles')
if not os.path.exists(dir_output):
    os.mkdir(dir_output)
os.chdir(dir_output)


# High level acceptance tests, these are the ones who are only meant to generate output files
# for the testers to verify (in Teamcity) whether the runs generate the expected files or not.
""" Run hatyan_main.py with test-configfiles as input """
list_configfiles = glob.glob(os.path.join(dir_tests,'configfiles','*.py'))
list_configfiles = [x for x in list_configfiles if '_interactive' not in x]
#list_configfiles = ['predictie_2019_b02ex2_19Ycomp4Ydia_CUXHVN_test.py']


@pytest.mark.acceptance
@pytest.mark.parametrize("file_config", [pytest.param(file_config, id=os.path.basename(file_config).replace('.py','')) for file_config in list_configfiles])
def test_configfiles(file_config):
    """
    file_config = os.path.join(dir_tests,'configfiles','predictie_2019_b02ex2_19Ycomp4Ydia_CUXHVN_test.py')
    """
    # 1. Set up test data
    dir_output = hatyan.get_outputfoldername(file_config)
    if not os.path.exists(dir_output):
        os.mkdir(dir_output)

    os.system("python {0} {1} > {1}/FILE_DIAGNOSTICS.txt 2>&1".format(file_config, dir_output))#+ " & pause")
    
    if os.path.exists(os.path.join(dir_output,'__NOT_FINISHED__')):
        file_diag = os.path.join(dir_output,'FILE_DIAGNOSTICS.txt')
        with open(file_diag,'r') as f:
            for i, l in enumerate(f):
                pass
        with open(file_diag,'r') as f:
            diag_contents = f.read().split('\n')[i-39:i+2]
        diag_contents_1str = '\n'.join(diag_contents)
        raise Exception('calculation did not finish properly since __NOT_FINISHED__ file is still present in output directory: %s\nlast 40 lines of diagnostics file:\n\n%s'%(dir_output, diag_contents_1str))


