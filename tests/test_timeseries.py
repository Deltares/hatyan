# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 12:14:02 2023

@author: veenstra
"""

import os
import pytest
import glob
import pandas as pd
import hatyan

dir_tests = os.path.dirname(__file__) #F9 doesnt work, only F5 (F5 also only method to reload external definition scripts)
dir_testdata = os.path.join(dir_tests,'data_unitsystemtests')


@pytest.mark.unittest
def test_ts_from_singlefile_equidistant_dia_hasfreq():
    """
    When reading asingle equidistant diafile,
    there should be a freq attribute that is not None
    """
    file_ts = os.path.join(dir_testdata,'VLISSGN_obs1.txt')
    
    ts_pd = hatyan.readts_dia(filename=file_ts)
    
    # check ts length (all four files are added)
    assert len(ts_pd) == 8760
    
    # assert on freq attribute
    assert hasattr(ts_pd.index,'freq')
    assert isinstance(ts_pd.index.freq,pd.offsets.Minute)
    assert ts_pd.index.freq is not None
    assert ts_pd.index.freq.nanos/1e9 == 3600


@pytest.mark.unittest
def test_ts_from_multifile_equidistant_dia_hasfreq():
    """
    When reading multiple equidistant diafiles that combine into a continuous timeseries,
    there should be a freq attribute that is not None
    """
    file_ts_pat = os.path.join(dir_testdata,'VLISSGN_obs?.txt')
    file_ts = glob.glob(file_ts_pat)
    
    ts_pd = hatyan.readts_dia(filename=file_ts)
    
    # check ts length (all four files are added)
    assert len(ts_pd) == 35064
    
    # assert on freq attribute
    assert hasattr(ts_pd.index,'freq')
    assert isinstance(ts_pd.index.freq,pd.offsets.Minute)
    assert ts_pd.index.freq is not None
    assert ts_pd.index.freq.nanos/1e9 == 3600


@pytest.mark.unittest
def test_ts_from_multifile_equidistant_dia_correctglob():
    """
    When providing a file pattern for reading multiple equidistant diafiles,
    glob is supported, this test checks if it results in the correct data
    """
    file_ts = os.path.join(dir_testdata,'VLISSGN_obs?.txt')
    
    ts_pd = hatyan.readts_dia(filename=file_ts)
    
    assert len(ts_pd) == 35064
