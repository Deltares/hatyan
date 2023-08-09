# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 12:14:02 2023

@author: veenstra
"""

import os
import glob
import hatyan

dir_tests = os.path.dirname(__file__) #F9 doesnt work, only F5 (F5 also only method to reload external definition scripts)
dir_testdata = os.path.join(dir_tests,'data_unitsystemtests')

def test_ts_from_multifile_equidistant_dia_hasfreq():
    """
    When reading multiple equidistant diafiles that combine into a continuous timeseries,
    there should be a freq attribute that is not None
    """
    file_ts_pat = os.path.join(dir_testdata,'VLISSGN_obs?.txt')
    file_ts = glob.glob(file_ts_pat)
    
    ts_pd = hatyan.readts_dia(filename=file_ts)
    
    ts_freq_min = ts_pd.index.freq.nanos/1e9/60
    
    assert hasattr(ts_pd.index,'freq')
    assert ts_pd.index.freq is not None
    assert ts_freq_min==60
