# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 14:33:40 2023

@author: veenstra
"""

import os
import pytest
import hatyan
from hatyan.metadata import metadata_from_obj

dir_tests = os.path.dirname(__file__) #F9 doesnt work, only F5 (F5 also only method to reload external definition scripts)
dir_testdata = os.path.join(dir_tests,'data_unitsystemtests')


@pytest.mark.systemtest
def test_components_timeshift():
    
    timeshift_hr = 1
    current_station = 'HOEKVHLD'
    file_data_comp0 = os.path.join(dir_testdata,f'{current_station}_ana.txt')
    
    comp_merged = hatyan.read_components(filename=file_data_comp0)
    
    comp_shift = hatyan.components_timeshift(comp_merged,hours=timeshift_hr)
    
    # assert shift
    assert comp_shift.loc['SA','phi_deg'] == 221.54106863959504
    
    # assert metadata
    comp_merged_meta = metadata_from_obj(comp_merged)
    comp_shift_meta = metadata_from_obj(comp_shift)
    
    comp_merged_tz_min = comp_merged_meta['tzone']._minutes
    comp_shift_tz_min = comp_shift_meta['tzone']._minutes
    
    assert comp_merged_tz_min + timeshift_hr*60 == comp_shift_tz_min
    