# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 21:47:33 2023

@author: veenstra
"""

import os
import pytest
import hatyan
from hatyan.metadata import metadata_from_diablocks, metadata_add_to_obj, metadata_from_obj, metadata_compare

dir_tests = os.path.dirname(__file__) #F9 doesnt work, only F5 (F5 also only method to reload external definition scripts)
dir_testdata = os.path.join(dir_tests,'data_unitsystemtests')



@pytest.mark.unittest
def readts_dia_metadata_multifile():
    
    current_station = 'VLISSGN'
    file_ts = [os.path.join(dir_testdata, f'{current_station}_obs{i}.txt') for i in [1,2,3,4]]
    ts_measurements_group0 = hatyan.readts_dia(filename=file_ts, station=current_station)
    meta_fromts = metadata_from_obj(ts_measurements_group0)
    
    meta_expected = {
        'vertref': 'NAP',
        'station': 'VLISSGN',
        'TYP': 'TE',
        'groepering': 'NVT',
        'grootheid': 'WATHTE',
        'eenheid': 'cm',
        'timestep_min': 60.0,
        'timestep_unit': 'min'}
    
    assert meta_fromts == meta_expected

@pytest.mark.unittest
def readts_dia_metadata_multiblock():
    file_ts = r'c:\DATA\hatyan_github\tests\data_unitsystemtests\hoek_har.dia'
    ts_measurements_group0 = hatyan.readts_dia(filename=file_ts, station='HOEKVHLD', block_ids='allstation')
    meta_fromts = metadata_from_obj(ts_measurements_group0)
    
    meta_expected = {
        'vertref': 'NAP',
        'station': 'HOEKVHLD',
        'TYP': 'TN',
        'groepering': 'GETETM2',
        'grootheid': 'WATHTE',
        'eenheid': 'cm',
        'timestep_min': None,
        'timestep_unit': None}
    
    assert meta_fromts == meta_expected
