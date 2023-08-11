# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 21:47:33 2023

@author: veenstra
"""

import os
import pytest
import pytz
import pandas as pd
import hatyan
from hatyan.metadata import metadata_from_obj, metadata_compare, wns_from_metadata #, metadata_from_diablocks, metadata_add_to_obj

dir_tests = os.path.dirname(__file__) #F9 doesnt work, only F5 (F5 also only method to reload external definition scripts)
dir_testdata = os.path.join(dir_tests,'data_unitsystemtests')


@pytest.mark.unittest
def test_readts_dia_metadata_multifile():
    
    current_station = 'VLISSGN'
    file_ts = [os.path.join(dir_testdata, f'{current_station}_obs{i}.txt') for i in [1,2,3,4]]
    ts_measurements_group0 = hatyan.readts_dia(filename=file_ts, station=current_station)
    meta_fromts = metadata_from_obj(ts_measurements_group0)
    
    meta_expected = {'station': 'VLISSGN',
     'grootheid': 'WATHTE',
     'eenheid': 'cm',
     'vertref': 'NAP',
     'tstart': pd.Timestamp('2009-01-01 00:00:00'),
     'tstop': pd.Timestamp('2012-12-31 23:00:00'),
     'timestep_min': 60.0,
     'timestep_unit': 'min',
     'TYP': 'TE',
     'groepering': 'NVT',
     'tzone': pytz.FixedOffset(60),
     'origin': 'from timeseries dia file'}
    
    assert meta_fromts == meta_expected


@pytest.mark.unittest
def test_readts_dia_metadata_multiblock():
    file_ts = os.path.join(dir_testdata, 'hoek_har.dia')
    ts_measurements_group0 = hatyan.readts_dia(filename=file_ts, station='HOEKVHLD', block_ids='allstation')
    meta_fromts = metadata_from_obj(ts_measurements_group0)
    
    meta_expected = {'station': 'HOEKVHLD',
     'grootheid': 'WATHTE',
     'eenheid': 'cm',
     'vertref': 'NAP',
     'tstart': pd.Timestamp('1980-01-01 01:32:00'),
     'tstop': pd.Timestamp('1991-12-31 23:45:00'),
     'timestep_min': None,
     'timestep_unit': None,
     'TYP': 'TN',
     'groepering': 'GETETM2',
     'tzone': pytz.FixedOffset(60),
     'origin': 'from timeseries dia file'}
    
    assert meta_fromts == meta_expected


@pytest.mark.unittest
def test_metadata_compare():
    metadata = {
        'vertref': 'NAP',
        'station': 'VLISSGN',
        'TYP': 'TE',
        'groepering': 'NVT',
        'grootheid': 'WATHTE',
        'eenheid': 'cm',
        'timestep_min': 60.0,
        'timestep_unit': 'min',
        'tstart': None,
        'tstop': None}
    
    metadata_compare([metadata,metadata,metadata])


@pytest.mark.unittest
def test_wns_from_metadata():
    metadata_1 = {'grootheid':'WATHTE', 'eenheid':'cm', 'vertref':'NAP'}
    metadata_54 = {'grootheid':'WATHTE', 'eenheid':'cm', 'vertref':'MSL'}
    metadata_18 = {'grootheid':'WATHTBRKD', 'eenheid':'cm', 'vertref':'NAP'}
    metadata_55 = {'grootheid':'WATHTBRKD', 'eenheid':'cm', 'vertref':'MSL'}
    
    wns_1 = wns_from_metadata(metadata_1)
    wns_54 = wns_from_metadata(metadata_54)
    wns_18 = wns_from_metadata(metadata_18)
    wns_55 = wns_from_metadata(metadata_55)
    
    assert wns_1 == 1
    assert wns_54 == 54
    assert wns_18 == 18
    assert wns_55 == 55
    
    # metadata_1 = {'grootheid':'GETETM2', 'eenheid':'cm', 'vertref':'NAP'}
    # metadata_54 = {'grootheid':'GETETM2', 'eenheid':'cm', 'vertref':'MSL'}
    # metadata_18 = {'grootheid':'GETETBRKD2', 'eenheid':'cm', 'vertref':'NAP'}
    # metadata_55 = {'grootheid':'GETETBRKD2', 'eenheid':'cm', 'vertref':'MSL'}
    
    # wns_1 = wns_from_metadata(metadata_1)
    # wns_54 = wns_from_metadata(metadata_54)
    # wns_18 = wns_from_metadata(metadata_18)
    # wns_55 = wns_from_metadata(metadata_55)
    
    # assert wns_1 == 1
    # assert wns_54 == 54
    # assert wns_18 == 18
    # assert wns_55 == 55
