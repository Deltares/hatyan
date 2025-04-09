# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 21:47:33 2023

@author: veenstra
"""

import os
import pytest
import pytz
import hatyan
from hatyan.metadata import metadata_from_obj, metadata_compare, wns_from_metadata

dir_tests = os.path.dirname(__file__) #F9 doesnt work, only F5 (F5 also only method to reload external definition scripts)
dir_testdata = os.path.join(dir_tests,'data_unitsystemtests')


@pytest.mark.unittest
def test_readts_dia_metadata_multifile():
    
    current_station = 'VLISSGN'
    file_ts = [os.path.join(dir_testdata, f'{current_station}_obs{i}.txt') for i in [1,2,3,4]]
    ts_measurements_group0 = hatyan.read_dia(filename=file_ts, station=current_station)
    meta_fromts = metadata_from_obj(ts_measurements_group0)
    
    meta_expected = {'station': 'VLISSGN',
     'grootheid': 'WATHTE',
     'eenheid': 'm',
     'vertref': 'NAP',
     'TYP': 'TE',
     'groepering': 'NVT',
     'origin': 'from timeseries dia file'}
    
    assert meta_fromts == meta_expected
    assert ts_measurements_group0.index.tz == pytz.FixedOffset(60)


@pytest.mark.unittest
def test_metadata_compare_valueerror():
    file_data_comp = os.path.join(dir_testdata,'VLISSGN_ana.txt')
    comp_fromfile = hatyan.read_components(filename=file_data_comp)
    comp_fromfile_fake = comp_fromfile.copy()
    comp_fromfile_fake.attrs["nodalfactors"] = False
    
    meta1 = metadata_from_obj(comp_fromfile)
    meta2 = metadata_from_obj(comp_fromfile_fake)
    with pytest.raises(ValueError) as e:
        metadata_compare([meta1,meta2])
    assert "equal" in str(e.value)


@pytest.mark.unittest
def test_anapred_metadata():
    current_station = 'VLISSGN'
    file_ts = os.path.join(dir_testdata, f'{current_station}_obs1.txt')
    ts_measurements_group0 = hatyan.read_dia(filename=file_ts, station=current_station)
    
    comp = hatyan.analysis(ts_measurements_group0,const_list='month')
    
    pred_xfac0 = hatyan.prediction(comp, times=ts_measurements_group0.index)
    # we also test if metadata is correctly passed if e.g. xfac is not in line with xfac in components file
    comp.attrs["xfac"] = True
    pred_xfac1 = hatyan.prediction(comp, times=ts_measurements_group0.index)
    
    meta_fromts_xfac0 = metadata_from_obj(pred_xfac0)
    meta_fromts_xfac1 = metadata_from_obj(pred_xfac1)
    
    meta_expected_xfac0 = {'station': 'VLISSGN',
     'grootheid': 'WATHTBRKD',
     'eenheid': 'm',
     'vertref': 'NAP',
     'TYP': 'TE',
     'groepering': 'NVT',
     'origin': 'from timeseries dia file',
     'nodalfactors': True,
     'xfac': False,
     'fu_alltimes': True,
     'source': 'schureman'}
    
    meta_expected_xfac1 = meta_expected_xfac0.copy()
    meta_expected_xfac1["xfac"] = True
    
    # compare metadata (raises ValueError if not equal)
    metadata_compare([meta_fromts_xfac0, meta_expected_xfac0])
    metadata_compare([meta_fromts_xfac1, meta_expected_xfac1])

    assert pred_xfac0.index.tz == pytz.FixedOffset(60)
    assert pred_xfac1.index.tz == pytz.FixedOffset(60)


@pytest.mark.unittest
def test_hwlw_metadata():
    current_station = 'VLISSGN'
    file_ts = os.path.join(dir_testdata, f'{current_station}_obs1.txt')
    ts_measurements_group0 = hatyan.read_dia(filename=file_ts, station=current_station)
    comp = hatyan.analysis(ts_measurements_group0,const_list='month')
    pred = hatyan.prediction(comp, times=ts_measurements_group0.index)
    
    meas_ext = hatyan.calc_HWLW(ts_measurements_group0)
    meas_ext = hatyan.calc_HWLWnumbering(meas_ext)
    pred_ext = hatyan.calc_HWLW(pred)
    pred_ext = hatyan.calc_HWLWnumbering(pred_ext)
    
    meas_ext_meta = metadata_from_obj(meas_ext)
    pred_ext_meta = metadata_from_obj(pred_ext)
    
    meas_ext_meta_expected = {'station': 'VLISSGN',
     'grootheid': 'WATHTE',
     'eenheid': 'm',
     'vertref': 'NAP',
     'TYP': 'TE',
     'groepering': 'NVT',
     'origin': 'from timeseries dia file'}
    
    pred_ext_meta_expected = {'station': 'VLISSGN',
     'grootheid': 'WATHTBRKD',
     'eenheid': 'm',
     'vertref': 'NAP',
     'TYP': 'TE',
     'groepering': 'NVT',
     'origin': 'from timeseries dia file',
     'nodalfactors': True,
     'xfac': False,
     'fu_alltimes': True,
     'source': 'schureman'}
    
    # compare metadata (raises ValueError if not equal)
    metadata_compare([pred_ext_meta, pred_ext_meta_expected])
    metadata_compare([meas_ext_meta, meas_ext_meta_expected])
    
    assert ts_measurements_group0.index.tz == pytz.FixedOffset(60)
    assert meas_ext.index.tz == pytz.FixedOffset(60)
    assert pred.index.tz == pytz.FixedOffset(60)
    assert pred_ext.index.tz == pytz.FixedOffset(60)


@pytest.mark.unittest
def test_readts_dia_metadata_multiblock():
    file_ts = os.path.join(dir_testdata, 'hoek_har.dia')
    ts_measurements_group0 = hatyan.read_dia(filename=file_ts, station='HOEKVHLD', block_ids='allstation')
    meta_fromts = metadata_from_obj(ts_measurements_group0)
    
    meta_expected = {'station': 'HOEKVHLD',
     'grootheid': 'WATHTE',
     'eenheid': 'm',
     'vertref': 'NAP',
     'TYP': 'TN',
     'groepering': 'GETETM2',
     'origin': 'from timeseries dia file'}
    
    assert meta_fromts == meta_expected
    assert ts_measurements_group0.index.tz == pytz.FixedOffset(60)


@pytest.mark.unittest
def test_metadata_compare():
    metadata = {
        'vertref': 'NAP',
        'station': 'VLISSGN',
        'TYP': 'TE',
        'groepering': 'NVT',
        'grootheid': 'WATHTE',
        'eenheid': 'm',
        }
    
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


@pytest.mark.unittest
def test_wns_from_metadata_invalid():
    metadata_invalid = {'grootheid':'WATHTBRKD', 'eenheid':'m', 'vertref':'MSL'}
    with pytest.raises(ValueError) as e:
        _ = wns_from_metadata(metadata_invalid)
    assert "combination of quantity/unit/vertref not defined" in str(e.value)


@pytest.mark.unittest
def test_metadata_persistence():
    current_station = 'VLISSGN'
    file_comp = os.path.join(dir_testdata, f'{current_station}_ana.txt')
    comp = hatyan.read_components(filename=file_comp)
    assert len(comp.attrs) > 0
    
    assert len(comp.iloc[:10].attrs) > 0
    assert len(comp["A"].attrs) > 0
    assert len(comp.max().attrs) > 0
