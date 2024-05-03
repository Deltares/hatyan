# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 14:33:40 2023

@author: veenstra
"""

import os
import pytest
import numpy as np
import hatyan
from hatyan.metadata import metadata_from_obj

dir_tests = os.path.dirname(__file__) #F9 doesnt work, only F5 (F5 also only method to reload external definition scripts)
dir_testdata = os.path.join(dir_tests,'data_unitsystemtests')


@pytest.mark.unittest
def test_read_write_components():
    
    current_station = 'HOEKVHLD'
    file_orig = os.path.join(dir_testdata,f'{current_station}_ana.txt')
    file_new = 'temp_components.txt'
    
    comp_orig = hatyan.read_components(filename=file_orig)
    hatyan.write_components(comp_orig, filename=file_new)
    comp_new = hatyan.read_components(filename=file_new)
    
    meta_orig = metadata_from_obj(comp_orig)
    meta_new = metadata_from_obj(comp_new)
    
    assert np.allclose(comp_orig, comp_new)
    assert meta_orig == meta_new
    os.remove(file_new)


@pytest.mark.unittest
def test_read_write_components_nometadata():
    """
    Test for component files written with hatyan 2.7.0 or lower,
    these files lack essential metadata for STAT, PERD and CODE lines.
    This tests tests whether the correct exception is raised.
    """
    
    file_orig = os.path.join(dir_testdata,'DENHDR_ana.txt')
    file_nometa = "temp_comp_dummymetadata.txt"
    with open(file_orig, "r") as f:
        data = f.readlines()
    with open(file_nometa, "w") as f:
        for line in data:
            if line.startswith("STAT"):
                continue
            if line.startswith("PERD"):
                continue
            if line.startswith("CODE"):
                continue
            f.write(line)
    
    try:
        _ = hatyan.read_components(file_nometa)
    except KeyError as e:
        assert "No STAT/PERD metadata available in component file" in str(e)

    os.remove(file_nometa)


@pytest.mark.unittest
def test_read_write_components_nondefaultsettings():
    
    current_station = 'HOEKVHLD'
    file_orig = os.path.join(dir_testdata,f'{current_station}_ana.txt')
    file_new = 'temp_components.txt'
    
    comp_orig = hatyan.read_components(filename=file_orig)
    comp_orig.attrs["nodalfactors"] = False
    with pytest.raises(ValueError) as e:
        hatyan.write_components(comp_orig, filename=file_new)
    assert "nodalfactors" in str(e.value)
    
    comp_orig = hatyan.read_components(filename=file_orig)
    comp_orig.attrs["fu_alltimes"] = True
    with pytest.raises(ValueError) as e:
        hatyan.write_components(comp_orig, filename=file_new)
    assert "fu_alltimes" in str(e.value)
    
    comp_orig = hatyan.read_components(filename=file_orig)
    comp_orig.attrs["source"] = "notschureman"
    with pytest.raises(ValueError) as e:
        hatyan.write_components(comp_orig, filename=file_new)
    assert "source" in str(e.value)


@pytest.mark.unittest
def test_plot_components_validation():
    
    current_station = 'HOEKVHLD'
    file_data_comp0 = os.path.join(dir_testdata,f'{current_station}_ana.txt')
    
    comp_merged = hatyan.read_components(filename=file_data_comp0)
    _ = hatyan.plot_components(comp=comp_merged, comp_validation=comp_merged)


@pytest.mark.unittest
def test_plot_components_allyears():
    current_station = 'VLISSGN'
    file_data_comp0 = os.path.join(dir_testdata,f'{current_station}_obs?.txt')
    ts_measurements = hatyan.read_dia(filename=file_data_comp0, station=current_station)
    
    ts_comp, ts_comp_all = hatyan.analysis(ts=ts_measurements, const_list='month', analysis_perperiod="Y", return_allperiods=True)
    _ = hatyan.plot_components(comp=ts_comp, comp_allperiods=ts_comp_all)


@pytest.mark.systemtest
def test_components_timeshift():
    
    timeshift_hr = 1
    current_station = 'HOEKVHLD'
    file_data_comp0 = os.path.join(dir_testdata,f'{current_station}_ana.txt')
    
    comp_merged = hatyan.read_components(filename=file_data_comp0)
    
    comp_shift = hatyan.components_timeshift(comp_merged,hours=timeshift_hr)
    
    # check timeshift shift
    assert np.abs(comp_shift.loc['SA','phi_deg'] - 221.54106863959504) < 1e-9
    
    # check metadata contents
    comp_merged_meta = metadata_from_obj(comp_merged)
    comp_shift_meta = metadata_from_obj(comp_shift)
    
    comp_merged_tz_min = comp_merged_meta['tzone']._minutes
    comp_shift_tz_min = comp_shift_meta['tzone']._minutes
    
    assert comp_merged_tz_min + timeshift_hr*60 == comp_shift_tz_min


@pytest.mark.unittest
def test_merge_componentgroups_comparesettings():
    file_data_comp = os.path.join(dir_testdata,'VLISSGN_ana.txt')
    comp_fromfile = hatyan.read_components(filename=file_data_comp)
    
    with pytest.raises(ValueError):
        comp_fromfile_fake = comp_fromfile.copy()
        comp_fromfile_fake.attrs["nodalfactors"] = False
        _ = hatyan.merge_componentgroups(comp_main=comp_fromfile, comp_sec=comp_fromfile_fake, comp_sec_list=['SA','SM'])

    with pytest.raises(ValueError):
        comp_fromfile_fake = comp_fromfile.copy()
        comp_fromfile_fake.attrs["xfac"] = False
        _ = hatyan.merge_componentgroups(comp_main=comp_fromfile, comp_sec=comp_fromfile_fake, comp_sec_list=['SA','SM'])

    with pytest.raises(ValueError):
        comp_fromfile_fake = comp_fromfile.copy()
        comp_fromfile_fake.attrs["fu_alltimes"] = True
        _ = hatyan.merge_componentgroups(comp_main=comp_fromfile, comp_sec=comp_fromfile_fake, comp_sec_list=['SA','SM'])

    with pytest.raises(ValueError):
        comp_fromfile_fake = comp_fromfile.copy()
        comp_fromfile_fake.attrs["source"] = 'foreman'
        _ = hatyan.merge_componentgroups(comp_main=comp_fromfile, comp_sec=comp_fromfile_fake, comp_sec_list=['SA','SM'])
