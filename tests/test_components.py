# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 14:33:40 2023

@author: veenstra
"""

import os
import pytest
import pandas as pd
import pytz
import datetime as dt
import numpy as np
import hatyan
from hatyan.metadata import metadata_from_obj
from hatyan.components import _get_tzone_minutes

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

    comp_orig = hatyan.read_components(filename=file_orig)
    comp_orig.attrs['tzone'] = None
    with pytest.raises(ValueError) as e:
        hatyan.write_components(comp_orig, filename=file_new)
    assert "tzone=None" in str(e.value)


@pytest.mark.unittest
def test_writecomponents_fromanalysis(tmp_path):
    """
    this is of added value to check if all required metadata is present from an analysis
    """
    current_station = 'VLISSGN'
    file_data_comp0 = os.path.join(dir_testdata,f'{current_station}_obs1.txt')
    ts = hatyan.read_dia(filename=file_data_comp0, station=current_station)
    
    comp = hatyan.analysis(ts=ts, const_list='month', fu_alltimes=False)
    file_comp = os.path.join(tmp_path, "temp_comp.txt")
    hatyan.write_components(comp, file_comp)


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


@pytest.mark.unittest
def test_merge_componentgroups():
    file_comp1 = os.path.join(dir_testdata,'VLISSGN_ana.txt')
    comp_fromfile1 = hatyan.read_components(filename=file_comp1)
    file_comp2 = os.path.join(dir_testdata,'HOEKVHLD_ana.txt')
    comp_fromfile2 = hatyan.read_components(filename=file_comp2)
    comp_fromfile2.attrs["station"] = "VLISSGN"
    
    comp_merged = hatyan.merge_componentgroups(comp_main=comp_fromfile1, comp_sec=comp_fromfile2.loc[['SA','SM']])
    
    assert (comp_fromfile1!=comp_merged).sum().sum() == 4


@pytest.mark.unittest
def test_merge_componentgroups_comparesettings():
    file_data_comp = os.path.join(dir_testdata,'VLISSGN_ana.txt')
    comp_fromfile = hatyan.read_components(filename=file_data_comp)
    
    with pytest.raises(ValueError):
        comp_fromfile_fake = comp_fromfile.copy()
        comp_fromfile_fake.attrs["nodalfactors"] = False
        _ = hatyan.merge_componentgroups(comp_main=comp_fromfile, comp_sec=comp_fromfile_fake.loc[['SA','SM']])

    with pytest.raises(ValueError):
        comp_fromfile_fake = comp_fromfile.copy()
        comp_fromfile_fake.attrs["xfac"] = False
        _ = hatyan.merge_componentgroups(comp_main=comp_fromfile, comp_sec=comp_fromfile_fake.loc[['SA','SM']])

    with pytest.raises(ValueError):
        comp_fromfile_fake = comp_fromfile.copy()
        comp_fromfile_fake.attrs["fu_alltimes"] = True
        _ = hatyan.merge_componentgroups(comp_main=comp_fromfile, comp_sec=comp_fromfile_fake.loc[['SA','SM']])

    with pytest.raises(ValueError):
        comp_fromfile_fake = comp_fromfile.copy()
        comp_fromfile_fake.attrs["source"] = 'foreman'
        _ = hatyan.merge_componentgroups(comp_main=comp_fromfile, comp_sec=comp_fromfile_fake.loc[['SA','SM']])


@pytest.mark.unittest
def test_get_tzone_minutes():
    # relevant for ddlpy timeseries
    tstart = pd.Timestamp("2020-01-01 00:00:00 +01:00")
    tzone_min = _get_tzone_minutes(tstart.tz)
    assert isinstance(tstart.tz, dt.timezone)
    assert tzone_min == 60
    
    # hatyan dia timeseries
    tstart = pd.Timestamp("2020-01-01 00:00:00", tz=pytz.FixedOffset(60))
    tzone_min = _get_tzone_minutes(tstart.tz)
    assert tzone_min == 60
    
    tstart = pd.Timestamp("2020-01-01 00:00:00")
    with pytest.raises(NotImplementedError) as e:
        _get_tzone_minutes(tstart.tz)
    assert "tzone of type <class 'NoneType'> is not yet supported" in str(e.value)
