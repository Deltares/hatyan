# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 13:13:31 2023

@author: veenstra
"""

import os
import pytest
import datetime as dt
import hatyan

dir_tests = os.path.dirname(__file__) #F9 doesnt work, only F5 (F5 also only method to reload external definition scripts)
dir_testdata = os.path.join(dir_tests,'data_unitsystemtests')


@pytest.mark.unittest
def test_write_tsdia_rounding():
    """
    rounding error occurred in older versions of hatyan2
    Therefore we test whether writing and reading a timeseries results in the same data (accuracy of 1cm)
    """
    
    current_station = 'HOEKVHLD'
    file_data_comp0 = os.path.join(dir_testdata,f'{current_station}_ana.txt')
    
    times_ext_pred = [dt.datetime(2019,1,1),dt.datetime(2020,1,1)]
    times_step_pred = 10
    
    comp_merged = hatyan.read_components(filename=file_data_comp0)
    
    #prediction and validation
    ts_prediction = hatyan.prediction(comp=comp_merged, nodalfactors=True, xfac=True, fu_alltimes=False, 
                                      times_ext=times_ext_pred, timestep_min=times_step_pred)
    
    #write to file
    fname_pred = 'prediction_%im_%s.dia'%(times_step_pred,current_station)
    hatyan.write_tsdia(ts=ts_prediction, filename=fname_pred)
    
    #read from file
    ts_prediction_fromfile = hatyan.readts_dia(filename=fname_pred, station=current_station)
    
    # assert max differences
    ts_diff = ts_prediction_fromfile -ts_prediction
    assert (ts_diff['values'].abs()<=0.005).all()
