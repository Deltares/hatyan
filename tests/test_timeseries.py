# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 12:14:02 2023

@author: veenstra
"""

import os
import pytest
import pandas as pd
import datetime as dt
import numpy as np
import hatyan
from hatyan.metadata import metadata_from_obj, metadata_compare

dir_tests = os.path.dirname(__file__) #F9 doesnt work, only F5 (F5 also only method to reload external definition scripts)
dir_testdata = os.path.join(dir_tests,'data_unitsystemtests')


@pytest.mark.unittest
def test_readts_dia_multifile():
    file_data_comp0 = [os.path.join(dir_testdata,f'VLISSGN_obs{file_id}.txt') for file_id in [1,2,3,4]]
    ts_measurements_group0 = hatyan.readts_dia(filename=file_data_comp0, station='VLISSGN')
    
    assert len(ts_measurements_group0) == 35064
    assert np.isclose(ts_measurements_group0['values'].iloc[0], -1.24)
    assert np.isclose(ts_measurements_group0['values'].iloc[-1], -1.5)
    assert (ts_measurements_group0['qualitycode'] != 0).sum() == 82
    assert list(np.unique(ts_measurements_group0['qualitycode'])) == [ 0, 25]


@pytest.mark.unittest
def test_readts_dia_multiblock():
    
    file1 = os.path.join(dir_testdata,'hoek_har.dia')
    #ts_measurements_group0_extno = hatyan.readts_dia(filename=file1, station='HOEKVHLD')
    ts_measurements_group0_ext0 = hatyan.readts_dia(filename=file1, station='HOEKVHLD', block_ids=0)
    ts_measurements_group0_ext1 = hatyan.readts_dia(filename=file1, station='HOEKVHLD', block_ids=1)
    ts_measurements_group0_ext2 = hatyan.readts_dia(filename=file1, station='HOEKVHLD', block_ids=2)
    ts_measurements_group0_ext012 = hatyan.readts_dia(filename=file1, station='HOEKVHLD', block_ids=[0,1,2])
    ts_measurements_group0_extall = hatyan.readts_dia(filename=file1, station='HOEKVHLD', block_ids='allstation')
    
    assert len(ts_measurements_group0_ext0) == 3977
    assert len(ts_measurements_group0_ext1) == 9913
    assert len(ts_measurements_group0_ext2) == 9403
    assert len(ts_measurements_group0_ext012) == 23293
    assert len(ts_measurements_group0_extall) == 23293


@pytest.mark.unittest
def test_readts_noos_resamplecrop():

    file_data_comp0 = os.path.join(dir_testdata,'VLISSGN_waterlevel_20180101_20180401.noos')
    times_ext_comp0 = slice(dt.datetime(2018,1,1),dt.datetime(2018,4,1))

    #component groups
    ts_measurements_group0 = hatyan.readts_noos(filename=file_data_comp0)
    ts_measurements_group0_res = hatyan.resample_timeseries(ts_measurements_group0, timestep_min=10)
    ts_measurements_group0_rescrop = hatyan.crop_timeseries(ts_measurements_group0_res, times=times_ext_comp0)
    
    assert len(ts_measurements_group0) == 12752
    assert len(ts_measurements_group0_res) == 12961
    assert len(ts_measurements_group0_rescrop) == 12961
    assert ts_measurements_group0_rescrop.index[0] == pd.Timestamp('2018-01-01')
    assert ts_measurements_group0_rescrop.index[-1] == pd.Timestamp('2018-04-01')
    assert np.isclose(ts_measurements_group0_rescrop['values'][0], 2.5)
    assert np.isclose(ts_measurements_group0_rescrop['values'][-1], 1.05)


@pytest.mark.unittest
def test_readts_dia_equidistant_singlefile_hasfreq():
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
def test_pandas_concat_hasfreq():
    """
    freq is None in case of index with non-constant freq (eg 2020+2022 or if one timestep is skipped/duplicated)
    freq is pd.offsets.Minute if index is equidistant
    """

    def df_index(year):
        pd_year = pd.date_range(f'{year}-01-01',f'{year}-12-31 23:50', freq='10min')
        df_year = pd.DataFrame(index=pd_year)
        return df_year
    
    df_2020 = df_index(2020)
    df_2021 = df_index(2021)
    df_2022 = df_index(2022)
    
    ts_pd = pd.concat([df_2020,df_2021])
    ts_pd_nonequi = pd.concat([df_2020,df_2022])
    
    # assert on freq attribute
    assert hasattr(ts_pd.index,'freq')
    assert isinstance(ts_pd.index.freq,pd.offsets.Minute)
    assert ts_pd.index.freq is not None
    assert ts_pd.index.freq.nanos/1e9 == 600    

    # assert on freq attribute
    assert hasattr(ts_pd_nonequi.index,'freq')
    assert isinstance(ts_pd_nonequi.index.freq,type(None))
    assert ts_pd_nonequi.index.freq is None


@pytest.mark.unittest
def test_readts_dia_equidistant_multifile_hasfreq():
    """
    When reading multiple equidistant diafiles that combine into a continuous timeseries,
    there should be a freq attribute that is not None
    """
    file_ts = [os.path.join(dir_testdata,f'VLISSGN_obs{i}.txt') for i in [1,2,3,4]]
    
    ts_pd = hatyan.readts_dia(filename=file_ts)
    
    # check ts length (all four files are added)
    assert len(ts_pd) == 35064
    
    # check index dtype
    assert ts_pd.index.dtype == '<M8[ns]'
    
    # checks for freq attribute
    assert hasattr(ts_pd.index,'freq')
    assert isinstance(ts_pd.index.freq,pd.offsets.Minute)
    assert ts_pd.index.freq is not None
    assert ts_pd.index.freq.nanos/1e9 == 3600


@pytest.mark.unittest
def test_readts_dia_equidistant_multifile_glob_hasfreq():
    """
    When providing a file pattern for reading multiple equidistant diafiles,
    glob is supported, this test checks if it results in the correct data
    """
    file_ts = os.path.join(dir_testdata,'VLISSGN_obs?.txt')
    
    ts_pd = hatyan.readts_dia(filename=file_ts)
    
    assert len(ts_pd) == 35064

    # check index dtype
    assert ts_pd.index.dtype == '<M8[ns]'
    
    # checks for freq attribute
    assert hasattr(ts_pd.index,'freq')
    assert isinstance(ts_pd.index.freq,pd.offsets.Minute)
    assert ts_pd.index.freq is not None
    assert ts_pd.index.freq.nanos/1e9 == 3600


@pytest.mark.unittest
def test_crop_timeseries():
    
    current_station = 'VLISSGN'
    times_ext = slice(dt.datetime(2019,1,1),dt.datetime(2019,6,1))
    
    file_pred = os.path.join(dir_testdata,f'{current_station}_pre.txt')
    ts_prediction = hatyan.readts_dia(filename=file_pred, station=current_station)
    ts_prediction_cropped = hatyan.crop_timeseries(ts_prediction, times=times_ext)
    
    assert len(ts_prediction_cropped) == 21745
    assert ts_prediction_cropped.index[0] == pd.Timestamp(times_ext.start)
    assert ts_prediction_cropped.index[-1] == pd.Timestamp(times_ext.stop)
    
    pred_meta = metadata_from_obj(ts_prediction)
    pred_cropped_meta = metadata_from_obj(ts_prediction_cropped)
    metadata_compare([pred_meta,pred_cropped_meta])
    
    assert pred_cropped_meta['tstart'] == pd.Timestamp(times_ext.start)
    assert pred_cropped_meta['tstop'] == pd.Timestamp(times_ext.stop)



@pytest.mark.unittest
def test_resample_timeseries():
    
    current_station = 'VLISSGN'
    timestep_min = 120
    
    file_pred = os.path.join(dir_testdata,f'{current_station}_pre.txt')
    ts_prediction = hatyan.readts_dia(filename=file_pred, station=current_station)
    ts_prediction_res = hatyan.resample_timeseries(ts_prediction, timestep_min=timestep_min)
    
    assert len(ts_prediction_res) == 4380
    assert ts_prediction_res.index[0] == pd.Timestamp(ts_prediction.index[0])
    assert ts_prediction_res.index[-1] == pd.Timestamp("2019-12-31 22:00")
    
    pred_meta = metadata_from_obj(ts_prediction)
    pred_meta.pop('timestep_min')
    pred_res_meta = metadata_from_obj(ts_prediction_res)
    pred_res_meta.pop('timestep_min')
    metadata_compare([pred_meta,pred_res_meta])
    
    pred_res_meta = metadata_from_obj(ts_prediction_res)
    assert pred_res_meta['timestep_min'] == timestep_min


@pytest.mark.unittest
def test_write_tsdia_rounding():
    """
    rounding error occurred in older versions of hatyan2
    Therefore we test whether writing and reading a timeseries results in the same data (accuracy of 1cm)
    """
    
    current_station = 'HOEKVHLD'
    file_data_comp0 = os.path.join(dir_testdata,f'{current_station}_ana.txt')
    
    times_pred = slice(dt.datetime(2019,1,1),dt.datetime(2020,1,1), 10)
    
    comp_merged = hatyan.read_components(filename=file_data_comp0)
    
    #prediction and validation
    ts_prediction = hatyan.prediction(comp=comp_merged, nodalfactors=True, xfac=True, fu_alltimes=False, 
                                      times=times_pred)
    
    #write to file
    fname_pred = 'prediction_%im_%s.dia'%(times_pred.step,current_station)
    hatyan.write_tsdia(ts=ts_prediction, filename=fname_pred)
    
    #read from file
    ts_prediction_fromfile = hatyan.readts_dia(filename=fname_pred, station=current_station)
    os.remove(fname_pred)
    
    # assert max differences
    ts_diff = ts_prediction_fromfile - ts_prediction
    assert (ts_diff['values'].abs()<=0.005).all()
