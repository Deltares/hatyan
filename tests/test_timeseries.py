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
from netCDF4 import Dataset, num2date
import hatyan
from hatyan.metadata import metadata_from_obj, metadata_compare

dir_tests = os.path.dirname(__file__) #F9 doesnt work, only F5 (F5 also only method to reload external definition scripts)
dir_testdata = os.path.join(dir_tests,'data_unitsystemtests')


@pytest.mark.unittest
def test_get_diaxycoords():
    file_pred = os.path.join(dir_testdata,"VLISSGN_pre.txt")
    stat_x, stat_y = hatyan.get_diaxycoords(filename=file_pred, crs=28992)
    assert np.isclose(stat_x, 30480.0)
    assert np.isclose(stat_y, 385220.0)
    
    file_pred = os.path.join(dir_testdata,"hoek_har.dia")
    stat_x, stat_y = hatyan.get_diaxycoords(filename=file_pred, crs=28992)
    assert np.allclose(stat_x, [67930., 67930., 67930., 49862.])
    assert np.allclose(stat_y, [444000., 444000., 444000., 431612.])


@pytest.mark.unittest
def test_readwrite_tsdia_rounding():
    file_pred = os.path.join(dir_testdata, "VLISSGN_pre.txt")
    file_new = 'temp_dia.txt'
    ts_pred = hatyan.read_dia(file_pred)
    for diff in [0,0.004, -0.004]:
        ts_pred_rounddiff = hatyan.read_dia(file_pred)
        ts_pred_rounddiff['values'] = ts_pred['values'] + diff
        hatyan.write_dia(ts=ts_pred_rounddiff, filename=file_new)
        ts_new = hatyan.read_dia(file_new)
        
        assert np.allclose(ts_pred, ts_new)
    os.remove(file_new)


@pytest.mark.unittest
def test_readwrite_tsdia_ext_rounding():
    file_pred = os.path.join(dir_testdata, "VLISSGN_ext.txt")
    file_new = 'temp_dia_ext.txt'
    ts_pred = hatyan.read_dia(file_pred)
    for diff in [0, 0.004, -0.004]:
        ts_pred_rounddiff = hatyan.read_dia(file_pred)
        ts_pred_rounddiff['values'] = ts_pred['values'] + diff
        hatyan.write_dia(ts=ts_pred_rounddiff, filename=file_new)
        ts_new = hatyan.read_dia(file_new)
        
        assert np.allclose(ts_pred, ts_new)
    os.remove(file_new)


@pytest.mark.unittest
def test_readts_dia_multifile():
    file_data_comp0 = [os.path.join(dir_testdata,f'VLISSGN_obs{file_id}.txt') for file_id in [1,2,3,4]]
    ts_measurements_group0 = hatyan.read_dia(filename=file_data_comp0, station='VLISSGN')
    
    assert len(ts_measurements_group0) == 35064
    assert np.isclose(ts_measurements_group0['values'].iloc[0], -1.24)
    assert np.isclose(ts_measurements_group0['values'].iloc[-1], -1.5)
    assert (ts_measurements_group0['qualitycode'] != 0).sum() == 82
    assert list(np.unique(ts_measurements_group0['qualitycode'])) == [ 0, 25]


@pytest.mark.unittest
def test_readts_dia_multiblock():
    
    file1 = os.path.join(dir_testdata,'hoek_har.dia')
    #ts_measurements_group0_extno = hatyan.read_dia(filename=file1, station='HOEKVHLD')
    ts_measurements_group0_ext0 = hatyan.read_dia(filename=file1, station='HOEKVHLD', block_ids=0)
    ts_measurements_group0_ext1 = hatyan.read_dia(filename=file1, station='HOEKVHLD', block_ids=1)
    ts_measurements_group0_ext2 = hatyan.read_dia(filename=file1, station='HOEKVHLD', block_ids=2)
    ts_measurements_group0_ext012 = hatyan.read_dia(filename=file1, station='HOEKVHLD', block_ids=[0,1,2])
    ts_measurements_group0_extall = hatyan.read_dia(filename=file1, station='HOEKVHLD', block_ids='allstation')
    
    assert len(ts_measurements_group0_ext0) == 3977
    assert len(ts_measurements_group0_ext1) == 9913
    assert len(ts_measurements_group0_ext2) == 9403
    assert len(ts_measurements_group0_ext012) == 23293
    assert len(ts_measurements_group0_extall) == 23293


@pytest.mark.unittest
def test_readwrite_noos():

    file_data_comp0 = os.path.join(dir_testdata,'VLISSGN_waterlevel_20180101_20180401.noos')
    
    #component groups
    filename_out = "noos_test.txt"
    ts_noos1 = hatyan.read_noos(filename=file_data_comp0)
    hatyan.write_noos(ts_noos1, filename=filename_out)
    ts_noos2 = hatyan.read_noos(filename=filename_out)
    
    assert (ts_noos1 == ts_noos2).all().all()
    
    # compare metadata, raises ValueError if not equal
    meta1 = metadata_from_obj(ts_noos1)
    meta2 = metadata_from_obj(ts_noos2)
    metadata_compare([meta1, meta2])
    os.remove(filename_out)


@pytest.mark.unittest
def test_readts_noos_resamplecrop():

    file_data_comp0 = os.path.join(dir_testdata,'VLISSGN_waterlevel_20180101_20180401.noos')
    times_ext_comp0 = slice(dt.datetime(2018,1,1),dt.datetime(2018,4,1))

    #component groups
    ts_measurements_group0 = hatyan.read_noos(filename=file_data_comp0)
    ts_measurements_group0_res = hatyan.resample_timeseries(ts_measurements_group0, timestep_min=10)
    ts_measurements_group0_rescrop = hatyan.crop_timeseries(ts_measurements_group0_res, times=times_ext_comp0)
    
    assert len(ts_measurements_group0) == 12752
    assert len(ts_measurements_group0_res) == 12961
    assert len(ts_measurements_group0_rescrop) == 12961
    assert ts_measurements_group0_rescrop.index[0] == pd.Timestamp('2018-01-01')
    assert ts_measurements_group0_rescrop.index[-1] == pd.Timestamp('2018-04-01')
    assert np.isclose(ts_measurements_group0_rescrop['values'].iloc[0], 2.5)
    assert np.isclose(ts_measurements_group0_rescrop['values'].iloc[-1], 1.05)


@pytest.mark.unittest
def test_readts_dia_equidistant_singlefile_hasfreq():
    """
    When reading asingle equidistant diafile,
    there should be a freq attribute that is not None
    """
    file_ts = os.path.join(dir_testdata,'VLISSGN_obs1.txt')
    
    ts_pd = hatyan.read_dia(filename=file_ts)
    
    # check ts length (all four files are added)
    assert len(ts_pd) == 8760
    
    # assert on freq attribute
    assert hasattr(ts_pd.index,'freq')
    assert isinstance(ts_pd.index.freq,pd.offsets.Minute)
    assert ts_pd.index.freq is not None
    assert np.isclose(ts_pd.index.freq.nanos/1e9, 3600)


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
    assert np.isclose(ts_pd.index.freq.nanos/1e9, 600)

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
    
    ts_pd = hatyan.read_dia(filename=file_ts)
    
    # check ts length (all four files are added)
    assert len(ts_pd) == 35064
    
    # check index dtype
    assert ts_pd.index.dtype == "datetime64[ns, pytz.FixedOffset(60)]"
    
    # checks for freq attribute
    assert hasattr(ts_pd.index,'freq')
    assert isinstance(ts_pd.index.freq,pd.offsets.Minute)
    assert ts_pd.index.freq is not None
    assert np.isclose(ts_pd.index.freq.nanos/1e9, 3600)


@pytest.mark.unittest
def test_readts_dia_equidistant_multifile_glob_hasfreq():
    """
    When providing a file pattern for reading multiple equidistant diafiles,
    glob is supported, this test checks if it results in the correct data
    """
    file_ts = os.path.join(dir_testdata,'VLISSGN_obs?.txt')
    
    ts_pd = hatyan.read_dia(filename=file_ts)
    
    assert len(ts_pd) == 35064

    # check index dtype
    assert ts_pd.index.dtype == "datetime64[ns, pytz.FixedOffset(60)]"
    
    # checks for freq attribute
    assert hasattr(ts_pd.index,'freq')
    assert isinstance(ts_pd.index.freq,pd.offsets.Minute)
    assert ts_pd.index.freq is not None
    assert np.isclose(ts_pd.index.freq.nanos/1e9, 3600)


@pytest.mark.unittest
def test_readwrite_diawia():
    
    current_station = 'HOEKVHLD'
    
    file_dia_wl = os.path.join(dir_testdata,'diawia_%s_astro_tijdreeks.dia'%(current_station))
    file_dia_ext = os.path.join(dir_testdata,'diawia_%s_astro_extremen.dia'%(current_station))
    
    for file_dia in [file_dia_wl,file_dia_ext]:
        file_wia = file_dia.replace('.dia','.wia')
        file_dia_out = file_dia.replace('.dia','_out.dia')
        file_wia_out = file_dia.replace('.dia','_out.wia')
        ts_dia = hatyan.read_dia(filename=file_dia, station=current_station)
        ts_wia = hatyan.read_dia(filename=file_wia, station=current_station)
        assert (ts_dia==ts_wia).all().all() #check if wia and dia input is equal
        
        #write to files
        hatyan.write_dia(ts=ts_dia, filename=file_dia_out)
        hatyan.write_dia(ts=ts_wia, filename=file_wia_out, headerformat='wia')
        
        #read from new files
        ts_dia_new = hatyan.read_dia(filename=file_dia_out, station=current_station)
        ts_wia_new = hatyan.read_dia(filename=file_wia_out, station=current_station)
        assert np.allclose(ts_dia, ts_dia_new) #check if wia and dia input is equal
        assert np.allclose(ts_wia, ts_wia_new) #check if wia and dia input is equal
        
        #meta
        meta_dia = metadata_from_obj(ts_dia)
        meta_wia = metadata_from_obj(ts_wia)
        meta_dia_new = metadata_from_obj(ts_dia_new)
        meta_wia_new = metadata_from_obj(ts_wia_new)
        assert (meta_dia == meta_wia == meta_dia_new == meta_wia_new)
        
        #remove files
        os.remove(file_dia_out)
        os.remove(file_wia_out)


@pytest.mark.unittest
def test_crop_timeseries():
    
    current_station = 'VLISSGN'
    times_ext = slice("2019-01-01 00:00:00 +00:00", "2019-06-01 00:00:00 +00:00")
    
    file_pred = os.path.join(dir_testdata,f'{current_station}_pre.txt')
    ts_prediction = hatyan.read_dia(filename=file_pred, station=current_station)
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
    ts_prediction = hatyan.read_dia(filename=file_pred, station=current_station)
    ts_prediction_res = hatyan.resample_timeseries(ts_prediction, timestep_min=timestep_min)
    
    assert len(ts_prediction_res) == 4380
    assert ts_prediction_res.index[0] == pd.Timestamp(ts_prediction.index[0])
    assert ts_prediction_res.index[-1] == pd.Timestamp("2019-12-31 22:00 +01:00")
    
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
    
    times_pred = slice(dt.datetime(2019,1,1),dt.datetime(2020,1,1), "10min")
    
    comp_merged = hatyan.read_components(filename=file_data_comp0)
    
    #prediction and validation
    ts_prediction = hatyan.prediction(comp=comp_merged, times=times_pred)
    
    #write to file
    fname_pred = 'prediction_%s_%s.dia'%(times_pred.step,current_station)
    hatyan.write_dia(ts=ts_prediction, filename=fname_pred)
    
    #read from file
    ts_prediction_fromfile = hatyan.read_dia(filename=fname_pred, station=current_station)
    os.remove(fname_pred)
    
    # assert max differences
    ts_diff = ts_prediction_fromfile - ts_prediction
    assert (ts_diff['values'].abs()<=0.005).all()


@pytest.mark.unittest
def test_timeseries_fft():
    file_pred = os.path.join(dir_testdata, "VLISSGN_pre.txt")
    ts_pred = hatyan.read_dia(file_pred)
        
    hatyan_freqs_suggestions_schureman = hatyan.timeseries_fft(ts_pred, min_prominence=2000, plot_fft=True, source="schureman")
    hatyan_freqs_suggestions_foreman = hatyan.timeseries_fft(ts_pred, min_prominence=2000, plot_fft=True, source="foreman")
    const_list_schureman = hatyan_freqs_suggestions_schureman.index.unique().tolist()
    const_list_foreman = hatyan_freqs_suggestions_foreman.index.unique().tolist()
    
    assert const_list_schureman == ['O1', 'N2', 'M2', 'L2B', 'S2', 'K2', 'M4', 'M6']
    assert const_list_foreman ==   ['O1', 'N2', 'M2', 'L2', 'S2', 'K2', 'M4', 'M6']


@pytest.mark.unittest
def test_plot_timeseries():
    file_pred = os.path.join(dir_testdata, "VLISSGN_pre.txt")
    file_ext = os.path.join(dir_testdata, "VLISSGN_ext.txt")
    ts_pred = hatyan.read_dia(file_pred)
    ts_ext = hatyan.read_dia(file_ext)
    hatyan.plot_timeseries(ts=ts_pred, ts_validation=ts_pred, ts_ext=ts_ext, ts_ext_validation=ts_ext)


@pytest.mark.unittest
def test_plot_HWLW_validatestats():
    file_ext = os.path.join(dir_testdata, "VLISSGN_ext.txt")
    ts_ext = hatyan.read_dia(file_ext)
    ts_ext_nos = hatyan.calc_HWLWnumbering(ts_ext)

    hatyan.plot_HWLW_validatestats(ts_ext=ts_ext, ts_ext_validation=ts_ext)
    hatyan.plot_HWLW_validatestats(ts_ext=ts_ext_nos, ts_ext_validation=ts_ext_nos)
    
    # check if original dataframe was not altered
    assert "HWLWno" not in ts_ext.columns
    assert "HWLWno" in ts_ext_nos.columns


@pytest.mark.unittest
def test_nyquist_folding():
    const_list = hatyan.get_const_list_hatyan('year')
    
    # drop_list = ['S4','3M2S10','2SM6','4M2S12'] # overlappende frequenties na folding, en S4 is precies op Nyquist frequentie
    # for const in drop_list:
    #     const_list.remove(const)
    
    station = 'DENHDR'
    year = 1955
    url_dataraw = 'https://watersysteemdata.deltares.nl/thredds/fileServer/watersysteemdata/Wadden/ddl/raw/waterhoogte/'
    file_csv = url_dataraw+f'{station}_OW_WATHTE_NVT_NAP_{year}_ddl_wq.csv'
    
    data_pd = pd.read_csv(file_csv,sep=';',parse_dates=['tijdstip'])
    
    ts_meas_raw = pd.DataFrame({'values':data_pd['numeriekewaarde'].values/100},index=data_pd['tijdstip'].dt.tz_localize(None)) #read dataset and convert to DataFrame with correct format #TODO: tijdzone is MET (UTC+1), ook van de resulterende getijcomponenten. Is dat wenselijk?
    ts_meas_raw = ts_meas_raw.sort_index(axis='index') #sort dataset on times (not necessary but easy for debugging)
    
    # check if we get "Exception: there is a component on the Nyquist frequency (0.16666666666666666 [1/hr]), this not possible"
    # without doing the nyquist check, we get a MatrixConditionTooHigh error instead
    with pytest.raises(ValueError) as e:
        hatyan.analysis(ts=ts_meas_raw, const_list=const_list, nodalfactors=True, xfac=True, fu_alltimes=True)
    assert "nyquist" in str(e.value).lower()


@pytest.mark.unittest
def test_writenetcdf():
    
    current_station = 'VLISSGN'
    
    file_pred = os.path.join(dir_testdata,f'{current_station}_pre.txt')
    ts_prediction = hatyan.read_dia(filename=file_pred, station=current_station)
    ts_ext_prediction = hatyan.calc_HWLW(ts=ts_prediction)
    
    file_nc = 'prediction_10m_%s.nc'%(current_station)
    hatyan.write_netcdf(ts=ts_prediction, ts_ext=ts_ext_prediction, filename=file_nc)
    
    data_nc = Dataset(file_nc,'r')
    
    timevar = data_nc.variables['time']
    timevar_dt = num2date(timevar[:],units=timevar.units, only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    
    #put netcdf file contents in pandas DataFrame for usage in hatyan
    ts_pd = pd.DataFrame({'values':data_nc.variables['waterlevel_astro'][:,0]}, index=timevar_dt)
    print(ts_pd)
    
    assert list(data_nc.dimensions.keys()) == ['stations', 'statname_len', 'time', 'analysis_time', 'time_HW', 'time_LW']
    assert list(data_nc.variables.keys()) == ['stations', 'analysis_time', 'time', 'waterlevel_astro', 'time_HW', 'waterlevel_astro_HW', 'time_LW', 'waterlevel_astro_LW']
    assert timevar_dt[0] == ts_prediction.index[0].tz_convert(None)
    assert timevar_dt[-1] == ts_prediction.index[-1].tz_convert(None)
    assert 'title' in data_nc.__dict__.keys()
    assert 'institution' in data_nc.__dict__.keys()
    assert 'source' in data_nc.__dict__.keys()
    assert timevar.units == 'minutes since 1900-01-01 00:00:00 +0100'

    data_nc.close()
    os.remove(file_nc)


@pytest.mark.unittest
def test_writenetcdf_nosidx():
    
    current_station = 'VLISSGN'
    
    file_pred = os.path.join(dir_testdata,f'{current_station}_pre.txt')
    ts_prediction = hatyan.read_dia(filename=file_pred, station=current_station)
    ts_ext_prediction = hatyan.calc_HWLW(ts=ts_prediction)
    ts_ext_prediction = hatyan.calc_HWLWnumbering(ts_ext=ts_ext_prediction)
    
    file_nc = 'prediction_10m_%s.nc'%(current_station)
    hatyan.write_netcdf(ts=ts_prediction, ts_ext=ts_ext_prediction, filename=file_nc, nosidx=True)
    
    data_nc = Dataset(file_nc,'r')
    
    timevar = data_nc.variables['time']
    timevar_dt = num2date(timevar[:],units=timevar.units, only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    
    #put netcdf file contents in pandas DataFrame for usage in hatyan
    ts_pd = pd.DataFrame({'values':data_nc.variables['waterlevel_astro'][:,0]}, index=timevar_dt)
    print(ts_pd)
    
    assert list(data_nc.dimensions.keys()) == ['stations', 'statname_len', 'time', 'analysis_time', 'HWLWno']
    assert list(data_nc.variables.keys()) == ['stations', 'analysis_time', 'time', 'waterlevel_astro', 
                                              'HWLWno', 'times_astro_HW', 'waterlevel_astro_HW', 'times_astro_LW', 'waterlevel_astro_LW']
    assert timevar_dt[0] == ts_prediction.index[0].tz_convert(None)
    assert timevar_dt[-1] == ts_prediction.index[-1].tz_convert(None)
    assert 'title' in data_nc.__dict__.keys()
    assert 'institution' in data_nc.__dict__.keys()
    assert 'source' in data_nc.__dict__.keys()
    assert timevar.units == 'minutes since 1900-01-01 00:00:00 +0100'

    data_nc.close()
    os.remove(file_nc)

