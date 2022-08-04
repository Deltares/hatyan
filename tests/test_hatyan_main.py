# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 17:21:13 2020

@author: veenstra

"""

import pytest
import os
import numpy as np
import pandas as pd
import datetime as dt
from netCDF4 import Dataset, num2date
import hatyan

dir_tests = os.path.dirname(__file__) #F9 doesnt work, only F5 (F5 also only method to reload external definition scripts)
dir_testdata = os.path.join(dir_tests,'data_unitsystemtests')
modulename_list = ['os','sys','glob','shutil','scipy','numpy','datetime','pandas','pyproj','matplotlib','netCDF4','hatyan']


@pytest.mark.parametrize("modulename", [pytest.param('%s'%(stat), id='%s'%(stat)) for stat in modulename_list])
@pytest.mark.unittest
def test_import_libraries(modulename):
    exec('import %s'%(modulename))
    if modulename == 'pyproj':
        exec('from pyproj import Transformer')


@pytest.mark.unittest
def test_readts_dia_multifile():
    file_data_comp0 = [os.path.join(dir_testdata,'%s_obs%i.txt'%('VLISSGN', file_id)) for file_id in [1,2,3,4]]
    ts_measurements_group0 = hatyan.readts_dia(filename=file_data_comp0, station='VLISSGN')
    
    assert len(ts_measurements_group0) == 35064


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
    times_ext_comp0 = [dt.datetime(2018,1,1),dt.datetime(2018,4,1)]

    #component groups
    ts_measurements_group0 = hatyan.readts_noos(filename=file_data_comp0)
    ts_measurements_group0_res = hatyan.resample_timeseries(ts_measurements_group0, timestep_min=10)
    ts_measurements_group0_rescrop = hatyan.crop_timeseries(ts_measurements_group0_res, times_ext=times_ext_comp0)
    
    assert len(ts_measurements_group0) == 12752
    assert len(ts_measurements_group0_res) == 12961
    assert len(ts_measurements_group0_rescrop) == 12961
    assert ts_measurements_group0_rescrop.index[0] == pd.Timestamp('2018-01-01')
    assert ts_measurements_group0_rescrop.index[-1] == pd.Timestamp('2018-04-01')
    assert ts_measurements_group0_rescrop['values'][0] == 2.5
    assert ts_measurements_group0_rescrop['values'][-1] == 1.05


@pytest.mark.unittest
def test_writenetcdf():
    
    current_station = 'VLISSGN'
    times_ext=[dt.datetime(2019,1,1),dt.datetime(2019,6,1)]
    
    file_pred = os.path.join(dir_testdata,f'{current_station}_pre.txt')
    ts_prediction = hatyan.readts_dia(filename=file_pred, station=current_station)
    ts_prediction = hatyan.crop_timeseries(ts_prediction, times_ext=times_ext)
    ts_ext_prediction = hatyan.calc_HWLW(ts=ts_prediction)
    
    file_nc = 'prediction_10m_%s.nc'%(current_station)
    hatyan.write_tsnetcdf(ts=ts_prediction, ts_ext=ts_ext_prediction, station=current_station, vertref='NAP', filename=file_nc, tzone_hr=1)
    
    data_nc = Dataset(file_nc,'r')
    
    timevar = data_nc.variables['time']
    timevar_dt = num2date(timevar[:],units=timevar.units, only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    
    #put netcdf file contents in pandas DataFrame for usage in hatyan
    ts_pd = pd.DataFrame({'values':data_nc.variables['waterlevel_astro'][:,0]}, index=timevar_dt)
    print(ts_pd)
    
    assert list(data_nc.dimensions.keys()) == ['stations', 'statname_len', 'time', 'analysis_time', 'time_HW', 'time_LW']
    assert list(data_nc.variables.keys()) == ['stations', 'analysis_time', 'time', 'waterlevel_astro', 'time_HW', 'waterlevel_astro_HW', 'time_LW', 'waterlevel_astro_LW']
    assert timevar_dt[0] == times_ext[0]
    assert timevar_dt[-1] == times_ext[-1]
    assert 'title' in data_nc.__dict__.keys()
    assert 'institution' in data_nc.__dict__.keys()
    assert 'source' in data_nc.__dict__.keys()
    assert timevar.units == 'minutes since 1900-01-01 00:00:00 +0100'

    data_nc.close()
    os.remove(file_nc)


@pytest.mark.unittest
def test_analysis_settings():
    
    current_station = 'VLISSGN'
    
    #comp0
    file_data_comp0 = [os.path.join(dir_testdata,'%s_obs%i.txt'%(current_station, file_id)) for file_id in [1]]
    ts_measurements_group0 = hatyan.readts_dia(filename=file_data_comp0, station=current_station)
    
    ts_comp_nfac1_fualltimes0_xfac1_peryear0 = hatyan.get_components_from_ts(ts=ts_measurements_group0, const_list='month', nodalfactors=True, fu_alltimes=False, xfac=True, analysis_perperiod=False)
    
    ts_comp_nfac1_fualltimes1_xfac1 = hatyan.analysis(ts=ts_measurements_group0, const_list='month', nodalfactors=True, fu_alltimes=True, xfac=True)
    ts_comp_nfac1_fualltimes0_xfac1 = hatyan.analysis(ts=ts_measurements_group0, const_list='month', nodalfactors=True, fu_alltimes=False, xfac=True)
    ts_comp_nfac1_fualltimes1_xfac0 = hatyan.analysis(ts=ts_measurements_group0, const_list='month', nodalfactors=True, fu_alltimes=True, xfac=False)
    ts_comp_nfac1_fualltimes0_xfac0 = hatyan.analysis(ts=ts_measurements_group0, const_list='month', nodalfactors=True, fu_alltimes=False, xfac=False)
    ts_comp_nfac0_fualltimes0_xfac0 = hatyan.analysis(ts=ts_measurements_group0, const_list='month', nodalfactors=False, fu_alltimes=False, xfac=False)
    
    comp_A_expected = np.array([9.92687566e-04, 3.11039747e-02, 9.73652407e-02, 6.68645368e-02,
                                4.02714135e-02, 2.44247279e-02, 1.34401343e-01, 2.69848016e-01,
                                1.74658516e+00, 1.33037114e-01, 4.78882284e-01, 4.42634321e-02,
                                1.93477077e-02, 4.29425377e-02, 1.28103067e-01, 8.89165480e-02,
                                4.66780541e-02, 8.50788597e-02, 9.01287840e-02, 3.06546626e-02,
                                4.72786278e-02, 1.47588322e-02])
    
    comp_phi_expected = np.array([  0.        , 141.78046951, 188.22394202,   7.35166063,
                                  279.19176166, 122.65613222, 159.77053294,  33.54516691,
                                   59.02752661, 256.04742485, 117.45538945, 345.62388999,
                                  200.90891452,  92.80114988, 115.31391482, 176.58434752,
                                   77.32631928, 103.4830287 , 153.94793646, 110.76924584,
                                  157.3933541 , 220.41306556])
    
    assert (np.abs(ts_comp_nfac1_fualltimes0_xfac1_peryear0-ts_comp_nfac1_fualltimes0_xfac1) < 10E-9).all().all()
    assert (np.abs(ts_comp_nfac1_fualltimes0_xfac1_peryear0['phi_deg'].values-comp_phi_expected) < 10E-9).all()
    assert (np.abs(ts_comp_nfac1_fualltimes0_xfac1_peryear0['A'].values-comp_A_expected) < 10E-9).all()


@pytest.mark.unittest
def test_analysis_foreman():
    current_station = 'VLISSGN'
    
    #comp0
    file_data_comp0 = [os.path.join(dir_testdata,'%s_obs%i.txt'%(current_station, file_id)) for file_id in [1]]
    ts_measurements_group0 = hatyan.readts_dia(filename=file_data_comp0, station=current_station)
    
    comp_frommeas_SCHU = hatyan.get_components_from_ts(ts=ts_measurements_group0, const_list='month', nodalfactors=True, fu_alltimes=True, xfac=True, source='schureman')
    comp_frommeas_FOR = hatyan.get_components_from_ts(ts=ts_measurements_group0, const_list='month', nodalfactors=True, fu_alltimes=True, xfac=True, source='foreman')
    SCHU_comp_phi_expected =  np.array([  0.        , 141.68386391, 188.27101759,   7.3992028 ,
                                       279.17807514, 122.81257343, 159.77061119,  33.51788732,
                                        59.01827634, 256.05796043, 117.47528038, 345.63438447,
                                       201.05307505,  92.81986396, 115.2820895 , 176.56952991,
                                        77.30925229, 103.43837741, 153.90952097, 110.64017054,
                                       157.30038107, 220.24021917])
    SCHU_comp_A_expected =    np.array([9.93677030e-04, 3.11102518e-02, 9.74940999e-02, 6.69089179e-02,
                                       4.01204828e-02, 2.43192049e-02, 1.34386805e-01, 2.69878623e-01,
                                       1.74649803e+00, 1.32976762e-01, 4.79019453e-01, 4.42086557e-02,
                                       1.92695651e-02, 4.29786476e-02, 1.28078709e-01, 8.89513755e-02,
                                       4.67055728e-02, 8.50832246e-02, 9.01742533e-02, 3.06526796e-02,
                                       4.72907179e-02, 1.47366661e-02])
    FOR_comp_phi_expected =  np.array([  0.        , 321.7403689 ,   8.15729584, 187.40999172,
                                       279.42669414, 122.65017113, 160.05967765,  33.32750808,
                                        59.01572119, 256.18495318, 117.35826025, 345.41402817,
                                       201.17575979,  92.60563462, 115.27544986, 176.43679913,
                                        77.08771066, 103.42351654, 153.76456299, 110.600143  ,
                                       157.16101443, 220.13261581])
    FOR_comp_A_expected =    np.array([9.96058525e-04, 3.01459094e-02, 9.72007804e-02, 6.69287134e-02,
                                       3.99816712e-02, 2.42471615e-02, 1.33778086e-01, 2.75247289e-01,
                                       1.76155491e+00, 1.39007149e-01, 4.85825897e-01, 4.40858825e-02,
                                       1.92198515e-02, 4.29817009e-02, 1.29452749e-01, 9.04975848e-02,
                                       4.66816274e-02, 8.62181973e-02, 9.27486375e-02, 3.13117938e-02,
                                       4.82570962e-02, 1.46874202e-02])
    
    assert (np.abs(comp_frommeas_SCHU['phi_deg'].values-SCHU_comp_phi_expected) < 10E-9).all()
    assert (np.abs(comp_frommeas_SCHU['A'].values-SCHU_comp_A_expected) < 10E-9).all()
    assert (np.abs(comp_frommeas_FOR['phi_deg'].values-FOR_comp_phi_expected) < 10E-9).all()
    assert (np.abs(comp_frommeas_FOR['A'].values-FOR_comp_A_expected) < 10E-9).all()


@pytest.mark.unittest
def test_getcomponentsfromts_settings():
    
    current_station = 'VLISSGN'
    
    #comp0
    file_data_comp0 = [os.path.join(dir_testdata,'%s_obs%i.txt'%(current_station, file_id)) for file_id in [1,2]]
    ts_measurements_group0 = hatyan.readts_dia(filename=file_data_comp0, station=current_station)
    
    ts_comp_nfac1_fualltimes1_xfac1_peryear0 = hatyan.get_components_from_ts(ts=ts_measurements_group0, const_list='month', nodalfactors=True, fu_alltimes=True, xfac=True, analysis_perperiod=False)
    
    ts_comp_nfac1_fualltimes1_xfac1_permonth0 = hatyan.get_components_from_ts(ts=ts_measurements_group0, const_list='month', nodalfactors=True, fu_alltimes=True, xfac=True, analysis_perperiod='M')
    
    ts_comp_nfac1_fualltimes1_xfac1 = hatyan.get_components_from_ts(ts=ts_measurements_group0, const_list='month', nodalfactors=True, fu_alltimes=True, xfac=True, analysis_perperiod='Y')
    ts_comp_nfac1_fualltimes0_xfac1 = hatyan.get_components_from_ts(ts=ts_measurements_group0, const_list='month', nodalfactors=True, fu_alltimes=False, xfac=True, analysis_perperiod='Y')
    ts_comp_nfac1_fualltimes1_xfac0 = hatyan.get_components_from_ts(ts=ts_measurements_group0, const_list='month', nodalfactors=True, fu_alltimes=True, xfac=False, analysis_perperiod='Y')
    ts_comp_nfac1_fualltimes0_xfac0 = hatyan.get_components_from_ts(ts=ts_measurements_group0, const_list='month', nodalfactors=True, fu_alltimes=False, xfac=False, analysis_perperiod='Y')
    ts_comp_nfac0_fualltimes0_xfac0 = hatyan.get_components_from_ts(ts=ts_measurements_group0, const_list='month', nodalfactors=False, fu_alltimes=False, xfac=False, analysis_perperiod='Y')


@pytest.mark.unittest
def test_predictionsettings():
    times_ext_pred_HWLWno = [dt.datetime(2009,12,31,14),dt.datetime(2010,1,2,12)]
    times_step_pred = 1
    current_station = 'DENHDR'
    
    file_data_comp0 = os.path.join(dir_testdata,'%s_ana.txt'%(current_station))
    COMP_merged = hatyan.read_components(filename=file_data_comp0)
    ts_prediction_nfac1_fualltimes1_xfac1 = hatyan.prediction(comp=COMP_merged, nodalfactors=True, fu_alltimes=True, xfac=True, times_ext=times_ext_pred_HWLWno, timestep_min=times_step_pred)
    ts_prediction_nfac1_fualltimes0_xfac1 = hatyan.prediction(comp=COMP_merged, nodalfactors=True, fu_alltimes=False, xfac=True, times_ext=times_ext_pred_HWLWno, timestep_min=times_step_pred)
    ts_prediction_nfac1_fualltimes1_xfac0 = hatyan.prediction(comp=COMP_merged, nodalfactors=True, fu_alltimes=True, xfac=False, times_ext=times_ext_pred_HWLWno, timestep_min=times_step_pred)
    ts_prediction_nfac1_fualltimes0_xfac0 = hatyan.prediction(comp=COMP_merged, nodalfactors=True, fu_alltimes=False, xfac=False, times_ext=times_ext_pred_HWLWno, timestep_min=times_step_pred)
    ts_prediction_nfac0_fualltimes0_xfac0 = hatyan.prediction(comp=COMP_merged, nodalfactors=False, fu_alltimes=False, xfac=False, times_ext=times_ext_pred_HWLWno, timestep_min=times_step_pred)


@pytest.mark.unittest
def test_prediction_1018():
    current_station = 'DENHDR'
    nodalfactors = True
    xfac=True
    
    times_ext_pred = [dt.datetime(1018,7,21),dt.datetime(1018,7,21,3)]
    times_step_pred = 10
    
    file_data_comp0 = os.path.join(dir_testdata,'%s_ana.txt'%(current_station))
    COMP_merged = hatyan.read_components(filename=file_data_comp0)
    
    #prediction
    ts_prediction = hatyan.prediction(comp=COMP_merged, nodalfactors=nodalfactors, xfac=xfac, fu_alltimes=False, times_ext=times_ext_pred, timestep_min=times_step_pred)
    
    ts_prediction_times = np.array([dt.datetime(1018, 7, 21, 0, 0),
                                    dt.datetime(1018, 7, 21, 0, 10),
                                    dt.datetime(1018, 7, 21, 0, 20),
                                    dt.datetime(1018, 7, 21, 0, 30),
                                    dt.datetime(1018, 7, 21, 0, 40),
                                    dt.datetime(1018, 7, 21, 0, 50),
                                    dt.datetime(1018, 7, 21, 1, 0),
                                    dt.datetime(1018, 7, 21, 1, 10),
                                    dt.datetime(1018, 7, 21, 1, 20),
                                    dt.datetime(1018, 7, 21, 1, 30),
                                    dt.datetime(1018, 7, 21, 1, 40),
                                    dt.datetime(1018, 7, 21, 1, 50),
                                    dt.datetime(1018, 7, 21, 2, 0),
                                    dt.datetime(1018, 7, 21, 2, 10),
                                    dt.datetime(1018, 7, 21, 2, 20),
                                    dt.datetime(1018, 7, 21, 2, 30),
                                    dt.datetime(1018, 7, 21, 2, 40),
                                    dt.datetime(1018, 7, 21, 2, 50),
                                    dt.datetime(1018, 7, 21, 3, 0)], dtype=object)
    ts_prediction_vals = np.array([-0.77416669, -0.78821831, -0.79718123, -0.80126649, -0.80055319,
                                   -0.79466499, -0.78252257, -0.76226094, -0.73136835, -0.68705389,
                                   -0.62679467, -0.54896408, -0.4534113 , -0.34185962, -0.218017  ,
                                   -0.08734459,  0.04350389,  0.16749106,  0.27811846])
    
    assert (np.abs(ts_prediction['values'].values-ts_prediction_vals) < 10E-9).all()
    assert (np.abs(ts_prediction.index.values-ts_prediction_times) < dt.timedelta(days=10E-9)).all()


@pytest.mark.systemtest
def test_frommergedcomp():
    # 1. define test data
    nodalfactors = True
    xfac = True
    current_station = 'VLISSGN'
    
    #comp
    file_data_comp = os.path.join(dir_testdata,'%s_ana.txt'%(current_station))
    
    #pred
    times_ext_pred = [dt.datetime(2019,1,1),dt.datetime(2019,1,1,12)]
    timestep_pred = 10
    
    # 2. define initial expectations
    expected_ts_prediction_data_pd_values = np.array([1.00809295,  0.89565827,  0.77557688,  0.64880508,  0.51669927,
                                                      0.38083644,  0.24282757,  0.1041679 , -0.03384444, -0.17011025,
                                                     -0.3036618 , -0.43355561, -0.55876447, -0.67810779, -0.79024821,
                                                     -0.89376351, -0.98728197, -1.06965037, -1.14009357, -1.19832413,
                                                     -1.2445716 , -1.27952004, -1.30416467, -1.31961912, -1.32691704,
                                                     -1.32685401, -1.31990568, -1.30623924, -1.28581172, -1.25852648,
                                                     -1.22440367, -1.18371671, -1.1370541 , -1.08528461, -1.02942801,
                                                     -0.97045777, -0.90908055, -0.84554483, -0.77952619, -0.71012042,
                                                     -0.63595176, -0.55537797, -0.46675305, -0.36869694, -0.26032191,
                                                     -0.14137807, -0.01230058,  0.12583519,  0.27142446,  0.4224715 ,
                                                      0.57671752,  0.73172934,  0.88493925,  1.03364646,  1.17500835,
                                                      1.30605751,  1.42377775,  1.52525728,  1.60791471,  1.66976848,
                                                      1.70970034,  1.72765351,  1.72471102,  1.70301805,  1.66554205,
                                                      1.61569821,  1.55689834,  1.49210117,  1.42344542,  1.35203362,
                                                      1.27790547,  1.20020227,  1.11748521])
    
    # 3. run test
    COMP_mergedfromfile = hatyan.read_components(filename=file_data_comp)
    ts_prediction_direct = hatyan.prediction(COMP_mergedfromfile, nodalfactors=nodalfactors, xfac=xfac, fu_alltimes=False, times_ext=times_ext_pred, timestep_min=timestep_pred)
    ts_prediction_direct_values = ts_prediction_direct['values'].values
    
    # 4. Vefiry final expectations
    assert type(ts_prediction_direct) == pd.core.frame.DataFrame
    assert ts_prediction_direct.index[0].to_pydatetime() == times_ext_pred[0]
    assert ts_prediction_direct.index[-1].to_pydatetime() == times_ext_pred[-1]
    assert len(ts_prediction_direct_values) == len(expected_ts_prediction_data_pd_values)
    assert type(ts_prediction_direct_values) == type(expected_ts_prediction_data_pd_values)
    assert (np.abs(ts_prediction_direct_values - expected_ts_prediction_data_pd_values) < 10E-9).all()


@pytest.mark.systemtest
def test_frommergedcomp_HWLW_toomuch():
    """
    This test produces a very short prediction for DENHDR, based on an imported component list. It then calculates extremes (HW/LW) and numbers them both.
    This test is meant to tweak the prominence parameter of scipy.signal.find_peaks() in hatyan.calc_HWLW() and it fails when the prominence is not set or is too small.
    A prominence value of None works for most stations, but fails for DENHDR in this period since a local dip is interpreted as a low water.
    This incorrect LW has a prominence of ~0.03, other peaks in this timeseries have a prominence of >1 and even when lowering M2 to an unrealistic value (0.5*M2) the prominence of peaks is >0.2
    A prominence value of 0.1 makes sure this test succeeds.
    
    """
    
    #for testing occurance of invalid low water at start of denhelder timeseries
    times_ext_pred_HWLWno = [dt.datetime(2009,12,31,14),dt.datetime(2010,1,2,12)]
    times_step_pred = 1
    current_station = 'DENHDR'

    file_data_comp0 = os.path.join(dir_testdata,'%s_ana.txt'%(current_station))
    COMP_merged = hatyan.read_components(filename=file_data_comp0)
    COMP_merged_temp = COMP_merged.copy()
    #COMP_merged_temp.loc['M2','A']=0.05
    ts_prediction_HWLWno = hatyan.prediction(comp=COMP_merged_temp, nodalfactors=True, xfac=True, fu_alltimes=True, times_ext=times_ext_pred_HWLWno, timestep_min=times_step_pred)
    ts_ext_prediction_HWLWno_pre = hatyan.calc_HWLW(ts=ts_prediction_HWLWno, debug=True)
    
    ts_ext_prediction_HWLWno = hatyan.calc_HWLWnumbering(ts_ext=ts_ext_prediction_HWLWno_pre, station=current_station)
    #fig,(ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction_HWLWno, ts_ext=ts_ext_prediction_HWLWno_pre)
    #for irow, pdrow in ts_ext_prediction_HWLWno.iterrows():
    #    ax1.text(pdrow.index,pdrow['values'],pdrow['HWLWno'], color='k')

    assert len(ts_ext_prediction_HWLWno_pre) == 5
    assert (np.abs(ts_ext_prediction_HWLWno_pre['HWLWcode'].values-np.array([2, 1, 2, 1, 2])) < 10E-9).all()
    assert (np.abs(ts_ext_prediction_HWLWno_pre['values'].values-np.array([-0.80339484,  0.61842428, -0.76465391,  0.79758734, -0.87872493 ])) < 10E-9).all()

    assert len(ts_ext_prediction_HWLWno['HWLWno']) == 5
    assert (np.abs(ts_ext_prediction_HWLWno['HWLWno'].values -np.array([7057, 7058, 7058, 7059, 7059])) < 10E-9).all()


@pytest.mark.systemtest
#@pytest.mark.parametrize("current_station, yr", [pytest.param(x, y, id='%s %d'%(x,y))  for y in range(1999,2022) for x in ['WICK','ABDN','LEITH','WHITBY','IMMHM','CROMR','FELSWE','CADZD','VLISSGN','TERNZN','ROOMPBTN','HARVT10','HOEKVHLD','ROTTDM','DORDT','SCHEVNGN','IJMDBTHVN','PETTZD','DENHDR','DENOVBTN','HARLGN','HOLWD','SCHIERMNOG','LAUWOG','EEMSHVN','DELFZL','CUXHVN']])
@pytest.mark.parametrize("current_station, yr", [pytest.param(x, y, id='%s_%d'%(x,y)) for y in [2018,2022] for x in ['HOEKVHLD','ROTTDM','DENHDR','LITHDP']])
def test_frommergedcomp_HWLW_toolittle(current_station, yr):
    """
    This test produces a prediction for an entire year for a specific station, based on an imported component list. It then calculates extremes (HW/LW) and numbers them both. 
    The test fails if the numbers for HW and/or LW do not always increase with 1 and prints parts of the HW/LW-arrays with these gaps.
    This test is meant to tweak the distance parameter of scipy.signal.find_peaks() in hatyan.calc_HWLW() and it fails when the distance is to small.
    A distance value of M2period/1.4 works for most stations, but fails for three stations of the selection in several years (of which the year 2000)
    A distance value of M2period/1.5 for HW and M2period/1.7 for LW works for all tested stations and years.
    
    #selected stations which previously resulted in missing HWLW values in year 2000 and 2018 (and more)
    current_station = 'HOEKVHLD'
    current_station = 'ROTTDM'
    current_station = 'DENHDR'
    current_station = 'LITHDP' #still fails, but is accepted
    yr=2000
    yr=2018
    """
    
    #stats_all = ['ABDN','AMLAHVN','BAALHK','BATH','BERGSDSWT','BORSSLE','BOURNMH','BRESKS','BROUWHVSGT02','BROUWHVSGT08','CADZD','CROMR','CUXHVN','DELFZL','DENHDR','DENOVBTN','DEVPT','DORDT','DOVR','EEMHVN','EEMSHVN','EURPFM','EURPHVN','FELSWE','FISHGD','GEULHVN','GOIDSOD','GOUDBG','HAGSBNDN','HANSWT','HARLGN','HARMSBG','HARTBG','HARTHVN','HARVT10','HEESBN','HELLVSS','HOEKVHLD','HOLWD','HUIBGT','IJMDBTHVN','IJMDSMPL','IMMHM','KATSBTN','KEIZVR','KINLBVE','KORNWDZBTN','KRAMMSZWT','KRIMPADIJSL','KRIMPADLK','K13APFM','LAUWOG','LEITH','LICHTELGRE','LITHDP','LLANDNO','LOWST','MAASMSMPL','MAASSS','MAESLKRZZDE','MARLGT','MOERDK','NES','NEWHVN','NEWLN','NIEUWSTZL','NORTHSS','OOSTSDE04','OOSTSDE11','OOSTSDE14','OUDSD','OVLVHWT','PARKSS','PETTZD','PORTSMH','RAKND','ROOMPBNN','ROOMPBTN','ROTTDM','ROZBSSNZDE','ROZBSSZZDE','SCHAARVDND','SCHEVNGN','SCHIERMNOG','SCHOONHVN','SHEERNS','SINTANLHVSGR','SPIJKNSE','STAVNSE','STELLDBTN','STORNWY','SUURHBNZDE','TENNSHVN','TERNZN','TERSLNZE','TEXNZE','VLAARDGN','VLAKTVDRN','VLIELHVN','VLISSGN','VURN','WALSODN','WERKDBTN','WESTKPLE','WESTTSLG','WEYMH','WHITBY','WICK','WIERMGDN','YERSKE','ZALTBML','A12','AUKFPFM','AWGPFM','D15','F16','F3PFM','J6','K14PFM','L9PFM','NORTHCMRT','Q1']
    stats_xfac0 = ['A12','ABDN','AUKFPFM','BOURNMH','CROMR','CUXHVN','D15','DEVPT','DOVR','F16','F3PFM','FELSWE','FISHGD','IMMHM','J6','K13APFM','K14PFM','KINLBVE','LEITH','LLANDNO','LOWST','NEWHVN','NEWLN','NORTHCMRT','NORTHSS','PORTSMH','SHEERNS','STORNWY','WEYMH','WHITBY','WICK']
    #stats_anaperyear0 = ['A12','ABDN','AUKFPFM','BOURNMH','CROMR','D15','DEVPT','DOVR','F16','F3PFM','FELSWE','FISHGD','IMMHM','J6','K14PFM','KINLBVE','LEITH','LLANDNO','LOWST','NEWHVN','NEWLN','NORTHCMRT','NORTHSS','PORTSMH','SHEERNS','STORNWY','WEYMH','WHITBY','WICK']
    #stats_MSL = ['EURPFM','K13APFM','LICHTELGRE','A12','AUKFPFM','AWGPFM','D15','F16','F3PFM','J6','K14PFM','L9PFM','NORTHCMRT','Q1']
    
    #selected_stations = stats_all
    #selected_stations = ['WICK','ABDN','LEITH','WHITBY','IMMHM','CROMR','FELSWE','CADZD','VLISSGN','TERNZN','ROOMPBTN','HARVT10','HOEKVHLD','ROTTDM','DORDT','SCHEVNGN','IJMDBTHVN','PETTZD','DENHDR','DENOVBTN','HARLGN','HOLWD','SCHIERMNOG','LAUWOG','EEMSHVN','DELFZL','CUXHVN']
    #selected_stations = ['CROMR','CADZD','HOEKVHLD','DENHDR','CUXHVN']
    #selected_stations = ['HOEKVHLD','ROTTDM','DENHDR'] #selected stations which resulted in missing HWLW values 
    #for yr in [2000]: #range(1999,2022):
    #fig, (ax1) = plt.subplots(1,1,figsize=(15,6))
    #n_colors = len(selected_stations)
    #colors = plt.cm.jet(np.linspace(0,1,n_colors))
    #for i_stat, current_station in enumerate(selected_stations):
    print('-'*100)
    print('%-45s = %s'%('station_name',current_station))
    print('-'*5)
    
    #START OF STATION SETTINGS
    #xfactor
    if current_station in stats_xfac0:
        xfac=False
    else:
        xfac=True
    #analysis_peryear
    #if current_station in stats_anaperyear0:
    #    analysis_perperiod=False
    #else:
    #    analysis_perperiod='Y'
    #constituent list
    #const_list = hatyan.get_const_list_hatyan('year') #94 const
    #vertical reference
    #vertref='NAP'
    #END OF STATION SETTINGS
    
    file_data_comp0 = os.path.join(dir_testdata,'%s_ana.txt'%(current_station))
    times_ext_pred = [dt.datetime(yr,1,1),dt.datetime(yr+1,1,1)]
    times_step_pred = 1
    
    #component groups
    COMP_merged = hatyan.read_components(filename=file_data_comp0)
    
    #prediction and validation
    ts_prediction = hatyan.prediction(comp=COMP_merged, nodalfactors=True, xfac=xfac, fu_alltimes=False, times_ext=times_ext_pred, timestep_min=times_step_pred)
    #ts_validation = hatyan.readts_dia(filename=file_data_predvali, station=current_station)
    #ts_ext_validation = hatyan.readts_dia(filename=file_data_predvaliHWLW, station=current_station)
    ts_ext_prediction = hatyan.calc_HWLW(ts=ts_prediction, debug=True)
    
    #calculate tidal wave number
    ts_ext_prediction_HWLWno = hatyan.calc_HWLWnumbering(ts_ext=ts_ext_prediction, station=current_station)
    #print(ts_ext_prediction_HWLWno)
    #for irow, pdrow in ts_ext_prediction_HWLWno.iterrows():
    #    ax1.text(pdrow.index,pdrow['values'],pdrow['HWLWno'], color=colors[i_stat])
    
    #ax1.plot(ts_prediction.index, ts_prediction['values'], label=current_station, color=colors[i_stat])
    #figa, (ax1a) = hatyan.plot_timeseries(ts=ts_prediction, ts_ext=ts_ext_prediction_HWLWno)
    
    #str_combination = '%i_%s'%(yr, current_station)
    #print(str_combination)
    print(current_station, yr)
    print('all HWLW values:')
    print(ts_ext_prediction_HWLWno)
    
    HW_data = ts_ext_prediction_HWLWno[ts_ext_prediction_HWLWno['HWLWcode']==1]
    HW_data_diff1bool = (HW_data['HWLWno'].diff().iloc[1:].values==1)
    print('%i parts of HW array containing gaps:'%((~HW_data_diff1bool).sum()))
    if not HW_data_diff1bool.all():
        ids_false = np.where(~HW_data_diff1bool)[0]
        for id_false in ids_false: 
            print(HW_data.iloc[id_false-1:id_false+3])
    
    LW_data = ts_ext_prediction_HWLWno[ts_ext_prediction_HWLWno['HWLWcode']==2]
    LW_data_diff1bool = (LW_data['HWLWno'].diff().iloc[1:].values==1)
    print('%i parts of LW array containing gaps:'%((~LW_data_diff1bool).sum()))
    if not LW_data_diff1bool.all():
        ids_false = np.where(~LW_data_diff1bool)[0]
        for id_false in ids_false: 
            print(LW_data.iloc[id_false-1:id_false+3])
    
    assert HW_data_diff1bool.all()
    assert LW_data_diff1bool.all()
    
    #ax1.set_xlim(times_ext_pred)
    #fig.tight_layout()
    #ax1.legend(loc=2, fontsize=7)#bbox_to_anchor=(1,1))
    #import matplotlib.dates as mdates
    #ax1.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d %H:%M"))
    #fig.savefig('tide_numbering_%i.png'%(yr), dpi=250)


@pytest.mark.parametrize("current_station", [pytest.param(x, id=x) for x in ['HOEKVHLD','DENHDR']])
@pytest.mark.systemtest
def test_frommergedcomp_HWLW_345(current_station):
    """
    This test produces a prediction for a period for HOEKVHLD, based on an imported component list. It then calculates extremes (HW/LW) with several settings and numbers them both (including 3/4/5 HWLWcodes around aggers, excluding 11/22 HWLWcodes). 
    #HOEKVHLD has alternating aggers, DENHDR has double HW's (which should not be numbered as aggers)
    
    current_station = 'HOEKVHLD'
    current_station = 'DENHDR'
    """
    
    times_step_pred = 1
    
    file_data_comp0 = os.path.join(dir_testdata,'%s_ana.txt'%(current_station))
    COMP_merged = hatyan.read_components(filename=file_data_comp0)
    COMP_merged_temp = COMP_merged.copy()
    
    times_ext_pred_HWLWno = [dt.datetime(2010,1,31,3),dt.datetime(2010,2,17,12)] #longer period with alternating aggers and no aggers, also eerste HW wordt als lokaal ipv primair HW gezien, also extra agger outside of 1stLW/agger/2ndLW sequence
    #times_ext_pred_HWLWno = [dt.datetime(2010,1,31),dt.datetime(2010,2,3)] #extra agger outside of 1stLW/agger/2ndLW sequence
    #times_ext_pred_HWLWno = [dt.datetime(2010,6,26),dt.datetime(2010,6,28)] #lokale extremen tussen twee laagste LWs
    #times_ext_pred_HWLWno = [dt.datetime(2019,2,1),dt.datetime(2019,2,2)] #eerste HW wordt als lokaal ipv primair HW gezien (lage prominence door dicht op begin tijdserie) >> warning
    
    ts_prediction_HWLWno = hatyan.prediction(comp=COMP_merged_temp, nodalfactors=True, xfac=True, fu_alltimes=True, times_ext=times_ext_pred_HWLWno, timestep_min=times_step_pred)
    ts_ext_prediction_main = hatyan.calc_HWLW(ts=ts_prediction_HWLWno)#, debug=True)
    #ts_ext_prediction_all = hatyan.calc_HWLW(ts=ts_prediction_HWLWno, calc_HWLW345=True, calc_HWLW345_cleanup1122=False)#, debug=True)
    ts_ext_prediction_clean = hatyan.calc_HWLW(ts=ts_prediction_HWLWno, calc_HWLW345=True)#, calc_HWLW345_cleanup1122=True) #for numbering, cannot cope with 11/22 HWLWcodes
    
    ts_ext_prediction_main_HWLWno = hatyan.calc_HWLWnumbering(ts_ext=ts_ext_prediction_main, station=current_station)
    ts_ext_prediction_clean_HWLWno = hatyan.calc_HWLWnumbering(ts_ext=ts_ext_prediction_clean, station=current_station)
    
    #fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction_HWLWno, ts_ext=ts_ext_prediction_all)#, ts_ext_validation=ts_ext_validation)
    #for irow, pdrow in ts_ext_prediction_clean_HWLWno.iterrows():
    #    ax1.text(pdrow.index,pdrow['values'],pdrow['HWLWno'])
    
    #expected values
    if current_station == 'HOEKVHLD':
        expected_ts_ext_prediction_main_HWLWno_HWLWno = np.array([7117, 7117, 7118, 7118, 7119, 7119, 7120, 7120, 7121,
                                                                  7121, 7122, 7122, 7123, 7123, 7124, 7124, 7125, 7125, 7126, 7126,
                                                                  7127, 7127, 7128, 7128, 7129, 7129, 7130, 7130, 7131, 7131, 7132,
                                                                  7132, 7133, 7133, 7134, 7134, 7135, 7135, 7136, 7136, 7137, 7137,
                                                                  7138, 7138, 7139, 7139, 7140, 7140, 7141, 7141, 7142, 7142, 7143,
                                                                  7143, 7144, 7144, 7145, 7145, 7146, 7146, 7147, 7147, 7148, 7148,
                                                                  7149])
        expected_ts_ext_prediction_main_HWLWcode = np.array([1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2,
                                                             1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2,
                                                             1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2,
                                                             1])
        """
        expected_ts_ext_prediction_all_HWLWcode = np.array([11, 22, 11, 22,  1,  3,  4,  5,  1,  3,  4,  5,  1,  3,  4,  5,  1,
                                                             3,  4,  5,  1,  3,  4,  5,  1,  3,  4,  5,  1,  3,  4,  5,  1,  3,
                                                             4,  5,  1,  3,  4,  5,  1,  3,  4,  5,  1,  2,  1,  2,  1,  2,  1,
                                                             2,  1,  2,  1,  2,  1,  2,  1,  2,  1,  2,  1,  2,  1,  3,  4,  5,
                                                             1,  2,  1,  3,  4,  5,  1,  2,  1,  3,  4,  5,  1,  2,  1,  3,  4,
                                                             5,  1,  2,  1,  2,  1,  3,  4,  5,  1,  3,  4,  5,  1,  3,  4,  5,
                                                             1, 22, 11])
        """
        expected_ts_ext_prediction_clean_HWLWno_HWLWcode = np.array([1, 3, 4, 5, 1, 3, 4, 5, 1, 3, 4, 5, 1, 3, 4, 5, 1, 3, 4, 5, 1, 3,
                                                                     4, 5, 1, 3, 4, 5, 1, 3, 4, 5, 1, 3, 4, 5, 1, 3, 4, 5, 1, 2, 1, 2,
                                                                     1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 3, 4, 5, 1, 2,
                                                                     1, 3, 4, 5, 1, 2, 1, 3, 4, 5, 1, 2, 1, 3, 4, 5, 1, 2, 1, 2, 1, 3,
                                                                     4, 5, 1, 3, 4, 5, 1, 3, 4, 5, 1])
        expected_ts_ext_prediction_clean_HWLWno_HWLWno = np.array([7117, 7117, 7117, 7117, 7118, 7118, 7118, 7118, 7119, 7119, 7119,
                                                                   7119, 7120, 7120, 7120, 7120, 7121, 7121, 7121, 7121, 7122, 7122,
                                                                   7122, 7122, 7123, 7123, 7123, 7123, 7124, 7124, 7124, 7124, 7125,
                                                                   7125, 7125, 7125, 7126, 7126, 7126, 7126, 7127, 7127, 7128, 7128,
                                                                   7129, 7129, 7130, 7130, 7131, 7131, 7132, 7132, 7133, 7133, 7134,
                                                                   7134, 7135, 7135, 7136, 7136, 7137, 7137, 7137, 7137, 7138, 7138,
                                                                   7139, 7139, 7139, 7139, 7140, 7140, 7141, 7141, 7141, 7141, 7142,
                                                                   7142, 7143, 7143, 7143, 7143, 7144, 7144, 7145, 7145, 7146, 7146,
                                                                   7146, 7146, 7147, 7147, 7147, 7147, 7148, 7148, 7148, 7148, 7149])
    elif current_station == 'DENHDR':
        expected_ts_ext_prediction_main_HWLWno_HWLWno = np.array([7116, 7116, 7117, 7117, 7118, 7118, 7119, 7119, 7120, 7120, 7121,
                                                                  7121, 7122, 7122, 7123, 7123, 7124, 7124, 7125, 7125, 7126, 7126,
                                                                  7127, 7127, 7128, 7128, 7129, 7129, 7130, 7130, 7131, 7131, 7132,
                                                                  7132, 7133, 7133, 7134, 7134, 7135, 7135, 7136, 7136, 7137, 7137,
                                                                  7138, 7138, 7139, 7139, 7140, 7140, 7141, 7141, 7142, 7142, 7143,
                                                                  7143, 7144, 7144, 7145, 7145, 7146, 7146, 7147, 7147, 7148, 7148])
        expected_ts_ext_prediction_main_HWLWcode = np.array([1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2,
                                                             1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2,
                                                             1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2])
        """
        expected_ts_ext_prediction_all_HWLWcode = np.array([11, 22,  1,  2, 11, 22,  1,  2, 11, 22,  1,  2, 11, 22,  1,  2, 11,
                                                            22,  1,  2, 11, 22,  1,  2, 11, 22,  1,  2, 11, 22,  1,  2,  1, 22,
                                                            11,  2,  1,  2,  1,  2,  1,  2,  1,  2,  1,  2,  1,  2,  1,  2,  1,
                                                             2,  1,  2,  1,  2,  1,  2,  1,  2,  1, 22, 11,  2,  1,  2, 11, 22,
                                                             1,  2,  1,  2, 11, 22,  1,  2,  1,  2,  1, 22, 11,  2,  1,  2,  1,
                                                             2,  1,  2,  1,  2,  1,  2, 11])
        """
        expected_ts_ext_prediction_clean_HWLWno_HWLWcode = np.array([1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2,
                                                                     1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2,
                                                                     1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2])
        expected_ts_ext_prediction_clean_HWLWno_HWLWno = np.array([7116, 7116, 7117, 7117, 7118, 7118, 7119, 7119, 7120, 7120, 7121,
                                                                   7121, 7122, 7122, 7123, 7123, 7124, 7124, 7125, 7125, 7126, 7126,
                                                                   7127, 7127, 7128, 7128, 7129, 7129, 7130, 7130, 7131, 7131, 7132,
                                                                   7132, 7133, 7133, 7134, 7134, 7135, 7135, 7136, 7136, 7137, 7137,
                                                                   7138, 7138, 7139, 7139, 7140, 7140, 7141, 7141, 7142, 7142, 7143,
                                                                   7143, 7144, 7144, 7145, 7145, 7146, 7146, 7147, 7147, 7148, 7148])
    
    assert (np.abs(ts_ext_prediction_main_HWLWno['HWLWno'].values - expected_ts_ext_prediction_main_HWLWno_HWLWno) < 10E-9).all()
    assert (np.abs(ts_ext_prediction_main['HWLWcode'].values - expected_ts_ext_prediction_main_HWLWcode) < 10E-9).all()
    #assert (np.abs(ts_ext_prediction_all['HWLWcode'].values - expected_ts_ext_prediction_all_HWLWcode) < 10E-9).all()
    assert (np.abs(ts_ext_prediction_clean_HWLWno['HWLWcode'].values - expected_ts_ext_prediction_clean_HWLWno_HWLWcode) < 10E-9).all()
    assert (np.abs(ts_ext_prediction_clean_HWLWno['HWLWno'].values - expected_ts_ext_prediction_clean_HWLWno_HWLWno) < 10E-9).all()


@pytest.mark.systemtest
def test_19Ycomp4Ydia():
    # 1. define test data
    nodalfactors = True
    xfac = True
    const_list = hatyan.get_const_list_hatyan('year')
    current_station = 'VLISSGN'
    
    #comp0
    file_data_comp0 = [os.path.join(dir_testdata,'%s_obs%i.txt'%(current_station, file_id)) for file_id in [1,2,3,4]]
    ts_measurements_group0 = hatyan.readts_dia(filename=file_data_comp0, station=current_station)
    comp_frommeasurements_avg_group0 = hatyan.get_components_from_ts(ts=ts_measurements_group0, const_list=const_list, nodalfactors=nodalfactors, xfac=xfac, fu_alltimes=False, analysis_perperiod='Y')
    comp_frommeasurements_avg_group0.station = current_station
    
    #comp1
    file_data_comp1 = os.path.join(dir_testdata,'%s_ana.txt'%(current_station))
    comp_fromfile_group1 = hatyan.read_components(filename=file_data_comp1)
    
    #merge component groups (SA/SM from 19Y, rest from 4Y)
    COMP_merged = hatyan.merge_componentgroups(comp_main=comp_frommeasurements_avg_group0, comp_sec=comp_fromfile_group1, comp_sec_list=['SA','SM'])
    
    #prediction and validation
    times_ext_pred = [dt.datetime(2019,1,1),dt.datetime(2019,1,1,12)]
    times_step_pred = 10
    ts_prediction = hatyan.prediction(COMP_merged, nodalfactors=nodalfactors, xfac=xfac, fu_alltimes=False, times_ext=times_ext_pred, timestep_min=times_step_pred)
    ts_prediction_values = ts_prediction['values'].values
    
    # 2. define initial expectations (VLISSGN)
    expected_ts_prediction_data_pd_values = np.array([1.00071654,  0.88827903,  0.76819314,  0.64141558,  0.5093034 ,
                                                      0.37343441,  0.23542045,  0.09675753, -0.04125568, -0.17751967,
                                                     -0.31106673, -0.44095367, -0.56615373, -0.68548689, -0.79761632,
                                                     -0.9011202 , -0.99462701, -1.07698355, -1.14741452, -1.2056322 ,
                                                     -1.25186587, -1.28679942, -1.31142807, -1.32686565, -1.33414624,
                                                     -1.33406597, -1.3271011 , -1.31341934, -1.29297812, -1.2656809 ,
                                                     -1.23154773, -1.19085167, -1.14418076, -1.09240329, -1.03653868,
                                                     -0.97756028, -0.91617494, -0.85263162, -0.78660665, -0.71719668,
                                                     -0.64302677, -0.56245529, -0.47383656, -0.37579035, -0.26742836,
                                                     -0.14849976, -0.01943851,  0.11868132,  0.26425605,  0.41529083,
                                                      0.56952725,  0.72453204,  0.8777369 ,  1.02644009,  1.16779776,
                                                      1.29884123,  1.4165532 ,  1.51802108,  1.60066312,  1.6624979 ,
                                                      1.70240777,  1.720337  ,  1.71736989,  1.69565303,  1.65815518,
                                                      1.60829263,  1.54947792,  1.48467006,  1.41600766,  1.3445927 ,
                                                      1.27046406,  1.19276203,  1.11004681])

    # 4. Vefiry final expectations
    assert (np.abs(ts_prediction_values - expected_ts_prediction_data_pd_values) < 10E-9).all()


@pytest.mark.systemtest
def test_19Ycomp4Ydia_compsplitsing():
    # 1. define test data
    nodalfactors = True
    xfac = False
    const_list = hatyan.get_const_list_hatyan('month')
    current_station = 'D15'
    file_data_comp0 = os.path.join(dir_testdata,'%s_obs1.txt'%(current_station))
    CS_comps = pd.DataFrame({'CS_comps_from':['K1','N2','2MN2','S2','S2'],
                             'CS_comps_derive':['P1','NU2','LABDA2','K2','T2'],
                             'CS_ampfacs':[0.33,0.22,0.48,0.29,0.05],
                             'CS_degincrs':[-11,-24,174,1,-24]})
    file_data_comp1 = os.path.join(dir_testdata,'%s_ana.txt'%(current_station))
    times_ext_pred = [dt.datetime(2019,1,1),dt.datetime(2019,1,1,12)]
    times_step_pred = 10
    
    # 3. run test
    #component groups
    ts_measurements_group0 = hatyan.readts_dia(filename=file_data_comp0, station=current_station)
    comp_frommeasurements_avg_group0 = hatyan.get_components_from_ts(ts=ts_measurements_group0, nodalfactors=nodalfactors, xfac=xfac, const_list=const_list, analysis_perperiod=False, fu_alltimes=False, CS_comps=CS_comps)
    comp_fromfile_group1 = hatyan.read_components(filename=file_data_comp1)
    
    #merge component groups (SA/SM from 19Y, rest from 4Y)
    COMP_merged = hatyan.merge_componentgroups(comp_main=comp_frommeasurements_avg_group0, comp_sec=comp_fromfile_group1, comp_sec_list=['SA','SM'])
    
    #prediction and validation
    ts_prediction = hatyan.prediction(COMP_merged, nodalfactors=nodalfactors, xfac=xfac, fu_alltimes=False, times_ext=times_ext_pred, timestep_min=times_step_pred)
    ts_prediction_values = ts_prediction['values'].values
    
    # 2. define initial expectations (D15)
    expected_ts_prediction_data_pd_values = np.array([0.21208738,  0.26853271,  0.32342168,  0.37563977,  0.4242852 ,
                                                      0.46874611,  0.50873038,  0.54424431,  0.57552477,  0.60293749,
                                                      0.6268598 ,  0.64756857,  0.66515311,  0.67946822,  0.69013541,
                                                      0.69659205,  0.6981796 ,  0.69425554,  0.68430933,  0.66806236,
                                                      0.64553448,  0.61706564,  0.58328909,  0.54506078,  0.50335705,
                                                      0.45915817,  0.41333717,  0.36657245,  0.31929793,  0.27169799,
                                                      0.22374628,  0.17527983,  0.12609415,  0.0760412 ,  0.02511229,
                                                     -0.0265089 , -0.07843269, -0.1300906 , -0.18079823, -0.22984136,
                                                     -0.27656929, -0.32047779, -0.36126658, -0.39886054, -0.43339088,
                                                     -0.46513999, -0.49446005, -0.52168068, -0.54702266, -0.57053377,
                                                     -0.59205869, -0.61124872, -0.62760946, -0.64057802, -0.64961524,
                                                     -0.65429615, -0.65438144, -0.64985669, -0.64093152, -0.62799831,
                                                     -0.61155776, -0.59212467, -0.57013128, -0.54584644, -0.5193262 ,
                                                     -0.49040627, -0.45873912, -0.42387076, -0.38534492, -0.3428169 ,
                                                     -0.29615745, -0.24552815, -0.19141437])
    
    # 4. Vefiry final expectations
    assert (np.abs(ts_prediction_values - expected_ts_prediction_data_pd_values) < 10E-9).all()        


@pytest.mark.systemtest
def test_allfromdia():
    # 1. define test data
    nodalfactors = True
    xfac = True
    const_list = hatyan.get_const_list_hatyan('year')
    current_station = 'VLISSGN'
    file_data_comp0 = [os.path.join(dir_testdata,'%s_obs%i.txt'%(current_station, file_id)) for file_id in [1,2,3,4]]
    file_data_comp1 = os.path.join(dir_testdata,'%s_obs19.txt'%(current_station))
    times_ext_pred = [dt.datetime(2019,1,1),dt.datetime(2019,1,1,12)]
    times_step_pred = 10

    #component groups
    ts_measurements_group0 = hatyan.readts_dia(filename=file_data_comp0, station=current_station)
    comp_frommeasurements_avg_group0 = hatyan.get_components_from_ts(ts=ts_measurements_group0, nodalfactors=nodalfactors, xfac=xfac, const_list=const_list, analysis_perperiod='Y', fu_alltimes=False)
    ts_measurements_group1 = hatyan.readts_dia(filename=file_data_comp1, station=current_station)
    comp_frommeasurements_avg_group1 = hatyan.get_components_from_ts(ts=ts_measurements_group1, nodalfactors=nodalfactors, xfac=xfac, const_list=const_list, analysis_perperiod=False, fu_alltimes=False)
    
    #merge component groups (SA/SM from 19Y, rest from 4Y)
    COMP_merged = hatyan.merge_componentgroups(comp_main=comp_frommeasurements_avg_group0, comp_sec=comp_frommeasurements_avg_group1, comp_sec_list=['SA','SM'])
    
    #prediction and validation
    ts_prediction = hatyan.prediction(COMP_merged, nodalfactors=nodalfactors, xfac=xfac, fu_alltimes=False, times_ext=times_ext_pred, timestep_min=times_step_pred)
    ts_prediction_values = ts_prediction['values'].values
    
    # 2. define initial expectations
    """expected_ts_prediction_data_pd_values_CADZD = np.array([  0.69689822,  0.58122683,  0.46063218,  0.3366989 ,  0.21084484,
                                                        0.08416485, -0.04259883, -0.16895802, -0.29448491, -0.41858223,
                                                       -0.54028668, -0.65816757, -0.77035852, -0.87472565, -0.96913982,
                                                       -1.05179199, -1.12147718, -1.17777783, -1.22109997, -1.25255039,
                                                       -1.27368164, -1.28616356, -1.29145829, -1.29057458, -1.28395729,
                                                       -1.27153413, -1.25290186, -1.22759891, -1.19538913, -1.15647792,
                                                       -1.11159891, -1.06194149, -1.00893047, -0.95390715, -0.89778804,
                                                       -0.84078547, -0.78226211, -0.72076111, -0.65421307, -0.5802802 ,
                                                       -0.49676665, -0.40200983, -0.29517405, -0.17639266, -0.04674223,
                                                        0.09192805,  0.23725639,  0.38664551,  0.53748768,  0.68728945,
                                                        0.83368193,  0.97433907,  1.10685673,  1.22865962,  1.33699851,
                                                        1.42907551,  1.50229865,  1.55462575,  1.58492517,  1.59326408,
                                                        1.58104106,  1.55090766,  1.50646777,  1.45179364,  1.39084215,
                                                        1.32688186,  1.26204443,  1.1970914 ,  1.13144293,  1.06345986,
                                                        0.99091634,  0.91155995,  0.82363968])
    """
    expected_ts_prediction_data_pd_values = np.array([1.00064431,  0.8882068 ,  0.76812092,  0.64134336,  0.50923119,
                                                      0.3733622 ,  0.23534825,  0.09668533, -0.04132787, -0.17759185,
                                                     -0.31113891, -0.44102584, -0.5662259 , -0.68555906, -0.79768847,
                                                     -0.90119235, -0.99469915, -1.07705569, -1.14748665, -1.20570432,
                                                     -1.25193799, -1.28687153, -1.31150017, -1.32693775, -1.33421834,
                                                     -1.33413806, -1.32717318, -1.31349142, -1.2930502 , -1.26575297,
                                                     -1.23161979, -1.19092373, -1.14425281, -1.09247534, -1.03661071,
                                                     -0.97763231, -0.91624696, -0.85270364, -0.78667866, -0.71726869,
                                                     -0.64309876, -0.56252728, -0.47390854, -0.37586232, -0.26750033,
                                                     -0.14857172, -0.01951046,  0.11860937,  0.26418411,  0.4152189 ,
                                                      0.56945533,  0.72446013,  0.87766499,  1.02636819,  1.16772586,
                                                      1.29876934,  1.41648132,  1.51794921,  1.60059126,  1.66242605,
                                                      1.70233593,  1.72026516,  1.71729806,  1.6955812 ,  1.65808337,
                                                      1.60822082,  1.54940611,  1.48459827,  1.41593587,  1.34452092,
                                                      1.27039229,  1.19269026,  1.10997505])
    
    # 4. Vefiry final expectations
    assert (np.abs(ts_prediction_values - expected_ts_prediction_data_pd_values) < 10E-9).all()


@pytest.mark.systemtest
def test_allfromdia_2008xfac0():
    # 1. define test data
    nodalfactors = True
    xfac = False
    const_list = hatyan.get_const_list_hatyan('year')
    current_station = 'DOVR'
    file_data_comp0 = os.path.join(dir_testdata,'%s_obs1.txt'%(current_station))
    times_ext_pred = [dt.datetime(2019,1,1),dt.datetime(2019,1,1,12)]
    times_step_pred = 10
    
    #component groups
    ts_measurements_group0 = hatyan.readts_dia(filename=file_data_comp0, station=current_station)
    comp_frommeasurements_avg_group0 = hatyan.get_components_from_ts(ts=ts_measurements_group0, nodalfactors=nodalfactors, xfac=xfac, const_list=const_list, analysis_perperiod='Y', fu_alltimes=False)
    
    #prediction and validation
    ts_prediction = hatyan.prediction(comp=comp_frommeasurements_avg_group0, nodalfactors=nodalfactors, xfac=xfac, fu_alltimes=False, times_ext=times_ext_pred, timestep_min=times_step_pred)
    ts_prediction_values = ts_prediction['values'].values
    
    # 2. define initial expectations
    expected_ts_prediction_data_pd_values = np.array([-0.60440198, -0.72142681, -0.83581976, -0.94743404, -1.05609766,
                                                      -1.1615773 , -1.2635405 , -1.36152022, -1.45488748, -1.54283832,
                                                      -1.62440113, -1.69846786, -1.76384908, -1.81934804, -1.86384401,
                                                      -1.89637159, -1.91618096, -1.92276574, -1.91584961, -1.89533042,
                                                      -1.86118924, -1.81338063, -1.75172686, -1.67584147, -1.58510561,
                                                      -1.4787132 , -1.35578995, -1.21557736, -1.05765978, -0.88220159,
                                                      -0.69015614, -0.48340872, -0.26482336, -0.03817622,  0.19202489,
                                                       0.42081806,  0.64312991,  0.8541411 ,  1.04962903,  1.22624559,
                                                       1.38169767,  1.514813  ,  1.62549027,  1.7145481 ,  1.78350005,
                                                       1.83429032,  1.86902593,  1.88973797,  1.89819605,  1.89578968,
                                                       1.88347976,  1.86181361,  1.83099019,  1.79095846,  1.74153131,
                                                       1.68249919,  1.61373091,  1.53525334,  1.44730507,  1.35036266,
                                                       1.24514015,  1.13256443,  1.01373016,  0.8898392 ,  0.76213036,
                                                       0.63180616,  0.49996392,  0.36753791,  0.23525891,  0.10363497,
                                                      -0.02704478, -0.15668462, -0.28533886])
    
    # 4. Vefiry final expectations
    assert (np.abs(ts_prediction_values - expected_ts_prediction_data_pd_values) < 10E-9).all()


@pytest.mark.acceptance
def test_DDL_QCvalues():
    tstart_dt = dt.datetime(2019,10,1)
    tstop_dt = dt.datetime(2019,10,10)
    tzone = 'UTC+00:00' #'UTC+00:00' for GMT and 'UTC+01:00' for MET
    
    catalog_dict = hatyan.get_DDL_catalog(catalog_extrainfo=['WaardeBepalingsmethoden','MeetApparaten','Typeringen'])
    
    #HARVT10
    cat_locatielijst = catalog_dict['LocatieLijst']#.set_index('Locatie_MessageID',drop=True)
    station_dict = cat_locatielijst[cat_locatielijst['Code']=='HARVT10'].iloc[0]
    ts_meas_pd, metadata, stationdata = hatyan.get_DDL_data(station_dict=station_dict,tstart_dt=tstart_dt,tstop_dt=tstop_dt,tzone=tzone,
                                                            meta_dict={'Grootheid.Code':'WATHTE','Groepering.Code':'NVT','WaardeBewerkingsmethode.Code':'NVT'})
    uniqueQC = ts_meas_pd['QC'].unique() #array([ 0., 99., nan, 25.]), but should be integers without nan array([0, 99, 25])
    assert uniqueQC.dtype=='float64' #this should be int in the future, if None/nan is not in QC list anymore
    assert np.isnan(uniqueQC).sum() == 1 #this one should become 0 in the future and then the second assertion should be valid without indexing
    assert (uniqueQC[~np.isnan(uniqueQC)] == np.array([ 0, 99, 25])).all()

    #STELLDBTN
    station_dict = cat_locatielijst[cat_locatielijst['Code']=='STELLDBTN'].iloc[0]
    ts_meas_pd, metadata, stationdata = hatyan.get_DDL_data(station_dict=station_dict,tstart_dt=tstart_dt,tstop_dt=tstop_dt,tzone=tzone,
                                                            meta_dict={'Grootheid.Code':'WATHTE','Groepering.Code':'NVT','WaardeBewerkingsmethode.Code':'NVT'})
    uniqueQC = ts_meas_pd['QC'].unique() #array([ 0, 25], dtype=int8)
    assert uniqueQC.dtype=='int8'
    assert (uniqueQC == np.array([ 0, 25])).all()

