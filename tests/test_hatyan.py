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
import hatyan
import pytz

dir_tests = os.path.dirname(__file__) #F9 doesnt work, only F5 (F5 also only method to reload external definition scripts)
dir_testdata = os.path.join(dir_tests,'data_unitsystemtests')


@pytest.mark.systemtest
def test_frommergedcomp():
    # 1. define test data
    current_station = 'VLISSGN'
    
    #comp
    file_data_comp = os.path.join(dir_testdata,'%s_ana.txt'%(current_station))
    
    #pred
    times_pred = slice(dt.datetime(2019,1,1),dt.datetime(2019,1,1,12), "10min")
    
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
    ts_prediction_direct = hatyan.prediction(COMP_mergedfromfile, times=times_pred)
    ts_prediction_direct_values = ts_prediction_direct['values'].values
    
    # 4. Vefiry final expectations
    assert type(ts_prediction_direct) == pd.core.frame.DataFrame
    assert ts_prediction_direct.index.tz==pytz.FixedOffset(60)
    assert ts_prediction_direct.index[0].tz_localize(None) == times_pred.start
    assert ts_prediction_direct.index[-1].tz_localize(None) == times_pred.stop
    assert len(ts_prediction_direct_values) == len(expected_ts_prediction_data_pd_values)
    assert type(ts_prediction_direct_values) == type(expected_ts_prediction_data_pd_values)
    assert (np.abs(ts_prediction_direct_values - expected_ts_prediction_data_pd_values) < 10E-9).all()


@pytest.mark.systemtest
def test_meas_HWLW_toomuch():
    """
    this test will fail if the minimal prominence is set to 0.01 or lower, or of the minimal width=None
    Then there will be an additional LW found right after the HW. With a slightly higher minimal prominence (>0.02) this is avoided.
    This testcase originates from an issue with a FEWS timeseries: https://github.com/Deltares/hatyan/issues/85
    """
    wl_vals = np.array([ 0.4437,  0.563 ,  0.6677,  0.7498,  0.8125,  0.8732,  0.9181,
            0.9384,  0.944 ,  0.94  ,  0.9254,  0.9033,  0.8602,  0.8258,
            0.8177,  0.7867,  0.7184,  0.6672,  0.6314,  0.5722,  0.5028,
            0.428 ,  0.3539,  0.2822,  0.2003,  0.1141,  0.0317, -0.0457,
           -0.1233, -0.2187, -0.3101, -0.3945, -0.476 , -0.5503, -0.6137,
           -0.671 , -0.7262, -0.7742, -0.8025, -0.8145, -0.8271, -0.829 ,
           -0.8131, -0.7856, -0.7616, -0.7289, -0.7176, -0.7114, -0.7088,
           -0.7111, -0.7113, -0.7106, -0.715 , -0.7282, -0.7392, -0.7434,
           -0.7384, -0.7405, -0.7321, -0.7202, -0.7054, -0.6847, -0.6514,
           -0.6082, -0.5656, -0.5155, -0.4549, -0.3814, -0.2934, -0.1994,
           -0.0949,  0.0291,  0.1794,  0.3537,  0.5279,  0.6908,  0.8318,
            0.9402,  1.024 ,  1.0951,  1.1474,  1.1639,  1.1582,  1.147 ,
            1.1298,  1.0997,  1.0717,  1.0828,  1.078 ,  1.0232,  0.9849,
            0.969 ,  0.9424,  0.9058,  0.8601,  0.8045,  0.7569,  0.6917,
            0.6198,  0.5531,  0.4802,  0.4109,  0.3291,  0.2442,  0.1692,
            0.0838, -0.0073, -0.088 , -0.1659, -0.2412, -0.3055, -0.3647,
           -0.405 , -0.4294, -0.4459, -0.4414, -0.4359, -0.4387, -0.4455,
           -0.4544, -0.4641, -0.4793, -0.4984, -0.5001, -0.4996, -0.5065,
           -0.5022, -0.4981, -0.5046, -0.5136, -0.5243, -0.526 , -0.5254,
           -0.5239, -0.5131, -0.4995, -0.4831, -0.4595, -0.4368, -0.4089,
           -0.3758, -0.3415, -0.299 , -0.2475, -0.1889, -0.11  , -0.0254,
            0.0702,  0.1958,  0.3341,  0.4824,  0.6129,  0.709 ,  0.7898,
            0.8648,  0.9185,  0.955 ,  0.9721,  0.971 ,  0.9652,  0.9486,
            0.9194,  0.8885,  0.8618,  0.8491,  0.8154,  0.7478,  0.7017,
            0.6724,  0.6105,  0.5356,  0.4669,  0.3918,  0.323 ,  0.2405])
    wl_times = pd.date_range('2023-06-13 21:00','2023-06-15 02:00',freq='10min')

    # Then store it as a pandas dataframe and compute extremes
    wl_pd = pd.DataFrame({'values':wl_vals},index=wl_times)
    wl_pd_ext = hatyan.calc_HWLW(wl_pd)

    assert len(wl_pd_ext) == 3
    assert (wl_pd_ext['HWLWcode'] == [2,1,2]).all()


@pytest.mark.systemtest
def test_meas_HWLW_toomuch_19y():
    """
    Additional testcase for https://github.com/Deltares/hatyan/issues/85
    The min_width parameter for calc_HWLW() has to be approx 2hrs to properly compute and number almost all of the 19y timeseries
    Stations HELLVSS/KRIMPADLK/RAKND still fail, but it seems arbitrary and sometimes computing+numbering extremes per year does work
    """
    current_station = 'HOEKVHLD'
    file_dia = os.path.join(dir_testdata,f'{current_station}_obs19.txt')
    
    wl_pd = hatyan.read_dia(file_dia)
    
    wl_pd_ext = hatyan.calc_HWLW(wl_pd)
    
    # For testing: plot water level timeseries with peaks identified and labeled. 
    #hatyan.plot_timeseries(ts=wl_pd,ts_ext=wl_pd_ext)
    
    # RUN HATYAN: Assign numbers to the extremes
    wl_pd_ext = hatyan.calc_HWLWnumbering(wl_pd_ext,station=current_station)

    
@pytest.mark.systemtest
def test_calc_HWLW_ts_with_gap():
    """
    testcase that checks if no extreme is computed at the start/end of a gap
    like was documented in https://github.com/Deltares/hatyan/issues/97
    """
    
    current_station = 'VLISSGN'
    
    file_pred = os.path.join(dir_testdata,f'{current_station}_pre.txt')
    ts_prediction = hatyan.read_dia(filename=file_pred, station=current_station)
    ts_prediction = ts_prediction.loc[slice("2019-01","2019-01")]
    ts_prediction.loc[slice("2019-01-12 04:00", "2019-01-19 04:00")] = np.nan
    
    ts_ext_prediction = hatyan.calc_HWLW(ts=ts_prediction)
    
    # assert the number of extremes, was 91 before the fix of the issue
    assert len(ts_ext_prediction) == 90


@pytest.mark.systemtest
def test_calc_HWLW_unsorted():
    """
    testcase that checks if no extreme is computed at the start/end of a gap
    like was documented in https://github.com/Deltares/hatyan/issues/97
    """
    
    current_station = 'VLISSGN'
    
    file_pred = os.path.join(dir_testdata,f'{current_station}_pre.txt')
    ts_prediction = hatyan.read_dia(filename=file_pred, station=current_station)
    ts_prediction = ts_prediction.loc[slice("2019-01","2019-01")]
    ts_prediction = ts_prediction.sort_values("values")
    
    with pytest.raises(ValueError) as e:
        _ = hatyan.calc_HWLW(ts=ts_prediction)
    assert "timeseries is not monotonic increasing" in str(e.value)


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
    times_pred = slice(dt.datetime(2009,12,31,14),dt.datetime(2010,1,2,12), "1min")
    current_station = 'DENHDR'

    file_data_comp0 = os.path.join(dir_testdata,'%s_ana.txt'%(current_station))
    comp_merged = hatyan.read_components(filename=file_data_comp0)
    
    ts_prediction_HWLWno = hatyan.prediction(comp=comp_merged, times=times_pred)
    ts_ext_prediction_HWLWno_pre = hatyan.calc_HWLW(ts=ts_prediction_HWLWno)
    
    ts_ext_prediction_HWLWno = hatyan.calc_HWLWnumbering(ts_ext=ts_ext_prediction_HWLWno_pre, station=current_station)
    
    values_actual = ts_ext_prediction_HWLWno_pre['values'].values
    values_expected = np.array([-0.80339617,  0.61842454, -0.76465349,  0.79758517, -0.87872312])
    hwlwcodes_actual = ts_ext_prediction_HWLWno_pre['HWLWcode'].values
    hwlwcodes_expected = np.array([2, 1, 2, 1, 2])
    hwlwnos_actual = ts_ext_prediction_HWLWno['HWLWno'].values
    hwlwnos_expected = np.array([7057, 7058, 7058, 7059, 7059])
    
    assert len(ts_ext_prediction_HWLWno_pre) == 5
    assert np.allclose(values_actual, values_expected)
    assert np.allclose(hwlwcodes_actual, hwlwcodes_expected)
    assert np.allclose(hwlwnos_actual, hwlwnos_expected)


@pytest.mark.systemtest
@pytest.mark.parametrize("current_station, yr", [pytest.param(x, y, id='%s_%d'%(x,y)) for y in [2018,2022,2026] for x in ['HOEKVHLD','ROTTDM','DENHDR','LITHDP','EURPHVN']])
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
    current_station = 'LITHDP'
    yr=2000
    yr=2018
    """
    
    print(current_station)
    file_data_comp0 = os.path.join(dir_testdata,'%s_ana.txt'%(current_station))
    times_pred = slice(dt.datetime(yr,1,1),dt.datetime(yr+1,1,1), "5min")
    
    #component groups
    COMP_merged = hatyan.read_components(filename=file_data_comp0)
    
    #prediction and validation
    ts_prediction = hatyan.prediction(comp=COMP_merged, times=times_pred)
    ts_ext_prediction = hatyan.calc_HWLW(ts=ts_prediction)
    
    #calculate tidal wave number
    ts_ext_prediction_HWLWno = hatyan.calc_HWLWnumbering(ts_ext=ts_ext_prediction, station=current_station)
    print(current_station, yr)
    print('all HWLW values:')
    print(ts_ext_prediction_HWLWno)
    
    HW_data = ts_ext_prediction_HWLWno[ts_ext_prediction_HWLWno['HWLWcode']==1]
    HW_data_diff1bool = (HW_data['HWLWno'].diff().iloc[1:].values==1)
    print('%i parts of HW array containing gaps:'%((~HW_data_diff1bool).sum()))
    if not HW_data_diff1bool.all():
        ids_false = np.nonzero(~HW_data_diff1bool)[0]
        for id_false in ids_false: 
            print(HW_data.iloc[id_false-1:id_false+3])
    
    LW_data = ts_ext_prediction_HWLWno[ts_ext_prediction_HWLWno['HWLWcode']==2]
    LW_data_diff1bool = (LW_data['HWLWno'].diff().iloc[1:].values==1)
    print('%i parts of LW array containing gaps:'%((~LW_data_diff1bool).sum()))
    if not LW_data_diff1bool.all():
        ids_false = np.nonzero(~LW_data_diff1bool)[0]
        for id_false in ids_false: 
            print(LW_data.iloc[id_false-1:id_false+3])
    
    assert HW_data_diff1bool.all()
    assert LW_data_diff1bool.all()


@pytest.mark.parametrize("current_station", [pytest.param(x, id=x) for x in ['HOEKVHLD','DENHDR']])
@pytest.mark.systemtest
def test_frommergedcomp_HWLW_345(current_station):
    """
    This test produces a prediction for a period for HOEKVHLD, based on an imported component list. It then calculates extremes (HW/LW) with several settings and numbers them both (including 3/4/5 HWLWcodes around aggers, excluding 11/22 HWLWcodes). 
    #HOEKVHLD has alternating aggers, DENHDR has double HW's (which should not be numbered as aggers)
    
    current_station = 'HOEKVHLD'
    current_station = 'DENHDR'
    """
    
    file_data_comp0 = os.path.join(dir_testdata,'%s_ana.txt'%(current_station))
    comp_merged = hatyan.read_components(filename=file_data_comp0)
    
    times_pred = slice(dt.datetime(2010,1,31,3),dt.datetime(2010,2,17,12), "1min") #longer period with alternating aggers and no aggers, also eerste HW wordt als lokaal ipv primair HW gezien, also extra agger outside of 1stLW/agger/2ndLW sequence
    # times_pred = slice(dt.datetime(2010,1,31),dt.datetime(2010,2,3), "1min") #extra agger outside of 1stLW/agger/2ndLW sequence
    # times_pred = slice(dt.datetime(2010,6,26),dt.datetime(2010,6,28), "1min") #lokale extremen tussen twee laagste LWs
    # times_pred = slice(dt.datetime(2019,2,1),dt.datetime(2019,2,2), "1min") #eerste HW wordt als lokaal ipv primair HW gezien (lage prominence door dicht op begin tijdserie) >> warning
    
    ts_prediction_HWLWno = hatyan.prediction(comp=comp_merged, times=times_pred)
    ts_ext_prediction_main = hatyan.calc_HWLW(ts=ts_prediction_HWLWno)
    #ts_ext_prediction_all = hatyan.calc_HWLW(ts=ts_prediction_HWLWno, calc_HWLW345=True, calc_HWLW345_cleanup1122=False)
    ts_ext_prediction_clean = hatyan.calc_HWLW(ts=ts_prediction_HWLWno, calc_HWLW345=True)#, calc_HWLW345_cleanup1122=True) #for numbering, cannot cope with 11/22 HWLWcodes
    ts_ext_prediction_clean_tomain = hatyan.calc_HWLW12345to12(ts_ext_prediction_clean)
    
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
    
    assert np.allclose(ts_ext_prediction_main_HWLWno['HWLWno'].values, expected_ts_ext_prediction_main_HWLWno_HWLWno)
    assert np.allclose(ts_ext_prediction_main['HWLWcode'].values, expected_ts_ext_prediction_main_HWLWcode)
    assert np.allclose(ts_ext_prediction_clean_tomain['HWLWcode'].values, expected_ts_ext_prediction_main_HWLWcode)
    assert np.allclose(ts_ext_prediction_clean_HWLWno['HWLWcode'].values, expected_ts_ext_prediction_clean_HWLWno_HWLWcode)
    assert np.allclose(ts_ext_prediction_clean_HWLWno['HWLWno'].values, expected_ts_ext_prediction_clean_HWLWno_HWLWno)


@pytest.mark.systemtest
def test_hwlw_numbering_aggers():
    """
    This is probably already tested with other hwlw tests
    """
    
    times_pred = slice(dt.datetime(2010,1,31,3),dt.datetime(2010,2,17,12), "1min") #longer period with alternating aggers and no aggers, also eerste HW wordt als lokaal ipv primair HW gezien, also extra agger outside of 1stLW/agger/2ndLW sequence
    #HOEKVHLD has alternating aggers, DENHDR has double HW's
    expect_345_len = {'HOEKVHLD':99,'DENHDR':66}
    expect_345_code_sum = {'HOEKVHLD':267,'DENHDR':99}
    expect_345_nos_sum = {'HOEKVHLD':706061,'DENHDR':470712}
    
    for current_station in expect_345_len.keys():
    
        file_data_comp0 = os.path.join(dir_testdata,'%s_ana.txt'%(current_station))
        comp_merged = hatyan.read_components(filename=file_data_comp0)
        
        ts_prediction_HWLWno = hatyan.prediction(comp=comp_merged, times=times_pred)
        ts_ext_prediction_345 = hatyan.calc_HWLW(ts=ts_prediction_HWLWno, calc_HWLW345=True)
        ts_ext_prediction_345_HWLWno = hatyan.calc_HWLWnumbering(ts_ext=ts_ext_prediction_345, station=current_station)
        
        assert len(ts_ext_prediction_345_HWLWno) == expect_345_len[current_station]
        assert ts_ext_prediction_345_HWLWno["HWLWcode"].sum() == expect_345_code_sum[current_station]
        assert ts_ext_prediction_345_HWLWno["HWLWno"].sum() == expect_345_nos_sum[current_station]
        
        fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction_HWLWno, ts_ext=ts_ext_prediction_345)
        for irow, pdrow in ts_ext_prediction_345_HWLWno.iterrows():
            ax1.text(pdrow.name,pdrow['values'],pdrow['HWLWno'].astype(int))
        ax1.set_ylim(-1.2,1.7)
        hatyan.close()
        

@pytest.mark.systemtest
def test_19Ycomp4Ydia():
    # 1. define test data
    nodalfactors = True
    xfac = True
    const_list = hatyan.get_const_list_hatyan('year')
    current_station = 'VLISSGN'
    
    #comp0
    file_data_comp0 = os.path.join(dir_testdata,f'{current_station}_obs?.txt')
    ts_measurements_group0 = hatyan.read_dia(filename=file_data_comp0, station=current_station)
    comp_frommeasurements_avg_group0 = hatyan.analysis(ts=ts_measurements_group0, const_list=const_list, nodalfactors=nodalfactors, xfac=xfac, fu_alltimes=False, analysis_perperiod='Y')
    comp_frommeasurements_avg_group0.station = current_station
    
    #comp1
    file_data_comp1 = os.path.join(dir_testdata,'%s_ana.txt'%(current_station))
    comp_fromfile_group1 = hatyan.read_components(filename=file_data_comp1)
    
    #merge component groups (SA/SM from 19Y, rest from 4Y)
    COMP_merged = hatyan.merge_componentgroups(comp_main=comp_frommeasurements_avg_group0, comp_sec=comp_fromfile_group1.loc[['SA','SM']])
    
    #prediction and validation
    times_pred = slice(dt.datetime(2019,1,1),dt.datetime(2019,1,1,12), "10min")
    ts_prediction = hatyan.prediction(COMP_merged, times=times_pred)
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
    cs_comps = pd.DataFrame({'CS_comps_from':['K1','N2','2MN2','S2','S2'],
                             'CS_comps_derive':['P1','NU2','LABDA2','K2','T2'],
                             'CS_ampfacs':[0.33,0.22,0.48,0.29,0.05],
                             'CS_degincrs':[-11,-24,174,1,-24]})
    file_data_comp1 = os.path.join(dir_testdata,'%s_ana.txt'%(current_station))
    times_pred = slice(dt.datetime(2019,1,1),dt.datetime(2019,1,1,12), "10min")
    
    # 3. run test
    #component groups
    ts_measurements_group0 = hatyan.read_dia(filename=file_data_comp0, station=current_station)
    comp_frommeasurements_avg_group0 = hatyan.analysis(ts=ts_measurements_group0, nodalfactors=nodalfactors, xfac=xfac, 
                                                       const_list=const_list, analysis_perperiod=False, fu_alltimes=False, 
                                                       cs_comps=cs_comps)
    comp_fromfile_group1 = hatyan.read_components(filename=file_data_comp1)
    
    #merge component groups (SA/SM from 19Y, rest from 4Y)
    COMP_merged = hatyan.merge_componentgroups(comp_main=comp_frommeasurements_avg_group0, comp_sec=comp_fromfile_group1.loc[['SA','SM']])
    
    #prediction and validation
    ts_prediction = hatyan.prediction(COMP_merged, times=times_pred)
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
    times_pred = slice(dt.datetime(2019,1,1),dt.datetime(2019,1,1,12), "10min")

    #component groups
    ts_measurements_group0 = hatyan.read_dia(filename=file_data_comp0, station=current_station)
    comp_frommeasurements_avg_group0 = hatyan.analysis(ts=ts_measurements_group0, nodalfactors=nodalfactors, xfac=xfac, const_list=const_list, analysis_perperiod='Y', fu_alltimes=False)
    ts_measurements_group1 = hatyan.read_dia(filename=file_data_comp1, station=current_station)
    comp_frommeasurements_avg_group1 = hatyan.analysis(ts=ts_measurements_group1, nodalfactors=nodalfactors, xfac=xfac, const_list=const_list, analysis_perperiod=False, fu_alltimes=False)
    
    #merge component groups (SA/SM from 19Y, rest from 4Y)
    COMP_merged = hatyan.merge_componentgroups(comp_main=comp_frommeasurements_avg_group0, comp_sec=comp_frommeasurements_avg_group1.loc[['SA','SM']])
    
    #prediction and validation
    ts_prediction = hatyan.prediction(COMP_merged, times=times_pred)
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
    times_pred = slice(dt.datetime(2019,1,1),dt.datetime(2019,1,1,12), "10min")
    
    #component groups
    ts_measurements_group0 = hatyan.read_dia(filename=file_data_comp0, station=current_station)
    comp_frommeasurements_avg_group0 = hatyan.analysis(ts=ts_measurements_group0, nodalfactors=nodalfactors, xfac=xfac, const_list=const_list, analysis_perperiod='Y', fu_alltimes=False)
    
    #prediction and validation
    ts_prediction = hatyan.prediction(comp=comp_frommeasurements_avg_group0, times=times_pred)
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


