# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 11:17:44 2021

@author: laan_st

"""


import pytest
import numpy as np
import pandas as pd
import datetime as dt
import hatyan


@pytest.mark.unittest
def test_astrog_deprecated_tzone_argument():
    func_list = [
        hatyan.astrog.astrog_culminations,
        hatyan.astrog.astrog_phases,
        hatyan.astrog.astrog_sunriseset,
        hatyan.astrog.astrog_moonriseset,
        hatyan.astrog.astrog_anomalies,
        hatyan.astrog.astrog_seasons,
    ]
    for func in func_list:
        with pytest.raises(DeprecationWarning) as e:
            func(tFirst="2000-01-01", tLast="2001-01-01", tzone="UTC+01:00")
        assert "Argument 'tzone' has been deprecated for" in str(e.value)


@pytest.mark.unittest
def test_astrog_convert_str2datetime_str_tzdifferent():
    with pytest.raises(AssertionError):
        hatyan.astrog.convert_str2datetime("1887-01-01 00:00 +00:00",
                                           "2022-12-31 00:00 +01:00")


@pytest.mark.unittest
def test_astrog_convert_str2datetime_enddate_toosmall():
    with pytest.raises(ValueError) as e:
        hatyan.astrog.convert_str2datetime("2022-01-01 00:00 +00:00",
                                           "2010-12-31 00:00 +00:00")
    assert "start_date 2022-01-01 00:00:00 is larger than end_date 2010-12-31 00:00:00" in str(e.value)


@pytest.mark.unittest
def test_astrog_convert_str2datetime_str_naive():
    tstart, tstop, tzone = hatyan.astrog.convert_str2datetime("1887-01-01","2022-12-31")
    assert tstart == pd.Timestamp("1887-01-01")
    assert tstart.tz is None
    assert tstop == pd.Timestamp("2022-12-31")
    assert tstop.tz is None
    assert tzone == "UTC"


@pytest.mark.unittest
def test_astrog_convert_str2datetime_pd_naive():
    tstart, tstop, tzone = hatyan.astrog.convert_str2datetime(pd.Timestamp("1887-01-01"),
                                                              pd.Timestamp("2022-12-31"))
    assert tstart == pd.Timestamp("1887-01-01")
    assert tstart.tz is None
    assert tstop == pd.Timestamp("2022-12-31")
    assert tstop.tz is None
    assert tzone == "UTC"


@pytest.mark.unittest
def test_astrog_convert_str2datetime_str_met():
    tstart, tstop, tzone = hatyan.astrog.convert_str2datetime("1887-01-01 00:00 +01:00",
                                                              "2022-12-31 00:00 +01:00")
    assert tstart == pd.Timestamp("1886-12-31 23:00:00")
    assert tstart.tz is None
    assert tstop == pd.Timestamp("2022-12-30 23:00:00")
    assert tstop.tz is None
    assert tzone == dt.timezone(dt.timedelta(seconds=3600))


@pytest.mark.unittest
def test_astrog_convert_str2datetime_pd_met():
    tstart, tstop, tzone = hatyan.astrog.convert_str2datetime(pd.Timestamp("1887-01-01 00:00 +01:00"),
                                                              pd.Timestamp("2022-12-31 00:00 +01:00"))
    assert tstart == pd.Timestamp("1886-12-31 23:00:00")
    assert tstart.tz is None
    assert tstop == pd.Timestamp("2022-12-30 23:00:00")
    assert tstop.tz is None
    assert tzone == dt.timezone(dt.timedelta(seconds=3600))


@pytest.mark.unittest
def test_astrog_convert_str2datetime_pd_utc():
    tstart, tstop, tzone = hatyan.astrog.convert_str2datetime(pd.Timestamp("1887-01-01 00:00 +00:00"),
                                                              pd.Timestamp("2022-12-31 00:00 +00:00"))
    assert tstart == pd.Timestamp("1887-01-01")
    assert tstart.tz is None
    assert tstop == pd.Timestamp("2022-12-31")
    assert tstop.tz is None
    assert tzone == dt.timezone.utc
    

@pytest.mark.systemtest
def test_astrog_dT():
    # 1. Input
    timeInput = pd.date_range(start=dt.datetime(1970,12,31,12,31),end=dt.datetime(2020,12,31,12,31),freq='YS')
    
    # 2. Expectations
    dTLastExp = np.array([39.468, 40.144, 40.82 , 41.496, 42.172, 42.848, 43.524, 44.2  ,
           44.876, 45.552, 46.228, 46.904, 47.58 , 48.256, 48.932, 49.608,
           50.284, 50.96 , 51.636, 52.312, 52.988, 53.664, 54.34 , 55.016,
           55.692, 56.368, 57.044, 57.72 , 58.396, 59.072, 59.748, 60.424,
           61.1  , 61.776, 62.452, 63.128, 63.804, 64.48 , 65.156, 65.832,
           66.508, 67.184, 67.86 , 68.536, 69.212, 69.888, 70.564, 71.24 ,
           71.916, 72.592])
    dTExctExp = np.array([42.184, 42.184, 44.184, 45.184, 46.184, 47.184, 48.184, 49.184,
           50.184, 51.184, 51.184, 52.184, 53.184, 54.184, 54.184, 55.184,
           55.184, 56.184, 56.184, 57.184, 58.184, 58.184, 59.184, 60.184,
           61.184, 62.184, 62.184, 63.184, 64.184, 64.184, 64.184, 64.184,
           64.184, 64.184, 64.184, 65.184, 65.184, 65.184, 66.184, 66.184,
           66.184, 66.184, 67.184, 67.184, 67.184, 68.184, 69.184, 69.184,
           69.184, 69.184])
    
    # 3. Run test
    dTLastOut = hatyan.astrog.dT(timeInput,dT_fortran=True)
    dTExctOut = hatyan.astrog.dT(timeInput,dT_fortran=False)
    
    # 4. Vefication
    assert type(dTLastOut) == np.ndarray
    assert type(dTExctOut) == np.ndarray
    assert np.allclose(dTLastExp, dTLastOut)
    assert np.allclose(dTExctExp, dTExctOut)


@pytest.mark.systemtest
def test_astrog_astrab():
    
    # 1. Input
    timeInput = dt.datetime(2008,1,1,12,31)
    
    # 2. Expectations
    parExpect = {'EHMOON': np.array([448.08691901]),
                 'DECMOO': np.array([-13.02683311]),
                 'PARLAX': np.array([3262.092635]),
                 'DPAXDT': np.array([-18.05431966]),
                 'ALTMOO': np.array([-13.09124181]),
                 'ELONG':  np.array([283.08899401]),
                 'ALTSUN': np.array([14.08083711]),
                 'LONSUN': np.array([280.45809545]),
                 'EQELON': np.array([278.83322595]),
                 'DECSUN': np.array([-23.02828019]),
                 'DISSUN': np.array([0.98330485]),
                 'EHARI':  np.array([5.0317035]),
                 'RASUN':  np.array([281.37522918]),
                 'EHSUN':  np.array([6.92014497]),
                 'LONMOO': np.array([203.54708946]),
                 'LATMOO': np.array([-4.18226297]),
                 'RAMOON': np.array([200.20845513]),
                 'ANM':    np.array([151.15606465])}
    
    # 3. Run test
    parOutput = hatyan.astrog.astrab(timeInput,hatyan.astrog.dT(timeInput,dT_fortran=True))
    
    # 4. Vefication
    assert type(parOutput) == dict
    assert len(parOutput)  == 18
    for parName in parExpect:
        assert len(parOutput[parName]) == 1
        assert abs(parExpect[parName][0]-parOutput[parName][0]) < 10E-9
    
    
@pytest.mark.systemtest
def test_astrog_astrac():
    
    # 1. Input
    timeInput  = dt.datetime(2008,1,1,6,31)
    
    # 2. Expectations
    timeExpect = [pd.DatetimeIndex(['2007-12-31 18:08:06.634848001']),
                  pd.DatetimeIndex(['2008-01-01 06:28:32.856047999']),
                  pd.DatetimeIndex(['2007-12-17 10:18:26.216188801']),
                  pd.DatetimeIndex(['2007-12-24 01:16:42.818976000']),
                  pd.DatetimeIndex(['2007-12-31 07:51:54.693312000']),
                  pd.DatetimeIndex(['2008-01-08 11:38:04.174895999']),
                  pd.DatetimeIndex(['2008-01-01 00:54:56.645894400']),
                  pd.DatetimeIndex(['2008-01-01 11:07:33.697843200']),
                  pd.DatetimeIndex(['2008-01-01 07:48:34.139385600']),
                  pd.DatetimeIndex(['2008-01-01 07:48:34.139471999']),
                  pd.DatetimeIndex(['2008-03-20 05:49:10.054166400']),
                  pd.DatetimeIndex(['2007-06-21 18:07:03.793017600']),
                  pd.DatetimeIndex(['2007-09-23 09:52:02.142432000']),
                  pd.DatetimeIndex(['2007-12-22 06:08:46.055884800']),
                  pd.DatetimeIndex(['2008-01-03 08:22:45.903254400']),
                  pd.DatetimeIndex(['2008-01-03 08:22:45.696239999'])]
     
    # 3. Run test
    timeOutput = []
    for iMode in range(1,17):
        timeOutput.append(hatyan.astrog.astrac(timeInput,dT_fortran=False,mode=np.array(iMode)))
    
    # 4. Vefication
    for iMode in range(1,17):
        print(iMode)
        assert type(timeOutput[iMode-1]) == pd.DatetimeIndex
        assert abs(timeExpect[iMode-1]-timeOutput[iMode-1]).total_seconds() < 10E-5


@pytest.mark.systemtest
def test_astrog_leapsecondslist():
    leap_seconds_pd, expirydate = hatyan.astrog.get_leapsecondslist_fromurlorfile()


@pytest.mark.systemtest
def test_astrog_culminations():
    start_date = "2020-01-01"
    end_date = "2021-01-01"
    # moon culminations
    culminations_python = hatyan.astrog_culminations(tFirst=start_date, tLast=end_date)
    
    datetimes = culminations_python.index
    assert datetimes.tz == dt.timezone.utc
    assert datetimes[0] == pd.Timestamp('2020-01-01 04:44:14.732116514+0000')
    assert datetimes[-1] == pd.Timestamp('2020-12-31 13:16:46.464269856+0000')

    subset = culminations_python.iloc[:10]
    assert subset["type"].tolist() == [1, 2, 1, 2, 1, 2, 1, 2, 1, 2]
    expected_type_str = [
        'lowerculmination', 'upperculmination',
        'lowerculmination', 'upperculmination',
        'lowerculmination', 'upperculmination',
        'lowerculmination', 'upperculmination',
        'lowerculmination', 'upperculmination',
        ]
    assert subset["type_str"].tolist() == expected_type_str
    expected_parallax = np.array([
        0.90436235, 0.90347198, 0.90333631, 0.90397446, 0.90539423,
           0.90759186, 0.9105516 , 0.91424508, 0.91863045, 0.92365127])
    assert np.allclose(subset["parallax"].values, expected_parallax)
    expected_declination = np.array([
        -9.12903352, -6.86053844, -4.53030521, -2.15683332,  0.24204067,
        2.64853148,  5.04424777,  7.40950843,  9.72265897, 11.95942792])
    assert np.allclose(subset["declination"].values, expected_declination)


@pytest.mark.systemtest
def test_astrog_culminations_met():
    start_date = "2020-01-01 00:00:00 +01:00"
    end_date = "2021-01-01 00:00:00 +01:00"
    # moon culminations
    culminations_python = hatyan.astrog_culminations(tFirst=start_date, tLast=end_date)
    
    datetimes = culminations_python.index
    assert datetimes.tz == dt.timezone(dt.timedelta(seconds=3600))
    assert datetimes[0] == pd.Timestamp('2020-01-01 05:44:14.732115815+0100')
    assert datetimes[-1] == pd.Timestamp('2020-12-31 14:16:46.464269693+0100')


@pytest.mark.systemtest
def test_astrog_phases():
    start_date = "2020-01-01"
    end_date = "2021-01-01"
    # lunar phases
    phases_python = hatyan.astrog_phases(tFirst=start_date, tLast=end_date)
    
    datetimes = phases_python.index
    assert datetimes.tz == dt.timezone.utc
    assert datetimes[0] == pd.Timestamp('2020-01-03 04:45:20+0000')
    assert datetimes[-1] == pd.Timestamp('2020-12-30 03:28:07+0000')
    
    expected_type = [1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 
                     1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 
                     1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 
                     1, 2]
    assert phases_python["type"].tolist() == expected_type
    expected_type_str = ['FQ', 'FM', 'LQ', 'NM', 'FQ', 'FM', 'LQ', 'NM', 
                         'FQ', 'FM', 'LQ', 'NM', 'FQ', 'FM', 'LQ', 'NM', 
                         'FQ', 'FM', 'LQ', 'NM', 'FQ', 'FM', 'LQ', 'NM', 
                         'FQ', 'FM', 'LQ', 'NM', 'FQ', 'FM', 'LQ', 'NM', 
                         'FQ', 'FM', 'LQ', 'NM', 'FQ', 'FM', 'LQ', 'NM', 
                         'FQ', 'FM', 'LQ', 'NM', 'FQ', 'FM', 'LQ', 'NM', 
                         'FQ', 'FM']
    assert phases_python["type_str"].tolist() == expected_type_str


@pytest.mark.systemtest
def test_astrog_moonriseset():
    start_date = "2020-01-01"
    end_date = "2021-01-01"
    # moonrise and -set
    moonriseset_python = hatyan.astrog_moonriseset(tFirst=start_date, tLast=end_date)
    moonriseset_python_perday = hatyan.convert2perday(moonriseset_python)
    
    datetimes = moonriseset_python.index
    assert datetimes.tz == dt.timezone.utc
    assert datetimes[0] == pd.Timestamp('2020-01-01 11:15:43+0000')
    assert datetimes[-1] == pd.Timestamp('2020-12-31 16:46:52+0000')
    
    datetimes_perday = pd.to_datetime(moonriseset_python_perday.index)
    assert datetimes_perday.tz is None
    assert datetimes_perday[0] == pd.Timestamp('2020-01-01 00:00:00')
    assert datetimes_perday[-1] == pd.Timestamp('2020-12-31 00:00:00')
    
    subset = moonriseset_python.iloc[:10]
    assert subset["type"].tolist() == [1, 2, 1, 2, 1, 2, 1, 2, 1, 2]
    expected_type_str = [
        'moonrise', 'moonset', 'moonrise', 'moonset',
        'moonrise', 'moonset', 'moonrise', 'moonset',
        'moonrise', 'moonset']
    assert subset["type_str"].tolist() == expected_type_str
    
    subset_perday = moonriseset_python_perday.iloc[:10]
    assert subset_perday["moonrise"].isnull().sum() == 0
    assert subset_perday["moonset"].isnull().sum() == 1


@pytest.mark.systemtest
def test_astrog_moonriseset_posinf():
    """
    This part of code resulted in timesteps being dropped in newer pandas versions.
    This was due to 0 RATES in astract, resulting in inf values in 'addtime'.
    This is solved by supplying posinf=0 to np.nan_to_num().
    This test therefore checks if the length of the resulting dataframe is correct.
    """
    
    # script settings
    tstart = dt.datetime(2003,1,1)
    tstop = dt.datetime(2003,6,1)
    
    # moonrise and -set
    moonriseset_python = hatyan.astrog_moonriseset(tFirst=tstart, tLast=tstop)
    
    assert len(moonriseset_python) == 292


@pytest.mark.systemtest
def test_astrog_sunriseset():
    start_date = "2020-01-01"
    end_date = "2021-01-01"
    # sunrise and -set
    sunriseset_python = hatyan.astrog_sunriseset(tFirst=start_date, tLast=end_date)
    sunriseset_python_perday = hatyan.convert2perday(sunriseset_python)
    
    datetimes = sunriseset_python.index
    assert datetimes.tz == dt.timezone.utc
    assert datetimes[0] == pd.Timestamp('2020-01-01 07:47:25+0000')
    assert datetimes[-1] == pd.Timestamp('2020-12-31 15:36:02+0000')
    
    datetimes_perday = pd.to_datetime(sunriseset_python_perday.index)
    assert datetimes_perday.tz is None
    assert datetimes_perday[0] == pd.Timestamp('2020-01-01 00:00:00')
    assert datetimes_perday[-1] == pd.Timestamp('2020-12-31 00:00:00')

    subset = sunriseset_python.iloc[:10]
    assert subset["type"].tolist() == [1, 2, 1, 2, 1, 2, 1, 2, 1, 2]
    expected_type_str = [
        'sunrise', 'sunset', 'sunrise', 'sunset', 
        'sunrise', 'sunset', 'sunrise', 'sunset', 
        'sunrise', 'sunset']
    assert subset["type_str"].tolist() == expected_type_str
    
    subset_perday = sunriseset_python_perday.iloc[:10]
    assert subset_perday["sunrise"].isnull().sum() == 0
    assert subset_perday["sunset"].isnull().sum() == 0
    

@pytest.mark.systemtest
def test_astrog_anomalies():
    start_date = "2020-01-01"
    end_date = "2021-01-01"
    # lunar anomalies
    anomalies_python = hatyan.astrog_anomalies(tFirst=start_date, tLast=end_date)
    
    datetimes = anomalies_python.index
    assert datetimes.tz == dt.timezone.utc
    assert datetimes[0] == pd.Timestamp('2020-01-02 01:25:44+0000')
    assert datetimes[-1] == pd.Timestamp('2020-12-24 16:48:18+0000')
    
    expected_type = [2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 
                     1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2]
    expected_type_str = ['apogeum', 'perigeum', 'apogeum', 'perigeum', 'apogeum', 'perigeum', 
                         'apogeum', 'perigeum', 'apogeum', 'perigeum', 'apogeum', 'perigeum', 
                         'apogeum', 'perigeum', 'apogeum', 'perigeum', 'apogeum', 'perigeum', 
                         'apogeum', 'perigeum', 'apogeum', 'perigeum', 'apogeum', 'perigeum', 
                         'apogeum', 'perigeum', 'apogeum']
    assert anomalies_python['type'].tolist() == expected_type
    assert anomalies_python['type_str'].tolist() == expected_type_str


@pytest.mark.systemtest
def test_astrog_seasons():
    start_date = "2020-01-01"
    end_date = "2021-01-01"
    # astronomical seasons
    seasons_python = hatyan.astrog_seasons(tFirst=start_date, tLast=end_date)
    
    datetimes = seasons_python.index
    assert datetimes.tz == dt.timezone.utc
    assert datetimes[0] == pd.Timestamp('2020-03-20 03:49:46+0000')
    assert datetimes[-1] == pd.Timestamp('2020-12-21 10:02:39+0000')
    
    assert seasons_python['type'].tolist() == [1, 2, 3, 4]
    assert seasons_python['type_str'].tolist() == ['spring', 'summer', 'autumn', 'winter']


@pytest.mark.unittest
def test_plot_astrog_diff():
    start_date_utc = pd.Timestamp(2000, 1, 1, tz="UTC")
    end_date_utc = pd.Timestamp(2000, 4, 1, tz="UTC")
    
    culminations_python = hatyan.astrog_culminations(tFirst=start_date_utc, tLast=end_date_utc)
    culminations_python_naive = culminations_python.tz_localize(None).reset_index()
    
    hatyan.plot_astrog_diff(culminations_python, culminations_python_naive, typeLab=['lower','upper'], timeBand=[-.18,.18])
