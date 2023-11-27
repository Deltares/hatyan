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


@pytest.mark.systemtest
def test_astrog_dT():
    # 1. Input
    timeInput = pd.date_range(start=dt.datetime(1970,12,31,12,31),end=dt.datetime(2020,12,31,12,31),freq='Y')
    
    # 2. Expectations
    dTLastExp = np.array([38.792, 39.468, 40.144, 40.82 , 41.496, 42.172, 42.848, 43.524,
                           44.2  , 44.876, 45.552, 46.228, 46.904, 47.58 , 48.256, 48.932,
                           49.608, 50.284, 50.96 , 51.636, 52.312, 52.988, 53.664, 54.34 ,
                           55.016, 55.692, 56.368, 57.044, 57.72 , 58.396, 59.072, 59.748,
                           60.424, 61.1  , 61.776, 62.452, 63.128, 63.804, 64.48 , 65.156,
                           65.832, 66.508, 67.184, 67.86 , 68.536, 69.212, 69.888, 70.564,
                           71.24 , 71.916, 72.592])
    dTExctExp = np.array([42.184, 42.184, 43.184, 44.184, 45.184, 46.184, 47.184, 48.184,
                            49.184, 50.184, 51.184, 52.184, 53.184, 54.184, 54.184, 55.184,
                            55.184, 55.184, 56.184, 56.184, 57.184, 58.184, 59.184, 60.184,
                            61.184, 61.184, 62.184, 63.184, 63.184, 64.184, 64.184, 64.184,
                            64.184, 64.184, 64.184, 64.184, 65.184, 65.184, 65.184, 66.184,
                            66.184, 66.184, 67.184, 67.184, 67.184, 68.184, 68.184, 69.184,
                            69.184, 69.184, 69.184])
    
    # 3. Run test
    dTLastOut = hatyan.dT(timeInput,dT_fortran=True)
    dTExctOut = hatyan.dT(timeInput,dT_fortran=False)
    
    # 4. Vefication
    assert type(dTLastOut) == np.ndarray
    assert type(dTExctOut) == np.ndarray
    assert (abs(dTLastExp-dTLastOut) < 10E-9).all()
    assert (abs(dTExctExp-dTExctOut) < 10E-9).all()
    

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
    parOutput = hatyan.astrab(timeInput,hatyan.dT(timeInput,dT_fortran=True))
    
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
        timeOutput.append(hatyan.astrac(timeInput,dT_fortran=False,mode=np.array(iMode)))
    
    # 4. Vefication
    for iMode in range(1,17):
        print(iMode)
        assert type(timeOutput[iMode-1]) == pd.DatetimeIndex
        assert abs(timeExpect[iMode-1]-timeOutput[iMode-1]).total_seconds() < 10E-5


def test_astrog_astrac_inf_to_zero():
    
    # script settings
    tstart = dt.datetime(2003,1,1)
    tstop = dt.datetime(2003,6,1)
    
    # moonrise and -set
    moonriseset_python = hatyan.astrog_moonriseset(tFirst=tstart, tLast=tstop)
    
    # assert len, should be 292
    # but 0 RATE results in nan if we do not use posinf=0 in np.nan_to_num(), which gives len 291
    assert len(moonriseset_python) == 292


@pytest.mark.systemtest
def test_astrog_leapsecondslist():
    leap_seconds_pd, expirydate = hatyan.get_leapsecondslist_fromurlorfile()
    
