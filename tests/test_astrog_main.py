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
    dTLastExp = np.array([0.00044898, 0.00045681, 0.00046463, 0.00047245, 0.00048028,
                          0.0004881 , 0.00049593, 0.00050375, 0.00051157, 0.0005194 ,
                          0.00052722, 0.00053505, 0.00054287, 0.00055069, 0.00055852,
                          0.00056634, 0.00057417, 0.00058199, 0.00058981, 0.00059764,
                          0.00060546, 0.00061329, 0.00062111, 0.00062894, 0.00063676,
                          0.00064458, 0.00065241, 0.00066023, 0.00066806, 0.00067588,
                          0.0006837 , 0.00069153, 0.00069935, 0.00070718, 0.000715  ,
                          0.00072282, 0.00073065, 0.00073847, 0.0007463 , 0.00075412,
                          0.00076194, 0.00076977, 0.00077759, 0.00078542, 0.00079324,
                          0.00080106, 0.00080889, 0.00081671, 0.00082454, 0.00083236,
                          0.00084019])
    dTHistExp = np.array([0.00047442, 0.00048597, 0.00049752, 0.00050907, 0.00052063,
                          0.00053218, 0.00054373, 0.00055528, 0.00056683, 0.00057838,
                          0.00058993, 0.00060148, 0.00061303, 0.00062458, 0.00063613,
                          0.00064769, 0.00065924, 0.00067079, 0.00068234, 0.00069389,
                          0.00070544, 0.00071699, 0.00072854, 0.00074009, 0.00069502,
                          0.00070312, 0.00071123, 0.00071933, 0.00072743, 0.00073553,
                          0.00074363, 0.0007463 , 0.00075116, 0.00075602, 0.00076088,
                          0.00076574, 0.0007706 , 0.00077546, 0.00078032, 0.00078519,
                          0.00079005, 0.00079491, 0.00079977, 0.00078542, 0.00079324,
                          0.00080106, 0.00080889, 0.00081671, 0.00082454, 0.00083236,
                          0.00084019])
    dTExctExp = np.array([0.00048824, 0.00048824, 0.00051139, 0.00052296, 0.00053454,
                          0.00054611, 0.00055769, 0.00056926, 0.00058083, 0.00059241,
                          0.00060398, 0.00061556, 0.00062713, 0.0006387 , 0.0006387 ,
                          0.00065028, 0.00065028, 0.00065028, 0.00066185, 0.00066185,
                          0.00067343, 0.000685  , 0.00069657, 0.00070815, 0.00071972,
                          0.00071972, 0.0007313 , 0.00074287, 0.00074287, 0.00075444,
                          0.00075444, 0.00075444, 0.00075444, 0.00075444, 0.00075444,
                          0.00075444, 0.00076602, 0.00076602, 0.00076602, 0.00077759,
                          0.00077759, 0.00077759, 0.00078917, 0.00078917, 0.00078917,
                          0.00080074, 0.00080074, 0.00080074, 0.00080074, 0.00080074,
                          0.00080074])
    
    # 3. Run test
    dTLastOut = hatyan.dT(timeInput,mode_dT='last')
    dTHistOut = hatyan.dT(timeInput,mode_dT='historical')
    dTExctOut = hatyan.dT(timeInput,mode_dT='exact')
    
    # 4. Vefication
    assert type(dTLastOut) == np.ndarray
    assert type(dTHistOut) == np.ndarray
    assert type(dTExctOut) == np.ndarray
    assert (dTHistOut[0:-8] != dTLastOut[0:-8]).all()
    assert (abs(dTLastExp-dTLastOut) < 10E-9).all()
    assert (abs(dTHistExp-dTHistOut) < 10E-9).all()
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
    parOutput = hatyan.astrab(timeInput,hatyan.dT(timeInput,mode_dT='last'))
    
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
    timeExpect = ([pd.DatetimeIndex(['2007-12-31 18:08:06.634848001']),
                   pd.DatetimeIndex(['2008-01-01 06:28:32.856047999']),
                   pd.DatetimeIndex(['2007-12-17 10:18:26.216188801']),
                   pd.DatetimeIndex(['2007-12-24 01:16:42.818976000']),
                   pd.DatetimeIndex(['2007-12-31 07:51:54.693312000']),
                   pd.DatetimeIndex(['2008-01-08 11:38:04.174895999']),
                   pd.DatetimeIndex(['2008-01-01 00:54:57.698764800']),
                   pd.DatetimeIndex(['2008-01-01 11:07:34.708809600']),
                   pd.DatetimeIndex(['2008-01-01 07:48:35.142057600']),
                   pd.DatetimeIndex(['2008-01-01 07:48:35.142143999']),
                   pd.DatetimeIndex(['2008-03-20 05:49:10.054166400']),
                   pd.DatetimeIndex(['2007-06-21 18:07:03.793017600']),
                   pd.DatetimeIndex(['2007-09-23 09:52:02.142432000']),
                   pd.DatetimeIndex(['2007-12-22 06:08:46.055884800']),
                   pd.DatetimeIndex(['2008-01-03 08:22:45.903254400']),
                   pd.DatetimeIndex(['2008-01-03 08:22:45.696239999'])])
    
    # 3. Run test
    timeOutput = []
    for iMode in range(1,17):
        timeOutput.append(hatyan.astrac(timeInput,hatyan.dT(timeInput),np.array(iMode)))
    
    # 4. Vefication
    for iMode in range(0,16):
        assert type(timeOutput[iMode]) == pd.DatetimeIndex
        assert abs(timeExpect[iMode]-timeOutput[iMode]).total_seconds() < 10E-5
