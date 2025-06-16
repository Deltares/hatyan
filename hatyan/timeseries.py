# -*- coding: utf-8 -*-
"""
timeseries.py contains all definitions related to hatyan timeseries.
"""

import os
import io
import glob
import numpy as np
import pandas as pd
import datetime as dt
import pytz
import scipy.signal as ssig
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq
from netCDF4 import Dataset, date2num, stringtoarr
import logging
from pyproj import Transformer

from hatyan.foreman import get_foreman_v0_freq
from hatyan.schureman import get_schureman_freqs
from hatyan.hatyan_core import get_const_list_hatyan
from hatyan.metadata import (metadata_from_diablocks, metadata_add_to_obj, 
                             metadata_from_obj, metadata_compare, 
                             wns_from_metadata)

__all__ = ["get_diaxycoords",
           "read_dia",
           "read_noos",
           "write_dia",
           "write_noos",
           "write_netcdf",
           "plot_timeseries",
           "plot_HWLW_validatestats",
           "resample_timeseries",
           "crop_timeseries",
           "timeseries_fft",
           "Timeseries_Statistics",
           "calc_HWLW",
           "calc_HWLWnumbering",
           "calc_HWLW12345to12",
           ]

file_path = os.path.realpath(__file__)
logger = logging.getLogger(__name__)


def calc_HWLW(ts, calc_HWLW345=False, buffer_hr=6):
    """
    
    Calculates extremes (high and low waters) for the provided timeseries. 
    This definition uses scipy.signal.find_peaks() with arguments 'distance' and 'prominence'. 
    The minimal 'distance' between two high or low water peaks is based on the
    M2 period: 12.42/1.5=8.28 hours for HW and 12.42/1.7=7.30 hours for LW (larger because of aggers). 
    The prominence for local extremes is set to 0.01m, to filter out very minor dips in the timeseries.
    If there are two equal high or low water values, the first one is taken. 
    There are no main high/low waters calculated within 6 hours of the start/end of the
    timeseries (keyword buffer_hr), since these can be invalid.
    Since scipy.signal.find_peaks() warns about nan values, those are removed first.
    Nans/gaps can influence the results since find_peaks does not know about time
    registration. This is also tricky for input timeseries with varying time interval.
    
    Parameters
    ----------
    ts : pandas.DataFrame
        The DataFrame should contain a 'values' column and a pd.DatetimeIndex as index,
        it contains the timeseries with a tidal prediction or water level measurements.
    calc_HWLW345 : boolean, optional
        Whether to also calculate local extremes, first/second low waters and 'aggers'. 
        The default is False, in which case only extremes per tidal period are calculated.
        When first/second low waters and aggers are calculated, the local extremes
        around highwater (eg double highwaters and dips) are filtered out first.
    
    Returns
    -------
    data_pd_HWLW : pandas.DataFrame
        The DataFrame contains colums 'times', 'values' and 'HWLWcode', it contains the
        times, values and codes of the timeseries that are extremes.
        1 (high water) and 2 (low water). And if calc_HWLW345=True also
        3 (first low water), 4 (agger) and 5 (second low water).

    """
    
    if not ts.index.is_monotonic_increasing:
        #otherwise "ValueError: 'list' argument must have no negative elements"
        raise ValueError(
            'timeseries is not monotonic increasing, supply sorted timeseries '
            '(ts = ts.index.sort_index()',
            )
    
    #calculate the amount of steps in a M2 period, based on the most occurring timestep 
    M2_period_sec = get_schureman_freqs(['M2']).loc['M2','period [hr]']*3600
    
    ts_steps_sec_most = np.argmax(np.bincount(pd.Series(np.diff(ts.index)).dt.total_seconds().astype(int).values))
    if ts_steps_sec_most > 60:
        logger.warning(
            'the timestep of the series for which to calculate extremes/HWLW is'
            f'{ts_steps_sec_most/60:.2f} minutes, but 1 minute is recommended',
            )
    elif ts_steps_sec_most == 0:
        raise ValueError('ts_steps_sec_most=0, check rounding issue')
    M2period_numsteps = M2_period_sec/ts_steps_sec_most
    # minimal width of 2 hours makes peakfinding more robust: https://github.com/Deltares/hatyan/issues/85
    minwidth_numsteps = 2*3600/ts_steps_sec_most

    data_pd_HWLW = pd.DataFrame({'times':ts.index,'values':ts['values'],'HWLWcode':np.nan}).reset_index(drop=True)
    #create empty HWLW dataframe
    if data_pd_HWLW['values'].isnull().any():
        data_pd_HWLW = data_pd_HWLW[~data_pd_HWLW['values'].isnull()].reset_index(drop=True)
        len_diff = len(ts)-len(data_pd_HWLW)
        len_diff_pct = (len(ts)-len(data_pd_HWLW))/len(ts)*100
        logger.warning(
            'the provided ts for extreme/HWLW calculation contained NaN values. To '
            f'avoid unexpected results from scipy.signal.find_peaks(), the {len_diff} '
            f'NaN values ({len_diff_pct:.2f}%) were removed from the ts before '
            'calculating extremes/HWLW.'
            )

    # minimal prominence to exclude very minor dips/peaks from being seen as aggers.
    min_prominence = 0.01
    
    if calc_HWLW345:
        # get all local extremes, including aggers and second high waters (1/2/11/22)
        # takes first value of two equal peaks
        LWid_all, LWid_all_properties = ssig.find_peaks(-data_pd_HWLW['values'].values, prominence=(min_prominence,None), width=(None,None), distance=None)
        HWid_all, HWid_all_properties = ssig.find_peaks(data_pd_HWLW['values'].values, prominence=(min_prominence,None), width=(None,None), distance=None)
        data_pd_HWLW.loc[LWid_all,'HWLWcode'] = 22 #all LW
        data_pd_HWLW.loc[HWid_all,'HWLWcode'] = 11 #all HW

    #get HWLW (extremes per tidal period).
    # LW: most stations work with factor 1.4
    # 1.5 results in all LW values for HoekvanHolland for 2000 (and for 2018/2022)
    # 1.7 results in all LW values for Rotterdam for 2000 (and for 2018/2022)
    # 1.8 results in all LW values for EURPHVN for 2026
    LWid_main_raw,LWid_main_properties = ssig.find_peaks(
        -data_pd_HWLW['values'].values,
        prominence=(min_prominence,None),
        width=(minwidth_numsteps,None),
        distance=M2period_numsteps/1.8,
        )
    # HW: most stations work with factor 1.4.
    # 1.5 results in all HW values for DenHelder for year 2000 (also for 1999-2002)
    # 1.7 results in all HW values for LITHDP 2018 (but it fails with 1.8)
    # 1.8 results in all HW values for LITHDP 2026
    # 1.9 results in all HW values for LITHDP 2022
    HWid_main_raw,HWid_main_properties = ssig.find_peaks(
        data_pd_HWLW['values'].values,
        prominence=(min_prominence,None),
        width=(minwidth_numsteps,None),
        distance=M2period_numsteps/1.9,
        )
    # remove main extremes within 6 hours of start/end of timeseries, since they are
    # often missed or invalid.
    buffer_dt = dt.timedelta(hours=buffer_hr)
    bool_larger = (data_pd_HWLW['times'] >= data_pd_HWLW['times'].iloc[0] + buffer_dt)
    bool_smaller = (data_pd_HWLW['times'] <= data_pd_HWLW['times'].iloc[-1] - buffer_dt)
    validtimes_idx = data_pd_HWLW.loc[bool_larger & bool_smaller].index
    LWid_main = LWid_main_raw[np.isin(LWid_main_raw,validtimes_idx)]
    HWid_main = HWid_main_raw[np.isin(HWid_main_raw,validtimes_idx)]
    #use valid values to continue process
    data_pd_HWLW.loc[LWid_main,'HWLWcode'] = 2
    data_pd_HWLW.loc[HWid_main,'HWLWcode'] = 1
    
    #drop all non-(local)extreme timesteps and convert HWLWcode column to integers
    data_pd_HWLW = data_pd_HWLW.dropna(subset=['HWLWcode'])
    data_pd_HWLW['HWLWcode'] = data_pd_HWLW['HWLWcode'].astype(int)

    # extensive debugging statistics
    prop_list = ['prominences','widths']
    data_pd_HWLW.loc[data_pd_HWLW['HWLWcode']==2,prop_list] = pd.DataFrame(LWid_main_properties,index=LWid_main_raw).loc[LWid_main,prop_list]
    logger.debug('LW values:\n%s\n'%(data_pd_HWLW[data_pd_HWLW['HWLWcode']==2]))
    data_pd_HWLW.loc[data_pd_HWLW['HWLWcode']==1,prop_list] = pd.DataFrame(HWid_main_properties,index=HWid_main_raw).loc[HWid_main,prop_list]
    logger.debug('HW values:\n%s\n'%(data_pd_HWLW[data_pd_HWLW['HWLWcode']==1]))
    if calc_HWLW345:
        LW_local_bool = ~np.isin(LWid_all, LWid_main)
        data_pd_HWLW.loc[data_pd_HWLW['HWLWcode']==22,prop_list] = pd.DataFrame(LWid_all_properties,index=LWid_all).loc[LW_local_bool,prop_list]
        logger.debug('LW_local values:\n%s\n'%(data_pd_HWLW[data_pd_HWLW['HWLWcode']==22]))
        HW_local_bool = ~np.isin(HWid_all, HWid_main)
        data_pd_HWLW.loc[data_pd_HWLW['HWLWcode']==11,prop_list] = pd.DataFrame(HWid_all_properties,index=HWid_all).loc[HW_local_bool,prop_list]
        logger.debug('HW_local values:\n%s\n'%(data_pd_HWLW[data_pd_HWLW['HWLWcode']==11]))
    
    if calc_HWLW345: #recalculate local LW/HWs between two main HWs to firstLW/agger/secondLW
        data_pd_HWLW = calc_HWLWlocalto345(data_pd_HWLW,HWid_main)
    
    #return to normal time-index
    data_pd_HWLW = data_pd_HWLW.set_index('times')
    
    # add metadata #TODO: update metadata for extremes?
    metadata = metadata_from_obj(ts)
    data_pd_HWLW = metadata_add_to_obj(data_pd_HWLW,metadata)
    
    return data_pd_HWLW


def calc_HWLWlocalto345(data_pd_HWLW,HWid_main):
    """
    Recalculate local LW/HWs between two main HWs to firstLW/agger/secondLW

    Parameters
    ----------
    data_pd_HWLW : TYPE
        DESCRIPTION.
    HWid_main : TYPE
        DESCRIPTION.

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    data_pd_HWLW : TYPE
        DESCRIPTION.

    """
    
    logger.info('calculating 1stLW/agger/2ndLW for all tidalperiods (between two HW values)...')
    for iTide, dummy in enumerate(HWid_main[:-1]):
        data_pd_HWLW_1tide = data_pd_HWLW.loc[HWid_main[iTide]:HWid_main[iTide+1],:]
        
        #filter local extremes around HW (only interested in aggers, so LW), this is necessary for eg DENHDR and PETTZD, otherwise second HW is seen as first LW
        data_pd_HWLW_1tide_minHW = data_pd_HWLW_1tide.loc[data_pd_HWLW_1tide['HWLWcode']==1,'values'].min()
        data_pd_HWLW_1tide_min = data_pd_HWLW_1tide['values'].min()
        data_pd_HWLW_1tide_mid = np.mean([data_pd_HWLW_1tide_minHW,data_pd_HWLW_1tide_min])
        bool_LWs = data_pd_HWLW_1tide['values']<data_pd_HWLW_1tide_mid
        data_pd_HWLW_1tide_noHWs = data_pd_HWLW_1tide[bool_LWs]
        
        if len(data_pd_HWLW_1tide_noHWs) > 3: #(attempt to) reduce to three values between two HWs
            logger.warning('more than 3 values between HWs, removing part of them')
            agger35_prim = data_pd_HWLW_1tide_noHWs[data_pd_HWLW_1tide_noHWs['HWLWcode']==2]
            if len(agger35_prim)>1:
                raise Exception('should be only one HWLWcode=2 per tide period')
            agger35_prim_loc = agger35_prim.index[0]
            agger35_sec_loc = data_pd_HWLW_1tide_noHWs.loc[data_pd_HWLW_1tide_noHWs['HWLWcode']==22,'values'].idxmin()
            agger35_loc = np.sort([agger35_prim_loc,agger35_sec_loc])
            data_pd_HWLW_1tide_noHWs = data_pd_HWLW_1tide_noHWs.loc[agger35_loc.min():agger35_loc.max(),:]
            agger4_loc = data_pd_HWLW_1tide_noHWs['values'].idxmax()
            data_pd_HWLW_1tide_noHWs = data_pd_HWLW_1tide_noHWs.loc[[agger35_loc.min(),agger4_loc,agger35_loc.max()],:]
        
        if len(data_pd_HWLW_1tide_noHWs) == 1: #primary low water already has code 2
            if data_pd_HWLW_1tide_noHWs['HWLWcode'].iloc[0] != 2:
                raise Exception('Only 1 LW value but does not have HWLWcode 2')
        elif len(data_pd_HWLW_1tide_noHWs) == 3:
            if not data_pd_HWLW_1tide_noHWs['values'].argmax() == 1:
                raise Exception('3 values between two HW values, but center one is not the largest:\n%s'%(data_pd_HWLW_1tide))
            agger345_loc = data_pd_HWLW_1tide_noHWs.index
            if not (data_pd_HWLW.loc[agger345_loc[0],'HWLWcode'] in [2,22] and data_pd_HWLW.loc[agger345_loc[1],'HWLWcode'] in [11] and data_pd_HWLW.loc[agger345_loc[2],'HWLWcode'] in [2,22]):
                raise Exception('3 values between two HW values, but do not correspond to LW/agger/LW:\n%s'%(data_pd_HWLW_1tide))
            data_pd_HWLW.loc[agger345_loc,'HWLWcode'] = [3,4,5]
        elif len(data_pd_HWLW_1tide_noHWs) == 2:
            logger.warning('2 values left between two HWs (slightly unexpected):\n%s'%(data_pd_HWLW_1tide))
        else:
            raise Exception('unexpected number of values between two HWs (0 or more than 3):\n%s'%(data_pd_HWLW_1tide))
    
    #remove remaining 11 and 22 values from array
    data_pd_HWLW = data_pd_HWLW.drop(data_pd_HWLW[data_pd_HWLW['HWLWcode']==11].index)
    data_pd_HWLW = data_pd_HWLW.drop(data_pd_HWLW[data_pd_HWLW['HWLWcode']==22].index)
    
    logger.info('finished calculating 1stLW/agger/2ndLW for all tidalperiods')
    
    return data_pd_HWLW


def calc_HWLW12345to12(data_HWLW_12345):
    """
    
    
    Parameters
    ----------
    data_HWLW12345 : TYPE
        DESCRIPTION.
    
    Returns
    -------
    None.
    
    """
    if len(data_HWLW_12345['HWLWcode'].unique()) <= 2:
        logger.info('skipping HWLW 12345 to 12 correction since no 345 found, returning input df')
        return data_HWLW_12345.copy()
    
    logger.info('starting HWLW 12345 to 12 correction')
    times_LWmin = []
    data_HW1 = data_HWLW_12345.loc[data_HWLW_12345['HWLWcode']==1]
    # computing minimum waterlevels after each HW. Using hardcoded 12hour period 
    # instead of from one HW to next HW since then we can also assess last LW values
    for timeHW in data_HW1.index: #np.arange(0,len(data_HW1)-1):
        if timeHW==data_HWLW_12345.index[-1]: #if last HW is last time of input dataframe
            continue
        tide_afterHW = data_HWLW_12345.loc[timeHW:timeHW+dt.timedelta(hours=12)]
        # remove first HW to avoid issues if LW is higher than HW due to surge
        tide_afterHW = tide_afterHW.iloc[1:]
        if len(tide_afterHW)==0:
            # this happens if there is no LW defined between two HWs
            # for instance after SCHEVNGN HW at '1948-04-30 19:50:00'
            continue
        time_minimum = tide_afterHW['values'].idxmin()
        times_LWmin.append(time_minimum)
    data_LW2 = data_HWLW_12345.loc[times_LWmin]
    data_LW2['HWLWcode'] = 2
    
    #optionally also concat first value if HWLWcode=2 (this is possible since we know above code starts on first HWLWcode=1)
    firstlast = []
    if data_HWLW_12345.iloc[0]['HWLWcode']==2:
        firstlast.append(data_HWLW_12345.iloc[[0]])
    
    data_HWLW_12 = pd.concat(firstlast+[data_HW1,data_LW2]).sort_index()
    
    return data_HWLW_12


def filter_duplicate_hwlwnos(ts_ext):
    # find rows that have a duplicated HWLW code+number
    bool_hwno_duplicated = ts_ext[['HWLWcode','HWLWno']].duplicated(keep=False)
    # return this filtered dataframe
    ts_ext_hwlwno_duplicated = ts_ext.loc[bool_hwno_duplicated,
                                          ['values','HWLWcode','HWLWno']]
    return ts_ext_hwlwno_duplicated


def calc_HWLWnumbering(ts_ext, station=None):
    """
    For calculation of the extremes numbering, w.r.t. the first high water at Cadzand in
    2000 (occurred on 1-1-2000 at approximately 9:45). 
    The number of every high and low water is calculated by taking the time difference
    between itself and the first high water at Cadzand, correcting it with the station
    phase difference (M2phasediff). 
    Low waters are searched for half an M2 period from the high waters. 
    By adding a search window of half the period of M2 (searchwindow_hr), even strong
    time variance between consecutive high or low waters should be caputered. 
    
    Parameters
    ----------
    ts_ext : pandas.DataFrame
        The DataFrame should contain a 'values' and 'HWLWcode' column and a
        pd.DatetimeIndex as index, it contains the times, values and codes of the
        timeseries that are extremes.
    station: string, optional
        The station for which the M2 phase difference should be retrieved from
        data_M2phasediff_perstation.txt. This value is the phase difference in degrees
        of the occurrence of the high water generated by the same tidal wave as the
        first high water in 2000 at Cadzand (actually difference between M2 phases
        of stations). This value is used to correct the search window of high/low water
        numbering. The default is None.
        Providing a value will result in a proper HWLWno, corresponing to CADZD.
        Providing None will result in a HWLWno that is a multiple of
        360degrees/M2_period_hr off (positive or negative). This is only an issue when
        comparing different stations, not comparing e.g. measured and predicted HW
        values of one station.
    
    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    ts_ext : pandas.DataFrame
        The input DataFrame with the column 'HWLWno' added, which contains the numbers of the extremes.

    """
        
    M2_period_hr = get_schureman_freqs(['M2']).loc['M2','period [hr]']
    firstHWcadz_fixed = pd.Timestamp("2000-01-01 09:45:00 +01:00")
    firstHWcadz_fixed = firstHWcadz_fixed.tz_convert(ts_ext.index.tz)
    searchwindow_hr = M2_period_hr/2
    
    if len(ts_ext) == 0:
        raise Exception('length of provided ts_ext is zero')
    
    if not all((ts_ext['HWLWcode']==1) | (ts_ext['HWLWcode']==2) | (ts_ext['HWLWcode']==3) | (ts_ext['HWLWcode']==4) | (ts_ext['HWLWcode']==5)):
        raise ValueError(
            'calc_HWLWnumbering() not implemented for HWLWcode other than 1,2,3,4,5 '
            '(so no HWLWcode 11 or 22 supported), provide extreme timeseries derived '
            'with Timeseries.calc_HWLW(calc_HWLW345=False) or '
            'Timeseries.calc_HWLW(calc_HWLW345=True, calc_HWLW345_cleanup1122=True)'
            )
    
    # copy object but retain metadata
    metadata_ext = metadata_from_obj(ts_ext)
    ts_ext = ts_ext.copy()
    ts_ext = metadata_add_to_obj(ts_ext, metadata_ext) # otherwise M2-analysis fails
    
    HW_bool = ts_ext['HWLWcode']==1
    HW_tdiff_cadzdraw = (ts_ext.loc[HW_bool].index.to_series()-firstHWcadz_fixed).dt.total_seconds()/3600
    if station is None:
        # TODO: local import since Importerror: cannot import name 'analysis' from
        # partially initialized module 'hatyan.analysis_prediction' (most likely due to a circular import) 
        from hatyan.analysis_prediction import analysis
        M2phase_cadzd = 48.81 # from analyse waterlevels CADZD over 2009 t/m 2012
        # high max_matrix_condition value necessary for some english stations. Not a big
        # issue since it should provide a phasediff also in case of only one HW+LW.
        comp_M2 = analysis(ts_ext,const_list=['M2'],max_matrix_condition=250)
        logger.debug(f"M2phase: {comp_M2.loc['M2','phi_deg']}, comp_M2_cadzd: {M2phase_cadzd}")
        M2phasediff_deg = (comp_M2.loc['M2','phi_deg'] - M2phase_cadzd+90)%360-90
        message_prefix = 'no value or None for argument M2phasediff provided, automatically calculated correction w.r.t. Cadzand based on M2phase:'
    else:
        file_M2phasediff = os.path.join(os.path.dirname(file_path),'data','data_M2phasediff_perstation.txt')
        stations_M2phasediff = pd.read_csv(file_M2phasediff, names=['M2phasediff'], comment='#', sep="\\s+")
        if station not in stations_M2phasediff.index:
            raise Exception(f'ERROR: station "{station}" not in file_M2phasediff ({file_M2phasediff})')
        M2phasediff_deg = stations_M2phasediff.loc[station,'M2phasediff']
        message_prefix = 'M2phasediff retrieved from file, correction w.r.t. Cadzand:'
    M2phasediff_hr = M2phasediff_deg/360*M2_period_hr
    logger.info(f'{message_prefix} {M2phasediff_hr:.2f} hours ({M2phasediff_deg:.2f} degrees)')
    HW_tdiff_cadzd = HW_tdiff_cadzdraw - M2phasediff_hr + searchwindow_hr
    HW_tdiff_div, HW_tdiff_mod_searchwindow = np.divmod(HW_tdiff_cadzd.values, M2_period_hr)
    HW_tdiff_mod = HW_tdiff_mod_searchwindow - searchwindow_hr
    if not all(np.abs(HW_tdiff_mod)<searchwindow_hr):
        raise ValueError('tidal wave numbering: not all HW fall into hardcoded search window')
    ts_ext.loc[HW_bool,'HWLWno'] = HW_tdiff_div
    
    for LWcode_2345 in [2,3,4,5]:
        LW_bool = ts_ext['HWLWcode']==LWcode_2345
        LW_tdiff_cadzdraw = (ts_ext.loc[LW_bool].index.to_series()-firstHWcadz_fixed).dt.total_seconds()/3600
        LW_tdiff_cadzd = LW_tdiff_cadzdraw - M2phasediff_hr + searchwindow_hr - M2_period_hr/2
        LW_tdiff_div, LW_tdiff_mod_searchwindow = np.divmod(LW_tdiff_cadzd.values, M2_period_hr)
        LW_tdiff_mod = LW_tdiff_mod_searchwindow - searchwindow_hr
        if not all(np.abs(LW_tdiff_mod)<searchwindow_hr):
            raise ValueError('tidal wave numbering: not all LW fall into hardcoded search window')
        ts_ext.loc[LW_bool,'HWLWno'] = LW_tdiff_div
    
    # check if HWLW code+numbers are duplicated
    ts_ext_hwlwno_duplicated = filter_duplicate_hwlwnos(ts_ext)
    if not ts_ext_hwlwno_duplicated.empty:
        raise ValueError('tidal wave numbering: HWLW code+numbers not always unique:\n'
                         f'{ts_ext_hwlwno_duplicated}'
                         )

    ts_ext['HWLWno'] = ts_ext['HWLWno'].astype(int)
    
    return ts_ext


def timeseries_fft(ts_residue, min_prominence=10**3, max_freqdiff=None, plot_fft=True, source='schureman'):
    
    logger.info('analyzing timeseries with fft and fftfreq')
    
    if ts_residue['values'].isnull().sum() > 0:
        raise ValueError('supplied timeseries contains nan values, use pd.interpolate first (dropping them will result in non-constant timestep which is also not possible for fft)')
    y = ts_residue['values'].values
    N = len(y)
    T = np.unique((ts_residue.index[1:]-ts_residue.index[:-1])).astype(float)/1e9/3600 #timestep in hours.
    if len(T)!=1:
        raise ValueError('timestep of supplied timeseries should be constant for fourier analysis')
    yf = fft(y)
    power = np.abs(yf)
    freq = fftfreq(N, T[0])
    peaks, peaks_properties = ssig.find_peaks(power[freq >=0], prominence=min_prominence)
    peak_freq =  freq[peaks]
    peak_power = power[peaks]
    
    if plot_fft:
        fig,ax = plt.subplots()
        ax.plot(freq[:N//2], power[:N//2])
        ax.plot(peak_freq, peak_power, 'ro')
        ax.grid()
        ax.set_xlim(0,0.5)
    
    if source=='schureman':
        const_list_all = get_const_list_hatyan(listtype='all_schureman')
        hatyan_freqs = get_schureman_freqs(const_list=const_list_all)
    elif source=='foreman':
        const_list_all = get_const_list_hatyan(listtype='all_foreman')
        dummy, hatyan_freqs = get_foreman_v0_freq(const_list=const_list_all)
        hatyan_freqs['period [hr]'] = 1/hatyan_freqs['freq']
    
    const_closest = []
    for peak_freq_one in peak_freq:
        hatyan_freqs_closest = hatyan_freqs.iloc[np.argmin(np.abs(hatyan_freqs['freq']-peak_freq_one)),:] #TODO: freq is not always close enough but is still added to list
        const_closest.append(hatyan_freqs_closest.name)
    hatyan_freqs_suggestions = hatyan_freqs.loc[const_closest,['freq','period [hr]']]
    hatyan_freqs_suggestions['peak_freq'] = peak_freq
    hatyan_freqs_suggestions['peak_freqdiff'] = (hatyan_freqs_suggestions['freq'] - hatyan_freqs_suggestions['peak_freq']).abs()
    hatyan_freqs_suggestions['peak_prominences'] = peaks_properties['prominences']
    if max_freqdiff is not None:
        #select below freqdiff treshold
        hatyan_freqs_suggestions = hatyan_freqs_suggestions.loc[hatyan_freqs_suggestions['peak_freqdiff']<max_freqdiff]
    logger.info('suggested constituents+freqs from hatyan:\n%s'%(hatyan_freqs_suggestions))
    
    return hatyan_freqs_suggestions


def plot_timeseries(ts, ts_validation=None, ts_ext=None, ts_ext_validation=None):
    """
    Creates a plot with the provided timeseries

    Parameters
    ----------
    ts : pandas.DataFrame
        The DataFrame should contain a 'values' column and a pd.DatetimeIndex as index,
        it contains the timeseries.
    ts_validation : pandas.DataFrame, optional
        The DataFrame should contain a 'values' column and a pd.DatetimeIndex as index,
        it contains the timeseries. The default is None.
    ts_ext : pandas.DataFrame, optional
        The DataFrame should contain a 'values' and 'HWLW_code' column and a
        pd.DatetimeIndex as index, it contains the times, values and codes of the
        timeseries that are extremes. The default is None.
    ts_ext_validation : pandas.DataFrame, optional
        The DataFrame should contain a 'values' and 'HWLW_code' column and a
        pd.DatetimeIndex as index, it contains the times, values and codes of the
        timeseries that are extremes. The default is None.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The generated figure handle, with which the figure can be adapted and saved.
    axs : (tuple of) matplotlib.axes._subplots.AxesSubplot
        The generated axis handle, whith which the figure can be adapted.

    """
        
    size_figure = (15,9)
    size_line_ts = 0.7
    size_marker_ts = 1
    figure_ylim_ts = [-3,3]
    figure_ylim_tsdiff = [-0.02,0.02]
    
    # get unit, assuming all other dataframes have the same unit
    eenheid = ts.attrs.get('eenheid', '-')
    
    if ts_validation is not None:
        times_predval_ext = [
            min(min(ts_validation.index),min(ts.index)),
            max(max(ts_validation.index),max(ts.index)),
            ]
    else:
        times_predval_ext = [min(ts.index), max(ts.index)]    

    fig, (ax1, ax2) = plt.subplots(2,1,figsize=size_figure, sharex=True, gridspec_kw={'height_ratios':[2,1]})
    
    ax1.set_title('hatyan timeseries')
    ax1.plot(ts.index, ts['values'],'o-',linewidth=size_line_ts,markersize=size_marker_ts, label='ts')
    if ts_validation is not None:
        # compute difference between timeseries for overlapping period
        ts_diff = ts['values'] - ts_validation['values']
        ts_diff = ts_diff.dropna()
        ax1.plot(ts_validation.index, ts_validation['values'],'o-',linewidth=size_line_ts,markersize=size_marker_ts, label='ts_validation', alpha=0.7)
        ax1.plot(ts_diff.index, ts_diff,'go-',linewidth=size_line_ts,markersize=size_marker_ts, label='difference', alpha=0.7)
    ax1.plot(times_predval_ext,[0,0],'-k',linewidth=size_line_ts)
    ts_mean = np.mean(ts['values'])
    ax1.plot(ts.index[[0,-1]],[ts_mean,ts_mean],'-r',linewidth=size_line_ts,label='mean of ts')
    if ts_ext is not None:
        HWLW_codesnames = {1:'HW (1)',
                           2:'LW (2)',
                           3:'LW1 (3)',
                           4:'topagger (4)',
                           5:'LW2 (5)',
                           11:'HW_local (11)',
                           22:'LW_local (22)'}
        for HWLW_code in HWLW_codesnames.keys():
            iExt = ts_ext['HWLWcode']==HWLW_code
            if iExt.any():
                HWLW_name = HWLW_codesnames[HWLW_code]
                HWLW_markersize=10
                if HWLW_code in [4,11,22]:
                    HWLW_markersize=5
                ax1.plot(ts_ext.index[iExt],ts_ext['values'][iExt],'x',markersize=HWLW_markersize,label=HWLW_name)
    if ts_ext_validation is not None:
        vali_codes = [1,2,3,4,5]
        vali_codenames = ['vali_HW','vali_LW','vali_LW1','vali_topagger','vali_LW2']
        for vali_code, vali_codename in zip(vali_codes,vali_codenames):
            vali_code_ids = ts_ext_validation['HWLWcode'].values==vali_code
            if any(vali_code_ids): #only plot vali_code in legend if present in HWLW_timeseries
                ax1.plot(ts_ext_validation.index[vali_code_ids],ts_ext_validation['values'][vali_code_ids],'1',markersize=10,label=vali_codename)
    ax1.set_ylim(figure_ylim_ts)
    ax2.set_xlabel('Time')
    ax1.set_ylabel(f'waterlevel [{eenheid}]')
    ax1.legend(loc='lower right')
    ax1.grid()
    if ts_validation is not None:
        ax2.plot(ts_diff.index, ts_diff, 'go-',linewidth=size_line_ts,markersize=size_marker_ts, label='difference')
        ax2.legend(loc='lower right') # create legend only if diff is plotted
    ax2.plot(times_predval_ext,[0,0],'-k',linewidth=size_line_ts)
    ax2.set_ylim(figure_ylim_tsdiff)
    rmse = np.nan
    if ts_validation is not None:
        if len(ts_diff) != 0:
            rmse = np.sqrt(np.nanmean(ts_diff ** 2))
    ax2.set_ylabel(f'timeseries difference [{eenheid}], RMSE = %.5f'%(rmse))
    ax2.grid()
    fig.tight_layout()
    
    axs = (ax1,ax2)
    return fig, axs


def plot_HWLW_validatestats(ts_ext, ts_ext_validation):
    """
    This definition calculates (and plots and prints) some statistics when comparing
    extreme values. This is done by calculating the extreme number (sort of relative to
    Cadzand 1jan2000, but see 'warning') and subtracting the ts_ext and
    ts_ext_validation dataframes based on these numbers (and HWLWcode).
    It will only result in values for the overlapping extremes, other values will be NaN
    and are not considered for the statistics.
    Warning: the calculated extreme numbers in this definition are not corrected for the
    real phase difference with the M2phasediff argument, the calculated extreme are fine
    for internal use (to match corresponding extremes) but the absolute number might be
    incorrect.

    Parameters
    ----------
    ts_ext : pandas.DataFrame
        The DataFrame should contain a 'values' and 'HWLW_code' column and a
        pd.DatetimeIndex as index, it contains the times, values and codes of the
        timeseries that are extremes.
    ts_ext_validation : pandas.DataFrame
        The DataFrame should contain a 'values' and 'HWLW_code' column and a
        pd.DatetimeIndex as index, values and codes of the timeseries that are extremes.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The generated figure handle, with which the figure can be adapted and saved.
    axs : (tuple of) matplotlib.axes._subplots.AxesSubplot
        The generated axis handle, whith which the figure can be adapted.

    """
    
    logger.info('Calculating comparison statistics for extremes')
    if 'HWLWno' not in ts_ext.columns or 'HWLWno' not in ts_ext_validation.columns:
        logger.info(
            'HWLWno is not present in ts_ext or ts_ext_validation, trying to '
            'automatically derive it without station argument (this might fail)'
            )
        ts_ext = calc_HWLWnumbering(ts_ext=ts_ext)
        ts_ext_validation = calc_HWLWnumbering(ts_ext=ts_ext_validation)
    
    #set HWLWcode and HWLWno as index, to make easy subtraction possible
    ts_ext['times'] = ts_ext.index
    ts_ext = ts_ext.set_index(['HWLWcode','HWLWno'],drop=False)
    ts_ext_validation['times'] = ts_ext_validation.index
    ts_ext_validation = ts_ext_validation.set_index(['HWLWcode','HWLWno'],drop=False)
    hwlw_diff = ts_ext.sub(ts_ext_validation[['times','values']])
    
    tdiff_minutes = hwlw_diff['times'].dt.total_seconds()/60
    assert hwlw_diff.attrs["eenheid"] == "m"
    vdiff_cm = hwlw_diff['values']*100
    message = ('Time differences [minutes]\n'
               f'  RMSE: {(np.sqrt(np.mean(tdiff_minutes**2))):.2f}\n'
               f'  std: {tdiff_minutes.std():.2f}\n'
               f'  abs max: {tdiff_minutes.abs().max():.2f}\n'
               f'  abs mean: {tdiff_minutes.abs().mean():.2f}\n'
               f'  #NaN: {tdiff_minutes.isnull().sum()} of {len(vdiff_cm)}')
    logger.info(message)
    message = ('Value differences [cm]\n'
               f'  RMSE: {(np.sqrt(np.mean(vdiff_cm**2))):.2f}\n'
               f'  std: {vdiff_cm.std():.2f}\n'
               f'  abs max: {vdiff_cm.abs().max():.2f}\n'
               f'  abs mean: {vdiff_cm.abs().mean():.2f}\n'
               f'  #NaN: {vdiff_cm.isnull().sum()} of {len(vdiff_cm)}')
    logger.info(message)
    
    fig, ax = plt.subplots()
    ax.plot(hwlw_diff.loc[1,'times'].dt.total_seconds()/60,hwlw_diff.loc[1,'values']*100,'+',label='HWdiff')
    ax.plot(hwlw_diff.loc[2,'times'].dt.total_seconds()/60,hwlw_diff.loc[2,'values']*100,'.',label='LWdiff')
    ax.set_xlabel('Time difference [minutes]')
    ax.set_ylabel('Value difference [cm]')
    ax.legend(loc=1)
    ax.grid()

    return fig, ax


def write_netcdf(ts, filename, ts_ext=None, nosidx=False, mode='w'):
    """
    Writes the timeseries to a netCDF file

    Parameters
    ----------
    ts : pandas.DataFrame
        The DataFrame should contain a 'values' column and a pd.DatetimeIndex as index,
        it contains the timeseries.
    station : str
        DESCRIPTION.
    vertref : str
        DESCRIPTION.
    filename : str
        The filename of the netCDF file that will be written.
    ts_ext : pandas.DataFrame, optional
        The DataFrame should contain a 'values' and 'HWLW_code' column and a
        pd.DatetimeIndex as index, it contains the times, values and codes of the
        timeseries that are extremes. The default is None.
    tzone_hr : int, optional
        The timezone (GMT+tzone_hr) that applies to the data. The default is 1 (MET).
    
    Returns
    -------
    None.
    
    """
    
    import hatyan
    version_no = hatyan.__version__
    
    # get metadata from ts dataframe
    metadata_ts = metadata_from_obj(ts)
    station = metadata_ts.pop("station")
    vertref = metadata_ts.pop("vertref")
    if ts_ext is not None:
        metadata_ext = metadata_from_obj(ts_ext)
        station_ext = metadata_ext.pop("station")
        vertref_ext = metadata_ext.pop("vertref")
        assert station==station_ext
        assert vertref==vertref_ext
    
    times_all = ts.index
    timeseries = ts['values']
    dt_analysistime = dt.datetime.now()
    data_nc = Dataset(filename, mode, format="NETCDF3_CLASSIC")
    attr_dict = {'title': 'tidal prediction for %s to %s'%(times_all[0].strftime('%Y-%m-%d %H:%M:%S'), times_all[-1].strftime('%Y-%m-%d %H:%M:%S')),
                 'institution': 'Rijkswaterstaat',
                 'source': 'hatyan-%s tidal analysis program of Rijkswaterstaat'%(version_no),
                 }
    data_nc.setncatts(attr_dict)
    
    ncvarlist = list(data_nc.variables.keys())
    ncdimlist = list(data_nc.dimensions.keys())
    statname_len = 64
    
    if 'stations' not in ncdimlist:
        data_nc.createDimension('stations',None)
    if 'statname_len' not in ncdimlist:
        data_nc.createDimension('statname_len',statname_len)
    if 'time' not in ncdimlist:
        data_nc.createDimension('time',len(times_all.tolist()))
    if 'analysis_time' not in ncdimlist:
        data_nc.createDimension('analysis_time',1)
    
    refdate_tz = dt.datetime(1900,1,1,tzinfo=ts.index.tz)
    time_units = 'minutes since %s'%(refdate_tz.strftime('%Y-%m-%d %H:%M:%S %z'))
    dict_statattr = {'cf_role': 'timeseries_id'}
    dict_anatimattr = {'units': time_units, 'standard_name':'forecast_reference_time', 'long_name':'forecast_reference_time'}
    dict_timattr = {'units': time_units}
    dict_wlattr = {'units':'m', 'vertical_reference': vertref,
                   'standard_name': 'sea_surface_height_above_geopotential_datum',
                   'long_name': 'astronomical prediction of water level above reference level',
                   }
    dict_HWattr = {'units':'m', 'vertical_reference': vertref,
                   'standard_name': 'sea_surface_height_above_geopotential_datum',
                   'long_name': 'astronomical prediction of high water extremes above reference level',
                   }
    dict_LWattr = {'units':'m', 'vertical_reference': vertref,
                   'standard_name': 'sea_surface_height_above_geopotential_datum',
                   'long_name': 'astronomical prediction of low water extremes above reference level',
                   }
    dict_HWLWnoattr = {'units':'n-th tidal wave since reference wave at Cadzand on 1-1-2000'} #, 'standard_name': '', 'long_name': ''}

    if 'stations' not in ncvarlist: #create empty variables if not yet present
        nc_newvar = data_nc.createVariable('stations','S1',('stations','statname_len',))
        nc_newvar.setncatts(dict_statattr)
    
    if 'analysis_time' not in ncvarlist:
        nc_newvar = data_nc.createVariable('analysis_time','f8',('analysis_time',))
        nc_newvar.setncatts(dict_anatimattr)
        data_nc.variables['analysis_time'][0] = date2num([dt_analysistime], units=data_nc.variables['analysis_time'].units)   
    
    #current length is used as index
    nstat = data_nc.variables['stations'].shape[0]
    #append current data to netcdf files
    data_nc.variables['stations'][nstat,:] = stringtoarr(station, statname_len, dtype='S')
    
    #general prediction
    if 'time' not in ncvarlist:
        nc_newvar = data_nc.createVariable('time','f8',('time',))
        nc_newvar.setncatts(dict_timattr)
        #set time contents upon creation of variable, is constant over loop
        data_nc.variables['time'][:] = date2num(times_all.tolist(),units=data_nc.variables['time'].units)
    if 'waterlevel_astro' not in ncvarlist:
        nc_newvar = data_nc.createVariable('waterlevel_astro','f8',('stations','time',))
        nc_newvar.setncatts(dict_wlattr)
    data_nc.variables['waterlevel_astro'][nstat,:] = timeseries
    
    if ts_ext is None:
        logger.info('no HWLW prediction written')
        data_nc.close()
        return #this skips the HWLW part of the definition
    
    #HWLW prediction
    if nosidx:
        #convert index from time to HWLWno
        data_HWLW_nosidx = ts_ext.copy()
        data_HWLW_nosidx['times'] = data_HWLW_nosidx.index
        data_HWLW_nosidx = data_HWLW_nosidx.set_index('HWLWno')
        HWLWno_all = data_HWLW_nosidx.index.unique()
        data_HW = data_HWLW_nosidx.loc[data_HWLW_nosidx['HWLWcode']==1]
        data_LW = data_HWLW_nosidx.loc[data_HWLW_nosidx['HWLWcode']==2]
        bool_HW = HWLWno_all.isin(data_HW.index)
        bool_LW = HWLWno_all.isin(data_LW.index)
        
        #HWLWno
        if 'HWLWno' not in ncdimlist:
            data_nc.createDimension('HWLWno',len(HWLWno_all))
        if 'HWLWno' not in ncvarlist:
            nc_newvar = data_nc.createVariable('HWLWno','i',('HWLWno',))
            nc_newvar.setncatts(dict_HWLWnoattr)
        data_nc.variables['HWLWno'][:] = HWLWno_all
        #HW
        if 'times_astro_HW' not in ncvarlist:
            nc_newvar = data_nc.createVariable('times_astro_HW','f8',('stations','HWLWno',))
            nc_newvar.setncatts(dict_timattr)
        data_nc.variables['times_astro_HW'][nstat,bool_HW] = date2num(data_HW['times'].tolist(),units=data_nc.variables['times_astro_HW'].units)
        if 'waterlevel_astro_HW' not in ncvarlist:
            nc_newvar = data_nc.createVariable('waterlevel_astro_HW','f8',('stations','HWLWno',))
            nc_newvar.setncatts(dict_HWattr)
        data_nc.variables['waterlevel_astro_HW'][nstat,bool_HW] = data_HW['values']
        #LW
        if 'times_astro_LW' not in ncvarlist:
            nc_newvar = data_nc.createVariable('times_astro_LW','f8',('stations','HWLWno',)) 
            nc_newvar.setncatts(dict_timattr)
        data_nc.variables['times_astro_LW'][nstat,bool_LW] = date2num(data_LW['times'].tolist(),units=data_nc.variables['times_astro_LW'].units)
        if 'waterlevel_astro_LW' not in ncvarlist:
            nc_newvar = data_nc.createVariable('waterlevel_astro_LW','f8',('stations','HWLWno',))
            nc_newvar.setncatts(dict_LWattr)
        data_nc.variables['waterlevel_astro_LW'][nstat,bool_LW] = data_LW['values']
    
    else: #use time as index and create array with gaps (not possible to combine multiple stations)
        if nstat>0:
            raise Exception(f'with nosidx={nosidx} it is not possible to write multiple stations per file')
        data_HWLW = ts_ext.copy()
        data_HWLW = data_HWLW.sort_index(axis=0)
        data_HW = data_HWLW[data_HWLW['HWLWcode']==1]
        data_LW = data_HWLW[data_HWLW['HWLWcode']==2]
        #create empty variables if not yet present
        
        #HW
        if 'time_HW' not in ncdimlist:
            data_nc.createDimension('time_HW',len(data_HW))
        if 'time_HW' not in ncvarlist:
            nc_newvar = data_nc.createVariable('time_HW','f8',('time_HW',))
            nc_newvar.setncatts(dict_timattr)
        data_nc.variables['time_HW'][:] = date2num(data_HW.index.tolist(),units=data_nc.variables['time_HW'].units)
        if 'waterlevel_astro_HW' not in ncvarlist:
            nc_newvar = data_nc.createVariable('waterlevel_astro_HW','f8',('stations','time_HW',))
            nc_newvar.setncatts(dict_HWattr)
        data_nc.variables['waterlevel_astro_HW'][nstat,:] = data_HW['values']
        
        #LW
        if 'time_LW' not in ncdimlist:
            data_nc.createDimension('time_LW',len(data_LW))
        if 'time_LW' not in ncvarlist:
            nc_newvar = data_nc.createVariable('time_LW','f8',('time_LW',))
            nc_newvar.setncatts(dict_timattr)
        data_nc.variables['time_LW'][:] = date2num(data_LW.index.tolist(),units=data_nc.variables['time_LW'].units)
        if 'waterlevel_astro_LW' not in ncvarlist:
            nc_newvar = data_nc.createVariable('waterlevel_astro_LW','f8',('stations','time_LW',))
            nc_newvar.setncatts(dict_LWattr)
        data_nc.variables['waterlevel_astro_LW'][nstat,:] = data_LW['values']
        
        #HWLW numbering
        if 'HWLWno' in ts_ext.columns:
            if 'waterlevel_astro_HW_numbers' not in ncvarlist:
                nc_newvar = data_nc.createVariable('waterlevel_astro_HW_numbers','i4',('stations','time_HW',))
                #nc_newvar.setncatts(dict_HWattr)
            data_nc.variables['waterlevel_astro_HW_numbers'][nstat,:] = data_HW['HWLWno']
            if 'waterlevel_astro_LW_numbers' not in ncvarlist:
                nc_newvar = data_nc.createVariable('waterlevel_astro_LW_numbers','i4',('stations','time_LW',))
                #nc_newvar.setncatts(dict_LWattr)
            data_nc.variables['waterlevel_astro_LW_numbers'][nstat,:] = data_LW['HWLWno']
    
    data_nc.close()
    return


def write_dia(ts, filename, headerformat='dia'):
    """
    Writes the timeseries to an equidistant dia file or the extremes to 
    a non-equidistant dia file. This is only supported 
    for timeseries with a UTC+1 timestamp, since DONAR (and therefore dia)
    data is always in UTC+1 (MET/CET), also during summertime periods.

    Parameters
    ----------
    ts : pandas.DataFrame
        The DataFrame should contain a 'values' column and a pd.DatetimeIndex as index, it contains the timeseries. 
        In case of extremes, the DataFrame should also contain a 'HWLW_code' column.
    filename : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    if "status" in ts.columns:
       logger.warning("status column is ignored by hatyan.write_dia(), all "
                      "status values in diafile will be 'Ongecontroleerd'")
    if "qualitycode" in ts.columns:
       logger.warning("qualitycode column is ignored by hatyan.write_dia(), all "
                      "qualitycode values in diafile will be 0")
    
    # optionally convert meters to centimeters
    # before asserting metadata in wns_from_metadata
    if ts.attrs["eenheid"] == "m":
        ts = ts.copy()
        ts["values"] *= 100
        ts.attrs["eenheid"] = "cm"
    
    metadata_pd = get_metadata_pd(ts, headerformat=headerformat)
    
    if "HWLWcode" in ts.columns:
        data_todia = (ts.index.strftime('%Y%m%d;%H%M')
                      + ';'
                      + ts['HWLWcode'].astype(str)
                      + '/0;'
                      + (ts['values']).round().astype(int).astype(str)
                      + ':')
    else:
        linestr_list = []
        linestr = ''
        for iV, ts_value in enumerate(ts['values']):
            linestr_add = "%i/0:"%(np.round(ts_value))
            linestr = linestr + linestr_add
            if (len(linestr) > 114) or (iV==len(ts)-1): # append linestr to linestr_list if linestr is longer than n characters or last item of ts_values was reached
                linestr_list.append(linestr)
                linestr = ''
        data_todia = pd.Series(linestr_list)
    
    # write metadata and data to file
    with io.open(filename,'w', newline='\n') as f: #open file and set linux newline style
        for metaline in metadata_pd:
            f.write('%s\n'%(metaline))
        data_todia.to_csv(f,index=False,header=False)


def get_metadata_pd(ts, headerformat):
    metadata = metadata_from_obj(ts)
    waarnemingssoort = wns_from_metadata(metadata)
    vertref = metadata['vertref']
    station = metadata['station']
    grootheid = metadata['grootheid']
    
    # check grootheid
    if grootheid == 'WATHTBRKD':
        grootheid = 'WATHTBRKD;Waterhoogte berekend'
        ana = 'F012;Waterhoogte astronomisch mbv harmonische analyse'
    else:
        raise ValueError(f'write_dia() expects quantity WATHTBRKD, but {grootheid} was provided.')
    
    # check tzone
    tzone = ts.index.tz
    if tzone not in [pytz.FixedOffset(60), dt.timezone(dt.timedelta(seconds=3600))]:
        raise ValueError('write_dia() expects tzone pytz.FixedOffset(60) (since tzone '
                         f'is not defined in dia-header), but {tzone} was provided.')
    
    # check vertref
    if vertref == 'NAP':
        vertreflong = 'T.o.v. Normaal Amsterdams Peil'
        parameterX = 'GETETBRKD2;Getijextreem berekend' # only for ext
    elif vertref == 'MSL':
        vertreflong = 'T.o.v. Mean Sea Level'
        parameterX = 'GETETBRKDMSL2;Getijextreem berekend t.o.v. MSL' # only for ext
    else:
        raise Exception('ERROR: currently only vertref="NAP" and vertref="MSL" are supported for writing diafiles')
    
    # get times
    time_today = dt.datetime.today().strftime('%Y%m%d')
    tstart_str = ts.index[0].strftime('%Y%m%d;%H%M')
    tstop_str = ts.index[-1].strftime('%Y%m%d;%H%M')
    
    if "HWLWcode" in ts.columns:
        if 11 in ts['HWLWcode'].values or 22 in ts['HWLWcode'].values:
            raise Exception('ERROR: invalid HWLWcodes in provided extreme timeseries (11 and/or 22)')
        
        metadata_pd = pd.Series(['[IDT;*DIF*;A;;%6s]'%(time_today), #identificatieblok (DIF voor dia en WIF voor wia, A voor ASCII) #TODO: kan *DIF*/*WIF* gebruikt worden voor identificatie dia/wia file?
                                 '[W3H]', #WIE, WAT, WAAR en HOE
                                 'MUX;%s'%(parameterX), #Mux
                                 'ANI;RIKZITSDHG;RIKZ - afdeling ZDI te Den Haag', #Analyserende-instantie
                                 'BHI;RIKZITSDHG;RIKZ - afdeling ZDI te Den Haag', #Beherende-instantie
                                 'BMI;NVT;Niet van toepassing', #Bemonsterende-instantie
                                 'OGI;RIKZMON_WAT;RIKZ - Landelijke monitoring waterhoogten gegevens', #Opdrachtgevende-instantie
                                 'LOC;%s'%(station), #Locatiecode;Omschrijving;Soort;Coördinaattype;X-coördinaat_GS;Y-coördinaat_GS, EPSG_code
                                 'ANA;%s'%(ana), #WBM in wia: Waardebepalingsmethode
                                 'TYP;TN', #Reekstype: niet-equidistant
                                 '[MUX]', #Multiplex administratieblok
                                 'MXW;1;15', #TODO: niet ondersteund in wia, wellicht niet essentieel voor dia?
                                 'MXP;1;GETETCDE;Getijextreem code', #MXG in wia: grootheid muxkanaal (TODO: wia GETETTPE in welke groep?)
                                 'MXC;1;10;Oppervlaktewater', #Compartiment Muxkanaal
                                 'MXE;1;T;DIMSLS', #domein (T: integer met waarde 1 tot 127, vertaaltabel TYPERING met dezelfde parameter/compartiment-combinatie is dan verplicht) en eenheid
                                 'MXH;1;NVT;Niet van toepassing', #Hoedanigheid Muxkanaal
                                 'MXW;2;%i'%(waarnemingssoort), #TODO: niet ondersteund in wia, wellicht niet essentieel voor dia?
                                 'MXP;2;%s'%(grootheid), #MXG in wia: grootheid muxkanaal
                                 'MXC;2;10;Oppervlaktewater', #Compartiment Muxkanaal
                                 'MXE;2;I;cm', #domein (I: integer) en eenheid
                                 'MXH;2;%s;%s'%(vertref, vertreflong), #Hoedanigheid Muxkanaal
                                 '[TYP]',
                                 'TVL;1;1;hoogwater',
                                 'TVL;1;2;laagwater',
                                 'TVL;1;3;laagwater 1',
                                 'TVL;1;4;topagger',
                                 'TVL;1;5;laagwater 2',
                                 '[RKS]',
                                 'TYD;%10s;%10s'%(tstart_str,tstop_str),
                                 '[TPS]',
                                 'STA;%10s;%10s;O'%(tstart_str,tstop_str),
                                 '[WRD]'])
    else:
        assert ts.index.freq is not None
        timestep_min = pd.Timedelta(ts.index.freq).total_seconds()/60
        
        #informatie in comments komt veelal uit "IDD-WIA-v0.9.2.docx"
        metadata_pd = pd.Series(['[IDT;*DIF*;A;;%6s]'%(time_today), #identificatieblok (DIF voor dia en WIF voor wia, A voor ASCII) #TODO: kan *DIF*/*WIF* gebruikt worden voor identificatie dia/wia file?
                                 '[W3H]', #WIE, WAT, WAAR en HOE
                                 'WNS;%i'%(waarnemingssoort), #TODO: niet ondersteund in wia, wellicht niet essentieel voor dia dus geheel weglaten?
                                 'PAR;%s'%(grootheid), #parameter/grootheid, gelijk voor waarnemingssoorten 18 en 55. GHD in wia (PAR is daar parameter, maar betekent wat anders)
                                 'CPM;10;Oppervlaktewater', #compartiment, gelijk voor waarnemingssoorten 18 en 55
                                 'EHD;I;cm', #domein (I: integer) en eenheid, gelijk voor waarnemingssoorten 18 en 55
                                 'HDH;%s;%s'%(vertref,vertreflong),
                                 'ANI;RIKZITSDHG;RIKZ - afdeling ZDI te Den Haag', #Analyserende-instantie
                                 'BHI;RIKZITSDHG;RIKZ - afdeling ZDI te Den Haag', #Beherende-instantie
                                 'BMI;NVT;Niet van toepassing', #Bemonsterende-instantie
                                 'OGI;RIKZMON_WAT;RIKZ - Landelijke monitoring waterhoogten gegevens', #Opdrachtgevende-instantie
                                 'LOC;%s'%(station), #Locatiecode;Omschrijving;Soort;Coördinaattype;X-coördinaat_GS;Y-coördinaat_GS, EPSG_code #;Hoek van Holland;P;RD;6793000;44400000
                                 'ANA;%s'%(ana), #WBM in wia: Waardebepalingsmethode
                                 'TYP;TE', #Reekstype: equidistant
                                 '[RKS]',
                                 'TYD;%10s;%10s;%i;min'%(tstart_str,tstop_str,timestep_min), #(Begin)datum;(Begin)tijdstip;Einddatum;Eindtijdstip;Tijdstap;Tijdstapeenheid
                                  '[TPS]',
                                 'STA;%10s;%10s;O'%(tstart_str,tstop_str), #Statuscode;Begindatum;(Begin)tijdstip;Einddatum;Eindtijdstip;Tijdstap;
                                 '[WRD]'])
    
    if headerformat=='wia':
        metadata_pd = dia_metadata_to_wia_metadata(metadata_pd, time_today, ana, grootheid)
    return metadata_pd
    

def dia_metadata_to_wia_metadata(metadata_pd, time_today, ana, grootheid):
    # drop metadata
    for metalinestart in ['WNS','MXW']:
        bool_drop = metadata_pd.str.startswith(metalinestart)
        metadata_pd = metadata_pd[~bool_drop]
    
    # rename metadata for ts and ext
    metadata_pd[metadata_pd.str.startswith('[IDT;')] = '[IDT;*WIF*;A;;%6s]'%(time_today)
    metadata_pd[metadata_pd.str.startswith('ANA;')] = 'WBM;other:%s'%(ana)
    # rename metadata for ts
    metadata_pd[metadata_pd.str.startswith('PAR;')] = 'GHD;%s'%(grootheid)
    metadata_pd[metadata_pd.str.startswith('CPM;')] = 'CPM;OW;Oppervlaktewater'
    # rename metadata for ext
    metadata_pd[metadata_pd.str.startswith('MXP;1')] = 'MXT;1;GETETTPE' # GETETCDE;Getijextreem code naar GETETTPE #TODO: MXT wordt niet ondersteund door wia, toch MXG?
    metadata_pd[metadata_pd.str.startswith('MXC;1')] = 'MXC;1;OW;Oppervlaktewater'
    metadata_pd[metadata_pd.str.startswith('MXP;2')] = 'MXG;2;%s'%(grootheid)
    metadata_pd[metadata_pd.str.startswith('MXC;2')] = 'MXC;2;OW;Oppervlaktewater'
    return metadata_pd


def write_noos(ts, filename):
    """
    Writes the timeseries to a noos file

    Parameters
    ----------
    ts : pandas.DataFrame
        The DataFrame should contain a 'values' column and a pd.DatetimeIndex as index, it contains the timeseries.
    filename : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    metadata = ts.attrs
    
    timestamp = dt.datetime.now().strftime('%c')
    ts_out = pd.DataFrame({'times':ts.index.strftime('%Y%m%d%H%M'),'values':ts['values']})
    
    if metadata is None:
        header_txt = f"""------------------------------------------------------
        Timeseries written by hatyan
        Created at {timestamp}
        ------------------------------------------------------
        Location    : None
        Position    : None
        Source      : None
        Unit        : waterlevel
        Analyse time: 000000000000
        Timezone    : None
        ------------------------------------------------------"""
    else:
        header_txt = f"""------------------------------------------------------
        Timeseries written by hatyan\nCreated at {timestamp}
        ------------------------------------------------------
        """
        for key in metadata.keys():
            header_txt = header_txt+('%-12s: %s\n'%(key, metadata[key]))
        header_txt = header_txt+'------------------------------------------------------'
    np.savetxt(filename,ts_out,fmt='%s %7.4f',header=header_txt)


def crop_timeseries(ts, times, onlyfull=True):
    """
    Crops the provided timeseries

    Parameters
    ----------
    ts : pandas.DataFrame
        The DataFrame should contain a 'values' column and a pd.DatetimeIndex as index, it contains the timeseries.
    times : slice
        slice(tstart,tstop).

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    ts_pd_out : TYPE
        DESCRIPTION.

    """
    ts_pd_in = ts
    tstart = pd.Timestamp(times.start)
    tstop = pd.Timestamp(times.stop)
    
    logger.info('cropping timeseries')
    if tstart >= tstop:
        raise ValueError(f'the tstart and tstop should be increasing, but they are not: {times}.')
    if (tstart < ts_pd_in.index.min()) or (tstop > ts_pd_in.index.max()):
        message = (
            'imported timeseries is not available within entire requested period:\n'
            f'requested period:    {tstart} to {tstop}\n'
            f'imported timeseries: {ts_pd_in.index[0]} to {ts_pd_in.index[-1]}'
            )
        if onlyfull:
            raise ValueError(message)
        else:
            logger.warning(message)
            
    ts_pd_out = ts_pd_in.loc[tstart:tstop]
    
    # add metadata
    metadata = metadata_from_obj(ts)
    ts_pd_out = metadata_add_to_obj(ts_pd_out,metadata)
    return ts_pd_out


def resample_timeseries(ts, timestep_min, tstart=None, tstop=None):
    """
    resamples the provided timeseries, only overlapping timesteps are selected, so no
    interpolation. with tstart/tstop it is possible to extend the timeseries with NaN
    values.

    Parameters
    ----------
    ts : pandas.DataFrame
        The DataFrame should contain a 'values' and 'HWLW_code' column and a
        pd.DatetimeIndex as index, it contains the timeseries to be resampled.
    timestep_min : int
        the amount of minutes with which to resample the timeseries.
    tstart : dt.datetime, optional
        the start date for the resampled timeseries, the default is None which results
        in using the start date of the input ts.
    tstop : dt.datetime, optional
        the stop date for the resampled timeseries, the default is None which results
        in using the stop date of the input ts.

    Returns
    -------
    data_pd_resample : pandas.DataFrame with a 'values' column and a pd.DatetimeIndex as index
        the resampled timeseries.

    """
    
    logger.info(f'resampling timeseries to {timestep_min} minutes')
    
    bool_duplicated_index = ts.index.duplicated()
    if bool_duplicated_index.sum()>0:
        raise ValueError(
            'there are duplicated values in the ts DatetimeIndex, this is not '
            'supported by Timeseries.resample_timeseries(). Try '
            '"ts_nodupl = ts[~ts.index.duplicated()]"'
            )
    
    if tstart is None:
        tstart = ts.index[0]
    if tstop is None:
        tstop = ts.index[-1]
    # generate timeseries with correct tstart/tstop and interval
    data_pd_resample = pd.DataFrame({},index=pd.date_range(tstart,tstop,freq='%dmin'%(timestep_min)))
    # put measurements into this timeseries, matches to correct index automatically
    data_pd_resample['values'] = ts['values']
    
    # add metadata
    metadata = metadata_from_obj(ts)
    data_pd_resample = metadata_add_to_obj(data_pd_resample,metadata)
    
    return data_pd_resample


def nyquist_folding(ts_pd,t_const_freq_pd):
    #deriving dominant timestep, sampling frequency and nyquist frequency
    timestep_hr_all = ((ts_pd.index[1:]-ts_pd.index[:-1]).total_seconds()/3600)
    uniq_vals, uniq_counts = np.unique(timestep_hr_all,return_counts=True)
    timestep_hr_dominant = uniq_vals[np.argmax(uniq_counts)]
    fs_phr = 1/timestep_hr_dominant #sampling freq [1/hr]
    fn_phr = fs_phr/2 #nyquist freq [1/hr]
    
    # check if there is a component with exactly the same frequency as the nyquist
    # frequency #TODO: or its multiples (but with remainder, A0 also gets flagged)
    bool_isnyquist = t_const_freq_pd['freq']==fn_phr
    if bool_isnyquist.any():
        raise ValueError(f'there is a component on the Nyquist frequency ({fn_phr} [1/hr]), '
                        f'this not possible:\n{t_const_freq_pd.loc[bool_isnyquist]}')
    
    logger.info('folding frequencies over Nyquist frequency, which is half '
                f'of the dominant timestep ({timestep_hr_dominant} hour), '
                f'there are {len(uniq_counts)} unique timesteps)')
    #folding frequencies over nyquist frequency: https://users.encs.concordia.ca/~kadem/CHAPTER%20V.pdf (fig 5.6)
    # remainder gives folded frequencies for even divisions (freqs for odd divisions are not valid yet)
    freq_div,freq_rem = np.divmod(t_const_freq_pd[['freq']],fn_phr)
    freq_div_isodd = (freq_div%2).astype(bool)
    # fn-remainder gives folded frequencies for odd divisions
    freq_rem[freq_div_isodd] = fn_phr-freq_rem[freq_div_isodd]
    return freq_rem


def check_rayleigh(ts_pd,t_const_freq_pd):
    """
    The Rayleigh criterion: |freq1-freq2| * T > R, where T is the period length of the timeseries, R is the Rayleigh number.
    This comes down to: ts_period / period_difference > R
    
    Parameters
    ----------
    ts_pd : TYPE
        DESCRIPTION.
    t_const_freq_pd : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    t_const_freq = t_const_freq_pd.sort_values('freq')['freq'].drop('A0',errors='ignore')
    if not len(t_const_freq)>1:
        return #Rayleigh check is only relevant (and possible) if there more than one non-A0 component, otherwise stop.
    freq_diffs = np.diff(t_const_freq)
    ts_period_hr = (ts_pd.index.max()-ts_pd.index.min()).total_seconds()/3600
    # Koos Doekes: "Bij het algoritme dat HATYAN gebruikt mag men in de praktijk het
    # Rayleigh-criterium enigszins schenden, tot zo'n 0,7 van de theoretisch vereiste reekslengte. "
    rayleigh_tresh = 0.7
    # TODO: might be better to drop timeseries nan-values first, especially from
    # start/end of series since it decreases ts_period_hr
    rayleigh = ts_period_hr*freq_diffs
    freq_diff_phr_minimum = rayleigh_tresh/ts_period_hr
    rayleigh_bool = rayleigh>rayleigh_tresh
    rayleigh_bool_id = np.nonzero(~rayleigh_bool)[0]
    
    if rayleigh_bool.all():
        logger.info('Rayleigh criterion OK (always>%.2f, minimum is %.2f)'%(rayleigh_tresh, np.min(rayleigh)))
        logger.info('Frequencies are far enough apart (always >%.6f, minimum is %.6f)'%(freq_diff_phr_minimum,np.min(freq_diffs)))
    else:
        logger.info('Rayleigh criterion vandalised (not always>%.2f, minimum is %.2f)'%(rayleigh_tresh, np.min(rayleigh)))
        logger.info('Frequencies with not enough difference (not always >%.6f, minimum is %.6f)'%(freq_diff_phr_minimum,np.min(freq_diffs)))
        for ray_id in rayleigh_bool_id:
            t_const_freq_sel = t_const_freq.iloc[[ray_id,ray_id+1]]
            t_const_freq_sel['diff'] = np.diff(t_const_freq_sel.values)[0]
            t_const_freq_sel['ndays min'] = rayleigh_tresh/np.diff(t_const_freq_sel.values)[0]/24
            logger.info(t_const_freq_sel)
            if t_const_freq_sel['diff'] < 1e-9:
                logger.warning(f'frequency difference between {t_const_freq_sel.index[0]} and {t_const_freq_sel.index[1]} almost zero, will result in ill conditioned matrix')


class Timeseries_Statistics:
    """
    returns several statistics of the provided timeseries as a Timeseries_Statistics class, which is a like a dict that pretty prints automatically.

    Parameters
    ----------
    ts : pandas.DataFrame
        The DataFrame should contain a 'values' column and a pd.DatetimeIndex as index, it contains the timeseries to be checked.

    Returns
    -------
    stats: class Timeseries_Statistics
        Timeseries_Statistics is a like a dict that pretty prints automatically.

    """
    #TODO: make like a dict with different __str__ method, instead of this mess https://stackoverflow.com/questions/4014621/a-python-class-that-acts-like-dict
    #TODO: improve output dict, keys are now not convenient to use. Maybe make keys and longname?
    def __init__(self,ts):
        timesteps_min_all = ts.index.to_series().diff()[1:].dt.total_seconds()/60
        bool_int = np.abs(timesteps_min_all-timesteps_min_all.round(0))<1e-9
        if bool_int.all():
            timesteps_min_all = timesteps_min_all.astype(int)
        
        #uniq timesteps and counts
        uniq_vals, uniq_counts = np.unique(timesteps_min_all,return_counts=True)
        timestep_min_valcounts_pd = pd.DataFrame({'timestep_min':uniq_vals},index=uniq_counts)
        timestep_min_valcounts_pd.index.name = 'counts'
        
        #ordering of timesteps
        if (timesteps_min_all>0).all():
            timesteps_incr_print = 'all time intervals are in increasing order and are never equal'
        else:
            timesteps_incr_print = 'the times-order of ts is not always increasing (duplicate values or wrong order)'
        
        ntimes_nonan = ts['values'].count()
        ntimes = len(ts)
        ntimesteps_uniq = len(timestep_min_valcounts_pd)
        
        if len(ts)==0:
            self.stats = {'timeseries contents':ts}
        else:
            self.stats = {'timeseries contents':ts,
                        'timeseries # unique timesteps': ntimesteps_uniq,
                        'timeseries dominant timesteps':timestep_min_valcounts_pd.sort_index(ascending=False),
                        'timeseries validity': timesteps_incr_print,
                        'timeseries length': ntimes,
                        'timeseries # nonan': ntimes_nonan,
                        'timeseries % nonan': ntimes_nonan/ntimes*100,#%.1f %
                        'timeseries # nan': ntimes-ntimes_nonan,
                        'timeseries % nan': (ntimes-ntimes_nonan)/ntimes*100, #%.1f %
                        }
    def __str__(self):
        pd_max_rows_backup = pd.options.display.max_rows #backup max_rows default setting (60 or so) and set to 15
        pd.options.display.max_rows = 10
        print_statement = ''
        for key in self.stats.keys():
            if key in ['timeseries contents','timeseries dominant timesteps']:
                print_statement += f'{key}:\n{self.stats[key]}\n'
            else:
                print_statement += f'{key}: {self.stats[key]}\n'
        pd.options.display.max_rows = pd_max_rows_backup #revert max_rows default setting after printing
        return print_statement
    def __repr__(self): #avoid printing the class name
        #return dict.__repr__
        return str(self.stats)
    """
    @classmethod
    def keys(self):
        return self.stats.keys()
    """
        
    
    
###############################
################# READING FILES
###############################


def get_diaxycoords(filename, crs):
    diablocks_pd_extra = get_diablocks(filename=filename)
    dia_x = diablocks_pd_extra['x'].values
    dia_y = diablocks_pd_extra['y'].values
    dia_epsg = diablocks_pd_extra['epsg'].astype(int).values
    
    # check if all epsg codes are the same
    if not (dia_epsg==dia_epsg[0]).all():
        raise ValueError(f"The diafile contains multiple EPSG codes, not supported yet: {dia_epsg}")
    dia_epsg = dia_epsg[0]
    
    if len(dia_x)==1:
        # to avoid DeprecationWarning: Conversion of an array with ndim > 0 to a scalar is deprecated, and will error in future.
        dia_x = dia_x[0]
        dia_y = dia_y[0]
    
    transformer = Transformer.from_crs(f'epsg:{dia_epsg}', f'epsg:{crs}', always_xy=True)
    stat_x, stat_y = transformer.transform(dia_x, dia_y)
    
    return stat_x, stat_y


def get_diablocks_startstopstation(filename):
    """
    Gets information about the data blocks present in a dia file

    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    diablocks_pd_startstopstation: pd.DataFrame
        Pandas DataFrame with 'block_starts','data_starts','data_ends','station'

    """
    
    #get list of starts/ends of datasets in diafile
    linenum_colnames = ['block_starts','data_starts','data_ends']
    diablocks_pd_startstopstation = pd.DataFrame({},columns=linenum_colnames)
    
    # add str type columns to avoid "FutureWarning: Setting an item of incompatible
    # dtype is deprecated and will raise in a future error of pandas. Value 'min' has
    # dtype incompatible with float64, please explicitly cast to a compatible dtype first."
    diablocks_pd_startstopstation['station'] = ""
    
    # 'latin1 is nodig om predictie diafile die rechtstreeks uit hatyan komen in te
    # lezen (validatietijdserie met op regel 4 (PAR) ongeldige tekens aan het einde)
    with open(filename, encoding='latin1') as f:
        block_id = -1
        for linenum, line in enumerate(f, 1):
            if '[W3H]' in line:
                block_id += 1
                diablocks_pd_startstopstation.loc[block_id,'block_starts'] = linenum
            elif '[WRD]' in line:
                diablocks_pd_startstopstation.loc[block_id,'data_starts'] = linenum
            elif 'LOC' in line:
                diablocks_pd_startstopstation.loc[block_id,'station'] = line.rstrip().split(';')[1]
    diablocks_pd_startstopstation['data_ends'] = (diablocks_pd_startstopstation['block_starts']-1).tolist()[1:]+[linenum]
    if block_id == -1:
        raise Exception('ERROR: empty dia file')
    if diablocks_pd_startstopstation.isnull().any().any():
        raise Exception('ERROR: multiple blocks in diafile, but unequal amount of start/end/datastart/stationnames')
    
    #convert columns with line numbers to integers
    diablocks_pd_startstopstation[linenum_colnames] = diablocks_pd_startstopstation[linenum_colnames].astype(int)
    
    return diablocks_pd_startstopstation


def get_diablocks(filename):
    
    logger.info('reading file: %s'%(filename))
    diablocks_pd = get_diablocks_startstopstation(filename)
    # add str type columns to avoid "FutureWarning: Setting an item of incompatible
    # dtype is deprecated and will raise in a future error of pandas. Value 'min' has
    # dtype incompatible with float64, please explicitly cast to a compatible dtype first."
    columns_str = ['TYP','groepering','grootheid','coordsys',
                   'eenheid','vertref',
                   'timestep_unit',
                   'STA']
    diablocks_pd[columns_str] = ""
    for block_id in diablocks_pd.index.tolist():
        # read diafile metadata as pandas series, prevent splitting of combined paramater names like MXH;2 by replacing ; with !
        data_meta_nrows = diablocks_pd.loc[block_id,'data_starts'] - diablocks_pd.loc[block_id,'block_starts']
        # series of metadata
        data_meta_series = pd.read_table(filename,skiprows=diablocks_pd.loc[block_id,'block_starts'],nrows=data_meta_nrows,header=None)[0]
        # wia files contain these parameters, dia files don't. Replace dia names with
        # wia names (wia files also contain PAR and MXP;2, but they should not be replaced)
        if not data_meta_series.str.contains('GHD|MXG;2').any():
            data_meta_series = data_meta_series.str.replace('PAR','GHD').str.replace('MXP;2','MXG;2')
        bool_combinedparname = (data_meta_series.str[3:6]==';1;') | (data_meta_series.str[3:6]==';2;')
        data_meta_series.loc[bool_combinedparname] = data_meta_series.loc[bool_combinedparname].str.slice_replace(3,4,'!')
        
        #get groepering and whether dia/wia is equidistant or non-equidistant
        bool_startswithmux = data_meta_series.str.startswith('MUX')
        row_TYP = data_meta_series.loc[data_meta_series.str.startswith('TYP')].iloc[0].split(';')[1]
        diablocks_pd.loc[block_id,'TYP'] = row_TYP
        if row_TYP=='TN': #bool_startswithmux.any(): #extreme waterlevel timeseries (non-equidistant)
            mincontent = ['MXG;2','LOC','MXH;2','MXE;2','TYD','STA']
            if bool_startswithmux.sum()==0:
                raise Exception(
                    f'ERROR: block_id={block_id} is of TYP={row_TYP} (non-equidistant, '
                    'extreme waterlevels), but no MUX is available in metadata header '
                    'so the file cannot be read:\n{diablocks_pd}'
                    )
            diablocks_pd.loc[block_id,'groepering'] = data_meta_series.loc[bool_startswithmux].iloc[0].split(';')[1]
        elif row_TYP=='TE': #normal waterlevel timeseries (equidistant)
            mincontent = ['GHD',  'LOC','HDH',  'EHD',  'TYD','STA'] #,WNSCPM,HDH,ANA
            diablocks_pd.loc[block_id,'groepering'] = 'NVT'
        else:
            raise Exception(f'TYP "{row_TYP}" not implemented in hatyan.read_dia()')
        
        #read all required metadata
        for get_content_sel in mincontent:
            bool_mincontent = data_meta_series.str.replace('!',';').str.startswith(get_content_sel)
            if bool_mincontent.sum()!=1:
                if get_content_sel!='STA':
                    raise Exception(f'unexpected amount of matched metadatalines ({bool_mincontent.sum()}) for {get_content_sel}')
            data_meta_mincontent = data_meta_series.loc[bool_mincontent].iloc[0].split(';') #list type
            if get_content_sel in ['GHD','MXG;2']: # Grootheid (dia/wia, dia/wia equidistant). Originally dia contains PAR and MXP;2, but they are replaced
                file_grootheidname = data_meta_mincontent[1]
                valid_grootheidnames = ['WATHTE','WATHTBRKD','NVT'] #NVT in wia files
                if file_grootheidname not in valid_grootheidnames:
                    raise Exception('ERROR: grootheid name (%s) should be in %s but is %s'%(get_content_sel, valid_grootheidnames, file_grootheidname))
                diablocks_pd.loc[block_id,'grootheid'] = file_grootheidname
            elif get_content_sel in ['LOC']: # Locatie. same in all files
                coords_pd = pd.DataFrame({'epsg_in':[28992,4326,4230], 'factor':[100,1000000,1000000]}, index=['RD','W84','E50'])
                if len(data_meta_mincontent)<7:
                    logger.warning('no coordinate data available in LOC line of dia file')
                    continue
                coordsys_str, coord_x, coord_y = data_meta_mincontent[4:]
                if coordsys_str not in coords_pd.index:
                    raise Exception('unknown coordinate system string in diafile ({coordsys_str})')
                diablocks_pd.loc[block_id,'x'] = int(coord_x)/coords_pd.loc[coordsys_str,'factor']
                diablocks_pd.loc[block_id,'y'] = int(coord_y)/coords_pd.loc[coordsys_str,'factor']
                diablocks_pd.loc[block_id,'coordsys'] = coordsys_str
                diablocks_pd.loc[block_id,'epsg'] = coords_pd.loc[coordsys_str,'epsg_in']
            elif get_content_sel in ['EHD','MXE;2']: # Eenheid. equidistant dia/wia, non-equidistant dia/wia
                file_eenheid = data_meta_mincontent[2]
                if file_eenheid != 'cm':
                    raise Exception('unknown eenheid in diafile: %s'%(file_eenheid))
                diablocks_pd.loc[block_id,'eenheid'] = file_eenheid
            elif get_content_sel in ['HDH','MXH;2']: # Hoedanigheid (NAP/MSL). equidistant dia/wia, non-equidistant dia/wia
                diablocks_pd.loc[block_id,'vertref'] = data_meta_mincontent[1]
            elif get_content_sel in ['TYD']: #Tijdstip. same in all files
                datestart = dt.datetime.strptime(data_meta_mincontent[1]+data_meta_mincontent[2], "%Y%m%d%H%M")
                datestop = dt.datetime.strptime(data_meta_mincontent[3]+data_meta_mincontent[4], "%Y%m%d%H%M")
                if len(data_meta_mincontent)==5: #nonequidistant timeseries
                    timestep_value = np.nan #TODO: None is supported by pandas 2.1.2 and maybe earlier versions, but not 2.0.3 (py39 only)
                    timestep_unit = np.nan
                elif len(data_meta_mincontent)==7: #equidistant timeseries contains also timeunit and timestep
                    timestep_unit = data_meta_mincontent[6]
                    if timestep_unit not in ['min','cs']: #minutes and 1/100 sec
                        raise Exception(f'ERROR: time unit from TYD is in unknown format (not "min" or "cs"): {timestep_unit}')
                    timestep_value = int(data_meta_mincontent[5])
                else:
                    raise Exception(f'ERROR: time metadata is not understood: {data_meta_mincontent}')
                diablocks_pd.loc[block_id,'tstart'] = datestart
                diablocks_pd.loc[block_id,'tstop'] = datestop
                diablocks_pd.loc[block_id,'timestep_min'] = timestep_value
                diablocks_pd.loc[block_id,'timestep_unit'] = timestep_unit
            elif get_content_sel in ['STA']: #Status. same in all files
                diablocks_pd.loc[block_id,'STA'] = '!'.join(data_meta_series.loc[bool_mincontent].tolist())
    return diablocks_pd


def read_dia_nonequidistant(filename, diablocks_pd, block_id):

    data_nrows = diablocks_pd.loc[block_id,'data_ends'] - diablocks_pd.loc[block_id,'data_starts']
    skiprows = diablocks_pd.loc[block_id,'data_starts']
    data_pd_HWLW = pd.read_csv(filename, skiprows=skiprows,nrows=data_nrows, header=None, sep=';',
                               names=['date','time','HWLWcode/qualitycode','valuecm:'], 
                               dtype={'date':str,'time':str},
                               )
    dates_pd = data_pd_HWLW.pop('date')
    times_pd = data_pd_HWLW.pop('time')
    data_pd_HWLW['times'] = pd.to_datetime(dates_pd+times_pd)
    
    #convert HWLW+quality code to separate columns
    data_pd_HWLWtemp = data_pd_HWLW.loc[:,'HWLWcode/qualitycode'].str.split('/', expand=True)
    data_pd_HWLW['HWLWcode'] = data_pd_HWLWtemp.iloc[:,0].astype('int')
    data_pd_HWLW['qualitycode'] = data_pd_HWLWtemp.iloc[:,1].astype('int')
    data_pd_HWLW = data_pd_HWLW.drop('HWLWcode/qualitycode',axis='columns')

    #construct df
    data_pd_HWLW['values'] = data_pd_HWLW['valuecm:'].str.strip(':').astype('int')
    data_pd_HWLW = data_pd_HWLW.drop('valuecm:',axis='columns')
    
    bool_hiaat = data_pd_HWLW['qualitycode'] == 99
    data_pd_HWLW.loc[bool_hiaat,'values'] = np.nan
    
    data_pd = data_pd_HWLW
    data_pd = data_pd.set_index('times')
    data_pd.index = data_pd.index.tz_localize(pytz.FixedOffset(60))
    
    # add metadata
    metadata = metadata_from_diablocks(diablocks_pd, block_id)
    data_pd = metadata_add_to_obj(data_pd,metadata)
    
    return data_pd


def read_dia_equidistant(filename, diablocks_pd, block_id):
    
    datestart = diablocks_pd.loc[block_id,'tstart']
    datestop = diablocks_pd.loc[block_id,'tstop']
    timestep_min = diablocks_pd.loc[block_id,'timestep_min']
    timestep_unit = diablocks_pd.loc[block_id,'timestep_unit']
    tzone = pytz.FixedOffset(60)
    if timestep_unit=='min':
        times_fromfile = pd.date_range(start=datestart,end=datestop,freq='%dmin'%(timestep_min), tz=tzone)
    else:
        times_fromfile = pd.date_range(start=datestart,end=datestop,freq=f'{timestep_min*10000000} ns', tz=tzone)
    
    # get data for station
    data_nrows = diablocks_pd.loc[block_id,'data_ends'] - diablocks_pd.loc[block_id,'data_starts']
    data_pd = pd.read_csv(filename,skiprows=diablocks_pd.loc[block_id,'data_starts'],nrows=data_nrows, header=None)
    data_pdser = data_pd[0].str.strip()
    data = data_pdser.str.cat()
    data = data.strip(':') #remove first and/or last colon if present
    data = data.split(':')
    
    if len(times_fromfile) != len(data):
        raise Exception(f'ERROR: times and values for block_id={block_id} are not of equal length\nlen(times_fromfile): %d\nlen(data): %d'%(len(times_fromfile),len(data)))
    
    # construct pandas dataframe (has equidistant index, so has freq property)
    data_pd = pd.DataFrame({'valuecm/qualitycode':data},index=times_fromfile)
    data_pd.index.name = 'times'
    
    # convert HWLW+quality code to separate columns
    data_pd_temp = data_pd.loc[:,'valuecm/qualitycode'].str.split('/', expand=True)
    data_pd['values'] = data_pd_temp.iloc[:,0].astype('int')
    data_pd['qualitycode'] = data_pd_temp.iloc[:,1].astype('int')
    data_pd = data_pd.drop('valuecm/qualitycode',axis='columns')
    
    # replace missing values with nan
    bool_hiaat = data_pd['qualitycode'] == 99
    data_pd.loc[bool_hiaat,'values'] = np.nan
    
    # add metadata
    metadata = metadata_from_diablocks(diablocks_pd, block_id)
    data_pd = metadata_add_to_obj(data_pd,metadata)
    return data_pd


def read_dia(filename, station=None, block_ids=None, allow_duplicates=False):
    """
    Reads an equidistant or non-equidistant dia file, or a list of dia files. Also works for diafiles containing multiple blocks for one station.

    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.
    station : TYPE
        DESCRIPTION. The default is None.
    block_ids : int, list of int or 'allstation', optional
        DESCRIPTION. The default is None.

    Returns
    -------
    data_pd : pandas.core.frame.DataFrame
        DataFrame with a 'values' column and a pd.DatetimeIndex as index in case of an equidistant file, or more columns in case of a non-equidistant file.

    """
    
    if not isinstance(filename,list):
        # solve wildcards and convert to list
        filename = glob.glob(filename)
        filename.sort()
    
    if len(filename)==0:
        raise FileNotFoundError('Filename list is empty')
    
    data_pd_list = []
    metadata_list = []
    for filename_one in filename:    
        diablocks_pd = get_diablocks(filename_one)
        pd.set_option('display.max_columns', 6) #default was 0, but need more to display groepering
        pd.set_option('display.width', 200) #default was 80, but need more to display groepering
        print_cols = ['block_starts', 'station', 'grootheid', 'groepering', 'tstart', 'tstop']
        logger.info('blocks in diafile:\n%s'%(diablocks_pd[print_cols]))
        
        #get equidistant timeseries from metadata
        if block_ids is None or block_ids=='allstation':
            if station is None:
                if len(diablocks_pd)==1:
                    station = diablocks_pd.loc[0,'station']
                else:
                    raise ValueError('If block_ids=None or block_ids="allstation", station argument should be provided. '
                                      f'Available blocks:\n{diablocks_pd[print_cols]}')
            bool_station = diablocks_pd['station']==station
            ids_station = diablocks_pd[bool_station].index.tolist()
            if len(ids_station)<1:
                raise ValueError(f"No data block with requested station ({station}) present in dia file. "
                                 f"Available blocks:\n{diablocks_pd[print_cols]}")
            elif len(ids_station)>1 and block_ids is None:
                raise ValueError(f"More than one data block with requested station ({station}) "
                                 "present in dia file. Provide block_ids argument to read_dia() (int, list of int or 'allstation'). "
                                 f"Available blocks:\n{diablocks_pd[print_cols]}")
            else: #exactly one occurrence or block_ids is provided or block_ids='allstation'
                block_ids_one = ids_station
        elif isinstance(block_ids,int):
            block_ids_one = [block_ids]
        else:
            # prevent overwriting of block_ids in this file loop
            block_ids_one = block_ids
        
        #check validity of blockids of type listlist
        if not isinstance(block_ids_one,list):
            raise TypeError('Invalid type for block_ids (should be int, list of int or "allstation")')
        if not pd.Series(block_ids_one).isin(diablocks_pd.index).all():
            raise ValueError(f"Invalid values in block_ids list ({block_ids_one}), "
                             f"possible are {diablocks_pd.index.tolist()} (all integers)")
        
        if station is not None:
            if not isinstance(station,str):
                raise TypeError('Station argument should be of type string')
            bool_samestation = diablocks_pd.loc[block_ids_one,'station']==station
            if not bool_samestation.all():
                raise ValueError("Both the arguments station and block_ids are provided, "
                                 "but at least one of the requested block_ids corresponds to a different station. "
                                 f"Available blocks:\n{diablocks_pd[print_cols]}")
            
        for block_id in block_ids_one:
            if np.isnan(diablocks_pd.loc[block_id,'timestep_min']): # non-equidistant
                data_pd_oneblock = read_dia_nonequidistant(filename_one, diablocks_pd, block_id)
            else: # equidistant
                data_pd_oneblock = read_dia_equidistant(filename_one, diablocks_pd, block_id)
            
            # add status columnm
            # first set dtype via empty column to avoid FutureWarning
            data_pd_oneblock['status'] = ''
            block_status_list = diablocks_pd.loc[block_id,'STA'].split('!')
            for block_status_one in block_status_list:
                status_tstart = pd.to_datetime(block_status_one[4:17],format='%Y%m%d;%H%M').tz_localize("UTC+01:00")
                status_tstop = pd.to_datetime(block_status_one[18:31],format='%Y%m%d;%H%M').tz_localize("UTC+01:00")
                status_val = block_status_one[-1]
                data_pd_oneblock.loc[status_tstart:status_tstop,'status'] = status_val
            
            data_pd_list.append(data_pd_oneblock)
            metadata = metadata_from_obj(data_pd_oneblock)
            metadata_list.append(metadata)
    
    #concat allyears dataset
    data_pd_all = pd.concat(data_pd_list)
    metadata_compare(metadata_list)
    metadata = metadata_list[0].copy()
    
    # convert cm to m
    assert metadata['eenheid'] == 'cm'
    data_pd_all['values'] /= 100
    metadata['eenheid'] = 'm'

    data_pd_all = metadata_add_to_obj(data_pd_all,metadata)
    
    if allow_duplicates:
        return data_pd_all
    
    #check overlapping timesteps, sort values on time
    if data_pd_all.index.duplicated().any():
        raise ValueError("Merged datasets have duplicate/overlapping timesteps, "
                         "clean up your input data or provide one file instead of a list. "
                         "Or pass `allow_duplicates=True`")
    if not data_pd_all.index.is_monotonic_increasing:
        data_pd_all = data_pd_all.sort_index()
    
    return data_pd_all


def read_noos(filename, datetime_format='%Y%m%d%H%M', na_values=None):
    """
    Reads a noos file

    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.
    datetime_format : TYPE, optional
        DESCRIPTION. The default is '%Y%m%d%H%M'.
    na_values : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    data_pd : TYPE
        DESCRIPTION.

    """
    
    logger.info('reading file: %s'%(filename))
    noosheader = []
    noosheader_dict = {} #TODO: this is not returned, could be valuable to do so
    with open(filename) as f:
        for linenum, line in enumerate(f, 0):
            if '#' in line:
                noosheader.append(line)
                comment_stripped = line.strip('#').strip().split(': ')
                if len(comment_stripped) == 1:
                    if comment_stripped[0] != '':
                        noosheader_dict[comment_stripped[0]] = ''
                else:
                    noosheader_dict[comment_stripped[0].strip()] = comment_stripped[1].strip()
            else:
                startdata = linenum
                break
        
    content_pd = pd.read_csv(filename,header=startdata-1, sep="\\s+",names=['times_str','values'],
                             na_values=na_values, dtype = {'times_str': str, 'values' : 'float'})
    noos_datetime = pd.to_datetime(content_pd['times_str'],format=datetime_format)
    data_pd = pd.DataFrame({'values':content_pd['values'].values},index=noos_datetime)
    
    # clean noos metadata and add as attrs
    noosheader_dict_clean = {k:v for k,v in noosheader_dict.items() if v!=''}
    data_pd.attrs = noosheader_dict_clean
    
    return data_pd

