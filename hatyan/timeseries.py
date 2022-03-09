# -*- coding: utf-8 -*-
"""
timeseries.py contains all definitions related to hatyan timeseries.

hatyan is a Python program for tidal analysis and prediction, based on the FORTRAN version. 
Copyright (C) 2019-2021 Rijkswaterstaat.  Maintained by Deltares, contact: Jelmer Veenstra (jelmer.veenstra@deltares.nl). 
Source code available at: https://github.com/Deltares/hatyan

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""

import os
import io
import numpy as np
import pandas as pd
import datetime as dt
import scipy.signal as ssig
file_path = os.path.realpath(__file__)
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq
from netCDF4 import Dataset, date2num, stringtoarr#, num2date
from hatyan.schureman import get_schureman_freqs #TODO: not generic
from hatyan.hatyan_core import get_const_list_hatyan


def calc_HWLW(ts, calc_HWLW345=False, calc_HWLW1122=False, debug=False):
    """
    
    Calculates extremes (high and low waters) for the provided timeseries. 
    This definition uses scipy.signal.find_peaks() with arguments 'distance' and 'prominence'. 
    The minimal 'distance' between two high or low water peaks is based on the M2 period: 12.42/1.5=8.28 hours for HW and 12.42/1.7=7.30 hours for LW (larger because of aggers). 
    The prominence for local extremes is set to 0.01m, to filter out very minor dips in the timeseries.
    If there are two equal high or low water values, the first one is taken. 
    There are no main high/low waters calculated within 6 hours of the start/end of the timeseries (keyword buffer_hr), since these can be invalid.
    This function can deal with gaps. Since scipy.signal.find_peaks() warns about nan values, those are removed first.
    This does influence the results since find_peaks does not know about time registration. This is also tricky for input timeseries with varying time interval.
    
    Parameters
    ----------
    ts : pandas.DataFrame
        The DataFrame should contain a 'values' column and a pd.DatetimeIndex as index, it contains the timeseries with a tidal prediction or water level measurements.
    calc_HWLW345 : boolean, optional
        Whether to also calculate local extremes, first/second low waters and 'aggers'. 
        The default is False, in which case only extremes per tidal period are calculated.
        When first/second low waters and aggers are calculated, the local extremes around highwater (eg double highwaters and dips) are filtered out first.
    calc_HWLW345_cleanup1122 : boolean, optional
        Whether to remove HWLWcodes 11 and 22 from DataFrame. The default is True.
    debug : boolean, optional
        Whether to print debug information. The default is False.
    
    Raises
    ------
    Exception
        DESCRIPTION.
    
    Returns
    -------
    data_pd_HWLW : pandas.DataFrame
        The DataFrame contains colums 'times', 'values' and 'HWLWcode', it contains the times, values and codes of the timeseries that are extremes.
        1 (high water) and 2 (low water). And if calc_HWLW345=True also 3 (first low water), 4 (agger) and 5 (second low water).

    """

    #calculate the amount of steps in a M2 period, based on the most occurring timestep 
    M2_period_min = get_schureman_freqs(['M2']).loc['M2','period [hr]']*60
    ts_steps_min_most = np.argmax(np.bincount((ts.index.to_series().diff().iloc[1:].dt.total_seconds()/60).astype(int).values))
    if ts_steps_min_most > 1:
        print('WARNING: the timestep of the series for which to calculate extremes/HWLW is %i minutes, but 1 minute is recommended'%(ts_steps_min_most))
    M2period_numsteps = M2_period_min/ts_steps_min_most
    
    data_pd_HWLW = pd.DataFrame({'times':ts.index,'values':ts['values'],'HWLWcode':np.nan}).reset_index(drop=True)
    #create empty HWLW dataframe
    if data_pd_HWLW['values'].isnull().any():
        data_pd_HWLW = data_pd_HWLW[~data_pd_HWLW['values'].isnull()].reset_index(drop=True)
        print('WARNING: the provided ts for extreme/HWLW calculation contained NaN values. To avoid unexpected results from scipy.signal.find_peaks(), the %i NaN values were removed from the ts (%.2f%%) before calculating extremes/HWLW.'%(len(ts)-len(data_pd_HWLW), (len(ts)-len(data_pd_HWLW))/len(ts)*100))

    if calc_HWLW345 or calc_HWLW1122:
        #get all local extremes, including aggers and second high waters (1/2/11/22) #takes first value of two equal peaks, prominence naar 0.01 om matige aggers uit te sluiten
        LWid_all, LWid_all_properties = ssig.find_peaks(-data_pd_HWLW['values'].values, prominence=(0.01,None), width=(None,None), distance=None)
        HWid_all, HWid_all_properties = ssig.find_peaks(data_pd_HWLW['values'].values, prominence=(0.01,None), width=(None,None), distance=None)
        data_pd_HWLW.loc[LWid_all,'HWLWcode'] = 22 #all LW
        data_pd_HWLW.loc[HWid_all,'HWLWcode'] = 11 #all HW

    #get HWLW (extremes per tidal period). 
    LWid_main_raw,LWid_main_properties = ssig.find_peaks(-data_pd_HWLW['values'].values, prominence=(0.01,None), width=(None,None), distance=M2period_numsteps/1.7) #most stations work with factor 1.4. 1.5 results in all LW values for HoekvanHolland for 2000, 1.7 results in all LW values for Rotterdam for 2000 (also for 1999-2002).
    HWid_main_raw,HWid_main_properties = ssig.find_peaks(data_pd_HWLW['values'].values, prominence=(0.01,None), width=(None,None), distance=M2period_numsteps/1.9) #most stations work with factor 1.4. 1.5 value results in all HW values for DenHelder for year 2000 (also for 1999-2002). 1.7 results in all HW values for LITHDP 2018. 1.9 results in all correct values for LITHDP 2022
    # remove main extremes within 6 hours of start/end of timeseries, since they are often missed or invalid.
    validtimes_idx = data_pd_HWLW.loc[(data_pd_HWLW['times']>=ts.index[0]+dt.timedelta(hours=6)) & (data_pd_HWLW['times']<=ts.index[-1]-dt.timedelta(hours=6))].index
    LWid_main = LWid_main_raw[np.in1d(LWid_main_raw,validtimes_idx)]
    HWid_main = HWid_main_raw[np.in1d(HWid_main_raw,validtimes_idx)]
    #use valid values to continue process
    data_pd_HWLW.loc[LWid_main,'HWLWcode'] = 2
    data_pd_HWLW.loc[HWid_main,'HWLWcode'] = 1
    
    #drop all non-(local)extreme timesteps and convert HWLWcode column to integers
    data_pd_HWLW = data_pd_HWLW.dropna(subset=['HWLWcode'])
    data_pd_HWLW['HWLWcode'] = data_pd_HWLW['HWLWcode'].astype(int)

    if debug: #debug statistics
        prop_list = ['prominences','widths']
        data_pd_HWLW.loc[data_pd_HWLW['HWLWcode']==2,prop_list] = pd.DataFrame(LWid_main_properties,index=LWid_main_raw).loc[LWid_main,prop_list]
        print('LW values:\n%s\n'%(data_pd_HWLW[data_pd_HWLW['HWLWcode']==2]))
        data_pd_HWLW.loc[data_pd_HWLW['HWLWcode']==1,prop_list] = pd.DataFrame(HWid_main_properties,index=HWid_main_raw).loc[HWid_main,prop_list]
        print('HW values:\n%s\n'%(data_pd_HWLW[data_pd_HWLW['HWLWcode']==1]))
        if calc_HWLW345 or calc_HWLW1122:
            LW_local_bool = ~np.in1d(LWid_all, LWid_main)
            data_pd_HWLW.loc[data_pd_HWLW['HWLWcode']==22,prop_list] = pd.DataFrame(LWid_all_properties,index=LWid_all).loc[LW_local_bool,prop_list]
            print('LW_local values:\n%s\n'%(data_pd_HWLW[data_pd_HWLW['HWLWcode']==22]))
            HW_local_bool = ~np.in1d(HWid_all, HWid_main)
            data_pd_HWLW.loc[data_pd_HWLW['HWLWcode']==11,prop_list] = pd.DataFrame(HWid_all_properties,index=HWid_all).loc[HW_local_bool,prop_list]
            print('HW_local values:\n%s\n'%(data_pd_HWLW[data_pd_HWLW['HWLWcode']==11]))
    
    if calc_HWLW345: #recalculate local LW/HWs between two main HWs to firstLW/agger/secondLW
        data_pd_HWLW = calc_HWLWlocalto345(data_pd_HWLW,HWid_main)
    
    #return to normal time-index
    data_pd_HWLW = data_pd_HWLW.set_index('times')
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
    
    print('calculating 1stLW/agger/2ndLW for all tidalperiods (between two HW values)...')
    for iTide, dummy in enumerate(HWid_main[:-1]):
        data_pd_HWLW_1tide = data_pd_HWLW.loc[HWid_main[iTide]:HWid_main[iTide+1],:]
        
        #filter local extremes around HW (only interested in aggers, so LW), this is necessary for eg DENHDR and PETTZD, otherwise second HW is seen as first LW
        data_pd_HWLW_1tide_minHW = data_pd_HWLW_1tide.loc[data_pd_HWLW_1tide['HWLWcode']==1,['values']].min()[0]
        data_pd_HWLW_1tide_min = data_pd_HWLW_1tide['values'].min()
        data_pd_HWLW_1tide_mid = np.mean([data_pd_HWLW_1tide_minHW,data_pd_HWLW_1tide_min])
        bool_LWs = data_pd_HWLW_1tide['values']<data_pd_HWLW_1tide_mid
        data_pd_HWLW_1tide_noHWs = data_pd_HWLW_1tide[bool_LWs]
        
        if len(data_pd_HWLW_1tide_noHWs) > 3: #(attempt to) reduce to three values between two HWs
            print('WARNING: more than 3 values between HWs, removing part of them')
            #print(data_pd_HWLW_1tide)
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
            print('WARNING: 2 values left between two HWs (slightly unexpected):\n%s'%(data_pd_HWLW_1tide))
        else:
            raise Exception('unexpected number of values between two HWs (0 or more than 3):\n%s'%(data_pd_HWLW_1tide))
    
    #remove remaining 11 and 22 values from array
    #if calc_HWLW345_cleanup1122:
    data_pd_HWLW = data_pd_HWLW.drop(data_pd_HWLW[data_pd_HWLW['HWLWcode']==11].index)
    data_pd_HWLW = data_pd_HWLW.drop(data_pd_HWLW[data_pd_HWLW['HWLWcode']==22].index)
    
    print('finished calculating 1stLW/agger/2ndLW for all tidalperiods')
    
    return data_pd_HWLW


def calc_HWLWnumbering(ts_ext, station=None, corr_tideperiods=None):
    """
    For calculation of the extremes numbering, w.r.t. the first high water at Cadzand in 2000 (occurred on 1-1-2000 at approximately 9:45). 
    The number of every high and low water is calculated by taking the time difference between itself and the first high water at Cadzand, correcting it with the station phase difference (M2phasediff). 
    Low waters are searched for half an M2 period from the high waters. 
    By adding a search window of half the period of M2 (searchwindow_hr), even strong time variance between consecutive high or low waters should be caputered. 
    
    Parameters
    ----------
    ts_ext : pandas.DataFrame
        The DataFrame should contain a 'values' and 'HWLWcode' column and a pd.DatetimeIndex as index, it contains the times, values and codes of the timeseries that are extremes.
    station: string, optional
        The station for which the M2 phase difference should be retrieved from data_M2phasediff_perstation.txt.
        This value is the phase difference in degrees of the occurrence of the high water generated by the same tidal wave as the first high water in 2000 at Cadzand (actually difference between M2 phases of stations).
        This value is used to correct the search window of high/low water numbering. The default is None.
    corr_tideperiods : integer, optional
        Test keyword to derive HWLWnumbering with a n*360 degrees offset only, but this does not work properly. The default is None.

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
    firstHWcadz_fixed = dt.datetime(2000, 1, 1, 9, 45)
    searchwindow_hr = M2_period_hr/2
    
    if not all((ts_ext['HWLWcode']==1) | (ts_ext['HWLWcode']==2) | (ts_ext['HWLWcode']==3) | (ts_ext['HWLWcode']==4) | (ts_ext['HWLWcode']==5)):
        raise Exception('calc_HWLWnumbering() not implemented for HWLWcode other than 1,2,3,4,5 (so no HWLWcode 11 or 22 supported), provide extreme timeseries derived with Timeseries.calc_HWLW(calc_HWLW345=False) or Timeseries.calc_HWLW(calc_HWLW345=True, calc_HWLW345_cleanup1122=True)')
    ts_ext = ts_ext.copy()
    
    HW_bool = ts_ext['HWLWcode']==1
    HW_tdiff_cadzdraw = (ts_ext.loc[HW_bool].index.to_series()-firstHWcadz_fixed).dt.total_seconds()/3600
    if station is None:
        HW_tdiff_cadzdraw_M2remainders = (HW_tdiff_cadzdraw)%M2_period_hr
        M2phasediff_hr = (HW_tdiff_cadzdraw_M2remainders).mean()
        M2phasediff_deg = M2phasediff_hr/M2_period_hr*360
        print('no value or None for argument M2phasediff provided, automatically calculated correction w.r.t. Cadzand is %.2f hours (%.2f degrees)'%(M2phasediff_hr, M2phasediff_deg))
        if corr_tideperiods is not None:
            M2phasediff_deg = M2phasediff_deg+corr_tideperiods
            M2phasediff_hr = M2phasediff_deg/360*M2_period_hr
            print('additional tideperiod correction provided via corr_tideperiods of %.1f degrees, new correction w.r.t. Cadzand is %.2f hours (%.2f degrees)'%(corr_tideperiods, M2phasediff_hr, M2phasediff_deg))
    else:
        file_M2phasediff = os.path.join(os.path.dirname(file_path),'data','data_M2phasediff_perstation.txt')
        stations_M2phasediff = pd.read_csv(file_M2phasediff, names=['M2phasediff'], comment='#', delim_whitespace=True)
        stat_M2phasediff = stations_M2phasediff.loc[station,'M2phasediff']
        M2phasediff_hr = stat_M2phasediff/360*M2_period_hr
    HW_tdiff_cadzd = HW_tdiff_cadzdraw - M2phasediff_hr + searchwindow_hr
    HW_tdiff_div, HW_tdiff_mod_searchwindow = np.divmod(HW_tdiff_cadzd.values, M2_period_hr)
    HW_tdiff_mod = HW_tdiff_mod_searchwindow - searchwindow_hr
    if not all(np.diff(HW_tdiff_div) > 0):
        raise Exception('tidal wave numbering: HW numbers not always increasing')
    if not all(np.abs(HW_tdiff_mod)<searchwindow_hr):
        raise Exception('tidal wave numbering: not all HW fall into hardcoded search window')
    ts_ext.loc[HW_bool,'HWLWno'] = HW_tdiff_div
    
    for LWcode_2345 in [2,3,4,5]:
        LW_bool = ts_ext['HWLWcode']==LWcode_2345
        LW_tdiff_cadzdraw = (ts_ext.loc[LW_bool].index.to_series()-firstHWcadz_fixed).dt.total_seconds()/3600
        LW_tdiff_cadzd = LW_tdiff_cadzdraw - M2phasediff_hr + searchwindow_hr - M2_period_hr/2
        LW_tdiff_div, LW_tdiff_mod_searchwindow = np.divmod(LW_tdiff_cadzd.values, M2_period_hr)
        LW_tdiff_mod = LW_tdiff_mod_searchwindow - searchwindow_hr
        if not all(np.diff(LW_tdiff_div) > 0):
            raise Exception('tidal wave numbering: LW numbers not always increasing')
        if not all(np.abs(LW_tdiff_mod)<searchwindow_hr):
            raise Exception('tidal wave numbering: not all LW fall into defined search window')
        ts_ext.loc[LW_bool,'HWLWno'] = LW_tdiff_div
    
    #check if LW is after HW
    ts_ext_checkfirst = ts_ext[ts_ext['HWLWno']==np.min(HW_tdiff_div)]
    tdiff_firstHWLW = (ts_ext_checkfirst.index.to_series().diff().dt.total_seconds()/3600).values[1]
    if (tdiff_firstHWLW<0) or (tdiff_firstHWLW>M2_period_hr):
        raise Exception('tidal wave numbering: first LW does not match first HW')
    
    ts_ext['HWLWno'] = ts_ext['HWLWno'].astype(int)
    
    return ts_ext


def timeseries_fft(ts_residue, prominence=10**3, plot_fft=True):
    
    print('analyzing timeseries with fft and fftfreq')
    
    y = ts_residue['values'].values
    N = len(y)
    T = np.unique((ts_residue.index[1:]-ts_residue.index[:-1])).astype(float)/1e9/3600 #timestep in hours.
    if len(T)!=1:
        raise Exception('timestep of supplied timeseries should be constant for fourier analysis')
    yf = fft(y)
    power = np.abs(yf)
    freq = fftfreq(N, T[0])
    peaks = ssig.find_peaks(power[freq >=0], prominence=prominence)[0]
    peak_freq =  freq[peaks]
    peak_power = power[peaks]
    
    if plot_fft:
        fig,ax = plt.subplots()
        ax.plot(freq[:N//2], power[:N//2])
        ax.plot(peak_freq, peak_power, 'ro')
        ax.grid()
        ax.set_xlim(0,0.5)
    
    const_list_all = get_const_list_hatyan(listtype='all_schureman') #TODO: not generic
    hatyan_freqs = get_schureman_freqs(const_list=const_list_all)[['freq']]
    const_match = []
    const_closest = []
    for peak_freq_one in peak_freq:
        hatyan_freqs_match = hatyan_freqs[np.abs(hatyan_freqs['freq']-peak_freq_one)<4e-5]
        #print(peak_freq_one)
        #print(hatyan_freqs_match)
        hatyan_freqs_match_list = [x for x in hatyan_freqs_match.index if '_IHO' not in x]
        const_match = const_match+hatyan_freqs_match_list
        hatyan_freqs_closest = hatyan_freqs.iloc[np.argmin(np.abs(hatyan_freqs-peak_freq_one)),:]
        const_closest.append(hatyan_freqs_closest.name)
    hatyan_freqs_matches = get_schureman_freqs(const_list=const_match)[['freq','period [hr]']]
    hatyan_freqs_suggestions = get_schureman_freqs(const_list=const_closest)[['freq','period [hr]']]
    hatyan_freqs_suggestions['peak_freq'] = peak_freq
    hatyan_freqs_suggestions['peak_power'] = peak_power
    print('dominant freqs from fft:\n%s'%(peak_freq))
    print('suggested constituents+freqs from hatyan:\n%s'%(hatyan_freqs_suggestions))
    
    return peak_freq, hatyan_freqs_suggestions, hatyan_freqs_matches


def plot_timeseries(ts, ts_validation=None, ts_ext=None, ts_ext_validation=None):
    """
    Creates a plot with the provided timeseries

    Parameters
    ----------
    ts : pandas.DataFrame
        The DataFrame should contain a 'values' column and a pd.DatetimeIndex as index, it contains the timeseries.
    ts_validation : pandas.DataFrame, optional
        The DataFrame should contain a 'values' column and a pd.DatetimeIndex as index, it contains the timeseries. The default is None.
    ts_ext : pandas.DataFrame, optional
        The DataFrame should contain a 'values' and 'HWLW_code' column and a pd.DatetimeIndex as index, it contains the times, values and codes of the timeseries that are extremes. The default is None.
    ts_ext_validation : pandas.DataFrame, optional
        The DataFrame should contain a 'values' and 'HWLW_code' column and a pd.DatetimeIndex as index, it contains the times, values and codes of the timeseries that are extremes. The default is None.

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
    
    if ts_validation is not None:
        times_predval_ext = [min(min(ts_validation.index),min(ts.index)), max(max(ts_validation.index),max(ts.index))]

    else:
        times_predval_ext = [min(ts.index), max(ts.index)]    

    fig, (ax1, ax2) = plt.subplots(2,1,figsize=size_figure, sharex=True, gridspec_kw={'height_ratios':[2,1]})
    
    ax1.set_title('hatyan timeseries')
    ax1.plot(ts.index, ts['values'],'o-',linewidth=size_line_ts,markersize=size_marker_ts, label='ts')
    if ts_validation is not None:
        #overlap between timeseries for difference plots
        times_id_validationinpred = np.where(ts_validation.index.isin(ts.index))[0]
        times_id_predinvalidation = np.where(ts.index.isin(ts_validation.index))[0]
        ax1.plot(ts_validation.index, ts_validation['values'],'o-',linewidth=size_line_ts,markersize=size_marker_ts, label='ts_validation', alpha=0.7)
        ax1.plot(ts.index[times_id_predinvalidation], ts['values'].iloc[times_id_predinvalidation].values-ts_validation['values'].iloc[times_id_validationinpred].values,'go-',linewidth=size_line_ts,markersize=size_marker_ts, label='difference', alpha=0.7)
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
        #print HWLW statistics
        try:
            plot_HWLW_validatestats(ts_ext=ts_ext, ts_ext_validation=ts_ext_validation, create_plot=False)        
        except:
            print('WARNING: plot_HWLW_validatestats() failed, probably due to missing HWLWno where autocalculation failed. Consider adding HWLWno to ts_ext and ts_ext_validation with calc_HWLWnumbering() before plotting.')
    ax1.set_ylim(figure_ylim_ts)
    ax2.set_xlabel('Time')
    ax1.set_ylabel('waterlevel [m]')
    ax1.legend(loc='lower right')
    ax1.grid()
    if ts_validation is not None:
        ax2.plot(ts.index[times_id_predinvalidation], ts['values'].iloc[times_id_predinvalidation].values-ts_validation['values'].iloc[times_id_validationinpred].values,'go-',linewidth=size_line_ts,markersize=size_marker_ts, label='difference')
    ax2.plot(times_predval_ext,[0,0],'-k',linewidth=size_line_ts)
    ax2.set_ylim(figure_ylim_tsdiff)
    rmse = np.nan
    if ts_validation is not None:
        overlapdiff = ts['values'].iloc[times_id_predinvalidation].values-ts_validation['values'].iloc[times_id_validationinpred].values
        if len(overlapdiff) != 0:
            rmse = np.sqrt(np.nanmean(overlapdiff ** 2))
    ax2.set_ylabel('timeseries difference [m], RMSE = %.5f'%(rmse))
    ax2.legend(loc='lower right')
    ax2.grid()
    fig.tight_layout()
    
    axs = (ax1,ax2)
    return fig, axs


def plot_HWLW_validatestats(ts_ext, ts_ext_validation, create_plot=True):
    """
    This definition calculates (and plots and prints) some statistics when comparing extreme values.
    This is done by calculating the extreme number (sort of relative to Cadzand 1jan2000, but see 'warning') and subtracting the ts_ext and ts_ext_validation dataframes based on these numbers (and HWLWcode).
    It will only result in values for the overlapping extremes, other values will be NaN and are not considered for the statistics.
    Warning: the calculated extreme numbers in this definition are not corrected for the real phase difference with the M2phasediff argument, the calculated extreme are fine for internal use (to match corresponding extremes) but the absolute number might be incorrect.

    Parameters
    ----------
    ts_ext : pandas.DataFrame
        The DataFrame should contain a 'values' and 'HWLW_code' column and a pd.DatetimeIndex as index, it contains the times, values and codes of the timeseries that are extremes.
    ts_ext_validation : pandas.DataFrame
        The DataFrame should contain a 'values' and 'HWLW_code' column and a pd.DatetimeIndex as index, values and codes of the timeseries that are extremes.
    create_plot : boolean, optional
        Whether to plot the time/value differences or only print the statistics. The default is True.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The generated figure handle, with which the figure can be adapted and saved.
    axs : (tuple of) matplotlib.axes._subplots.AxesSubplot
        The generated axis handle, whith which the figure can be adapted.

    """
    
    print('Calculating comparison statistics for extremes')
    if 'HWLWno' not in ts_ext.columns or 'HWLWno' not in ts_ext_validation.columns:
        print('HWLWno is not present in ts_ext or ts_ext_validation, trying to automatically derive it without M2phasediff argument (this might fail)')
        try:
            ts_ext_nrs = calc_HWLWnumbering(ts_ext=ts_ext)
            ts_ext_validation_nrs = calc_HWLWnumbering(ts_ext=ts_ext_validation)
        except:
            raise Exception('ERROR: deriving HWLWno failed, so HWLW statistics cannot be calculated. Add HWLWno with calc_HWLWnumbering() before calling plot_HWLW_validatestats().')
    else:
        ts_ext_nrs = ts_ext.copy()
        ts_ext_validation_nrs = ts_ext_validation.copy()
    
    #set HWLWcode and HWLWno as index, to make easy subtraction possible
    ts_ext_nrs['times'] = ts_ext_nrs.index
    ts_ext_nrs = ts_ext_nrs.set_index(['HWLWcode','HWLWno'],drop=False)
    ts_ext_validation_nrs['times'] = ts_ext_validation_nrs.index
    ts_ext_validation_nrs = ts_ext_validation_nrs.set_index(['HWLWcode','HWLWno'],drop=False)
    HWLW_diff = ts_ext_nrs.sub(ts_ext_validation_nrs)
    
    tdiff_minutes = HWLW_diff['times'].dt.total_seconds()/60
    vdiff_cm = HWLW_diff['values']*100
    print('Time differences [minutes]')
    print('    RMSE: %.2f'%(np.sqrt(np.mean(tdiff_minutes**2))))
    print('    std: %.2f'%(tdiff_minutes.std()))
    print('    abs max: %.2f'%(tdiff_minutes.abs().max()))
    print('    abs mean: %.2f'%(tdiff_minutes.abs().mean()))
    print('    #NaN: %i of %i'%(tdiff_minutes.isnull().sum(),len(vdiff_cm)))
    print('Value differences [cm]')
    print('    RMSE: %.2f'%(np.sqrt(np.mean(vdiff_cm**2))))
    print('    std: %.2f'%(vdiff_cm.std()))
    print('    abs max: %.2f'%(vdiff_cm.abs().max()))
    print('    abs mean: %.2f'%(vdiff_cm.abs().mean()))
    print('    #NaN: %i of %i'%(vdiff_cm.isnull().sum(),len(vdiff_cm)))
    
    if create_plot:
        fig, ax1 = plt.subplots()
        ax1.plot(HWLW_diff.loc[1,'times'].dt.total_seconds()/60,HWLW_diff.loc[1,'values']*100,'+',label='HWdiff')
        ax1.plot(HWLW_diff.loc[2,'times'].dt.total_seconds()/60,HWLW_diff.loc[2,'values']*100,'.',label='LWdiff')
        ax1.set_xlabel('Time difference [minutes]')
        ax1.set_ylabel('Value difference [cm]')
        ax1.legend(loc=1)
        ax1.grid()
    
        axs = (ax1)
        return fig, axs


def write_tsnetcdf(ts, station, vertref, filename, ts_ext=None, tzone_hr=1, nosidx=False, mode='w'):
    """
    Writes the timeseries to a netCDF file

    Parameters
    ----------
    ts : pandas.DataFrame
        The DataFrame should contain a 'values' column and a pd.DatetimeIndex as index, it contains the timeseries.
    station : str
        DESCRIPTION.
    vertref : str
        DESCRIPTION.
    filename : str
        The filename of the netCDF file that will be written.
    ts_ext : pandas.DataFrame, optional
        The DataFrame should contain a 'values' and 'HWLW_code' column and a pd.DatetimeIndex as index, it contains the times, values and codes of the timeseries that are extremes. The default is None.
    tzone_hr : int, optional
        The timezone (GMT+tzone_hr) that applies to the data. The default is 1 (MET).
    
    Returns
    -------
    None.
    
    """
    
    import hatyan
    version_no = hatyan.__version__
    
    times_all = ts.index
    timeseries = ts['values']
    times_stepmin = (ts.index[1]-ts.index[0]).total_seconds()/60
    dt_analysistime = dt.datetime.now()
    data_nc = Dataset(filename, mode, format="NETCDF3_CLASSIC")
    attr_dict = {'title': 'tidal prediction for %s to %s'%(times_all[0].strftime('%Y-%m-%d %H:%M:%S'), times_all[-1].strftime('%Y-%m-%d %H:%M:%S')),
                 'institution': 'Rijkswaterstaat',
                 'source': 'hatyan-%s tidal analysis program of Rijkswaterstaat'%(version_no),
                 'timestep_min': times_stepmin}
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
    
    refdate_tz = dt.datetime(1900,1,1,tzinfo=dt.timezone(dt.timedelta(hours=tzone_hr)))
    dict_statattr = {'cf_role': 'timeseries_id'}
    dict_anatimattr = {'units': 'minutes since %s'%(refdate_tz.strftime('%Y-%m-%d %H:%M:%S %z')), 'standard_name':'forecast_reference_time', 'long_name':'forecast_reference_time'}
    dict_timattr = {'units': 'minutes since %s'%(refdate_tz.strftime('%Y-%m-%d %H:%M:%S %z'))}
    dict_wlattr = {'units':'m', 'vertical_reference': vertref, 'standard_name': 'sea_surface_height_above_geopotential_datum', 'long_name': 'astronomical prediction of water level above reference level'}
    dict_HWattr = {'units':'m', 'vertical_reference': vertref, 'standard_name': 'sea_surface_height_above_geopotential_datum', 'long_name': 'astronomical prediction of high water extremes above reference level'}
    dict_LWattr = {'units':'m', 'vertical_reference': vertref, 'standard_name': 'sea_surface_height_above_geopotential_datum', 'long_name': 'astronomical prediction of low water extremes above reference level'}
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
        print('no HWLW prediction written')
        data_nc.close()
        return #this skips the HWLW part of the definition
    
    #HWLW prediction
    if nosidx:
        #convert index from time to HWLWno
        data_HWLW_nosidx = ts_ext.copy()
        data_HWLW_nosidx['times'] = data_HWLW_nosidx.index
        data_HWLW_nosidx = data_HWLW_nosidx.set_index('HWLWno')
        HWLWno_all = data_HWLW_nosidx.index.unique()
        data_HW = pd.DataFrame(data_HWLW_nosidx.loc[data_HWLW_nosidx['HWLWcode']==1],index=HWLWno_all)
        data_LW = pd.DataFrame(data_HWLW_nosidx.loc[data_HWLW_nosidx['HWLWcode']==2],index=HWLWno_all)
        
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
        data_nc.variables['times_astro_HW'][nstat,:] = date2num(data_HW['times'].tolist(),units=data_nc.variables['times_astro_HW'].units)
        if 'waterlevel_astro_HW' not in ncvarlist:
            nc_newvar = data_nc.createVariable('waterlevel_astro_HW','f8',('stations','HWLWno',))
            nc_newvar.setncatts(dict_HWattr)
        data_nc.variables['waterlevel_astro_HW'][nstat,:] = data_HW['values']
        #LW
        if 'times_astro_LW' not in ncvarlist:
            nc_newvar = data_nc.createVariable('times_astro_LW','f8',('stations','HWLWno',)) 
            nc_newvar.setncatts(dict_timattr)
        data_nc.variables['times_astro_LW'][nstat,:] = date2num(data_LW['times'].tolist(),units=data_nc.variables['times_astro_LW'].units)
        if 'waterlevel_astro_LW' not in ncvarlist:
            nc_newvar = data_nc.createVariable('waterlevel_astro_LW','f8',('stations','HWLWno',))
            nc_newvar.setncatts(dict_LWattr)
        data_nc.variables['waterlevel_astro_LW'][nstat,:] = data_LW['values']
    
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


def write_tsdia(ts, station, vertref, filename, headerformat='dia'):
    """
    Writes the timeseries to an equidistant dia file

    Parameters
    ----------
    ts : pandas.DataFrame
        The DataFrame should contain a 'values' column and a pd.DatetimeIndex as index, it contains the timeseries.
    station : TYPE
        DESCRIPTION.
    vertref : TYPE
        DESCRIPTION.
    filename : TYPE
        DESCRIPTION.

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    if vertref == 'NAP':
        waarnemingssoort = 18
        vertreflong = 'T.o.v. Normaal Amsterdams Peil'
    elif vertref == 'MSL':
        waarnemingssoort = 55
        vertreflong = 'T.o.v. Mean Sea Level'
    else:
        raise Exception('ERROR: currently only vertref="NAP" and vertref="MSL" are supported for writing diafiles')
    grootheid = 'WATHTBRKD;Waterhoogte berekend;J'
    ana = 'F012;Waterhoogte astronomisch mbv harmonische analyse'
    
    time_today = dt.datetime.today().strftime('%Y%m%d')
    tstart_str = ts.index[0].strftime('%Y%m%d;%H%M')
    tstop_str = ts.index[-1].strftime('%Y%m%d;%H%M')
    timestep_min = (ts.index[1]-ts.index[0]).total_seconds()/60
    
    ts_values = ts['values']
    metadata_pd = pd.Series(['[IDT;*DIF*;A;;%6s]'%(time_today),
                             '[W3H]',
                             'WNS;%i'%(waarnemingssoort),
                             'PAR;%s'%(grootheid), #parameter/grootheid, gelijk voor waarnemingssoorten 18 en 55
                             'CPM;10;Oppervlaktewater', #compartiment, gelijk voor waarnemingssoorten 18 en 55
                             'EHD;I;cm', #eenheid, gelijk voor waarnemingssoorten 18 en 55
                             'HDH;%s;%s'%(vertref,vertreflong),
                             ##'ORG;NVT;Niet van toepassing',
                             ##'SGK;NVT',
                             ##'IVS;NVT;Niet van toepassing',
                             ##'BTX;NVT;NVT;Niet van toepassing',
                             ##'BTN;Niet van toepassing',
                             'ANI;RIKZITSDHG;RIKZ - afdeling ZDI te Den Haag', #niet_essentieel?
                             'BHI;RIKZITSDHG;RIKZ - afdeling ZDI te Den Haag', #niet_essentieel?
                             'BMI;NVT;Niet van toepassing', #niet_essentieel?
                             'OGI;RIKZMON_WAT;RIKZ - Landelijke monitoring waterhoogten gegevens', #niet_essentieel?
                             ##'GBD;NIEUWWTWG;Nieuwe Waterweg',
                             'LOC;%s'%(station), #;Hoek van Holland;P;RD;6793000;44400000
                             'ANA;%s'%(ana),
                             'BEM;NVT', #niet_essentieel?
                             'BEW;NVT', #niet_essentieel?
                             'VAT;NVT', #niet_essentieel?
                             'TYP;TE',
                             '[RKS]',
                             'TYD;%10s;%10s;%i;min'%(tstart_str,tstop_str,timestep_min),
                             ##'PLT;NVT;-999999999;6793000;44400000',
                             ##'SYS;CENT',
                             '[TPS]',
                             'STA;%10s;%10s;O'%(tstart_str,tstop_str),
                             '[WRD]'])
    if headerformat=='wia':
        for metalinestart in ['[IDT;','WNS']:
            bool_drop = metadata_pd.str.startswith(metalinestart)
            metadata_pd = metadata_pd[~bool_drop]
        metadata_pd[metadata_pd.str.startswith('PAR')] = 'GHD;%s'%(grootheid)#.split(';')[0])
        metadata_pd[metadata_pd.str.startswith('CPM')] = 'CPM;OW;Oppervlaktewater'
        metadata_pd[metadata_pd.str.startswith('ANA')] = 'WBM;other:%s'%(ana)#.split(';')[0])
    
    linestr_list = []
    linestr = ''
    for iV, ts_value in enumerate(ts_values): # iterate over ts_values
        linestr_add = "%i/0:"%(np.round(ts_value*100))
        linestr = linestr + linestr_add
        if (len(linestr) > 114) or (iV==len(ts_values)-1): # append linestr to linestr_list if linestr is longer than n characters or last item of ts_values was reached
            linestr_list.append(linestr)
            linestr = ''
    data_todia = pd.Series(linestr_list)
    
    with io.open(filename,'w', newline='\n') as f: #open file and set linux newline style
        for metaline in metadata_pd:
            f.write('%s\n'%(metaline))
        data_todia.to_csv(f,index=False,header=False)


def write_tsdia_HWLW(ts_ext, station, vertref, filename, headerformat='dia'):
    """
    writes the extremes timeseries to a non-equidistant dia file

    Parameters
    ----------
    ts_ext : pandas.DataFrame
        The DataFrame should contain a 'values' and 'HWLW_code' column and a pd.DatetimeIndex as index, it contains the times, values and codes of the timeseries that are extremes.
    station : TYPE
        DESCRIPTION.
    vertref : TYPE
        DESCRIPTION.
    filename : TYPE
        DESCRIPTION.

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    if vertref == 'NAP':
        waarnemingssoort = 18
        vertreflong = 'T.o.v. Normaal Amsterdams Peil'
        parameterX = 'GETETBRKD2;Getijextreem berekend'
    elif vertref == 'MSL':
        waarnemingssoort = 55
        vertreflong = 'T.o.v. Mean Sea Level'
        parameterX = 'GETETBRKDMSL2;Getijextreem berekend t.o.v. MSL'
    else:
        raise Exception('ERROR: currently only vertref="NAP" and vertref="MSL" are supported for writing diafiles')
    grootheid = 'WATHTBRKD;Waterhoogte berekend;J'
    ana = 'F012;Waterhoogte astronomisch mbv harmonische analyse' #HW en LW uit 1 min. waterhoogten gefilterd uit 10 min. gem.
    time_today = dt.datetime.today().strftime('%Y%m%d')
    tstart_str = ts_ext.index[0].strftime('%Y%m%d;%H%M')
    tstop_str = ts_ext.index[-1].strftime('%Y%m%d;%H%M')
    
    if 11 in ts_ext['HWLWcode'].values or 22 in ts_ext['HWLWcode'].values:
        raise Exception('ERROR: invalid HWLWcodes in provided extreme timeseries (11 and/or 22)')
    
    metadata_pd = pd.Series(['[IDT;*DIF*;A;;%6s]'%(time_today),
                             '[W3H]',
                             'MUX;%s'%(parameterX),
                             ##IVS;NVT;Niet van toepassing
                             ##BTX;NVT;NVT;Niet van toepassing
                             ##BTN;Niet van toepassing
                             'ANI;RIKZITSDHG;RIKZ - afdeling ZDI te Den Haag', #niet_essentieel?
                             'BHI;RIKZITSDHG;RIKZ - afdeling ZDI te Den Haag', #niet_essentieel?
                             'BMI;NVT;Niet van toepassing', #niet_essentieel?
                             'OGI;RIKZMON_WAT;RIKZ - Landelijke monitoring waterhoogten gegevens', #niet_essentieel?
                             ##GBD;NIEUWWTWG;Nieuwe Waterweg
                             'LOC;%s'%(station),
                             'ANA;%s'%(ana),
                             'BEM;NVT;Niet van toepassing', #niet_essentieel?
                             'BEW;NVT;Niet van toepassing', #niet_essentieel?
                             'VAT;NVT;Niet van toepassing', #niet_essentieel?
                             'TYP;TN',
                             '[MUX]',
                             'MXW;1;15',
                             'MXP;1;GETETCDE;Getijextreem code;J',
                             'MXC;1;10;Oppervlaktewater',
                             'MXE;1;T;DIMSLS',
                             'MXH;1;NVT;Niet van toepassing', #niet_essentieel?
                             'MXO;1;NVT;Niet van toepassing', #niet_essentieel?
                             'MXS;1;NVT', #niet_essentieel?
                             'MXW;2;%i'%(waarnemingssoort),
                             'MXP;2;%s'%(grootheid),
                             'MXC;2;10;Oppervlaktewater',
                             'MXE;2;I;cm',
                             'MXH;2;%s;%s'%(vertref, vertreflong),
                             'MXO;2;NVT;Niet van toepassing', #niet_essentieel?
                             'MXS;2;NVT', #niet_essentieel?
                             '[TYP]',
                             'TVL;1;1;hoogwater',
                             'TVL;1;2;laagwater',
                             'TVL;1;3;laagwater 1',
                             'TVL;1;4;topagger',
                             'TVL;1;5;laagwater 2',
                             '[RKS]',
                             'TYD;%10s;%10s'%(tstart_str,tstop_str),
                             ##'PLT;NVT;-999999999;6793000;44400000',
                             'SYS;CENT', #niet_essentieel?
                             '[TPS]',
                             'STA;%10s;%10s;O'%(tstart_str,tstop_str),
                             '[WRD]'])
    if headerformat=='wia':
        for metalinestart in ['[IDT;','MXW']:
            bool_drop = metadata_pd.str.startswith(metalinestart)
            metadata_pd = metadata_pd[~bool_drop]
        metadata_pd[metadata_pd.str.startswith('ANA')] = 'WBM;other:%s'%(ana)#.split(';')[0])
        metadata_pd[metadata_pd.str.startswith('MXP;1')] = 'MXT;1;GETETTPE' # GETETCDE;Getijextreem code naar GETETTPE
        metadata_pd[metadata_pd.str.startswith('MXC;1')] = 'MXC;1;OW;Oppervlaktewater'
        metadata_pd[metadata_pd.str.startswith('MXP;2')] = 'MXG;2;%s'%(grootheid)
        metadata_pd[metadata_pd.str.startswith('MXC;2')] = 'MXC;2;OW;Oppervlaktewater'

    data_todia = ts_ext.index.strftime('%Y%m%d;%H%M')+';'+ts_ext['HWLWcode'].astype(str)+'/0;'+(ts_ext['values']*100).round().astype(int).astype(str)+':'

    with io.open(filename,'w', newline='\n') as f: #open file and set linux newline style
        for metaline in metadata_pd:
            f.write('%s\n'%(metaline))
        data_todia.to_csv(f,index=False,header=False)


def writets_noos(ts, filename, metadata=None):
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
        header_txt = f"""------------------------------------------------------\nTimeseries written by hatyan\nCreated at {timestamp}\n------------------------------------------------------\n"""
        for key in metadata.keys():
            header_txt = header_txt+('%-12s: %s\n'%(key, metadata[key]))
        header_txt = header_txt+'------------------------------------------------------'
    np.savetxt(filename,ts_out,fmt='%s %7.4f',header=header_txt)


def crop_timeseries(ts, times_ext, onlyfull=True):
    """
    Crops the provided timeseries

    Parameters
    ----------
    ts : pandas.DataFrame
        The DataFrame should contain a 'values' column and a pd.DatetimeIndex as index, it contains the timeseries.
    times_ext : TYPE
        DESCRIPTION.

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
    
    print('-'*50)
    print('cropping timeseries')
    if not times_ext[0]<times_ext[1]:
        raise Exception('ERROR: the two times times_ext should be increasing, but they are not: %s.'%(times_ext))
    if (times_ext[0] < ts_pd_in.index[0]) or (times_ext[-1] > ts_pd_in.index[-1]):
        message = 'imported timeseries is not available within entire requested period:\nrequested period:    %s to %s\nimported timeseries: %s to %s'%(times_ext[0],times_ext[-1],ts_pd_in.index[0],ts_pd_in.index[-1])
        if onlyfull:
            raise Exception('ERROR: %s'%(message))
        else:
            print('WARNING: %s'%(message))
            
    times_selected_bool = (ts_pd_in.index >= times_ext[0]) & (ts_pd_in.index <= times_ext[-1])
    ts_pd_out = ts_pd_in.loc[times_selected_bool]
    
    check_ts(ts_pd_out)
    return ts_pd_out


def resample_timeseries(ts, timestep_min, tstart=None, tstop=None):
    """
    resamples the provided timeseries, only overlapping timesteps are selected, so no interpolation. with tstart/tstop it is possible to extend the timeseries with NaN values.

    Parameters
    ----------
    ts : pandas.DataFrame
        The DataFrame should contain a 'values' and 'HWLW_code' column and a pd.DatetimeIndex as index, it contains the timeseries to be resampled.
    timestep_min : int
        the amount of minutes with which to resample the timeseries.
    tstart : dt.datetime, optional
        the start date for the resampled timeseries, the default is None which results in using the start date of the input ts.
    tstop : dt.datetime, optional
        the stop date for the resampled timeseries, the default is None which results in using the stop date of the input ts.

    Returns
    -------
    data_pd_resample : pandas.DataFrame with a 'values' column and a pd.DatetimeIndex as index
        the resampled timeseries.

    """
    
    print('-'*50)
    print('resampling timeseries to %i minutes'%(timestep_min))
    
    bool_duplicated_index = ts.index.duplicated()
    if bool_duplicated_index.sum()>0:
        raise Exception('there are duplicated values in the ts DatetimeIndex, this is not supported by Timeseries.resample_timeseries(). Try "ts_nodupl = ts[~ts.index.duplicated()]"')
    
    if tstart is None:
        tstart = ts.index[0]
    if tstop is None:
        tstop = ts.index[-1]
    data_pd_resample = pd.DataFrame({},index=pd.date_range(tstart,tstop,freq='%dmin'%(timestep_min))) #generate timeseries with correct tstart/tstop and interval
    data_pd_resample['values'] = ts['values'] #put measurements into this timeseries, matches to correct index automatically
    
    check_ts(data_pd_resample)
    return data_pd_resample


def check_rayleigh(ts_pd,t_const_freq_pd):
    
    t_const_freq = t_const_freq_pd['freq']
    freq_diffs = np.diff(t_const_freq)
    rayleigh_tresh = 0.99
    rayleigh = len(ts_pd['values'])*freq_diffs
    freq_diff_min = rayleigh_tresh/len(ts_pd['values'])
    rayleigh_bool = rayleigh>rayleigh_tresh
    rayleigh_bool_id = np.where(~rayleigh_bool)[0]
    
    if rayleigh_bool.all():
        print('Rayleigh criterion OK (always>%.2f, minimum is %.2f)'%(rayleigh_tresh, np.min(rayleigh)))
        print('Frequencies are far enough apart (always >%.6f, minimum is %.6f)'%(freq_diff_min,np.min(freq_diffs)))
    else:
        print('Rayleigh criterion vandalised (not always>%.2f, minimum is %.2f)'%(rayleigh_tresh, np.min(rayleigh)))
        print('Frequencies with not enough difference (not always >%.6f, minimum is %.6f)'%(freq_diff_min,np.min(freq_diffs)))
        for ray_id in rayleigh_bool_id:
            print(t_const_freq.iloc[[ray_id,ray_id+1]])


def check_ts(ts):
    """
    prints several statistics of the provided timeseries

    Parameters
    ----------
    ts : pandas.DataFrame
        The DataFrame should contain a 'values' column and a pd.DatetimeIndex as index, it contains the timeseries to be checked.

    Returns
    -------
    print_statement: str
        For printing as a substring of another string.

    """
    
    timesteps_min_all = ts.index.to_series().diff()[1:].dt.total_seconds()/60
    bool_int = (timesteps_min_all-timesteps_min_all.round(0))<1e-9
    if bool_int.all():
        timesteps_min_all = timesteps_min_all.astype(int)
    else: #in case of non integer minute timesteps (eg seconds)
        timesteps_min_all[bool_int] = timesteps_min_all[bool_int].round(0)
    timesteps_min = set(timesteps_min_all)
    #print(timesteps_min)
    if len(timesteps_min)<=100:
        timesteps_min_print = timesteps_min
    else:
        timesteps_min_print = 'too much unique time intervals (>100) to display all of them, %i intervals ranging from %i to %i minutes'%(len(timesteps_min),np.min(list(timesteps_min)),np.max(list(timesteps_min)))
    if (timesteps_min_all>0).all():
        timesteps_incr_print = 'all time intervals are in increasing order and are never equal'
    else:
        timesteps_incr_print = 'the times-order of ts is not always increasing (duplicate values or wrong order)'
    
    ntimes_nonan = ts['values'].count()
    ntimes = len(ts)
    ntimesteps_uniq = len(timesteps_min)
    if len(ts)==0:
        print_statement = 'timeseries contents:\n%s'%(ts)
        print(print_statement)
    else:
        list_statements = ['timeseries contents:\n%s'%(ts),
                           'timeseries # unique timesteps: %i'%(ntimesteps_uniq),
                           'timeseries unique timesteps (minutes):\n%s'%(timesteps_min_print),
                           'timeseries validity: %s'%(timesteps_incr_print),
                           'timeseries length: %i'%(ntimes),
                           'timeseries # nonan: %i'%(ntimes_nonan),
                           'timeseries %% nonan: %.1f%%'%(ntimes_nonan/ntimes*100),
                           'timeseries # nan: %i'%(ntimes-ntimes_nonan),
                           'timeseries %% nan: %.1f%%'%((ntimes-ntimes_nonan)/ntimes*100)]
        print_statement = '\n'.join(list_statements)
        print(print_statement)
    
    return print_statement
    
    
###############################
################# READING FILES
###############################


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
    
    with open(filename, encoding='latin1') as f: #'latin1 is nodig om predictie diafile die rechtstreeks uit hatyan komen in te lezen (validatietijdserie met op regel 4 (PAR) ongeldige tekens aan het einde)
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
    if diablocks_pd_startstopstation.isnull().any().any():
        raise Exception('ERROR: multiple blocks in diafile, but unequal amount of start/end/datastart/stationnames')
    
    #convert columns with line numbers to integers
    diablocks_pd_startstopstation[linenum_colnames] = diablocks_pd_startstopstation[linenum_colnames].astype(int)
    
    return diablocks_pd_startstopstation


def get_diablocks(filename):
    
    print('reading file: %s'%(filename))
    diablocks_pd = get_diablocks_startstopstation(filename)
    for block_id in diablocks_pd.index.tolist():
        #read diafile metadata as pandas series, prevent splitting of combined paramater names like MXH;2 by replacing ; with !
        data_meta_nrows = diablocks_pd.loc[block_id,'data_starts'] - diablocks_pd.loc[block_id,'block_starts']
        data_meta_series = pd.read_table(filename,skiprows=diablocks_pd.loc[block_id,'block_starts'],nrows=data_meta_nrows,header=None)[0] #series of metadata
        if not data_meta_series.str.contains('GHD|MXG;2').any(): #wia files contain these parameters, dia files don't. Replace dia names with wia names (wia files also contain PAR and MXP;2, but they should not be replaced)
            data_meta_series = data_meta_series.str.replace('PAR','GHD').str.replace('MXP;2','MXG;2')
        bool_combinedparname = (data_meta_series.str[3:6]==';1;') | (data_meta_series.str[3:6]==';2;')
        data_meta_series.loc[bool_combinedparname] = data_meta_series.loc[bool_combinedparname].str.slice_replace(3,4,'!')
        
        #get groepering and whether dia/wia is equidistant or non-equidistant
        bool_startswithmux = data_meta_series.str.startswith('MUX')
        if bool_startswithmux.any(): #extreme waterlevel timeseries (non-equidistant)
            mincontent = ['MXG;2','LOC','MXH;2','MXE;2','TYD']
            diablocks_pd.loc[block_id,'groepering'] = data_meta_series.loc[bool_startswithmux].iloc[0].split(';')[1]
        else: #normal waterlevel timeseries (equidistant)
            mincontent = ['GHD',  'LOC','HDH',  'EHD',  'TYD']
            diablocks_pd.loc[block_id,'groepering'] = 'NVT'
        
        #read all required metadata
        for get_content_sel in mincontent:
            bool_mincontent = data_meta_series.str.replace('!',';').str.startswith(get_content_sel)
            if bool_mincontent.sum()!=1:
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
                    print('no coordinate data available in LOC line of dia file')
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
                    timestep_value = None
                elif len(data_meta_mincontent)==7: #equidistant timeseries contains also timeunit and timestep
                    timestep_unit = data_meta_mincontent[6]
                    if timestep_unit != 'min':
                        raise Exception('ERROR: time unit from TYD is in unknown format (not "min")')
                    timestep_value = int(data_meta_mincontent[5]) #int(timestep_value_raw)
                else:
                    raise Exception(f'ERROR: time metadata is not understood: {data_meta_mincontent}')
                diablocks_pd.loc[block_id,'tstart'] = datestart
                diablocks_pd.loc[block_id,'tstop'] = datestop
                diablocks_pd.loc[block_id,'timestep_min'] = timestep_value

    return diablocks_pd


def readts_dia_nonequidistant(filename, diablocks_pd, block_id):

    data_nrows = diablocks_pd.loc[block_id,'data_ends'] - diablocks_pd.loc[block_id,'data_starts']
    data_pd_HWLW = pd.read_csv(filename,skiprows=diablocks_pd.loc[block_id,'data_starts'],nrows=data_nrows, header=None, names=['date','time','HWLWcode/qualitycode','valuecm:'], sep=';', parse_dates={'times':[0,1]})
    
    #convert HWLW+quality code to separate columns
    data_pd_HWLWtemp = data_pd_HWLW.loc[:,'HWLWcode/qualitycode'].str.split('/', expand=True)
    data_pd_HWLW['HWLWcode'] = data_pd_HWLWtemp.iloc[:,0].astype('int')
    data_pd_HWLW['qualitycode'] = data_pd_HWLWtemp.iloc[:,1].astype('int')
    data_pd_HWLW = data_pd_HWLW.drop('HWLWcode/qualitycode',axis='columns')

    #convert value from cm to m
    data_pd_HWLW['values'] = data_pd_HWLW['valuecm:'].str.strip(':').astype('int')/100
    data_pd_HWLW = data_pd_HWLW.drop('valuecm:',axis='columns')
    
    bool_hiaat = data_pd_HWLW['qualitycode'] == 99
    data_pd_HWLW.loc[bool_hiaat,'values'] = np.nan
    
    data_pd = data_pd_HWLW
    data_pd = data_pd.set_index('times')
    
    return data_pd


def readts_dia_equidistant(filename, diablocks_pd, block_id):
    
    datestart = diablocks_pd.loc[block_id,'tstart']
    datestop = diablocks_pd.loc[block_id,'tstop']
    timestep_min = diablocks_pd.loc[block_id,'timestep_min']
    times_fromfile = pd.date_range(start=datestart,end=datestop,freq='%dmin'%(timestep_min))
    
    #get data for station
    data_nrows = diablocks_pd.loc[block_id,'data_ends'] - diablocks_pd.loc[block_id,'data_starts']
    data_pd = pd.read_csv(filename,skiprows=diablocks_pd.loc[block_id,'data_starts'],nrows=data_nrows, header=None)
    data_pdser = data_pd[0].str.strip()
    data = data_pdser.str.cat()
    data = data.strip(':') #remove first and/or last colon if present
    data = data.split(':')
    
    if len(times_fromfile) != len(data):
        raise Exception('ERROR: times and values ts are not of equal length\nlen(times_fromfile): %d\nlen(data): %d'%(len(times_fromfile),len(data)))
    data_pd = pd.DataFrame({'times':times_fromfile,'valuecm/qualitycode':data})
    
    #convert HWLW+quality code to separate columns
    data_pd_temp = data_pd.loc[:,'valuecm/qualitycode'].str.split('/', expand=True)
    data_pd['values'] = data_pd_temp.iloc[:,0].astype('int')/100
    data_pd['qualitycode'] = data_pd_temp.iloc[:,1].astype('int')
    data_pd = data_pd.drop('valuecm/qualitycode',axis='columns')

    bool_hiaat = data_pd['qualitycode'] == 99
    data_pd.loc[bool_hiaat,'values'] = np.nan
    data_pd = data_pd.set_index('times')
    
    return data_pd


def readts_dia(filename, station=None, block_ids=None):
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

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    data_pd : pandas.core.frame.DataFrame
        DataFrame with a 'values' column and a pd.DatetimeIndex as index in case of an equidistant file, or more columns in case of a non-equidistant file.

    """
    
    if not isinstance(filename,list):
        filename = [filename]
        
    data_pd_all = pd.DataFrame()
    for iF, filename_one in enumerate(filename):    
        diablocks_pd = get_diablocks(filename_one)
        pd.set_option('display.max_columns', 6) #default was 0, but need more to display groepering
        pd.set_option('display.width', 200) #default was 80, but need more to display groepering
        print_cols = ['block_starts', 'station', 'grootheid', 'groepering', 'tstart', 'tstop']
        print('blocks in diafile:\n%s'%(diablocks_pd[print_cols]))
        str_getdiablockspd = 'A summary of the available blocks is printed above, obtain a full DataFrame of available diablocks with "diablocks_pd=Timeseries.get_diablocks(filename)"'
        
        #get equidistant timeseries from metadata
        if block_ids is None or block_ids=='allstation':
            if station is None:
                raise Exception('ERROR: if block_ids argument is not provided (or None) or is "allstation", station argument should be provided.')
            bool_station = diablocks_pd['station']==station
            ids_station = diablocks_pd[bool_station].index.tolist()
            if len(ids_station)<1:
                raise Exception('ERROR: no data block with requested station (%s) present in dia file. %s'%(station, str_getdiablockspd))
            elif len(ids_station)>1 and block_ids is None:
                raise Exception('ERROR: more than one data block with requested station (%s) present in dia file. Provide block_ids argument to readts_dia() (int, list of int or "allstation"). %s'%(station, str_getdiablockspd))
            else: #exactly one occurrence or block_ids is provided or block_ids='allstation'
                block_ids = ids_station
        
        #check validity of blockids of type listlist
        if isinstance(block_ids,int):
            block_ids = [block_ids]
        if not isinstance(block_ids,list):
            raise Exception('ERROR: invalid type for block_ids (should be int, list of int or "allstation")')
        if not pd.Series(block_ids).isin(diablocks_pd.index).all():
            raise Exception(f'ERROR: invalid values in block_ids list ({block_ids}), possible are {diablocks_pd.index.tolist()} (all integers)')
            
        if station is not None:
            if not isinstance(station,str):
                raise Exception('ERROR: station argument should be of type string')
            bool_samestation = diablocks_pd.loc[block_ids,'station']==station
            if not bool_samestation.all():
                raise Exception('ERROR: both the arguments station and block_ids are provided, but at least one of the requested block_ids corresponds to a different station. %s'%(str_getdiablockspd))
            
        data_pd_allblocks = pd.DataFrame()
        for block_id in block_ids:
            if np.isnan(diablocks_pd.loc[block_id,'timestep_min']): #non-equidistant
                data_pd_oneblock = readts_dia_nonequidistant(filename_one, diablocks_pd, block_id)
            else: #equidistant
                data_pd_oneblock = readts_dia_equidistant(filename_one, diablocks_pd, block_id)
            #check_ts(data_pd_oneblock)
            data_pd_allblocks = data_pd_allblocks.append(data_pd_oneblock, ignore_index=False)
        
        #append to allyears dataset
        data_pd_all = data_pd_all.append(data_pd_allblocks, ignore_index=False)

    #check overlapping timesteps, sort values on time and check_ts
    if len(data_pd_all) != len(data_pd_all.index.unique()):
        raise Exception('ERROR: merged datasets have duplicate/overlapping timesteps, clean up your input data or provide one file instead of a list')
    data_pd_all = data_pd_all.sort_index(axis=0)
    check_ts(data_pd_all)
    
    return data_pd_all


def readts_noos(filename, datetime_format='%Y%m%d%H%M', na_values=None):
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
    
    print('-'*50)
    print('reading file: %s'%(filename))
    noosheader = []
    noosheader_dict = {}
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
    
    content_pd = pd.read_csv(filename,header=startdata-1,delim_whitespace=True,names=['times_str','values'], na_values=na_values)
    noos_datetime = pd.to_datetime(content_pd['times_str'],format=datetime_format)
    data_pd = pd.DataFrame({'values':content_pd['values'].values},index=noos_datetime)
    
    check_ts(data_pd)
    return data_pd

