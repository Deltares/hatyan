# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 14:12:54 2022

@author: veenstra
"""
import matplotlib.pyplot as plt
import numpy as np
import datetime as dt
import pandas as pd

def reshape_signal(ts, ts_ext, HW_goal, LW_goal, tP_goal=None):
    """
    scales tidal signal to provided HW/LW value and up/down going time
    tP_goal (tidal period time) is used to fix tidalperiod to 12h25m (for BOI timeseries)
    
    time_down was scaled with havengetallen before, but not anymore to avoid issues with aggers
    """
    TR_goal = HW_goal-LW_goal
    
    #selecteer alle hoogwaters en opvolgende laagwaters
    bool_HW = ts_ext['HWLWcode']==1
    idx_HW = np.where(bool_HW)[0]
    timesHW = ts_ext.index[idx_HW]
    timesLW = ts_ext.index[idx_HW[:-1]+1] #assuming alternating 1,2,1 or 1,3,1, this is always valid in this workflow
    
    #crop from first to last HW (rest is not scaled anyway)
    ts_time_firstHW = ts_ext[bool_HW].index[0]
    ts_time_lastHW = ts_ext[bool_HW].index[-1]
    ts_corr = ts.copy().loc[ts_time_firstHW:ts_time_lastHW]

    ts_corr['times'] = ts_corr.index #this is necessary since datetimeindex with freq is not editable, and Series is editable
    ts_corr['values_new'] = np.nan #necessary since HW is read and overwitten twice
    for i in np.arange(0,len(timesHW)-1):
        HW1_val = ts_corr.loc[timesHW[i],'values']
        HW2_val = ts_corr.loc[timesHW[i+1],'values']
        LW_val = ts_corr.loc[timesLW[i],'values']
        TR1_val = HW1_val-LW_val
        TR2_val = HW2_val-LW_val
        tP_val = timesHW[i+1]-timesHW[i]
        if tP_goal is None:
            tP_goal = tP_val
        
        tide_HWtoHW = ts_corr.loc[timesHW[i]:timesHW[i+1]]
        
        ts_corr.loc[timesHW[i]:timesLW[i],'values_new'] = (ts_corr.loc[timesHW[i]:timesLW[i],'values']-LW_val)/TR1_val*TR_goal+LW_goal
        ts_corr.loc[timesLW[i]:timesHW[i+1],'values_new'] = (ts_corr.loc[timesLW[i]:timesHW[i+1],'values']-LW_val)/TR2_val*TR_goal+LW_goal
        ts_corr.loc[timesHW[i]:timesHW[i+1],'times'] = pd.date_range(start=ts_corr.loc[timesHW[i],'times'],end=ts_corr.loc[timesHW[i],'times']+tP_goal,periods=len(tide_HWtoHW))

    ts_corr = ts_corr.set_index('times',drop=True)
    ts_corr['values'] = ts_corr['values_new']
    ts_corr = ts_corr.drop(['values_new'],axis=1)
    return ts_corr


def ts_to_trefHW(ts,HWreftime):
    """
    converts to hours relative to HWreftime, to plot av/sp/np tidal signals in one plot
    """
    ts.index.name = 'times' #just to be sure
    ts_trefHW = ts.reset_index()
    ts_trefHW.index = (ts_trefHW['times']-HWreftime).dt.total_seconds()/3600
    return ts_trefHW


def repeat_signal(ts_one_HWtoHW, nb, na):
    """
    repeat tidal signal, necessary for sp/np, since they are computed as single tidal signal first
    """
    tidalperiod = ts_one_HWtoHW.index[-1] - ts_one_HWtoHW.index[0]
    ts_rep = pd.DataFrame()
    for iAdd in np.arange(-nb,na+1):
        ts_add = pd.DataFrame({'values':ts_one_HWtoHW['values'].values},
                              index=ts_one_HWtoHW.index + iAdd*tidalperiod)
        ts_rep = pd.concat([ts_rep,ts_add])
    ts_rep = ts_rep.loc[~ts_rep.index.duplicated()]
    return ts_rep


def plot_aardappelgrafiek(HWLW_culmhr_summary):
    def timeTicks(x, pos):
        d = dt.timedelta(hours=np.abs(x))
        if np.sign(x)>0:
            d_str = str(d)
        else:
            d_str = '-'+str(d)
        return d_str
    
    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(7.5,4), sharex=False)
    ax1.set_title('HW')
    ax1.set_xlabel('maansverloop in uu:mm:ss' )
    ax1.set_ylabel('waterstand in m t.o.v. NAP')
    ax1.plot(HWLW_culmhr_summary['HW_delay_median'].dt.total_seconds()/3600,HWLW_culmhr_summary['HW_values_median'],'.-')
    ax1.xaxis.set_major_formatter(timeTicks)
    ax1.grid()
    ax2.set_title('LW')
    ax2.set_xlabel('maansverloop in uu:mm:ss' )
    ax2.set_ylabel('waterstand in m t.o.v. NAP')
    ax2.plot(HWLW_culmhr_summary['LW_delay_median'].dt.total_seconds()/3600,HWLW_culmhr_summary['LW_values_median'],'.-')
    ax2.xaxis.set_major_formatter(timeTicks)
    ax2.grid()
    for iH,row in HWLW_culmhr_summary.iterrows():
        ax1.text(row['HW_delay_median'].total_seconds()/3600,row['HW_values_median'], str(int(iH)))
        ax2.text(row['LW_delay_median'].total_seconds()/3600,row['LW_values_median'], str(int(iH)))
    #set equal ylims
    ax1_xlimmean = np.mean(ax1.get_xlim())
    ax2_xlimmean = np.mean(ax2.get_xlim())
    ax1_ylimmean = np.mean(ax1.get_ylim())
    ax2_ylimmean = np.mean(ax2.get_ylim())
    xlimrange = 2
    ylimrange = 1
    ax1.set_xlim([ax1_xlimmean-xlimrange/2,ax1_xlimmean+xlimrange/2])
    ax2.set_xlim([ax2_xlimmean-xlimrange/2,ax2_xlimmean+xlimrange/2])
    ax1.set_ylim([ax1_ylimmean-ylimrange/2,ax1_ylimmean+ylimrange/2])
    ax2.set_ylim([ax2_ylimmean-ylimrange/2,ax2_ylimmean+ylimrange/2])
    #plot gemtij dotted lines
    ax1.plot(ax1.get_xlim(),[HWLW_culmhr_summary['HW_values_median'].mean(),HWLW_culmhr_summary['HW_values_median'].mean()],'k--')
    ax1.plot([HWLW_culmhr_summary['HW_delay_median'].mean().total_seconds()/3600,HWLW_culmhr_summary['HW_delay_median'].mean().total_seconds()/3600],ax1.get_ylim(),'k--')
    ax2.plot(ax2.get_xlim(),[HWLW_culmhr_summary['LW_values_median'].mean(),HWLW_culmhr_summary['LW_values_median'].mean()],'k--')
    ax2.plot([HWLW_culmhr_summary['LW_delay_median'].mean().total_seconds()/3600,HWLW_culmhr_summary['LW_delay_median'].mean().total_seconds()/3600],ax2.get_ylim(),'k--')
    fig.tight_layout()
    
    axs = (ax1,ax2)
    return fig, axs
