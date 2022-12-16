# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 14:12:54 2022

@author: veenstra
"""
import matplotlib.pyplot as plt
import numpy as np
import datetime as dt


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
