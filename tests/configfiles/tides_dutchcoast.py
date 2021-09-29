# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 15:11:07 2019

@author: veenstra
"""
import os, sys
#sys.path.append(r'c:\DATA\hatyan_python')
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
plt.close('all')
import matplotlib.dates as mdates

from hatyan.analysis_prediction import prediction
import hatyan.components as Components
from hatyan.hatyan_core import get_const_list_hatyan
from hatyan.wrapper_RWS import init_RWS, exit_RWS

file_config = os.path.realpath(__file__) #F9 doesnt work, only F5 (F5 also only method to reload external definition scripts)
dir_output, timer_start = init_RWS(file_config, sys.argv, interactive_plots=False) #provides header to outputfile when calling this script with python
#dir_testdata = 'P:\\1209447-kpp-hydraulicaprogrammatuur\\hatyan\\hatyan_data_acceptancetests'
dir_testdata = 'C:\\DATA\\hatyan_data_acceptancetests'

polar_fig = False

const_list = get_const_list_hatyan('year')
selected_stations = ['CADZD','BATH','VLISSGN','HOEKVHLD','IJMDBTHVN','DENHDR','TERSLNZE','SCHIERMNOG','DELFZL']
selected_stations_names = ['Cadzand','Bath','Vlissingen','Hoek van Holland','IJmuiden Buitenhaven','Den Helder','Terschelling','Schiermonnikoog','Delfzijl']

times_ext = [dt.datetime(2009,1,1),dt.datetime(2012,12,31,23,0)]

tstart = dt.datetime(2019,1,6,2,30) #nieuwe maan op 6 jan
tstart = dt.datetime(2019,1,6,13,1) #midden tussen opkomst en ondergang van de maan, dus maansdoorgang?
times_ext_pred = [tstart, tstart+dt.timedelta(hours=24, minutes=49)]
timestep_pred = 1


if polar_fig:
    fig_polar = plt.figure(figsize=(9,9))
    ax1_polar = fig_polar.add_subplot(111, polar=True)
fig, ax1 = plt.subplots(1,1,figsize=(14,5))
n_colors = len(selected_stations)
colors = plt.cm.jet(np.linspace(0,1,n_colors))
for i_stat, current_station in enumerate(selected_stations):
    
    comp_frommeasurements_avg_group = Components.read_components(filename=os.path.join(dir_testdata,'predictie2019','%s_ana.txt'%(current_station)))
    
    ts_prediction = prediction(comp=comp_frommeasurements_avg_group, times_ext=times_ext_pred, timestep_min=timestep_pred)
    
    vals_real = ts_prediction['values']
    times_real = ts_prediction.index
    if polar_fig:
        t_from0_sec = ((ts_prediction.index-ts_prediction.index[0])/np.timedelta64(1, 's')).values
        t_from0_2pi = t_from0_sec/np.max(t_from0_sec)*2*np.pi
        vals_norm = (vals_real-np.min(vals_real))/(-np.min(vals_real)+np.max(vals_real))*3.5
        ax1_polar.plot(t_from0_2pi, vals_norm, label=current_station, color=colors[i_stat])
    ax1.plot(times_real, vals_real, label=current_station, color=colors[i_stat])


if polar_fig:
    ax1_polar.plot([0],[1.75],'ok',markersize=15)
    ax1_polar.plot([np.pi],[1.75],'ok',markersize=15, mfc='none')
    ax1_polar.set_rorigin(-51.1/2) #binnendiameter 51.1cm, buitendiameter 58.2cm
    ax1_polar.legend()
    ax1_polar.set_rlim([-0.1,3.6])
    #ax1_polar.set_rlim([0,1])
    ax1_polar.set_rticks([0,1.75,3.5])
    ax1_polar.set_yticklabels([])
    ax1_polar.set_theta_direction(1)
    ax1_polar.set_theta_zero_location('N')
    ax1_polar.set_thetagrids([0,45,90,135,180,225,270,315], labels=['NM','','EK','','VM','','LK',''])
    #ax1_polar.set_thetagrids([0,45,90,135,180,225,270,315], labels=['','','','','','','',''])
    ax1_polar.set_frame_on(False)
    fig_polar.tight_layout()
    #fig_polar.savefig('tide_clock_low.png', dpi=100)
    fig_polar.savefig('tide_clock.png', dpi=250)
    #fig_polar.savefig('tide_clock.eps', format='eps')


ax1.grid()
#ax1.plot([0],[1.75],'ok',markersize=15)
#ax1.plot([2*np.pi],[1.75],'ok',markersize=15)
#ax1.plot([np.pi],[1.75],'ok',markersize=15, mfc='none')
ax1.legend(bbox_to_anchor=(1,1))
ax1.set_xlim(tstart,times_ext_pred[1])
ax1.xaxis.set_major_formatter(mdates.DateFormatter('%d %b %H:%M'))
#ax1.set_ylim([-0.1,3.6])
#ax1.set_yticks([0,1.75,3.5])
#ax1.set_yticklabels(['LOW','','HIGH'])
#ax1.set_xlim([-0.1,2*np.pi+0.1])
#ax1.set_xticks([np.deg2rad(x) for x in [0,45,90,135,180,225,270,315,360]])
#ax1.set_xticklabels(['NM','','EK','','VM','','LK','','NM'])
#ax1.set_frame_on(False)
fig.tight_layout()
fig.savefig('tide_clock_nonpolar_copy.png', dpi=250)

exit_RWS(timer_start)


