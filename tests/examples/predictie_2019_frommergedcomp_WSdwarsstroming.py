# -*- coding: utf-8 -*-
"""
Created on Thu Dec 24 10:33:54 2020

@author: veenstra
"""
import os
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt
plt.close('all')
import hatyan

dir_testdata = 'C:\\DATA\\hatyan_data_acceptancetests'

file_slotgemiddelden = os.path.join(dir_testdata,'predictie2019','_slotgemiddelden_predictie2019.txt')
stations_slotgem = pd.read_csv(file_slotgemiddelden, names=['slotgemiddelde'], comment='#', sep="\\s+")

fig, (ax1, ax2) = plt.subplots(2,1,figsize=(15,9))#, sharex=True, gridspec_kw={'height_ratios':[2,1]})
list_matig = []
list_sterk = []
for yr, ax in zip([2020,2021],[ax1,ax2]):

    times_pred = slice(dt.datetime(yr,1,1),dt.datetime(yr+1,1,1),1) #longer period with alternating aggers and no aggers, also eerste HW wordt als lokaal ipv primair HW gezien, also extra agger outside of 1stLW/agger/2ndLW sequence
    
    if yr == 2020:
        file_data_comp_HANSWT = os.path.join(dir_testdata,'predictie2019','HANSWT_ana.txt') #2009-2012
        file_data_comp_TERNZN = os.path.join(dir_testdata,'predictie2019','TERNZN_ana.txt') #2009-2012
    elif yr == 2021:
        file_data_comp_HANSWT = os.path.join(dir_testdata,'other','HANSWT_ana_new.txt') #newer period
        file_data_comp_TERNZN = os.path.join(dir_testdata,'other','TERNZN_ana_new.txt') #newer period
    COMP_merged_HANSWT = hatyan.read_components(filename=file_data_comp_HANSWT)
    COMP_merged_TERNZN = hatyan.read_components(filename=file_data_comp_TERNZN)
    
    ts_prediction_HANSWT = hatyan.prediction(comp=COMP_merged_HANSWT, times=times_pred)
    ts_prediction_TERNZN = hatyan.prediction(comp=COMP_merged_TERNZN, times=times_pred)
    
    ts_prediction_diff = ts_prediction_TERNZN-ts_prediction_HANSWT # TODO: metadata is dropped
    ts_prediction_diff['values'] = ts_prediction_diff['values'].round(2) #round to cm
    ts_prediction_diff_HWLW = hatyan.calc_HWLW(ts=ts_prediction_diff)
    
    value_matig = .75
    value_sterk = .80
    ts_prediction_diff_matig = ts_prediction_diff_HWLW[ts_prediction_diff_HWLW['values'] >= value_matig]
    ts_prediction_diff_sterk = ts_prediction_diff_HWLW[ts_prediction_diff_HWLW['values'] >= value_sterk]
    list_matig.append(ts_prediction_diff_matig)
    list_sterk.append(ts_prediction_diff_sterk)
    
    #fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction_diff)
    ax.plot(ts_prediction_diff.index, ts_prediction_diff['values'],'-',linewidth=0.7,markersize=1, label='verval')
    ax.plot(ts_prediction_diff_matig.index, ts_prediction_diff_matig['values'],'oy',markersize=5, label='matig (>%.2fm)'%(value_matig))
    ax.plot(ts_prediction_diff_sterk.index, ts_prediction_diff_sterk['values'],'or',markersize=5, label='sterk (>%.2fm)'%(value_sterk))
    ax.plot([ts_prediction_diff.index[0],ts_prediction_diff.index[-1]],[value_matig,value_matig],'-y')
    ax.plot([ts_prediction_diff.index[0],ts_prediction_diff.index[-1]],[value_sterk,value_sterk],'-r')
    #ax1.set_ylim(-1.2,1.7)
    ax.legend(loc=1)
    ax.grid()
    ax.set_xlabel('Tijd')
    ax.set_ylabel('astro verval %d [m] (TERNZN-HANSWT)'%yr)
    #hatyan.write_dia(ts=ts_ext_prediction_main, station=current_station, vertref='NAP', filename='prediction_HWLW_%im_%s_main.dia'%(times_step_pred, current_station))
    #hatyan.write_dia(ts=ts_ext_prediction_clean, station=current_station, vertref='NAP', filename='prediction_HWLW_%im_%s_agger345.dia'%(times_step_pred, current_station))
fig.tight_layout()
fig.savefig('prediction_WSdwarsstroming')

for yr, pd_matig, pd_sterk in zip([2020,2021],list_matig,list_sterk):
    print('matige vervaloverschrijdingen in', yr)
    print(pd_matig)
    print('sterke vervaloverschrijdingen in', yr)
    print(pd_sterk)

for yr, pd_matig, pd_sterk in zip([2020,2021],list_matig,list_sterk):
    print('%i matige vervaloverschrijdingen in %i'%(len(pd_matig),yr))
    print('%i sterke vervaloverschrijdingen in %i'%(len(pd_sterk),yr))
