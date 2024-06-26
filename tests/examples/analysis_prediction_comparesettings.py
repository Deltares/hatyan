# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 09:53:03 2021

@author: veenstra

analyse three years and predict for forth year, compare analysis/prediction settings and schureman/foreman
"""

import os
import pandas as pd
import numpy as np
import hatyan

#dir_testdata = 'P:\\1209447-kpp-hydraulicaprogrammatuur\\hatyan\\hatyan_data_acceptancetests'
dir_testdata = 'C:\\DATA\\hatyan_data_acceptancetests'

selected_stations = ['HOEKVHLD']#,'VLISSGN','CUXHVN']

#file_slotgemiddelden = os.path.join(dir_testdata,'predictie2019','_slotgemiddelden_predictie2019.txt')
#stations_slotgem = pd.read_csv(file_slotgemiddelden, names=['slotgemiddelde'], comment='#', sep="\\s+")

stats_list = []

for current_station in selected_stations:
    const_list = hatyan.get_const_list_hatyan('year')
    file_data_comp0_lastyear = os.path.join(dir_testdata,'predictie2019','%s_obs4.txt'%(current_station))
    file_data_comp0_raw = [os.path.join(dir_testdata,'predictie2019','%s_obs%i.txt'%(current_station, file_id)) for file_id in [1,2,3]]
    file_data_comp0 = [x for x in file_data_comp0_raw if os.path.exists(x)] #slim filename list down to available files/years
    
    ts_measurements_group0_lastyear = hatyan.read_dia(filename=file_data_comp0_lastyear, station=current_station)
    ts_measurements_group0 = hatyan.read_dia(filename=file_data_comp0, station=current_station)
    #ts_measurements_group0 = hatyan.crop_timeseries(ts_measurements_group0, times_ext=[dt.datetime(2012,1,1),dt.datetime(2013,1,1)])

    stats_row = pd.DataFrame(index=[current_station])
    
    #prediction and comparison to measurements
    for fu_alltimes in [True,False]:
        for xfac in [True, False]:
            for source in ['schureman','foreman']:
                COMP_merged = hatyan.analysis(ts=ts_measurements_group0, const_list=const_list, analysis_perperiod='Y', fu_alltimes=fu_alltimes, xfac=xfac, source=source)
                #print(COMP_merged)
                ts_prediction = hatyan.prediction(comp=COMP_merged, times=ts_measurements_group0_lastyear.index)
                overlapdiff = ts_prediction['values']-ts_measurements_group0_lastyear['values']
                rmse = np.sqrt(np.nanmean(overlapdiff ** 2))
                #fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction, ts_validation=ts_measurements_group0)
                #ax2.set_ylim(-0.3,0.3)
                stats_row['fu_alltimes=%s_xfac=%s_%s'%(fu_alltimes,xfac,source)] = [rmse]
    #print(stats_row)
    stats_list.append(stats_row)

stats = pd.concat(stats_list,axis=0)
statsT = stats.T
print('RMSE values [cm] for several settings and stations:')
print((statsT*100).round(3))
