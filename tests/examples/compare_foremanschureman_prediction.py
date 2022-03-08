# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 09:53:03 2021

@author: veenstra

analyse three years and predict for forth year, compare analysis/prediction settings and schureman/foreman
"""

import os, sys
import pandas as pd
import hatyan

file_config = os.path.realpath(__file__)
dir_output, timer_start = hatyan.init_RWS(file_config, sys.argv, interactive_plots=False)
#dir_testdata = 'P:\\1209447-kpp-hydraulicaprogrammatuur\\hatyan\\hatyan_data_acceptancetests'
dir_testdata = 'C:\\DATA\\hatyan_data_acceptancetests'

selected_stations = ['HOEKVHLD']#,'VLISSGN','CUXHVN']

for current_station in selected_stations:
    const_list = hatyan.get_const_list_hatyan('year')
    file_data_comp0_lastyear = os.path.join(dir_testdata,'predictie2019','%s_obs1.txt'%(current_station))
    #file_data_comp0_raw = [os.path.join(dir_testdata,'predictie2019','%s_obs%i.txt'%(current_station, file_id)) for file_id in [1,2,3]]
    #file_data_comp0 = [x for x in file_data_comp0_raw if os.path.exists(x)] #slim filename list down to available files/years
    
    ts_measurements_group0_lastyear = hatyan.readts_dia(filename=file_data_comp0_lastyear, station=current_station)
    #ts_measurements_group0 = hatyan.readts_dia(filename=file_data_comp0, station=current_station)
    #ts_measurements_group0 = hatyan.crop_timeseries(ts_measurements_group0, times_ext=[dt.datetime(2012,1,1),dt.datetime(2013,1,1)])

    stats_row = pd.DataFrame(index=[current_station])
    for fu_alltimes in [True]:#,False]:
        xfac = False
        return_prediction = True
        #prediction and comparison to measurements
        settings_schu = hatyan.HatyanSettings(source='schureman',fu_alltimes=fu_alltimes,xfac=xfac,return_prediction=return_prediction)
        settings_for = hatyan.HatyanSettings(source='foreman',fu_alltimes=fu_alltimes,xfac=xfac,return_prediction=return_prediction)
        COMP_schu, ts_prediction_schu = hatyan.analysis(ts=ts_measurements_group0_lastyear, const_list=const_list, hatyan_settings=settings_schu)
        COMP_for, ts_prediction_for = hatyan.analysis(ts=ts_measurements_group0_lastyear, const_list=const_list, hatyan_settings=settings_for)
        
        fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction_schu, ts_validation=ts_prediction_for)
        ax1.set_title(f'{current_station} fualltimes={fu_alltimes}')
        ax1.set_ylim(-2.5,2.5)
        ax2.set_ylim(-0.02,0.02)
        ax1.legend(['schureman','foreman','difference'],loc=4)
        fig.savefig(f'compare_foreman_schureman_{current_station}_fualltimes={fu_alltimes}')

hatyan.exit_RWS(timer_start)


