# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 09:53:03 2021

@author: veenstra
"""

import os
import sys
import datetime as dt
import hatyan

file_config = os.path.realpath(__file__)
dir_output, timer_start = hatyan.init_RWS(file_config, sys.argv, interactive_plots=False)
#dir_testdata = 'P:\\1209447-kpp-hydraulicaprogrammatuur\\hatyan\\hatyan_data_acceptancetests'
dir_testdata = 'C:\\DATA\\hatyan_data_acceptancetests'

current_station = 'VLISSGN'

#file_slotgemiddelden = os.path.join(dir_testdata,'predictie2019','_slotgemiddelden_predictie2019.txt')
#stations_slotgem = pd.read_csv(file_slotgemiddelden, names=['slotgemiddelde'], comment='#', delim_whitespace=True)

for fu_alltimes in [True,False]:
    const_list = hatyan.get_const_list_hatyan('springneap')
    #file_data_comp0 = os.path.join(dir_testdata,'predictie2019','%s_pre.txt'%(current_station))    
    file_data_comp0 = os.path.join(dir_testdata,'predictie2019','%s_obs4.txt'%(current_station))
    
    file_data_compvali = os.path.join(dir_testdata,'predictie2019','%s_ana.txt'%(current_station))
    
    times_ext_pred = [dt.datetime(2019,1,1),dt.datetime(2020,1,1)]
    times_step_pred = 10

    file_data_predvali = os.path.join(dir_testdata,'predictie2019','%s_pre.txt'%(current_station))
    #file_data_predvaliHWLW = os.path.join(dir_testdata,'predictie2019','%s_ext.txt'%(current_station))
    
    ts_measurements_group0 = hatyan.readts_dia(filename=file_data_comp0, station=current_station)
    ts_measurements_group0['values'] -=.1
    #ts_measurements_group0 = hatyan.crop_timeseries(ts_measurements_group0, times_ext=[dt.datetime(2012,1,1),dt.datetime(2012,2,1)])
    times_ext_comp0 = [ts_measurements_group0.index[0],ts_measurements_group0.index[-1]]
    times_step_comp0 = (ts_measurements_group0.index[1]-ts_measurements_group0.index[0]).total_seconds()/60
    
    #COMP_merged = hatyan.get_components_from_ts(ts=ts_measurements_group0, const_list=const_list, fu_alltimes=True)
    COMP_merged_xfac0 = hatyan.analysis(ts=ts_measurements_group0, const_list=const_list, fu_alltimes=fu_alltimes, xfac=False)#, return_prediction=True)
    COMP_merged_xfac1 = hatyan.analysis(ts=ts_measurements_group0, const_list=const_list, fu_alltimes=fu_alltimes, xfac=True)#, return_prediction=True)
    print(COMP_merged_xfac0)
    print(COMP_merged_xfac1)
    #fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_predCOMPyear, ts_validation=ts_measurements_group0)

    #print(COMP_merged['A'].max())
    #COMP_validation = hatyan.read_components(filename=file_data_compvali)
    #fig, (ax1,ax2) = hatyan.plot_components(COMP_merged, comp_validation=COMP_validation)
    
    #prediction and validation
    #ts_prediction = hatyan.prediction(comp=COMP_merged, fu_alltimes=True, times_ext=times_ext_pred, timestep_min=times_step_pred)
    ts_prediction_xfac0_xfac0 = hatyan.prediction(comp=COMP_merged_xfac0, fu_alltimes=fu_alltimes, xfac=False, times_ext=times_ext_pred, timestep_min=times_step_pred)
    ts_prediction_xfac0_xfac1 = hatyan.prediction(comp=COMP_merged_xfac0, fu_alltimes=fu_alltimes, xfac=True, times_ext=times_ext_pred, timestep_min=times_step_pred)
    ts_prediction_xfac1_xfac0 = hatyan.prediction(comp=COMP_merged_xfac1, fu_alltimes=fu_alltimes, xfac=False, times_ext=times_ext_pred, timestep_min=times_step_pred)
    ts_prediction_xfac1_xfac1 = hatyan.prediction(comp=COMP_merged_xfac1, fu_alltimes=fu_alltimes, xfac=True, times_ext=times_ext_pred, timestep_min=times_step_pred)

    #ts_prediction1min = hatyan.prediction(comp=COMP_merged, nodalfactors=nodalfactors, xfac=xfac, fu_alltimes=False, times_ext=times_ext_pred, timestep_min=1)
    #ts_validation = hatyan.readts_dia(filename=file_data_predvali, station=current_station)
    fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction_xfac0_xfac0, ts_validation=ts_prediction_xfac1_xfac1)
    fig.savefig(os.path.join(dir_output,'fualltimes%d_anaxfac0_predxfac1.png'%(fu_alltimes)))
    fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction_xfac1_xfac1, ts_validation=ts_prediction_xfac0_xfac0)
    fig.savefig(os.path.join(dir_output,'fualltimes%d_anaxfac1_predxfac0.png'%(fu_alltimes)))
    fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction_xfac0_xfac0, ts_validation=ts_prediction_xfac0_xfac1)
    fig.savefig(os.path.join(dir_output,'fualltimes%d_anaxfac0_predxfacdiff.png'%(fu_alltimes)))
    fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction_xfac1_xfac0, ts_validation=ts_prediction_xfac1_xfac1)
    fig.savefig(os.path.join(dir_output,'fualltimes%d_anaxfac1_predxfacdiff.png'%(fu_alltimes)))
    fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction_xfac0_xfac0, ts_validation=ts_prediction_xfac1_xfac0)
    fig.savefig(os.path.join(dir_output,'fualltimes%d_anaxfacdiff_predxfac0.png'%(fu_alltimes)))
    fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction_xfac0_xfac1, ts_validation=ts_prediction_xfac1_xfac1)
    fig.savefig(os.path.join(dir_output,'fualltimes%d_anaxfacdiff_predxfac1.png'%(fu_alltimes)))
    #fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction, ts_validation=ts_measurements_group0)

    
hatyan.exit_RWS(timer_start)


