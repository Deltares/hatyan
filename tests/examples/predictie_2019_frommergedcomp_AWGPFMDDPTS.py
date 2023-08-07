# -*- coding: utf-8 -*-
"""
predictie_2019_frommergedcomp_AWGPFM_testDDPTS_test
shows the difference between calculating the nodal factors:
    - in the middle of the prediction period (fu_alltimes=False, hatyan default)
    - for all timesteps in the prediction period (fu_alltimes=True)
The results show that with the latter method, overlapping periods match perfectly

"""

import os
import datetime as dt
import hatyan

#dir_testdata = 'P:\\1209447-kpp-hydraulicaprogrammatuur\\hatyan\\hatyan_data_acceptancetests'
dir_testdata = 'C:\\DATA\\hatyan_data_acceptancetests'

selected_stations = ['AWGPFM']

for current_station in selected_stations:
    
    file_data_comp = os.path.join(dir_testdata,'predictie2019','%s_ana.txt'%(current_station))
    
    times_ext_pred13 = [dt.datetime(2019,1,1),dt.datetime(2020,1,1)]
    times_ext_pred24 = [dt.datetime(2019,11,1),dt.datetime(2020,3,1)]
    times_step_pred = 10
    
    file_data_predvali = os.path.join(dir_testdata,'predictie2019','%s_pre.txt'%(current_station))
    
    #components
    COMP_merged = hatyan.read_components(filename=file_data_comp)
    COMP_validation = hatyan.read_components(filename=file_data_comp)
    
    #prediction and validation
    ts_prediction1 = hatyan.prediction(comp=COMP_merged, nodalfactors=True, xfac=True, fu_alltimes=False, times_ext=times_ext_pred13, timestep_min=times_step_pred)
    ts_prediction2 = hatyan.prediction(comp=COMP_merged, nodalfactors=True, xfac=True, fu_alltimes=False, times_ext=times_ext_pred24, timestep_min=times_step_pred)
    ts_prediction3 = hatyan.prediction(comp=COMP_merged, nodalfactors=True, xfac=True, fu_alltimes=True, times_ext=times_ext_pred13, timestep_min=times_step_pred)
    ts_prediction4 = hatyan.prediction(comp=COMP_merged, nodalfactors=True, xfac=True, fu_alltimes=True, times_ext=times_ext_pred24, timestep_min=times_step_pred)
    ts_validation = hatyan.readts_dia(filename=file_data_predvali, station=current_station)

    fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction1, ts_validation=ts_validation)
    fig.savefig('prediction_%im_%s_validation_default'%(times_step_pred, current_station))
    fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction1, ts_validation=ts_prediction2)
    fig.savefig('prediction_%im_%s_validation_default_twoperiods'%(times_step_pred, current_station))
    fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction3, ts_validation=ts_validation)
    fig.savefig('prediction_%im_%s_validation_DDPTST'%(times_step_pred, current_station))
    fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction3, ts_validation=ts_prediction4)
    fig.savefig('prediction_%im_%s_validation_DDPTST_twoperiods'%(times_step_pred, current_station))
