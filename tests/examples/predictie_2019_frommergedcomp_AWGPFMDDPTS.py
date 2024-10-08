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

dir_testdata = 'C:\\DATA\\hatyan_data_acceptancetests'

selected_stations = ['AWGPFM']

for current_station in selected_stations:
    print(f'current_station: {current_station}')
    
    file_data_comp = os.path.join(dir_testdata,'predictie2019','%s_ana.txt'%(current_station))
    
    times_pred13 = slice(dt.datetime(2019,1,1),dt.datetime(2020,1,1), "10min")
    times_pred24 = slice(dt.datetime(2019,11,1),dt.datetime(2020,3,1), "10min")
    
    file_data_predvali = os.path.join(dir_testdata,'predictie2019','%s_pre.txt'%(current_station))
    
    #components
    comp_merged = hatyan.read_components(filename=file_data_comp)
    
    #prediction and validation
    ts_prediction1 = hatyan.prediction(comp=comp_merged, times=times_pred13)
    ts_prediction2 = hatyan.prediction(comp=comp_merged, times=times_pred24)
    comp_merged.attrs["fu_alltimes"] = True
    ts_prediction3 = hatyan.prediction(comp=comp_merged, times=times_pred13)
    ts_prediction4 = hatyan.prediction(comp=comp_merged, times=times_pred24)
    ts_validation = hatyan.read_dia(filename=file_data_predvali, station=current_station)

    fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction1, ts_validation=ts_validation)
    fig.savefig('prediction_%s_%s_validation_default'%(times_pred13.step, current_station))
    fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction1, ts_validation=ts_prediction2)
    fig.savefig('prediction_%s_%s_validation_default_twoperiods'%(times_pred13.step, current_station))
    fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction3, ts_validation=ts_validation)
    fig.savefig('prediction_%s_%s_validation_DDPTST'%(times_pred13.step, current_station))
    fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction3, ts_validation=ts_prediction4)
    fig.savefig('prediction_%s_%s_validation_DDPTST_twoperiods'%(times_pred13.step, current_station))
