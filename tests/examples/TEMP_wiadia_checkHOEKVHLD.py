# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 18:33:22 2021

@author: veenstra
"""

import os, sys
import datetime as dt
import hatyan
hatyan.close('all')

file_config = os.path.realpath(__file__)
#dir_output, timer_start = hatyan.init_RWS(file_config, sys.argv, interactive_plots=False)
dir_output = 'TEMP_output'

dir_testdata = 'C:\\DATA\\hatyan_data_acceptancetests'

selected_stations = ['HOEKVHLD']

for current_station in selected_stations:
    nodalfactors = True
    xfac=True
    analysis_peryear=True
    const_list = hatyan.get_const_list_hatyan('year') #94 const
    vertref='NAP'
    times_ext_pred = [dt.datetime(2019,1,1),dt.datetime(2020,1,1)]
    times_step_pred = 10
    
    """
    file_data_comp0_raw = [os.path.join(dir_testdata,'','%s_obs%i.txt'%(current_station, file_id)) for file_id in [1,2,3,4]]
    file_data_comp0 = [x for x in file_data_comp0_raw if os.path.exists(x)] #slim filename list down to available files/years
    file_data_comp1 = os.path.join(dir_testdata,'','%s_ana.txt'%(current_station))
    file_data_compvali = os.path.join(dir_testdata,'','%s_ana.txt'%(current_station))
    file_data_predvali = os.path.join(dir_testdata,'','%s_pre.txt'%(current_station))
    file_data_predvaliHWLW = os.path.join(dir_testdata,'','%s_ext.txt'%(current_station))
    
    ts_measurements_group0 = hatyan.readts_dia(filename=file_data_comp0, station=current_station)
    times_ext_comp0 = [ts_measurements_group0.index[0],ts_measurements_group0.index[-1]]
    times_step_comp0 = (ts_measurements_group0.index[1]-ts_measurements_group0.index[0]).total_seconds()/60

    comp_frommeasurements_avg_group0, comp_frommeasurements_all_group0 = hatyan.get_components_from_ts(ts=ts_measurements_group0, const_list=const_list, nodalfactors=nodalfactors, xfac=xfac, fu_alltimes=False, analysis_peryear=analysis_peryear, return_allyears=True)

    #prediction and validation
    ts_prediction = hatyan.prediction(comp=comp_frommeasurements_avg_group0, nodalfactors=nodalfactors, xfac=xfac, fu_alltimes=False, times_ext=times_ext_pred, timestep_min=times_step_pred)
    ts_ext_prediction = hatyan.calc_HWLW(ts=ts_prediction)
    ts_validation = hatyan.readts_dia(filename=file_data_predvali, station=current_station)
    ts_ext_validation = hatyan.readts_dia(filename=file_data_predvaliHWLW, station=current_station)
    hatyan.write_tsdia(ts=ts_prediction, station=current_station, vertref=vertref, filename='prediction_%im_%s.dia'%(times_step_pred,current_station))
    hatyan.write_tsdia_HWLW(ts_ext=ts_ext_prediction, station=current_station, vertref=vertref, filename='prediction_HWLW_%im_%s.dia'%(times_step_pred, current_station))
    fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction, ts_ext=ts_ext_prediction, ts_validation=ts_validation, ts_ext_validation=ts_ext_validation)
    """
    file_dia_wl = os.path.join(dir_testdata,'other','diawia_HOEKVHLD_astro_tijdreeks.dia')
    file_dia_ext = os.path.join(dir_testdata,'other','diawia_HOEKVHLD_astro_extremen.dia')
    file_wia_wl = os.path.join(dir_testdata,'other','diawia_HOEKVHLD_astro_tijdreeks.wia')
    file_wia_ext = os.path.join(dir_testdata,'other','diawia_HOEKVHLD_astro_extremen.wia')
    
    ts_dia_wl = hatyan.readts_dia(filename=file_dia_wl, station=current_station)
    ts_dia_ext = hatyan.readts_dia(filename=file_dia_ext, station=current_station)
    ts_wia_wl = hatyan.readts_dia(filename=file_wia_wl, station=current_station)
    ts_wia_ext = hatyan.readts_dia(filename=file_wia_ext, station=current_station)
    
    hatyan.write_tsdia(ts=ts_dia_wl, station=current_station, vertref=vertref, filename='diawia_HOEKVHLD_astro_tijdreeks_out.dia')
    hatyan.write_tsdia_HWLW(ts_ext=ts_dia_ext, station=current_station, vertref=vertref, filename='diawia_HOEKVHLD_astro_extremen_out.dia')
    #hatyan.write_tsdia(ts=ts_wia_wl, station=current_station, vertref=vertref, filename='diawia_HOEKVHLD_astro_tijdreeks_out.wia')
    #hatyan.write_tsdia_HWLW(ts_ext=ts_wia_ext, station=current_station, vertref=vertref, filename='diawia_HOEKVHLD_astro_extremen_out.wia')
    
    assert (ts_dia_wl==ts_wia_wl).all().all()
    assert (ts_dia_ext==ts_wia_ext).all().all()
#hatyan.exit_RWS(timer_start) #provides footer to outputfile when calling this script with python






