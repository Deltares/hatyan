# -*- coding: utf-8 -*-
"""
predictie_2019_19Ycomp4Ydia_VLISSGN_interactive.py
hatyan master configfile
voor alle stations indien mogelijk:
    - analyse van 4 jaar aan data
    - combineren met SA+SM uit analyseresultatenbestand
    - predictie maken

"""

import os
import datetime as dt
import pandas as pd
import hatyan
hatyan.close('all')

dir_testdata = 'C:\\DATA\\hatyan_github\\tests'

selected_stations = ['VLISSGN']

file_slotgemiddelden = os.path.join(dir_testdata,'data_unitsystemtests','_slotgemiddelden_predictie2019.txt')
stations_slotgem = pd.read_csv(file_slotgemiddelden, names=['slotgemiddelde'], comment='#', delim_whitespace=True)

for current_station in selected_stations:
    
    #START OF STATION SETTINGS
    #nodalfactors
    nodalfactors = True
    #xfactor
    xfac=True
    #analysis_perperiod
    analysis_perperiod='Y'
    #constituent list
    const_list = hatyan.get_const_list_hatyan('year') #94 const
    #component splitting
    CS_comps = None
    #vertical reference
    vertref='NAP'
    #END OF STATION SETTINGS
    

    file_data_comp0_raw = [os.path.join(dir_testdata,'data_unitsystemtests','%s_obs%i.txt'%(current_station, file_id)) for file_id in [1,2,3,4]]
    file_data_comp0 = [x for x in file_data_comp0_raw if os.path.exists(x)] #slim filename list down to available files/years
        
    file_data_comp1 = os.path.join(dir_testdata,'data_unitsystemtests','%s_ana.txt'%(current_station))
    
    file_data_compvali = os.path.join(dir_testdata,'data_unitsystemtests','%s_ana.txt'%(current_station))
    
    times_ext_pred = [dt.datetime(2019,1,1),dt.datetime(2020,1,1)]
    times_step_pred = 10

    file_data_predvali = os.path.join(dir_testdata,'data_unitsystemtests','%s_pre.txt'%(current_station))
    file_data_predvaliHWLW = os.path.join(dir_testdata,'data_unitsystemtests','%s_ext.txt'%(current_station))
    
    ts_measurements_group0 = hatyan.readts_dia(filename=file_data_comp0, station=current_station)
    times_ext_comp0 = [ts_measurements_group0.index[0],ts_measurements_group0.index[-1]]
    times_step_comp0 = (ts_measurements_group0.index[1]-ts_measurements_group0.index[0]).total_seconds()/60

    comp_frommeasurements_avg_group0, comp_frommeasurements_all_group0 = hatyan.get_components_from_ts(ts=ts_measurements_group0, const_list=const_list, nodalfactors=nodalfactors, xfac=xfac, fu_alltimes=False, analysis_perperiod=analysis_perperiod, return_allperiods=True)

    fig,(ax1,ax2) = hatyan.plot_components(comp_frommeasurements_avg_group0, comp_allyears=comp_frommeasurements_all_group0)
    fig.savefig('components_%s_4Y.png'%(current_station))
    comp_metadata = {'station':current_station, 'vertref':vertref, 'times_ext':[x.strftime('%Y%m%d%H%M') for x in times_ext_comp0], 'times_step':times_step_comp0, 'xfac':xfac}
    hatyan.write_components(comp_frommeasurements_avg_group0, filename='components_%s_4Y.txt'%(current_station), metadata=comp_metadata)
    
    comp_fromfile_group1 = hatyan.read_components(filename=file_data_comp1)
    
    #merge component groups (SA/SM from 19Y, rest from 4Y)
    COMP_merged = hatyan.merge_componentgroups(comp_main=comp_frommeasurements_avg_group0, comp_sec=comp_fromfile_group1, comp_sec_list=['SA','SM'])
    #replace A0 amplitude (middenstand) by slotgemiddelde
    if current_station in stations_slotgem.index.tolist():
        COMP_merged.loc['A0','A'] = stations_slotgem.loc[current_station,'slotgemiddelde']

    COMP_validation = hatyan.read_components(filename=file_data_compvali)
    fig, (ax1,ax2) = hatyan.plot_components(COMP_merged, comp_validation=COMP_validation)
    fig.savefig('components_%s_merged.png'%(current_station))
    comp_metadata = {'station':current_station, 'vertref':vertref, 'times_ext':times_ext_comp0, 'times_ext2':'SA and SM imported from analyseresultatenbestand', 'times_step':times_step_comp0, 'xfac':xfac}
    hatyan.write_components(COMP_merged, filename='components_%s_merged.txt'%(current_station), metadata=comp_metadata)
    
    #prediction and validation
    ts_prediction = hatyan.prediction(comp=COMP_merged, nodalfactors=nodalfactors, xfac=xfac, fu_alltimes=False, times_ext=times_ext_pred, timestep_min=times_step_pred)
    ts_prediction1min = hatyan.prediction(comp=COMP_merged, nodalfactors=nodalfactors, xfac=xfac, fu_alltimes=False, times_ext=times_ext_pred, timestep_min=1)
    ts_validation = hatyan.readts_dia(filename=file_data_predvali, station=current_station)
    ts_ext_validation = hatyan.readts_dia(filename=file_data_predvaliHWLW, station=current_station)
    hatyan.write_tsdia(ts=ts_prediction, station=current_station, vertref=vertref, filename='prediction_%im_%s.dia'%(times_step_pred,current_station))
    ts_ext_prediction1min = hatyan.calc_HWLW(ts=ts_prediction1min)
    hatyan.write_tsdia_HWLW(ts_ext=ts_ext_prediction1min, station=current_station, vertref=vertref, filename='prediction_HWLW_%im_%s.dia'%(times_step_pred, current_station))
    fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction, ts_ext=ts_ext_prediction1min, ts_ext_validation=ts_ext_validation)
    fig.savefig('prediction_%im_%s_HWLW'%(times_step_pred, current_station))
    fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction, ts_validation=ts_validation)
    fig.savefig('prediction_%im_%s_validation'%(times_step_pred, current_station))
    fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction, ts_validation=ts_measurements_group0)
    fig.savefig('prediction_%im_%s_measurements'%(times_step_pred, current_station))
    
    #plot and print HWLW statistics
    fig, ax = hatyan.plot_HWLW_validatestats(ts_ext=ts_ext_prediction1min, ts_ext_validation=ts_ext_validation)
    fig.savefig('HWLWstats_%im_%s_extvalidation'%(times_step_pred, current_station))
