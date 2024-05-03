# -*- coding: utf-8 -*-
"""
predictie_2019_19Ycomp4Ydia.py
hatyan master configfile
voor alle stations indien mogelijk:
    - lees opgegeven data in
    - analyse van de ingelezen data, eventueel met componentensplitsing
    - combineren met SA+SM uit analyseresultatenbestand
    - predictie maken

"""

import os
import pandas as pd
import hatyan
hatyan.close('all')

dir_testdata = 'C:\\DATA\\hatyan\\tests'

selected_stations = ['VLISSGN']

file_slotgemiddelden = os.path.join(dir_testdata,'data_unitsystemtests','_slotgemiddelden_predictie2019.txt')
stations_slotgem = pd.read_csv(file_slotgemiddelden, names=['slotgemiddelde'], comment='#', sep="\\s+")

for current_station in selected_stations:
    print(f'current_station: {current_station}')
    
    # station settings
    nodalfactors = True
    xfac=True
    analysis_perperiod='Y'
    const_list = hatyan.get_const_list_hatyan('year') # 94 constituents
    CS_comps = None # component splitting
    
    # file pattern for multiple diafiles. Use ? instead of * to avoid matching of obs19.txt
    file_data_comp0 = os.path.join(dir_testdata,'data_unitsystemtests',f'{current_station}_obs?.txt')
        
    file_data_comp1 = os.path.join(dir_testdata,'data_unitsystemtests',f'{current_station}_ana.txt')
    
    file_data_compvali = os.path.join(dir_testdata,'data_unitsystemtests',f'{current_station}_ana.txt')
    
    times_pred = slice("2019-01-01","2020-01-01", 10)
    times_pred_1min = slice("2019-01-01","2020-01-01", 1)

    file_data_predvali = os.path.join(dir_testdata,'data_unitsystemtests',f'{current_station}_pre.txt')
    file_data_predvaliHWLW = os.path.join(dir_testdata,'data_unitsystemtests',f'{current_station}_ext.txt')
    
    ts_measurements_group0 = hatyan.read_dia(filename=file_data_comp0, station=current_station)
    
    comp_frommeasurements_avg_group0, comp_frommeasurements_all_group0 = hatyan.analysis(ts=ts_measurements_group0, const_list=const_list, nodalfactors=nodalfactors, xfac=xfac, fu_alltimes=False, analysis_perperiod=analysis_perperiod, return_allperiods=True)

    fig,(ax1,ax2) = hatyan.plot_components(comp_frommeasurements_avg_group0, comp_allperiods=comp_frommeasurements_all_group0)
    fig.savefig('components_%s_4Y.png'%(current_station))
    hatyan.write_components(comp_frommeasurements_avg_group0, filename=f'components_{current_station}_4Y.txt')
    
    comp_fromfile_group1 = hatyan.read_components(filename=file_data_comp1)
    assert comp_fromfile_group1.attrs["xfac"] == xfac

    #merge component groups (SA/SM from 19Y, rest from 4Y)
    COMP_merged = hatyan.merge_componentgroups(comp_main=comp_frommeasurements_avg_group0, comp_sec=comp_fromfile_group1, comp_sec_list=['SA','SM'])
    #replace A0 amplitude (middenstand) by slotgemiddelde
    if current_station in stations_slotgem.index.tolist():
        COMP_merged.loc['A0','A'] = stations_slotgem.loc[current_station,'slotgemiddelde']

    COMP_validation = hatyan.read_components(filename=file_data_compvali)
    assert COMP_validation.attrs["xfac"] == xfac
    fig, (ax1,ax2) = hatyan.plot_components(COMP_merged, comp_validation=COMP_validation)
    fig.savefig('components_%s_merged.png'%(current_station))
    hatyan.write_components(COMP_merged, filename='components_%s_merged.txt'%(current_station))
    
    #prediction and validation
    ts_prediction = hatyan.prediction(comp=COMP_merged, times=times_pred)
    ts_prediction1min = hatyan.prediction(comp=COMP_merged, times=times_pred_1min)
    ts_validation = hatyan.read_dia(filename=file_data_predvali, station=current_station)
    ts_ext_validation = hatyan.read_dia(filename=file_data_predvaliHWLW, station=current_station)
    hatyan.write_dia(ts=ts_prediction, filename='prediction_%im_%s.dia'%(times_pred.step,current_station))
    ts_ext_prediction1min = hatyan.calc_HWLW(ts=ts_prediction1min)
    hatyan.write_dia(ts=ts_ext_prediction1min, filename='prediction_HWLW_%im_%s.dia'%(times_pred.step, current_station))
    fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction, ts_ext=ts_ext_prediction1min, ts_ext_validation=ts_ext_validation)
    fig.savefig('prediction_%im_%s_HWLW'%(times_pred.step, current_station))
    fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction, ts_validation=ts_validation)
    fig.savefig('prediction_%im_%s_validation'%(times_pred.step, current_station))
    fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction, ts_validation=ts_measurements_group0)
    fig.savefig('prediction_%im_%s_measurements'%(times_pred.step, current_station))
    
    #plot and print HWLW statistics
    fig, ax = hatyan.plot_HWLW_validatestats(ts_ext=ts_ext_prediction1min, ts_ext_validation=ts_ext_validation)
    fig.savefig('HWLWstats_%im_%s_extvalidation'%(times_pred.step, current_station))
