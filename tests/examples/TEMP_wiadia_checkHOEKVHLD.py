# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 18:33:22 2021

@author: veenstra
"""

import os
import hatyan
hatyan.close('all')

file_config = os.path.realpath(__file__)
dir_output, timer_start = hatyan.init_RWS(file_config, interactive_plots=False)

dir_testdata = 'C:\\DATA\\hatyan_data_acceptancetests'

current_station = 'HOEKVHLD'

vertref='NAP'

file_dia_wl = os.path.join(dir_testdata,'other','diawia_%s_astro_tijdreeks.dia'%(current_station))
file_dia_ext = os.path.join(dir_testdata,'other','diawia_%s_astro_extremen.dia'%(current_station))

for file_dia in [file_dia_wl,file_dia_ext]:
    file_wia = file_dia.replace('.dia','.wia')
    file_dia_out = file_dia.replace('.dia','_out.dia')
    file_wia_out = file_dia.replace('.dia','_out.wia')
    ts_dia = hatyan.readts_dia(filename=file_dia, station=current_station)
    ts_wia = hatyan.readts_dia(filename=file_wia, station=current_station)
    assert (ts_dia==ts_wia).all().all() #check if wia and dia input is equal
    
    #write to files
    if 'HWLWcode' in ts_dia.columns:
        hatyan.write_tsdia_HWLW(ts_ext=ts_dia, station=current_station, vertref=vertref, filename=file_dia_out)
        hatyan.write_tsdia_HWLW(ts_ext=ts_wia, station=current_station, vertref=vertref, filename=file_wia_out, headerformat='wia')
    else:
        hatyan.write_tsdia(ts=ts_dia, station=current_station, vertref=vertref, filename=file_dia_out)
        hatyan.write_tsdia(ts=ts_wia, station=current_station, vertref=vertref, filename=file_wia_out, headerformat='wia')
    
    #read from new files
    ts_dia_new = hatyan.readts_dia(filename=file_dia_out, station=current_station)
    ts_wia_new = hatyan.readts_dia(filename=file_wia_out, station=current_station)
    assert (ts_dia==ts_dia_new).all().all() #check if wia and dia input is equal
    assert (ts_wia==ts_wia_new).all().all() #check if wia and dia input is equal
    
    #remove files
    os.remove(file_dia_out)
    os.remove(file_wia_out)
    
hatyan.exit_RWS(timer_start) #provides footer to outputfile when calling this script with python






