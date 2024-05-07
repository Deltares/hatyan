# -*- coding: utf-8 -*-
"""
predictie_2019_frommergedcomp.py
hatyan master configfile
voor alle stations indien mogelijk:
    - inlezen analyseresultatenbestand
    - predictie maken

"""

import os
import datetime as dt
import hatyan

dir_testdata = 'C:\\DATA\\hatyan_data_acceptancetests'

stats_all = ['ABDN','AMLAHVN','BAALHK','BATH','BERGSDSWT','BORSSLE','BOURNMH','BRESKS','BROUWHVSGT02','BROUWHVSGT08','CADZD','CROMR','CUXHVN','DELFZL','DENHDR','DENOVBTN','DEVPT','DORDT','DOVR','EEMHVN','EEMSHVN','EURPFM','EURPHVN','FELSWE','FISHGD','GEULHVN','GOIDSOD','GOUDBG','HAGSBNDN','HANSWT','HARLGN','HARMSBG','HARTBG','HARTHVN','HARVT10','HEESBN','HELLVSS','HOEKVHLD','HOLWD','HUIBGT','IJMDBTHVN','IJMDSMPL','IMMHM','KATSBTN','KEIZVR','KINLBVE','KORNWDZBTN','KRAMMSZWT','KRIMPADIJSL','KRIMPADLK','K13APFM','LAUWOG','LEITH','LICHTELGRE','LITHDP','LLANDNO','LOWST','MAASMSMPL','MAASSS','MAESLKRZZDE','MARLGT','MOERDK','NES','NEWHVN','NEWLN','NIEUWSTZL','NORTHSS','OOSTSDE04','OOSTSDE11','OOSTSDE14','OUDSD','OVLVHWT','PARKSS','PETTZD','PORTSMH','RAKND','ROOMPBNN','ROOMPBTN','ROTTDM','ROZBSSNZDE','ROZBSSZZDE','SCHAARVDND','SCHEVNGN','SCHIERMNOG','SCHOONHVN','SHEERNS','SINTANLHVSGR','SPIJKNSE','STAVNSE','STELLDBTN','STORNWY','SUURHBNZDE','TENNSHVN','TERNZN','TERSLNZE','TEXNZE','VLAARDGN','VLAKTVDRN','VLIELHVN','VLISSGN','VURN','WALSODN','WERKDBTN','WESTKPLE','WESTTSLG','WEYMH','WHITBY','WICK','WIERMGDN','YERSKE','ZALTBML','A12','AUKFPFM','AWGPFM','D15','F16','F3PFM','J6','K14PFM','L9PFM','NORTHCMRT','Q1']
stats_xfac0 = ['A12','ABDN','AUKFPFM','BOURNMH','CROMR','CUXHVN','D15','DEVPT','DOVR','F16','F3PFM','FELSWE','FISHGD','IMMHM','J6','K13APFM','K14PFM','KINLBVE','L9PFM','LEITH','LLANDNO','LOWST','NEWHVN','NEWLN','NORTHCMRT','NORTHSS','PORTSMH','Q1','SHEERNS','STORNWY','WEYMH','WHITBY','WICK']

#selected_stations = stats_all
# TODO: temporarily disabled AUKFPFM because of https://github.com/Deltares/hatyan/issues/197
selected_stations = ['CADZD','DORDT','ABDN']#,'AUKFPFM']

stats_noana = []

for current_station in selected_stations:
    print(f'current_station: {current_station}')
    
    # station settings
    nodalfactors = True
    if current_station in stats_xfac0:
        xfac=False
    else:
        xfac=True
    analysis_perperiod='Y'
    if current_station in ['D15','F3PFM','K14PFM','MAESLKRZZDE','Q1','A12','AWGPFM','F16','J6','L9PFM']:
        const_list = hatyan.get_const_list_hatyan('month') #21 const, potentially extended with component splitting (5 components) and SA+SM
    elif current_station in ['AMLAHVN']:
        const_list = hatyan.get_const_list_hatyan('halfyear') #88 const
    else:
        const_list = hatyan.get_const_list_hatyan('year') #94 const

    file_data_comp0 = os.path.join(dir_testdata,'predictie2019','%s_ana.txt'%(current_station))

    file_data_compvali = os.path.join(dir_testdata,'predictie2019','%s_ana.txt'%(current_station))
    
    times_pred = slice(dt.datetime(2019,1,1),dt.datetime(2020,1,1), "10min")
    
    file_data_predvali = os.path.join(dir_testdata,'predictie2019','%s_pre.txt'%(current_station))
    file_data_predvaliHWLW = os.path.join(dir_testdata,'predictie2019','%s_ext.txt'%(current_station))
    
    if not os.path.exists(file_data_compvali):
        stats_noana.append(current_station)
        continue
    
    #component groups
    COMP_merged = hatyan.read_components(filename=file_data_comp0)
    assert COMP_merged.attrs["xfac"] == xfac
    
    #prediction and validation
    ts_prediction = hatyan.prediction(comp=COMP_merged, times=times_pred)
    ts_validation = hatyan.read_dia(filename=file_data_predvali, station=current_station)
    #ts_ext_validation = hatyan.read_dia(filename=file_data_predvaliHWLW, station=current_station)
    hatyan.write_dia(ts=ts_prediction, filename='prediction_%s_%s.dia'%(times_pred.step,current_station))
    #ts_ext_prediction = hatyan.calc_HWLW(ts=ts_prediction)
    #hatyan.write_dia(ts=ts_ext_prediction, filename='prediction_HWLW_%s_%s.dia'%(times_pred.step, current_station))
    #fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction, ts_ext=ts_ext_prediction, ts_ext_validation=ts_ext_validation)
    #fig.savefig('prediction_%s_%s_HWLW'%(times_pred.step, current_station))
    fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction, ts_validation=ts_validation)
    fig.savefig('prediction_%s_%s_validation'%(times_pred.step, current_station))

print('\nthese %i stations were requested for processing:\n%s'%(len(selected_stations),selected_stations))
print('\nthese %i stations were not processed because there is no ana/comp dataset available:\n%s'%(len(stats_noana),stats_noana))
