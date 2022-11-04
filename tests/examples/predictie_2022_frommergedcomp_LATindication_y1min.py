# -*- coding: utf-8 -*-

"""
#Opm Anneke
15mrt21
de predictie_2022_frommergedcomp_all_1min.py is omgezet naar een predictie file om de LAT over een willekeurige 19 jaar te berekenen. Door Jelmer Veenstra is hiervoor in zijn predictiefiles een aanpassing gedaan. Deze heb ik omgezet in mijn predictiefile.

soort disclaimer van Jelmer:
Deze file maakt per station voor 19 jaar een astro predictie (per jaar, met knoopfactor in het midden van ieder jaar), op basis van een componentenset afgeleid over meestal 4 jaar (SA/SM van meestal 19 jaar).
Van deze predictie wordt per station en per jaar de minimum waterstand (en bijbehorende tijdstip) verzameld in de tabel min_vallist.
Van deze tabel wordt de minimumwaarde (en bijbehorende stationsnaam) aan de tabel min_vallist_allstats toegevoegd. Deze tabel wordt geprint naar LAT_indication.csv
Deze waardes geven een ruwe indicatie van LAT, maar kan niet zonder meer als 'de waarheid' worden beschouwd.
De (lengte van) de geanalyseerde periode is hierbij zeer bepalend (nu veelal 4 jaar maar dit verschilt per station, 19 jaar zou misschien beter zijn), SA/SM komen van een andere periode. Ook is de predictietijdstap (times_step_pred) bepalend voor de precisie.
 
"""

import os, sys
import datetime as dt
import pandas as pd
import hatyan

file_config = os.path.realpath(__file__)

dir_output, timer_start = hatyan.init_RWS(file_config, sys.argv, interactive_plots=True)

dir_base = os.path.abspath(os.path.join(file_config,os.pardir,os.pardir,os.pardir))
#dir_base = '/home/rikz/baak/hatyan2/'
#dir_base = '/home/ad.rws.nl/baaka01/hatyan/'

stats_all = ['AMLAHVN','BAALHK','BATH','BERGSDSWT','GATVBSLE','BRESKVHVN','BROUWHVSGT02','BROUWHVSGT08','CADZD','DELFZL','DENHDR','DENOVBTN','DORDT','EEMHVN','EEMSHVN','EURPFM','EURPHVN','GEULHVN','GOIDSOD','GOUDBG','HAGSBNDN','HANSWT','HARLGN','HARMSBG','HARTBG','HARTHVN','HARVT10','HEESBN','HELLVSS','HOEKVHLD','HOLWD','HUIBGT','IJMDBTHVN','IJMDSMPL','KATSBTN','KEIZVR','KORNWDZBTN','KRAMMSZWT','KRIMPADIJSL','KRIMPADLK','K13APFM','LAUWOG','LICHTELGRE','LITHDP','MAASSS','MAESLKRZZDE','MARLGT','MOERDK','NES','NIEUWSTZL','OOSTSDE04','OOSTSDE11','OOSTSDE14','OUDSD','OVLVHWT','PARKSS','RAKND','ROOMPBNN','ROOMPBTN','ROTTDM','ROZBSSNZDE','ROZBSSZZDE','SCHAARVDND','SCHEVNGN','SCHIERMNOG','SCHOONHVN','SINTANLHVSGR','SPIJKNSE','STAVNSE','STELLDBTN','SUURHBNZDE','TENNSHVN','TERNZN','TERSLNZE','TEXNZE','VLAARDGN','VLAKTVDRN','VLIELHVN','VLISSGN','VURN','WALSODN','WERKDBTN','WESTKPLE','WESTTSLG','WIERMGDN','YERSKE','ZALTBML','A12','AWGPFM','D15','F16','F3PFM','J6','K14PFM','L9PFM','NORTHCMRT','Q1']
stats_xfac0 = ['A12','D15','F16','F3PFM','J6','K13APFM','K14PFM','L9PFM','NORTHCMRT','Q1']
stats_anaperyear0 = ['A12','D15','F16','F3PFM','J6','K14PFM','NORTHCMRT']
stats_MSL = ['EURPFM','K13APFM','LICHTELGRE','A12','AWGPFM','D15','F16','F3PFM','J6','K14PFM','L9PFM','NORTHCMRT','Q1']
stats_all = ['A12','AWGPFM','D15','F16','F3PFM','J6','K14PFM','L9PFM','NORTHCMRT','Q1','IJMDSMPL','MAESLKRZZDE','K13APFM']
stats_all = ['HOLWD','HUIBGT','IJMDBTHVN','IJMDSMPL','KATSBTN','KEIZVR','KORNWDZBTN','KRAMMSZWT','KRIMPADIJSL','KRIMPADLK','K13APFM','LAUWOG','LICHTELGRE','LITHDP','MAASSS','MAESLKRZZDE','MARLGT','MOERDK','NES','NIEUWSTZL','OOSTSDE04','OOSTSDE11','OOSTSDE14','OUDSD','OVLVHWT','PARKSS','RAKND','ROOMPBNN','ROOMPBTN','ROTTDM','ROZBSSNZDE','ROZBSSZZDE','SCHAARVDND','SCHEVNGN','SCHIERMNOG','SCHOONHVN','SINTANLHVSGR','SPIJKNSE','STAVNSE','STELLDBTN','SUURHBNZDE','TENNSHVN','TERNZN','TERSLNZE','TEXNZE','VLAARDGN','VLAKTVDRN','VLIELHVN','VLISSGN','VURN','WALSODN','WERKDBTN','WESTKPLE','WESTTSLG','WIERMGDN','YERSKE','ZALTBML','A12','AWGPFM','D15','F16','F3PFM','J6','K14PFM','L9PFM','NORTHCMRT','Q1']

selected_stations = ['EURPFM']#stats_all

stats_noana = []
min_vallist_allstats = pd.DataFrame()

for current_station in selected_stations:
    print(f'processing station: {current_station}')
    
    #START OF STATION SETTINGS
    nodalfactors = True
    if current_station in stats_xfac0:
        xfac=False
    else:
        xfac=True
    if current_station in stats_anaperyear0:
        analysis_peryear=False
    else:
        analysis_peryear=True
    if current_station in ['F16','D15','F3PFM','K14PFM','MAESLKRZZDE','Q1','A12','AWGPFM','J6','L9PFM']:
        const_list = hatyan.get_const_list_hatyan('month') #21 const, potentially extended with component splitting (5 components) and SA+SM
    elif current_station in ['AMLAHVN']:
        const_list = hatyan.get_const_list_hatyan('halfyear') #88 const
    else:
        const_list = hatyan.get_const_list_hatyan('year') #94 const
    CS_comps = None
    if current_station in stats_MSL:
        vertref='MSL'
    else:
        vertref='NAP'

    file_data_comp0 = os.path.join(dir_base,'data','%s_ana.txt'%(current_station))
    file_data_compvali = os.path.join(dir_base,'data','%s_ana.txt'%(current_station))
    file_data_predvali = os.path.join(dir_base,'data','%s_pre.txt'%(current_station))
    file_data_predvaliHWLW = os.path.join(dir_base,'data','%s_ext.txt'%(current_station))
    
    if not os.path.exists(file_data_compvali):
        stats_noana.append(current_station)
        continue
    
    COMP_merged = hatyan.read_components(filename=file_data_comp0)
    min_vallist = pd.DataFrame()
    
    for year in range(2020,2039):
        times_ext_pred = [dt.datetime(year,1,1),dt.datetime(year+1,1,1)]
        times_step_pred = 10
   
        #ts_prediction = hatyan.prediction(comp=COMP_merged, nodalfactors=nodalfactors, xfac=xfac, fu_alltimes=False, times_ext=times_ext_pred, timestep_min=times_step_pred)
        ts_prediction = hatyan.prediction(comp=COMP_merged, nodalfactors=nodalfactors, xfac=xfac, fu_alltimes=False, times_ext=times_ext_pred, timestep_min=1)

        id_minvalue = ts_prediction['values'].argmin()
        ts_prediction_min = ts_prediction.iloc[[id_minvalue]]
        min_vallist = min_vallist.append(ts_prediction_min,ignore_index=True)
    print(min_vallist)
    min_vallist.to_csv('LAT_indication_19Y_%s.csv'%(current_station))

    ts_prediction_19Ymin = min_vallist['values'].min()
    min_1stat = pd.DataFrame({'station':[current_station],'values':[ts_prediction_19Ymin]})
    min_vallist_allstats = min_vallist_allstats.append(min_1stat)
   

print(min_vallist_allstats)
min_vallist_allstats.to_csv('LAT_indication.csv')
        
hatyan.exit_RWS(timer_start)

print('\nthese %i stations were requested for processing:\n%s'%(len(selected_stations),selected_stations))
print('\nthese %i stations were not processed because there is no ana/comp dataset available:\n%s'%(len(stats_noana),stats_noana))

