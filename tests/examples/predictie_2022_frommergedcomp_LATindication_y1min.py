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
dir_base = r'c:\DATA\hatyan_data_acceptancetests\predictie2019'

stats_all = ['AMLAHVN','BAALHK','BATH','BERGSDSWT','GATVBSLE','BRESKVHVN','BROUWHVSGT02','BROUWHVSGT08','CADZD','DELFZL','DENHDR','DENOVBTN','DORDT','EEMHVN','EEMSHVN','EURPFM','EURPHVN','GEULHVN','GOIDSOD','GOUDBG','HAGSBNDN','HANSWT','HARLGN','HARMSBG','HARTBG','HARTHVN','HARVT10','HEESBN','HELLVSS','HOEKVHLD','HOLWD','HUIBGT','IJMDBTHVN','IJMDSMPL','KATSBTN','KEIZVR','KORNWDZBTN','KRAMMSZWT','KRIMPADIJSL','KRIMPADLK','K13APFM','LAUWOG','LICHTELGRE','LITHDP','MAASSS','MAESLKRZZDE','MARLGT','MOERDK','NES','NIEUWSTZL','OOSTSDE04','OOSTSDE11','OOSTSDE14','OUDSD','OVLVHWT','PARKSS','RAKND','ROOMPBNN','ROOMPBTN','ROTTDM','ROZBSSNZDE','ROZBSSZZDE','SCHAARVDND','SCHEVNGN','SCHIERMNOG','SCHOONHVN','SINTANLHVSGR','SPIJKNSE','STAVNSE','STELLDBTN','SUURHBNZDE','TENNSHVN','TERNZN','TERSLNZE','TEXNZE','VLAARDGN','VLAKTVDRN','VLIELHVN','VLISSGN','VURN','WALSODN','WERKDBTN','WESTKPLE','WESTTSLG','WIERMGDN','YERSKE','ZALTBML','A12','AWGPFM','D15','F16','F3PFM','J6','K14PFM','L9PFM','NORTHCMRT','Q1']
stats_xfac0 = ['A12','D15','F16','F3PFM','J6','K13APFM','K14PFM','L9PFM','NORTHCMRT','Q1']
stats_all = ['A12','AWGPFM','D15','F16','F3PFM','J6','K14PFM','L9PFM','NORTHCMRT','Q1','IJMDSMPL','MAESLKRZZDE','K13APFM']
stats_all = ['HOLWD','HUIBGT','IJMDBTHVN','IJMDSMPL','KATSBTN','KEIZVR','KORNWDZBTN','KRAMMSZWT','KRIMPADIJSL','KRIMPADLK','K13APFM','LAUWOG','LICHTELGRE','LITHDP','MAASSS','MAESLKRZZDE','MARLGT','MOERDK','NES','NIEUWSTZL','OOSTSDE04','OOSTSDE11','OOSTSDE14','OUDSD','OVLVHWT','PARKSS','RAKND','ROOMPBNN','ROOMPBTN','ROTTDM','ROZBSSNZDE','ROZBSSZZDE','SCHAARVDND','SCHEVNGN','SCHIERMNOG','SCHOONHVN','SINTANLHVSGR','SPIJKNSE','STAVNSE','STELLDBTN','SUURHBNZDE','TENNSHVN','TERNZN','TERSLNZE','TEXNZE','VLAARDGN','VLAKTVDRN','VLIELHVN','VLISSGN','VURN','WALSODN','WERKDBTN','WESTKPLE','WESTTSLG','WIERMGDN','YERSKE','ZALTBML','A12','AWGPFM','D15','F16','F3PFM','J6','K14PFM','L9PFM','NORTHCMRT','Q1']

selected_stations = ['HOEKVHLD']#stats_all

stats_noana = []
vallist_allstats = pd.DataFrame(columns=['LAT','HAT'])


for current_station in selected_stations:
    print(f'processing station: {current_station}')
    
    #START OF STATION SETTINGS
    nodalfactors = True
    if current_station in stats_xfac0:
        xfac=False
    else:
        xfac=True

    file_data_comp0 = os.path.join(dir_base,'%s_ana.txt'%(current_station))
    if not os.path.exists(file_data_comp0):
        stats_noana.append(current_station)
        continue
    
    COMP_merged = hatyan.read_components(filename=file_data_comp0)
    vallist_allyears = pd.DataFrame()
    for year in range(2020,2039):
        times_pred_all = pd.date_range(start=dt.datetime(year,1,1), end=dt.datetime(year+1,1,1), freq='1min')
        ts_prediction = hatyan.prediction(comp=COMP_merged, nodalfactors=nodalfactors, xfac=xfac, fu_alltimes=False, times_pred_all=times_pred_all)
        
        vallist_allyears.loc[year,'min'] = ts_prediction['values'].min()
        vallist_allyears.loc[year,'max'] = ts_prediction['values'].max()
    #vallist_allyears.plot()
    #print(vallist_allyears)
    #vallist_allyears.to_csv('LAT_HAT_indication_19Y_%s.csv'%(current_station))
    
    vallist_allstats.loc[current_station,'LAT'] = vallist_allyears['min'].min()
    vallist_allstats.loc[current_station,'HAT'] = vallist_allyears['max'].max()


print(vallist_allstats)
vallist_allstats.to_csv('LAT_HAT_indication.csv')
        
hatyan.exit_RWS(timer_start)

print('\nthese %i stations were requested for processing:\n%s'%(len(selected_stations),selected_stations))
print('\nthese %i stations were not processed because there is no ana/comp dataset available:\n%s'%(len(stats_noana),stats_noana))

