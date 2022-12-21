# -*- coding: utf-8 -*-

"""
Gebruik van componentenfile om LAT/HAT te berekenen.

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
    hatyan_settings = hatyan.HatyanSettings(nodalfactors=nodalfactors, xfac=xfac, fu_alltimes=True)

    file_data_comp0 = os.path.join(dir_base,'%s_ana.txt'%(current_station))
    if not os.path.exists(file_data_comp0):
        stats_noana.append(current_station)
        continue
    
    COMP_merged = hatyan.read_components(filename=file_data_comp0) # make sure A0 is replaced by the slotgemiddelde (or do it afterwards)
    
    LAT, HAT = hatyan.calc_LAT_HAT_fromcomponents(comp=COMP_merged,hatyan_settings=hatyan_settings)
    
    vallist_allstats.loc[current_station,'LAT'] = LAT
    vallist_allstats.loc[current_station,'HAT'] = HAT


print(vallist_allstats)
vallist_allstats.to_csv('LAT_HAT_indication.csv')
        
hatyan.exit_RWS(timer_start)

print('\nthese %i stations were requested for processing:\n%s'%(len(selected_stations),selected_stations))
print('\nthese %i stations were not processed because there is no ana/comp dataset available:\n%s'%(len(stats_noana),stats_noana))

