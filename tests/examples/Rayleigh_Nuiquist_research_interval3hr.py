# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 11:36:14 2022

@author: veenstra
"""


import os
import datetime as dt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')
from matplotlib import cm
import hatyan


#defining a list of the components to be analysed (can also be 'half_year' and others, 'year' contains 94 components and the mean H0)
const_list = hatyan.get_const_list_hatyan('year') #['A0','M2','S2','M4'] # 
#const_list = ['3MS2','3M2S10']
station_list = ['WIERMWD1','WIERMGDN','WESTTSLG','UITHZWD1','TEXNZE','TERSLNZE','SCHIERMNOG','NES','LAUWOG','HUIBGT','HOLWD','HARLGN','EEMSHVN','DENHDR','DELFZL']
station_list = ['UITHZWD1','WIERMWD1','DENHDR','WESTTSLG','TEXNZE','HARLGN','TERSLNZE','HOLWD','WIERMGDN','NES','HUIBGT','SCHIERMNOG','LAUWOG','EEMSHVN','DELFZL'][::-1] #sorted on M2 amplitude
station_list = ['DENHDR']
year_list = list(range(1891,2021)) #list(range(1891,2021))
year_list = [1955]

#const_list = const_list[:48] #max 4maal daags en ook geen S4
drop_list = ['S4','3M2S10','2SM6','4M2S12']
for const in drop_list:
    const_list.remove(const)
#comp_close = ['3M2S10','3MS2','2SM6','2SM2','4MS4','4M2S12']

metadata = pd.read_pickle(r'p:\11208031-010-kenmerkende-waarden-k\work\measurements_wl_18700101_20220101\meta_DENHDR_measwl.pkl')

for station in station_list:
    
    for year in year_list:
        print(f'processing {station} for {year}')
        url_dataraw = 'http://watersysteemdata.deltares.nl/thredds/fileServer/watersysteemdata/Wadden/ddl/raw/waterhoogte/'
        file_csv = url_dataraw+f'{station}_OW_WATHTE_NVT_NAP_{year}_ddl_wq.csv'
        
        data_pd = pd.read_csv(file_csv,sep=';',parse_dates=['tijdstip'])
        nyears = 1#1 #TODO: 1 year gives a lot of rayleigh failures, but 3 are clearly outliers (almost 0 difference), is criterion computed correctly?
        for addyear in range(1,nyears):
            data_pd_add = pd.read_csv(file_csv.replace(f'_{year}_',f'_{year+addyear}_'),sep=';',parse_dates=['tijdstip'])
            data_pd = data_pd.append(data_pd_add).sort_values('tijdstip')
        #print(data_pd)
        
        bool_invalid = (data_pd['kwaliteitswaarde.code']!=0) | (data_pd['numeriekewaarde']>340) | (data_pd['numeriekewaarde']<-340) #TODO: even kijken of threshold nog nodig is -> hij werkt voor sommige gevallen niet zonder
        data_pd.loc[bool_invalid,'numeriekewaarde'] = np.nan
        ts_meas_raw = pd.DataFrame({'values':data_pd['numeriekewaarde'].values/100},index=data_pd['tijdstip'].dt.tz_localize(None)) #read dataset and convert to DataFrame with correct format #TODO: tijdzone is MET (UTC+1), ook van de resulterende getijcomponenten. Is dat wenselijk?
        ts_meas_raw = ts_meas_raw.sort_index(axis='index') #sort dataset on times (not necessary but easy for debugging)
        bool_duplicated_index = ts_meas_raw.index.duplicated()
        if bool_duplicated_index.sum()>0:
            ts_meas_raw = ts_meas_raw[~bool_duplicated_index] #remove duplicate timesteps if they are present
        
        stats = hatyan.Timeseries_Statistics(ts=ts_meas_raw)
        
        print(hatyan.check_ts(ts_meas_raw))
        
        comp_frommeas = hatyan.get_components_from_ts(ts=ts_meas_raw, const_list=const_list, nodalfactors=True, xfac=True, return_allyears=False, fu_alltimes=True, analysis_peryear=False)
        ts_pred = hatyan.prediction(comp=comp_frommeas,times_ext=[ts_meas_raw.index.min(),ts_meas_raw.index.max()],timestep_min=20) 
        fig,(ax1,ax2) = hatyan.plot_timeseries(ts=ts_meas_raw,ts_validation=ts_pred)
        ax2.set_ylim(-1,1)
        
