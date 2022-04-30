# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 14:17:13 2022

@author: veenstra
"""

import sys
import os
import sys
import pandas as pd
import datetime as dt
import numpy as np
import hatyan # beschikbaar via https://github.com/Deltares/hatyan
import matplotlib.pyplot as plt
plt.close('all')

#TODO: apply to all measurements: remove QC==99 (always, or maybe make nans?), crop_timeseries (when applicable), NAP2005 correction?, SLR trend correctie voor overschrijdingsfrequenties en evt ook voor andere KW?
#TODO: when to deliver data for Anneke and Robert and for which stations? (stationslijst opvragen)
#TODO: move all parts to hatyan.kenmerkendewaarden.*, maybe also the stuff in hatyan/overschrijding.py (and include license header) >> indeed put it in hatyan or not?
#TODO: add tidal indicators (LAT etc)


tstart_dt_DDL = dt.datetime(1870,1,1) #1870,1,1 for measall folder
tstop_dt_DDL = dt.datetime(2022,1,1)
tzone_DLL = 'UTC+01:00' #'UTC+00:00' for GMT and 'UTC+01:00' for MET
tstart_dt = dt.datetime(2001,1,1)
tstop_dt = dt.datetime(2011,1,1)
reproduce_2011_olddata = False #TODO: difference in gemgetijkromme (summary figure) for HARVT10 2011.0 (spnp/sp/np lines), probably because of duplicate values in measwl DDL
NAP2005correction = True
if ((tstop_dt.year-tstart_dt.year)==10) & (tstop_dt.month==tstop_dt.day==tstart_dt.month==tstart_dt.day==1):
    year_slotgem = tstop_dt.year
    if reproduce_2011_olddata:
        if not ((tstart_dt==dt.datetime(2001,1,1)) & (tstop_dt==dt.datetime(2011,1,1))):
            raise Exception('INVALID DATES WITH reproduce_2011_olddata')
        else:
            year_slotgem = '2011_olddata'
else:
    year_slotgem = 'invalid'
print(f'year_slotgem: {year_slotgem}')

dir_base = r'p:\11208031-010-kenmerkende-waarden-k\work'
if reproduce_2011_olddata:
    dir_meas = r'p:\11208031-010-kenmerkende-waarden-k\work\measurements_wl_20010101_20110101_olddata'
else:
    dir_meas = os.path.join(dir_base,'measurements_wl_20000101_20220101')
    dir_meas_alldata = os.path.join(dir_base,'measurements_wl_18700101_20220101')
dir_meas_DDL = os.path.join(dir_base,f"measurements_wl_{tstart_dt_DDL.strftime('%Y%m%d')}_{tstop_dt_DDL.strftime('%Y%m%d')}")
if not os.path.exists(dir_meas_DDL):
    os.mkdir(dir_meas_DDL)
dir_havget = os.path.join(dir_base,f'out_havengetallen_{year_slotgem}')
if not os.path.exists(dir_havget):
    os.mkdir(dir_havget)
dir_slotgem = os.path.join(dir_base,'out_slotgem')
if not os.path.exists(dir_slotgem):
    os.mkdir(dir_slotgem)
dir_gemgetij = os.path.join(dir_base,f'out_gemgetij_{year_slotgem}')
if not os.path.exists(dir_gemgetij):
    os.mkdir(dir_gemgetij)
dir_overschrijding = os.path.join(dir_base,'out_overschrijding')
if not os.path.exists(dir_overschrijding):
    os.mkdir(dir_overschrijding)

catalog_dict = hatyan.get_DDL_catalog(catalog_extrainfo=['WaardeBepalingsmethoden','MeetApparaten','Typeringen'])
cat_locatielijst = catalog_dict['LocatieLijst']#.set_index('Locatie_MessageID',drop=True)
cat_locatielijst.to_pickle(os.path.join(dir_meas_DDL,'catalog_lokatielijst.pkl'))
cat_aquometadatalijst_ext, cat_locatielijst_ext = hatyan.get_DDL_stationmetasubset(catalog_dict=catalog_dict,station_dict=None,meta_dict ={'Grootheid.Code':'WATHTE','Groepering.Code':'GETETM2'})
K13APFM_entry = cat_locatielijst.loc[cat_locatielijst['Code']=='K13APFM'].set_index('Locatie_MessageID',drop=True) #K13a does not have extremes, so is manually added to the interest-list
cat_locatielijst_ext = cat_locatielijst_ext.append(K13APFM_entry)
cat_locatielijst_ext['RDx'],cat_locatielijst_ext['RDy'] = hatyan.convert_coordinates(coordx_in=cat_locatielijst_ext['X'].values, coordy_in=cat_locatielijst_ext['Y'].values, epsg_in=int(cat_locatielijst_ext['Coordinatenstelsel'].iloc[0]),epsg_out=28992)
cat_locatielijst_ext_codeidx = cat_locatielijst_ext.reset_index(drop=False).set_index('Code',drop=False)

""" #TODO: maybe add these coastal stations
Holwerd
Stortemelk
Eierland
Texel Westgat (Texel Noordzee is er wel)
Oostoever
Petten (Petten zuid is er wel)
IJmuiden zuidelijk havenhoofd / IJmuiden semafoor / IJmuiden Noordersluis (IJmuiden buitenhaven is er wel)
Noordwijk meetpost
Katwijk / Katwijk paal
Brouwershavensche Gat 05 / Brouwershavensche Gat, punt 02
Oosterschelde 04 / 09 / 10 / 11 / 12 / 13 / 14 / 15
Oranjezon
Oostkapelle
Breskens
Perkpolder Walsoorden
"""

#stat_name_list = ['BATH','DELFZIJL','DEN HELDER','DORDRECHT','EEMSHAVEN','EURO PLATFORM','HANSWEERT','HARINGVLIETSLUIZEN','HARLINGEN','HOEK VAN HOLLAND','HUIBERTGAT','IJMUIDEN','KORNWERDERZAND','LAUWERSOOG','ROOMPOT BUITEN','ROTTERDAM','SCHEVENINGEN','STAVENISSE','TERNEUZEN','VLISSINGEN','WEST-TERSCHELLING'] # lijst AB
stat_name_list = ['Terneuzen','Bath','Hansweert','Vlissingen','Bergse Diepsluis west','Krammersluizen west','Stavenisse','Roompot binnen','Cadzand','Westkapelle','Roompot buiten','Brouwershavensche Gat 08','Haringvliet 10','Hoek van Holland','Scheveningen','IJmuiden buitenhaven','Petten zuid','Den Helder','Texel Noordzee','Terschelling Noordzee','Wierumergronden','Huibertgat','Oudeschild','Vlieland haven','West-Terschelling','Nes','Schiermonnikoog','Den Oever buiten','Kornwerderzand buiten','Harlingen','Lauwersoog','Eemshaven','Delfzijl','Nieuwe Statenzijl','Lichteiland Goeree','Euro platform','K13a platform'] + ['Dordrecht','Stellendam Buiten','Rotterdam']#Dillingh 2013, aangevuld met 3 stations AB
stat_list = []
for stat_name in stat_name_list:
    bool_isstation = cat_locatielijst_ext_codeidx['Naam'].str.contains(stat_name,case=False)
    if not bool_isstation.sum()==1:
        raise Exception(f'station name {stat_name} found {bool_isstation.sum()} times, should be 1.')
    stat_list.append(cat_locatielijst_ext_codeidx.loc[bool_isstation,'Code'].iloc[0])
    #print(f'{stat_name:30s}: {bool_isstation.sum()}')
#stat_list = ['BATH','DELFZL','DENHDR','DORDT','EEMSHVN','EURPFM','HANSWT','STELLDBTN','HARLGN','HOEKVHLD','HUIBGT','IJMDBTHVN','KORNWDZBTN','LAUWOG','ROOMPBTN','ROTTDM','SCHEVNGN','STAVNSE','TERNZN','VLISSGN','WESTTSLG'] # lijst AB vertaald naar DONAR namen
#stat_list = ['HOEKVHLD','HARVT10','VLISSGN']

M2_period_timedelta = pd.Timedelta(hours=hatyan.get_schureman_freqs(['M2']).loc['M2','period [hr]'])


def nap2005_correction(data_pd,current_station):
    #NAP correction for dates before 1-1-2005
    #TODO: make this flexible per station, where to get the data or is the RWS data already corrected for it? Also does it matter? for havengetallen it makes a slight difference so yes. For gemgetijkromme it only makes a difference for spring/doodtij, because A0 is in componentlist (which it should not be probably?) (now only applied at gemgetij en havengetallen)
    print('applying NAP2005 correction')
    data_pd_corr = data_pd.copy()
    before2005bool = data_pd_corr.index<dt.datetime(2005,1,1)
    dict_correct_nap2005 = {'HOEKVHLD':-0.0277,
                            'HARVT10':-0.0210,
                            'VLISSGN':-0.0297}
    if not current_station in dict_correct_nap2005.keys():
        raise Exception('ERROR nap2005 correction not implemented for this station')

    correct_value = dict_correct_nap2005[current_station]
    data_pd_corr.loc[before2005bool,'values'] = data_pd_corr.loc[before2005bool,'values']+correct_value
    
    return data_pd_corr


### RETRIEVE DATA FROM DDL AND WRITE TO PICKLE
for current_station in []:#stat_list:
    file_wl_pkl = os.path.join(dir_meas_DDL,f"{current_station}_measwl.pkl")
    file_wlmeta_pkl = os.path.join(dir_meas_DDL,f"meta_{current_station}_measwl.pkl")
    
    station_dict = cat_locatielijst_ext_codeidx.loc[current_station,['Locatie_MessageID','X','Y','Naam','Code']]
    
    allow_multipleresultsfor = ['WaardeBepalingsmethode'] # necessary for retrieving very long timeseries
    
    #retrieving waterlevels
    if os.path.exists(file_wl_pkl):
        print(f'measwl data for {current_station} already available in {os.path.basename(dir_meas_DDL)}')
    elif 0: #tstart_dt_DDL < dt.datetime(2000,1,1):
        print('skipping this station since the period is too long')
    else:
        print(f'retrieving measwl data from DDL for {current_station} to {os.path.basename(dir_meas_DDL)}')
        request_output = hatyan.get_DDL_data(station_dict=station_dict,tstart_dt=tstart_dt_DDL,tstop_dt=tstop_dt_DDL,tzone=tzone_DLL, allow_multipleresultsfor=allow_multipleresultsfor,
                                             meta_dict={'Grootheid.Code':'WATHTE','Groepering.Code':'NVT',
                                                        'Hoedanigheid.Code':'NAP',  # Hoedanigheid is necessary for eg EURPFM, where NAP and MSL values are available.
                                                        'MeetApparaat.Code':'127'}) # MeetApparaat.Code is necessary for IJMDBTHVN/ROOMPBTN, where also radar measurements are available (all other stations are vlotter and these stations also have all important data in vlotter)
                                            #meta_dict={'Grootheid.Code':'WATHTE','Groepering.Code':'NVT'} #ts_measwl
                                            #meta_dict={'Grootheid.Code':'WATHTE','Groepering.Code':'GETETM2'} #ts_measwlHWLW
                                            #meta_dict={'Groepering.Code':'GETETM2','Typering.Code':'GETETTPE'} #ts_measwlHWLWtype
                                            #meta_dict={'Grootheid.Code':'WATHTBRKD','Groepering.Code':'NVT'} #ts_astro
                                            #meta_dict={'Grootheid.Code':'WATHTBRKD','Groepering.Code':'GETETBRKD2'} #ts_astroHWLW
                                            #meta_dict={'Groepering.Code':'GETETBRKD2','Typering.Code':'GETETTPE'} #ts_astroHWLWtype
        if request_output is None:
            continue
        ts_meas_pd, metadata, stationdata = request_output
        ts_meas_pd['values'] = ts_meas_pd['values']/100 #convert from cm to m
        ts_meas_pd.to_pickle(file_wl_pkl)
        metadata.to_pickle(file_wlmeta_pkl)
    
    #retrieving measured extremes
    if os.path.exists(file_wl_pkl.replace('_measwl','_measext')):
        print(f'measext data for {current_station} already available in {os.path.basename(dir_meas_DDL)}')
    else:
        print(f'retrieving measext data from DDL for {current_station} to {os.path.basename(dir_meas_DDL)}')
        request_output_extval = hatyan.get_DDL_data(station_dict=station_dict,tstart_dt=tstart_dt_DDL,tstop_dt=tstop_dt_DDL,tzone=tzone_DLL, allow_multipleresultsfor=allow_multipleresultsfor,
                                                    meta_dict={'Grootheid.Code':'WATHTE','Groepering.Code':'GETETM2'})#,'MeetApparaat.Code':'127'}) #ts_measwlHWLW # TODO: MeetApparaat is necessary for IJMBTHVN, maybe remove if servicedesk has resolved this probable Vlotter/Radar issue (gemeld op 28-4-2022) (or add Hoedanigheid.Code, alles is toch NAP)
        request_output_exttyp = hatyan.get_DDL_data(station_dict=station_dict,tstart_dt=tstart_dt_DDL,tstop_dt=tstop_dt_DDL,tzone=tzone_DLL, allow_multipleresultsfor=allow_multipleresultsfor,
                                                    meta_dict={'Groepering.Code':'GETETM2','Typering.Code':'GETETTPE'})#,'MeetApparaat.Code':'127'}) #ts_measwlHWLWtype
        if request_output_extval is None:
            continue
        ts_meas_ext_pd, metadata, stationdata = request_output_extval
        ts_meas_exttyp_pd, metadata2, dummy = request_output_exttyp
        ts_meas_ext_pd['values'] = ts_meas_ext_pd['values']/100 #convert from cm to m
        if ts_meas_exttyp_pd['values'].isnull().any(): #TODO: remove this exception for SCHEVNGN after DDL exttype data is fixed
            continue
        ts_meas_ext_pd = hatyan.convert_HWLWstr2num(ts_meas_ext_pd,ts_meas_exttyp_pd)
        ts_meas_ext_pd.to_pickle(file_wl_pkl.replace('_measwl','_measext'))
        metadata.to_pickle(file_wlmeta_pkl.replace('_measwl','_measext'))






### LOAD DATA FROM PICKLE plot and do checks
#TODO: getijanalyse+predictie?
#TODO: visually check availability (start/stop/gaps/aggers) of wl/ext, monthmean wl, outliers (nog niet gedaan voor hele periode, wel voor 2000-2022 (listAB+HARVT10):
#   IJMDBTHVN extremen missen vanaf 2018 want Radar ipv Vlotter (al gemeld op 28-4-2022)
#   Missende data vanaf 2000 (gemeld op 26-4):
#       BATH (2000-2020, measwl en measext, komt doordat er twee stations zijn genaamd Bath/BATH, dit is waarschijnlijk de realtime versie)
#       EURPFM (2000-2001, measwl en measext)
#       HOEKVHLD (2000-2006, 2013, 2019, measext)
#       ROOMPBTN (2016, measext)
#       ROTTDM (2013, 2019, measext)
#       SCHEVNGN (extrementype bevatten nans ipv strings als 'hoogwater', data is dus ongeldig)
#       STELLDBTN (geen data beschikbaar) >> geen ext data vanaf 2000
#   sterke outliers in tijdreeksen (na filtering QC=99, gemeld op 26-4): IJMDBTHVN/ROOMPBTN (2001) >> is niet meer zo na verwijderen Radar metingen
#TODO: report dubbelingen HARVT10 (2000-2022, al gedaan?) en andere stations (1900-2000), en EURPFM ext, zie data_summary.csv (er zijn ook dubbelingen met nan-waardes)
#TODO: report wl/ext missings in recent period 2000-2021 (vanuit data_summary)
#TODO: request gem hw/lw/wl bij Anneke voor alle jaren (gedaan op 28-04-2022)
#TODO: vergelijking yearmean wl/HW/LW met validatiedata Anneke (nu alleen beschikbaar voor HOEKVHLD en HARVT10, sowieso wl is nodig voor slotgemiddelde), it is clear in the HARVT10 figures that something is off for meanwl, dit gebeurt misschien ook bij andere stations met duplicate times in data_summary_filtered.xlsx (also check on nanvalues that are not nan in validationdata, this points to missing data in DDL)
#TODO: als extremen missen evt zelf afleiden, maar is misschien niet zomaar gedaan? (wanneer met/zonder aggers?, calcHWLW verwacht redelijk constante tijdstap) >> HWLW numbering werkt sowieso niet heel goed met metingen blijkt nu
"""
#TODO: melden servicedesk data: zes duplicate timesteps in extremen aanwezig met gelijke waarden
                     values  QC         Status  HWLWcode
Tijdstip                                                
2012-12-31 09:35:00   -0.70   0  Gecontroleerd         2
2012-12-31 09:35:00   -0.70   0  Gecontroleerd         2
2012-12-31 15:44:00    0.95   0  Gecontroleerd         1
2012-12-31 15:44:00    0.95   0  Gecontroleerd         1
2012-12-31 20:50:00   -0.63   0  Gecontroleerd         2
2012-12-31 20:50:00   -0.63   0  Gecontroleerd         2
2013-01-01 04:04:00    1.40   0  Gecontroleerd         1
2013-01-01 04:04:00    1.40   0  Gecontroleerd         1
2013-01-01 09:34:00   -0.36   0  Gecontroleerd         2
2013-01-01 09:34:00   -0.36   0  Gecontroleerd         2
2013-01-01 16:26:00    1.39   0  Gecontroleerd         1
2013-01-01 16:26:00    1.39   0  Gecontroleerd         1
"""
"""
#gemeld op 28-4-2022 bij servicedesk data: Radar extremen IJMDBTHVN vanaf 2018 (dus missings)
import hatyan # "pip install hatyan"
station_dict_IJMDBTHVN = {'Locatie_MessageID': 20503,
                          'Coordinatenstelsel': '25831',
                          'X': 605633.035699228,
                          'Y': 5813598.03897256,
                          'Naam': 'IJmuiden buitenhaven',
                          'Code': 'IJMDBTHVN'}
request_output_extval = hatyan.get_DDL_data(station_dict=station_dict_IJMDBTHVN,tstart_dt=dt.datetime(2018,1,1),tstop_dt=dt.datetime(2018,12,31,23,50),tzone='UTC+01:00', allow_multipleresultsfor=['WaardeBepalingsmethode'],
                                            meta_dict={'Grootheid.Code':'WATHTE','Groepering.Code':'GETETM2'}) #HWLWvalues
                                            #meta_dict={'Groepering.Code':'GETETM2','Typering.Code':'GETETTPE'}) #HWLWtypes
#IMPROVEMENT (data): results in duplicate MeetApparaat, states that 'Radar' is used to derive waterlevels for extremes from somewhere in 2018 for IJMDBTHVN. This is probably not true since the measured waterlevels for that period are measured with 'Vlotter', just like extremes in all other years for this station and all other stations I checked (for period 1900-2022). This also goes for the GETETTPE
Result 1:
              'MeetApparaat': {'MeetApparaat.Code': '127', 'MeetApparaat.Omschrijving': 'Vlotter'},
Result 2:
              'MeetApparaat': {'MeetApparaat.Code': '109', 'MeetApparaat.Omschrijving': 'Radar'},

request_output_extval = hatyan.get_DDL_data(station_dict=station_dict_IJMDBTHVN,tstart_dt=dt.datetime(2018,1,1),tstop_dt=dt.datetime(2018,12,31,23,50),tzone='UTC+01:00', allow_multipleresultsfor=['WaardeBepalingsmethode'],
                                            meta_dict={'Grootheid.Code':'WATHTE','Groepering.Code':'GETETM2','MeetApparaat.Omschrijving':'Vlotter'}) #HWLWvalues
#IMPROVEMENT (DDL): now with 'MeetApparaat.Omschrijving':'Vlotter' included in query, but this seems not to be read by the DDL (same error, this should not happen I think)  (ook ROOMPBTN heeft alleen 'Vlotter' extremen, terwijl daar in 2000 ook radarmetingen aanwezig zijn in de gemeten waterstanden)

request_output_extval = hatyan.get_DDL_data(station_dict=station_dict_IJMDBTHVN,tstart_dt=dt.datetime(2018,1,1),tstop_dt=dt.datetime(2018,12,31,23,50),tzone='UTC+01:00', allow_multipleresultsfor=['WaardeBepalingsmethode'],
                                            meta_dict={'Grootheid.Code':'WATHTE','Groepering.Code':'GETETM2','MeetApparaat.Code':'127'}) #HWLWvalues
#now with 'MeetApparaat.Code':'127' included in query, this does the trick.
"""
#TODO: add to data_summary: std/mean monthavg/yearavg wl/HW/LW, stats splitsen voor HW/LW?
data_summary = pd.DataFrame(index=stat_list)
for current_station in []:#stat_list:
    print(f'checking data for {current_station}')
    
    #add coordinates to data_summary
    data_summary.loc[current_station,['RDx','RDy']] = cat_locatielijst_ext_codeidx.loc[current_station,['RDx','RDy']]
    
    #load measwl data
    file_wl_pkl = os.path.join(dir_meas_alldata,f"{current_station}_measwl.pkl")
    if not os.path.exists(file_wl_pkl):
        data_summary.loc[current_station,'data_wl'] = False
        data_summary.loc[current_station,'data_ext'] = False
        continue
    data_summary.loc[current_station,'data_wl'] = True
    ts_meas_pd = pd.read_pickle(file_wl_pkl)
    ts_meas_pd = ts_meas_pd[['values','QC']] # reduces the memory consumption significantly
    if str(ts_meas_pd.index[0].tz) != 'Etc/GMT-1': #this means UTC+1
        raise Exception(f'measwl data for {current_station} is not in expected timezone (Etc/GMT-1): {ts_meas_pd.index[0].tz}')
    ts_meas_pd.index = ts_meas_pd.index.tz_localize(None)
    bool_99 = ts_meas_pd['QC']==99
    if bool_99.any(): #ts contains invalid values
        ts_meas_pd[bool_99] = np.nan
    data_summary.loc[current_station,'tstart_wl'] = ts_meas_pd.index[0]#.tz_localize(None)
    data_summary.loc[current_station,'tstop_wl'] = ts_meas_pd.index[-1]#.tz_localize(None)
    data_summary.loc[current_station,'nvals_wl'] = len(ts_meas_pd['values'])
    data_summary.loc[current_station,'#nans_wl'] = bool_99.sum()
    data_summary.loc[current_station,'min_wl'] = ts_meas_pd['values'].min()
    data_summary.loc[current_station,'max_wl'] = ts_meas_pd['values'].max()
    data_summary.loc[current_station,'std_wl'] = ts_meas_pd['values'].std()
    data_summary.loc[current_station,'mean_wl'] = ts_meas_pd['values'].mean()
    ts_meas_dupltimes = ts_meas_pd.index.duplicated()
    data_summary.loc[current_station,'dupltimes_wl'] = ts_meas_dupltimes.sum()
    #count #nans for duplicated times, happens at HARVT10/HUIBGT/STELLDBTN
    data_summary.loc[current_station,'#nans_dupltimes_wl'] = ts_meas_pd.loc[ts_meas_pd.index.duplicated(keep=False),'values'].isnull().sum()
    
    #calc #nan-values in recent period
    ts_meas_2000to202102a = ts_meas_pd.loc[~ts_meas_dupltimes,['values']].loc[dt.datetime(2000,1,1):dt.datetime(2021,2,1)]
    ts_meas_2000to202102b = pd.DataFrame({'values':ts_meas_pd.loc[~ts_meas_dupltimes,'values']},index=pd.date_range(start=dt.datetime(2000,1,1),end=dt.datetime(2021,2,1),freq='10min'))
    data_summary.loc[current_station,'#nans_2000to202102a_wl'] = ts_meas_2000to202102a['values'].isnull().sum()
    data_summary.loc[current_station,'#nans_2000to202102b_wl'] = ts_meas_2000to202102b['values'].isnull().sum()

    #load measext data
    file_ext_pkl = os.path.join(dir_meas_alldata,f"{current_station}_measext.pkl")
    if not os.path.exists(file_ext_pkl):
        data_summary.loc[current_station,'data_ext'] = False
    else:
        data_summary.loc[current_station,'data_ext'] = True
        ts_meas_ext_pd = pd.read_pickle(file_ext_pkl)
        if str(ts_meas_ext_pd.index[0].tz) != 'Etc/GMT-1': #this means UTC+1
            raise Exception(f'measext data for {current_station} is not in expected timezone (Etc/GMT-1): {ts_meas_ext_pd.index[0].tz}')
        ts_meas_ext_pd.index = ts_meas_ext_pd.index.tz_localize(None)
        ts_meas_ext_dupltimes = ts_meas_ext_pd.index.duplicated()
        data_summary.loc[current_station,'dupltimes_ext'] = ts_meas_ext_dupltimes.sum()
        data_summary.loc[current_station,'tstart_ext'] = ts_meas_ext_pd.index[0]#.tz_localize(None)
        data_summary.loc[current_station,'tstop_ext'] = ts_meas_ext_pd.index[-1]#.tz_localize(None)
        data_summary.loc[current_station,'nvals_ext'] = len(ts_meas_ext_pd['values'])
        data_summary.loc[current_station,'min_ext'] = ts_meas_ext_pd['values'].min()
        data_summary.loc[current_station,'max_ext'] = ts_meas_ext_pd['values'].max()
        data_summary.loc[current_station,'std_ext'] = ts_meas_ext_pd['values'].std()
        data_summary.loc[current_station,'mean_ext'] = ts_meas_ext_pd['values'].mean()
        if len(ts_meas_ext_pd['HWLWcode'].unique()) > 2:
            data_summary.loc[current_station,'aggers_ext'] = True
        else:
            data_summary.loc[current_station,'aggers_ext'] = False
        try:
            ts_meas_ext_2000to202102 = ts_meas_ext_pd.loc[(ts_meas_ext_pd.index>=dt.datetime(2000,1,1)) & (ts_meas_ext_pd.index<=dt.datetime(2021,2,1))]
            ts_meas_ext_2000to202102 = hatyan.calc_HWLWnumbering(ts_meas_ext_2000to202102, station=current_station) #station argument helpt bij 3 extra stations
            HWmissings = (ts_meas_ext_2000to202102.loc[ts_meas_ext_pd['HWLWcode']==1,'HWLWno'].diff().dropna()!=1).sum()
            data_summary.loc[current_station,'#HWgaps_2000to202102_ext'] = HWmissings
        except Exception as e: #"tidal wave numbering: HW/LW numbers not always increasing" and "zero-size array to reduction operation minimum which has no identity" #TODO: fix by calulate and providing station or corr_tideperiods argument? Or fix otherwise in hatyan (maybe under different project)
            print(f'ERROR: {e}')
    
    #calculate monthly mean
    mean_peryearmonth_long = ts_meas_pd.groupby(pd.PeriodIndex(ts_meas_pd.index, freq="M"))['values'].mean()
    mean_peryear_long = ts_meas_pd.groupby(pd.PeriodIndex(ts_meas_pd.index, freq="Y"))['values'].mean()/100
    """#TODO: move to hatyan.timeseries.* or hatyan.kenmerkendewaarden.*. Add minimum # values to calculate monthmean? Make long2array edit simpler with pandas smart stuff?
    numvals_peryearmonth_long = ts_meas_pd.groupby(pd.PeriodIndex(ts_meas_pd.index, freq="M"))['values'].count()
    mean_peryearmonth_array = pd.DataFrame(index=range(1,13))
    for year in mean_peryearmonth_long.index.year.unique():
        bool_year = mean_peryearmonth_long.index.year==year
        mean_peryearmonth_long_oneyr = mean_peryearmonth_long.loc[bool_year]
        mean_peryearmonth_array.loc[mean_peryearmonth_long_oneyr.index.month,year] = mean_peryearmonth_long_oneyr.values
    mean_permonth = mean_peryearmonth_array.mean(axis=1)
    """
    
    if current_station == stat_list[-1]: #indication of last station
        #data_summary
        #print(data_summary[['RDx','RDy']])
        print(data_summary[['data_wl','tstart_wl','tstop_wl','nvals_wl','dupltimes_wl','#nans_wl','#nans_2000to202102a_wl']])
        print(data_summary[['data_ext','dupltimes_ext','#HWgaps_2000to202102_ext']])
        data_summary.to_csv(os.path.join(dir_meas_alldata,'data_summary.csv'))
        
        file_ldb = r'p:\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\computations\model_setup\run_205\20101209-06.ldb' #TODO: make ldb available in code?
        ldb_pd = pd.read_csv(file_ldb, delim_whitespace=True,skiprows=4,names=['RDx','RDy'],na_values=[999.999])
        
        fig_map,ax_map = plt.subplots(figsize=(8,7))
        ax_map.plot(ldb_pd['RDx']/1000,ldb_pd['RDy']/1000,'-k',linewidth=0.4)
        ax_map.plot(cat_locatielijst_ext['RDx']/1000,cat_locatielijst_ext['RDy']/1000,'xk',alpha=0.4)
        ax_map.plot(cat_locatielijst_ext_codeidx.loc[stat_list,'RDx']/1000,cat_locatielijst_ext_codeidx.loc[stat_list,'RDy']/1000,'xr')
        ax_map.plot(data_summary.loc[data_summary['data_ext'],'RDx']/1000,data_summary.loc[data_summary['data_ext'],'RDy']/1000,'xm')#,color='orange')
        """
        for iR, row in cat_locatielijst_ext.iterrows():
            ax_map.text(row['RDx']/1000,row['RDy']/1000,row['Code'])
        """
        ax_map.set_xlim(-50,300)
        ax_map.set_ylim(350,650)
        ax_map.set_title('overview of stations with GETETM2 data')
        ax_map.set_aspect('equal')
        ax_map.grid()
        fig_map.tight_layout()
        fig_map.savefig(os.path.join(dir_meas_alldata,'stations_map.png'))

    #plotting
    file_wl_png = os.path.join(dir_meas_alldata,f'ts_{current_station}.png')
    if os.path.exists(file_wl_png):
        continue #skip the plotting if there is already a png available
    if os.path.exists(file_ext_pkl):
        fig,(ax1,ax2) = hatyan.plot_timeseries(ts=ts_meas_pd, ts_ext=ts_meas_ext_pd)
    else:
        fig,(ax1,ax2) = hatyan.plot_timeseries(ts=ts_meas_pd)
    ax1.plot(mean_peryearmonth_long)
    ax1.plot(mean_peryear_long)
    ax1.set_ylim(-4,4)
    fig_times_ext = [dt.datetime.strptime(x,'%Y%m%d') for x in os.path.basename(dir_meas_alldata).split('_')[2:]]
    ax1.set_xlim(fig_times_ext) #to make clear that e.g. Bath contains too few data
    #ax1.grid()
    fig.savefig(file_wl_png.replace('.png','_alldata.png'))
    ax1.set_xlim(dt.datetime(2000,1,1),dt.datetime(2022,1,1)) #to make clear that e.g. Bath contains too few data
    fig.savefig(file_wl_png)
    plt.close(fig)
    






#### SLOTGEMIDDELDEN
for current_station in []: #stat_list:
    
    ####################
    #READ HWLW
    print(f'slotgemiddelden for {current_station}')
    
    #get validation timeseries
    add_validation = True
    station_name_dict = {'HOEKVHLD':'hoek',#TODO: request data for all stations with DONARcode in filename
                         'HARVT10':'ha10'}
    if current_station in station_name_dict.keys():
        dir_meas_gemHWLWwlAB = r'p:\11208031-010-kenmerkende-waarden-k\work\data_KW-RMM'
        file_yearmeanHW = os.path.join(dir_meas_gemHWLWwlAB,f'{station_name_dict[current_station]}_hw.txt')
        file_yearmeanLW = os.path.join(dir_meas_gemHWLWwlAB,f'{station_name_dict[current_station]}_lw.txt')
        file_yearmeanwl = os.path.join(dir_meas_gemHWLWwlAB,f'{station_name_dict[current_station]}_Z.txt')
        yearmeanHW = pd.read_csv(file_yearmeanHW, delim_whitespace=True, skiprows=1, names=['datetime','values'], parse_dates=['datetime'], na_values=-999.9, index_col='datetime')/100
        yearmeanLW = pd.read_csv(file_yearmeanLW, delim_whitespace=True, skiprows=1, names=['datetime','values'], parse_dates=['datetime'], na_values=-999.9, index_col='datetime')/100
        yearmeanwl = pd.read_csv(file_yearmeanwl, delim_whitespace=True, skiprows=1, names=['datetime','values'], parse_dates=['datetime'], na_values=-999.9, index_col='datetime')/100
    else:
        add_validation = False
    
    #derive yearmean wl from wl values
    file_wl_pkl = os.path.join(dir_meas_alldata,f"{current_station}_measwl.pkl")
    data_pd_meas = pd.read_pickle(file_wl_pkl)
    data_pd_meas = data_pd_meas[['values','QC']] # reduces the memory consumption significantly
    data_pd_meas.index = data_pd_meas.index.tz_localize(None)
    data_pd_meas = data_pd_meas.loc[~(data_pd_meas['QC']==99)]
    wl_mean_peryear_long = data_pd_meas.groupby(pd.PeriodIndex(data_pd_meas.index, freq="y"))['values'].mean()
    wl_mean_peryear_long.index = wl_mean_peryear_long.index.to_timestamp()

    #derive yearmean HWLW from HWLW values
    file_ext_pkl = os.path.join(dir_meas_alldata,f"{current_station}_measext.pkl")
    #if not os.path.exists(file_ext_pkl):
    #    continue
    data_pd_HWLW_all = pd.read_pickle(file_ext_pkl)
    data_pd_HWLW_all.index = data_pd_HWLW_all.index.tz_localize(None)
    if len(data_pd_HWLW_all['HWLWcode'].unique()) > 2:
        data_pd_HWLW_12 = hatyan.calc_HWLW12345to21(data_pd_HWLW_all) #convert 12345 to 12 by taking minimum of 345 as 2 (laagste laagwater) #TODO: this drops first/last value if it is a LW, should be fixed
    else:
        data_pd_HWLW_12 = data_pd_HWLW_all.copy()
    data_pd_HW = data_pd_HWLW_12.loc[data_pd_HWLW_12['HWLWcode']==1]
    data_pd_LW = data_pd_HWLW_12.loc[data_pd_HWLW_12['HWLWcode']==2]
    HW_mean_peryear_long = data_pd_HW.groupby(pd.PeriodIndex(data_pd_HW.index, freq="y"))['values'].mean()
    HW_mean_peryear_long.index = HW_mean_peryear_long.index.to_timestamp()
    LW_mean_peryear_long = data_pd_LW.groupby(pd.PeriodIndex(data_pd_LW.index, freq="y"))['values'].mean()
    LW_mean_peryear_long.index = LW_mean_peryear_long.index.to_timestamp()
    
    #plotting (yearly averages are plotted on 1jan, would be better on 1jul)
    fig,ax1 = plt.subplots(figsize=(14,7))
    if add_validation:
        ax1.plot(yearmeanHW['values'],'+g')
        ax1.plot(yearmeanLW['values'],'+g')
        ax1.plot(yearmeanwl['values'],'+g')
        yearmeanHW_diff = (yearmeanHW['values']-HW_mean_peryear_long).dropna() #TODO: move to data check part, when validationdata for more stations is available
        yearmeanLW_diff = (yearmeanLW['values']-LW_mean_peryear_long).dropna()
        yearmeanwl_diff = (yearmeanwl['values']-wl_mean_peryear_long).dropna()

    ax1.plot(HW_mean_peryear_long,'xr')
    ax1.plot(LW_mean_peryear_long,'xr')
    ax1.plot(wl_mean_peryear_long,'xr')
    ax1.grid()
    ax1.set_ylabel('waterstand [m]')
    ax1.set_title(f'yearly mean HW/wl/LW {current_station}')
    fig.tight_layout()
    fig.savefig(os.path.join(dir_slotgem,f'yearly_values_{current_station}'))
    plt.close()
    
    #TODO: fit with sm.OLS: https://github.com/openearth/sealevel/blob/master/notebooks/analysis/gtsm/nodal-tide.ipynb (does this not work better with monthly values?, PLSS because of potential trendbreuk?, make sure yearly values are on 1jul instead of 1jan)
    







### HAVENGETALLEN
"""
LWaggercode uitleg
TVL;1;1;hoogwater
TVL;1;2;laagwater
TVL;1;3;laagwater 1
TVL;1;4;topagger
TVL;1;5;laagwater 2
"""
culm_addtime = 2*dt.timedelta(hours=24,minutes=50)-dt.timedelta(minutes=20)+dt.timedelta(hours=1) # link with moonculmination (or M2) two days before, 24h rotates entire graph. # furthermore: 2u20min correction, this shifts the x-axis: HW is 2 days after culmination (so 4x25min difference between length of avg moonculm and length of 2 days), 20 minutes (0 to 5 meridian), 1 hour (GMT to MET, alternatively derive moonculm in different tzone)
#data_pd_moonculm = hatyan.astrog_culminations(tFirst=dt.datetime(1979,1,1),tLast=dt.datetime(2022,1,1))
data_pd_moonculm = hatyan.astrog_culminations(tFirst=tstart_dt-culm_addtime,tLast=tstop_dt-culm_addtime)#,tzone='UTC+01:00')
if str(data_pd_moonculm.loc[0,'datetime'].tz) != 'UTC':
    raise Exception(f'culmination data is not in expected timezone (UTC): {data_pd_moonculm.loc[0,"datetime"].tz}')

for current_station in []:#['HARVT10', 'VLISSGN']:#stat_list:

    ####################
    #READ HWLW
    print(f'havengetallen for {current_station}')
    
    file_ext_pkl = os.path.join(dir_meas,f"{current_station}_measext.pkl")
    if not os.path.exists(file_ext_pkl):
        continue
    data_pd_HWLW_all = pd.read_pickle(file_ext_pkl)
    
    #remove timezone-awareness, crop timeseries and apply NAP correction
    data_pd_HWLW_all.index = data_pd_HWLW_all.index.tz_localize(None)
    data_pd_HWLW_all = hatyan.crop_timeseries(data_pd_HWLW_all, times_ext=[tstart_dt,tstop_dt+dt.timedelta(days=30)],onlyfull=False) #TODO: should be possible to crop timeseries to tstop, but results in mising last LWs for HOEKVHLD 2011.0 and possibly others (then for 2021.0 we need also jan2022, but that is not always valid data). Still decide whether to select on extremes or moonculminations in period of interest
    if NAP2005correction:
        data_pd_HWLW_all = nap2005_correction(data_pd_HWLW_all,current_station=current_station)
    
    print('SELECT/CALC HWLW VALUES')
    LWaggercode = 3 # timings LW aardappelgrafiek kloppen voor 1991.0 het best bij LWaggercode=3, misschien doordat eerste laagwater dominant is voor HvH. TODO: delays should then also be used to scale with first LW in gemgetijkromme but now dominant one is used (which depends per station/period, how to automate?). Or simpler: getijkromme1991.0 "Bij meetpunten waar zich aggers voordoen, is, afgezien van de dominantie, de vorm bepaald door de ruwe krommen; dit in tegenstelling tot vroegere bepa-lingen. Bij spring- en doodtij is bovendien de differentiele getijduur, en daarmee de duur rijzing, afgeleid uit de ruwe krommen."
    data_pd_HWLW = data_pd_HWLW_all.loc[(data_pd_HWLW_all['HWLWcode']==1) | (data_pd_HWLW_all['HWLWcode']==2) | (data_pd_HWLW_all['HWLWcode']==LWaggercode)]
    """
    if (LWaggercode == 4) and (4 in data_pd_HWLW_all['HWLWcode'].values): #TODO: this is not used so can be removed?
        # select dominant LW, time of first or dominant LW, LWaggercode must be 4 because of iR
        for iR, HWLWrow in data_pd_HWLW.loc[data_pd_HWLW['HWLWcode'] == 4].iterrows():
            print('%d of %d'%(iR,data_pd_HWLW.index.max()))
            id_min = data_pd_HWLW_all.loc[iR-1:iR+1, 'values'].idxmin()
            data_pd_HWLW.loc[iR,'times'] = data_pd_HWLW_all.loc[id_min,'times'] #on time of lowest LW
            data_pd_HWLW.loc[iR,'values'] = data_pd_HWLW_all.loc[id_min,'values']
            data_pd_HWLW.loc[iR,'HWLWcode'] = LWaggercode
    """
    data_pd_HWLW.index.name = 'times' #index is called 'Tijdstip' if retrieved from DDL.
    data_pd_HWLW = data_pd_HWLW.reset_index() # needed since we need numbered HWLW, HW is a value and LW is value+1
    
    #add duur getijperiode
    data_pd_HWLW['duur getyper hr'] = np.nan
    HW_bool = data_pd_HWLW['HWLWcode']==1
    gtyper_pds_hr = (data_pd_HWLW.loc[HW_bool,'times'].iloc[1:].values - data_pd_HWLW.loc[HW_bool,'times'].iloc[:-1]).dt.total_seconds()/3600
    data_pd_HWLW.loc[HW_bool,'duur getyper hr'] = np.concatenate([gtyper_pds_hr.values,[12+25/60]]) #TODO: concatenate should be in opposite order?
    
    ##### CULMINATIEBEREKENING/HAVENGETALLEN
    print('calculate HW/LWs corresponding to each culmintation') #TODO: maybe make more efficient by working with HWLWno? (might also be easier with gaps)
    data_pd_HWLW['culm'] = np.nan
    data_pd_HWLW['culm_hr'] = np.nan
    data_pd_HWLW['HWLW_delay_hours'] = np.nan
    for iC,culmrow in data_pd_moonculm.iterrows():
        culm_time = culmrow.datetime.tz_localize(None)
        data_pd_HW_selid = (data_pd_HWLW.loc[data_pd_HWLW['HWLWcode']==1,'times']-(culm_time+culm_addtime)).abs().idxmin()
        data_pd_LW_selid = data_pd_HW_selid+1
        if np.abs((data_pd_HWLW.loc[data_pd_HW_selid,'times']-(culm_time+culm_addtime)).total_seconds()/3600) > 8:
            raise Exception('ERROR: no high waters found within 8 hours of culmination at %s +culm_addtime(=%.1f), range HW: \n%s)'%(culm_time, culm_addtime.total_seconds()/3600, data_pd_HWLW.loc[[data_pd_HWLW.index.min(),data_pd_HWLW.index.max()],'times']))
        data_pd_HWLW.loc[[data_pd_HW_selid,data_pd_LW_selid],'culm'] = culm_time
        data_pd_HWLW.loc[[data_pd_HW_selid,data_pd_LW_selid],'culm_hr'] = culm_time.round('h').hour%12

        HW_delay_h = (data_pd_HWLW.loc[data_pd_HW_selid,'times']-(culm_time+culm_addtime))/np.timedelta64(1,'h')
        data_pd_HWLW.loc[data_pd_HW_selid,'HWLW_delay_hours'] = HW_delay_h
        LW_delay_h = (data_pd_HWLW.loc[data_pd_LW_selid,'times']-(culm_time+culm_addtime))/np.timedelta64(1,'h')
        data_pd_HWLW.loc[data_pd_LW_selid,'HWLW_delay_hours'] = LW_delay_h
    
    HW_culmhr_valavg = []
    HW_culmhr_timavg = []
    HW_culmhr_gtyperavg = []
    LW_culmhr_valavg = []
    LW_culmhr_timavg = []
    
    print('calculate medians per hour group for LW and HW (instead of 1991 method: average of subgroups with removal of outliers)')
    for iH in range(0,12): #TODO: havenget ipv culm_hr loop ook pd median op uren doen zoals overschrfreq?
        data_pd_HW_hour = data_pd_HWLW[(data_pd_HWLW['HWLWcode']==1) & (data_pd_HWLW['culm_hr']==iH)]
        HW_culmhr_valavg.append(np.median(data_pd_HW_hour['values']))
        HW_culmhr_timavg.append(np.median(data_pd_HW_hour['HWLW_delay_hours']))
        HW_culmhr_gtyperavg.append(np.mean(data_pd_HW_hour['duur getyper hr']))
        data_pd_LW_hour = data_pd_HWLW[(data_pd_HWLW['HWLWcode']!=1) & (data_pd_HWLW['culm_hr']==iH)]
        LW_culmhr_valavg.append(np.median(data_pd_LW_hour['values']))
        LW_culmhr_timavg.append(np.median(data_pd_LW_hour['HWLW_delay_hours']))
    HW_culmhr_valavg_gemtij = np.mean(HW_culmhr_valavg)
    HW_culmhr_timavg_gemtij = np.mean(HW_culmhr_timavg)
    HW_culmhr_gtyperavg_gemtij = np.mean(HW_culmhr_gtyperavg)
    LW_culmhr_valavg_gemtij = np.mean(LW_culmhr_valavg)
    LW_culmhr_timavg_gemtij = np.mean(LW_culmhr_timavg)

    print('HWLW FIGUREN PER TIJDSKLASSE, INCLUSIEF MEDIAN LINE')
    HW_data = data_pd_HWLW[(data_pd_HWLW['HWLWcode']==1)]
    LW_data = data_pd_HWLW[(data_pd_HWLW['HWLWcode']!=1)] #HWLWcode==2 or HWLWcode==LWaggercode (=3)
    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(18,8), sharex=True)
    ax1.set_title('HW values %s'%(current_station))
    ax1.plot(HW_data['culm_hr'],HW_data['values'],'.')
    ax1.plot(range(0,12),HW_culmhr_valavg,'.-')
    ax2.set_title('LW values %s'%(current_station))
    ax2.plot(LW_data['culm_hr'],LW_data['values'],'.')
    ax2.plot(range(0,12),LW_culmhr_valavg,'.-')
    ax3.set_title('HW time delays %s'%(current_station))
    ax3.plot(HW_data['culm_hr'],HW_data['HWLW_delay_hours'],'.')
    ax3.plot(range(0,12),HW_culmhr_timavg,'.-')
    ax4.set_title('LW time delays %s'%(current_station))
    ax4.plot(LW_data['culm_hr'],LW_data['HWLW_delay_hours'],'.')
    ax4.plot(range(0,12),LW_culmhr_timavg,'.-')
    ax4.set_xlim([0-0.5,12-0.5])
    fig.tight_layout()
    fig.savefig(os.path.join(dir_havget,f'HWLW_pertijdsklasse_inclmedianline_{current_station}'))


    file_outname = os.path.join(dir_havget, 'aardappelgrafiek_%s_%s_aggercode%s'%(year_slotgem, current_station, LWaggercode))

    with open('%s.txt'%(file_outname),'w') as f:
        f.write('### HIGH WATERS ###\n')
        #for iHR, time, val in zip(range(len(HW_culmhr_valavg)), HW_culmhr_timavg, HW_culmhr_valavg):
        #    f.write('%2i: %6s %6.2f\n'%(iHR, str(dt.timedelta(hours=time)), val))
        f.write('%2i: %-14s %6.2f  gtyper %6s\n'%(0,str(dt.timedelta(hours=HW_culmhr_timavg[0])), HW_culmhr_valavg[0], str(dt.timedelta(hours=HW_culmhr_gtyperavg[0]))))
        f.write('av: %-14s %6.2f  gtyper %6s\n'%(str(dt.timedelta(hours=HW_culmhr_timavg_gemtij)), HW_culmhr_valavg_gemtij, str(dt.timedelta(hours=HW_culmhr_gtyperavg_gemtij))))
        f.write('%2i: %-14s %6.2f  gtyper %6s\n'%(6,str(dt.timedelta(hours=HW_culmhr_timavg[6])), HW_culmhr_valavg[6], str(dt.timedelta(hours=HW_culmhr_gtyperavg[6]))))
        f.write('### LOW WATERS ###\n')
        #for iHR, time, val in zip(range(len(LW_culmhr_valavg)), LW_culmhr_timavg, LW_culmhr_valavg):
        #    f.write('%2i: %6s %6.2f\n'%(iHR, str(dt.timedelta(hours=time)), val))
        f.write('%2i: %-14s %6.2f  daling %6s\n'%(0,str(dt.timedelta(hours=LW_culmhr_timavg[0])), LW_culmhr_valavg[0], str(dt.timedelta(hours=LW_culmhr_timavg[0]-HW_culmhr_timavg[0]))))
        f.write('av: %-14s %6.2f  daling %6s\n'%(str(dt.timedelta(hours=LW_culmhr_timavg_gemtij)), LW_culmhr_valavg_gemtij, str(dt.timedelta(hours=LW_culmhr_timavg_gemtij-HW_culmhr_timavg_gemtij))))
        f.write('%2i: %-14s %6.2f  daling %6s\n'%(6,str(dt.timedelta(hours=LW_culmhr_timavg[6])), LW_culmhr_valavg[6], str(dt.timedelta(hours=LW_culmhr_timavg[6]-HW_culmhr_timavg[6]))))

    print('AARDAPPELGRAFIEK')
    def timeTicks(x, pos):
        d = dt.timedelta(hours=np.abs(x))
        if np.sign(x)>0:
            d_str = str(d)
        else:
            d_str = '-'+str(d)
        return d_str

    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(7.5,4), sharex=False)
    ax1.set_title(f'HW {current_station} {year_slotgem}')
    ax1.set_xlabel('maansverloop in uu:mm:ss' )
    ax1.set_ylabel('waterstand in m t.o.v. NAP')
    ax1.plot(HW_culmhr_timavg,HW_culmhr_valavg,'.-',label=current_station)
    for iHR in range(0,12):
        ax1.text(HW_culmhr_timavg[iHR],HW_culmhr_valavg[iHR], str(iHR))
    ax1.xaxis.set_major_formatter(timeTicks)
    ax1.grid()
    ax2.set_title(f'LW {current_station} {year_slotgem}')
    ax2.set_xlabel('maansverloop in uu:mm:ss' )
    ax2.set_ylabel('waterstand in m t.o.v. NAP')
    ax2.plot(LW_culmhr_timavg,LW_culmhr_valavg,'.-',label=current_station)
    for iHR in range(0,12):
        ax2.text(LW_culmhr_timavg[iHR],LW_culmhr_valavg[iHR], str(iHR))
    ax2.xaxis.set_major_formatter(timeTicks)
    ax2.grid()
    #set equal ylims
    ax1_xlimold = ax1.get_xlim()
    ax2_xlimold = ax2.get_xlim()
    ax1_ylimold = ax1.get_ylim()
    ax2_ylimold = ax2.get_ylim()
    xlimrange = 2
    ylimrange = 1
    ax1.set_xlim([np.mean(ax1_xlimold)-xlimrange/2,np.mean(ax1_xlimold)+xlimrange/2])
    ax2.set_xlim([np.mean(ax2_xlimold)-xlimrange/2,np.mean(ax2_xlimold)+xlimrange/2])
    ax1.set_ylim([np.mean(ax1_ylimold)-ylimrange/2,np.mean(ax1_ylimold)+ylimrange/2])
    ax2.set_ylim([np.mean(ax2_ylimold)-ylimrange/2,np.mean(ax2_ylimold)+ylimrange/2])
    #plot gemtij
    ax1.plot(ax1.get_xlim(),[HW_culmhr_valavg_gemtij,HW_culmhr_valavg_gemtij],'k--')
    ax1.plot([HW_culmhr_timavg_gemtij,HW_culmhr_timavg_gemtij],ax1.get_ylim(),'k--')
    ax2.plot(ax2.get_xlim(),[LW_culmhr_valavg_gemtij,LW_culmhr_valavg_gemtij],'k--')
    ax2.plot([LW_culmhr_timavg_gemtij,LW_culmhr_timavg_gemtij],ax2.get_ylim(),'k--')
    fig.tight_layout()
    fig.savefig(file_outname)







##### gemiddelde getijkrommen
# slotgemiddelden uit:
# =============================================================================
# slotGem  = 'rapportRWS'
# slotGem  = 'havengetallen2011'
slotGem  = 'havengetallen2011_PLSS'

fig_sum,ax_sum = plt.subplots(figsize=(14,7))
for current_station in []:#['HOEKVHLD','HARVT10']:#stat_list:
    """
    uit: gemiddelde getijkrommen 1991.0
        
    Voor meetpunten in het onbeinvloed gebied is per getijfase eerst een "ruwe kromme" berekend met de resultaten van de harmonische analyse, 
    welke daarna een weinig is bijgesteld aan de hand van de volgende slotgemiddelden:
    gemiddeld hoog- en laagwater, duur daling. Deze bijstelling bestaat uit een eenvoudige vermenigvuldiging.
    
    Voor de ruwe krommen voor springtij en doodtij is het getij voorspeld voor een jaar met gemiddelde helling maansbaan 
    met uitsluitend zuivere combinaties van de componenten M2 en S2:
    tabel: Gebruikte componenten voor de spring- en doodtijkromme
    SM, 3MS2, u2, M2, S2, 2SM2, 3MS4, M4, MS4, 
    4MS6, M6, 2MS6, M8, 3MS8, M10, 4MS10, M12, 5MS12
       
    In het aldus gemodelleerde getij is de vorm van iedere getijslag, gegeven de getijfase, identiek. 
    Vervolgens is aan de hand van de havengetallen een springtij- en een doodtijkromme geselecteerd.
    
    Voor de ruwe krommen voor gemiddeld tij zijn uitsluitend zuivere harmonischen van M2 gebruikt: M2, M4, M6, M8, M10, M12, 
    waarbij de amplituden per component zijn vervangen door de wortel uit de kwadraatsom van de amplituden 
    van alle componenten in de betreffende band, voor zover voorkomend in de standaardset van 94 componenten. 
    Zoals te verwachten is de verhouding per component tussen deze wortel en de oorspronkelijke amplitude voor alle plaatsen gelijk.
    tabel: Verhouding tussen amplitude en oorspronkelijke amplitude
    M2 (tweemaaldaagse band) 1,06
    M4 1,28
    M6 1,65
    M8 2,18
    M10 2,86
    M12 3,46
    
    In het aldus gemodelleerde getij is de vorm van iedere getijslag identiek, met een getijduur van 12 h 25 min.
    Bij meetpunten waar zich aggers voordoen, is, afgezien van de dominantie, de vorm bepaald door de ruwe krommen; 
    dit in tegenstelling tot vroegere bepalingen. Bij spring- en doodtij is bovendien de differen-tiele getijduur, 
    en daarmee de duur rijzing, afgeleid uit de ruwe krommen.
    
    De aanpassing aan de slotgemiddelden kwam bij springtij en gemiddeld tij steeds op zeer geringe wijzigingen neer, 
    bij doodtij echter werd het tijverschil door de ruwe krommen ca. 10% onderschat. 
    Het verschil tussen spring- en gemiddeld tij enerzijds en doodtij anderzijds is vermoedelijk 
    toe te schrijven aan enkele algemene eigenschappen van harmonische componenten.
    
    Voor wat betreft de juistheid van de vorm van de krommen kan slechts worden afgegaan op de ervaring, 
    aangezien geen der vroegere bepalingen consis-tente en reproduceerbare resultaten opleverde.
    
    De boven omschreven methode voor gemiddeld tij biedt de mogelijkheid op geheel analoge wijze krommen te berekenen, 
    zodra maar enkele componenten (eventueel alleen M2, M4 en M6) bekend zijn. 
    Te denken valt aan: buitenlandse meetpunten; gedurende korte tijd bemeten locaties; 
    modelresultaten; hypothetische hydrologische omstandigheden.
    """

    print(f'gem getijkrommen for {current_station}')

    dir_vali_krommen = r'p:\archivedprojects\11205258-005-kpp2020_rmm-g5\C_Work\00_KenmerkendeWaarden\07_Figuren\figures_ppSCL_2\final20201211'
    file_vali_doodtijkromme = os.path.join(dir_vali_krommen,f'doodtijkromme_{current_station}_{slotGem}.csv')
    file_vali_gemtijkromme = os.path.join(dir_vali_krommen,f'gemGetijkromme_{current_station}_{slotGem}.csv')
    file_vali_springtijkromme = os.path.join(dir_vali_krommen,f'springtijkromme_{current_station}_{slotGem}.csv')
    
    if year_slotgem not in [2011,'2011_olddata']:
        raise Exception(f'gemiddelde getijkromme only possible for 2011: {year_slotgem}')
        
    #HvH
    if current_station == 'HOEKVHLD':
        #TODO: make this automatic, also apply to retrieved data at top of script if necessary, or remove entirely. Also ask stendert where these values come from (delay values seem a bit weird)
        if slotGem == 'rapportRWS': #2011.0 rapport van Douwe Dillingh
            HW_sp = 1.32
            HW_av = 1.15
            HW_np = 0.90
            LW_sp = -0.63
            LW_av = -0.60
            LW_np = -0.55
            # tijdsduur voor daling
            tD_sp = dt.timedelta(hours=6-1,minutes=51-32)
            tD_av = dt.timedelta(hours=7-1,minutes=17-34)
            tD_np = dt.timedelta(hours=7-1,minutes=39-36)
        elif slotGem == 'havengetallen2011': #KW-RMM havengetallen programma
            HW_sp = 1.28
            HW_av = 1.11
            HW_np = 0.87
            LW_sp = -0.64
            LW_av = -0.61
            LW_np = -0.57
            # tijdsduur voor daling
            tD_sp = dt.timedelta(hours=6-1,minutes=51-32,seconds=20.5-58)
            tD_av = dt.timedelta(hours=7-1,minutes=21-33,seconds=29-13.875)
            tD_np = dt.timedelta(hours=7-1,minutes=45-34,seconds=19-38)
        elif slotGem == 'havengetallen2011_PLSS': #KW-RMM havengetallen programma, bewerkt met PLSS correctie van Douwe (in excelsheet)
            HW_sp = 1.32
            HW_av = 1.15
            HW_np = 0.91
            LW_sp = -0.63
            LW_av = -0.60
            LW_np = -0.56
            # tijdsduur voor daling
            tD_sp = dt.timedelta(hours=6-1,minutes=51-32,seconds=20.5-58)
            tD_av = dt.timedelta(hours=7-1,minutes=21-33,seconds=29-13.875)
            tD_np = dt.timedelta(hours=7-1,minutes=45-34,seconds=19-38)

        # tijdsverschil voor verplaatsing HvH-->Maasmond
        tDiff_sp = dt.timedelta(minutes=-5)
        tDiff_av = dt.timedelta(minutes=-5)
        tDiff_np = dt.timedelta(minutes=-5)

    #HV10
    elif current_station == 'HARVT10':

        if slotGem == 'rapportRWS':
            HW_sp = 1.45
            HW_av = 1.24
            HW_np = 0.93
            LW_sp = -0.92
            LW_av = -0.86
            LW_np = -0.77
            # tijdsduur voor daling
            tD_sp = dt.timedelta(hours=6-1,minutes=44-14)
            tD_av = dt.timedelta(hours=7-1,minutes=5-14)
            tD_np = dt.timedelta(hours=7-1,minutes=32-13)
            # tijdsverschil met HvH voor start van tijdserie van kromme (KW kustwater en grote rivieren, G, verschil havengetal LW)
            tDiff_sp = dt.timedelta(minutes=14-32)
            tDiff_av = dt.timedelta(minutes=14-34)
            tDiff_np = dt.timedelta(minutes=13-36)

        elif slotGem == 'havengetallen2011':
            HW_sp = 1.41
            HW_av = 1.21
            HW_np = 0.92
            LW_sp = -0.94
            LW_av = -0.86
            LW_np = -0.77
            # tijdsduur voor daling
            tD_sp = dt.timedelta(hours=6-1,minutes=45-14,seconds=53-45.5)
            tD_av = dt.timedelta(hours=7-1,minutes=5-13,seconds=41.625-33.416667)
            tD_np = dt.timedelta(hours=7-1,minutes=31-11,seconds=17-12)
            # tijdsverschil met HvH voor start van tijdserie van kromme (verschil in havengetal HW)
            tDiff_sp = dt.timedelta(minutes=14-32,seconds=45.5-58)
            tDiff_av = dt.timedelta(minutes=13-33,seconds=33.416667-13.875)
            tDiff_np = dt.timedelta(minutes=11-34,seconds=12-38)
        elif slotGem == 'havengetallen2011_PLSS':
            HW_sp = 1.44
            HW_av = 1.24
            HW_np = 0.95
            LW_sp = -0.94
            LW_av = -0.86
            LW_np = -0.77
            # tijdsduur voor daling
            tD_sp = dt.timedelta(hours=6-1,minutes=45-14,seconds=53-45.5)
            tD_av = dt.timedelta(hours=7-1,minutes=5-13,seconds=41.625-33.416667)
            tD_np = dt.timedelta(hours=7-1,minutes=31-11,seconds=17-12)
            # tijdsverschil met HvH voor start van tijdserie van kromme (verschil in havengetal HW)
            tDiff_sp = dt.timedelta(minutes=14-32,seconds=45.5-58)
            tDiff_av = dt.timedelta(minutes=13-33,seconds=33.416667-13.875)
            tDiff_np = dt.timedelta(minutes=11-34,seconds=12-38)

    else:
        raise Exception(f'station {current_station} not implemented for gemiddelde getijkrommen')
        
    #load measurement data
    file_wl_pkl = os.path.join(dir_meas,f"{current_station}_measwl.pkl")
    ts_meas_pd = pd.read_pickle(file_wl_pkl)
    ts_meas_pd = ts_meas_pd[['values','QC']] # reduces the memory consumption significantly
    ts_meas_pd.index = ts_meas_pd.index.tz_localize(None)
    ts_meas_pd = ts_meas_pd.loc[~(ts_meas_pd['QC']==99)]
    ts_meas_pd = hatyan.crop_timeseries(ts_meas_pd, times_ext=[tstart_dt,tstop_dt-dt.timedelta(minutes=10)],onlyfull=False) #TODO: moet onlyfull wel False zijn?
    if NAP2005correction:
        ts_meas_pd = nap2005_correction(ts_meas_pd,current_station)
    
    # =============================================================================
    # Hatyan voor 10 jaar (alle componenten voor gemiddelde getijcyclus)
    # =============================================================================
    
    const_list = hatyan.get_const_list_hatyan('year') #this should not be changed, since higher harmonics are necessary
    hatyan_settings = hatyan.HatyanSettings(nodalfactors = True,
                                            fu_alltimes = False, # False is RWS-default
                                            xfac = True, #wordt niet besproken, moet die wel aan?
                                            analysis_peryear = True,
                                            xTxmat_condition_max=15, #TODO: for some reason this is necessary for HOEKVHLD 2006 (default=10), also strong difference in springneap ts when using smaller component set, what is happening?
                                            return_allyears=True)
    comp_frommeasurements_avg, comp_frommeasurements_allyears = hatyan.get_components_from_ts(ts_meas_pd, const_list=const_list, hatyan_settings=hatyan_settings)
    
    # =============================================================================
    # gemiddelde getijkromme
    # =============================================================================
    #kwadraatsommen voor M2 tot M12
    components_av = ['M2','M4','M6','M8','M10','M12']
    comp_av = comp_frommeasurements_avg.loc[components_av]
    for comp_higherharmonics in components_av:
        iM = int(comp_higherharmonics[1:])
        bool_endswithiM = comp_frommeasurements_avg.index.str.endswith(str(iM)) & comp_frommeasurements_avg.index.str.replace(str(iM),'').str[-1].str.isalpha()
        comp_iM = comp_frommeasurements_avg.loc[bool_endswithiM]
        comp_av.loc[comp_higherharmonics,'A'] = np.sqrt((comp_iM['A']**2).sum()) #kwadraatsom
    
    print('verhouding tussen originele en kwadratensom componenten:')
    print(comp_av/comp_frommeasurements_avg.loc[components_av]) #TODO: values are different than 1991.0 document, but could be because of different year so check with 1981-1991 data
    
    comp_av.loc['A0'] = comp_frommeasurements_avg.loc['A0']
    times_pred_1mnth = pd.date_range(start=dt.datetime(tstop_dt.year, 1, 1, 0, 0), end=dt.datetime(tstop_dt.year, 2, 1, 0, 0), freq='60 S') # TODO hatyan: when using <60sec, hatyan.calc_HWLW() goes wrong, probably since there is a minute-rounding thing somewhere, fix this
    prediction_av = hatyan.prediction(comp_av, times_pred_all=times_pred_1mnth, hatyan_settings=hatyan_settings)
    prediction_av_ext = hatyan.calc_HWLW(ts=prediction_av)#,calc_HWLWlocal=False)


    # karateristieken uit ruwe gemiddelde getijkromme >> schalingsratio
    idHW_av = prediction_av_ext.index[prediction_av_ext.HWLWcode==1][:-1]
    # selecteer eerste laagwater
    idLW_av = prediction_av_ext.iloc[np.where(prediction_av_ext.HWLWcode==1)[0][:-1]+1].index
    
    def get_tide_meanext_valstimes(ts_ext):
        #TODO: vorm van iedere getijslag is in principe identiek (maar niet als deze is afgerond op 1min), dus onderstaande is eigenlijk niet nodig hoewel er nu het risico is op 12:24 of 12:26 getijduur >> 10sec voorspelling maken (maar kan nu niet). Afronding prediction_av op 1min zorgt voor timeUp/timeDown die meestal 0 maar soms 60 seconden van elkaar verschillen. Geldt ook voor spring en doodtij?
        #TODO: ongetwijfeld gaat er iets in dit script uit van 1/2/1/2 alternerende HWLW, bouw hier een check voor in (eg identify potential gaps)
        HW_val_mean = ts_ext.loc[ts_ext['HWLWcode']==1,'values'].mean() # np.mean(HW)
        LW_val_mean = ts_ext.loc[ts_ext['HWLWcode']==2,'values'].mean() # np.mean(LW)
        timediff = pd.Series(ts_ext.index,index=ts_ext.index).diff()
        time_up = timediff.loc[ts_ext['HWLWcode']==1].mean() # np.mean(timeUp)
        time_down = timediff.loc[ts_ext['HWLWcode']==2].mean() # np.mean(timeDown)
        return HW_val_mean, LW_val_mean, time_up, time_down
    HW_cav, LW_cav, tU_cav, tD_cav = get_tide_meanext_valstimes(prediction_av_ext)
    
    # tijd daling uit metingen
    #tD_av = tD_cav
    tU_av = M2_period_timedelta - tD_av #dt.timedelta(hours=12,minutes=25)-tD_av
    
    # bereken schalingsratio's voor kromme #TODO: _av etc komen nu uit hardgecodeerd start van script, maar moeten komen uit live afgeleidde havengetallen/slotgemiddelden
    rHW_av = HW_av/HW_cav
    rLW_av = LW_av/LW_cav
    rtU_av = tU_av/tU_cav
    rtD_av = tD_av/tD_cav

    
    def vermenigvuldig_kromme(ts, timesHW, timesLW, ratioHW, ratioLW, ratioDown, ratioUp):
        #TODO: is boven/onder nul goede indicator aangezien A0 ook wordt gebruikt? >> misschien boven/onder A0 of A0 weglaten?
        ts_corr = ts.copy()
        ts_corr['values'][ts_corr['values']>0] = ts_corr['values'][ts_corr['values']>0]*ratioHW
        ts_corr['values'][ts_corr['values']<0] = ts_corr['values'][ts_corr['values']<0]*ratioLW
        ts_corr['times'] = ts_corr.index #this is necessary since datetimeindex with freq is not editable, and Series is editable
        for i in np.arange(0,len(timesHW)): #TODO: what happens if ts starts with LW instead of HW? TODO: prevent timesHW/timesLW possible?
            #iHW = ts_corr.index.get_loc(timesHW[i])
            #iLW = ts_corr.index.get_loc(timesLW[i])
            tide_HWtoLW = ts_corr.loc[timesHW[i]:timesLW[i]]
            ts_corr.loc[timesHW[i]:timesLW[i],'times'] = pd.date_range(start=ts_corr.loc[timesHW[i],'times'],freq=f'{int(ratioDown*1e9*60)} N',periods=len(tide_HWtoLW))
            if i == len(timesHW)-1: #not for last HW
                continue
            tide_LWtoHW = ts_corr.loc[timesLW[i]:timesHW[i+1]]
            ts_corr.loc[timesLW[i]:timesHW[i+1],'times'] = pd.date_range(start=ts_corr.loc[timesLW[i],'times'],freq=f'{int(ratioUp*1e9*60)} N',periods=len(tide_LWtoHW))
        ts_corr = ts_corr.set_index('times',drop=True)
        return ts_corr
    
    #vermenigvuldiging van kromme met ratio's
    prediction_av_corr = vermenigvuldig_kromme(prediction_av, idHW_av, idLW_av, rHW_av, rLW_av, rtD_av, rtU_av)
    
    fig,(ax1,ax2) = hatyan.plot_timeseries(ts=prediction_av,ts_ext=prediction_av_ext)
    ax1.plot(prediction_av_corr['values'],'r',label='gecorrigeerde kromme')
    ax1.legend(labels=['ruwe kromme','0m+NAP','gemiddelde waterstand','hoogwater','laagwater','gecorrigeerde kromme'],loc=4)
    ax1.set_ylabel('waterstand [m]')
    ax1.set_title('gemiddelde getijkromme')
    fig.savefig(os.path.join(dir_gemgetij,"gemGetijkromme_%s_%s.png"%(current_station,slotGem)))
    
    def shift_HW_tostart(ts, timesHW, tstop_dt, tDiff):
        #TODO: this timeshift derived from old csv writing should be eliminated, maybe not necessary to correct HvH/HARVT10 if normal times are accounted for
        bool_av = ts.index>=timesHW[0]
        ts_shift = ts.loc[bool_av]
        ts_shift.index -= timesHW[0]-tstop_dt-tDiff
        return ts_shift
    
    prediction_av_corr_timeshift = shift_HW_tostart(prediction_av_corr, idHW_av, tstop_dt, tDiff_av)
    
    prediction_av_corr_timeshift.to_csv(os.path.join(dir_gemgetij,"gemGetijkromme_%s_%s.csv"%(current_station,slotGem)),float_format='%.3f',date_format='%Y-%m-%d %H:%M:%S')
    fig,(ax1,ax2) = hatyan.plot_timeseries(ts=prediction_av_corr_timeshift,ts_ext=prediction_av_ext)
    if file_vali_gemtijkromme is not None:
        data_vali_gemtij = pd.read_csv(file_vali_gemtijkromme,index_col=0,parse_dates=True)
        ax1.plot(data_vali_gemtij)
    fig.savefig(os.path.join(dir_gemgetij,"springdoodtijkromme_gemiddeld_%s_%s.png"%(current_station,slotGem)))
    ax1.set_xlim(tstop_dt-dt.timedelta(days=0.5),tstop_dt+dt.timedelta(days=4))
    fig.savefig(os.path.join(dir_gemgetij,"springdoodtijkromme_gemiddeld_%s_%s_zoom.png"%(current_station,slotGem)))
    
    
    
    
    # =============================================================================
    # Hatyan voor 1 jaar met gemiddelde helling maansbaan (voor afleiden spring-doodtijcyclus)
    #DONE: 2001 heeft gemiddelde nodalfactor f voor M2 (getijkrommen 1991.0 spreekt van "Voor de ruwe krommen voor springtij en doodtij is het getij voorspeld voor een jaar met gemiddelde helling maansbaan")
    #DONE: de analyse wordt op die metingen gedaan, maar de predictie vervolgens op een andere periode (times_ext) Moet het niet times_ext_pred zijn? (script is dan veel trager omdat het een jaar ipv een maand is, vooral door write_csv, resultaten zijn vrijwel gelijk)
    # =============================================================================
    components_sn = ['SM','3MS2','MU2','M2','S2','2SM2','3MS4','M4','MS4','4MS6','M6','2MS6','M8','3MS8','M10','4MS10','M12','5MS12','A0'] #TODO: should A0 be added since we look at zerocrossings eventually? >> makes a difference with A0 far from 0 since tD/tU scaling is then different?
    
    # derive f values for M2 and select year where the value is closest to 0. TODO: it seems to not matter too much what year is chosen, but maybe for the scaling factors?
    yearcenters_time = pd.date_range(start=tstart_dt, end=tstop_dt, freq='Y') - dt.timedelta(days=364/2)
    yearcenters_ffactor = hatyan.get_schureman_f(const_list=['M2'], dood_date=yearcenters_time, xfac=hatyan_settings.xfac)
    yearcenters_ffactor.columns = yearcenters_time
    year_neutralffactor = (yearcenters_ffactor.T['M2']-1).abs().idxmin().year
    
    ts_measurements_oneyear = hatyan.crop_timeseries(ts_meas_pd, times_ext=[dt.datetime(year_neutralffactor,1,1),dt.datetime(year_neutralffactor,12,31,23,50,00)])
    comp_oneyear_sncomp, dummy = hatyan.get_components_from_ts(ts_measurements_oneyear, const_list=components_sn, hatyan_settings=hatyan_settings)
    comp_oneyear_sncomp = comp_oneyear_sncomp.loc[components_sn]
    
    #TODO: use this to automatically select year with neutrale f voor M2 (helling maansbaan) >> slighly different values, mainly for SM phase, due to componentset difference (const_list vs components_sn)
    #comp_oneyear_minffactor = comp_frommeasurements_allyears.loc[components_sn,(slice(None),year_minffactor)]
    #comp_oneyear_minffactor.columns = comp_oneyear_minffactor.columns.droplevel(1)
    #print(comp_oneyear_sncomp-comp_oneyear_minffactor)
    
    #hatyan_settings_sn = hatyan.HatyanSettings(nodalfactors = False) #TODO: year does not matter too much (maybe it does for scaling), but nodalfactors=False does matter a bit for doodtij duration
    prediction_sn     = hatyan.prediction(comp_oneyear_sncomp, times_pred_all=times_pred_1mnth, hatyan_settings=hatyan_settings)
    
    #TODO KW-RMM2020: "In het geval van aggers is het eerste laagwater gebruikt." >> laagste laagwater wordt genomen
    prediction_sn_ext = hatyan.calc_HWLW(ts=prediction_sn)#, calc_HWLW345=True)
    
    #karateristieken uit ruwe spring-/doodtijkromme >> schalingsratio
    #selecteer alle hoogwaters en opvolgende laagwaters
    idHW_sn = prediction_sn_ext.index[prediction_sn_ext.HWLWcode==1][:-1]
    idLW_sn = prediction_sn_ext.iloc[np.where(prediction_sn_ext.HWLWcode==1)[0][:-1]+1].index
    
    zero_crossings_bool = np.sign(prediction_sn['values']).diff()>0
    zero_crossings_times = prediction_sn.loc[zero_crossings_bool].index #list of zero crossings timestamps

    # maak ruwe springtijkromme (selecteer getijslag na maximale HW)
    time_highestHW = prediction_sn_ext['values'][idHW_sn].idxmax()
    zerocrossing_highestHW_idx = (np.abs(zero_crossings_times-time_highestHW)).argmin()
    is1 = zero_crossings_times[zerocrossing_highestHW_idx]
    is2 = zero_crossings_times[zerocrossing_highestHW_idx+1]
    tC_sp = is2-is1
    tU_sp = tC_sp - tD_sp
    # repeat ruwe springtijkromme in time. #TODO: this shifts the getijkromme in time, which should probably not happen (also for doodtij)
    prediction_sp_one = prediction_sn.loc[is1:is2].iloc[:-3] #TODO: -3 is nodig om reproductie oude lijnen te krijgen, maar dat is niet goed (moet -1 zijn) en je ziet ook een hickup bij ieder begin/eind (also for doodtij)
    prediction_sp = pd.DataFrame(index=prediction_sn.index)
    prediction_sp['values'] = np.tile(prediction_sp_one['values'].values,int(np.ceil(len(prediction_sn.index)/len(prediction_sp_one))))[0:len(prediction_sn.index)]
    prediction_sp_ext = hatyan.calc_HWLW(ts=prediction_sp)
    
    # maak ruwe doodtijkromme (selecteer getijslag na minimale HW)
    time_lowestHW = prediction_sn_ext['values'][idHW_sn].idxmin()
    zerocrossing_lowestHW_idx = (np.abs(zero_crossings_times-time_lowestHW)).argmin()
    in1 = zero_crossings_times[zerocrossing_lowestHW_idx]
    in2 = zero_crossings_times[zerocrossing_lowestHW_idx+1]
    tC_np = in2-in1
    tU_np = tC_np - tD_np
    # repeat ruwe doodtijkromme in time
    prediction_np_one = prediction_sn.loc[in1:in2].iloc[:-3]
    prediction_np = pd.DataFrame(index=prediction_sn.index)
    prediction_np['values'] = np.tile(prediction_np_one['values'].values,int(np.ceil(len(prediction_sn.index)/len(prediction_np_one))))[0:len(prediction_sn.index)] 
    prediction_np_ext = hatyan.calc_HWLW(ts=prediction_np)
        
    # plot selection of neap/spring
    fig, (ax1,ax2) = hatyan.plot_timeseries(ts=prediction_sn,ts_ext=prediction_sn_ext)
    ax1.plot(prediction_sp_one['values'],'r')
    ax1.plot(prediction_np_one['values'],'r')
    ax1.legend(labels=['ruwe kromme','0m+NAP','gemiddelde waterstand','hoogwater','laagwater','kromme spring','kromme neap'],loc=4)
    ax1.set_ylabel('waterstand [m]')
    ax1.set_title('spring- en doodtijkromme')
    fig.savefig(os.path.join(dir_gemgetij,"springdoodtijkromme_%s_%s.png"%(current_station,slotGem)))
        
    print('SPRINGTIJ')
    #karakteristieken springtij
    #selecteer alle hoogwaters en opvolgende laagwaters
    idHW_sp = prediction_sp_ext.index[prediction_sp_ext.HWLWcode==1][:-1]
    idLW_sp = prediction_sp_ext.iloc[np.where(prediction_sp_ext.HWLWcode==1)[0][:-1]+1].index
    
    HW_csp, LW_csp, tU_csp, tD_csp = get_tide_meanext_valstimes(prediction_sp_ext)
    tU_csp = tC_sp-tD_csp
    rHW_sp = HW_sp/HW_csp
    rLW_sp = LW_sp/LW_csp
    rtU_sp = tU_sp/tU_csp
    rtD_sp = tD_sp/tD_csp
    ratio_tideduration_sp = M2_period_timedelta/(tU_csp+tD_csp) #TODO: use this to convert to 12u25m timeseries? beware on minute/second rounding
    
    #vermenigvuldiging van kromme met ratios
    print('vermenigvuldig_kromme')
    prediction_sp_corr = vermenigvuldig_kromme(prediction_sp, idHW_sp, idLW_sp, rHW_sp, rLW_sp, rtD_sp, rtU_sp)
    print('timeshift')
    prediction_sp_corr_timeshift = shift_HW_tostart(prediction_sp_corr, idHW_sp, tstop_dt, tDiff_sp)
    
    print('write to csv')
    prediction_sp_corr_timeshift.to_csv(os.path.join(dir_gemgetij,"springtijkromme_%s_%s.csv"%(current_station,slotGem)),float_format='%.3f',date_format='%Y-%m-%d %H:%M:%S')
    print('plot figure')
    fig,(ax1,ax2) = hatyan.plot_timeseries(ts=prediction_sp_corr_timeshift,ts_ext=prediction_sp_ext)
    if file_vali_springtijkromme is not None:
        data_vali_springtij = pd.read_csv(file_vali_springtijkromme,index_col=0,parse_dates=True)
        ax1.plot(data_vali_springtij)
    print('save figure')
    fig.savefig(os.path.join(dir_gemgetij,"springdoodtijkromme_spring_%s_%s.png"%(current_station,slotGem)))
    ax1.set_xlim(tstop_dt-dt.timedelta(days=0.5),tstop_dt+dt.timedelta(days=4))
    fig.savefig(os.path.join(dir_gemgetij,"springdoodtijkromme_spring_%s_%s_zoom.png"%(current_station,slotGem)))
    print('SPRINGTIJ finished')


    print('DOODTIJ')
    #karakteristieken doodtij
    #selecteer alle hoogwaters en opvolgende laagwaters
    idHW_np = prediction_np_ext.index[prediction_np_ext.HWLWcode==1][:-1]
    idLW_np = prediction_np_ext.iloc[np.where(prediction_np_ext.HWLWcode==1)[0][:-1]+1].index

    HW_cnp, LW_cnp, tU_cnp, tD_cnp = get_tide_meanext_valstimes(prediction_np_ext)
    tU_cnp = tC_np-tD_cnp
    rHW_np = HW_np/HW_cnp
    rLW_np = LW_np/LW_cnp
    rtU_np = tU_np/tU_cnp
    rtD_np = tD_np/tD_cnp
    ratio_tideduration_np = M2_period_timedelta/(tU_cnp+tD_cnp)

    #vermenigvuldiging van kromme met ratios
    print('vermenigvuldig_kromme')
    prediction_np_corr = vermenigvuldig_kromme(prediction_np, idHW_np, idLW_np, rHW_np, rLW_np, rtD_np, rtU_np)
    print('timeshift')
    prediction_np_corr_timeshift = shift_HW_tostart(prediction_np_corr, idHW_np, tstop_dt, tDiff_np)
    
    print('write to csv')
    prediction_np_corr_timeshift.to_csv(os.path.join(dir_gemgetij,"doodtijkromme_%s_%s.csv"%(current_station,slotGem)),float_format='%.3f',date_format='%Y-%m-%d %H:%M:%S')
    print('plot figure')
    fig,(ax1,ax2) = hatyan.plot_timeseries(ts=prediction_np_corr_timeshift,ts_ext=prediction_np_ext)
    if file_vali_doodtijkromme is not None:
        data_vali_doodtij = pd.read_csv(file_vali_doodtijkromme,index_col=0,parse_dates=True)
        ax1.plot(data_vali_doodtij)
    print('save figure')
    fig.savefig(os.path.join(dir_gemgetij,"springdoodtijkromme_doodtij_%s_%s.png"%(current_station,slotGem)))
    ax1.set_xlim(tstop_dt-dt.timedelta(days=0.5),tstop_dt+dt.timedelta(days=4))
    fig.savefig(os.path.join(dir_gemgetij,"springdoodtijkromme_doodtij_%s_%s_zoom.png"%(current_station,slotGem)))
    print('DOODTIJ finished')
    
    ax_sum.plot(prediction_sp_one['values'],'k',label=f'sp kromme {current_station}, origtiming')
    ax_sum.plot(prediction_np_one['values'],'k',label=f'np kromme {current_station}, origtiming')
    ax_sum.plot(prediction_av_corr['values'],'-',label=f'gem kromme {current_station}, origtiming')
    ax_sum.plot(prediction_sn['values'],':',label=f'spnp kromme {current_station}, origtiming')
    ax_sum.plot(prediction_sp_corr['values'],'--',label=f'sp kromme {current_station}, start0')
    ax_sum.plot(prediction_np_corr['values'],'-.',label=f'np kromme {current_station}, start0')

    #TODO: this is useful to calculate delays between stations, so write to summary DataFrame? This is ruwe sprinneap cycle, so no exact values
    print(f'HWLW gemiddeld tij {current_station}:\n',prediction_av_ext.loc[is1:is2])
    print(f'HWLW springtij {current_station} (tide duration is {tU_csp+tD_csp}):\n',prediction_sn_ext.loc[is1:is2])
    print(f'HWLW doodtij {current_station} (tide duration is {tU_cnp+tD_cnp}):\n',prediction_sn_ext.loc[in1:in2])
    

ax_sum.legend()
ax_sum.grid()
ax_sum.set_xlim(tstop_dt,tstop_dt+dt.timedelta(days=1))
fig_sum.tight_layout()
fig_sum.savefig(os.path.join(dir_gemgetij,'gemgetij_allstations_noshift'))






###OVERSCHRIJDINGSFREQUENTIES
#TODO: discuss edits with RWS:
#       included data up to 2021-1-1 (was 1-1-2012 for all stations) >> trekt weibull krommer omhoog en dichter bij hydraNL
#       included data from max(last3hrint/first1hrint) (HOEKVHLD: 1970-12-31 23:00:00, was 1971-1-1, so almost the same but tiny yaxis offset)
#       break: included data from first10minint in trendanalyis (HOEKVHLD: 1987-01-01 00:10:00, was 1-1-1998 for almost all RMM stations, unless linear trend) >> trekt weibull lijn krommer omlaag en verder van hydraNL
#       max ipv mean in reasampling to hours >> trekt weibull lijn omhoog en dichter bij hydraNL >> combinatie met data tot 2021 is niet per se beter
#       OR: use extremes from RWS dataset >> this also works but there are quite some gaps at some stations (it does sove the need for resampling and mean/max decision, but different tstart/tstop/break_value should still be tested, although it makes sense to include all data and use no break_value since extremes are now independent of meas storing interval) >> then SLR trend removal is even more important
#       OR: derive extremes instead of sampling to 12H (TODO: does hatyan.calc_HWLW work when time interval is not constant? it already crashes on WESTTSLG)
#TODO: Is het niet nodig om een correctie toe te passen voor zeespiegelstijging? >> zorgt misschien voor minder afhankelijkheid van break_value
#       je zou zeespiegelstijging moeten verdisconteren in de (linearie) trend
#       in rapport van HKV wordt gerefereerd naar Goederen(RWS,2003) en daar staat trend in uitgewerkt, die heeft Boyan overgenomen uit HKV rapport.
#       is voor de statistiek wel belangrijk dat die trend wordt verwijderd (anders wordt extreme freqs onderschat), referentievlak is laatste deel meetperiode.
#       je kunt automatisch lineair corrigeren, maar RWS moet akkoord gaan over methodiek want heeft invloed op ontwerpcriteria etc (onder welke condities wel/niet corrigeren).
#TODO: in principe kun je zoveel mogelijk data meenemen, maar je bent hier alleen geinteresseerd in absolute maximum per piek. Je moet bij 3u interval dus wel weten dat dat het maximum is (en dat is het niet), interval van 1u is eigenlijk ook al te weinig, maar soort van acceptabel. >> evt is het beter om alleen de HW's mee te nemen, mits dat echt het maximum representeert (is dat zo: F010, HW en LW uit 1 min. waterhoogten gefilterd uit 10 min. gem.)
#TODO: Je hebt break = '1-1-1998' ingesteld voor dit station, ik heb een maand/jaar gemiddelde bijgevoegd en zie daarin geen reden om dat te doen. Waar komt deze waarde vandaan?
#       break 1998: in buurt van spui was er een trendbreuk, zie rapport boyan, ook voor HOEKVHLD toegepast.
#       die break wordt voor trendanalyse toegepast, daardoor is de trendlijn korter dan de (on)gefilterd lijnen
#weibull lijn begint pas bij hogere freq want die begint pas bij n-de waarde, want die gaat niet met het staartje naar beneden. is niet van toepassing voor die hoogfrequente situaties, is ontwikkeld voor extremen. te voorspellen freqs wordt met np.logspace() opgegeven.
#TODO: zie vragen in script
#TODO: hoe plots beoordelen? >> rode lijn moet soort van verlengde zijn van groene, als die ineens omhoog piekt komt dat door hele extreme wardes die je dan vrmoedelijk ook al ziet in je groene lijn

dir_meas_overschr = os.path.join(dir_base,'data_overschrijding')

#station_break_dict = {'HOEKVHLD':'01-01-1998'} #TODO: possible to make generic?
station_name_dict = {'HOEKVHLD':'Hoek_van_Holland'}

Tfreqs_interested = [5, 2, 1, 1/2, 1/5, 1/10, 1/20, 1/50, 1/100, 1/200,
                     1/500, 1/1000, 1/2000, 1/4000, 1/5000, 1/10000] #TODO: which frequencies are realistic with n years of data? probably remove this entire row >> met 40 jaar data kun je in principe tot 1/40 gaan, maar met weibull kun je extrapoleren en in theorie >> dit is voor tabel die je eruit wil hebben

color_map = {'Ongefilterd':  'b', 'Gefilterd': 'orange', 'Trendanalyse': 'g',
             'Weibull': 'r', 'Hydra-NL': 'm', 'Hydra-NL met modelonzekerheid': 'cyan',
             'Gecombineerd': 'k'}

reproduce_oldoverschr = False

temp = {}
tstarts = pd.DataFrame()
for current_station in []:#stat_list:
    print(f'overschrijdingsfrequenties for {current_station}')

    file_wl_pkl = os.path.join(dir_meas_alldata,f"{current_station}_measwl.pkl")
    data_pd_meas = pd.read_pickle(file_wl_pkl)
    data_pd_meas.index = data_pd_meas.index.tz_localize(None)
    if not reproduce_oldoverschr: #saves time and is not used in old method
        file_ext_pkl = os.path.join(dir_meas_alldata,f"{current_station}_measext.pkl")
        if not os.path.exists(file_ext_pkl):
            continue
        data_pd_measext = pd.read_pickle(file_ext_pkl)
        data_pd_measext.index = data_pd_measext.index.tz_localize(None)
        if len(data_pd_measext['HWLWcode'].unique()) > 2:
            data_pd_HWLW_12 = hatyan.calc_HWLW12345to21(data_pd_measext) #convert 12345 to 12 by taking minimum of 345 as 2 (laagste laagwater) #TODO: this drops first/last value if it is a LW, should be fixed
        else:
            data_pd_HWLW_12 = data_pd_measext.copy()
        data_pd_HW = data_pd_HWLW_12.loc[data_pd_HWLW_12['HWLWcode']==1]
        data_pd_LW = data_pd_HWLW_12.loc[data_pd_HWLW_12['HWLWcode']==2]
        
        #TODO: move this to data-check part (first/last occurrences of WaardeBepalingsmethode)
        data_pd_measext_WBM_tstart = data_pd_measext[['WaardeBepalingsmethode.Code','WaardeBepalingsmethode.Omschrijving']].drop_duplicates(keep='first')
        data_pd_measext_WBM_tstop = data_pd_measext[['WaardeBepalingsmethode.Code','WaardeBepalingsmethode.Omschrijving']].drop_duplicates(keep='last')
        data_pd_measext_WBM_times = pd.concat([data_pd_measext_WBM_tstart,data_pd_measext_WBM_tstop]).sort_index()
    
    #TODO: move this to data-check part (first/last occurrences of WaardeBepalingsmethode)
    data_pd_meas_WBM_tstart = data_pd_meas[['WaardeBepalingsmethode.Code','WaardeBepalingsmethode.Omschrijving']].drop_duplicates(keep='first')
    data_pd_meas_WBM_tstop = data_pd_meas[['WaardeBepalingsmethode.Code','WaardeBepalingsmethode.Omschrijving']].drop_duplicates(keep='last')
    data_pd_meas_WBM_times = pd.concat([data_pd_meas_WBM_tstart,data_pd_meas_WBM_tstop]).sort_index()
    
    #TODO: move this to data-check part (first/last occurrences of unique time intervals), can be use to identify where time duplciates and gaps are
    data_pd_meas['interval'] = pd.Series(data_pd_meas.index,index=data_pd_meas.index).diff()
    data_pd_meas_int_tstart = data_pd_meas[['interval']].drop_duplicates(keep='first')
    data_pd_meas_int_tstop = data_pd_meas[['interval']].drop_duplicates(keep='last')
    data_pd_meas_int_times = pd.concat([data_pd_meas_int_tstart,data_pd_meas_int_tstop]).sort_values(by='interval',ascending=False)
    datetime_last3hrint = data_pd_meas_int_times.loc[data_pd_meas_int_times['interval']==dt.timedelta(hours=3)].index.max() #last 3hr timestep is often also the first 1hr time interval and this prevents accidental 1hr timestaps in the past to beincluded
    datetime_first1hrint = data_pd_meas_int_times.loc[data_pd_meas_int_times['interval']==dt.timedelta(hours=1)].index.min() #first 1hr timestep, sometimes also occurs somewhere in 3hr interval period so is not always a good indicator of where 1hr-interval measurements start
    datetime_first10minint = data_pd_meas_int_times.loc[data_pd_meas_int_times['interval']==dt.timedelta(minutes=10)].index.min() #first 10min timestep
    
    tstart_usefuldata = datetime_last3hrint
    if pd.isnull(tstart_usefuldata): #if not available, revert to actual first 1hr timestep
        tstart_usefuldata = datetime_first1hrint
    if pd.isnull(tstart_usefuldata): #if not available, revert to actual first 1hr timestep
        tstart_usefuldata = datetime_first10minint
    tstarts.loc[current_station,['last3hrint','first1hrint','first10minint']] = [datetime_last3hrint,datetime_first1hrint,datetime_first10minint]
    
    data_pd_meas = data_pd_meas[['values','QC']] # reduces the memory consumption significantly
    data_pd_meas = data_pd_meas.loc[~(data_pd_meas['QC']==99)]
    #data_pd_meas = data_pd_meas*100 #TODO: maybe remove this conversion to cm, but then also fix the fine gridlines in the plot

    station_rule_type = 'break' #TODO: compare results to the ones withouth this break or break on different date
    if reproduce_oldoverschr:
        station_break_value = dt.datetime(1998,1,1)#'01-01-1998' #station_break_dict[current_station] 
        df_alldata = data_pd_meas.resample('H').mean() #TODO: "Rekenkundig gemiddelde waarde over vorige 5 en volgende 5 minuten" >> resampling method moet .max() zijn #TODO: is this resampling method ok (probably means in hour class) or should it be 30min before/after? (tijdcomponent maakt voor fit niet uit)
        df = hatyan.crop_timeseries(ts=df_alldata,times_ext=[tstart_usefuldata,dt.datetime(2012,1,1)]) #available data HOEKVHLD was 1971-1-1 to 2011-12-31 23:50 #TODO: discuss with RWS of deze automatische tstart bepaling acceptabel is
    else:
        #all different
        station_break_value = datetime_first10minint
        df_alldata = data_pd_meas.resample('H').max()
        df = hatyan.crop_timeseries(ts=df_alldata,times_ext=[tstart_usefuldata,dt.datetime(2021,1,1)]) #crop data to data that is uesful for deriving frequencies
        #select yes/no change
        station_break_value = dt.datetime(1998,1,1)
        df_alldata = data_pd_meas.resample('H').mean()
        df = hatyan.crop_timeseries(ts=df_alldata,times_ext=[tstart_usefuldata,dt.datetime(2012,1,1)]) #crop data to data that is uesful for deriving frequencies
    
    """
    statname_overschr = station_name_dict[current_station]
    #metadata_station = dict(metadata.loc[statname_overschr])
    df_old_raw = pd.read_csv(os.path.join(dir_meas_overschr, 'Processed_RWS', f'{statname_overschr}.csv'),
                             sep=';', header=[0], index_col=[0], parse_dates=True)
    df_old = df_old_raw.resample('H').mean()
    diff_array = (df_old/100-df[['values']]).dropna()
    fig,(ax1,ax2) = hatyan.plot_timeseries(ts=df,ts_validation=df_old/100)
    """
    
    # 1. Exceedance
    print('Exceedance')
    dist = {}
    
    print('Calculate unfiltered distribution')
    
    if reproduce_oldoverschr:
        try:
            df_extrema = df.loc[df.resample('12H')[['values']].idxmax().dropna()['values'].values].copy()
        except TypeError as e: #TypeError: The DTypes <class 'numpy.dtype[float64]'> and <class 'numpy.dtype[datetime64]'> do not have a common DType. For example they cannot be stored in a single array unless the dtype is `object`.
            print(f'FAILED: {e}')
            continue #TODO: .idxmax() does not work with non-dropped HOEKVHLD timeseries, why? Also not with cropped HARVT10 timeseries. 
    else:
        df_extrema = hatyan.crop_timeseries(ts=data_pd_HW, times_ext=[tstart_usefuldata,dt.datetime(2012,1,1)],onlyfull=False) #TODO: decide on tstart/tstop with this method
        #hatyan.plot_timeseries(ts=df_extrema, ts_ext=df_extrema)
        if df_extrema.index.max() < station_break_value: #to catch STELLDBTN since extreme data stops after 1996
            continue
    
    dist['Ongefilterd'] = hatyan.distribution(df_extrema.copy())

    """# filtering is only applicable for stations with high river discharge influence, so disabled
    print('Calculate filtered distribution')
    df_peaks, threshold, _ = hatyan.detect_peaks(df_extrema.copy())
    if metadata_station['apply_treshold']: 
        temp[metadata_station['id']] = threshold
        df_extrema_filt = hatyan.filter_with_threshold(df_extrema.copy(), df_peaks, threshold)
    else:
        df_extrema_filt = df_extrema.copy()
    dist['Gefilterd'] = hatyan.distribution(df_extrema_filt.copy())
    """

    print('Calculate filtered distribution with trendanalysis')
    df_trend = hatyan.apply_trendanalysis(df_extrema.copy(),#df_extrema_filt.copy(), #TODO: only starttime 1998-1-1 applied to HOEKVHLD, where to find that information for all coastal stations?
                                          rule_type=station_rule_type,# metadata_station['rule_type'],
                                          rule_value=station_break_value)# metadata_station['rule_value_high'])
    dist['Trendanalyse'] = hatyan.distribution(df_trend.copy())

    print('Fit Weibull to filtered distribution with trendanalysis')
    # Last 100 datapoints from distribution (assuming it is sorted with Tfreqs from large to small)
    dist['Weibull'] = hatyan.get_weibull(dist['Trendanalyse'].copy(),
                                         threshold=dist['Trendanalyse']['values'].iloc[-101],
                                         Tfreqs=np.logspace(-5, np.log10(dist['Trendanalyse']['values_Tfreq'].iloc[-101]), 5000))

    if current_station in station_name_dict.keys(): #TODO: useful validation data, asked Ferdinand if and where this is also available for other kuststations
        stat_name = station_name_dict[current_station]
        print('Load Hydra-NL dstribution data')
        dist['Hydra-NL'] = pd.read_csv(os.path.join(dir_meas_overschr,'Processed_HydraNL','Without_model_uncertainty',f'{stat_name}.csv'), sep=';', header=[0])
        dist['Hydra-NL']['values'] /= 100 # cm to m
        dist['Hydra-NL met modelonzekerheid'] = pd.read_csv(os.path.join(dir_meas_overschr,'Processed_HydraNL','With_model_uncertainty',f'{stat_name}_with_model_uncertainty.csv'), sep=';', header=[0])
        dist['Hydra-NL met modelonzekerheid']['values'] /= 100 # cm to m
    """
    print('Blend trend, weibull and Hydra-NL together')
    dist['Gecombineerd'] = hatyan.blend_distributions(dist['Trendanalyse'].copy(),
                                                      dist['Weibull'].copy(),
                                                      dist['Hydra-NL'].copy())
    if row['apply_treshold']:
        keys = list(dist.keys())
    else:
        keys = [x for x in list(dist.keys()) if x != 'Gefilterd']
    """
    fig, ax = hatyan.plot_distributions(dist, name=current_station,
                                        keys=None, color_map=color_map, legend_loc='lower right',
                                        xlabel='Frequentie [1/jaar]', ylabel='Hoogwater [m+NAP]')
    ax.set_ylim(0,5.5)
    fig.savefig(os.path.join(dir_overschrijding, f'Exceedance_lines_{current_station}.png')) #.svg
    plt.close(fig)
    """
    hatyan.interpolate_interested_Tfreqs_to_csv(dist['Gecombineerd'], Tfreqs=Tfreqs_interested, id=current_station,
                                              csv_dir=dir_overschrijding, prefix='Exceedance_lines')
    """
    continue
    # 2. Deceedance
    print('Deceedance')
    dist = {}

    print('Calculate unfiltered distribution')
    df_extrema = df.loc[df.resample('12H')[['values']].idxmin().dropna()['values'].values]
    dist['Ongefilterd'] = hatyan.distribution(df_extrema.copy(), inverse=True)

    #print('Calculate filtered distribution (direct copy of unfiltered')
    #dist['Gefilterd'] = hatyan.distribution(df_extrema.copy(), inverse=True)

    print('Calculate filtered distribution with trendanalysis')
    df_trend = hatyan.apply_trendanalysis(df_extrema.copy(),
                                          rule_type=station_rule_type,# metadata_station['rule_type'],
                                          rule_value=station_break_value)# metadata_station['rule_value_high'])
    dist['Trendanalyse'] = hatyan.distribution(df_trend.copy(), inverse=True)

    print('Fit Weibull to filtered distribution with trendanalysis')
    # Last 100 datapoints from distribution (assuming it is sorted with Tfreqs from large to small)
    dist['Weibull'] = hatyan.get_weibull(dist['Trendanalyse'].copy(),
                                         threshold=dist['Trendanalyse']['values'].iloc[-100],
                                         Tfreqs=np.logspace(-5, np.log10(dist['Trendanalyse']['values_Tfreq'].iloc[-100]), 5000),
                                         inverse=True)

    """
    print('Blend trend and weibull together')
    dist['Gecombineerd'] = hatyan.blend_distributions(dist['Trendanalyse'].copy(), dist['Weibull'].copy())
    """
    fig, ax = hatyan.plot_distributions(dist, name=current_station,
                                        keys=None,#['Ongefilterd', 'Trendanalyse', 'Weibull', 'Gecombineerd'],
                                        color_map=color_map, legend_loc='upper right',
                                        xlabel='Frequentie [1/jaar]', ylabel='Laagwater [m+NAP]')
    fig.savefig(os.path.join(dir_overschrijding, f'Deceedance_lines_{current_station}.png')) #.svg
    plt.close(fig)
    """
    hatyan.interpolate_interested_Tfreqs_to_csv(dist['Gecombineerd'], Tfreqs=Tfreqs_interested, id=current_station,
                                                csv_dir=dir_overschrijding, prefix='Deceedance_lines')
    """





#report on memory usage
print('memory usage')
var_list = locals().copy() #dir() globals() locals()
var_keys = var_list.keys()
max_size = 0
sum_size = 0
for varname in var_keys:
    var_size = sys.getsizeof(var_list[varname])/1024**2 #in MegaBytes
    max_size = np.maximum(max_size,var_size)
    sum_size += var_size
    if var_size > 5:
        print(f'{varname}: {var_size:.2f}MB')
print(f'max_size: {max_size:.2f}MB')
print(f'sum_size: {sum_size:.2f}MB')




