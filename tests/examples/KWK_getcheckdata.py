# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 11:35:58 2022

@author: veenstra
"""

import os
import pandas as pd
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')
from matplotlib import ticker
import hatyan # available via `pip install hatyan` or at https://github.com/Deltares/hatyan
#import contextily as ctx #`conda install -c conda-forge contextily -y` #commented since not part of hatyan_hmcenv

#TODO: convert to netcdf instead of pkl, think of convenient netcdf format (align with GTSM and DCSM)

get_catalog = False
dataTKdia = True #TODO: communicate data issues to TK (wl and ext): p:\11208031-010-kenmerkende-waarden-k\work\data_vanRWS_20220805\convert_dia2pickle_dataTK.py

tstart_dt_DDL = dt.datetime(1870,1,1) #1870,1,1 for measall folder
tstop_dt_DDL = dt.datetime(2022,1,1)
tzone_DLL = 'UTC+01:00' #'UTC+00:00' for GMT and 'UTC+01:00' for MET
tstart_dt = dt.datetime(2001,1,1)
tstop_dt = dt.datetime(2011,1,1)
NAP2005correction = False #True #TODO: define for all stations
if ((tstop_dt.year-tstart_dt.year)==10) & (tstop_dt.month==tstop_dt.day==tstart_dt.month==tstart_dt.day==1):
    year_slotgem = tstop_dt.year
else:
    year_slotgem = 'invalid'
print(f'year_slotgem: {year_slotgem}')

dir_base = r'p:\11208031-010-kenmerkende-waarden-k\work'
if dataTKdia:
    dir_meas = os.path.join(dir_base,'measurements_wl_18700101_20220101_dataTKdia')
    dir_meas_alldata = os.path.join(dir_base,'measurements_wl_18700101_20220101_dataTKdia')
else:
    dir_meas = os.path.join(dir_base,'measurements_wl_18700101_20220101')
    dir_meas_alldata = os.path.join(dir_base,'measurements_wl_18700101_20220101')
    
dir_meas_DDL = os.path.join(dir_base,f"measurements_wl_{tstart_dt_DDL.strftime('%Y%m%d')}_{tstop_dt_DDL.strftime('%Y%m%d')}")
if not os.path.exists(dir_meas_DDL):
    os.mkdir(dir_meas_DDL)

fig_alltimes_ext = [dt.datetime.strptime(x,'%Y%m%d') for x in os.path.basename(dir_meas_alldata).split('_')[2:4]]

if get_catalog:
    print('retrieving DDL catalog')
    catalog_dict = hatyan.get_DDL_catalog(catalog_extrainfo=['WaardeBepalingsmethoden','MeetApparaten','Typeringen'])
    pd.to_pickle(catalog_dict,os.path.join(dir_base,'DDL_catalog.pkl'))
    print('...done')
else:
    catalog_dict = pd.read_pickle(os.path.join(dir_base,'DDL_catalog.pkl'))
cat_locatielijst = catalog_dict['LocatieLijst']#.set_index('Locatie_MessageID',drop=True)
cat_locatielijst.to_pickle(os.path.join(dir_meas_DDL,'catalog_lokatielijst.pkl'))

#get list of stations with extremes #TODO: before, stations K13A, MAASMSMPL did not have extremes (many kenmerkende waarden are not possible then so skipping is fine)
if dataTKdia:
    cat_aquometadatalijst_sel, cat_locatielijst_sel = hatyan.get_DDL_stationmetasubset(catalog_dict=catalog_dict,station_dict=None, meta_dict={'Grootheid.Code':'WATHTE$','Groepering.Code':'NVT','Hoedanigheid.Code':'NAP'})
    bool_duplicatestatcodes = cat_locatielijst_sel['Code'].duplicated(keep='first')
    cat_locatielijst_sel = cat_locatielijst_sel.loc[~bool_duplicatestatcodes] #drop duplicate station Codes (keep first), just to let it run but this is not desireable. Not a big issue since with TK data, only RDx and RDy are used and probably safe to assume these are equal
else:
    cat_aquometadatalijst_sel, cat_locatielijst_sel = hatyan.get_DDL_stationmetasubset(catalog_dict=catalog_dict,station_dict=None, meta_dict={'Grootheid.Code':'WATHTE$','Groepering.Code':'GETETM2'})
cat_locatielijst_sel['RDx'],cat_locatielijst_sel['RDy'] = hatyan.convert_coordinates(coordx_in=cat_locatielijst_sel['X'].values, coordy_in=cat_locatielijst_sel['Y'].values, epsg_in=int(cat_locatielijst_sel['Coordinatenstelsel'].iloc[0]),epsg_out=28992)
cat_locatielijst_sel_codeidx = cat_locatielijst_sel.reset_index(drop=False).set_index('Code',drop=False)

#stat_name_list = ['BATH','DELFZIJL','DEN HELDER','DORDRECHT','EEMSHAVEN','EURO PLATFORM','HANSWEERT','HARINGVLIETSLUIZEN','HARLINGEN','HOEK VAN HOLLAND','HUIBERTGAT','IJMUIDEN','KORNWERDERZAND','LAUWERSOOG','ROOMPOT BUITEN','ROTTERDAM','SCHEVENINGEN','STAVENISSE','TERNEUZEN','VLISSINGEN','WEST-TERSCHELLING'] # lijst AB
stat_name_list = ['Terneuzen','Bath','HANSWT','Vlissingen','Bergse Diepsluis west','Krammersluizen west','Stavenisse','Roompot binnen','Cadzand','Westkapelle','Roompot buiten','Brouwershavensche Gat 08','Haringvliet 10','Hoek van Holland','Scheveningen','IJmuiden buitenhaven','Petten zuid','Den Helder','Texel Noordzee','Terschelling Noordzee','Wierumergronden','Huibertgat','Oudeschild','Vlieland haven','West-Terschelling','Nes','Schiermonnikoog','Den Oever buiten','Kornwerderzand buiten','Harlingen','Lauwersoog','Eemshaven','Delfzijl','Nieuwe Statenzijl','Lichteiland Goeree','Euro platform','K13a platform'] + ['Dordrecht','Stellendam Buiten','Rotterdam'] + ['Maasmond','Oosterschelde 11'] #+ stat_list_addnonext[2:] #"KW kust en GR Dillingh 2013" en "KW getijgebied RWS 2011.0", aangevuld met 3 stations AB, aangevuld met BOI wensen, aangevuld met dialijst ABCT
stat_list = []
for stat_name in stat_name_list:
    bool_isstation = cat_locatielijst_sel_codeidx['Naam'].str.contains(stat_name,case=False) | cat_locatielijst_sel_codeidx['Code'].str.contains(stat_name,case=False)
    if bool_isstation.sum()!=1:
        print(f'station name {stat_name} found {bool_isstation.sum()} times, should be 1.:\n{cat_locatielijst_sel_codeidx.loc[bool_isstation,["Naam"]]}')
    if bool_isstation.sum()==0: #skip if none found
        continue
    stat_list.append(cat_locatielijst_sel_codeidx.loc[bool_isstation,'Code'].iloc[0])
    #print(f'{stat_name:30s}: {bool_isstation.sum()}')
#stat_list = ['BATH','DELFZL','DENHDR','DORDT','EEMSHVN','EURPFM','HANSWT','STELLDBTN','HARLGN','HOEKVHLD','HUIBGT','IJMDBTHVN','KORNWDZBTN','LAUWOG','ROOMPBTN','ROTTDM','SCHEVNGN','STAVNSE','TERNZN','VLISSGN','WESTTSLG'] # lijst AB vertaald naar DONAR namen
#stat_list = ['HOEKVHLD','HARVT10','VLISSGN']

if dataTKdia:
    stat_list = ['A12','AWGPFM','BAALHK','BATH','BERGSDSWT','BROUWHVSGT02','BROUWHVSGT08','GATVBSLE','BRESKVHVN','CADZD','D15','DELFZL','DENHDR','EEMSHVN','EURPFM','F16','F3PFM','HARVT10','HANSWT','HARLGN','HOEKVHLD','HOLWD','HUIBGT','IJMDBTHVN','IJMDSMPL','J6','K13APFM','K14PFM','KATSBTN','KORNWDZBTN','KRAMMSZWT','L9PFM','LAUWOG','LICHTELGRE','MARLGT','NES','NIEUWSTZL','NORTHCMRT','DENOVBTN','OOSTSDE04','OOSTSDE11','OOSTSDE14','OUDSD','OVLVHWT','Q1','ROOMPBNN','ROOMPBTN','SCHAARVDND','SCHEVNGN','SCHIERMNOG','SINTANLHVSGR','STAVNSE','STELLDBTN','TERNZN','TERSLNZE','TEXNZE','VLAKTVDRN','VLIELHVN','VLISSGN','WALSODN','WESTKPLE','WESTTSLG','WIERMGDN','YERSKE'] #all stations from TK
    stat_list = ['BAALHK','BATH','BERGSDSWT','BRESKVHVN','CADZD','DELFZL','DENHDR','DENOVBTN','EEMSHVN','GATVBSLE','HANSWT','HARLGN','HARVT10','HOEKVHLD','IJMDBTHVN','KATSBTN','KORNWDZBTN','KRAMMSZWT','LAUWOG','OUDSD','ROOMPBNN','ROOMPBTN','SCHAARVDND','SCHEVNGN','SCHIERMNOG','STAVNSE','STELLDBTN','TERNZN','VLAKTVDRN','VLIELHVN','VLISSGN','WALSODN','WESTKPLE','WESTTSLG','WIERMGDN'] #all files with valid data for 2010 to 2021
    #stat_list = stat_list[stat_list.index('STELLDBTN'):]
M2_period_timedelta = pd.Timedelta(hours=hatyan.get_schureman_freqs(['M2']).loc['M2','period [hr]'])


def clean_data(ts_meas_pd,current_station):
    if 'HWLWcode' in ts_meas_pd.columns:
        keep_columns = ['values','QC','HWLWcode']
    else:
        keep_columns = ['values','QC']
    ts_meas_pd = ts_meas_pd[keep_columns] # reduces the memory consumption significantly
    ts_meas_pd.index = ts_meas_pd.index.tz_localize(None)
    ts_meas_pd = ts_meas_pd.loc[~(ts_meas_pd['QC']==99)] #TODO: remove or make nans?
    
    #optional nap correction
    if NAP2005correction:
        ts_meas_pd = nap2005_correction(ts_meas_pd,current_station)
    return ts_meas_pd


def nap2005_correction(data_pd,current_station):
    #NAP correction for dates before 1-1-2005
    #TODO: make this flexible per station, where to get the data or is the RWS data already corrected for it? Also does it matter? for havengetallen it makes a slight difference so yes. For gemgetijkromme it only makes a difference for spring/doodtij. (now only applied at gemgetij en havengetallen)
    #herdefinitie van NAP (2 tot 5 mm, relevant?): https://puc.overheid.nl/PUC/Handlers/DownloadDocument.ashx?identifier=PUC_113484_31&versienummer=1
    #Dit is de rapportage waar het gebruik voor PSMSL data voor het eerst beschreven is: https://puc.overheid.nl/PUC/Handlers/DownloadDocument.ashx?identifier=PUC_137204_31&versienummer=1
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
    
    station_dict = cat_locatielijst_sel_codeidx.loc[current_station,['Locatie_MessageID','X','Y','Naam','Code']]
    
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
                                                        'Hoedanigheid.Code':'NAP',  # Hoedanigheid is necessary for eg EURPFM/LICHTELGRE, where NAP and MSL values are available. #TODO: also look at MSL data? (then duplicate MeetApparaat must be allowed and that is inconvenient as default) >>NAP data begint vanaf 2001 (en bevat ext) en afwezig voor K13APFM, MSL data begint veel eerder (ook tot later?, maar bevat geen ext)
                                                        'MeetApparaat.Code':'127'}) # MeetApparaat.Code is necessary for IJMDBTHVN/ROOMPBTN, where also radar measurements are available (all other stations are vlotter and these stations also have all important data in vlotter) TODO: Except LICHTELGRE/K13APFM which have Radar/Stappenbaak en Radar as MeetApparaat
                                            #Hoedanigheid en MeetApparaat zijn anders voor LICHTELGRE en K13APFM (MSL en variabel)
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
                                                    meta_dict={'Grootheid.Code':'WATHTE','Groepering.Code':'GETETM2'})#,'MeetApparaat.Code':'127'}) #ts_measwlHWLW # TODO: MeetApparaat is necessary for IJMBTHVN/NIEUWSTZL/HOLWD, maybe remove if servicedesk has resolved this probable Vlotter/Radar issue (gemeld op 28-4-2022 voor IJMBTHVN) (or keep and also add Hoedanigheid.Code, alle ext data is toch NAP)
        request_output_exttyp = hatyan.get_DDL_data(station_dict=station_dict,tstart_dt=tstart_dt_DDL,tstop_dt=tstop_dt_DDL,tzone=tzone_DLL, allow_multipleresultsfor=allow_multipleresultsfor,
                                                    meta_dict={'Groepering.Code':'GETETM2','Typering.Code':'GETETTPE'})#,'MeetApparaat.Code':'127'}) #ts_measwlHWLWtype
        if request_output_extval is None:
            continue
        ts_meas_ext_pd, metadata, stationdata = request_output_extval
        ts_meas_exttyp_pd, metadata2, dummy = request_output_exttyp
        ts_meas_ext_pd['values'] = ts_meas_ext_pd['values']/100 #convert from cm to m
        if ts_meas_exttyp_pd['values'].isnull().any(): #TODO: remove this exception for SCHEVNGN after DDL exttype data is fixed
            print(f'WARNING: invalid ext type values for {current_station}, skipping station')    
            continue
        ts_meas_ext_pd = hatyan.convert_HWLWstr2num(ts_meas_ext_pd,ts_meas_exttyp_pd)
        ts_meas_ext_pd.to_pickle(file_wl_pkl.replace('_measwl','_measext'))
        metadata.to_pickle(file_wlmeta_pkl.replace('_measwl','_measext'))






### LOAD DATA FROM PICKLE plot and do checks
#TODO: visually check availability (start/stop/gaps/aggers) of wl/ext, monthmean wl, outliers (nog niet gedaan voor hele periode, wel voor 2000-2022 (listAB+HARVT10):
#   IJMDBTHVN extremen missen vanaf 2018 want Radar ipv Vlotter (al gemeld op 28-4-2022). HOLWD ook vanaf 2012, terwijl measwl allemaal Vlotter is.
#   Missende data vanaf 2000 (gemeld op 26-4):
#       BATH (2000-2020, measwl en measext, komt doordat er twee stations zijn genaamd Bath/BATH) >> andere station bevat wel een goede dataset
#       EURPFM (2000-2001, measwl en measext)
#       HOEKVHLD (2000-2006, 2013, 2019, measext)
#       ROOMPBTN (2016, measext)
#       ROTTDM (2013, 2019, measext)
#       SCHEVNGN (extrementype bevatten nans ipv strings als 'hoogwater', data is dus ongeldig)
#       STELLDBTN (geen data beschikbaar) >> geen ext data vanaf 2000
#   sterke outliers in tijdreeksen (na filtering QC=99, gemeld op 26-4): IJMDBTHVN/ROOMPBTN (2001) >> is niet meer zo na verwijderen Radar metingen
#   >> NIEUWSTZL extremen missen vanaf 2012/2013 want Radar ipv Vlotter, maar hier missen ook de metingen dus is het waarschijnlijker dat ze op Radar over zijn gegaan? (IJMDBTHVN heeft nog wel lang Vlotter measwl data na stoppen van Vlotter extremen)
#TODO: wadden data opvallendheden melden:
#   Outliers HUIBGT 1979/1985/1987 en WIERMGDN 1985/1987 >> en meer jaren
#   TEXNZE 2007/2012/2015: veel grote gaps (2007 heeft 10 maanden gap) >> 2007 is nog steeds het geval
#   HUIBGT 1982: veel ongeldige waardes rond 1.7m >> nog steeds?
#   HUIBGT 2017: veel missende waardes >> nog steeds?
#   UITHZWD1/WIERMWD1 2008 tm 2012 negatieve outliers. Sowieso alle ts vlak aan onderkant door droogval >> niet bij KWK meegenomen?
#TODO: report dubbelingen HARVT10 (2000-2022, al gedaan?) en andere stations (1900-2000), en EURPFM ext, zie data_summary.csv (er zijn ook dubbelingen met nan-waardes)
#TODO: report wl/ext missings in recent period 2000-2021 (vanuit data_summary)
#TODO: vergelijking yearmean wl/HW/LW met validatiedata Anneke (opgevraagd op 28-04-2022) (nu alleen beschikbaar voor HOEKVHLD en HARVT10, sowieso wl is nodig voor slotgemiddelde), it is clear in the HARVT10 figures that something is off for meanwl, dit gebeurt misschien ook bij andere stations met duplicate times in data_summary_filtered.xlsx (also check on nanvalues that are not nan in validationdata, this points to missing data in DDL)
#TODO: DORDT getijslag lijkt in 1970 ineens kleiner te worden, is dit geen foute dataset? >> deltawerken?

"""
#TODO: controleren of andere datasets nuttige data bevatten nadat gemiddelde HW/LW/wl uitwijzen dat ze meer data bevatten? (onderstaande komt uit data_summary.T)
	BATH (Vlotter en NAP voor ext)
DDL_MeetApparaat.Code_wl	127|155
DDL_MeetApparaat.Omschrijving_wl	Vlotter|Druksensor
DDL_Hoedanigheid.Code_wl	NAP
	EURPFM (Vlotter en NAP voor ext)
DDL_MeetApparaat.Code_wl	125|127
DDL_MeetApparaat.Omschrijving_wl	Stappenbaak|Vlotter
DDL_Hoedanigheid.Code_wl	MSL|NAP
	IJMDBTHVN (Radar/Vlotter en NAP voor ext, Radar is tijdelijke meting voor wl maar lijkt of ext per ongeluk ook zo zijn geregistreerd) >> gemeld
DDL_MeetApparaat.Code_wl	109|127
DDL_MeetApparaat.Omschrijving_wl	Radar|Vlotter
DDL_Hoedanigheid.Code_wl	NAP
	K13APFM (ext ontbreekt)
DDL_MeetApparaat.Code_wl	109
DDL_MeetApparaat.Omschrijving_wl	Radar
DDL_Hoedanigheid.Code_wl	MSL
	LICHTELGRE (Radar en NAP voor ext)
DDL_MeetApparaat.Code_wl	109|125
DDL_MeetApparaat.Omschrijving_wl	Radar|Stappenbaak
DDL_Hoedanigheid.Code_wl	MSL|NAP
	NES (Vlotter en NAP voor ext)
DDL_MeetApparaat.Code_wl	109|127
DDL_MeetApparaat.Omschrijving_wl	Radar|Vlotter
DDL_Hoedanigheid.Code_wl	NAP
	NIEUWSTZL (Radar/Vlotter en NAP voor ext)
DDL_MeetApparaat.Code_wl	109|127
DDL_MeetApparaat.Omschrijving_wl	Radar|Vlotter
DDL_Hoedanigheid.Code_wl	NAP
	ROOMPBTN (Vlotter en NAP voor ext, Radar is tijdelijke meting) >> goed zo dus
DDL_MeetApparaat.Code_wl	109|127
DDL_MeetApparaat.Omschrijving_wl	Radar|Vlotter
DDL_Hoedanigheid.Code_wl	NAP
	VLISSGN (Vlotter en NAP voor ext)
DDL_MeetApparaat.Code_wl	124|127
DDL_MeetApparaat.Omschrijving_wl	Peilschaal|Vlotter
DDL_Hoedanigheid.Code_wl	NAP
	WESTKPLE (Vlotter en NAP voor ext)
DDL_MeetApparaat.Code_wl	109|127
DDL_MeetApparaat.Omschrijving_wl	Radar|Vlotter
DDL_Hoedanigheid.Code_wl	NAP

MAASMSMPL wl >> geen vlotter (geen ext data beschikbaar)
"""

"""
#TODO: melden servicedesk data: zes duplicate timesteps in extremen aanwezig met gelijke waarden EURPFM en NIEUWSTZL (laatste van ander MeetApparaat)
ts_meas_ext_pd.loc[ts_meas_ext_pd.index.duplicated(keep=False),['values','QC','Status','HWLWcode']].sort_index()
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
                           values  QC  ... MeetApparaat.Omschrijving HWLWcode
Tijdstip                               ...                                   
2012-12-31 09:20:00+01:00   -0.15   0  ...                   Vlotter        2
2012-12-31 09:20:00+01:00   -0.15   0  ...                     Radar        2
2012-12-31 14:04:00+01:00    1.47   0  ...                   Vlotter        1
2012-12-31 14:04:00+01:00    1.47   0  ...                     Radar        1
2012-12-31 21:28:00+01:00   -0.89   0  ...                   Vlotter        2
2012-12-31 21:28:00+01:00   -0.89   0  ...                     Radar        2
2013-01-01 02:30:00+01:00    1.79   0  ...                   Vlotter        1
2013-01-01 02:30:00+01:00    1.79   0  ...                     Radar        1
2013-01-01 09:50:00+01:00   -0.78   0  ...                   Vlotter        2
2013-01-01 09:50:00+01:00   -0.78   0  ...                     Radar        2
2013-01-01 14:31:00+01:00    2.17   0  ...                   Vlotter        1
2013-01-01 14:31:00+01:00    2.17   0  ...                     Radar        1
[12 rows x 8 columns]
"""
"""
#gemeld op 28-4-2022 bij servicedesk data: Radar extremen IJMDBTHVN vanaf 2018 (dus missings) #TODO: is ook het geval voor NIEUWSTZL
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
"""
not in M2phasediff document: ['LICHTELGRE','EURPFM']
HW/LW numbers not always increasing: ['HANSWT','BROUWHVSGT08','PETTZD','DORDT']
no extremes in requested time frame: ['STELLDBTN','OOSTSDE11']
Catalog query yielded no results (no ext available like K13APFM): A12
"""
data_summary = pd.DataFrame(index=stat_list).sort_index()
for current_station in []:#stat_list:
    print(f'checking data for {current_station}')
    list_relevantmetadata = ['WaardeBepalingsmethode.Code','WaardeBepalingsmethode.Omschrijving','MeetApparaat.Code','MeetApparaat.Omschrijving','Hoedanigheid.Code','Grootheid.Code','Groepering.Code','Typering.Code']
    list_relevantDDLdata = ['WaardeBepalingsmethode.Code','MeetApparaat.Code','MeetApparaat.Omschrijving','Hoedanigheid.Code']
    
    if not dataTKdia:
        station_dict = dict(cat_locatielijst_sel_codeidx.loc[current_station,['Naam','Code']]) #TODO: put comment in hatyan.getonlinedata.py: get_DDL_stationmetasubset() does not work if 'X','Y','Locatie_MessageID' is added, since there is no column with that name (is index) and if it is, it is an int and not a str
        cat_aquometadatalijst_temp, cat_locatielijst_temp = hatyan.get_DDL_stationmetasubset(catalog_dict=catalog_dict,station_dict=station_dict,meta_dict={'Grootheid.Code':'WATHTE','Groepering.Code':'NVT'})
        for metakey in list_relevantDDLdata:
            data_summary.loc[current_station,f'DDL_{metakey}_wl'] = '|'.join(cat_aquometadatalijst_temp[metakey].unique())
        if not current_station in ['K13APFM','MAASMSMPL']:# no ext available for these stations
            cat_aquometadatalijst_temp, cat_locatielijst_temp = hatyan.get_DDL_stationmetasubset(catalog_dict=catalog_dict,station_dict=station_dict,meta_dict={'Grootheid.Code':'WATHTE','Groepering.Code':'GETETM2'})
            for metakey in list_relevantDDLdata:
                data_summary.loc[current_station,f'DDL_{metakey}_ext'] = '|'.join(cat_aquometadatalijst_temp[metakey].unique())
    
    #add coordinates to data_summary
    data_summary.loc[current_station,['RDx','RDy']] = cat_locatielijst_sel_codeidx.loc[current_station,['RDx','RDy']]
    time_interest_start = dt.datetime(2000,1,1)
    time_interest_stop = dt.datetime(2021,2,1)
    
    #load measwl data
    file_wl_pkl = os.path.join(dir_meas_alldata,f"{current_station}_measwl.pkl")
    file_wlmeta_pkl = os.path.join(dir_meas_alldata,f"meta_{current_station}_measwl.pkl")
    if not os.path.exists(file_wl_pkl):
        data_summary.loc[current_station,'data_wl'] = False
        data_summary.loc[current_station,'data_ext'] = False
        continue
    data_summary.loc[current_station,'data_wl'] = True
    ts_meas_pd = pd.read_pickle(file_wl_pkl)
    if not dataTKdia:
        metawl = pd.read_pickle(file_wlmeta_pkl)
        for metakey in list_relevantmetadata:
            data_summary.loc[current_station,f'{metakey}_wl'] = '|'.join(metawl[metakey].unique())
    ts_meas_pd = ts_meas_pd[['values','QC']] # reduces the memory consumption significantly
    if str(ts_meas_pd.index[0].tz) != 'Etc/GMT-1': #this means UTC+1
        raise Exception(f'measwl data for {current_station} is not in expected timezone (Etc/GMT-1): {ts_meas_pd.index[0].tz}')
    ts_meas_pd.index = ts_meas_pd.index.tz_localize(None)
    bool_99 = ts_meas_pd['QC']==99
    if bool_99.any(): #ts contains invalid values
        ts_meas_pd[bool_99] = np.nan
    data_summary.loc[current_station,'tstart_wl'] = ts_meas_pd.index[0]
    data_summary.loc[current_station,'tstop_wl'] = ts_meas_pd.index[-1]
    data_summary.loc[current_station,'tstart2000_wl'] = ts_meas_pd.index[0]<=time_interest_start
    data_summary.loc[current_station,'tstop202102_wl'] = ts_meas_pd.index[-1]>=time_interest_stop
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
    ts_meas_2000to202102a = ts_meas_pd.loc[~ts_meas_dupltimes,['values']].loc[time_interest_start:min(ts_meas_pd.index[-1],time_interest_stop)]
    ts_meas_2000to202102b = pd.DataFrame({'values':ts_meas_pd.loc[~ts_meas_dupltimes,'values']},index=pd.date_range(start=time_interest_start,end=time_interest_stop,freq='10min'))
    data_summary.loc[current_station,'#nans_2000to202102a_wl'] = ts_meas_2000to202102a['values'].isnull().sum()
    data_summary.loc[current_station,'#nans_2000to202102b_wl'] = ts_meas_2000to202102b['values'].isnull().sum()
    
    #calculate monthly/yearly mean for meas wl data #TODO: use hatyan.calc_wltidalindicators() instead (with threshold of eg 2900 like slotgem)
    mean_peryearmonth_long = ts_meas_pd.groupby(pd.PeriodIndex(ts_meas_pd.index, freq="M"))['values'].mean()
    data_summary.loc[current_station,'monthmean_mean_wl'] = mean_peryearmonth_long.mean()
    data_summary.loc[current_station,'monthmean_std_wl'] = mean_peryearmonth_long.std()
    mean_peryear_long = ts_meas_pd.groupby(pd.PeriodIndex(ts_meas_pd.index, freq="Y"))['values'].mean()
    data_summary.loc[current_station,'yearmean_mean_wl'] = mean_peryear_long.mean()
    data_summary.loc[current_station,'yearmean_std_wl'] = mean_peryear_long.std()
    """#TODO: move to hatyan.timeseries.* or hatyan.kenmerkendewaarden.*. Add minimum # values to calculate monthmean? Make long2array edit simpler with pandas smart stuff?
    numvals_peryearmonth_long = ts_meas_pd.groupby(pd.PeriodIndex(ts_meas_pd.index, freq="M"))['values'].count()
    mean_peryearmonth_array = pd.DataFrame(index=range(1,13))
    for year in mean_peryearmonth_long.index.year.unique():
        bool_year = mean_peryearmonth_long.index.year==year
        mean_peryearmonth_long_oneyr = mean_peryearmonth_long.loc[bool_year]
        mean_peryearmonth_array.loc[mean_peryearmonth_long_oneyr.index.month,year] = mean_peryearmonth_long_oneyr.values
    mean_permonth = mean_peryearmonth_array.mean(axis=1)
    """

    #load measext data
    file_ext_pkl = os.path.join(dir_meas_alldata,f"{current_station}_measext.pkl")
    file_extmeta_pkl = os.path.join(dir_meas_alldata,f"meta_{current_station}_measext.pkl")
    if not os.path.exists(file_ext_pkl):
        data_summary.loc[current_station,'data_ext'] = False
    else:
        data_summary.loc[current_station,'data_ext'] = True
        ts_meas_ext_pd = pd.read_pickle(file_ext_pkl)
        timediff_ext = ts_meas_ext_pd.index[1:]-ts_meas_ext_pd.index[:-1]
        if timediff_ext.min() < dt.timedelta(hours=4): #TODO: min timediff for e.g. BROUWHVSGT08 is 3 minutes: ts_meas_ext_pd.loc[dt.datetime(2015,1,1):dt.datetime(2015,1,2),['values', 'QC', 'Status']]. This should not happen and with new dataset should be converted to an error
            print(f'WARNING: extreme data contains values that are too close ({timediff_ext.min()}), should be at least 4 hours difference')
        if not dataTKdia:
            metaext = pd.read_pickle(file_extmeta_pkl)
            for metakey in list_relevantmetadata:
                data_summary.loc[current_station,f'{metakey}_ext'] = '|'.join(metaext[metakey].unique())
            if str(ts_meas_ext_pd.index[0].tz) != 'Etc/GMT-1': #this means UTC+1
                raise Exception(f'measext data for {current_station} is not in expected timezone (Etc/GMT-1): {ts_meas_ext_pd.index[0].tz}')
        ts_meas_ext_pd.index = ts_meas_ext_pd.index.tz_localize(None)
        ts_meas_ext_dupltimes = ts_meas_ext_pd.index.duplicated()
        data_summary.loc[current_station,'mintimediff_ext'] = timediff_ext.min()
        data_summary.loc[current_station,'dupltimes_ext'] = ts_meas_ext_dupltimes.sum()
        data_summary.loc[current_station,'tstart_ext'] = ts_meas_ext_pd.index[0]
        data_summary.loc[current_station,'tstop_ext'] = ts_meas_ext_pd.index[-1]
        data_summary.loc[current_station,'tstart2000_ext'] = ts_meas_ext_pd.index[0]<=(time_interest_start+M2_period_timedelta)
        data_summary.loc[current_station,'tstop202102_ext'] = ts_meas_ext_pd.index[-1]>=(time_interest_stop-M2_period_timedelta)
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
            ts_meas_ext_2000to202102 = ts_meas_ext_pd.loc[(ts_meas_ext_pd.index>=time_interest_start) & (ts_meas_ext_pd.index<=time_interest_stop)]
            ts_meas_ext_pd.loc[dt.datetime(2015,1,1):dt.datetime(2015,1,2)]
            ts_meas_ext_2000to202102 = hatyan.calc_HWLWnumbering(ts_meas_ext_2000to202102)
            HWmissings = (ts_meas_ext_2000to202102.loc[ts_meas_ext_pd['HWLWcode']==1,'HWLWno'].diff().dropna()!=1).sum()
            data_summary.loc[current_station,'#HWgaps_2000to202102_ext'] = HWmissings
        except Exception as e: #"tidal wave numbering: HW/LW numbers not always increasing" and "zero-size array to reduction operation minimum which has no identity" #TODO: fix by calulate and providing station or corr_tideperiods argument? Or fix otherwise in hatyan (maybe under different project)
            print(f'ERROR: {e}')
        
        #calculate monthly/yearly mean for meas ext data
        if len(ts_meas_ext_pd['HWLWcode'].unique()) > 2:
            data_pd_HWLW_12 = hatyan.calc_HWLW12345to12(ts_meas_ext_pd) #convert 12345 to 12 by taking minimum of 345 as 2 (laagste laagwater). TODO: currently, first/last values are skipped if LW
        else:
            data_pd_HWLW_12 = ts_meas_ext_pd.copy()
        data_pd_HW = data_pd_HWLW_12.loc[data_pd_HWLW_12['HWLWcode']==1]
        data_pd_LW = data_pd_HWLW_12.loc[data_pd_HWLW_12['HWLWcode']==2]
        HW_mean_peryear_long = data_pd_HW.groupby(pd.PeriodIndex(data_pd_HW.index, freq="y"))['values'].mean() #TODO: use hatyan.calc_HWLWtidalindicators() instead (with threshold of eg 1400 like slotgem)
        LW_mean_peryear_long = data_pd_LW.groupby(pd.PeriodIndex(data_pd_LW.index, freq="y"))['values'].mean()
        
    if data_summary['data_ext'].isnull().sum() == 0: #if all stat_list stations were processed (only True/False in this array, no nans)
        #print and save data_summary
        print(data_summary[['data_wl','tstart_wl','tstop_wl','nvals_wl','dupltimes_wl','#nans_wl','#nans_2000to202102a_wl']])
        print(data_summary[['data_ext','dupltimes_ext','#HWgaps_2000to202102_ext']])
        data_summary.to_csv(os.path.join(dir_meas_alldata,'data_summary.csv'))
        
        #make spatial plot of available/retrieved stations
        fig_map,ax_map = plt.subplots(figsize=(8,7))
        file_ldb = r'p:\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\computations\model_setup\run_205\20101209-06.ldb' #TODO: make ldb available in code or at least KWK project drive
        if os.path.exists(file_ldb):
            ldb_pd = pd.read_csv(file_ldb, delim_whitespace=True,skiprows=4,names=['RDx','RDy'],na_values=[999.999])
            ax_map.plot(ldb_pd['RDx'],ldb_pd['RDy'],'-k',linewidth=0.4)
        ax_map.plot(cat_locatielijst_sel['RDx'],cat_locatielijst_sel['RDy'],'xk')#,alpha=0.4) #all ext stations
        ax_map.plot(cat_locatielijst_sel_codeidx.loc[stat_list,'RDx'],cat_locatielijst_sel_codeidx.loc[stat_list,'RDy'],'xr') # selected ext stations (stat_list)
        ax_map.plot(data_summary.loc[data_summary['data_ext'],'RDx'],data_summary.loc[data_summary['data_ext'],'RDy'],'xm') # data retrieved
        """
        for iR, row in cat_locatielijst_sel.iterrows():
            ax_map.text(row['RDx'],row['RDy'],row['Code'])
        """
        ax_map.set_xlim(-50000,300000)
        ax_map.set_ylim(350000,650000)
        ax_map.set_title('overview of stations with GETETM2 data')
        ax_map.set_aspect('equal')
        def div1000(x,pos): return f'{int(x//1000)}'
        ax_map.xaxis.set_major_formatter(ticker.FuncFormatter(div1000))
        ax_map.yaxis.set_major_formatter(ticker.FuncFormatter(div1000))
        ax_map.set_xlabel('RDx [km]')
        ax_map.set_ylabel('RDy [km]')
        ax_map.grid(alpha=0.5)
        fig_map.tight_layout()
        #ctx.add_basemap(ax_map, source=ctx.providers.Esri.WorldImagery, crs="EPSG:28992", attribution=False)
        fig_map.savefig(os.path.join(dir_meas_alldata,'stations_map.png'))
    
    #plotting
    file_wl_png = os.path.join(dir_meas_alldata,f'ts_{current_station}.png')
    if 0:#os.path.exists(file_wl_png):
        continue #skip the plotting if there is already a png available
    if os.path.exists(file_ext_pkl):
        fig,(ax1,ax2) = hatyan.plot_timeseries(ts=ts_meas_pd, ts_ext=ts_meas_ext_pd)
    else:
        fig,(ax1,ax2) = hatyan.plot_timeseries(ts=ts_meas_pd)
    ax1.set_title(f'timeseries for {current_station}')
    ax1_legendlabels = ax1.get_legend_handles_labels()[1]
    ax2_legendlabels = ['zero']
    ax1_legendlabels.insert(1,'zero') #legend for zero line was not displayed but will be now so it needs to be added
    ax1_legendlabels[0] = 'measured waterlevels'
    ax1_legendlabels[2] = 'mean'
    ax1.plot(mean_peryearmonth_long,'c',linewidth=0.7); ax1_legendlabels.append('monthly mean')
    ax1.plot(mean_peryear_long,'m',linewidth=0.7); ax1_legendlabels.append('yearly mean')
    ax2.plot(mean_peryearmonth_long,'c',linewidth=0.7); ax2_legendlabels.append('monthly mean')
    ax2.plot(mean_peryear_long,'m',linewidth=0.7); ax2_legendlabels.append('yearly mean')
    ax1.set_ylim(-4,4)
    ax1.legend(ax1_legendlabels,loc=4)
    ax2.legend(ax2_legendlabels,loc=1)
    if os.path.exists(file_ext_pkl): #plot after legend creation, so these entries are not included
        ax1.plot(HW_mean_peryear_long,'m',linewidth=0.7)#; ax1_legendlabels.append('yearly mean')
        ax1.plot(LW_mean_peryear_long,'m',linewidth=0.7)#; ax1_legendlabels.append('yearly mean')
    ax2.set_ylim(-0.5,0.5)
    ax1.set_xlim(fig_alltimes_ext) # entire period
    fig.savefig(file_wl_png.replace('.png','_alldata.png'))
    ax1.set_xlim(dt.datetime(2000,1,1),dt.datetime(2022,1,1)) # period of interest
    fig.savefig(file_wl_png)
    plt.close(fig)
