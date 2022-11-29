# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 14:17:13 2022

@author: veenstra
"""

import os
import sys
import pandas as pd
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')
from matplotlib import ticker
import hatyan # available via `pip install hatyan` or at https://github.com/Deltares/hatyan
import contextily as ctx #`conda install -c conda-forge contextily -y`

#TODO: SLR trend correctie voor overschrijdingsfrequenties en evt ook voor andere KW?
#TODO: move all parts to hatyan.kenmerkendewaarden.*, maybe also the stuff in hatyan/overschrijding.py (and include license header)
#TODO: add LAT/HAT (AB needs this for RWS work)
#TODO: add tidal coefficient?: The tidal coefficient is the size of the tide in relation to its mean. It usually varies between 20 and 120. The higher the tidal coefficient, the larger the tidal range – i.e. the difference in water height between high and low tide. This means that the sea level rises and falls back a long way. The mean value is 70. We talk of strong tides – called spring tides – from coefficient 95.  Conversely, weak tides are called neap tides. https://escales.ponant.com/en/high-low-tide/ en https://www.manche-toerisme.com/springtij
get_catalog = False
dataTKdia = True
closefigatstart = True

tstart_dt_DDL = dt.datetime(1870,1,1) #1870,1,1 for measall folder
tstop_dt_DDL = dt.datetime(2022,1,1)
tzone_DLL = 'UTC+01:00' #'UTC+00:00' for GMT and 'UTC+01:00' for MET
tstart_dt = dt.datetime(2011,1,1)
tstop_dt = dt.datetime(2021,1,1)
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
    dir_meas = os.path.join(dir_base,'measurements_wl_20000101_20220101')
    dir_meas_alldata = os.path.join(dir_base,'measurements_wl_18700101_20220101')
    
dir_meas_DDL = os.path.join(dir_base,f"measurements_wl_{tstart_dt_DDL.strftime('%Y%m%d')}_{tstop_dt_DDL.strftime('%Y%m%d')}")
if not os.path.exists(dir_meas_DDL):
    os.mkdir(dir_meas_DDL)
dir_havget = os.path.join(dir_base,f'out_havengetallen_{year_slotgem}')
if not os.path.exists(dir_havget):
    os.mkdir(dir_havget)
dir_slotgem = os.path.join(dir_base,f'out_slotgem_{year_slotgem}')
if not os.path.exists(dir_slotgem):
    os.mkdir(dir_slotgem)
dir_gemgetij = os.path.join(dir_base,f'out_gemgetij_{year_slotgem}')
if not os.path.exists(dir_gemgetij):
    os.mkdir(dir_gemgetij)
dir_overschrijding = os.path.join(dir_base,f'out_overschrijding_{year_slotgem}')
if not os.path.exists(dir_overschrijding):
    os.mkdir(dir_overschrijding)

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


def clean_data(ts_meas_pd):
    if 'HWLWcode' in ts_meas_pd.columns:
        keep_columns = ['values','QC','HWLWcode']
    else:
        keep_columns = ['values','QC']
    ts_meas_pd = ts_meas_pd[keep_columns] # reduces the memory consumption significantly
    ts_meas_pd.index = ts_meas_pd.index.tz_localize(None)
    ts_meas_pd = ts_meas_pd.loc[~(ts_meas_pd['QC']==99)] #TODO: remove or make nans?
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
    
    #calculate monthly/yearly mean for meas wl data #TODO: use hatyan.calc_wltidalindicators() instead (with threshold of eg 2900 like slotgem
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
        ctx.add_basemap(ax_map, source=ctx.providers.Esri.WorldImagery, crs="EPSG:28992", attribution=False)
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
    





#TODO: more data is needed for proper working of models (2011: BAALHK, BRESKVHVN, GATVBSLE, SCHAARVDND)
#### SLOTGEMIDDELDEN
physical_break_dict = {'DENOVBTN':1933, #laatste sluitgat afsluitdijk in 1932
                       'HARLGN':1933, #laatste sluitgat afsluitdijk in 1932
                       'VLIELHVN':1933, #laatste sluitgat afsluitdijk in 1932
                       } #TODO: add physical_break for STAVNSE and KATSBTN? (Oosterscheldekering)
for current_station in []:#stat_list:#
    
    if closefigatstart:
        plt.close('all')
    print(f'slotgemiddelden for {current_station}')
    
    #derive yearmean wl from wl values
    file_wl_pkl = os.path.join(dir_meas_alldata,f"{current_station}_measwl.pkl")
    if not os.path.exists(file_wl_pkl):
        continue
    data_pd_meas = pd.read_pickle(file_wl_pkl)
    data_pd_meas = clean_data(data_pd_meas)
    
    #calculate yearly mean
    dict_wltidalindicators = hatyan.calc_wltidalindicators(data_pd_meas)
    wl_mean_peryear = dict_wltidalindicators['wl_mean_peryear']
    dict_wltidalindicators_valid = hatyan.calc_wltidalindicators(data_pd_meas, tresh_yearlywlcount=2900) #24*365=8760 (hourly interval), 24/3*365=2920 (3-hourly interval)
    wl_mean_peryear_valid = dict_wltidalindicators_valid['wl_mean_peryear']

    #derive tidal indicators like yearmean HWLW from HWLW values
    file_ext_pkl = os.path.join(dir_meas_alldata,f"{current_station}_measext.pkl")
    if os.path.exists(file_ext_pkl):
        data_pd_HWLW_all = pd.read_pickle(file_ext_pkl)
        dict_HWLWtidalindicators = hatyan.calc_HWLWtidalindicators(data_pd_HWLW_all)
        HW_mean_peryear = dict_HWLWtidalindicators['HW_mean_peryear']
        LW_mean_peryear = dict_HWLWtidalindicators['LW_mean_peryear']
        dict_HWLWtidalindicators_valid = hatyan.calc_HWLWtidalindicators(data_pd_HWLW_all, tresh_yearlyHWLWcount=1400) #2*24*365/12.42=1410.6 (12.42 hourly extreme)
        HW_mean_peryear_valid = dict_HWLWtidalindicators_valid['HW_mean_peryear']
        LW_mean_peryear_valid = dict_HWLWtidalindicators_valid['LW_mean_peryear']
    
    #plotting (yearly averages are plotted on 1jan, would be better on 1jul)
    fig,ax1 = plt.subplots(figsize=(14,7))
    
    #get validation timeseries (yearly mean wl/HW/LW)
    station_name_dict = {'HOEKVHLD':'hoek', #TODO: request data for all stations with DONARcode in filename (is not a standard RWS procedure)
                         'HARVT10':'ha10'}
    if current_station in station_name_dict.keys():
        dir_meas_gemHWLWwlAB = r'p:\11208031-010-kenmerkende-waarden-k\work\data_KW-RMM'
        file_yearmeanHW = os.path.join(dir_meas_gemHWLWwlAB,f'{station_name_dict[current_station]}_hw.txt')
        file_yearmeanLW = os.path.join(dir_meas_gemHWLWwlAB,f'{station_name_dict[current_station]}_lw.txt')
        file_yearmeanwl = os.path.join(dir_meas_gemHWLWwlAB,f'{station_name_dict[current_station]}_Z.txt')
        yearmeanHW = pd.read_csv(file_yearmeanHW, delim_whitespace=True, skiprows=1, names=['datetime','values'], parse_dates=['datetime'], na_values=-999.9, index_col='datetime')/100
        yearmeanLW = pd.read_csv(file_yearmeanLW, delim_whitespace=True, skiprows=1, names=['datetime','values'], parse_dates=['datetime'], na_values=-999.9, index_col='datetime')/100
        yearmeanwl = pd.read_csv(file_yearmeanwl, delim_whitespace=True, skiprows=1, names=['datetime','values'], parse_dates=['datetime'], na_values=-999.9, index_col='datetime')/100
        ax1.plot(yearmeanHW['values'],'+g')
        ax1.plot(yearmeanLW['values'],'+g')
        ax1.plot(yearmeanwl['values'],'+g')

    if os.path.exists(file_ext_pkl):
        ax1.plot(HW_mean_peryear,'x',color='grey')
        ax1.plot(LW_mean_peryear,'x',color='grey')
        ax1.plot(HW_mean_peryear_valid,'xr')
        ax1.plot(LW_mean_peryear_valid,'xr')
    ax1.plot(wl_mean_peryear,'x',color='grey')
    ax1.plot(wl_mean_peryear_valid,'xr')
    ax1.grid()
    ax1.set_xlim(fig_alltimes_ext) # entire period
    ax1.set_ylabel('waterstand [m]')
    ax1.set_title(f'yearly mean HW/wl/LW {current_station}')
    fig.tight_layout()
    
    if os.path.exists(file_ext_pkl):
        mean_list = [wl_mean_peryear_valid,HW_mean_peryear_valid,LW_mean_peryear_valid]
    else:
        mean_list = [wl_mean_peryear]
    for iM, mean_array in enumerate(mean_list):
        if current_station in physical_break_dict.keys():
            tstart_year_trend = physical_break_dict[current_station]
            tstart_dt_trend = dt.datetime(tstart_year_trend,1,1) #TODO: instead of removing pre-break data, maybe use broken_linear_model and set some years to nan? (more complex method might not be beneficial)
        else:
            tstart_year_trend = None
            tstart_dt_trend = None
        tstop_dt_trend = tstop_dt-dt.timedelta(days=1)
        mean_array_todate = mean_array.loc[tstart_dt_trend:tstop_dt_trend] #remove all values after tstop_dt (is year_slotgem)
        
        # We'll just use the years. This assumes that annual waterlevels are used that are stored left-padded, the mean waterlevel for 2020 is stored as 2020-1-1. This is not logical, but common practice.
        allyears_DTI = pd.date_range(mean_array_todate.index.min(),mean_array_todate.index.max()+dt.timedelta(days=5*360),freq='AS')
        mean_array_allyears = pd.Series(mean_array_todate,index=allyears_DTI)
        
        df = pd.DataFrame({'year':mean_array_allyears.index.year, 'height':mean_array_allyears.values}) #TODO: make functions accept mean_array instead of df as argument?
        
        # below methods are copied from https://github.com/openearth/sealevel/blob/master/slr/slr/models.py #TODO: install slr package as dependency or keep separate?
        #fit, names, X = hatyan.broken_linear_model(df, with_wind=False,start_acceleration=tstart_year_trend)
        #pred_broken_linear = fit.predict(X)
        fit, names, X = hatyan.linear_model(df, with_wind=False, with_nodal=False)
        pred_linear_nonodal = fit.predict(X)
        fit, names, X = hatyan.linear_model(df, with_wind=False)
        pred_linear_winodal = fit.predict(X)
        
        pred_pd = pd.DataFrame({#'pred_linear_acceleration':pred_linear_acceleration,
                                #'pred_quadratic':pred_quadratic,
                                #'pred_broken_linear':pred_broken_linear,
                                'pred_linear_nonodal':pred_linear_nonodal,
                                'pred_linear_winodal':pred_linear_winodal},
                                index=allyears_DTI)
        ax1.plot(pred_pd, ".-", label=pred_pd.columns)
        ax1.set_prop_cycle(None) #reset matplotlib colors
        
        #2021.0 value
        if iM==0: #only for 
            pred_slotgem = pred_pd.loc[[tstop_dt]]
            pred_slotgem.to_csv(os.path.join(dir_slotgem,f'slotgem_value_{current_station}.txt'))
        pred_future = pred_pd.loc[tstop_dt:,'pred_linear_winodal']
        ax1.plot(pred_future, ".k", label=f'pred_linear from {year_slotgem}')
        
        
    ax1.legend(loc=2)
    fig.savefig(os.path.join(dir_slotgem,f'yearly_values_{current_station}'))





#TODO IMPORTANT: provide feedback on incorrect values in extreme timeseries (include in generic data edits?)
#TODO IMPORTANT: check culm_addtime and HWLWno+4 offsets. culm_addtime could also be 2 days or 2days +1h GMT-MET correction. 20 minutes seems odd since moonculm is about tidal wave from ocean
### HAVENGETALLEN
culm_addtime = 2*dt.timedelta(hours=24,minutes=50)-dt.timedelta(minutes=20)+dt.timedelta(hours=1) # 2d and 2u20min correction, this shifts the x-axis of aardappelgrafiek: HW is 2 days after culmination (so 4x25min difference between length of avg moonculm and length of 2 days), 20 minutes (0 to 5 meridian), 1 hour (GMT to MET) #TODO: do we really want to correct for all this now moonculm and HWLW are matched via HWLWno?
data_pd_moonculm = hatyan.astrog_culminations(tFirst=tstart_dt-culm_addtime-dt.timedelta(hours=2*24),tLast=tstop_dt,dT_fortran=True) #TODO: dT_fortran since connection was forcibly closed, revert. #,tzone='UTC+01:00') #timezone rotates entire aardappelgrafiek, so decides what is neap and springtide
if str(data_pd_moonculm.loc[0,'datetime'].tz) != 'UTC': # important since data_pd_HWLW['culm_hr']=range(12) hourvalues should be in UTC since that relates to the relation dateline/sun
    raise Exception(f'culmination data is not in expected timezone (UTC): {data_pd_moonculm.loc[0,"datetime"].tz}')
data_pd_moonculm['datetime'] = data_pd_moonculm['datetime'].dt.tz_localize(None)
data_pd_moonculm = data_pd_moonculm.set_index('datetime',drop=False)
data_pd_moonculm['values'] = data_pd_moonculm['type'] #dummy values for TA in hatyan.calc_HWLWnumbering()
data_pd_moonculm['HWLWcode'] = 1 #all HW values since one every ~12h25m
data_pd_moonculm = hatyan.calc_HWLWnumbering(data_pd_moonculm,doHWLWcheck=False) #TODO: currently w.r.t. cadzd, is that an issue? With DELFZL the matched culmination is incorrect (since far away), but that might not be a big issue
data_pd_moonculm['HWLWno_offset'] = data_pd_moonculm['HWLWno']+4 #correlate HWLW to moonculmination 2 days before. TODO: check this offset in relation to culm_addtime.
moonculm_idxHWLWno = data_pd_moonculm.set_index('HWLWno_offset')

for current_station in ['STELLDBTN']:# stat_list:#['HOEKVHLD']:#['CADZD','VLISSGN','HARVT10','HOEKVHLD','IJMDBTHVN','DENOVBTN','KATSBTN','KORNWDZBTN','OUDSD','SCHEVNGN']:#stat_list:
    if closefigatstart:
        plt.close('all')
    print(f'havengetallen for {current_station}')
    
    #read HWLW data
    file_ext_pkl = os.path.join(dir_meas,f"{current_station}_measext.pkl")
    if not os.path.exists(file_ext_pkl):
        continue
    data_pd_HWLW = pd.read_pickle(file_ext_pkl)
    data_pd_HWLW = clean_data(data_pd_HWLW)
    if NAP2005correction: # apply NAP correction #TODO: move to clean_data?
        data_pd_HWLW = nap2005_correction(data_pd_HWLW,current_station=current_station)
    #crop timeseries
    data_pd_HWLW = hatyan.crop_timeseries(data_pd_HWLW, times_ext=[tstart_dt,tstop_dt],onlyfull=False)
    
    #check if amount of HWs is enough
    numHWs_expected = (tstop_dt-tstart_dt).total_seconds()/M2_period_timedelta.total_seconds()
    numHWs = (data_pd_HWLW['HWLWcode']==1).sum()
    if numHWs < 0.95*numHWs_expected:
        raise Exception(f'ERROR: not enough high waters present in period, {numHWs} instead of >=0.95*{int(numHWs_expected):d}')
    
    print('SELECT/CALC HWLW VALUES')
    if len(data_pd_HWLW['HWLWcode'].unique()) > 2:
        data_pd_HWLW = hatyan.calc_HWLW12345to12(data_pd_HWLW) #convert 12345 to 12 by taking minimum of 345 as 2 (laagste laagwater)
    
    if current_station in ['KATSBTN','GATVBSLE','HANSWT']:
        #TODO: move this to data check part?
        #TODO: this removes extreme values that are 1/5/28 min from each other, but they should not be present to begin with
        timediff = data_pd_HWLW.index[1:]-data_pd_HWLW.index[:-1]
        data_pd_HWLW['timediff'] = pd.TimedeltaIndex([pd.NaT]).append(timediff)
        bool_tooclose = data_pd_HWLW['timediff']<dt.timedelta(minutes=30)
        print('unique small timestep_min:',data_pd_HWLW.loc[bool_tooclose,'timediff'].unique()/1e9/60)
        data_pd_HWLW = data_pd_HWLW.loc[~bool_tooclose]
    if current_station in ['STELLDBTN']: #TODO: manual removal of invalid HW value from STELLDBTN (is flat line in wl timeseries)
        """
        file_wl_pkl = os.path.join(dir_meas,f"{current_station}_measwl.pkl")
        data_pd_wl_all = pd.read_pickle(file_wl_pkl)        
        data_pd_wl_all.index = data_pd_wl_all.index.tz_localize(None)
        data_pd_wl_all = hatyan.crop_timeseries(data_pd_wl_all, times_ext=[tstart_dt,tstop_dt],onlyfull=False)
        data_pd_wl_all_ext = hatyan.calc_HWLW(data_pd_wl_all)
        fig,(ax1,ax2) = hatyan.plot_timeseries(ts=data_pd_wl_all, ts_ext=data_pd_HWLW)
        fig,(ax1,ax2) = hatyan.plot_timeseries(ts=data_pd_wl_all, ts_ext=data_pd_wl_all_ext)
        """
        drop_time_STELLDBTN = '2012-02-09 09:36:00'
        if drop_time_STELLDBTN in data_pd_HWLW.index:
            data_pd_HWLW = data_pd_HWLW.drop(drop_time_STELLDBTN)
    
    data_pd_HWLW_idxHWLWno = hatyan.calc_HWLWnumbering(data_pd_HWLW)
    data_pd_HWLW_idxHWLWno['times'] = data_pd_HWLW_idxHWLWno.index
    data_pd_HWLW_idxHWLWno = data_pd_HWLW_idxHWLWno.set_index('HWLWno',drop=False)
    
    HW_bool = data_pd_HWLW_idxHWLWno['HWLWcode']==1
    data_pd_HWLW_idxHWLWno.loc[HW_bool,'getijperiod'] = data_pd_HWLW_idxHWLWno.loc[HW_bool,'times'].iloc[1:].values - data_pd_HWLW_idxHWLWno.loc[HW_bool,'times'].iloc[:-1] #this works properly since index is HWLW
    data_pd_HWLW_idxHWLWno.loc[HW_bool,'duurdaling'] = data_pd_HWLW_idxHWLWno.loc[~HW_bool,'times'] - data_pd_HWLW_idxHWLWno.loc[HW_bool,'times']
    data_pd_HWLW_idxHWLWno['culm_time'] = moonculm_idxHWLWno['datetime'] #couple HWLW to moonculminations two days earlier (this works since index is HWLWno)
    data_pd_HWLW_idxHWLWno['culm_hr'] = (data_pd_HWLW_idxHWLWno['culm_time'].round('h').dt.hour)%12
    data_pd_HWLW_idxHWLWno['HWLW_delay'] = (data_pd_HWLW_idxHWLWno['times']-(data_pd_HWLW_idxHWLWno['culm_time']+culm_addtime))
    data_pd_HWLW = data_pd_HWLW_idxHWLWno.set_index('times')
    
    print('calculate medians per hour group for LW and HW (instead of 1991 method: average of subgroups with removal of outliers)')
    data_pd_HW = data_pd_HWLW.loc[data_pd_HWLW['HWLWcode']==1]
    data_pd_LW = data_pd_HWLW.loc[data_pd_HWLW['HWLWcode']==2]
    HWLW_culmhr_summary = pd.DataFrame()
    HWLW_culmhr_summary['HW_values_median'] = data_pd_HW.groupby(data_pd_HW['culm_hr'])['values'].median()
    HWLW_culmhr_summary['HW_delay_median'] = data_pd_HW.groupby(data_pd_HW['culm_hr'])['HWLW_delay'].median()
    HWLW_culmhr_summary['LW_values_median'] = data_pd_LW.groupby(data_pd_LW['culm_hr'])['values'].median()
    HWLW_culmhr_summary['LW_delay_median'] = data_pd_LW.groupby(data_pd_LW['culm_hr'])['HWLW_delay'].median()
    HWLW_culmhr_summary['getijperiod_median'] = data_pd_HW.groupby(data_pd_HW['culm_hr'])['getijperiod'].median()
    HWLW_culmhr_summary['duurdaling_median'] = data_pd_HW.groupby(data_pd_HW['culm_hr'])['duurdaling'].median()
        
    print('HWLW FIGUREN PER TIJDSKLASSE, INCLUSIEF MEDIAN LINE')
    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(18,8), sharex=True)
    ax1.set_title('HW values %s'%(current_station))
    ax1.plot(data_pd_HW['culm_hr'],data_pd_HW['values'],'.')
    ax1.plot(HWLW_culmhr_summary['HW_values_median'],'.-')
    ax2.set_title('LW values %s'%(current_station))
    ax2.plot(data_pd_LW['culm_hr'],data_pd_LW['values'],'.')
    ax2.plot(HWLW_culmhr_summary['LW_values_median'],'.-')
    ax3.set_title('HW time delays %s'%(current_station))
    ax3.plot(data_pd_HW['culm_hr'],data_pd_HW['HWLW_delay'].dt.total_seconds()/3600,'.')
    ax3.plot(HWLW_culmhr_summary['HW_delay_median'].dt.total_seconds()/3600,'.-')
    ax4.set_title('LW time delays %s'%(current_station))
    ax4.plot(data_pd_LW['culm_hr'],data_pd_LW['HWLW_delay'].dt.total_seconds()/3600,'.')
    ax4.plot(HWLW_culmhr_summary['LW_delay_median'].dt.total_seconds()/3600,'.-')
    ax4.set_xlim([0-0.5,12-0.5])
    fig.tight_layout()
    fig.savefig(os.path.join(dir_havget,f'HWLW_pertijdsklasse_inclmedianline_{current_station}'))
    
    file_outname = os.path.join(dir_havget, f'aardappelgrafiek_{year_slotgem}_{current_station}')
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
    ax1.plot(HWLW_culmhr_summary['HW_delay_median'].dt.total_seconds()/3600,HWLW_culmhr_summary['HW_values_median'],'.-',label=current_station)
    ax1.xaxis.set_major_formatter(timeTicks)
    ax1.grid()
    ax2.set_title(f'LW {current_station} {year_slotgem}')
    ax2.set_xlabel('maansverloop in uu:mm:ss' )
    ax2.set_ylabel('waterstand in m t.o.v. NAP')
    ax2.plot(HWLW_culmhr_summary['LW_delay_median'].dt.total_seconds()/3600,HWLW_culmhr_summary['LW_values_median'],'.-',label=current_station)
    ax2.xaxis.set_major_formatter(timeTicks)
    ax2.grid()
    for iH,row in HWLW_culmhr_summary.iterrows():
        ax1.text(row['HW_delay_median'].total_seconds()/3600,row['HW_values_median'], str(int(iH)))
        ax2.text(row['LW_delay_median'].total_seconds()/3600,row['LW_values_median'], str(int(iH)))
    #set equal ylims
    ax1_xlimmean = np.mean(ax1.get_xlim())
    ax2_xlimmean = np.mean(ax2.get_xlim())
    ax1_ylimmean = np.mean(ax1.get_ylim())
    ax2_ylimmean = np.mean(ax2.get_ylim())
    xlimrange = 2
    ylimrange = 1
    ax1.set_xlim([ax1_xlimmean-xlimrange/2,ax1_xlimmean+xlimrange/2])
    ax2.set_xlim([ax2_xlimmean-xlimrange/2,ax2_xlimmean+xlimrange/2])
    ax1.set_ylim([ax1_ylimmean-ylimrange/2,ax1_ylimmean+ylimrange/2])
    ax2.set_ylim([ax2_ylimmean-ylimrange/2,ax2_ylimmean+ylimrange/2])
    #plot gemtij dotted lines
    ax1.plot(ax1.get_xlim(),[HWLW_culmhr_summary['HW_values_median'].mean(),HWLW_culmhr_summary['HW_values_median'].mean()],'k--')
    ax1.plot([HWLW_culmhr_summary['HW_delay_median'].mean().total_seconds()/3600,HWLW_culmhr_summary['HW_delay_median'].mean().total_seconds()/3600],ax1.get_ylim(),'k--')
    ax2.plot(ax2.get_xlim(),[HWLW_culmhr_summary['LW_values_median'].mean(),HWLW_culmhr_summary['LW_values_median'].mean()],'k--')
    ax2.plot([HWLW_culmhr_summary['LW_delay_median'].mean().total_seconds()/3600,HWLW_culmhr_summary['LW_delay_median'].mean().total_seconds()/3600],ax2.get_ylim(),'k--')
    fig.tight_layout()
    fig.savefig(file_outname)
    
    #write to csv
    HWLW_culmhr_summary_out = HWLW_culmhr_summary.copy()
    HWLW_culmhr_summary_out.loc['mean',:] = HWLW_culmhr_summary_out.mean() #add mean row to dataframe (not convenient to add immediately due to plotting with index 0-11)
    for colname in HWLW_culmhr_summary_out.columns: #round timedelta to make outputformat nicer
        if HWLW_culmhr_summary_out[colname].dtype == 'timedelta64[ns]':
            HWLW_culmhr_summary_out[colname] = HWLW_culmhr_summary_out[colname].round('S')
    HWLW_culmhr_summary_out.to_csv(file_outname+'.csv',float_format='%.3f')
    



#TODO IMPORTANT: uncertainty about length of analysis period (and SA/SM origin)
#TODO IMPORTANT: correct havengetallen with slotgemiddelden before using them for gemiddelde getijkromme
#TODO IMPORTANT: scaling is now max 18.2% but this is quite a lot, check values for all stations?
##### gemiddelde getijkrommen
for current_station in stat_list:#['SCHEVNGN']:#stat_list[stat_list.index('SCHEVNGN'):]:#stat_list:#['HOEKVHLD']:#['HOEKVHLD','HARVT10']: stat_list[stat_list.index('SCHEVNGN'):]
    """
    
    """
    if closefigatstart:
        plt.close('all')
    print(f'gem getijkrommen for {current_station}')
    
    dir_vali_krommen = r'p:\archivedprojects\11205258-005-kpp2020_rmm-g5\C_Work\00_KenmerkendeWaarden\07_Figuren\figures_ppSCL_2\final20201211'
    file_vali_doodtijkromme = os.path.join(dir_vali_krommen,f'doodtijkromme_{current_station}_havengetallen{year_slotgem}.csv')
    file_vali_gemtijkromme = os.path.join(dir_vali_krommen,f'gemGetijkromme_{current_station}_havengetallen{year_slotgem}.csv')
    file_vali_springtijkromme = os.path.join(dir_vali_krommen,f'springtijkromme_{current_station}_havengetallen{year_slotgem}.csv')        
    
    #TODO: add correctie havengetallen HW/LW av/sp/np met slotgemiddelde uit PLSS/modelfit (HW/LW av)
    file_havget = os.path.join(dir_havget,f'aardappelgrafiek_{year_slotgem}_{current_station}.csv')
    if not os.path.exists(file_havget):
        raise Exception(f'havengetallen file does not exist: {file_havget}')
    data_havget = pd.read_csv(file_havget)
    for colname in ['HW_delay_median','LW_delay_median','getijperiod_median','duurdaling_median']:
        data_havget[colname] = data_havget[colname].apply(lambda x: pd.Timedelta(x))
    HW_sp, LW_sp, tD_sp = data_havget.loc[0,['HW_values_median','LW_values_median','duurdaling_median']]
    HW_np, LW_np, tD_np = data_havget.loc[6,['HW_values_median','LW_values_median','duurdaling_median']]
    HW_av, LW_av, tD_av = data_havget.loc[12,['HW_values_median','LW_values_median','duurdaling_median']]
        
    
    def reshape_signal(ts, ts_ext, HW_goal, LW_goal, tP_goal=None):
        """
        scales tidal signal to provided HW/LW value and up/down going time
        tP_goal (tidal period time) is used to fix tidalperiod to 12h25m (for BOI timeseries)
        
        time_down was scaled with havengetallen before, but not anymore to avoid issues with aggers
        """
        TR_goal = HW_goal-LW_goal
        
        #selecteer alle hoogwaters en opvolgende laagwaters
        bool_HW = ts_ext['HWLWcode']==1
        idx_HW = np.where(bool_HW)[0]
        timesHW = ts_ext.index[idx_HW]
        timesLW = ts_ext.index[idx_HW[:-1]+1] #assuming alternating 1,2,1 or 1,3,1, this is always valid in this workflow
        
        #crop from first to last HW (rest is not scaled anyway)
        ts_time_firstHW = ts_ext[bool_HW].index[0]
        ts_time_lastHW = ts_ext[bool_HW].index[-1]
        ts_corr = ts.copy().loc[ts_time_firstHW:ts_time_lastHW]

        ts_corr['times'] = ts_corr.index #this is necessary since datetimeindex with freq is not editable, and Series is editable
        ts_corr['values_new'] = np.nan #necessary since HW is read and overwitten twice
        for i in np.arange(0,len(timesHW)-1):
            HW1_val = ts_corr.loc[timesHW[i],'values']
            HW2_val = ts_corr.loc[timesHW[i+1],'values']
            LW_val = ts_corr.loc[timesLW[i],'values']
            TR1_val = HW1_val-LW_val
            TR2_val = HW2_val-LW_val
            tP_val = timesHW[i+1]-timesHW[i]
            if tP_goal is None:
                tP_goal = tP_val
            tD_val = timesLW[i]-timesHW[i]
            tD_goal = tD_val/tP_val*tP_goal #no change if tP_goal is None
            tU_goal = tP_goal-tD_val #equal to tP_val-tD_goal if tP_goal is None
            
            print(f'tidalrange factor: {TR_goal/TR1_val:.3f}')
            #print(f'timeDown factor: {tD_goal/tD_val:.3f}')
            factors = np.array([TR_goal/TR1_val,tD_goal/tD_val])
            allowed_perc = 18.2 #TODO: 14.6 necesary for 2021 DOODTIJ tidalrange factor BAALHK, 15.1 for 2021 DOODTIJ tidalrange factor BATH, 17.0 for 2021 DOODTIJ tidalrange factor HARLGN, 18.2 for 2021 DOODTIJ tidalrange factor STELLDBTN
            if (factors>(1+allowed_perc/100)).any() or (factors<(1-allowed_perc/100)).any():
                raise Exception(f'more than {allowed_perc}% decrease or increase')
            
            tide_HWtoLW = ts_corr.loc[timesHW[i]:timesLW[i]]
            tide_LWtoHW = ts_corr.loc[timesLW[i]:timesHW[i+1]]
            
            ts_corr.loc[timesHW[i]:timesLW[i],'times'] = pd.date_range(start=ts_corr.loc[timesHW[i],'times'],end=ts_corr.loc[timesHW[i],'times']+tD_goal,periods=len(tide_HWtoLW))
            ts_corr.loc[timesHW[i]:timesLW[i],'values_new'] = (ts_corr.loc[timesHW[i]:timesLW[i],'values']-LW_val)/TR1_val*TR_goal+LW_goal
            ts_corr.loc[timesLW[i]:timesHW[i+1],'times'] = pd.date_range(start=ts_corr.loc[timesLW[i],'times'],end=ts_corr.loc[timesLW[i],'times']+tU_goal,periods=len(tide_LWtoHW))
            ts_corr.loc[timesLW[i]:timesHW[i+1],'values_new'] = (ts_corr.loc[timesLW[i]:timesHW[i+1],'values']-LW_val)/TR2_val*TR_goal+LW_goal
        ts_corr = ts_corr.set_index('times',drop=True)
        ts_corr['values'] = ts_corr['values_new']
        ts_corr = ts_corr.drop(['values_new'],axis=1)
        return ts_corr

    def ts_to_trefHW(ts,HWreftime):
        """
        converts to hours relative to HWreftime, to plot av/sp/np tidal signals in one plot
        """
        ts.index.name = 'times' #just to be sure
        ts_trefHW = ts.reset_index()
        ts_trefHW.index = (ts_trefHW['times']-HWreftime).dt.total_seconds()/3600
        return ts_trefHW
    
    def repeat_signal(ts_one_HWtoHW, nb, na):
        """
        repeat tidal signal, necessary for sp/np, since they are computed as single tidal signal first
        """
        tidalperiod = ts_one_HWtoHW.index[-1] - ts_one_HWtoHW.index[0]
        ts_rep = pd.DataFrame()
        for iAdd in np.arange(-nb,na+1):
            ts_add = pd.DataFrame({'values':ts_one_HWtoHW['values'].values},
                                  index=ts_one_HWtoHW.index + iAdd*tidalperiod)
            ts_rep = pd.concat([ts_rep,ts_add])
        ts_rep = ts_rep.loc[~ts_rep.index.duplicated()]
        return ts_rep
    
    #load measurement data
    file_wl_pkl = os.path.join(dir_meas,f"{current_station}_measwl.pkl")
    ts_meas_pd = pd.read_pickle(file_wl_pkl)
    ts_meas_pd = clean_data(ts_meas_pd)
    ts_meas_pd = hatyan.crop_timeseries(ts_meas_pd, times_ext=[tstart_dt,tstop_dt-dt.timedelta(minutes=10)])#,onlyfull=False)
    if NAP2005correction:
        ts_meas_pd = nap2005_correction(ts_meas_pd,current_station)
    
    # =============================================================================
    # Hatyan analyse voor 10 jaar (alle componenten voor gemiddelde getijcyclus) #TODO: maybe use original 4y period instead? SA/SM should come from 19y analysis
    # =============================================================================
    const_list = hatyan.get_const_list_hatyan('year') #this should not be changed, since higher harmonics are necessary
    hatyan_settings_ana = hatyan.HatyanSettings(nodalfactors=True,
                                                fu_alltimes=False, # False is RWS-default
                                                xfac=True, # True is RWS-default
                                                analysis_perperiod='Y',
                                                xTxmat_condition_max=15, #TODO: for some reason this is necessary for HOEKVHLD 2006 (default=10)
                                                return_allperiods=True)
    comp_frommeasurements_avg, comp_frommeasurements_allyears = hatyan.get_components_from_ts(ts_meas_pd, const_list=const_list, hatyan_settings=hatyan_settings_ana)
    
    #check if all years are available
    comp_years = comp_frommeasurements_allyears['A'].columns
    expected_years = tstop_dt.year-tstart_dt.year
    if len(comp_years) < expected_years:
        raise Exception('ERROR: analysis result contains not all years')
        
    # =============================================================================
    # gemiddelde getijkromme
    # =============================================================================
    """
    uit: gemiddelde getijkrommen 1991.0
    Voor meetpunten in het onbeinvloed gebied is per getijfase eerst een "ruwe kromme" berekend met de resultaten van de harmonische analyse, 
    welke daarna een weinig is bijgesteld aan de hand van de volgende slotgemiddelden:
    gemiddeld hoog- en laagwater, duur daling. Deze bijstelling bestaat uit een eenvoudige vermenigvuldiging.    
    
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
    dit in tegenstelling tot vroegere bepalingen. Bij spring- en doodtij is bovendien de differentiele getijduur, 
    en daarmee de duur rijzing, afgeleid uit de ruwe krommen.
    
    """
    #kwadraatsommen voor M2 tot M12
    components_av = ['M2','M4','M6','M8','M10','M12']
    comp_av = comp_frommeasurements_avg.loc[components_av]
    for comp_higherharmonics in components_av:
        iM = int(comp_higherharmonics[1:])
        bool_endswithiM = comp_frommeasurements_avg.index.str.endswith(str(iM)) & comp_frommeasurements_avg.index.str.replace(str(iM),'').str[-1].str.isalpha()
        comp_iM = comp_frommeasurements_avg.loc[bool_endswithiM]
        #print(comp_iM)
        comp_av.loc[comp_higherharmonics,'A'] = np.sqrt((comp_iM['A']**2).sum()) #kwadraatsom
    
    print(f'verhouding tussen originele en kwadratensom componenten: {current_station}')
    print(comp_av/comp_frommeasurements_avg.loc[components_av]) #TODO: values are different than 1991.0 document, but could be because of different year so check with 1981-1991 data. Statement "Zoals te verwachten is de verhouding per component tussen deze wortel en de oorspronkelijke amplitude voor alle plaatsen gelijk." seems to be not true. Could also differ because of 10 instead of 4 analysis years?
    
    comp_av.loc['A0'] = comp_frommeasurements_avg.loc['A0']
    freq_sec = 10 #TODO: frequency decides accuracy of tU/tD and other timings (and is writing freq of BOI timeseries)
    times_pred_1mnth = pd.date_range(start=dt.datetime(tstop_dt.year, 1, 1, 0, 0)-dt.timedelta(hours=12), end=dt.datetime(tstop_dt.year, 2, 1, 0, 0), freq=f'{freq_sec} S') #start 12 hours in advance, to assure also corrected values on desired tstart
    prediction_av = hatyan.prediction(comp_av, times_pred_all=times_pred_1mnth, nodalfactors=False) #nodalfactors=False to guarantee repetative signal
    prediction_av_ext = hatyan.calc_HWLW(ts=prediction_av, calc_HWLW345=False)
        
    time_firstHW = prediction_av_ext.loc[prediction_av_ext['HWLWcode']==1].index[0] #time of first HW
    ia1 = prediction_av_ext.loc[time_firstHW:].index[0] #time of first HW
    ia2 = prediction_av_ext.loc[time_firstHW:].index[2] #time of second HW
    prediction_av_one = prediction_av.loc[ia1:ia2]
    prediction_av_ext_one = prediction_av_ext.loc[ia1:ia2]
    
    
    # =============================================================================
    # Hatyan predictie voor 1 jaar met gemiddelde helling maansbaan (voor afleiden spring-doodtijcyclus) >> predictie zonder nodalfactors instead
    # =============================================================================
    """
    uit: gemiddelde getijkrommen 1991.0
    Voor de ruwe krommen voor springtij en doodtij is het getij voorspeld voor een jaar met gemiddelde helling maansbaan 
    met uitsluitend zuivere combinaties van de componenten M2 en S2:
    tabel: Gebruikte componenten voor de spring- en doodtijkromme
    SM, 3MS2, mu2, M2, S2, 2SM2, 3MS4, M4, MS4, 
    4MS6, M6, 2MS6, M8, 3MS8, M10, 4MS10, M12, 5MS12
    
    In het aldus gemodelleerde getij is de vorm van iedere getijslag, gegeven de getijfase, identiek. 
    Vervolgens is aan de hand van de havengetallen een springtij- en een doodtijkromme geselecteerd.
    """
    
    """
    #TODO, below is different than provided list, these shallow ones are extra: ['S4','2SM6','M7','4MS4','2(MS)8','3M2S10','4M2S12']
    #shallow relations, derive 'zuivere harmonischen van M2 en S2' (this means averaging over eenmaaldaagse componenten, but why is that chosen?)
    #adding above extra components or oneday freqs, gives a modulation and therefore there is no repetative signal anymore. Apperantly all components in this list have an integer number of periods in one springneap cycle?
    dummy,shallowrel,dummy = hatyan.get_foreman_shallowrelations()
    bool_M2S2only = shallowrel[1].isin([1,2]) & shallowrel[3].isin(['M2','S2']) & shallowrel[5].isin(['M2','S2',np.nan]) & shallowrel.index.isin(const_list_year)
    shallowdeps_M2S2 = shallowrel.loc[bool_M2S2only,:5]
    print(shallowdeps_M2S2)
    """
    components_sn = ['A0','SM','3MS2','MU2','M2','S2','2SM2','3MS4','M4','MS4','4MS6','M6','2MS6','M8','3MS8','M10','4MS10','M12','5MS12']
    
    #make prediction with springneap components with nodalfactors=False (alternative for choosing a year with a neutral nodal factor). Using 1yr instead of 1month does not make a difference in min/max tidal range and shape, also because of nodalfactors=False. (when using more components, there is a slight difference)
    comp_frommeasurements_avg_sncomp = comp_frommeasurements_avg.loc[components_sn]
    prediction_sn = hatyan.prediction(comp_frommeasurements_avg_sncomp, times_pred_all=times_pred_1mnth, nodalfactors=False) #nodalfactors=False to make independent on chosen year
    
    prediction_sn_ext = hatyan.calc_HWLW(ts=prediction_sn, calc_HWLW345=False)
    
    #selecteer getijslag met minimale tidalrange en maximale tidalrange (werd geselecteerd adhv havengetallen in 1991.0 doc)
    #TODO: hatyan.calc_HWLW_tidalrange() creëeren (code staat hieronder), geeft column 'tidalrange' in ts_ext dataframe
    prediction_sn_ext = hatyan.calc_HWLWnumbering(ts_ext=prediction_sn_ext)
    prediction_sn_ext['times_backup'] = prediction_sn_ext.index
    prediction_sn_ext_idxHWLWno = prediction_sn_ext.set_index('HWLWno',drop=False)
    prediction_sn_ext_idxHWLWno['tidalrange'] = prediction_sn_ext_idxHWLWno.loc[prediction_sn_ext_idxHWLWno['HWLWcode']==1,'values'] - prediction_sn_ext_idxHWLWno.loc[prediction_sn_ext_idxHWLWno['HWLWcode']==2,'values']
    prediction_sn_ext = prediction_sn_ext_idxHWLWno.set_index('times_backup')
    
    time_TRmax = prediction_sn_ext.loc[prediction_sn_ext['HWLWcode']==1,'tidalrange'].idxmax()
    is1 = prediction_sn_ext.loc[time_TRmax:].index[0]
    is2 = prediction_sn_ext.loc[time_TRmax:].index[2]
    
    time_TRmin = prediction_sn_ext.loc[prediction_sn_ext['HWLWcode']==1,'tidalrange'].idxmin()
    in1 = prediction_sn_ext.loc[time_TRmin:].index[0]
    in2 = prediction_sn_ext.loc[time_TRmin:].index[2]
    
    #select one tideperiod for springtide and one for neaptide
    prediction_sp_one = prediction_sn.loc[is1:is2]
    prediction_sp_ext_one = prediction_sn_ext.loc[is1:is2]
    prediction_np_one = prediction_sn.loc[in1:in2]
    prediction_np_ext_one = prediction_sn_ext.loc[in1:in2]
    
    # plot selection of neap/spring
    fig, (ax1,ax2) = hatyan.plot_timeseries(ts=prediction_sn,ts_ext=prediction_sn_ext)
    ax1.plot(prediction_sp_one['values'],'r')
    ax1.plot(prediction_np_one['values'],'r')
    ax1.legend(labels=ax1.get_legend_handles_labels()[1]+['kromme spring','kromme neap'],loc=4)
    ax1.set_ylabel('waterstand [m]')
    ax1.set_title(f'spring- en doodtijkromme {current_station}')
    fig.savefig(os.path.join(dir_gemgetij,f'springdoodtijkromme_{current_station}_slotgem{year_slotgem}.png'))
    
    
    print(f'reshape_signal GEMGETIJ: {current_station}')
    prediction_av_one_trefHW = ts_to_trefHW(prediction_av_one,HWreftime=ia1) # repeating one is not necessary for av, but easier to do the same for av/sp/np
    prediction_av_corr_one = reshape_signal(prediction_av_one, prediction_av_ext_one, HW_goal=HW_av, LW_goal=LW_av, tP_goal=None)
    prediction_av_corr_rep5 = repeat_signal(prediction_av_corr_one, nb=2, na=2)
    prediction_av_corr_rep5_trefHW = ts_to_trefHW(prediction_av_corr_rep5,HWreftime=ia1)

    print(f'reshape_signal SPRINGTIJ: {current_station}')
    prediction_sp_one_trefHW = ts_to_trefHW(prediction_sp_one,HWreftime=is1)
    prediction_sp_corr_one = reshape_signal(prediction_sp_one, prediction_sp_ext_one, HW_goal=HW_sp, LW_goal=LW_sp, tP_goal=None)
    prediction_sp_corr_rep5 = repeat_signal(prediction_sp_corr_one, nb=2, na=2)
    prediction_sp_corr_rep5_trefHW = ts_to_trefHW(prediction_sp_corr_rep5,HWreftime=is1)
    
    print(f'reshape_signal DOODTIJ: {current_station}')
    prediction_np_one_trefHW = ts_to_trefHW(prediction_np_one,HWreftime=in1)
    prediction_np_corr_one = reshape_signal(prediction_np_one, prediction_np_ext_one, HW_goal=HW_np, LW_goal=LW_np, tP_goal=None)
    prediction_np_corr_rep5 = repeat_signal(prediction_np_corr_one, nb=2, na=2)
    prediction_np_corr_rep5_trefHW = ts_to_trefHW(prediction_np_corr_rep5,HWreftime=in1)
    
    
    #12u25m timeseries for BOI computations (no relation between HW and moon, HW has to come at same time for av/sp/np tide, HW timing does differ between stations)
    print(f'reshape_signal BOI GEMGETIJ and write to csv: {current_station}')
    prediction_av_corrBOI_one = reshape_signal(prediction_av_one, prediction_av_ext_one, HW_goal=HW_av, LW_goal=LW_av, tP_goal=pd.Timedelta(hours=12,minutes=25))
    prediction_av_corrBOI_one_roundtime = prediction_av_corrBOI_one.resample(f'{freq_sec}S').nearest()
    prediction_av_corrBOI_one_roundtime.to_csv(os.path.join(dir_gemgetij,f'gemGetijkromme_BOI_{current_station}_slotgem{year_slotgem}.csv'),float_format='%.3f',date_format='%Y-%m-%d %H:%M:%S')
    prediction_av_corrBOI_repn_roundtime = repeat_signal(prediction_av_corrBOI_one_roundtime, nb=0, na=10)
    
    print(f'reshape_signal BOI SPRINGTIJ and write to csv: {current_station}')
    prediction_sp_corrBOI_one = reshape_signal(prediction_sp_one, prediction_sp_ext_one, HW_goal=HW_sp, LW_goal=LW_sp, tP_goal=pd.Timedelta(hours=12,minutes=25))
    prediction_sp_corrBOI_one.index = prediction_sp_corrBOI_one.index - prediction_sp_corrBOI_one.index[0] + prediction_av_corrBOI_one.index[0] #shift times to first HW from gemgetij
    prediction_sp_corrBOI_one_roundtime = prediction_sp_corrBOI_one.resample(f'{freq_sec}S').nearest()
    prediction_sp_corrBOI_one_roundtime.to_csv(os.path.join(dir_gemgetij,f'springtijkromme_BOI_{current_station}_slotgem{year_slotgem}.csv'),float_format='%.3f',date_format='%Y-%m-%d %H:%M:%S')
    prediction_sp_corrBOI_repn_roundtime = repeat_signal(prediction_sp_corrBOI_one_roundtime, nb=0, na=10)

    print(f'reshape_signal BOI DOODTIJ and write to csv: {current_station}')
    prediction_np_corrBOI_one = reshape_signal(prediction_np_one, prediction_np_ext_one, HW_goal=HW_np, LW_goal=LW_np, tP_goal=pd.Timedelta(hours=12,minutes=25))
    prediction_np_corrBOI_one.index = prediction_np_corrBOI_one.index - prediction_np_corrBOI_one.index[0] + prediction_av_corrBOI_one.index[0] #shift times to first HW from gemgetij
    prediction_np_corrBOI_one_roundtime = prediction_np_corrBOI_one.resample(f'{freq_sec}S').nearest()
    prediction_np_corrBOI_one_roundtime.to_csv(os.path.join(dir_gemgetij,f'doodtijkromme_BOI_{current_station}_slotgem{year_slotgem}.csv'),float_format='%.3f',date_format='%Y-%m-%d %H:%M:%S')
    prediction_np_corrBOI_repn_roundtime = repeat_signal(prediction_np_corrBOI_one_roundtime, nb=0, na=10)
    
    
    cmap = plt.get_cmap("tab10")
        
    print(f'plot getijkromme trefHW: {current_station}')
    fig_sum,ax_sum = plt.subplots(figsize=(14,7))
    ax_sum.set_title(f'getijkromme trefHW {current_station}')
    ax_sum.plot(prediction_av_one_trefHW['values'],'--', color=cmap(0),linewidth=0.7, label='gem kromme, one')
    ax_sum.plot(prediction_av_corr_rep5_trefHW['values'], color=cmap(0), label='gem kromme, corr')
    ax_sum.plot(prediction_sp_one_trefHW['values'],'--', color=cmap(1),linewidth=0.7, label='sp kromme, one')
    ax_sum.plot(prediction_sp_corr_rep5_trefHW['values'], color=cmap(1), label='sp kromme, corr')
    ax_sum.plot(prediction_np_one_trefHW['values'],'--', color=cmap(2),linewidth=0.7, label='np kromme, one')
    ax_sum.plot(prediction_np_corr_rep5_trefHW['values'], color=cmap(2), label='np kromme, corr')
    ax_sum.legend(loc=4)
    ax_sum.grid()
    ax_sum.set_xlim(-15.5,15.5)
    ax_sum.set_xlabel('hours since HW (ts are shifted to this reference)')
    fig_sum.tight_layout()
    fig_sum.savefig(os.path.join(dir_gemgetij,f'gemgetij_trefHW_{current_station}'))
    
    print(f'plot BOI figure and compare to KW2020: {current_station}')
    fig_boi,ax1_boi = plt.subplots(figsize=(14,7))
    ax1_boi.set_title(f'getijkromme BOI {current_station}')
    #gemtij
    ax1_boi.plot(prediction_av_corrBOI_repn_roundtime['values'],color=cmap(0),label='prediction gemtij')
    if 0:#os.path.exists(file_vali_gemtijkromme):#TODO: for some reason, spyder/explorer/tcmd freezes when interacting with these files/folders
        data_vali_gemtij = pd.read_csv(file_vali_gemtijkromme,index_col=0,parse_dates=True)
        ax1_boi.plot(data_vali_gemtij['Water Level [m]'],'--',color=cmap(0),linewidth=0.7,label='validation KW2020 gemtij')
    #springtij
    ax1_boi.plot(prediction_sp_corrBOI_repn_roundtime['values'],color=cmap(1),label='prediction springtij')
    if 0:#os.path.exists(file_vali_springtijkromme):
        data_vali_springtij = pd.read_csv(file_vali_springtijkromme,index_col=0,parse_dates=True)
        ax1_boi.plot(data_vali_springtij['Water Level [m]'],'--',color=cmap(1),linewidth=0.7,label='validation KW2020 springtij')
    #doodtij
    ax1_boi.plot(prediction_np_corrBOI_repn_roundtime['values'],color=cmap(2),label='prediction doodtij')
    if 0:#os.path.exists(file_vali_doodtijkromme):
        data_vali_doodtij = pd.read_csv(file_vali_doodtijkromme,index_col=0,parse_dates=True)
        ax1_boi.plot(data_vali_doodtij['Water Level [m]'],'--',color=cmap(2),linewidth=0.7, label='validation KW2020 doodtij')
    ax1_boi.grid()
    ax1_boi.legend(loc=4)
    ax1_boi.set_xlabel('times since first av HW (start of ts)')
    ax1_boi.set_xlim(tstop_dt-dt.timedelta(hours=2),tstop_dt+dt.timedelta(hours=48))
    fig_boi.tight_layout()
    fig_boi.savefig(os.path.join(dir_gemgetij,f'gemspringdoodtijkromme_BOI_{current_station}_slotgem{year_slotgem}.png'))
    
    







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

#TODO: overschrijdingsfreqs met extremen: neem datareeks vanaf waar geen grote gaps meer voorkomen, of maakt dat niet uit?


"""
#TODO: aantekeningen gesprek Boyan
○ Je vertaalt niet x aantal datapunten naar frequentie, maar je zet de punten op volgorde en je rankt ze, daarvan maak je distributie, ranking en frequentie is niet 1 op 1
○ Max freq is 2 getij per dag, keer 365 dagen, maximale frequentie komt daarmee overeen. (on)gefilterd en trendanalys is datapunten op volgorde en frequentie, 
○ Lezen:
    o rapport boyan kw-rmm: n:\\Projects\\11205000\11205232\\C. Report - advise\\007 - Kenmerkende waarden RMM\\11205232-007-ZKS-0003_v0.1-Kenmerkende Waarden Rijn-Maasmonding - Over- en Onderschrijdingsfrequenties.docx
    o HKV rapport pag 5-102 = -97 113, "Methode II Conditionele Weibull fit en zichtduur": p:\\11208031-010-kenmerkende-waarden-k\\literatuur\\Waterstandsfrequenties in de RMM - 2006.pdf
    o Ook goederen/Fiole (oa trendbreuk 1998): https://puc.overheid.nl/rijkswaterstaat/doc/PUC_102024_31/ (tabel die Boyan heeft gebruikt, is in HKV overgenomen en ook door Boyan overgenomen)
○ Voor bepaalde locaties waar afvoergolf rivier werkte methode van HKV het beste, Boyan heeft dit in Python gezet en veel duidelijker. Conclusies zijn in zijn rapport gezet

"""

dir_meas_overschr = os.path.join(dir_base,'data_overschrijding')
dir_vali_overschr = r'p:\archivedprojects\11205258-005-kpp2020_rmm-g5\C_Work\00_KenmerkendeWaarden\Onder_overschrijdingslijnen_Boyan\Tables'

#station_break_dict = {'HOEKVHLD':'01-01-1998'} #TODO: possible to make generic?
station_name_dict = {'HOEKVHLD':'Hoek_van_Holland'}

Tfreqs_interested = [5, 2, 1, 1/2, 1/5, 1/10, 1/20, 1/50, 1/100, 1/200,
                     1/500, 1/1000, 1/2000, 1/4000, 1/5000, 1/10000] #TODO: which frequencies are realistic with n years of data? probably remove this entire row >> met 40 jaar data kun je in principe tot 1/40 gaan, maar met weibull kun je extrapoleren en in theorie >> dit is voor tabel die je eruit wil hebben

color_map = {'Ongefilterd':  'b', 'Gefilterd': 'orange', 'Trendanalyse': 'g',
             'Weibull': 'r', 'Hydra-NL': 'm', 'Hydra-NL met modelonzekerheid': 'cyan',
             'Gecombineerd': 'k'}

temp = {}
tstarts = pd.DataFrame()
for current_station in []:#stat_list:
    print(f'overschrijdingsfrequenties for {current_station}')
    if closefigatstart:
        plt.close('all')

    file_ext_pkl = os.path.join(dir_meas_alldata,f"{current_station}_measext.pkl")
    if not os.path.exists(file_ext_pkl):
        continue
    data_pd_measext = pd.read_pickle(file_ext_pkl)
    data_pd_measext = clean_data(data_pd_measext)
    
    data_pd_measext = data_pd_measext.loc[:tstop_dt] # only include data up to year_slotgem #TODO: add trendbreuk tstart for DENOVBTN etc
    
    if len(data_pd_measext['HWLWcode'].unique()) > 2:
        data_pd_measext = hatyan.calc_HWLW12345to12(data_pd_measext) #convert 12345 to 12 by taking minimum of 345 as 2 (laagste laagwater)
    data_pd_HW = data_pd_measext.loc[data_pd_measext['HWLWcode']==1]
    data_pd_LW = data_pd_measext.loc[data_pd_measext['HWLWcode']!=1]
    
    #TODO: move this to data-check part (first/last occurrences of WaardeBepalingsmethode)
    # data_pd_measext_WBM_tstart = data_pd_measext[['WaardeBepalingsmethode.Code','WaardeBepalingsmethode.Omschrijving']].drop_duplicates(keep='first')
    # data_pd_measext_WBM_tstop = data_pd_measext[['WaardeBepalingsmethode.Code','WaardeBepalingsmethode.Omschrijving']].drop_duplicates(keep='last')
    # data_pd_measext_WBM_times = pd.concat([data_pd_measext_WBM_tstart,data_pd_measext_WBM_tstop]).sort_index()
    tstart_usefuldata = None
    
    
    station_rule_type = 'break' #TODO: compare results to the ones withouth this break or break on different date
    station_break_value = dt.datetime(1998,1,1).strftime('%Y-%m-%d') #TODO: adjust? 
    
    # 1. Exceedance
    print('Exceedance')
    dist = {}
    
    print('Calculate unfiltered distribution')
    
    df_extrema = data_pd_HW
    
    dist['Ongefilterd'] = hatyan.distribution(df_extrema.copy())
    
    """# filtering is only applicable for stations with high river discharge influence, so disabled #TODO: ext is geschikt voor getij, maar bij hoge afvoergolf wil je alleen het echte extreem. Er is dan een treshold per station nodig, is nodig om de rivierafvoerpiek te kunnen duiden.
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
        print('Load Hydra-NL distribution data')
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
    if current_station in station_name_dict.keys():
        file_vali_exeed = os.path.join(dir_vali_overschr,'Exceedance_lines',f'Exceedance_lines_{stat_name}.csv')
        if 0:#os.path.exists(file_vali_exeed):#TODO: for some reason, "OSError: [Errno 22] Invalid argument" when accessing this file
            data_vali = pd.read_csv(file_vali_exeed,sep=';')
            ax.plot(data_vali['value_Tfreq'],data_vali['value']/100,'--',label='validation')
            ax.legend(loc=4)
    
    fig.savefig(os.path.join(dir_overschrijding, f'Exceedance_lines_{current_station}.png')) #.svg
    """
    hatyan.interpolate_interested_Tfreqs_to_csv(dist['Gecombineerd'], Tfreqs=Tfreqs_interested, id=current_station,
                                              csv_dir=dir_overschrijding, prefix='Exceedance_lines')
    """
    #continue #TODO: also check and reactivate deceedance part. For instance select wl/ext dataset here also
    # 2. Deceedance
    print('Deceedance')
    dist = {}
    
    print('Calculate unfiltered distribution')
    df_extrema = data_pd_LW
    
    dist['Ongefilterd'] = hatyan.distribution(df_extrema.copy(), inverse=True) #TODO: with ext, this line is different than trendanalyse
    
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
    if current_station in station_name_dict.keys():
        file_vali = os.path.join(dir_vali_overschr,'Deceedance_lines',f'Deceedance_lines_{stat_name}.csv')
        if 0:#os.path.exists(file_vali):#TODO: for some reason, "OSError: [Errno 22] Invalid argument" when accessing this file
            data_vali = pd.read_csv(file_vali,sep=';')
            ax.plot(data_vali['value_Tfreq'],data_vali['value']/100,'--',label='validation')
            ax.legend(loc=4)
    fig.savefig(os.path.join(dir_overschrijding, f'Deceedance_lines_{current_station}.png')) #.svg
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




