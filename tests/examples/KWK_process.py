# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 14:17:13 2022

@author: veenstra
"""

import os
import sys
import glob
import pandas as pd
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')
from matplotlib import ticker
import hatyan # available via `pip install hatyan` or at https://github.com/Deltares/hatyan
import contextily as ctx #`conda install -c conda-forge contextily -y`
import statsmodels.api as sm # `conda install -c conda-forge statsmodels -y`

#TODO: apply to all measurements: remove QC==99 (always, or maybe make nans?), crop_timeseries (when applicable), NAP2005 correction?, SLR trend correctie voor overschrijdingsfrequenties en evt ook voor andere KW?
#TODO: move all parts to hatyan.kenmerkendewaarden.*, maybe also the stuff in hatyan/overschrijding.py (and include license header) >> indeed put it in hatyan or not?
#TODO: add tidal indicators (LAT etc) >> done at slotgemiddelden part
get_catalog = False

tstart_dt_DDL = dt.datetime(1870,1,1) #1870,1,1 for measall folder #TODO: HOEKVHLD contains yearmeanwl data from 1864, so is not all inclusive
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

fig_alltimes_ext = [dt.datetime.strptime(x,'%Y%m%d') for x in os.path.basename(dir_meas_alldata).split('_')[2:]]

if get_catalog:
    print('retrieving DDL catalog')
    catalog_dict = hatyan.get_DDL_catalog(catalog_extrainfo=['WaardeBepalingsmethoden','MeetApparaten','Typeringen'])
    pd.to_pickle(catalog_dict,os.path.join(dir_base,'DDL_catalog.pkl'))
    print('...done')
else:
    catalog_dict = pd.read_pickle(os.path.join(dir_base,'DDL_catalog.pkl'))
cat_locatielijst = catalog_dict['LocatieLijst']#.set_index('Locatie_MessageID',drop=True)
cat_locatielijst.to_pickle(os.path.join(dir_meas_DDL,'catalog_lokatielijst.pkl'))

#get list of stations with extremes and add K13A
cat_aquometadatalijst_sel, cat_locatielijst_sel = hatyan.get_DDL_stationmetasubset(catalog_dict=catalog_dict,station_dict=None, meta_dict={'Grootheid.Code':'WATHTE','Groepering.Code':'GETETM2'}) #TODO: this was filtered for extremes, but was not possible with list AB (should be possible after DDL synchronization)
"""
stat_list_ABTC = []
list_dia_ABTC = glob.glob('p:\\11208031-010-kenmerkende-waarden-k\\work\\data_vanRWS\\wetransfer_waterstandsgegevens_2022-08-05_1306\\WATHTE_oud\\*.dia')
for file_dia_ABTC in list_dia_ABTC:
    diablocks = hatyan.get_diablocks_startstopstation(file_dia_ABTC)
    stat_code = diablocks.loc[0,'station']
    stat_list_ABTC.append(stat_code)
print(stat_list_ABTC)
"""
stat_list_ABTC = ['A12','AWGPFM','BAALHK','BATH','BERGSDSWT','BROUWHVSGT02','BROUWHVSGT08','GATVBSLE','BRESKVHVN','CADZD','D15','DELFZL','DENHDR','EEMSHVN','EURPFM','F16','F3PFM','HARVT10','HANSWT','HARLGN','HOEKVHLD','HOLWD','HUIBGT','IJMDBTHVN','IJMDSMPL','J6','K13APFM','K14PFM','KATSBTN','KORNWDZBTN','KRAMMSZWT','L9PFM','LAUWOG','LICHTELGRE','MARLGT','NES','NIEUWSTZL','NORTHCMRT','DENOVBTN','OOSTSDE04','OOSTSDE11','OOSTSDE14','OUDSD','OVLVHWT','Q1','ROOMPBNN','ROOMPBTN','SCHAARVDND','SCHEVNGN','SCHIERMNOG','SINTANLHVSGR','STAVNSE','STELLDBTN','TERNZN','TERSLNZE','TEXNZE','VLAKTVDRN','VLIELHVN','VLISSGN','WALSODN','WESTKPLE','WESTTSLG','WIERMGDN','YERSKE']
stat_list_addnonext=['K13APFM','MAASMSMPL'] + stat_list_ABTC #two stations + missing from dialist ABCT (OVLVHWT is left out since HANSWT is there)
for stat_addnonext in stat_list_addnonext:
    if stat_addnonext in cat_locatielijst_sel['Code'].tolist():
        continue
    addnonext_entry = cat_locatielijst.loc[cat_locatielijst['Code']==stat_addnonext].set_index('Locatie_MessageID',drop=True) #K13APFM/MAASMSMPL/etc do not have extremes (at least not in DDL), so is manually added to the interest-list
    if not len(addnonext_entry)==1:
        #print(f'station name {stat_addnonext} found {len(addnonext_entry)} times, should be 1.:\n{addnonext_entry}, using last one')
        addnonext_entry = addnonext_entry.iloc[[-1]] #TODO: stations ['A12', 'D15', 'J6', 'Q1'] are now duplicated in cat_locatielijst, using last entry which corresponds to "Platform *" instead of "* Platform", which is what is in the dia. This should be fixed.
    cat_locatielijst_sel = cat_locatielijst_sel.append(addnonext_entry)
cat_locatielijst_sel['RDx'],cat_locatielijst_sel['RDy'] = hatyan.convert_coordinates(coordx_in=cat_locatielijst_sel['X'].values, coordy_in=cat_locatielijst_sel['Y'].values, epsg_in=int(cat_locatielijst_sel['Coordinatenstelsel'].iloc[0]),epsg_out=28992)
cat_locatielijst_sel_codeidx = cat_locatielijst_sel.reset_index(drop=False).set_index('Code',drop=False)

#stat_name_list = ['BATH','DELFZIJL','DEN HELDER','DORDRECHT','EEMSHAVEN','EURO PLATFORM','HANSWEERT','HARINGVLIETSLUIZEN','HARLINGEN','HOEK VAN HOLLAND','HUIBERTGAT','IJMUIDEN','KORNWERDERZAND','LAUWERSOOG','ROOMPOT BUITEN','ROTTERDAM','SCHEVENINGEN','STAVENISSE','TERNEUZEN','VLISSINGEN','WEST-TERSCHELLING'] # lijst AB
stat_name_list = ['Terneuzen','Bath','HANSWT','Vlissingen','Bergse Diepsluis west','Krammersluizen west','Stavenisse','Roompot binnen','Cadzand','Westkapelle','Roompot buiten','Brouwershavensche Gat 08','Haringvliet 10','Hoek van Holland','Scheveningen','IJmuiden buitenhaven','Petten zuid','Den Helder','Texel Noordzee','Terschelling Noordzee','Wierumergronden','Huibertgat','Oudeschild','Vlieland haven','West-Terschelling','Nes','Schiermonnikoog','Den Oever buiten','Kornwerderzand buiten','Harlingen','Lauwersoog','Eemshaven','Delfzijl','Nieuwe Statenzijl','Lichteiland Goeree','Euro platform','K13a platform'] + ['Dordrecht','Stellendam Buiten','Rotterdam'] + ['Maasmond','Oosterschelde 11'] + stat_list_addnonext[2:] #"KW kust en GR Dillingh 2013" en "KW getijgebied RWS 2011.0", aangevuld met 3 stations AB, aangevuld met BOI wensen, aangevuld met dialijst ABCT
stat_list = []
for stat_name in stat_name_list:
    bool_isstation = cat_locatielijst_sel_codeidx['Naam'].str.contains(stat_name,case=False) | cat_locatielijst_sel_codeidx['Code'].str.contains(stat_name,case=False)
    if not bool_isstation.sum()==1:
        print(f'station name {stat_name} found {bool_isstation.sum()} times, should be 1.:\n{cat_locatielijst_sel_codeidx.loc[bool_isstation,["Naam"]]}')
    stat_list.append(cat_locatielijst_sel_codeidx.loc[bool_isstation,'Code'].iloc[0])
    #print(f'{stat_name:30s}: {bool_isstation.sum()}')
#stat_list = ['BATH','DELFZL','DENHDR','DORDT','EEMSHVN','EURPFM','HANSWT','STELLDBTN','HARLGN','HOEKVHLD','HUIBGT','IJMDBTHVN','KORNWDZBTN','LAUWOG','ROOMPBTN','ROTTDM','SCHEVNGN','STAVNSE','TERNZN','VLISSGN','WESTTSLG'] # lijst AB vertaald naar DONAR namen
#stat_list = ['HOEKVHLD','HARVT10','VLISSGN']
#stat_list2 = ['DELFZL','EEMSHVN','NIEUWSTZL','IJMDBTHVN','IJMDSMPL','SCHEVNGN','HOEKVHLD','A12','AWGPFM','D15','EURPFM','F16','J6','K13APFM','K14PFM','L9PFM','LICHTELGRE','Q1','F3PFM','BERGSDSWT','KATSBTN','KRAMMSZWT','MARLGT','ROOMPBNN','ROOMPBTN','STAVNSE','YERSKE','NORTHCMRT','SINTANLHVSGR','BROUWHVSGT02','BROUWHVSGT08','HARVT10','OOSTSDE04','OOSTSDE11','OOSTSDE14','STELLDBTN','VLAKTVDRN','HUIBGT','TERSLNZE','TEXNZE','WIERMGDN','LAUWOG','SCHIERMNOG','DENHDR','DENOVBTN','HARLGN','HOLWD','KORNWDZBTN','NES','OUDSD','VLIELHVN','WESTTSLG','BAALHK','BATH','BRESKVHVN','CADZD','GATVBSLE','HANSWT','OVLVHWT','SCHAARVDND','TERNZN','VLISSGN','WALSODN','WESTKPLE'] #list by AB 1-6-2022, ext often not available for new stations. Not in this list from previous list: PETTZD, DORDT, ROTTDM, MAASMSMPL

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
for current_station in stat_list:
    print(f'checking data for {current_station}')
    list_relevantmetadata = ['WaardeBepalingsmethode.Code','WaardeBepalingsmethode.Omschrijving','MeetApparaat.Code','MeetApparaat.Omschrijving','Hoedanigheid.Code','Grootheid.Code','Groepering.Code','Typering.Code']
    list_relevantDDLdata = ['WaardeBepalingsmethode.Code','MeetApparaat.Code','MeetApparaat.Omschrijving','Hoedanigheid.Code']
    
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
    
    #calculate monthly/yearly mean for meas wl data
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
        if timediff_ext.min() < dt.timedelta(hours=4): #TODO: fix this: for e.g. BROUWHVSGT08: ts_meas_ext_pd.loc[dt.datetime(2015,1,1):dt.datetime(2015,1,2)]
            raise Exception(f'extreme data contains values that are too close ({timediff_ext.min()}), should be at least 4 hours difference')
        metaext = pd.read_pickle(file_extmeta_pkl)
        for metakey in list_relevantmetadata:
            data_summary.loc[current_station,f'{metakey}_ext'] = '|'.join(metaext[metakey].unique())
        if str(ts_meas_ext_pd.index[0].tz) != 'Etc/GMT-1': #this means UTC+1
            raise Exception(f'measext data for {current_station} is not in expected timezone (Etc/GMT-1): {ts_meas_ext_pd.index[0].tz}')
        ts_meas_ext_pd.index = ts_meas_ext_pd.index.tz_localize(None)
        ts_meas_ext_dupltimes = ts_meas_ext_pd.index.duplicated()
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
            ts_meas_ext_2000to202102 = hatyan.calc_HWLWnumbering(ts_meas_ext_2000to202102, station=current_station) #station argument helpt bij 3 extra stations
            HWmissings = (ts_meas_ext_2000to202102.loc[ts_meas_ext_pd['HWLWcode']==1,'HWLWno'].diff().dropna()!=1).sum()
            data_summary.loc[current_station,'#HWgaps_2000to202102_ext'] = HWmissings
        except Exception as e: #"tidal wave numbering: HW/LW numbers not always increasing" and "zero-size array to reduction operation minimum which has no identity" #TODO: fix by calulate and providing station or corr_tideperiods argument? Or fix otherwise in hatyan (maybe under different project)
            print(f'ERROR: {e}')
        #continue
        
        #calculate monthly/yearly mean for meas ext data
        if len(ts_meas_ext_pd['HWLWcode'].unique()) > 2:
            data_pd_HWLW_12 = hatyan.calc_HWLW12345to12(ts_meas_ext_pd) #convert 12345 to 12 by taking minimum of 345 as 2 (laagste laagwater). first/last values are skipped if LW
        else:
            data_pd_HWLW_12 = ts_meas_ext_pd.copy()
        data_pd_HW = data_pd_HWLW_12.loc[data_pd_HWLW_12['HWLWcode']==1]
        data_pd_LW = data_pd_HWLW_12.loc[data_pd_HWLW_12['HWLWcode']==2]
        HW_mean_peryear_long = data_pd_HW.groupby(pd.PeriodIndex(data_pd_HW.index, freq="y"))['values'].mean()
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
    if os.path.exists(file_wl_png):
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
    






#### SLOTGEMIDDELDEN
for current_station in []:#stat_list:
    
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
    if not os.path.exists(file_wl_pkl):
        continue
    data_pd_meas = pd.read_pickle(file_wl_pkl)
    data_pd_meas = data_pd_meas[['values','QC']] # reduces the memory consumption significantly
    data_pd_meas.index = data_pd_meas.index.tz_localize(None)
    data_pd_meas = data_pd_meas.loc[~(data_pd_meas['QC']==99)]
    wl_mean_peryear = data_pd_meas.groupby(pd.PeriodIndex(data_pd_meas.index, freq="y"))['values'].mean()
    wl_mean_peryear.index = wl_mean_peryear.index.to_timestamp()

    #derive tidal indicators like yearmean HWLW from HWLW values
    file_ext_pkl = os.path.join(dir_meas_alldata,f"{current_station}_measext.pkl")
    if os.path.exists(file_ext_pkl):
        data_pd_HWLW_all = pd.read_pickle(file_ext_pkl)
        dict_HWLWtidalindicators = hatyan.calc_HWLWtidalindicators(data_pd_HWLW_all)
        
    #plotting (yearly averages are plotted on 1jan, would be better on 1jul)
    fig,ax1 = plt.subplots(figsize=(14,7))
    if add_validation:
        ax1.plot(yearmeanHW['values'],'+g')
        ax1.plot(yearmeanLW['values'],'+g')
        ax1.plot(yearmeanwl['values'],'+g')
        if os.path.exists(file_ext_pkl):
            yearmeanHW_diff = (yearmeanHW['values']-dict_HWLWtidalindicators['HW_mean_peryear']).dropna() #TODO: move to data check part, when validationdata for more stations is available
            yearmeanLW_diff = (yearmeanLW['values']-dict_HWLWtidalindicators['LW_mean_peryear']).dropna()
        yearmeanwl_diff = (yearmeanwl['values']-wl_mean_peryear).dropna()

    if os.path.exists(file_ext_pkl):
        ax1.plot(dict_HWLWtidalindicators['HW_mean_peryear'],'xr')
        ax1.plot(dict_HWLWtidalindicators['LW_mean_peryear'],'xr')
    ax1.plot(wl_mean_peryear,'xr')
    ax1.grid()
    ax1.set_xlim(fig_alltimes_ext) # entire period
    ax1.set_ylabel('waterstand [m]')
    ax1.set_title(f'yearly mean HW/wl/LW {current_station}')
    fig.tight_layout()
    
    # fit with sm.OLS, method from fbaart. https://github.com/openearth/sealevel/blob/master/notebooks/analysis/gtsm/nodal-tide.ipynb
    #TODO: maybe use nodal_epoch fit from same ipynb? PLSS because of potential trendbreuk? include SLR trend
    #TODO: need to include nodal tide or only trend?
    if 0:#os.path.exists(file_ext_pkl):
        mean_list = [wl_mean_peryear,dict_HWLWtidalindicators['HW_mean_peryear'],dict_HWLWtidalindicators['LW_mean_peryear']]
    else:
        mean_list = [wl_mean_peryear]
    for mean_array in mean_list:
        # We'll just use the years. This assumes that annual waterlevels are used that are stored left-padded, the mean waterlevel for 2020 is stored as 2020-1-1. This is not logical, but common practice.
        times_OLS = mean_array.index
        years_OLS = times_OLS.year.values
        if not np.allclose(np.diff(years_OLS), 1):
            print(f"SKIPPED: all years should be sequential: {years_OLS}")
            continue
        #define a simple linear model, with nodal tide and without wind and without sea-level rise
        X = np.c_[np.cos(2 * np.pi * (years_OLS - 1970) / 18.613), np.sin(2 * np.pi * (years_OLS - 1970) / 18.613),]
        # X is of shape n year x 3 parameters
        X = sm.add_constant(X)
        Y = wl_mean_peryear.values        
        # define and fit the model, fit the reanalysis nodal cycle through all stations
        model = sm.OLS(Y, X, missing="drop")
        fit = model.fit()
        OLS_wl_yearmean = fit.predict(X)
        ax1.plot(times_OLS, OLS_wl_yearmean, ".-", label="yearmean OLSpredict stat0")
            
        """
        # define the names
        names = ["Constant", "Nodal U", "Nodal V"]
        Constant = fit_meanwl.params[0]
        A = fit_meanwl.params[1]
        B = fit_meanwl.params[2]
        phase = np.arctan2(B, A)
        amplitude = np.sqrt(A ** 2 + B ** 2)
        mean = fit_meanwl.params[0]
        """
    fig.savefig(os.path.join(dir_slotgem,f'yearly_values_{current_station}'))
    plt.close()






### HAVENGETALLEN
"""
LWaggercode uitleg
TVL;1;1;hoogwater
TVL;1;2;laagwater
TVL;1;3;laagwater 1
TVL;1;4;topagger
TVL;1;5;laagwater 2
"""
culm_addtime = 2*dt.timedelta(hours=24,minutes=50)-dt.timedelta(minutes=20)+dt.timedelta(hours=1) # link with moonculmination (or M2) two days before, 24h rotates entire graph. # furthermore: 2u20min correction, this shifts the x-axis: HW is 2 days after culmination (so 4x25min difference between length of avg moonculm and length of 2 days), 20 minutes (0 to 5 meridian), 1 hour (GMT to MET)
data_pd_moonculm = hatyan.astrog_culminations(tFirst=tstart_dt-culm_addtime,tLast=tstop_dt)#,tzone='UTC+01:00')
if str(data_pd_moonculm.loc[0,'datetime'].tz) != 'UTC': # important since data_pd_HWLW['culm_hr']=range(12) hourvalues should be in UTC since that relates to the relation dateline/sun
    raise Exception(f'culmination data is not in expected timezone (UTC): {data_pd_moonculm.loc[0,"datetime"].tz}')
data_pd_moonculm['datetime'] = data_pd_moonculm['datetime'].dt.tz_localize(None)

for current_station in []:#['HARVT10', 'VLISSGN']:#stat_list:
    print(f'havengetallen for {current_station}')
    
    #read HWLW data
    file_ext_pkl = os.path.join(dir_meas,f"{current_station}_measext.pkl")
    if not os.path.exists(file_ext_pkl):
        continue
    data_pd_HWLW_all = pd.read_pickle(file_ext_pkl)
    data_pd_HWLW_all = data_pd_HWLW_all[['values','QC','HWLWcode']] #saves memory (only a bit, unless WaardeBepalingsmethode is included)
    
    #remove timezone-awareness, crop timeseries and apply NAP correction
    data_pd_HWLW_all.index = data_pd_HWLW_all.index.tz_localize(None)
    data_pd_HWLW_all = hatyan.crop_timeseries(data_pd_HWLW_all, times_ext=[tstart_dt,tstop_dt],onlyfull=False)
    if NAP2005correction:
        data_pd_HWLW_all = nap2005_correction(data_pd_HWLW_all,current_station=current_station)
    
    #check if amount of HWs is enough
    numdays = (tstop_dt-tstart_dt).total_seconds()/3600/24
    numHWs_expected = numdays*24*3600/M2_period_timedelta.total_seconds()
    numHWs = len(data_pd_HWLW_all[data_pd_HWLW_all['HWLWcode']==1])
    if numHWs < 0.95*numHWs_expected:
        raise Exception(f'ERROR: not enough high waters present in period, {numHWs} instead of >=0.95*{int(numHWs_expected):d}')
    
    print('SELECT/CALC HWLW VALUES')
    LWaggercode = 3 # timings LW aardappelgrafiek kloppen voor 1991.0 het best bij LWaggercode=3, misschien doordat eerste laagwater dominant is voor HvH. #TODO: delays should then also be used to scale with first LW in gemgetijkromme but now dominant one is used (which depends per station/period, how to automate?). Or simpler: getijkromme1991.0 "Bij meetpunten waar zich aggers voordoen, is, afgezien van de dominantie, de vorm bepaald door de ruwe krommen; dit in tegenstelling tot vroegere bepalingen. Bij spring- en doodtij is bovendien de differentiele getijduur, en daarmee de duur rijzing, afgeleid uit de ruwe krommen." 3 is sowieso niet generiek, evt ruwe kromme maken en daar dominantie uit bepalen?
    if LWaggercode == 2: #use time/value of lowest LW, 2 is actually not aggercode, but lowest LWs are converted to 2. #TODO: does not help for HOEKVHLD, what to do?
        if len(data_pd_HWLW_all['HWLWcode'].unique()) > 2:
            data_pd_HWLW = hatyan.calc_HWLW12345to12(data_pd_HWLW_all) #convert 12345 to 12 by taking minimum of 345 as 2 (laagste laagwater) #TODO: this drops first/last value if it is a LW, should be fixed
        else:
            data_pd_HWLW = data_pd_HWLW_all.copy()
    else:
        data_pd_HWLW = data_pd_HWLW_all.loc[(data_pd_HWLW_all['HWLWcode']==1) | (data_pd_HWLW_all['HWLWcode']==2) | (data_pd_HWLW_all['HWLWcode']==LWaggercode)]
    
    data_pd_HWLW.index.name = 'times' #index is called 'Tijdstip' if retrieved from DDL.
    data_pd_HWLW = data_pd_HWLW.reset_index() # needed since we need numbered HWLW, HW is a value and LW is value+1
    
    #add duur getijperiode
    HW_bool = data_pd_HWLW['HWLWcode']==1
    data_pd_HWLW['getijperiod'] = (data_pd_HWLW.loc[HW_bool,'times'].iloc[1:].values - data_pd_HWLW.loc[HW_bool,'times'].iloc[:-1])
    
    ##### CULMINATIEBEREKENING/HAVENGETALLEN
    print('select culminations corresponding to each HW/LW')
    data_pd_HWLW['culm_time'] = pd.NaT
    for iHWLW,HWLWrow in data_pd_HWLW.iterrows():
        if HWLWrow['HWLWcode']!=1: #skip non-HW rows
            continue
        #select culmination for this HW
        timediff_withculm = (HWLWrow['times']-(data_pd_moonculm['datetime']+culm_addtime)).abs()
        if timediff_withculm.min() > dt.timedelta(hours=8):
            raise Exception(f'ERROR: no culmination found within 8 hours of high water at {HWLWrow["times"]} +culm_addtime(={culm_addtime.total_seconds()/3600:.1f}hr), range culm: \n{data_pd_moonculm["datetime"].iloc[[0,-1]]}')#%s)'%(culm_time, culm_addtime.total_seconds()/3600, data_pd_HWLW.loc[[data_pd_HWLW.index.min(),data_pd_HWLW.index.max()],'times']))
        data_pd_HWLW.loc[iHWLW:iHWLW+1,'culm_time'] = data_pd_moonculm.loc[timediff_withculm.idxmin(),'datetime']
        #compute duur daling for this HW
        if iHWLW<data_pd_HWLW.index[-1]:
            data_pd_HWLW.loc[iHWLW,'duurdaling'] = data_pd_HWLW.loc[iHWLW+1,'times']-data_pd_HWLW.loc[iHWLW,'times']
    
    #compute the rest for all extremes at once
    data_pd_HWLW['culm_hr'] = (data_pd_HWLW['culm_time'].round('h').dt.hour)%12
    data_pd_HWLW['HWLW_delay'] = (data_pd_HWLW['times']-(data_pd_HWLW['culm_time']+culm_addtime))
    
    print('calculate medians per hour group for LW and HW (instead of 1991 method: average of subgroups with removal of outliers)')
    data_pd_HW = data_pd_HWLW.loc[data_pd_HWLW['HWLWcode']==1]
    data_pd_LW = data_pd_HWLW.loc[data_pd_HWLW['HWLWcode']!=1] #HWLWcode==2 or HWLWcode==LWaggercode (=3)
    HWLW_culmhr_summary = pd.DataFrame()
    HWLW_culmhr_summary['HW_values_median'] = data_pd_HW.groupby(data_pd_HW['culm_hr'])['values'].median()
    HWLW_culmhr_summary['HW_delay_median'] = data_pd_HW.groupby(data_pd_HW['culm_hr'])['HWLW_delay'].median()
    HWLW_culmhr_summary['LW_values_median'] = data_pd_LW.groupby(data_pd_LW['culm_hr'])['values'].median()
    HWLW_culmhr_summary['LW_delay_median'] = data_pd_LW.groupby(data_pd_LW['culm_hr'])['HWLW_delay'].median()
    HWLW_culmhr_summary['getijperiod_mean'] = data_pd_HW.groupby(data_pd_HW['culm_hr'])['getijperiod'].mean()
    HWLW_culmhr_summary['duurdaling_median'] = HWLW_culmhr_summary['LW_delay_median']-HWLW_culmhr_summary['HW_delay_median'] #data_pd_HW.groupby(data_pd_HW['culm_hr'])['duurdaling'].mean() gives very different result for spring/mean HARVT10
    
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
    
    file_outname = os.path.join(dir_havget, 'aardappelgrafiek_%s_%s_aggercode%s'%(year_slotgem, current_station, LWaggercode))
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






##### gemiddelde getijkrommen
# slotgemiddelden uit:
# =============================================================================
# slotGem  = 'rapportRWS'
# slotGem  = 'havengetallen2011'
slotGem  = 'havengetallen2011improved' #'rapportRWS' 'havengetallen2011' 'havengetallen2011_PLSS'
#TODO: evt schaling naar 12u25m om repetitief signaal te maken (voor boi), dan 1 plotperiode selecteren en weer terugschalen. Voorafgaand aan dit alles de ene kromme schalen met havengetallen? (Ext berekening is ingewikkelder van 1 kromme dan repetitief signaal)
#TODO: gemgetijkromme is maar 1x of 1.5x nodig voor figuur, dus verplaatsen naar 1 datum en ext afleiden (buffer_hr=0 keyword gebruiken). Voor boi av/sp/np eerst schalen naar 12h25m en interpoleren, dan repeteren, dan is alles precies even lang en makkelijk te repeteren.
fig_sum,ax_sum = plt.subplots(figsize=(14,7))
for current_station in []:#'HOEKVHLD']:#['HOEKVHLD','HARVT10']:#stat_list:
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
    
    if slotGem in ['rapportRWS','havengetallen2011','havengetallen2011_PLSS']:
        slotGem_file = slotGem
    elif slotGem=='havengetallen2011improved':
        slotGem_file = 'havengetallen2011'
    dir_vali_krommen = r'p:\archivedprojects\11205258-005-kpp2020_rmm-g5\C_Work\00_KenmerkendeWaarden\07_Figuren\figures_ppSCL_2\final20201211'
    file_vali_doodtijkromme = os.path.join(dir_vali_krommen,f'doodtijkromme_{current_station}_{slotGem_file}.csv')
    file_vali_gemtijkromme = os.path.join(dir_vali_krommen,f'gemGetijkromme_{current_station}_{slotGem_file}.csv')
    file_vali_springtijkromme = os.path.join(dir_vali_krommen,f'springtijkromme_{current_station}_{slotGem_file}.csv')        
    
    if year_slotgem not in [2011,'2011_olddata']:
        raise Exception(f'gemiddelde getijkromme only possible for 2011: {year_slotgem}') #TODO: almost not anymore
        
    #TODO: make this automatic >> also add PLSS
    if slotGem == 'havengetallen2011improved': #KW-RMM havengetallen programma (was hardcoded in script)
        file_havget = os.path.join(dir_havget,f'aardappelgrafiek_2011_{current_station}_aggercode3.csv') #TODO: aggercode 3? do not scale if aggers? (see comments at havengetallen section)
        if not os.path.exists(file_havget):
            raise Exception(f'havengetallen file does not exist: {file_havget}')
        data_havget = pd.read_csv(file_havget)
        for colname in ['HW_delay_median','LW_delay_median','getijperiod_mean','duurdaling_median']:
            data_havget[colname] = data_havget[colname].apply(lambda x: pd.Timedelta(x))
        HW_sp, LW_sp, tD_sp = data_havget.loc[0,['HW_values_median','LW_values_median','duurdaling_median']]
        HW_np, LW_np, tD_np = data_havget.loc[6,['HW_values_median','LW_values_median','duurdaling_median']]
        HW_av, LW_av, tD_av = data_havget.loc[12,['HW_values_median','LW_values_median','duurdaling_median']]
        tDiff_sp = tDiff_av = tDiff_np = None #timeshift def is disabled anyway
        
    elif slotGem == 'havengetallen2011': #KW-RMM havengetallen programma (was hardcoded in script)
        if current_station == 'HOEKVHLD':
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
            # tijdsverschil voor verplaatsing HvH-->Maasmond
            tDiff_sp = dt.timedelta(minutes=-5)
            tDiff_av = dt.timedelta(minutes=-5)
            tDiff_np = dt.timedelta(minutes=-5)
        elif current_station == 'HARVT10':
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
        else:
            raise Exception(f'station {current_station} not implemented for gemiddelde getijkromme method: {slotGem}')
    elif slotGem == 'rapportRWS': #2011.0 rapport van Douwe Dillingh
        if current_station == 'HOEKVHLD':
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
            # tijdsverschil voor verplaatsing HvH-->Maasmond
            tDiff_sp = dt.timedelta(minutes=-5)
            tDiff_av = dt.timedelta(minutes=-5)
            tDiff_np = dt.timedelta(minutes=-5)
        elif current_station == 'HARVT10':
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
        else:
            raise Exception(f'station {current_station} not implemented for gemiddelde getijkromme method: {slotGem}')
    elif slotGem == 'havengetallen2011_PLSS': #KW-RMM havengetallen programma, bewerkt met PLSS correctie van Douwe (in excelsheet)
        if current_station == 'HOEKVHLD': 
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
        elif current_station == 'HARVT10':
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
            raise Exception(f'station {current_station} not implemented for gemiddelde getijkromme method: {slotGem}')
    else:
        raise Exception(f'non-existent gemiddelde getijkromme method: {slotGem}')
        
    #load measurement data
    file_wl_pkl = os.path.join(dir_meas,f"{current_station}_measwl.pkl")
    ts_meas_pd = pd.read_pickle(file_wl_pkl)
    ts_meas_pd = ts_meas_pd[['values','QC']] # reduces the memory consumption significantly
    ts_meas_pd.index = ts_meas_pd.index.tz_localize(None)
    ts_meas_pd = ts_meas_pd.loc[~(ts_meas_pd['QC']==99)]
    ts_meas_pd = hatyan.crop_timeseries(ts_meas_pd, times_ext=[tstart_dt,tstop_dt-dt.timedelta(minutes=10)])#,onlyfull=False)
    if NAP2005correction:
        ts_meas_pd = nap2005_correction(ts_meas_pd,current_station)
    
    # =============================================================================
    # Hatyan voor 10 jaar (alle componenten voor gemiddelde getijcyclus)
    # =============================================================================
    
    const_list = hatyan.get_const_list_hatyan('year') #this should not be changed, since higher harmonics are necessary
    hatyan_settings = hatyan.HatyanSettings(nodalfactors = True,
                                            fu_alltimes = False, # False is RWS-default
                                            xfac = True, #wordt niet besproken, moet die wel aan?
                                            analysis_perperiod = 'Y',
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
    times_pred_1mnth = pd.date_range(start=dt.datetime(tstop_dt.year, 1, 1, 0, 0), end=dt.datetime(tstop_dt.year, 2, 1, 0, 0), freq='20 S') # TODO hatyan: when using <60sec, hatyan.calc_HWLW() goes wrong, probably since there is a minute-rounding thing somewhere, fix this
    prediction_av = hatyan.prediction(comp_av, times_pred_all=times_pred_1mnth, hatyan_settings=hatyan_settings)
    prediction_av_ext = hatyan.calc_HWLW(ts=prediction_av)#,calc_HWLWlocal=False)
    continue


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
    
    # bereken schalingsratio's voor kromme
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
        for i in np.arange(0,len(timesHW)):
            tide_HWtoLW = ts_corr.loc[timesHW[i]:timesLW[i]]
            ts_corr.loc[timesHW[i]:timesLW[i],'times'] = pd.date_range(start=ts_corr.loc[timesHW[i],'times'],freq=f'{int(ratioDown*1e9*60)} N',periods=len(tide_HWtoLW))
            if i == len(timesHW)-1: #not for last HW
                continue
            tide_LWtoHW = ts_corr.loc[timesLW[i]:timesHW[i+1]]
            ts_corr.loc[timesLW[i]:timesHW[i+1],'times'] = pd.date_range(start=ts_corr.loc[timesLW[i],'times'],freq=f'{int(ratioUp*1e9*60)} N',periods=len(tide_LWtoHW))
        ts_corr = ts_corr.set_index('times',drop=True)
        return ts_corr
    
    #vermenigvuldiging van kromme met ratio's
    print('vermenigvuldig_kromme gemgetij')
    prediction_av_corr = vermenigvuldig_kromme(prediction_av, idHW_av, idLW_av, rHW_av, rLW_av, rtD_av, rtU_av)
    
    fig,(ax1,ax2) = hatyan.plot_timeseries(ts=prediction_av,ts_ext=prediction_av_ext)
    ax1.plot(prediction_av_corr['values'],'r',label='gecorrigeerde kromme')
    ax1.legend(labels=['ruwe kromme','0m+NAP','gemiddelde waterstand','hoogwater','laagwater','gecorrigeerde kromme'],loc=4)
    ax1.set_ylabel('waterstand [m]')
    ax1.set_title('gemiddelde getijkromme')
    fig.savefig(os.path.join(dir_gemgetij,"gemGetijkromme_%s_%s.png"%(current_station,slotGem)))
    
    def shift_HW_tostart(ts, timesHW, tstop_dt, tDiff):
        if tDiff is not None: #this timeshift derived from old csv writing should be eliminated, None results in no change
            bool_av = ts.index>=timesHW[0]
            ts_shift = ts.loc[bool_av]
            ts_shift.index -= timesHW[0]-tstop_dt-tDiff
        else:
            ts_shift = ts.copy()
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
    
    #hatyan_settings_sn = hatyan.HatyanSettings(nodalfactors = False) #TODO: year does not matter too much (maybe it does for scaling), but nodalfactors=False does matter a bit for doodtij duration >> optionally check sensitivities?
    prediction_sn     = hatyan.prediction(comp_oneyear_sncomp, times_pred_all=times_pred_1mnth, hatyan_settings=hatyan_settings)
    
    #TODO KW-RMM2020: "In het geval van aggers is het eerste laagwater gebruikt." >> laagste laagwater wordt nu genomen, gaat niet goed met schaling
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
    
    # maak ruwe doodtijkromme (selecteer getijslag na minimale HW)
    time_lowestHW = prediction_sn_ext['values'][idHW_sn].idxmin()
    zerocrossing_lowestHW_idx = (np.abs(zero_crossings_times-time_lowestHW)).argmin()
    in1 = zero_crossings_times[zerocrossing_lowestHW_idx]
    in2 = zero_crossings_times[zerocrossing_lowestHW_idx+1]
    tC_np = in2-in1
    tU_np = tC_np - tD_np
    # repeat ruwe spring/doodtijkromme in time.
    if slotGem in ['rapportRWS','havengetallen2011','havengetallen2011_PLSS']:
        prediction_sp_one = prediction_sn.loc[is1:is2].iloc[:-3] #-3 is nodig om reproductie oude lijnen te krijgen, maar dat is niet goed (moet -1 zijn) en je ziet ook een hickup bij ieder begin/eind (also for doodtij)
        prediction_np_one = prediction_sn.loc[in1:in2].iloc[:-3]
        #this repeat-method shifts the getijkromme in time, which should probably not happen
        prediction_sp = pd.DataFrame(index=prediction_sn.index)
        prediction_sp['values'] = np.tile(prediction_sp_one['values'].values,int(np.ceil(len(prediction_sn.index)/len(prediction_sp_one))))[0:len(prediction_sn.index)]
        prediction_np = pd.DataFrame(index=prediction_sn.index)
        prediction_np['values'] = np.tile(prediction_np_one['values'].values,int(np.ceil(len(prediction_sn.index)/len(prediction_np_one))))[0:len(prediction_sn.index)] 
    else: #TODO: tijd op xas in 1991.0 was uren tov HW. Dan zou bovenstaande gelden, maar dan is het ongeschikt voor BOI. (maar wel belangrijk om sp/np/av getijduur anders te hebben in getallen)
        ntide_1month_av = int(np.ceil((prediction_av.index[-1]-prediction_av.index[0])/M2_period_timedelta)*1.1) #add 1.1 factor to just add more tideperiods to be sure
        prediction_sp_one = prediction_sn.loc[is1:is2].iloc[:-1]
        prediction_np_one = prediction_sn.loc[in1:in2].iloc[:-1]
        tideperiod_sp = prediction_sp_one.index[-1]-prediction_sp_one.index[0]
        tideperiod_np = prediction_np_one.index[-1]-prediction_np_one.index[0]
        prediction_sp_more = pd.DataFrame(index=prediction_av.index)
        prediction_np_more = pd.DataFrame(index=prediction_av.index)
        for iT in range(-ntide_1month_av,ntide_1month_av+1): #repeat n times #TODO: pre-scaling of sp/np ts to M2/av length or use original lenghts? (now the ts shift in time with tideperiod_sp/np or it has gaps with M2_period_timedelta)
            pred_sp_rep = pd.DataFrame({'values':prediction_sp_one['values'].values},index=prediction_sp_one.index+iT*M2_period_timedelta)#tideperiod_sp)
            prediction_sp_more = prediction_sp_more.append(pred_sp_rep)
            pred_np_rep = pd.DataFrame({'values':prediction_np_one['values'].values},index=prediction_np_one.index+iT*M2_period_timedelta)#tideperiod_np)
            prediction_np_more = prediction_np_more.append(pred_np_rep)
        #interpolate to 1min values and take original 1 month subset (sorting index first is required)
        prediction_sp = prediction_sp_more.sort_index().interpolate(method='index').loc[prediction_av.index]
        prediction_np = prediction_np_more.sort_index().interpolate(method='index').loc[prediction_av.index]
        #drop duplicate whole-minutes values
        prediction_sp = prediction_sp[~prediction_sp.index.duplicated(keep='first')]
        prediction_np = prediction_np[~prediction_np.index.duplicated(keep='first')]
            
    #calculating extremes
    prediction_sp_ext = hatyan.calc_HWLW(ts=prediction_sp)
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
    
    #vermenigvuldiging van kromme met ratios
    print('vermenigvuldig_kromme springtij')
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
    
    #vermenigvuldiging van kromme met ratios
    print('vermenigvuldig_kromme doodtij')
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
    ax_sum.plot(prediction_sp['values'],'--',label=f'sp kromme {current_station}')
    ax_sum.plot(prediction_np['values'],'-.',label=f'np kromme {current_station}')
    ax_sum.plot(prediction_sp_corr['values'],'--',label=f'sp kromme corr {current_station}')
    ax_sum.plot(prediction_np_corr['values'],'-.',label=f'np kromme corr {current_station}')

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
            data_pd_HWLW_12 = hatyan.calc_HWLW12345to12(data_pd_measext) #convert 12345 to 12 by taking minimum of 345 as 2 (laagste laagwater) #TODO: this drops first/last value if it is a LW, should be fixed
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
    data_pd_meas = data_pd_meas.loc[~(data_pd_meas['QC']==99)] #TODO: maybe also drop duplicate times

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
    continue #TODO: skips rest of function
    fig.savefig(os.path.join(dir_overschrijding, f'Exceedance_lines_{current_station}.png')) #.svg
    plt.close(fig)
    """
    hatyan.interpolate_interested_Tfreqs_to_csv(dist['Gecombineerd'], Tfreqs=Tfreqs_interested, id=current_station,
                                              csv_dir=dir_overschrijding, prefix='Exceedance_lines')
    """
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




