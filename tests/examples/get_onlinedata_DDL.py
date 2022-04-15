# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 17:03:52 2021

@author: veenstra
"""

import datetime as dt
import hatyan
hatyan.close('all')

# input parameters
tstart_dt = dt.datetime(2020,11,25,9,47,0) #quite recent period
tstop_dt = dt.datetime(2021,1,30,9,50,0)
#tstart_dt = dt.datetime(1993,8,25,9,47,0) #VLISSGN got new Waardebepalingsmethode in this year
#tstop_dt = dt.datetime(1994,11,30,9,50,0)
#tstart_dt = dt.datetime(2009,1,1) #common RWS retrieval period
#tstop_dt = dt.datetime(2012,12,31,23,50)

#file_config = os.path.realpath(__file__)
#dir_output, timer_start = hatyan.init_RWS(file_config, sys.argv, interactive_plots=False)
dir_testdata = 'C:\\DATA\\hatyan_data_acceptancetests'

catalog_dict = hatyan.get_DDL_catalog(catalog_extrainfo=['WaardeBepalingsmethoden','MeetApparaten','Typeringen'])

######### oneline waterlevel data retrieval for one station
if 0: #for RWS
    def convert_HWLWstr2num(ts_measwlHWLW,ts_measwlHWLWtype):
        """
        TVL;1;1;hoogwater
        TVL;1;2;laagwater
        TVL;1;3;laagwater 1
        TVL;1;4;topagger
        TVL;1;5;laagwater 2
        """
        ts_measwlHWLW.loc[ts_measwlHWLWtype['values']=='hoogwater','HWLWcode'] = 1
        ts_measwlHWLW.loc[ts_measwlHWLWtype['values']=='laagwater','HWLWcode'] = 2
        ts_measwlHWLW.loc[ts_measwlHWLWtype['values']=='laagwater 1','HWLWcode'] = 3
        ts_measwlHWLW.loc[ts_measwlHWLWtype['values']=='topagger','HWLWcode'] = 4
        ts_measwlHWLW.loc[ts_measwlHWLWtype['values']=='laagwater 2','HWLWcode'] = 5
        ts_measwlHWLW['HWLWcode'] = ts_measwlHWLW['HWLWcode'].astype(int)
        return ts_measwlHWLW
    
    include_extremes = True
    stationcode = 'HOEKVHLD'
    cat_locatielijst = catalog_dict['LocatieLijst']#.set_index('Locatie_MessageID',drop=True)
    station_dict = cat_locatielijst[cat_locatielijst['Code']==stationcode].iloc[0]
    
    ts_astro, metadata, stationdata = hatyan.get_DDL_data(station_dict=station_dict,tstart_dt=tstart_dt,tstop_dt=tstop_dt,tzone='UTC',
                                                          meta_dict={'Grootheid.Code':'WATHTBRKD','Groepering.Code':'NVT'},allow_multipleresultsfor=['WaardeBepalingsmethode'])
    ts_astro['values'] = ts_astro['values']/100 #convert from cm to m
    ts_measwl, metadata, stationdata = hatyan.get_DDL_data(station_dict=station_dict,tstart_dt=tstart_dt,tstop_dt=tstop_dt,
                                                           meta_dict={'Grootheid.Code':'WATHTE','Groepering.Code':'NVT'},allow_multipleresultsfor=['WaardeBepalingsmethode']) #default meta_dict selects waterlevels
    ts_measwl['values'] = ts_measwl['values']/100 #convert from cm to m
    if include_extremes:
        ts_measwlHWLW, metadata, stationdata = hatyan.get_DDL_data(station_dict=station_dict,tstart_dt=tstart_dt,tstop_dt=tstop_dt,allow_multipleresultsfor=['WaardeBepalingsmethode'],
                                                                   meta_dict={'Grootheid.Code':'WATHTE','Groepering.Code':'GETETM2'})
        ts_measwlHWLW['values'] = ts_measwlHWLW['values']/100 #convert from cm to m
        ts_measwlHWLWtype, metadata, stationdata = hatyan.get_DDL_data(station_dict=station_dict,tstart_dt=tstart_dt,tstop_dt=tstop_dt,allow_multipleresultsfor=['WaardeBepalingsmethode'],
                                                                       meta_dict={'Groepering.Code':'GETETM2','Typering.Code':'GETETTPE'})
        ts_measwlHWLW = convert_HWLWstr2num(ts_measwlHWLW,ts_measwlHWLWtype)
        
        ts_astroHWLW, metadata, stationdata = hatyan.get_DDL_data(station_dict=station_dict,tstart_dt=tstart_dt,tstop_dt=tstop_dt,allow_multipleresultsfor=['WaardeBepalingsmethode'],
                                                                  meta_dict={'Grootheid.Code':'WATHTBRKD','Groepering.Code':'GETETBRKD2'})
        ts_astroHWLW['values'] = ts_astroHWLW['values']/100 #convert from cm to m
        ts_astroHWLWtype, metadata, stationdata = hatyan.get_DDL_data(station_dict=station_dict,tstart_dt=tstart_dt,tstop_dt=tstop_dt,allow_multipleresultsfor=['WaardeBepalingsmethode'],
                                                                      meta_dict={'Groepering.Code':'GETETBRKD2','Typering.Code':'GETETTPE'})
        ts_astrolHWLW = convert_HWLWstr2num(ts_astroHWLW,ts_astroHWLWtype)

    if include_extremes:
        fig,(ax1,ax2) = hatyan.plot_timeseries(ts=ts_astro,ts_validation=ts_measwl,ts_ext=ts_astroHWLW,ts_ext_validation=ts_measwlHWLW)
        ax1.set_title('waterlevels for %s (%s)'%(stationcode, stationdata['Naam'][0]))
        ax2.set_ylim(-0.5,0.5)
    else:
        fig,(ax1,ax2) = hatyan.plot_timeseries(ts=ts_astro,ts_validation=ts_measwl)
        ax1.set_title('waterlevels for %s (%s)'%(stationcode, stationdata['Naam'][0]))
        ax2.set_ylim(-0.5,0.5)


######### simple waterlevel data retrieval for all waterlevel stations or all stations
if 0: #for CMEMS
    #list of all waterlevel stations
    #cat_aquometadatalijst_sel, cat_locatielijst_sel = hatyan.get_DDL_stationmetasubset(catalog_dict=catalog_dict,station=None,stationcolumn='Code',
    #                                                                                   meta_dict={'Grootheid.Omschrijving':'waterhoogte','Groepering.Code':'NVT'})
    cat_locatielijst_sel = catalog_dict['LocatieLijst'] #list of all stations
    #cat_locatielijst_sel = cat_locatielijst_sel[cat_locatielijst_sel['Code']=='VLISSGN']
    for iR, locatie_row in cat_locatielijst_sel.iterrows(): 
        request_output = hatyan.get_DDL_data(station_dict=locatie_row,meta_dict={'Grootheid.Code':'WATHTE','Groepering.Code':'NVT'},
                                             tstart_dt=tstart_dt,tstop_dt=tstop_dt,tzone='UTC',allow_multipleresultsfor=['WaardeBepalingsmethode'])
        if request_output is not None:
            ts_meas_pd, metadata, stationdata = request_output
            ts_meas_pd['values'] = ts_meas_pd['values']/100 #convert from cm to m
            print(stationdata['Naam'][0])
            print(ts_meas_pd)


######### more complex retrieval of selection of data from DDL from selection of stations
if 0:
    import os
    import pandas as pd
    import matplotlib.pyplot as plt
    plt.close('all')
    
    plot_stations = False
    write_measurement_files = False
    
    #cat_aquometadatalijst_sel, cat_locatielijst_sel = hatyan.get_DDL_stationmetasubset(catalog_dict=catalog_dict,station_dict={'Naam':'Hoek van holland'},meta_dict={'Grootheid.Omschrijving':'waterhoogte'})
    #cat_aquometadatalijst_sel, cat_locatielijst_sel = hatyan.get_DDL_stationmetasubset(catalog_dict=catalog_dict,station_dict={'Code':'HOEKVHLD'},meta_dict={'Grootheid.Omschrijving':'waterhoogte'})
    #cat_aquometadatalijst_sel, cat_locatielijst_sel = hatyan.get_DDL_stationmetasubset(catalog_dict=catalog_dict,station_dict={'Code':'VLISSGN'},meta_dict={'Grootheid.Omschrijving':'waterhoogte'})
    cat_aquometadatalijst_sel, cat_locatielijst_sel = hatyan.get_DDL_stationmetasubset(catalog_dict=catalog_dict,station_dict={'Code':'VLISSGN'},meta_dict={'Grootheid.Omschrijving':'waterhoogte','Groepering.Code':'NVT'})
    #cat_aquometadatalijst_sel, cat_locatielijst_sel = hatyan.get_DDL_stationmetasubset(catalog_dict=catalog_dict,station_dict=None,meta_dict={'Grootheid.Omschrijving':'waterhoogte'})
    #cat_aquometadatalijst_sel, cat_locatielijst_sel = hatyan.get_DDL_stationmetasubset(catalog_dict=catalog_dict,station_dict=None,meta_dict={'Grootheid.Omschrijving':'stroom'})
    #cat_aquometadatalijst_sel, cat_locatielijst_sel = hatyan.get_DDL_stationmetasubset(catalog_dict=catalog_dict,station_dict={'Code':'VLISSGN|HOEKVHLD'},meta_dict=None)
    
    print('Grootheid/Groepering Code/Omschrijving:\n', cat_aquometadatalijst_sel[['Grootheid.Code','Grootheid.Omschrijving','Groepering.Code']])
    print('station selection:\n', cat_locatielijst_sel[['Naam','Code']])
    
    # construct meta_dict dict, use cat_aquometadatalijst_sel to select the data you are interested in
    meta_dict = {#'Compartiment.Code':'OW', #OW (oppervlaktewater)
                 #'Eenheid.Code':'cm', #cm (centimeter)
                 #'MeetApparaat.Code':'109',
                 #'Hoedanigheid.Code':'NAP', #MSL, NAP, PLAATSLR, TAW, NVT (from cat_aquometadatalijst['Hoedanigheid.Code'])
                 'Grootheid.Code':'WATHTE', #WATHTBRKD (Waterhoogte berekend), WATHTE (Waterhoogte), WATHTEASTRO (Waterhoogte astronomisch), WATHTEVERWACHT (Waterhoogte verwacht), STROOMSHD (Stroomsnelheid), STROOMRTG (Stroomrichting) (from cat_aquometadatalijst['Grootheid.Code'])
                 'Groepering.Code':'NVT', #GETETBRKD2 (Getijextreem berekend), GETETBRKDMSL2 (Getijextreem berekend t.o.v. MSL), GETETM2 (Getijextremen), GETETMSL2 (Getijextremen t.o.v. MSL), NVT (entire waterlevel timeseries) (from cat_aquometadatalijst['Groepering.Code'])
                 #'WaardeBepalingsmethode.Code': 'other:F007', #other:F009 ('Visuele aflezing van blad'), other:F001 (Rekenkundig gemiddelde waarde over vorige 10 minuten), other:F007 (Rekenkundig gemiddelde waarde over vorige 5 en volgende 5 minuten)
                 }

    if plot_stations: #plot all stations
        file_ldb = os.path.join(dir_testdata,'other','wvs_coastline3.ldb') #WGS84 ldb is converted to RD, but does not change anything wrt to matlab converted ldb, which is good
        ldb_pd_wgs = pd.read_csv(file_ldb, delim_whitespace=True,skiprows=4,names=['x','y'],na_values=[999.999])
        x_out, y_out = hatyan.convert_coordinates(coordx_in=ldb_pd_wgs['x'].values, coordy_in=ldb_pd_wgs['y'].values, epsg_in=4326, epsg_out=28992)
        ldb_pd = pd.DataFrame({'RDx':x_out/1000, 'RDy':y_out/1000})
        
        x_out, y_out = hatyan.convert_coordinates(coordx_in=cat_locatielijst_sel['X'].values, coordy_in=cat_locatielijst_sel['Y'].values, epsg_in=25831, epsg_out=28992)
        fig,ax_map = plt.subplots()
        ax_map.plot(ldb_pd['RDx'],ldb_pd['RDy'],'-k',linewidth=0.4)
        ax_map.plot(x_out/1000,y_out/1000,'xr')
        #x_NZBN, y_NZBN = hatyan.convert_coordinates(coordx_in=4.07075, coordy_in=52.16182, epsg_in=4326, epsg_out=28992)
        #ax_map.plot(x_NZBN/1000,y_NZBN/1000,'xb')
        for statx,staty,statname in zip(x_out,y_out,cat_locatielijst_sel['Naam']):
            ax_map.text(statx/1000,staty/1000,statname)
    
    #cat_locatielijst_sel = cat_locatielijst_sel.iloc[[2]] # use only one station to speed up testing of script
    for iR, locatie_row in cat_locatielijst_sel.iterrows(): #loop over the each station
        #locatie_row['Locatie_MessageID'] = cat_locatielijst_sel.loc[iR].name
        station_dict = locatie_row
        request_output = hatyan.get_DDL_data(station_dict=station_dict,tstart_dt=tstart_dt,tstop_dt=tstop_dt,tzone='UTC+01:00',
                                             meta_dict=meta_dict,)#,allow_multipleresultsfor=['WaardeBepalingsmethode'])
        if request_output is not None:
            ts_meas_pd, metadata, stationdata = request_output
            if (metadata['Eenheid.Code']=='cm').all():
                ts_meas_pd['values'] = ts_meas_pd['values']/100 #convert from cm to m in case of waterstanden
            print(stationdata['Naam'][0])
            print(ts_meas_pd)
            if write_measurement_files and ts_meas_pd is not None:
                ts_meas_pd.to_csv('waterlevel_%s.txt'%(stationdata['Code']))
            
            # PLOTTING
            if plot_stations: #plot all stations
                x_out, y_out = hatyan.convert_coordinates(coordx_in=stationdata['X'], coordy_in=stationdata['Y'], epsg_in=25831, epsg_out=28992)
                ax_map.plot(x_out/1000,y_out/1000,'xg')
            
            if 0: #compare to local available measurements (Vlis 2009)
                ts_meas_pd.index = ts_meas_pd.index.tz_localize(None) #remove timezone info without changing the time
                ts_meas_validation = hatyan.readts_dia('c:\\Users\\veenstra\\Downloads\\VLISSGN_obs1.dia',station='VLISSGN')
                fig,(ax1,ax2) = hatyan.plot_timeseries(ts=ts_meas_pd,ts_validation=ts_meas_validation)
            elif 0: #compare to tidal prediciton of retrieved dataset
                ts_meas_pd.index = ts_meas_pd.index.tz_localize(None) #remove timezone info without changing the time
                comp, ts_pred = hatyan.analysis(ts=ts_meas_pd,const_list='month',return_prediction=True)
                fig,(ax1,ax2) = hatyan.plot_timeseries(ts=ts_meas_pd,ts_validation=ts_pred)
            elif 0: #RWS style tidal analysis (but SA/SM should come from other 1976...1994)
                ts_meas_pd.index = ts_meas_pd.index.tz_localize(None) #remove timezone info without changing the time
                ts_pred_validation = hatyan.readts_dia('c:\\Users\\veenstra\\Downloads\\VLISSGN_pre.dia',station='VLISSGN')
                comp = hatyan.get_components_from_ts(ts=ts_meas_pd,const_list='year',analysis_peryear=True,xfac=True,fu_alltimes=False)
                ts_pred = hatyan.prediction(comp=comp,times_pred_all=ts_pred_validation.index,xfac=True,fu_alltimes=False)
                fig,(ax1,ax2) = hatyan.plot_timeseries(ts=ts_pred,ts_validation=ts_pred_validation)
            else:
                fig,(ax1,ax2) = hatyan.plot_timeseries(ts=ts_meas_pd)
            ax1.set_title('%s (%s)'%(stationdata['Naam'][0],stationdata['Code'][0]))


#hatyan.exit_RWS(timer_start) #provides footer to outputfile when calling this script with python
            