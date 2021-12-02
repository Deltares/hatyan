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
tstop_dt = dt.datetime(2021,11,30,9,50,0)
#tstart_dt = dt.datetime(1993,8,25,9,47,0) #VLISSGN got new Waardebepalingsmethode in this year
#tstop_dt = dt.datetime(1994,11,30,9,50,0)
#tstart_dt = dt.datetime(2009,1,1) #common RWS retrieval period
#tstop_dt = dt.datetime(2012,12,31,23,50)

#file_config = os.path.realpath(__file__)
#dir_output, timer_start = hatyan.init_RWS(file_config, sys.argv, interactive_plots=False)
dir_testdata = 'C:\\DATA\\hatyan_data_acceptancetests'

catalog_dict = hatyan.get_DDL_catalog(catalog_extrainfo=['WaardeBepalingsmethoden','MeetApparaten'])


######### oneline waterlevel data retrieval for one station
if 0: #for RWS
    stationcode = 'HOEKVHLD'
    cat_locatielijst = catalog_dict['LocatieLijst'].set_index('Locatie_MessageID',drop=True) # TODO IMPROVEMENT: list of all stations, setting index is necessary since it 'Locatie_MessageID' is not allowed to be in the dict
    query_station = dict(cat_locatielijst[cat_locatielijst['Code']==stationcode].iloc[0])
    
    ts_astro, metadata = hatyan.get_DDL_data(query_station=query_station,query_tstart=tstart_dt,query_tstop=tstop_dt,allow_multipleresultsfor=['WaardeBepalingsmethode'],
                                             meta_dict={'Grootheid.Code':'WATHTBRKD','Groepering.Code':'NVT'})
    ts_astro['values'] = ts_astro['values']/100 #convert from cm to m
    ts_measwl, metadata = hatyan.get_DDL_data(query_station=query_station,query_tstart=tstart_dt,query_tstop=tstop_dt,allow_multipleresultsfor=['WaardeBepalingsmethode']) #default meta_dict selects waterlevels
    ts_measwl['values'] = ts_measwl['values']/100 #convert from cm to m
    """
    # TODO IMPROVEMENT: not possible to distinguish between HW and LW, since the codes are not available in the output
    ts_measwlHWLW, metadata = hatyan.get_DDL_data(query_station=query_station,query_tstart=tstart_dt,query_tstop=tstop_dt,allow_multipleresultsfor=['WaardeBepalingsmethode'],
                                                  meta_dict={'Grootheid.Code':'WATHTE','Groepering.Code':'GETETM2'})
    ts_measwlHWLW['values'] = ts_measwlHWLW['values']/100 #convert from cm to m
    ts_astroHWLW, metadata = hatyan.get_DDL_data(query_station=query_station,query_tstart=tstart_dt,query_tstop=tstop_dt,allow_multipleresultsfor=['WaardeBepalingsmethode'],
                                                 meta_dict={'Grootheid.Code':'WATHTBRKD','Groepering.Code':'GETETBRKD2'})
    ts_measwlHWLW['values'] = ts_measwlHWLW['values']/100 #convert from cm to m
    """
    fig,(ax1,ax2) = hatyan.plot_timeseries(ts=ts_astro,ts_validation=ts_measwl)
    ax1.set_title('waterlevels for %s'%(stationcode))
    ax2.set_ylim(-0.5,0.5)


######### simple waterlevel data retrieval for all waterlevel stations or all stations
if 0: #for CMEMS
    #list of all waterlevel stations
    #cat_aquometadatalijst_sel, cat_locatielijst_sel = hatyan.get_DDL_stationmetasubset(catalog_dict=catalog_dict,station=None,stationcolumn='Code',
    #                                                                                   meta_dict={'Grootheid.Omschrijving':'waterhoogte','Groepering.Code':'NVT'})
    cat_locatielijst_sel = catalog_dict['LocatieLijst'] #list of all stations
    #cat_locatielijst_sel = cat_locatielijst_sel[cat_locatielijst_sel['Code']=='VLISSGN']
    for iR, row in cat_locatielijst_sel.iterrows(): 
        query_station=dict(row)
        ts_meas_pd, metadata = hatyan.get_DDL_data(query_station=query_station,query_tstart=tstart_dt,query_tstop=tstop_dt,query_tzone='UTC',allow_multipleresultsfor=['WaardeBepalingsmethode'])
        if ts_meas_pd is not None:
            ts_meas_pd['values'] = ts_meas_pd['values']/100 #convert from cm to m
            print(ts_meas_pd)


######### more complex retrieval of selection of data from DDL from selection of stations
if 1:
    import os
    import pandas as pd
    import matplotlib.pyplot as plt
    plt.close('all')
    
    plot_stations = False
    write_measurement_files = True
    
    #cat_aquometadatalijst_sel, cat_locatielijst_sel = hatyan.get_DDL_stationmetasubset(catalog_dict=catalog_dict,station='Hoek van holland',stationcolumn='Naam',meta_dict={'Grootheid.Omschrijving':'waterhoogte'})
    #cat_aquometadatalijst_sel, cat_locatielijst_sel = hatyan.get_DDL_stationmetasubset(catalog_dict=catalog_dict,station='HOEKVHLD',stationcolumn='Code',meta_dict={'Grootheid.Omschrijving':'waterhoogte'})
    #cat_aquometadatalijst_sel, cat_locatielijst_sel = hatyan.get_DDL_stationmetasubset(catalog_dict=catalog_dict,station='VLISSGN',stationcolumn='Code',meta_dict={'Grootheid.Omschrijving':'waterhoogte'})
    cat_aquometadatalijst_sel, cat_locatielijst_sel = hatyan.get_DDL_stationmetasubset(catalog_dict=catalog_dict,station='VLISSGN',stationcolumn='Code',meta_dict={'Grootheid.Omschrijving':'waterhoogte','Groepering.Code':'NVT'})
    #cat_aquometadatalijst_sel, cat_locatielijst_sel = hatyan.get_DDL_stationmetasubset(catalog_dict=catalog_dict,station=None,stationcolumn='Code',meta_dict={'Grootheid.Omschrijving':'waterhoogte'})
    #cat_aquometadatalijst_sel, cat_locatielijst_sel = hatyan.get_DDL_stationmetasubset(catalog_dict=catalog_dict,station=None,stationcolumn='Code',meta_dict={'Grootheid.Omschrijving':'stroom'})
    #cat_aquometadatalijst_sel, cat_locatielijst_sel = hatyan.get_DDL_stationmetasubset(catalog_dict=catalog_dict,station='VLISSGN|HOEKVHLD',stationcolumn='Code',meta_dict=None)
    
    print('Grootheid/Groepering Code/Omschrijving:\n', cat_aquometadatalijst_sel[['Grootheid.Code','Grootheid.Omschrijving','Groepering.Code']])
    print('station selection:\n', cat_locatielijst_sel[['Naam','Code']])
    
    # construct query_metadata dict, use cat_aquometadatalijst_sel to select the data you are interested in
    query_metadata = {#'Compartiment.Code':'OW', #OW (oppervlaktewater)
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
        x_out, y_out = hatyan.convertcoordinates(coordx_in=ldb_pd_wgs['x'].values, coordy_in=ldb_pd_wgs['y'].values, epsg_in=4326, epsg_out=28992)
        ldb_pd = pd.DataFrame({'RDx':x_out/1000, 'RDy':y_out/1000})
        
        x_out, y_out = hatyan.convertcoordinates(coordx_in=cat_locatielijst_sel['X'].values, coordy_in=cat_locatielijst_sel['Y'].values, epsg_in=25831, epsg_out=28992)
        fig,ax_map = plt.subplots()
        ax_map.plot(ldb_pd['RDx'],ldb_pd['RDy'],'-k',linewidth=0.4)
        ax_map.plot(x_out/1000,y_out/1000,'xr')
        #x_NZBN, y_NZBN = hatyan.convertcoordinates(coordx_in=4.07075, coordy_in=52.16182, epsg_in=4326, epsg_out=28992)
        #ax_map.plot(x_NZBN/1000,y_NZBN/1000,'xb')
        for statx,staty,statname in zip(x_out,y_out,cat_locatielijst_sel['Naam']):
            ax_map.text(statx/1000,staty/1000,statname)
    
    #cat_locatielijst_sel = cat_locatielijst_sel.iloc[[2]] # use only one station to speed up testing of script
    for iR, row in cat_locatielijst_sel.iterrows(): #loop over the each station
        query_station=dict(row)
        #query_station = {cat_locatielijst_sel.loc[[iR]].index.name: str(cat_locatielijst_sel.loc[[iR]].index[0])}
        request_output = hatyan.get_DDL_data(query_station=query_station,meta_dict=query_metadata,
                                             query_tstart=tstart_dt,query_tstop=tstop_dt,query_tzone='UTC+01:00')#,allow_multipleresultsfor=['WaardeBepalingsmethode'])
        if request_output is not None:
            ts_meas_pd, metadata = request_output
            if (metadata['Eenheid.Code']=='cm').all():
                ts_meas_pd['values'] = ts_meas_pd['values']/100 #convert from cm to m in case of waterstanden
            print(ts_meas_pd)
            if write_measurement_files and ts_meas_pd is not None:
                ts_meas_pd.to_csv('waterlevel_%s.txt'%(row['Code']))
            
            # PLOTTING
            if plot_stations: #plot all stations
                x_out, y_out = hatyan.convertcoordinates(coordx_in=row['X'], coordy_in=row['Y'], epsg_in=25831, epsg_out=28992)
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
            ax1.set_title('%s (%s)'%(row['Naam'],row['Code']))


#hatyan.exit_RWS(timer_start) #provides footer to outputfile when calling this script with python
            