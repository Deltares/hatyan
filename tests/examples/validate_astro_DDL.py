# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 13:36:22 2022

@author: veenstra
"""

import datetime as dt
#import pandas as pd
import hatyan

selected_stations = ['VLISSGN','HOEKVHLD']

tstart_dt = dt.datetime(2019,12,24) #period begins with Gecontroleerd and ends with Ongecontroleerd for HOEKVHLD
tstop_dt = dt.datetime(2020,3,5)

print('retrieving DDL catalog')
catalog_dict = hatyan.get_DDL_catalog(catalog_extrainfo=['WaardeBepalingsmethoden','MeetApparaten','Typeringen'])
print('...done')

for current_station in selected_stations:
    
    #START OF STATION SETTINGS
    nodalfactors = True
    xfac=True
    analysis_perperiod=False#'Y'
    const_list = hatyan.get_const_list_hatyan('month')#year') #94 const
    CS_comps = None
    vertref='NAP'
    
    #DDL data retrieval
    print('retrieving data from DDL')
    stationcode = current_station
    cat_locatielijst = catalog_dict['LocatieLijst']#.set_index('Locatie_MessageID',drop=True)
    station_dict = cat_locatielijst[cat_locatielijst['Code']==stationcode].iloc[0]

    ts_astro, metadata, stationdata = hatyan.get_DDL_data(station_dict=station_dict,tstart_dt=tstart_dt,tstop_dt=tstop_dt,tzone='UTC+01:00',
                                                          meta_dict={'Grootheid.Code':'WATHTBRKD','Groepering.Code':'NVT'},allow_multipleresultsfor=['WaardeBepalingsmethode'])
    ts_astro['values'] = ts_astro['values']/100 #convert from cm to m
    ts_astro.index = ts_astro.index.tz_localize(None)
    ts_measwl, metadata, stationdata = hatyan.get_DDL_data(station_dict=station_dict,tstart_dt=tstart_dt,tstop_dt=tstop_dt,tzone='UTC+01:00',
                                                           meta_dict={'Grootheid.Code':'WATHTE','Groepering.Code':'NVT'},allow_multipleresultsfor=['WaardeBepalingsmethode']) #default meta_dict selects waterlevels
    ts_measwl['values'] = ts_measwl['values']/100 #convert from cm to m
    ts_measwl.index = ts_measwl.index.tz_localize(None)

    fig,(ax1,ax2) = hatyan.plot_timeseries(ts=ts_astro,ts_validation=ts_measwl)
    ax1.set_title('waterlevels for %s (%s)'%(stationcode, stationdata['Naam'][0]))
    leglabels = ax1.get_legend_handles_labels()[1]
    leglabels[0] = 'astro DDL'
    leglabels[1] = 'measured DDL'
    ax1.legend(leglabels)
    ax2.set_ylim(-0.5,0.5)
    
    comp_frommeas = hatyan.get_components_from_ts(ts=ts_measwl, const_list=const_list, nodalfactors=nodalfactors, xfac=xfac, fu_alltimes=False)
    #fig,(ax1,ax2) = hatyan.plot_components(comp_frommeas)
    
    #prediction and validation
    ts_prediction = hatyan.prediction(comp=comp_frommeas, nodalfactors=nodalfactors, xfac=xfac, fu_alltimes=False, times_pred_all=ts_measwl.index)

    fig,(ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction,ts_validation=ts_astro)
    ax1.set_title('waterlevels for %s (%s)'%(stationcode, stationdata['Naam'][0]))
    leglabels = ax1.get_legend_handles_labels()[1]
    leglabels[0] = 'astro hatyan'
    leglabels[1] = 'astro DDL'
    ax1.legend(leglabels)
    ax2.set_ylim(-0.5,0.5)
