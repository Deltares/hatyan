# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 11:55:02 2023

@author: veenstra
"""


import os
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt
plt.close('all')
import hatyan

dir_meas = r'p:\11208031-010-kenmerkende-waarden-k\work\measurements_wl_18700101_20220101_dataTKdia'

def clean_data(ts_meas_pd,current_station):
    if 'HWLWcode' in ts_meas_pd.columns:
        keep_columns = ['values','QC','HWLWcode']
    else:
        keep_columns = ['values','QC']
    ts_meas_pd = ts_meas_pd[keep_columns] # reduces the memory consumption significantly in case of DDL data with a lot of metadata
    ts_meas_pd.index = ts_meas_pd.index.tz_localize(None)
    ts_meas_pd = ts_meas_pd.loc[~(ts_meas_pd['QC']==99)] #remove invalid data
    
    #optional nap correction
    return ts_meas_pd

pred_year = 2022 #2019 or 2022
offset_19y = 0

for current_station in ['HOEKVHLD']:#stat_list: #stat_list[stat_list.index('SCHEVNGN'):]: #['HOEKVHLD','DENOVBTN']:#
    plt.close('all')
    
    tstart_pred = dt.datetime(pred_year,1,1)
    tstop_pred = dt.datetime(pred_year+1,1,1)
    if pred_year==2019:
        times_ext_4y = [dt.datetime(2009,1,1),dt.datetime(2012,12,31,23,50)]
    elif pred_year==2022:
        times_ext_4y = [dt.datetime(2015,1,1),dt.datetime(2018,12,31,23,50)]
        
    
    file_astro = os.path.join(r'c:\Users\veenstra\Downloads',f'astro_HOEKVHLD_{pred_year}.pkl')
    if not os.path.exists(file_astro):
        # input parameters
        print('retrieving DDL catalog')
        catalog_dict = hatyan.get_DDL_catalog(catalog_extrainfo=['WaardeBepalingsmethoden','MeetApparaten','Typeringen'])
        print('...done')
        
        stationcode = 'HOEKVHLD'
        cat_locatielijst = catalog_dict['LocatieLijst']
        station_dict = cat_locatielijst[cat_locatielijst['Code']==stationcode].iloc[0]
        
        ts_astro, metadata, stationdata = hatyan.get_DDL_data(station_dict=station_dict,tstart_dt=tstart_pred,tstop_dt=tstop_pred,tzone='UTC+01:00',
                                                              meta_dict={'Grootheid.Code':'WATHTBRKD','Groepering.Code':'NVT'},allow_multipleresultsfor=['WaardeBepalingsmethode'])
        ts_astro['values'] = ts_astro['values']/100 #convert from cm to m
        ts_astro.index = ts_astro.index.tz_localize(None)
        ts_astro.to_pickle(file_astro)
        
        fig,(ax1,ax2) = hatyan.plot_timeseries(ts=ts_astro)
    else:
        ts_astro = pd.read_pickle(file_astro)
        
    print(f'loading data for {current_station}')
    file_wl_pkl = os.path.join(dir_meas,f"{current_station}_measwl.pkl")
    if os.path.exists(file_wl_pkl): #for slotgemiddelden, gemgetijkrommen (needs slotgem+havget)
        data_pd_meas_all = pd.read_pickle(file_wl_pkl)
        data_pd_meas_all = clean_data(data_pd_meas_all,current_station)
        #crop measurement data
        #data_pd_meas_10y = hatyan.crop_timeseries(data_pd_meas_all, times_ext=[tstart_dt,tstop_dt-dt.timedelta(minutes=10)])#,onlyfull=False)
    data_pd_meas_all_H = data_pd_meas_all.loc[data_pd_meas_all.index.minute==0]
    
    ts_19y_H = hatyan.crop_timeseries(data_pd_meas_all_H, times_ext=[dt.datetime(1976+offset_19y,1,1),dt.datetime(1995+offset_19y,1,1)-dt.timedelta(minutes=10)])
    comp_19y = hatyan.analysis(ts_19y_H,const_list=['SA','SM'],fu_alltimes=False,xfac=True)
    
    ts_4y_H = hatyan.crop_timeseries(data_pd_meas_all_H, times_ext=times_ext_4y)
    comp_4y = hatyan.get_components_from_ts(ts_4y_H,const_list='year',analysis_perperiod='Y',fu_alltimes=False,xfac=True)
    
    comp_merged = hatyan.merge_componentgroups(comp_main=comp_4y, comp_sec=comp_19y, comp_sec_list=['SA','SM'])
    
    pred = hatyan.prediction(comp=comp_merged,times_ext=[tstart_pred,tstop_pred],timestep_min=10,fu_alltimes=False,xfac=True)
    hatyan.write_tsdia(pred, station=current_station,vertref='NAP',filename=file_astro.replace('.pkl','.dia'))
    pred_re = hatyan.readts_dia(file_astro.replace('.pkl','.dia'))
    #pred2019_vali = hatyan.readts_dia(r'c:\DATA\hatyan_data_acceptancetests\predictie2019\HOEKVHLD_pre.txt')
    
    fig,(ax1,ax2) = hatyan.plot_timeseries(ts=pred,ts_validation=ts_astro)
    #fig,(ax1,ax2) = hatyan.plot_timeseries(ts=pred_re,ts_validation=ts_astro)
    
    
    
    