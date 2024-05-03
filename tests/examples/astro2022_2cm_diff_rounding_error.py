# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 12:51:56 2023

@author: veenstra
"""

import os
import pandas as pd
import datetime as dt
import ddlpy
import hatyan
hatyan.close('all')

# TODO: difference with 2022 astro ts from ddl for 2022/2023, not with 2019.
# probably due to bug in older version of hatyan2 that was used for 2022/2023 prediction (2019 was computed with hatyan1)

dir_meas = r'p:\archivedprojects\11208031-010-kenmerkende-waarden-k\work\measurements_wl_18700101_20220101_dataTKdia'

# year_list = [2019, 2020, 2021, 2022, 2023]
year_list = [2019, 2023]
for pred_year in year_list:

    for current_station in ['HOEKVHLD']:#,'DORDT']:
        
        tstart_pred = dt.datetime(pred_year,1,1)
        tstop_pred = dt.datetime(pred_year+1,1,1)
        
        if pred_year in [2019, 2020]:
            times_ext_4y = slice(dt.datetime(2009,1,1),dt.datetime(2012,12,31,23,50))
        elif pred_year in [2021, 2022, 2023]:
            times_ext_4y = slice(dt.datetime(2015,1,1),dt.datetime(2018,12,31,23,50))
            
        
        file_astro = os.path.join(r'c:\Users\veenstra\Downloads',f'astro_{current_station}_{pred_year}.pkl')
        if not os.path.exists(file_astro):
            print('retrieving DDL catalog')
            locations = ddlpy.locations()
            bool_station = locations.index.isin([current_station])
            bool_grootheid = locations["Grootheid.Code"].isin(["WATHTBRKD"])
            bool_groepering = locations["Groepering.Code"].isin(["NVT"])
            selected = locations.loc[bool_station & bool_grootheid & bool_groepering]
            
            measurements = ddlpy.measurements(selected.iloc[0], start_date=tstart_pred, end_date=tstop_pred)
            
            ts_astro = hatyan.ddlpy_to_hatyan(measurements)
            ts_astro.index = ts_astro.index.tz_localize(None)
            ts_astro.to_pickle(file_astro)
            
            fig,(ax1,ax2) = hatyan.plot_timeseries(ts=ts_astro)
        else:
            ts_astro = pd.read_pickle(file_astro)

        print(f'loading data for {current_station}')
        file_wl_pkl = os.path.join(dir_meas,f"{current_station}_measwl.pkl")
        data_pd_meas_all = pd.read_pickle(file_wl_pkl)
        data_pd_meas_all.index = data_pd_meas_all.index.tz_localize(None)
        data_pd_meas_all_H = data_pd_meas_all.loc[data_pd_meas_all.index.minute==0]
        
        ts_19y_H = hatyan.crop_timeseries(data_pd_meas_all_H, times=slice(dt.datetime(1976,1,1),dt.datetime(1995,1,1)-dt.timedelta(minutes=10)))
        comp_19y = hatyan.analysis(ts_19y_H,const_list=['SA','SM'],fu_alltimes=False,xfac=True)
        
        ts_4y_H = hatyan.crop_timeseries(data_pd_meas_all_H, times=times_ext_4y)
        comp_4y = hatyan.analysis(ts_4y_H,const_list='year',analysis_perperiod='Y',fu_alltimes=False,xfac=True)
        
        comp_merged = hatyan.merge_componentgroups(comp_main=comp_4y, comp_sec=comp_19y, comp_sec_list=['SA','SM'])
        
        pred = hatyan.prediction(comp=comp_merged,times=slice(tstart_pred,tstop_pred,10),fu_alltimes=False,xfac=True)
        
        fig,(ax1,ax2) = hatyan.plot_timeseries(ts=pred,ts_validation=ts_astro)
        ax2.set_ylim(None,0.035)
        ax2.set_xlim(tstart_pred,tstart_pred+dt.timedelta(days=40))
        fig.savefig(f'roundingdiff_{current_station}_{pred_year}.png')
