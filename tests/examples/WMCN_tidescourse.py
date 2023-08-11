# -*- coding: utf-8 -*-
"""
predictie_2019_frommergedcomp.py
hatyan master configfile
voor alle stations indien mogelijk:
    - inlezen analyseresultatenbestand
    - predictie maken

"""

import os
import numpy as np
import datetime as dt
import hatyan
import pandas as pd
import matplotlib.pyplot as plt
plt.close('all')


# predictin M2 / S2, spring/neap cycle
if 1:
    dir_testdata = 'C:\\DATA\\hatyan_data_acceptancetests'
    
    stat_list = ['HOEKVHLD']#,'DENHDR','IJMDBTHVN'] #'K13APFM'
    
    times_ext_pred = [dt.datetime(2010,1,1),dt.datetime(2010,2,1)]
    times_ext_twoweeks = [dt.datetime(2010,1,1),dt.datetime(2010,1,16)]
    times_ext_somedays = [dt.datetime(2010,1,9),dt.datetime(2010,1,16)]
    times_step_pred = 10
    
    for current_station in stat_list:
    
        #component groups
        file_data_comp0 = os.path.join(dir_testdata,'predictie2019','%s_ana.txt'%(current_station))
        COMP_merged = hatyan.read_components(filename=file_data_comp0)
        
        #prediction and validation
        bool_end1 = COMP_merged.index.astype(str).str.endswith('1')
        bool_end2 = COMP_merged.index.astype(str).str.endswith('2')
        bool_end4 = COMP_merged.index.astype(str).str.endswith('4')
        bool_Mstar = COMP_merged.index.isin([f'M{num}' for num in [1,2,3,4,5,6,7,8,9,12,11,12]])
        bool_Sstar = COMP_merged.index.isin([f'S{num}' for num in [1,2,3,4,5,6,7,8,9,12,11,12]])
        bool_Msome = COMP_merged.index.isin([f'M{num}' for num in [1,2,3,5,6,7,8,9,12,11,12]])
        
        #ts_prediction_M2_nonodal = hatyan.prediction(comp=COMP_merged.loc[['M2']], nodalfactors=False, xfac=False, fu_alltimes=True, times_ext=times_ext_pred, timestep_min=times_step_pred)
        ts_prediction = hatyan.prediction(comp=COMP_merged, nodalfactors=True, xfac=False, fu_alltimes=True, times_ext=times_ext_pred, timestep_min=times_step_pred)
        ts_prediction_20 = hatyan.prediction(comp=COMP_merged.loc[['M2','S2','M4','N2','O1','MS4','A0','SA','MU2','K1','2MN2','MN4','K2','NU2','M6','Q1','2MS6','MK4','P1','3MS8']], nodalfactors=True, xfac=False, fu_alltimes=True, times_ext=times_ext_pred, timestep_min=times_step_pred)
        ts_prediction_M2 = hatyan.prediction(comp=COMP_merged.loc[['M2']], nodalfactors=True, xfac=False, fu_alltimes=True, times_ext=times_ext_pred, timestep_min=times_step_pred)
        ts_prediction_M4 = hatyan.prediction(comp=COMP_merged.loc[['M4']], nodalfactors=True, xfac=False, fu_alltimes=True, times_ext=times_ext_pred, timestep_min=times_step_pred)
        ts_prediction_end1 = hatyan.prediction(comp=COMP_merged.loc[bool_end1], nodalfactors=True, xfac=False, fu_alltimes=True, times_ext=times_ext_pred, timestep_min=times_step_pred)
        ts_prediction_end2 = hatyan.prediction(comp=COMP_merged.loc[bool_end2], nodalfactors=True, xfac=False, fu_alltimes=True, times_ext=times_ext_pred, timestep_min=times_step_pred)
        ts_prediction_end4 = hatyan.prediction(comp=COMP_merged.loc[bool_end4], nodalfactors=True, xfac=False, fu_alltimes=True, times_ext=times_ext_pred, timestep_min=times_step_pred)
        #ts_prediction_M1 = hatyan.prediction(comp=COMP_merged.loc[['M1']], nodalfactors=True, xfac=False, fu_alltimes=True, times_ext=times_ext_pred, timestep_min=times_step_pred)
        ts_prediction_S2 = hatyan.prediction(comp=COMP_merged.loc[['S2']], nodalfactors=True, xfac=False, fu_alltimes=True, times_ext=times_ext_pred, timestep_min=times_step_pred)
        ts_prediction_S4 = hatyan.prediction(comp=COMP_merged.loc[['S4']], nodalfactors=True, xfac=False, fu_alltimes=True, times_ext=times_ext_pred, timestep_min=times_step_pred)
        ts_prediction_Mstar = hatyan.prediction(comp=COMP_merged.loc[bool_Mstar], nodalfactors=True, xfac=False, fu_alltimes=True, times_ext=times_ext_pred, timestep_min=times_step_pred)
        ts_prediction_Sstar = hatyan.prediction(comp=COMP_merged.loc[bool_Sstar], nodalfactors=True, xfac=False, fu_alltimes=True, times_ext=times_ext_pred, timestep_min=times_step_pred)
        ts_prediction_Msome = hatyan.prediction(comp=COMP_merged.loc[bool_Msome], nodalfactors=True, xfac=False, fu_alltimes=True, times_ext=times_ext_pred, timestep_min=times_step_pred)
        
        #fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction, ts_validation=None)
        fig,(ax1) = plt.subplots(1,1,figsize=(10,4),sharex=True,sharey=True)
        ax1.set_title(f'maansgetij (M2) en zonsgetij (S2) {current_station}')
        ax1.plot(ts_prediction_M2,linewidth=1,label='M2')
        ax1.plot(ts_prediction_S2,linewidth=1,label='S2')
        ax1.plot(ts_prediction_Mstar,linewidth=1,label='Moon')
        ax1.plot(ts_prediction_Sstar,linewidth=1,label='Sun')
        ax1.legend(loc=1)
        #ax2.legend(loc=1)
        ax1.set_xlim(times_ext_pred)
        ax1.set_xlim(times_ext_somedays)
        fig.tight_layout()
        fig.savefig(f'moonsun_{current_station}')
        
        #fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction, ts_validation=None)
        fig,(ax1,ax2) = plt.subplots(2,1,figsize=(10,6),sharex=True,sharey=True)
        ax1.set_title(f'maansgetij (M2) en zonsgetij (S2) {current_station}')
        ax1.plot(ts_prediction_M2,linewidth=1,label='M2')
        ax1.plot(ts_prediction_S2,linewidth=1,label='S2')
        ax2.set_title(f'samengesteld getij (M2+S2) {current_station}')
        ax2.plot(ts_prediction_M2+ts_prediction_S2,linewidth=1,label='M2+S2')
        ax1.legend(loc=1)
        ax2.legend(loc=1)
        ax1.set_xlim(times_ext_pred)
        fig.tight_layout()
        fig.savefig(f'springneap_{current_station}')
        
        #fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction, ts_validation=None)
        fig,(ax1,ax2) = plt.subplots(2,1,figsize=(10,6),sharex=True,sharey=True)
        ax1.set_title(f'samengesteld (M2+S2) eenmaaldaags getij (*1) {current_station}')
        #ax1.plot(ts_prediction_M2_nonodal,linewidth=1,label='M2_nonodal')
        ax1.plot(ts_prediction_M2+ts_prediction_S2,linewidth=1,label='M2+S2')
        ax1.plot(ts_prediction_end1,linewidth=1,label='*1')
        ax2.set_title(f'combinatie (M2+S2+*1) {current_station}')
        ax2.plot(ts_prediction_M2+ts_prediction_S2+ts_prediction_end1,linewidth=1,label='M2+S2+*1')
        ax1.legend(loc=1)
        ax1.grid()
        ax2.legend(loc=1)
        ax2.grid()
        ax1.set_xlim(times_ext_somedays)
        fig.tight_layout()
        fig.savefig(f'dagelijkseongelijkheid_{current_station}')
        
        
        fig,(ax1) = plt.subplots(1,1,figsize=(10,6),sharex=True,sharey=True)
        ax1.set_title(f'fullset vs M2 {current_station}')
        ax1.plot(ts_prediction,linewidth=1,label='full set')
        ax1.plot(ts_prediction_M2,linewidth=1,label='M2')
        ax1.legend(loc=1)
        ax1.grid()
        ax1.set_xlim(times_ext_twoweeks)
        ax1.set_ylim(-1.1,1.6)
        fig.tight_layout()
        fig.savefig(f'fullset_vs_M2_{current_station}')
        
        fig,(ax1) = plt.subplots(1,1,figsize=(10,6),sharex=True,sharey=True)
        ax1.set_title(f'fullset vs M2+S2 {current_station}')
        ax1.plot(ts_prediction,linewidth=1,label='full set')
        ax1.plot(ts_prediction_M2+ts_prediction_S2,linewidth=1,label='M2+S2')
        ax1.legend(loc=1)
        ax1.grid()
        ax1.set_xlim(times_ext_twoweeks)
        ax1.set_ylim(-1.1,1.6)
        fig.tight_layout()
        fig.savefig(f'fullset_vs_M2S2_{current_station}')
        
        fig,(ax1) = plt.subplots(1,1,figsize=(10,6),sharex=True,sharey=True)
        ax1.set_title(f'fullset vs M2+S2+M4 {current_station}')
        ax1.plot(ts_prediction,linewidth=1,label='full set')
        ax1.plot(ts_prediction_M2+ts_prediction_S2+ts_prediction_M4,linewidth=1,label='M2+S2+M4')
        ax1.legend(loc=1)
        ax1.grid()
        ax1.set_xlim(times_ext_twoweeks)
        ax1.set_ylim(-1.1,1.6)
        fig.tight_layout()
        fig.savefig(f'fullset_vs_M2S2M4_{current_station}')
        
        fig,(ax1) = plt.subplots(1,1,figsize=(10,6),sharex=True,sharey=True)
        ax1.set_title(f'fullset vs 20 components {current_station}')
        ax1.plot(ts_prediction,linewidth=1,label='full set')
        ax1.plot(ts_prediction_20,linewidth=1,label='20 components')
        ax1.legend(loc=1)
        ax1.grid()
        ax1.set_xlim(times_ext_twoweeks)
        ax1.set_ylim(-1.1,1.6)
        fig.tight_layout()
        fig.savefig(f'fullset_vs_20comp_{current_station}')
        
        """ #with difference plot, but would be more valuable with actual measurement
        fig,(ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction_M2,ts_validation=ts_prediction)
        ax1.set_title(f'fullset vs M2 {current_station}')
        ax1.legend(['full_set','20 components','difference','mean'],loc=4)
        ax1.set_xlim(times_ext_twoweeks)
        ax1.set_ylim(-1.1,1.6)
        ax2.set_ylim(-0.5,0.5)
        fig.savefig(f'fullset_vs_M2_{current_station}')

        fig,(ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction_M2+ts_prediction_S2,ts_validation=ts_prediction)
        ax1.set_title(f'fullset vs M2+S2 {current_station}')
        ax1.legend(['full_set','M2+S2','difference','mean'],loc=4)
        ax1.set_xlim(times_ext_twoweeks)
        ax1.set_ylim(-1.1,1.6)
        ax2.set_ylim(-0.5,0.5)
        fig.savefig(f'fullset_vs_M2S2_{current_station}')

        fig,(ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction_M2+ts_prediction_S2+ts_prediction_M4,ts_validation=ts_prediction)
        ax1.set_title(f'fullset vs M2+S2+M4 {current_station}')
        ax1.legend(['full_set','M2+S2+M4','difference','mean'],loc=4)
        ax1.set_xlim(times_ext_twoweeks)
        ax1.set_ylim(-1.1,1.6)
        ax2.set_ylim(-0.5,0.5)
        fig.savefig(f'fullset_vs_M2S2M4_{current_station}')
        
        fig,(ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction_20,ts_validation=ts_prediction)
        ax1.set_title(f'fullset vs 20 components {current_station}')
        ax1.legend(['full_set','20 components','difference','mean'],loc=4)
        ax1.set_xlim(times_ext_twoweeks)
        ax1.set_ylim(-1.1,1.6)
        ax2.set_ylim(-0.5,0.5)
        fig.savefig(f'fullset_vs_20comp_{current_station}')
        """

########################
# prediction from comp, tried to show LAT
if 0:
    #dir_testdata = 'P:\\1209447-kpp-hydraulicaprogrammatuur\\hatyan\\hatyan_data_acceptancetests'
    dir_testdata = 'C:\\DATA\\hatyan_data_acceptancetests'
    
    current_station = 'HOEKVHLD'
    nodalfactors = True
    xfac=True
    analysis_perperiod='Y'
    #constituent list
    if current_station in ['D15','F3PFM','K14PFM','MAESLKRZZDE','Q1','A12','AWGPFM','F16','J6','L9PFM']:
        const_list = hatyan.get_const_list_hatyan('month') #21 const, potentially extended with component splitting (5 components) and SA+SM
    elif current_station in ['AMLAHVN']:
        const_list = hatyan.get_const_list_hatyan('halfyear') #88 const
    else:
        const_list = hatyan.get_const_list_hatyan('year') #94 const
    
    file_data_comp0 = os.path.join(dir_testdata,'predictie2019','%s_ana.txt'%(current_station))
    times_ext_pred = [dt.datetime(2000,1,1),dt.datetime(2019,12,31,23,50)]
    times_step_pred = 60
        
    #component groups
    COMP_merged = hatyan.read_components(filename=file_data_comp0)
    
    #prediction and validation
    ts_prediction = hatyan.prediction(comp=COMP_merged, nodalfactors=nodalfactors, xfac=xfac, fu_alltimes=False, times_ext=times_ext_pred, timestep_min=times_step_pred)
    A0_allyears = ts_prediction.groupby(by=pd.Grouper(freq='Y')).mean()
    
    #fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction, ts_validation=None)
    fig,(ax1,ax2) = plt.subplots(2,1,figsize=(15,7))
    ax1.plot(ts_prediction,linewidth=0.7)
    ax2.plot(A0_allyears.index,A0_allyears)
    fig.tight_layout()




########################
# analysis form long meas, tried to show Nodal cycle >> failed
if 0:
    current_station = 'HOEKVHLD'
    file_pkl = os.path.join(r'p:\11208031-010-kenmerkende-waarden-k\work\measurements_wl_18700101_20220101',f'{current_station}_measwl.pkl')
    data_pkl = pd.read_pickle(file_pkl)
    ts_meas = data_pkl[['values']]
    ts_meas.index = ts_meas.index.tz_localize(None)
    ts_meas = hatyan.crop_timeseries(ts_meas,times_ext=[dt.datetime(1972,1,1),dt.datetime(2019,12,31,23,50)])#,onlyfull=False)
    #bool_duplicatetimes = ts_meas.index.duplicated(keep='first')
    #ts_meas = ts_meas.loc[~bool_duplicatetimes]
        
    A0_allyears_meas = ts_meas['values'].groupby(by=pd.Grouper(freq='Y')).mean()
    min_allyears_meas = ts_meas['values'].groupby(by=pd.Grouper(freq='Y')).min()
    max_allyears_meas = ts_meas['values'].groupby(by=pd.Grouper(freq='Y')).max()

    comp_avg, comp_allperiods = hatyan.analysis(ts_meas,const_list='year',analysis_perperiod='Y',return_allperiods=True)
    ts_pred_py = hatyan.prediction_perperiod(comp_allperiods, timestep_min=20)
    ts_pred_ext = hatyan.calc_HWLW(ts_pred_py)
    if 0: #dropp close extremes (why do they occur with 5/10/30min interval??)
        idx_smalltdiff = np.where((ts_pred_ext.index[1:]-ts_pred_ext.index[:-1])<dt.timedelta(hours=2))[0]
        ts_pred_ext = ts_pred_ext.drop(ts_pred_ext.index[idx_smalltdiff],axis=0)
    ts_pred_ext = hatyan.calc_HWLWnumbering(ts_pred_ext, station=current_station)
    ts_pred_HW = ts_pred_ext.loc[ts_pred_ext['HWLWcode']==1].copy()
    ts_pred_HW['times'] = ts_pred_HW.index
    ts_pred_HW = ts_pred_HW.set_index('HWLWno')
    ts_pred_LW = ts_pred_ext.loc[ts_pred_ext['HWLWcode']==2].copy()
    ts_pred_LW['times'] = ts_pred_LW.index
    ts_pred_LW = ts_pred_LW.set_index('HWLWno')
    ts_pred_range = (ts_pred_HW-ts_pred_LW).loc[ts_pred_HW.index]
    ts_pred_range['times'] = ts_pred_HW.loc[ts_pred_range.index]['times']
    ts_pred_range = ts_pred_range.set_index('times')
    
    A0_allyears_pred = ts_pred_py['values'].groupby(by=pd.Grouper(freq='Y')).mean()
    min_allyears_pred = ts_pred_py['values'].groupby(by=pd.Grouper(freq='Y')).min()
    max_allyears_pred = ts_pred_py['values'].groupby(by=pd.Grouper(freq='Y')).max()
    max_allyears_predrange = ts_pred_range['values'].groupby(by=pd.Grouper(freq='Y')).max()
    

    fig, (ax2) = plt.subplots()
    ax2.plot(max_allyears_predrange-A0_allyears_pred,label='A0_allyears_predrange')
    #ax1.plot(ts_pred_py)
    ax2.plot(A0_allyears_pred,label='A0_allyears_meas')
    #ax2.plot(A0_allyears_pred,label='A0_allyears_pred')
    ax2.plot(min_allyears_pred-A0_allyears_pred,label='min_allyears_meas')
    #ax2.plot(min_allyears_pred,label='min_allyears_pred')
    ax2.plot(max_allyears_pred-A0_allyears_pred,label='max_allyears_meas')
    #ax2.plot(max_allyears_pred,label='max_allyears_pred')
    ax2.set_ylim(-0.2,0.2)
    ax2.legend()
    
    
    fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_meas, ts_validation=None)
    #ax1.plot(ts_pred_py)
    #ax2.plot(A0_allyears_meas,label='A0_allyears_meas')
    #ax2.plot(A0_allyears_pred,label='A0_allyears_pred')
    #ax2.plot(min_allyears_meas,label='min_allyears_meas')
    #ax2.plot(min_allyears_pred,label='min_allyears_pred')
    ax2.plot(max_allyears_meas,label='max_allyears_meas')
    #ax2.plot(max_allyears_pred,label='max_allyears_pred')
    ax2.set_ylim(max_allyears_meas.min()-0.02,max_allyears_meas.max()+0.02)
    ax2.legend()




if 0:
    import matplotlib.dates as mdates
    dir_testdata = 'C:\\DATA\\hatyan_data_acceptancetests'

    const_list = hatyan.get_const_list_hatyan('year')
    selected_stations = ['CADZD','BATH','VLISSGN','HOEKVHLD','IJMDBTHVN','DENHDR','TERSLNZE','SCHIERMNOG','DELFZL']
    selected_stations_names = ['Cadzand','Bath','Vlissingen','Hoek van Holland','IJmuiden Buitenhaven','Den Helder','Terschelling','Schiermonnikoog','Delfzijl']

    times_ext = [dt.datetime(2009,1,1),dt.datetime(2012,12,31,23,0)]

    tstart = dt.datetime(2019,1,6,2,30) #nieuwe maan op 6 jan
    tstart = dt.datetime(2019,1,6,13,1) #midden tussen opkomst en ondergang van de maan, dus maansdoorgang?
    times_ext_pred = [tstart, tstart+dt.timedelta(hours=1.2*24, minutes=49)]
    timestep_pred = 1
    
    fig, ax1 = plt.subplots(1,1,figsize=(10,5))
    n_colors = len(selected_stations)
    colors = plt.cm.jet(np.linspace(0,1,n_colors))
    for i_stat, current_station in enumerate(selected_stations):
        comp_frommeasurements_avg_group = hatyan.read_components(filename=os.path.join(dir_testdata,'predictie2019','%s_ana.txt'%(current_station)))
        ts_prediction = hatyan.prediction(comp=comp_frommeasurements_avg_group, times_ext=times_ext_pred, timestep_min=timestep_pred)
        vals_real = ts_prediction['values']
        times_real = ts_prediction.index
        ax1.plot(times_real, vals_real, label=current_station, color=colors[i_stat])
    
    ax1.grid()
    ax1.legend(bbox_to_anchor=(1,1))
    ax1.set_xlim(tstart,times_ext_pred[1])
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    fig.tight_layout()
    fig.savefig('tides_dutchcoast.png', dpi=250)




if 0:
    import xarray as xr
    data_nc_TS = xr.open_dataset(r'p:\1204257-dcsmzuno\1980-2017\DCSM-FM\A03_FES2014_GTSM_H1H2\DFM_OUTPUT_DCSM-FM_0_5nm\DCSM-FM_0_5nm_0000_his.nc')
    data_nc_tide = xr.open_dataset(r'p:\1204257-dcsmzuno\1980-2017\DCSM-FM\A03_FES2014_GTSM_H1H2_astro\DFM_OUTPUT_DCSM-FM_0_5nm\DCSM-FM_0_5nm_0000_his.nc')
    
    station_name_pd = pd.Series(data_nc_TS.station_name.astype(str)).str.strip()
    bool_HvH = station_name_pd.str.contains('hoekvh',case=False)
    idx_HvH = station_name_pd.loc[bool_HvH].index[0]
    
    time_slice = slice(dt.datetime(2013,1,1),dt.datetime(2013,12,31))
    
    wl_HvH_TS = data_nc_TS.waterlevel.sel(time=time_slice).isel(stations=idx_HvH)
    wl_HvH_tide = data_nc_tide.waterlevel.sel(time=time_slice).isel(stations=idx_HvH)
    
    ts_TS = pd.DataFrame({'values':wl_HvH_TS.to_numpy()},index=wl_HvH_TS.time)
    comp_set = hatyan.analysis(ts_TS,const_list='year',xfac=True)
    ts_pred = hatyan.prediction(comp_set,times_pred_all=ts_TS.index)
    
    fig,ax = plt.subplots()
    ax.plot(wl_HvH_TS.time, wl_HvH_TS, linewidth=1, label='model TS (tide+surge)')
    ax.plot(wl_HvH_tide.time, wl_HvH_tide, linewidth=1, label='model tide-only')
    ax.plot(ts_pred, linewidth=1, label='astro from TS')
    ax.plot(wl_HvH_tide.time, wl_HvH_TS-wl_HvH_tide, linewidth=1, label='model TS-tide')
    ax.plot(ts_TS-ts_pred, linewidth=1, label='TS - astro')
    
    #surgediff = wl_HvH_TS.data-wl_HvH_tide.data-(ts_TS['values'].values-ts_pred['values'].values)
    #ax.plot(wl_HvH_tide.time, surgediff, label='surgediff')
    ax.plot(wl_HvH_tide.time[[0,-1]],[0,0],'k', linewidth=0.7)
    ax.grid()
    ax.legend()
    ax.set_ylim(-1.2,2.9)
    ax.set_xlim(dt.datetime(2013,12,4),dt.datetime(2013,12,10))
    fig.tight_layout()
    fig.savefig('TS_model_astro.png', dpi=250)









