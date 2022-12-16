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
import hatyan # available via `pip install hatyan` or at https://github.com/Deltares/hatyan (pip will not work since not all KWK functions are already in release)

#TODO: convert to netcdf instead of pkl, think of convenient netcdf format
#TODO: SLR trend correctie voor overschrijdingsfrequenties en evt ook voor andere KW?
#TODO: move all parts to hatyan.kenmerkendewaarden.*, maybe also the stuff in hatyan/overschrijding.py (and include license header)
#TODO: add LAT/HAT (AB needs this for RWS work)
dataTKdia = True #TODO: communicate data issues to TK (wl and ext): p:\11208031-010-kenmerkende-waarden-k\work\data_vanRWS_20220805\convert_dia2pickle_dataTK.py

tstart_dt = dt.datetime(2001,1,1)
tstop_dt = dt.datetime(2011,1,1)
NAP2005correction = False #True #TODO: define for all stations
if ((tstop_dt.year-tstart_dt.year)==10) & (tstop_dt.month==tstop_dt.day==tstart_dt.month==tstart_dt.day==1):
    year_slotgem = tstop_dt.year
else:
    year_slotgem = 'invalid'
print(f'year_slotgem: {year_slotgem}')

dir_base = r'p:\11208031-010-kenmerkende-waarden-k\work'
dir_meas = os.path.join(dir_base,'measurements_wl_18700101_20220101')
if dataTKdia:
    dir_meas += '_dataTKdia'
    
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

fig_alltimes_ext = [dt.datetime.strptime(x,'%Y%m%d') for x in os.path.basename(dir_meas).split('_')[2:4]]

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
    ts_meas_pd = ts_meas_pd[keep_columns] # reduces the memory consumption significantly in case of DDL data with a lot of metadata
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


compute_slotgem = True
compute_havengetallen = False
compute_gemgetij = False
compute_overschrijding = False

#physical_break_dict for slotgemiddelden and overschrijdingsfrequenties (maybe use everywhere, e.g. in clean_data?)
physical_break_dict = {'DENOVBTN':'1933', #laatste sluitgat afsluitdijk in 1932 
                       'HARLGN':'1933', #laatste sluitgat afsluitdijk in 1932
                       'VLIELHVN':'1933', #laatste sluitgat afsluitdijk in 1932
                       } #TODO: add physical_break for STAVNSE and KATSBTN? (Oosterscheldekering)


for current_station in ['HOEKVHLD']:#stat_list:#
    
    print(f'loading data for {current_station}')
    file_wl_pkl = os.path.join(dir_meas,f"{current_station}_measwl.pkl")
    if os.path.exists(file_wl_pkl): #for slotgemiddelden, gemgetijkrommen (needs slotgem+havget)
        data_pd_meas_all = pd.read_pickle(file_wl_pkl)
        data_pd_meas_all = clean_data(data_pd_meas_all,current_station)
    
    file_ext_pkl = os.path.join(dir_meas,f"{current_station}_measext.pkl")
    if os.path.exists(file_ext_pkl): #for slotgemiddelden, havengetallen, overschrijding
        data_pd_HWLW_all = pd.read_pickle(file_ext_pkl)
        data_pd_HWLW_all = clean_data(data_pd_HWLW_all,current_station)

    

    #### SLOTGEMIDDELDEN
    #TODO: more data is needed for proper working of fitting for some stations (2011: BAALHK, BRESKVHVN, GATVBSLE, SCHAARVDND)
    if compute_slotgem and os.path.exists(file_wl_pkl):
        print(f'slotgemiddelden for {current_station}')
        
        #calculate yearly mean
        dict_wltidalindicators = hatyan.calc_wltidalindicators(data_pd_meas_all) #TODO: indeed use all data?
        wl_mean_peryear = dict_wltidalindicators['wl_mean_peryear']
        dict_wltidalindicators_valid = hatyan.calc_wltidalindicators(data_pd_meas_all, tresh_yearlywlcount=2900) #24*365=8760 (hourly interval), 24/3*365=2920 (3-hourly interval)
        wl_mean_peryear_valid = dict_wltidalindicators_valid['wl_mean_peryear']
    
        #derive tidal indicators like yearmean HWLW from HWLW values
        if os.path.exists(file_ext_pkl):
            dict_HWLWtidalindicators = hatyan.calc_HWLWtidalindicators(data_pd_HWLW_all)
            HW_mean_peryear = dict_HWLWtidalindicators['HW_mean_peryear']
            LW_mean_peryear = dict_HWLWtidalindicators['LW_mean_peryear']
            dict_HWLWtidalindicators_valid = hatyan.calc_HWLWtidalindicators(data_pd_HWLW_all, tresh_yearlyHWLWcount=1400) #2*24*365/12.42=1410.6 (12.42 hourly extreme)
            HW_mean_peryear_valid = dict_HWLWtidalindicators_valid['HW_mean_peryear']
            LW_mean_peryear_valid = dict_HWLWtidalindicators_valid['LW_mean_peryear']
        
        #plotting (yearly averages are plotted on 1jan, would be better on 1jul)
        fig,ax1 = plt.subplots(figsize=(14,7))
        
        #get validation timeseries (yearly mean wl/HW/LW)
        station_name_dict = {'HOEKVHLD':'hoek',
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
                tstart_dt_trend = pd.Timestamp(physical_break_dict[current_station])
            else:
                tstart_dt_trend = None
            tstop_dt_trend = tstop_dt-dt.timedelta(days=1)
            mean_array_todate = mean_array.loc[tstart_dt_trend:tstop_dt_trend] #remove all values after tstop_dt (is year_slotgem)
            
            # We'll just use the years. This assumes that annual waterlevels are used that are stored left-padded, the mean waterlevel for 2020 is stored as 2020-1-1. This is not logical, but common practice.
            allyears_DTI = pd.date_range(mean_array_todate.index.min(),mean_array_todate.index.max()+dt.timedelta(days=5*360),freq='AS')
            mean_array_allyears = pd.Series(mean_array_todate,index=allyears_DTI)
            
            df = pd.DataFrame({'year':mean_array_allyears.index.year, 'height':mean_array_allyears.values}) #TODO: make functions accept mean_array instead of df as argument?
            
            # below methods are copied from https://github.com/openearth/sealevel/blob/master/slr/slr/models.py #TODO: install slr package as dependency or keep separate?
            fit, names, X = hatyan.linear_model(df, with_wind=False, with_nodal=False)
            pred_linear_nonodal = fit.predict(X)
            fit, names, X = hatyan.linear_model(df, with_wind=False)
            pred_linear_winodal = fit.predict(X)
            
            pred_pd = pd.DataFrame({'pred_linear_nonodal':pred_linear_nonodal,
                                    'pred_linear_winodal':pred_linear_winodal},
                                    index=allyears_DTI)
            ax1.plot(pred_pd, ".-", label=pred_pd.columns)
            ax1.set_prop_cycle(None) #reset matplotlib colors
            
            #2021.0 value
            if iM==0: #only for mean wl?
                pred_slotgem = pred_pd.loc[[tstop_dt],['pred_linear_winodal']]
                pred_slotgem.to_csv(os.path.join(dir_slotgem,f'slotgem_value_{current_station}.txt'))
            pred_future = pred_pd.loc[tstop_dt:,'pred_linear_winodal']
            ax1.plot(pred_future, ".k", label=f'pred_linear from {year_slotgem}')
            
        ax1.legend(loc=2)
        fig.savefig(os.path.join(dir_slotgem,f'yearly_values_{current_station}'))




    #TODO IMPORTANT: check culm_addtime and HWLWno+4 offsets. culm_addtime could also be 2 days or 2days +1h GMT-MET correction. 20 minutes seems odd since moonculm is about tidal wave from ocean
    ### HAVENGETALLEN 
    if compute_havengetallen and os.path.exists(file_ext_pkl):
        #TODO: move this first part to top
        culm_addtime = 4*dt.timedelta(hours=12,minutes=25)+dt.timedelta(hours=1)#-dt.timedelta(minutes=20) # 2d and 2u20min correction, this shifts the x-axis of aardappelgrafiek: HW is 2 days after culmination (so 4x25min difference between length of avg moonculm and length of 2 days), 1 hour (GMT to MET), 20 minutes (0 to 5 meridian, is commented now)
        data_pd_moonculm = hatyan.astrog_culminations(tFirst=tstart_dt-culm_addtime-dt.timedelta(hours=2*24),tLast=tstop_dt,dT_fortran=False)
        if str(data_pd_moonculm.loc[0,'datetime'].tz) != 'UTC': # important since data_pd_HWLW['culm_hr']=range(12) hourvalues should be in UTC since that relates to the relation dateline/sun
            raise Exception(f'culmination data is not in expected timezone (UTC): {data_pd_moonculm.loc[0,"datetime"].tz}')
        data_pd_moonculm['datetime'] = data_pd_moonculm['datetime'].dt.tz_localize(None)
        data_pd_moonculm = data_pd_moonculm.set_index('datetime',drop=False)
        data_pd_moonculm['values'] = data_pd_moonculm['type'] #dummy values for TA in hatyan.calc_HWLWnumbering()
        data_pd_moonculm['HWLWcode'] = 1 #all HW values since one every ~12h25m
        data_pd_moonculm = hatyan.calc_HWLWnumbering(data_pd_moonculm,doHWLWcheck=False) #TODO: currently w.r.t. cadzd, is that an issue? With DELFZL the matched culmination is incorrect (since far away), but that might not be a big issue
        data_pd_moonculm['HWLWno_offset'] = data_pd_moonculm['HWLWno']+4 #correlate HWLW to moonculmination 2 days before. TODO: check this offset in relation to culm_addtime.
        moonculm_idxHWLWno = data_pd_moonculm.set_index('HWLWno_offset')
    
        print(f'havengetallen for {current_station}')
        
        #crop timeseries
        data_pd_HWLW = hatyan.crop_timeseries(data_pd_HWLW_all, times_ext=[tstart_dt,tstop_dt],onlyfull=False) #TODO: this should be done at the top (or different varname)
        #data_pd_HWLW_new = data_pd_HWLW.copy()
        
        file_outname = os.path.join(dir_havget, f'aardappelgrafiek_{year_slotgem}_{current_station}')
        #check if amount of HWs is enough
        numHWs_expected = (tstop_dt-tstart_dt).total_seconds()/M2_period_timedelta.total_seconds()
        numHWs = (data_pd_HWLW['HWLWcode']==1).sum()
        if numHWs < 0.95*numHWs_expected:
            raise Exception(f'ERROR: not enough high waters present in period, {numHWs} instead of >=0.95*{int(numHWs_expected):d}')
        
        print('SELECT/CALC HWLW VALUES')
        if len(data_pd_HWLW['HWLWcode'].unique()) > 2:
            data_pd_HWLW = hatyan.calc_HWLW12345to12(data_pd_HWLW) #convert 12345 to 12 by taking minimum of 345 as 2 (laagste laagwater)
        
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
        
        print('AARDAPPELGRAFIEK')
        fig, (ax1,ax2) = hatyan.plot_aardappelgrafiek(HWLW_culmhr_summary)
        ax1.set_title(f'HW {current_station} {year_slotgem}')
        ax2.set_title(f'LW {current_station} {year_slotgem}')
        fig.savefig(file_outname)
        
        HWLW_culmhr_summary.loc['mean',:] = HWLW_culmhr_summary.mean() #add mean row to dataframe (not convenient to add immediately due to plotting with index 0-11)
        HWLW_culmhr_summary = HWLW_culmhr_summary.loc[[6,'mean',0]] #select neap/mean/springtide
        HWLW_culmhr_summary.index = ['neap','mean','spring']
        
        
        # #TODO: use tidal coefficient instead?: The tidal coefficient is the size of the tide in relation to its mean. It usually varies between 20 and 120. The higher the tidal coefficient, the larger the tidal range – i.e. the difference in water height between high and low tide. This means that the sea level rises and falls back a long way. The mean value is 70. We talk of strong tides – called spring tides – from coefficient 95.  Conversely, weak tides are called neap tides. https://escales.ponant.com/en/high-low-tide/ en https://www.manche-toerisme.com/springtij
        # #for HOEKVHLD, sp=0 is approx tc=1.2, np=6 is approx tc=0.8, av=mean is approx tc=1.0 (for HW, for LW it is different)
        # data_pd_HWLW_new = hatyan.calc_HWLWtidalrange(data_pd_HWLW_new)
        # data_pd_HWLW_new['tidalcoeff'] = data_pd_HWLW_new['tidalrange']/data_pd_HWLW_new['tidalrange'].mean()
        # data_pd_HWLW_new['tidalcoeff_round'] = data_pd_HWLW_new['tidalcoeff'].round(1)
        # TR_groupby_median = data_pd_HWLW_new.groupby('tidalcoeff_round')['tidalrange'].median()
        # HW_groupby_median = data_pd_HWLW_new.loc[data_pd_HWLW_new['HWLWcode']==1].groupby('tidalcoeff_round')['values'].median()
        # LW_groupby_median = data_pd_HWLW_new.loc[data_pd_HWLW_new['HWLWcode']==2].groupby('tidalcoeff_round')['values'].median()
        # #HWLW_culmhr_summary.loc[[0,6,'mean']]
        # # HWLW_culmhr_summary = pd.DataFrame()
        # # HWLW_culmhr_summary['HW_values_median'] = HW_groupby_median
        # # HWLW_culmhr_summary['LW_values_median'] = LW_groupby_median
        # # HWLW_culmhr_summary['tidalrange_median'] = TR_groupby_median
        # # HWLW_culmhr_summary = HWLW_culmhr_summary.loc[[0.8,1.0,1.2]] #select neap/mean/springtide
        # # HWLW_culmhr_summary.index = ['neap','mean','spring']
        
        
        #write to csv
        for colname in HWLW_culmhr_summary.columns: #round timedelta to make outputformat nicer
            if HWLW_culmhr_summary[colname].dtype == 'timedelta64[ns]':
                HWLW_culmhr_summary[colname] = HWLW_culmhr_summary[colname].round('S')
        HWLW_culmhr_summary.to_csv(file_outname+'.csv',float_format='%.3f')
    




    ##### gemiddelde getijkrommen
    pred_freq_sec = 10 #for gemgetijkromme #TODO: frequency decides accuracy of tU/tD and other timings (and is writing freq of BOI timeseries)
    if compute_gemgetij and os.path.exists(file_wl_pkl):
        """
        
        """
        print(f'gem getijkrommen for {current_station}')
        
        #TODO: add correctie havengetallen HW/LW av/sp/np met slotgemiddelde uit PLSS/modelfit (HW/LW av)
        file_havget = os.path.join(dir_havget,f'aardappelgrafiek_{year_slotgem}_{current_station}.csv')
        if not os.path.exists(file_havget):
            raise Exception(f'havengetallen file does not exist: {file_havget}')
        data_havget = pd.read_csv(file_havget,index_col=0)
        for colname in ['HW_delay_median','LW_delay_median','getijperiod_median','duurdaling_median']:
            if colname in data_havget.columns:
                data_havget[colname] = data_havget[colname].apply(lambda x: pd.Timedelta(x))
        # HW_sp, LW_sp, tD_sp = data_havget.loc['spring',['HW_values_median','LW_values_median','duurdaling_median']]
        # HW_np, LW_np, tD_np = data_havget.loc['neap',['HW_values_median','LW_values_median','duurdaling_median']]
        # HW_av, LW_av, tD_av = data_havget.loc['mean',['HW_values_median','LW_values_median','duurdaling_median']]
        HW_sp, LW_sp = data_havget.loc['spring',['HW_values_median','LW_values_median']]
        HW_np, LW_np = data_havget.loc['neap',['HW_values_median','LW_values_median']]
        HW_av, LW_av = data_havget.loc['mean',['HW_values_median','LW_values_median']]
        
        
        #crop measurement data #TODO: do on top?
        ts_meas_pd = hatyan.crop_timeseries(data_pd_meas_all, times_ext=[tstart_dt,tstop_dt-dt.timedelta(minutes=10)])#,onlyfull=False)
        
        # =============================================================================
        # Hatyan analyse voor 10 jaar (alle componenten voor gemiddelde getijcyclus) #TODO: maybe use original 4y period/componentfile instead? SA/SM should come from 19y analysis
        # =============================================================================
        const_list = hatyan.get_const_list_hatyan('year') #this should not be changed, since higher harmonics are necessary
        hatyan_settings_ana = hatyan.HatyanSettings(nodalfactors=True,
                                                    fu_alltimes=False, # False is RWS-default
                                                    xfac=True, # True is RWS-default
                                                    analysis_perperiod='Y',
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
        print(comp_av/comp_frommeasurements_avg.loc[components_av]) # values are different than 1991.0 document and differs per station while the document states "Zoals te verwachten is de verhouding per component tussen deze wortel en de oorspronkelijke amplitude voor alle plaatsen gelijk"
        
        comp_av.loc['A0'] = comp_frommeasurements_avg.loc['A0']
        times_pred_1mnth = pd.date_range(start=dt.datetime(tstop_dt.year, 1, 1, 0, 0)-dt.timedelta(hours=12), end=dt.datetime(tstop_dt.year, 2, 1, 0, 0), freq=f'{pred_freq_sec} S') #start 12 hours in advance, to assure also corrected values on desired tstart
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
        #NOTE: background on choice of components
        #below is different than provided list, these shallow ones are extra: ['S4','2SM6','M7','4MS4','2(MS)8','3M2S10','4M2S12']
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
        prediction_sn_ext = hatyan.calc_HWLWtidalrange(prediction_sn_ext)
        
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
        prediction_av_one_trefHW = hatyan.ts_to_trefHW(prediction_av_one,HWreftime=ia1) # repeating one is not necessary for av, but easier to do the same for av/sp/np
        prediction_av_corr_one = hatyan.reshape_signal(prediction_av_one, prediction_av_ext_one, HW_goal=HW_av, LW_goal=LW_av, tP_goal=None)
        prediction_av_corr_rep5 = hatyan.repeat_signal(prediction_av_corr_one, nb=2, na=2)
        prediction_av_corr_rep5_trefHW = hatyan.ts_to_trefHW(prediction_av_corr_rep5,HWreftime=ia1)
    
        print(f'reshape_signal SPRINGTIJ: {current_station}')
        prediction_sp_one_trefHW = hatyan.ts_to_trefHW(prediction_sp_one,HWreftime=is1)
        prediction_sp_corr_one = hatyan.reshape_signal(prediction_sp_one, prediction_sp_ext_one, HW_goal=HW_sp, LW_goal=LW_sp, tP_goal=None)
        prediction_sp_corr_rep5 = hatyan.repeat_signal(prediction_sp_corr_one, nb=2, na=2)
        prediction_sp_corr_rep5_trefHW = hatyan.ts_to_trefHW(prediction_sp_corr_rep5,HWreftime=is1)
        
        print(f'reshape_signal DOODTIJ: {current_station}')
        prediction_np_one_trefHW = hatyan.ts_to_trefHW(prediction_np_one,HWreftime=in1)
        prediction_np_corr_one = hatyan.reshape_signal(prediction_np_one, prediction_np_ext_one, HW_goal=HW_np, LW_goal=LW_np, tP_goal=None)
        prediction_np_corr_rep5 = hatyan.repeat_signal(prediction_np_corr_one, nb=2, na=2)
        prediction_np_corr_rep5_trefHW = hatyan.ts_to_trefHW(prediction_np_corr_rep5,HWreftime=in1)
        
        
        #12u25m timeseries for BOI computations (no relation between HW and moon, HW has to come at same time for av/sp/np tide, HW timing does differ between stations)
        print(f'reshape_signal BOI GEMGETIJ and write to csv: {current_station}')
        prediction_av_corrBOI_one = hatyan.reshape_signal(prediction_av_one, prediction_av_ext_one, HW_goal=HW_av, LW_goal=LW_av, tP_goal=pd.Timedelta(hours=12,minutes=25))
        prediction_av_corrBOI_one_roundtime = prediction_av_corrBOI_one.resample(f'{pred_freq_sec}S').nearest()
        prediction_av_corrBOI_one_roundtime.to_csv(os.path.join(dir_gemgetij,f'gemGetijkromme_BOI_{current_station}_slotgem{year_slotgem}.csv'),float_format='%.3f',date_format='%Y-%m-%d %H:%M:%S')
        prediction_av_corrBOI_repn_roundtime = hatyan.repeat_signal(prediction_av_corrBOI_one_roundtime, nb=0, na=10)
        
        print(f'reshape_signal BOI SPRINGTIJ and write to csv: {current_station}')
        prediction_sp_corrBOI_one = hatyan.reshape_signal(prediction_sp_one, prediction_sp_ext_one, HW_goal=HW_sp, LW_goal=LW_sp, tP_goal=pd.Timedelta(hours=12,minutes=25))
        prediction_sp_corrBOI_one.index = prediction_sp_corrBOI_one.index - prediction_sp_corrBOI_one.index[0] + prediction_av_corrBOI_one.index[0] #shift times to first HW from gemgetij
        prediction_sp_corrBOI_one_roundtime = prediction_sp_corrBOI_one.resample(f'{pred_freq_sec}S').nearest()
        prediction_sp_corrBOI_one_roundtime.to_csv(os.path.join(dir_gemgetij,f'springtijkromme_BOI_{current_station}_slotgem{year_slotgem}.csv'),float_format='%.3f',date_format='%Y-%m-%d %H:%M:%S')
        prediction_sp_corrBOI_repn_roundtime = hatyan.repeat_signal(prediction_sp_corrBOI_one_roundtime, nb=0, na=10)
    
        print(f'reshape_signal BOI DOODTIJ and write to csv: {current_station}')
        prediction_np_corrBOI_one = hatyan.reshape_signal(prediction_np_one, prediction_np_ext_one, HW_goal=HW_np, LW_goal=LW_np, tP_goal=pd.Timedelta(hours=12,minutes=25))
        prediction_np_corrBOI_one.index = prediction_np_corrBOI_one.index - prediction_np_corrBOI_one.index[0] + prediction_av_corrBOI_one.index[0] #shift times to first HW from gemgetij
        prediction_np_corrBOI_one_roundtime = prediction_np_corrBOI_one.resample(f'{pred_freq_sec}S').nearest()
        prediction_np_corrBOI_one_roundtime.to_csv(os.path.join(dir_gemgetij,f'doodtijkromme_BOI_{current_station}_slotgem{year_slotgem}.csv'),float_format='%.3f',date_format='%Y-%m-%d %H:%M:%S')
        prediction_np_corrBOI_repn_roundtime = hatyan.repeat_signal(prediction_np_corrBOI_one_roundtime, nb=0, na=10)
        
        
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
        
        #plot gemtij/springtij/doodtij
        ax1_boi.plot(prediction_av_corrBOI_repn_roundtime['values'],color=cmap(0),label='prediction gemtij')
        ax1_boi.plot(prediction_sp_corrBOI_repn_roundtime['values'],color=cmap(1),label='prediction springtij')
        ax1_boi.plot(prediction_np_corrBOI_repn_roundtime['values'],color=cmap(2),label='prediction doodtij')
        
        #plot validation lines if available
        dir_vali_krommen = r'p:\archivedprojects\11205258-005-kpp2020_rmm-g5\C_Work\00_KenmerkendeWaarden\07_Figuren\figures_ppSCL_2\final20201211'
        file_vali_doodtijkromme = os.path.join(dir_vali_krommen,f'doodtijkromme_{current_station}_havengetallen{year_slotgem}.csv')
        file_vali_gemtijkromme = os.path.join(dir_vali_krommen,f'gemGetijkromme_{current_station}_havengetallen{year_slotgem}.csv')
        file_vali_springtijkromme = os.path.join(dir_vali_krommen,f'springtijkromme_{current_station}_havengetallen{year_slotgem}.csv')        
        if os.path.exists(file_vali_gemtijkromme):
            data_vali_gemtij = pd.read_csv(file_vali_gemtijkromme,index_col=0,parse_dates=True)
            ax1_boi.plot(data_vali_gemtij['Water Level [m]'],'--',color=cmap(0),linewidth=0.7,label='validation KW2020 gemtij')
        if os.path.exists(file_vali_springtijkromme):
            data_vali_springtij = pd.read_csv(file_vali_springtijkromme,index_col=0,parse_dates=True)
            ax1_boi.plot(data_vali_springtij['Water Level [m]'],'--',color=cmap(1),linewidth=0.7,label='validation KW2020 springtij')
        if os.path.exists(file_vali_doodtijkromme):
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
    #    included data up to 2021-1-1 (was 1-1-2012 for all stations) >> trekt weibull krommer omhoog en dichter bij hydraNL
    #    moved to ext so: varying sampling interval of wl-data is not relevant anymore
    #    moved to ext so: resampling to tidal extremes with max/mean is not necessary anymore
    #    moved to ext so: break 1-1-1998 is now replaced with beginning of time (Boyan: in buurt van spui was er een trendbreuk, zie rapport boyan, daar ook voor HOEKVHLD toegepast maar niet per se van toepassing op kuststations)
    #       die break wordt voor trendanalyse toegepast, daardoor is de trendlijn korter dan de (on)gefilterd lijnen
    #TODO: add correction for SLR
    #    je zou zeespiegelstijging moeten verdisconteren in de (lineare) trend
    #    in rapport van HKV wordt gerefereerd naar Goederen(RWS,2003) en daar staat trend in uitgewerkt, die heeft Boyan overgenomen uit HKV rapport. (koppelen aan SLM)
    #    is voor de statistiek wel belangrijk dat die trend wordt verwijderd (anders wordt extreme freqs onderschat), referentievlak is laatste deel meetperiode.
    #    je kunt automatisch lineair corrigeren, maar RWS moet akkoord gaan over methodiek want heeft invloed op ontwerpcriteria etc (onder welke condities wel/niet corrigeren).
    
    """
    #TODO: aantekeningen gesprek Boyan
    ○ Je vertaalt niet x aantal datapunten naar frequentie, maar je zet de punten op volgorde en je rankt ze, daarvan maak je distributie, ranking en frequentie is niet 1 op 1
    ○ Max freq is 2 getij per dag, keer 365 dagen, maximale frequentie komt daarmee overeen. (on)gefilterd en trendanalys is datapunten op volgorde en frequentie, 
    ○ Lezen:
        o rapport boyan kw-rmm: n:\\Projects\\11205000\11205232\\C. Report - advise\\007 - Kenmerkende waarden RMM\\11205232-007-ZKS-0003_v0.1-Kenmerkende Waarden Rijn-Maasmonding - Over- en Onderschrijdingsfrequenties.docx
        o HKV rapport pag 5-102 = -97 113, "Methode II Conditionele Weibull fit en zichtduur": p:\\11208031-010-kenmerkende-waarden-k\\literatuur\\Waterstandsfrequenties in de RMM - 2006.pdf
        o Ook goederen/Fiole (oa trendbreuk 1998): https://puc.overheid.nl/rijkswaterstaat/doc/PUC_102024_31/ (tabel die Boyan heeft gebruikt, is in HKV overgenomen en ook door Boyan overgenomen)
    ○ Voor bepaalde locaties waar afvoergolf rivier werkte methode van HKV het beste, Boyan heeft dit in Python gezet en veel duidelijker. Conclusies zijn in zijn rapport gezet
    ○ weibull lijn begint pas bij hogere freq want die begint pas bij n-de waarde, want die gaat niet met het staartje naar beneden. is niet van toepassing voor die hoogfrequente situaties, is ontwikkeld voor extremen. te voorspellen freqs wordt met np.logspace() opgegeven.
    ○ plots beoordelen: rode lijn moet soort van verlengde zijn van groene, als die ineens omhoog piekt komt dat door hele extreme wardes die je dan vrmoedelijk ook al ziet in je groene lijn
    """
    Tfreqs_interested = [5, 2, 1, 1/2, 1/5, 1/10, 1/20, 1/50, 1/100, 1/200, #overschrijdingsfreqs
                         1/500, 1/1000, 1/2000, 1/4000, 1/5000, 1/10000] #TODO: which frequencies are realistic with n years of data? probably remove this entire row >> met 40 jaar data kun je in principe tot 1/40 gaan, maar met weibull kun je extrapoleren en in theorie >> dit is voor tabel die je eruit wil hebben
    
    if compute_overschrijding and os.path.exists(file_ext_pkl):
    
        print(f'overschrijdingsfrequenties for {current_station}')
        
        #clip data #TODO: do at top?
        data_pd_measext = data_pd_HWLW_all.loc[:tstop_dt] # only include data up to year_slotgem
        
        if len(data_pd_measext['HWLWcode'].unique()) > 2:
            data_pd_measext = hatyan.calc_HWLW12345to12(data_pd_measext) #convert 12345 to 12 by taking minimum of 345 as 2 (laagste laagwater)
        data_pd_HW = data_pd_measext.loc[data_pd_measext['HWLWcode']==1]
        data_pd_LW = data_pd_measext.loc[data_pd_measext['HWLWcode']!=1]
        
        #get Hydra-NL and KWK-RMM validation data (only for HOEKVHLD)
        dist_vali_exc = {}
        dist_vali_dec = {}
        if current_station =='HOEKVHLD':
            dir_vali_overschr = os.path.join(dir_base,'data_overschrijding')
            stat_name = 'Hoek_van_Holland'
            print('Load Hydra-NL distribution data and other validation data')
            dist_vali_exc = {}
            dist_vali_exc['Hydra-NL'] = pd.read_csv(os.path.join(dir_vali_overschr,'Processed_HydraNL','Without_model_uncertainty',f'{stat_name}.csv'), sep=';', header=[0])
            dist_vali_exc['Hydra-NL']['values'] /= 100 # cm to m
            dist_vali_exc['Hydra-NL met modelonzekerheid'] = pd.read_csv(os.path.join(dir_vali_overschr,'Processed_HydraNL','With_model_uncertainty',f'{stat_name}_with_model_uncertainty.csv'), sep=';', header=[0])
            dist_vali_exc['Hydra-NL met modelonzekerheid']['values'] /= 100 # cm to m
            file_vali_exeed = os.path.join(dir_vali_overschr,'Tables','Exceedance_lines',f'Exceedance_lines_{stat_name}.csv')
            if os.path.exists(file_vali_exeed):
                dist_vali_exc['validation'] = pd.read_csv(file_vali_exeed,sep=';')
                dist_vali_exc['validation']['values'] /= 100
            file_vali_dec = os.path.join(dir_vali_overschr,'Tables','Deceedance_lines',f'Deceedance_lines_{stat_name}.csv')
            if os.path.exists(file_vali_dec):
                dist_vali_dec['validation'] = pd.read_csv(file_vali_dec,sep=';')
                dist_vali_dec['validation']['values'] /= 100
    
        #set station rules
        station_rule_type = 'break'
        if current_station in physical_break_dict.keys(): 
            station_break_value = pd.Timestamp(physical_break_dict[current_station]) #TODO: maybe better to just not select the data by doing data_pd_measext.loc[station_break_value:tstop] instead of data_pd_measext.loc[:tstop]
        else:
            station_break_value = data_pd_measext.index.min()
    
        # 1. Exceedance
        print('Exceedance')
        dist_exc = hatyan.compute_overschrijding(data_pd_HW, rule_type=station_rule_type, rule_value=station_break_value)
        dist_exc.update(dist_vali_exc)
        df_interp = hatyan.interpolate_interested_Tfreqs(dist_exc['Gecombineerd'], Tfreqs=Tfreqs_interested)
        df_interp.to_csv(os.path.join(dir_overschrijding, f'Exceedance_{current_station}.csv'), index=False, sep=';')
        
        fig, ax = hatyan.plot_distributions(dist_exc, name=current_station, color_map='default')
        ax.set_ylim(0,5.5)
        fig.savefig(os.path.join(dir_overschrijding, f'Exceedance_lines_{current_station}.png'))
        
        # 2. Deceedance
        print('Deceedance')
        dist_dec = hatyan.compute_overschrijding(data_pd_LW, rule_type=station_rule_type, rule_value=station_break_value, inverse=True)
        dist_dec.update(dist_vali_dec)
        df_interp = hatyan.interpolate_interested_Tfreqs(dist_dec['Gecombineerd'], Tfreqs=Tfreqs_interested)
        df_interp.to_csv(os.path.join(dir_overschrijding, f'Deceedance_{current_station}.csv'), index=False, sep=';')
        
        fig, ax = hatyan.plot_distributions(dist_dec, name=current_station, color_map='default')
        fig.savefig(os.path.join(dir_overschrijding, f'Deceedance_lines_{current_station}.png'))
    


# #report on memory usage
# print('memory usage')
# var_list = locals().copy() #dir() globals() locals()
# var_keys = var_list.keys()
# max_size = 0
# sum_size = 0
# for varname in var_keys:
#     var_size = sys.getsizeof(var_list[varname])/1024**2 #in MegaBytes
#     max_size = np.maximum(max_size,var_size)
#     sum_size += var_size
#     if var_size > 5:
#         print(f'{varname}: {var_size:.2f}MB')
# print(f'max_size: {max_size:.2f}MB')
# print(f'sum_size: {sum_size:.2f}MB')




