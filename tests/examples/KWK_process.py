# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 14:17:13 2022

@author: veenstra
"""

import os
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt
plt.close('all')
import hatyan # available via `pip install hatyan` or at https://github.com/Deltares/hatyan (pip will not work since not all KWK functions are already in release)

dataTKdia = True #TODO: communicate data issues to TK (wl and ext): p:\11208031-010-kenmerkende-waarden-k\work\data_vanRWS_20220805\convert_dia2pickle_dataTK.py
NAP2005correction = False #True #TODO: define for all stations

tstart_dt = dt.datetime(2011,1,1)
tstop_dt = dt.datetime(2021,1,1)
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
else:
    stat_list = None #TODO: get from DDL catalog (check KWK_download script)


def clean_data(ts_meas_pd,current_station):
    if 'HWLWcode' in ts_meas_pd.columns:
        keep_columns = ['values','QC','HWLWcode']
    else:
        keep_columns = ['values','QC']
    ts_meas_pd = ts_meas_pd[keep_columns] # reduces the memory consumption significantly in case of DDL data with a lot of metadata
    ts_meas_pd.index = ts_meas_pd.index.tz_localize(None)
    ts_meas_pd = ts_meas_pd.loc[~(ts_meas_pd['QC']==99)] #remove invalid data
    
    #optional nap correction
    if NAP2005correction:
        ts_meas_pd = nap2005_correction(ts_meas_pd,current_station)
    return ts_meas_pd


def nap2005_correction(data_pd,current_station):
    #NAP correction for dates before 1-1-2005
    #TODO: check if ths make a difference (for havengetallen it makes a slight difference so yes. For gemgetijkromme it only makes a difference for spring/doodtij. (now only applied at gemgetij en havengetallen)). If so, make this flexible per station, where to get the data or is the RWS data already corrected for it?
    #herdefinitie van NAP (~20mm voor HvH in fig2, relevant?): https://puc.overheid.nl/PUC/Handlers/DownloadDocument.ashx?identifier=PUC_113484_31&versienummer=1
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



#physical_break_dict for slotgemiddelden and overschrijdingsfrequenties TODO: maybe use everywhere to crop data?
physical_break_dict = {'DENOVBTN':'1933', #laatste sluitgat afsluitdijk in 1932 
                       'HARLGN':'1933', #laatste sluitgat afsluitdijk in 1932
                       'VLIELHVN':'1933', #laatste sluitgat afsluitdijk in 1932
                       } #TODO: add physical_break for STAVNSE and KATSBTN? (Oosterscheldekering)

compute_slotgem = True
compute_havengetallen = True
compute_gemgetij = True
compute_overschrijding = True


for current_station in ['HOEKVHLD','DENOVBTN']:#stat_list[stat_list.index('SCHEVNGN'):]:#
    plt.close('all')
    
    print(f'loading data for {current_station}')
    file_wl_pkl = os.path.join(dir_meas,f"{current_station}_measwl.pkl")
    if os.path.exists(file_wl_pkl): #for slotgemiddelden, gemgetijkrommen (needs slotgem+havget)
        data_pd_meas_all = pd.read_pickle(file_wl_pkl)
        data_pd_meas_all = clean_data(data_pd_meas_all,current_station)
        #crop measurement data
        data_pd_meas_10y = hatyan.crop_timeseries(data_pd_meas_all, times_ext=[tstart_dt,tstop_dt-dt.timedelta(minutes=10)])#,onlyfull=False)
    
    file_ext_pkl = os.path.join(dir_meas,f"{current_station}_measext.pkl")
    if os.path.exists(file_ext_pkl): #for slotgemiddelden, havengetallen, overschrijding
        data_pd_HWLW_all = pd.read_pickle(file_ext_pkl)
        data_pd_HWLW_all = clean_data(data_pd_HWLW_all,current_station)
        if compute_slotgem or compute_havengetallen or compute_overschrijding: #TODO: make calc_HWLW12345to12() faster
            data_pd_HWLW_all_12 = hatyan.calc_HWLW12345to12(data_pd_HWLW_all) #convert 12345 to 12 by taking minimum of 345 as 2 (laagste laagwater)
            #crop timeseries to 10y
            data_pd_HWLW_10y_12 = hatyan.crop_timeseries(data_pd_HWLW_all_12, times_ext=[tstart_dt,tstop_dt],onlyfull=False)
            
            #check if amount of HWs is enough
            M2_period_timedelta = pd.Timedelta(hours=hatyan.get_schureman_freqs(['M2']).loc['M2','period [hr]'])
            numHWs_expected = (tstop_dt-tstart_dt).total_seconds()/M2_period_timedelta.total_seconds()
            numHWs = (data_pd_HWLW_10y_12['HWLWcode']==1).sum()
            if numHWs < 0.95*numHWs_expected:
                raise Exception(f'ERROR: not enough high waters present in period, {numHWs} instead of >=0.95*{int(numHWs_expected):d}')
    
    
    
    
    
    #### SLOTGEMIDDELDEN
    #TODO: nodal cycle is not in same phase for all stations, this is not physically correct.
    #TODO: more data is needed for proper working of fitting for some stations (2011: BAALHK, BRESKVHVN, GATVBSLE, SCHAARVDND)
    if compute_slotgem and os.path.exists(file_wl_pkl):
        print(f'slotgemiddelden for {current_station}')
        
        #calculate yearly mean
        dict_wltidalindicators = hatyan.calc_wltidalindicators(data_pd_meas_all)
        wl_mean_peryear = dict_wltidalindicators['wl_mean_peryear']
        dict_wltidalindicators_valid = hatyan.calc_wltidalindicators(data_pd_meas_all, tresh_yearlywlcount=2900) #24*365=8760 (hourly interval), 24/3*365=2920 (3-hourly interval)
        wl_mean_peryear_valid = dict_wltidalindicators_valid['wl_mean_peryear']
        
        #derive tidal indicators like yearmean HWLW from HWLW values
        if os.path.exists(file_ext_pkl):
            dict_HWLWtidalindicators = hatyan.calc_HWLWtidalindicators(data_pd_HWLW_all_12)
            HW_mean_peryear = dict_HWLWtidalindicators['HW_mean_peryear']
            LW_mean_peryear = dict_HWLWtidalindicators['LW_mean_peryear']
            dict_HWLWtidalindicators_valid = hatyan.calc_HWLWtidalindicators(data_pd_HWLW_all_12, tresh_yearlyHWLWcount=1400) #2*24*365/12.42=1410.6 (12.42 hourly extreme)
            HW_mean_peryear_valid = dict_HWLWtidalindicators_valid['HW_mean_peryear']
            LW_mean_peryear_valid = dict_HWLWtidalindicators_valid['LW_mean_peryear']
        
        #plotting (yearly averages are plotted on 1jan)
        fig,ax1 = plt.subplots(figsize=(12,6))
        
        #get and plot validation timeseries (yearly mean wl/HW/LW)
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
            ax1.plot(yearmeanwl['values'],'+g',label='yearmean validation')
        
        #plot values
        if os.path.exists(file_ext_pkl):
            ax1.plot(HW_mean_peryear,'x',color='grey')
            ax1.plot(LW_mean_peryear,'x',color='grey')
            ax1.plot(HW_mean_peryear_valid,'xr')
            ax1.plot(LW_mean_peryear_valid,'xr')
        ax1.plot(wl_mean_peryear,'x',color='grey',label='yearmean')
        ax1.plot(wl_mean_peryear_valid,'xr',label='yearmean significant')
        ax1.grid()
        ax1.set_xlim(fig_alltimes_ext) # entire period
        ax1.set_ylabel('waterstand [m]')
        ax1.set_title(f'yearly mean HW/wl/LW {current_station}')
        fig.tight_layout()
        
        if current_station in physical_break_dict.keys():
            tstart_dt_trend = physical_break_dict[current_station]
        else:
            tstart_dt_trend = None
        
        #fit linear models over yearly mean values
        wl_mean_array_todate = wl_mean_peryear_valid.loc[tstart_dt_trend:tstop_dt] #remove all values after tstop_dt (is year_slotgem)
        pred_pd_wl = hatyan.fit_models(wl_mean_array_todate)
        ax1.plot(pred_pd_wl, ".-", label=pred_pd_wl.columns)
        ax1.set_prop_cycle(None) #reset matplotlib colors
        #2021.0 value (and future)
        ax1.plot(pred_pd_wl.loc[tstop_dt:,'pred_linear_winodal'], ".k", label=f'pred_linear from {year_slotgem}')
        pred_slotgem = pred_pd_wl.loc[[tstop_dt],['pred_linear_winodal']]
        pred_slotgem.to_csv(os.path.join(dir_slotgem,f'slotgem_value_{current_station}.txt'))
        ax1.legend(loc=2)
        
        if os.path.exists(file_ext_pkl):
            HW_mean_array_todate = HW_mean_peryear_valid.loc[tstart_dt_trend:tstop_dt] #remove all values after tstop_dt (is year_slotgem)
            pred_pd_HW = hatyan.fit_models(HW_mean_array_todate)
            ax1.plot(pred_pd_HW, ".-", label=pred_pd_HW.columns)
            ax1.set_prop_cycle(None) #reset matplotlib colors
            
            LW_mean_array_todate = LW_mean_peryear_valid.loc[tstart_dt_trend:tstop_dt] #remove all values after tstop_dt (is year_slotgem)
            pred_pd_LW = hatyan.fit_models(LW_mean_array_todate)
            ax1.plot(pred_pd_LW, ".-", label=pred_pd_LW.columns)
            ax1.set_prop_cycle(None) #reset matplotlib colors

        fig.savefig(os.path.join(dir_slotgem,f'yearly_values_{current_station}'))
    
    
    
    
    ### HAVENGETALLEN 
    if compute_havengetallen and os.path.exists(file_ext_pkl):
        
        print(f'havengetallen for {current_station}')
        
        #TODO: move calc_HWLW_moonculm_combi() to top since it is the same for all stations
        culm_addtime = 4*dt.timedelta(hours=12,minutes=25)+dt.timedelta(hours=1)-dt.timedelta(minutes=20) # 2d and 2u20min correction, this shifts the x-axis of aardappelgrafiek: HW is 2 days after culmination (so 4x25min difference between length of avg moonculm and length of 2 days), 1 hour (GMT to MET), 20 minutes (0 to 5 meridian, is commented now)
        #TODO: check culm_addtime and HWLWno+4 offsets. culm_addtime could also be 2 days or 2days +1h GMT-MET correction. 20 minutes seems odd since moonculm is about tidal wave from ocean
        data_pd_HWLW = hatyan.calc_HWLW_moonculm_combi(data_pd_HWLW_12=data_pd_HWLW_10y_12, culm_addtime=culm_addtime) #culm_addtime=None provides the same gemgetijkromme now delay is not used for scaling anymore
        HWLW_culmhr_summary = hatyan.calc_HWLW_culmhr_summary(data_pd_HWLW)
        
        print('HWLW FIGUREN PER TIJDSKLASSE, INCLUSIEF MEDIAN LINE')
        fig, axs = hatyan.plot_HWLW_pertimeclass(data_pd_HWLW, HWLW_culmhr_summary)
        for ax in axs.ravel():
            ax.set_title(f'{ax.get_title()} {current_station} {year_slotgem}')
        fig.savefig(os.path.join(dir_havget,f'HWLW_pertijdsklasse_inclmedianline_{current_station}'))
        
        print('AARDAPPELGRAFIEK')
        fig, (ax1,ax2) = hatyan.plot_aardappelgrafiek(HWLW_culmhr_summary)
        for ax in (ax1,ax2):
            ax.set_title(f'{ax.get_title()} {current_station} {year_slotgem}')
        fig.savefig(os.path.join(dir_havget, f'aardappelgrafiek_{year_slotgem}_{current_station}'))
        
        #write to csv
        HWLW_culmhr_summary_exp = HWLW_culmhr_summary.loc[[6,'mean',0]] #select neap/mean/springtide
        HWLW_culmhr_summary_exp.index = ['neap','mean','spring']
        HWLW_culmhr_summary_exp.to_csv(os.path.join(dir_havget, f'havengetallen_{year_slotgem}_{current_station}.csv'),float_format='%.3f')
    
    
    
    
    
    ##### GEMIDDELDE GETIJKROMMEN
    if compute_gemgetij and os.path.exists(file_wl_pkl):
        
        print(f'gem getijkrommen for {current_station}')
        pred_freq_sec = 10 #TODO: frequency decides accuracy of tU/tD and other timings (and is writing freq of BOI timeseries)
        
        #TODO: add correctie havengetallen HW/LW av/sp/np met slotgemiddelde uit PLSS/modelfit (HW/LW av)
        file_havget = os.path.join(dir_havget,f'havengetallen_{year_slotgem}_{current_station}.csv')
        if not os.path.exists(file_havget):
            raise Exception(f'havengetallen file does not exist: {file_havget}')
        data_havget = pd.read_csv(file_havget,index_col=0,parse_dates=['HW_delay_median','LW_delay_median','getijperiod_median','duurdaling_median'])
        HW_sp, LW_sp = data_havget.loc['spring',['HW_values_median','LW_values_median']]
        HW_np, LW_np = data_havget.loc['neap',['HW_values_median','LW_values_median']]
        HW_av, LW_av = data_havget.loc['mean',['HW_values_median','LW_values_median']]
        
        #derive components via TA on measured waterlevels
        comp_frommeasurements_avg, comp_av = hatyan.get_gemgetij_components(data_pd_meas_10y)
        
        times_pred_1mnth = pd.date_range(start=dt.datetime(tstop_dt.year,1,1,0,0)-dt.timedelta(hours=12), end=dt.datetime(tstop_dt.year,2,1,0,0), freq=f'{pred_freq_sec} S') #start 12 hours in advance, to assure also corrected values on desired tstart
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
        
        
        #timeseries for gele boekje (av/sp/np have different lengths, time is relative to HW of av and HW of sp/np are shifted there) #TODO: is this product still necessary?
        print(f'reshape_signal GEMGETIJ: {current_station}')
        prediction_av_one_trefHW = hatyan.ts_to_trefHW(prediction_av_one) # repeating one is not necessary for av, but easier to do the same for av/sp/np
        prediction_av_corr_one = hatyan.reshape_signal(prediction_av_one, prediction_av_ext_one, HW_goal=HW_av, LW_goal=LW_av, tP_goal=None)
        prediction_av_corr_rep5 = hatyan.repeat_signal(prediction_av_corr_one, nb=2, na=2)
        prediction_av_corr_rep5_trefHW = hatyan.ts_to_trefHW(prediction_av_corr_rep5,HWreftime=ia1)
        
        print(f'reshape_signal SPRINGTIJ: {current_station}')
        prediction_sp_one_trefHW = hatyan.ts_to_trefHW(prediction_sp_one)
        prediction_sp_corr_one = hatyan.reshape_signal(prediction_sp_one, prediction_sp_ext_one, HW_goal=HW_sp, LW_goal=LW_sp, tP_goal=None)
        prediction_sp_corr_rep5 = hatyan.repeat_signal(prediction_sp_corr_one, nb=2, na=2)
        prediction_sp_corr_rep5_trefHW = hatyan.ts_to_trefHW(prediction_sp_corr_rep5,HWreftime=is1)
        
        print(f'reshape_signal DOODTIJ: {current_station}')
        prediction_np_one_trefHW = hatyan.ts_to_trefHW(prediction_np_one)
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
        ax_sum.plot(prediction_av_one_trefHW['values'],'--', color=cmap(0), linewidth=0.7, label='gem kromme, one')
        ax_sum.plot(prediction_av_corr_rep5_trefHW['values'], color=cmap(0), label='gem kromme, corr')
        ax_sum.plot(prediction_sp_one_trefHW['values'],'--', color=cmap(1), linewidth=0.7, label='sp kromme, one')
        ax_sum.plot(prediction_sp_corr_rep5_trefHW['values'], color=cmap(1), label='sp kromme, corr')
        ax_sum.plot(prediction_np_one_trefHW['values'],'--', color=cmap(2), linewidth=0.7, label='np kromme, one')
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
    #TODO: SLR trend correctie voor overschrijdingsfrequenties en evt ook voor andere KW?
    #plots beoordelen: rode lijn moet ongeveer verlengde zijn van groene, als die ineens omhoog piekt komt dat door hele extreme waardes die je dan vermoedelijk ook al ziet in je groene lijn
    
    Tfreqs_interested = [5, 2, 1, 1/2, 1/5, 1/10, 1/20, 1/50, 1/100, 1/200, #overschrijdingsfreqs
                         1/500, 1/1000, 1/2000, 1/4000, 1/5000, 1/10000] #TODO: which frequencies are realistic with n years of data? probably remove this entire row >> met 40 jaar data kun je in principe tot 1/40 gaan, maar met weibull kun je extrapoleren en in theorie >> dit is voor tabel die je eruit wil hebben
    
    if compute_overschrijding and os.path.exists(file_ext_pkl):
    
        print(f'overschrijdingsfrequenties for {current_station}')
        
        #clip data #TODO: do at top?
        data_pd_measext = data_pd_HWLW_all_12.loc[:tstop_dt] # only include data up to year_slotgem
        
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
            station_break_value = physical_break_dict[current_station] #TODO: maybe better to just not select the data by doing data_pd_measext.loc[station_break_value:tstop] instead of data_pd_measext.loc[:tstop]
        else:
            station_break_value = data_pd_measext.index.min()
    
        # 1. Exceedance
        print('Exceedance') #TODO: hatyan.get_weibull.der_pfunc() throws "RuntimeWarning: invalid value encountered in double_scalars"
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



