# -*- coding: utf-8 -*-
"""
predictie_2019_frommergedcomp.py
hatyan master configfile
voor alle stations indien mogelijk:
    - inlezen analyseresultatenbestand
    - predictie maken

"""

import os, sys
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
        ax1.set_title(f'fullset vs M2+S2+M4 {current_station}')
        ax1.plot(ts_prediction,linewidth=1,label='full set')
        ax1.plot(ts_prediction_20,linewidth=1,label='20 components')
        ax1.legend(loc=1)
        ax1.grid()
        ax1.set_xlim(times_ext_twoweeks)
        ax1.set_ylim(-1.1,1.6)
        fig.tight_layout()
        fig.savefig(f'fullset_vs_20comp_{current_station}')


########################
# prediction from comp, tried to show LAT
if 0:
    #dir_testdata = 'P:\\1209447-kpp-hydraulicaprogrammatuur\\hatyan\\hatyan_data_acceptancetests'
    dir_testdata = 'C:\\DATA\\hatyan_data_acceptancetests'
    
    stats_all = ['ABDN','AMLAHVN','BAALHK','BATH','BERGSDSWT','BORSSLE','BOURNMH','BRESKS','BROUWHVSGT02','BROUWHVSGT08','CADZD','CROMR','CUXHVN','DELFZL','DENHDR','DENOVBTN','DEVPT','DORDT','DOVR','EEMHVN','EEMSHVN','EURPFM','EURPHVN','FELSWE','FISHGD','GEULHVN','GOIDSOD','GOUDBG','HAGSBNDN','HANSWT','HARLGN','HARMSBG','HARTBG','HARTHVN','HARVT10','HEESBN','HELLVSS','HOEKVHLD','HOLWD','HUIBGT','IJMDBTHVN','IJMDSMPL','IMMHM','KATSBTN','KEIZVR','KINLBVE','KORNWDZBTN','KRAMMSZWT','KRIMPADIJSL','KRIMPADLK','K13APFM','LAUWOG','LEITH','LICHTELGRE','LITHDP','LLANDNO','LOWST','MAASMSMPL','MAASSS','MAESLKRZZDE','MARLGT','MOERDK','NES','NEWHVN','NEWLN','NIEUWSTZL','NORTHSS','OOSTSDE04','OOSTSDE11','OOSTSDE14','OUDSD','OVLVHWT','PARKSS','PETTZD','PORTSMH','RAKND','ROOMPBNN','ROOMPBTN','ROTTDM','ROZBSSNZDE','ROZBSSZZDE','SCHAARVDND','SCHEVNGN','SCHIERMNOG','SCHOONHVN','SHEERNS','SINTANLHVSGR','SPIJKNSE','STAVNSE','STELLDBTN','STORNWY','SUURHBNZDE','TENNSHVN','TERNZN','TERSLNZE','TEXNZE','VLAARDGN','VLAKTVDRN','VLIELHVN','VLISSGN','VURN','WALSODN','WERKDBTN','WESTKPLE','WESTTSLG','WEYMH','WHITBY','WICK','WIERMGDN','YERSKE','ZALTBML','A12','AUKFPFM','AWGPFM','D15','F16','F3PFM','J6','K14PFM','L9PFM','NORTHCMRT','Q1']
    stats_xfac0 = ['A12','ABDN','AUKFPFM','BOURNMH','CROMR','CUXHVN','D15','DEVPT','DOVR','F16','F3PFM','FELSWE','FISHGD','IMMHM','J6','K13APFM','K14PFM','KINLBVE','LEITH','LLANDNO','LOWST','NEWHVN','NEWLN','NORTHCMRT','NORTHSS','PORTSMH','SHEERNS','STORNWY','WEYMH','WHITBY','WICK']
    selected_stations = []
    
    current_station = 'HOEKVHLD'
    nodalfactors = True
    if current_station in stats_xfac0:
        xfac=False
    else:
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
    data_pkl = pd.read_pickle(r'p:\11208031-010-kenmerkende-waarden-k\work\measurements_wl_18700101_20220101\DELFZL_measwl.pkl')
    ts_meas = data_pkl[['values']]
    ts_meas.index = ts_meas.index.tz_localize(None)
    ts_meas = hatyan.crop_timeseries(ts_meas,times_ext=[dt.datetime(2001,1,1),dt.datetime(2019,12,31,23,50)])#,onlyfull=False)
    
    comp_avg, comp_allperiods = hatyan.get_components_from_ts(ts_meas,const_list='year',analysis_perperiod='Y',return_allperiods=True)
    #comp_allyears = comp_allyears.drop('A0')
    ts_pred_py = hatyan.prediction_perperiod(comp_allperiods, timestep_min=10)
    
    A0_allyears_meas = ts_meas['values'].groupby(by=pd.Grouper(freq='Y')).mean()
    A0_allyears_pred = ts_pred_py['values'].groupby(by=pd.Grouper(freq='Y')).mean()
    min_allyears_meas = ts_meas['values'].groupby(by=pd.Grouper(freq='Y')).min()
    min_allyears_pred = ts_pred_py['values'].groupby(by=pd.Grouper(freq='Y')).min()
    max_allyears_meas = ts_meas['values'].groupby(by=pd.Grouper(freq='Y')).max()
    max_allyears_pred = ts_pred_py['values'].groupby(by=pd.Grouper(freq='Y')).max()
    
    fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_meas, ts_validation=None)
    ax1.plot(ts_pred_py)
    #ax2.plot(A0_allyears_meas,label='A0_allyears_meas')
    ax2.plot(A0_allyears_pred,label='A0_allyears_pred')
    #ax2.plot(min_allyears_meas,label='min_allyears_meas')
    #ax2.plot(min_allyears_pred,label='min_allyears_pred')
    #ax2.plot(max_allyears_meas,label='max_allyears_meas')
    #ax2.plot(max_allyears_pred,label='max_allyears_pred')
    ax2.set_ylim(A0_allyears_pred.min()-0.02,A0_allyears_pred.max()+0.02)
    ax2.legend()




