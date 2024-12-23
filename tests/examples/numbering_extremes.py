# -*- coding: utf-8 -*-
"""
numbering_extremes.py
Deze configfile kan gebruikt worden om de dataset data_M2phasediff_perstation.txt bij te werken

"""

import os
import datetime as dt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')
import hatyan
try:
    import dfm_tools as dfmt # pip install dfm_tools
    add_coastlines = True
except ModuleNotFoundError:
    add_coastlines = False

create_spatialplot = True
crs = 28992

#dir_testdata = 'P:\\1209447-kpp-hydraulicaprogrammatuur\\hatyan\\hatyan_data_acceptancetests'
dir_testdata = 'C:\\DATA\\hatyan_data_acceptancetests'

stats_all = ['ABDN','AMLAHVN','BAALHK','BATH','BERGSDSWT','BORSSLE','BOURNMH','BRESKS','BROUWHVSGT02','BROUWHVSGT08','CADZD','CROMR','CUXHVN','DELFZL','DENHDR','DENOVBTN','DEVPT','DORDT','DOVR','EEMHVN','EEMSHVN','EURPFM','EURPHVN','FELSWE','FISHGD','GEULHVN','GOIDSOD','GOUDBG','HAGSBNDN','HANSWT','HARLGN','HARMSBG','HARTBG','HARTHVN','HARVT10','HEESBN','HELLVSS','HOEKVHLD','HOLWD','HUIBGT','IJMDBTHVN','IJMDSMPL','IMMHM','KATSBTN','KEIZVR','KINLBVE','KORNWDZBTN','KRAMMSZWT','KRIMPADIJSL','KRIMPADLK','K13APFM','LAUWOG','LEITH','LICHTELGRE','LITHDP','LLANDNO','LOWST','MAASMSMPL','MAASSS','MAESLKRZZDE','MARLGT','MOERDK','NES','NEWHVN','NEWLN','NIEUWSTZL','NORTHSS','OOSTSDE04','OOSTSDE11','OOSTSDE14','OUDSD','OVLVHWT','PARKSS','PETTZD','PORTSMH','RAKND','ROOMPBNN','ROOMPBTN','ROTTDM','ROZBSSNZDE','ROZBSSZZDE','SCHAARVDND','SCHEVNGN','SCHIERMNOG','SCHOONHVN','SHEERNS','SINTANLHVSGR','SPIJKNSE','STAVNSE','STELLDBTN','STORNWY','SUURHBNZDE','TENNSHVN','TERNZN','TERSLNZE','TEXNZE','VLAARDGN','VLAKTVDRN','VLIELHVN','VLISSGN','VURN','WALSODN','WERKDBTN','WESTKPLE','WESTTSLG','WEYMH','WHITBY','WICK','WIERMGDN','YERSKE','ZALTBML','A12','AUKFPFM','AWGPFM','D15','F16','F3PFM','J6','K14PFM','L9PFM','NORTHCMRT','Q1']
stats_xfac0 = ['A12','ABDN','AUKFPFM','BOURNMH','CROMR','CUXHVN','D15','DEVPT','DOVR','F16','F3PFM','FELSWE','FISHGD','IMMHM','J6','K13APFM','K14PFM','KINLBVE','L9PFM','LEITH','LLANDNO','LOWST','NEWHVN','NEWLN','NORTHCMRT','NORTHSS','PORTSMH','Q1','SHEERNS','STORNWY','WEYMH','WHITBY','WICK']

stats_CADZDm2 = ['WICK']
stats_CADZDm1 = ['ABDN','CROMR','DOVR','EURPFM','FELSWE','IMMHM','LEITH','LOWST','NORTHSS','SHEERNS','WHITBY','VLAKTVDRN']
#stats_CADZD0 = ['CADZD','BATH','VLISSGN','ROOMPBTN','HARVT10','HOEKVHLD','ROTTDM','DORDT','SCHEVNGN','IJMDBTHVN','PETTZD','DENHDR','HARLGN','CUXHVN']
stats_CADZDp1 = []
stats_CADZDtoofar = ['NORTHCMRT','FISHGD','LLANDNO','NEWLN','DEVPT','WEYMH','PORTSMH','BOURNMH','NEWHVN','STORNWY','KINLBVE']

#selected_stations = stats_all
# TODO: temporarily disabled CUXHVN because of https://github.com/Deltares/hatyan/issues/197
selected_stations = ['AUKFPFM']
selected_stations = ['WICK','ABDN','LEITH','WHITBY','IMMHM','CROMR','FELSWE','CADZD','VLISSGN','TERNZN','ROOMPBTN','HARVT10','HOEKVHLD','ROTTDM','DORDT','SCHEVNGN','IJMDBTHVN','PETTZD','DENHDR','DENOVBTN','HARLGN','HOLWD','SCHIERMNOG','LAUWOG','EEMSHVN','DELFZL']#,'CUXHVN']
#selected_stations = ['CROMR','CADZD','HOEKVHLD','DENHDR','CUXHVN']
#selected_stations = ['CADZD','DENHDR']

if 0:
    """
    #all stations in NLD, including several from English coast for which it is know how many tidal waves before CADZD they are. ZALTBML excluded because there is no _ana file available
    stats_xfac1_ana4yr_NAP = [x for x in stats_all if x not in stats_xfac0+stats_MSL+['ZALTBML']]
    selected_stations = stats_xfac1_ana4yr_NAP+stats_CADZDm2+stats_CADZDm1+['CUXHVN']
    selected_stations = list(np.unique(selected_stations))
    selected_stations = ['WICK', 'ABDN', 'LEITH', 'WHITBY', 'IMMHM', 'CROMR', 'FELSWE',
           'VLAKTVDRN', 'CADZD', 'WESTKPLE', 'OOSTSDE11', 'BRESKS', 'VLISSGN',
           'BROUWHVSGT02', 'ROOMPBTN', 'OOSTSDE14', 'BORSSLE', 'OOSTSDE04',
           'BROUWHVSGT08', 'TERNZN', 'HARVT10', 'OVLVHWT', 'STELLDBTN', 'HANSWT',
           'WALSODN', 'MAASMSMPL', 'ROOMPBNN', 'TENNSHVN', 'BAALHK', 'HOEKVHLD',
           'EURPHVN', 'HARTHVN', 'ROZBSSNZDE', 'AMLAHVN', 'SCHAARVDND',
           'MAESLKRZZDE', 'SUURHBNZDE', 'BATH', 'STAVNSE', 'KRAMMSZWT', 'KATSBTN',
           'SINTANLHVSGR', 'YERSKE', 'BERGSDSWT', 'SCHEVNGN', 'MARLGT', 'MAASSS',
           'HARMSBG', 'ROZBSSZZDE', 'GEULHVN', 'VLAARDGN', 'HARTBG', 'SPIJKNSE',
           'EEMHVN', 'PARKSS', 'ROTTDM', 'GOIDSOD', 'IJMDSMPL', 'KRIMPADIJSL',
           'IJMDBTHVN', 'KRIMPADLK', 'GOUDBG', 'DORDT', 'PETTZD', 'SCHOONHVN',
           'DENHDR', 'TEXNZE', 'WERKDBTN', 'HAGSBNDN', 'VURN', 'MOERDK', 'OUDSD',
           'RAKND', 'HELLVSS', 'DENOVBTN', 'TERSLNZE', 'KEIZVR', 'VLIELHVN',
           'WESTTSLG', 'KORNWDZBTN', 'WIERMGDN', 'HARLGN', 'HUIBGT', 'NES',
           'HEESBN', 'LAUWOG', 'SCHIERMNOG', 'HOLWD', 'EEMSHVN', 'LITHDP',
           'DELFZL', 'NIEUWSTZL', 'CUXHVN'] #ordered version of above
    """
    #all stations and ones that are too far away (north, south and southwest of UK) removed
    selected_stations = stats_all
    for stat_remove in stats_CADZDtoofar:
        selected_stations.remove(stat_remove)
    selected_stations.remove('LITHDP') #problem station for at least 2018
    selected_stations.remove('ZALTBML') #no ana file available

yr=2000
yr_HWLWno = 2010
for yr_HWLWno in [2000,2010,2021]: #range(1999,2022):
    stats = pd.DataFrame()
    pd_firstlocalHW_list = pd.DataFrame()

    fig, (ax1,ax2) = plt.subplots(2,1,figsize=(15,8))
    n_colors = len(selected_stations)
    colors = plt.cm.jet(np.linspace(0,1,n_colors))
    for i_stat, current_station in enumerate(selected_stations):
        print('-'*50)
        print('%-45s = %s'%('station_name',current_station))
        print('-'*5)
        
        # station settings
        if current_station in stats_xfac0:
            xfac=False
        else:
            xfac=True
        const_list = hatyan.get_const_list_hatyan('year') #94 const
        
        
        file_data_comp0 = os.path.join(dir_testdata,'predictie2019','%s_ana.txt'%(current_station))
        times_pred = slice(dt.datetime(yr-1,12,31),dt.datetime(yr,1,2,12), "1min")
        file_data_predvali = os.path.join(dir_testdata,'predictie2019','%s_pre.txt'%(current_station))
    
        #component groups
        COMP_merged = hatyan.read_components(filename=file_data_comp0)
        
        #prediction and validation
        ts_prediction = hatyan.prediction(comp=COMP_merged, times=times_pred)

        ts_ext_prediction = hatyan.calc_HWLW(ts=ts_prediction)
        
        if i_stat == 0:
            if 'CADZD' in stats_xfac0:
                xfac_cadzd=False
            else:
                xfac_cadzd=True
            COMP_merged_CADZD = hatyan.read_components(filename=file_data_comp0.replace(current_station,'CADZD'))
            COMP_merged_CADZD_M2 = COMP_merged_CADZD.loc[['A0','M2']]
            ts_prediction_CADZD = hatyan.prediction(comp=COMP_merged_CADZD, times=times_pred)
            ts_prediction_CADZD_M2 = hatyan.prediction(comp=COMP_merged_CADZD_M2, times=times_pred)
            ax1.plot(ts_prediction_CADZD_M2.index, ts_prediction_CADZD_M2['values'], label='CADZD_M2', color='k')
            ts_ext_prediction_CADZD = hatyan.calc_HWLW(ts=ts_prediction_CADZD)
            ts_ext_prediction_CADZD_newyr = ts_ext_prediction_CADZD.loc[str(yr):]
            bool_newyr_hw = ts_ext_prediction_CADZD_newyr['HWLWcode']==1
            firstHWcadz = ts_ext_prediction_CADZD_newyr.loc[bool_newyr_hw].index[0]
        
        M2phasediff_raw = (COMP_merged.loc['M2','phi_deg']-COMP_merged_CADZD.loc['M2','phi_deg'])%360
        
        if current_station in stats_CADZDp1:
            bool_wrtHWcadz = (ts_ext_prediction.index >= firstHWcadz) & (ts_ext_prediction['HWLWcode']==1)
            ts_firstlocalHW = ts_ext_prediction.loc[bool_wrtHWcadz].iloc[1] #second HW after HWcadzd
            M2phasediff = M2phasediff_raw+360
        elif current_station in stats_CADZDm1:
            bool_wrtHWcadz = (ts_ext_prediction.index <= firstHWcadz) & (ts_ext_prediction['HWLWcode']==1)
            ts_firstlocalHW = ts_ext_prediction.loc[bool_wrtHWcadz].iloc[-1] #last HW before HWcadzd
            M2phasediff = M2phasediff_raw-360
        elif current_station in stats_CADZDm2:
            bool_wrtHWcadz = (ts_ext_prediction.index <= firstHWcadz) & (ts_ext_prediction['HWLWcode']==1)
            ts_firstlocalHW = ts_ext_prediction.loc[bool_wrtHWcadz].iloc[-2] #second to last HW before HWcadzd
            M2phasediff = M2phasediff_raw-360-360
        else:
            bool_wrtHWcadz = (ts_ext_prediction.index >= firstHWcadz) & (ts_ext_prediction['HWLWcode']==1)
            ts_firstlocalHW = ts_ext_prediction.loc[bool_wrtHWcadz].iloc[0] #first HW after HWcadzd
            M2phasediff = M2phasediff_raw
            
        #print('tdiff %s:'%(current_station), ts_firstlocalHW.index-firstHWcadz)
        pdrow = pd.DataFrame({'time': [ts_firstlocalHW.name], 'HWtdiff_hr': [(ts_firstlocalHW.name-firstHWcadz).total_seconds()/3600], 'M2phase':COMP_merged.loc['M2','phi_deg'], 'M2phasediff':M2phasediff}, index=[current_station])
        if create_spatialplot:
            RDx, RDy = hatyan.get_diaxycoords(filename=file_data_predvali, crs=crs)
            pdrow['RDx'] = RDx
            pdrow['RDy'] = RDy
        stats = pd.concat([stats,pdrow])

        ax1.plot(ts_prediction.index, ts_prediction['values'], label=current_station, color=colors[i_stat])
        ax1.plot([ts_firstlocalHW.name,ts_firstlocalHW.name], [ts_firstlocalHW['values'], 2.5], '--', linewidth=1.5, color=colors[i_stat])
        ax1.plot(ts_firstlocalHW.name,ts_firstlocalHW['values'],'x', color=colors[i_stat])
        
        if 1: #validation case
            #calculate tidal wave number
            times_pred_HWLWno = slice(dt.datetime(yr_HWLWno-1,12,31),dt.datetime(yr_HWLWno,1,2,12), "1min")
            ts_prediction_HWLWno = hatyan.prediction(comp=COMP_merged, times=times_pred_HWLWno)
            ts_ext_prediction_HWLWno_pre = hatyan.calc_HWLW(ts=ts_prediction_HWLWno)
            
            print(current_station)
            ts_ext_prediction_HWLWno = hatyan.calc_HWLWnumbering(ts_ext=ts_ext_prediction_HWLWno_pre, station=current_station)
            
            print(ts_ext_prediction_HWLWno)
            for irow, pdrow in ts_ext_prediction_HWLWno.iterrows():
                ax2.text(pdrow.name,pdrow['values'],pdrow['HWLWno'].astype(int), color=colors[i_stat])
            HWLWno_focus = int(np.round((dt.datetime(yr_HWLWno,1,1,9,45)-dt.datetime(2000,1,1,9,45)).total_seconds()/3600/12.420601))
            ts_firstlocalHW_fromcalc = ts_ext_prediction_HWLWno[(ts_ext_prediction_HWLWno['HWLWcode']==1) & (ts_ext_prediction_HWLWno['HWLWno']==HWLWno_focus)]
            
            pd_firstlocalHW_list = pd.concat([pd_firstlocalHW_list,ts_firstlocalHW_fromcalc])
            
            ax2.plot(ts_prediction_HWLWno.index, ts_prediction_HWLWno['values'], label=current_station, color=colors[i_stat])
            ax2.plot([ts_firstlocalHW_fromcalc.index[0],ts_firstlocalHW_fromcalc.index[0]], [ts_firstlocalHW_fromcalc['values'].iloc[0], 2.5], '--', linewidth=1.5, color=colors[i_stat])
            ax2.plot(ts_firstlocalHW_fromcalc.index,ts_firstlocalHW_fromcalc['values'],'x', color=colors[i_stat])
    
    ax2.plot(pd_firstlocalHW_list.index,pd_firstlocalHW_list['values'],'-ok')
    
    stats['M2phasediff_hr'] = stats['M2phasediff']/360*12.420601
    stats_M2phasediff_out = stats.sort_values('M2phasediff_hr')['M2phasediff']

    print(stats)
    print('')
    print(hatyan.schureman.get_schureman_freqs(['M2']))
    
    ax1.set_xlim(times_pred.start,times_pred.stop)
    ax2.set_xlim(times_pred_HWLWno.start,times_pred_HWLWno.stop)
    ax2.set_xlim([dt.datetime(yr_HWLWno-1,12,31),dt.datetime(yr_HWLWno,1,2,12)])
    ax1_ylim = ax1.get_ylim()
    ax1.plot([dt.datetime(yr,1,1),dt.datetime(yr,1,1)],ax1_ylim,'k--')
    fig.tight_layout()
    for ax in (ax1,ax2):
        ax.legend(loc=2, fontsize=7)#bbox_to_anchor=(1,1))
        ax.set_ylim(ax1_ylim)
        import matplotlib.dates as mdates
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d %H:%M"))
    fig.savefig('tide_numbering_%i.png'%(yr_HWLWno), dpi=250)

if create_spatialplot:
    fig2, (fig2_ax1) = plt.subplots(1,1,figsize=(10,9))
    fig2_ax1.plot(stats.loc['CADZD','RDx'], stats.loc['CADZD','RDy'],'xk')
    pc = fig2_ax1.scatter(stats['RDx'], stats['RDy'],10,stats['M2phasediff'], vmin=-360,vmax=360, cmap='hsv')
    fig2.colorbar(pc, ax=fig2_ax1)
    fig2_ax1.set_xlim(-400000,360000)
    fig2_ax1.set_ylim(200000,1200000)
    if add_coastlines:
        dfmt.plot_coastlines(ax=fig2_ax1, res='i', min_area=50, linewidth=0.5, crs=crs)
        dfmt.plot_borders(ax=fig2_ax1, crs=crs)
        
    for irow, pdrow in stats.iterrows():
        fig2_ax1.text(pdrow['RDx'], pdrow['RDy'], '%.1f'%(pdrow['M2phasediff']))
    fig2.tight_layout()
    fig2.savefig('tide_numbering_phasediff.png', dpi=250)
