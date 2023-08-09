# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 16:35:26 2022

@author: veenstra
"""

import numpy as np
import pandas as pd
from hatyan.hatyan_core import get_const_list_hatyan
from hatyan.analysis_prediction import analysis
from hatyan import HatyanSettings

def get_gemgetij_components(data_pd_meas):
    # =============================================================================
    # Hatyan analyse voor 10 jaar (alle componenten voor gemiddelde getijcyclus) #TODO: maybe use original 4y period/componentfile instead? SA/SM should come from 19y analysis
    # =============================================================================
    const_list = get_const_list_hatyan('year') #components should not be reduced, since higher harmonics are necessary
    hatyan_settings_ana = HatyanSettings(nodalfactors=True, fu_alltimes=False, xfac=True, analysis_perperiod='Y', return_allperiods=True) #RWS-default settings
    comp_frommeasurements_avg, comp_frommeasurements_allyears = analysis(data_pd_meas, const_list=const_list, hatyan_settings=hatyan_settings_ana)
    
    # #check if all years are available
    # comp_years = comp_frommeasurements_allyears['A'].columns
    # expected_years = tstop_dt.year-tstart_dt.year
    # if len(comp_years) < expected_years:
    #     raise Exception('ERROR: analysis result contains not all years')
    
    #check if nans in analysis
    if comp_frommeasurements_avg.isnull()['A'].any():
        raise Exception('ERROR: analysis result contains nan values')
    
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
        comp_av.loc[comp_higherharmonics,'A'] = np.sqrt((comp_iM['A']**2).sum()) #kwadraatsom
    
    comp_av.loc['A0'] = comp_frommeasurements_avg.loc['A0']
    
    print('verhouding tussen originele en kwadratensom componenten')
    print(comp_av/comp_frommeasurements_avg.loc[components_av]) # values are different than 1991.0 document and differs per station while the document states "Zoals te verwachten is de verhouding per component tussen deze wortel en de oorspronkelijke amplitude voor alle plaatsen gelijk"

    return comp_frommeasurements_avg, comp_av


def reshape_signal(ts, ts_ext, HW_goal, LW_goal, tP_goal=None):
    """
    scales tidal signal to provided HW/LW value and up/down going time
    tP_goal (tidal period time) is used to fix tidalperiod to 12h25m (for BOI timeseries)
    
    time_down was scaled with havengetallen before, but not anymore to avoid issues with aggers
    """
    TR_goal = HW_goal-LW_goal
    
    #selecteer alle hoogwaters en opvolgende laagwaters
    bool_HW = ts_ext['HWLWcode']==1
    idx_HW = np.where(bool_HW)[0]
    timesHW = ts_ext.index[idx_HW]
    timesLW = ts_ext.index[idx_HW[:-1]+1] #assuming alternating 1,2,1 or 1,3,1, this is always valid in this workflow
    
    #crop from first to last HW (rest is not scaled anyway)
    ts_time_firstHW = ts_ext[bool_HW].index[0]
    ts_time_lastHW = ts_ext[bool_HW].index[-1]
    ts_corr = ts.copy().loc[ts_time_firstHW:ts_time_lastHW]

    ts_corr['times'] = ts_corr.index #this is necessary since datetimeindex with freq is not editable, and Series is editable
    for i in np.arange(0,len(timesHW)-1):
        HW1_val = ts_corr.loc[timesHW[i],'values']
        HW2_val = ts_corr.loc[timesHW[i+1],'values']
        LW_val = ts_corr.loc[timesLW[i],'values']
        TR1_val = HW1_val-LW_val
        TR2_val = HW2_val-LW_val
        tP_val = timesHW[i+1]-timesHW[i]
        if tP_goal is None:
            tP_goal = tP_val
        
        temp1 = (ts_corr.loc[timesHW[i]:timesLW[i],'values']-LW_val)/TR1_val*TR_goal+LW_goal
        temp2 = (ts_corr.loc[timesLW[i]:timesHW[i+1],'values']-LW_val)/TR2_val*TR_goal+LW_goal
        temp = pd.concat([temp1,temp2.iloc[1:]]) #.iloc[1:] since timesLW[i] is in both timeseries (values are equal)
        ts_corr['values_new'] = temp
        
        tide_HWtoHW = ts_corr.loc[timesHW[i]:timesHW[i+1]]
        ts_corr['times'] = pd.date_range(start=ts_corr.loc[timesHW[i],'times'],end=ts_corr.loc[timesHW[i],'times']+tP_goal,periods=len(tide_HWtoHW))
        
    ts_corr = ts_corr.set_index('times',drop=True)
    ts_corr['values'] = ts_corr['values_new']
    ts_corr = ts_corr.drop(['values_new'],axis=1)
    return ts_corr


def ts_to_trefHW(ts,HWreftime=None):
    """
    converts to hours relative to HWreftime, to plot av/sp/np tidal signals in one plot
    """
    ts.index.name = 'times' #just to be sure
    ts_trefHW = ts.reset_index()
    if HWreftime is None:
        ts_trefHW.index = (ts_trefHW['times']-ts_trefHW['times'].iloc[0]).dt.total_seconds()/3600
    else:
        ts_trefHW.index = (ts_trefHW['times']-HWreftime).dt.total_seconds()/3600
    return ts_trefHW


def repeat_signal(ts_one_HWtoHW, nb, na):
    """
    repeat tidal signal, necessary for sp/np, since they are computed as single tidal signal first
    """
    tidalperiod = ts_one_HWtoHW.index[-1] - ts_one_HWtoHW.index[0]
    ts_rep = pd.DataFrame()
    for iAdd in np.arange(-nb,na+1):
        ts_add = pd.DataFrame({'values':ts_one_HWtoHW['values'].values},
                              index=ts_one_HWtoHW.index + iAdd*tidalperiod)
        ts_rep = pd.concat([ts_rep,ts_add])
    ts_rep = ts_rep.loc[~ts_rep.index.duplicated()]
    return ts_rep
