# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 17:12:42 2022

@author: veenstra
"""

import pandas as pd
from hatyan.timeseries import calc_HWLW12345to21


def calc_HWLWtidalindicators(data_pd_HWLW_all):
    """
    computes several tidal extreme indicators from tidal extreme dataset

    Parameters
    ----------
    data_pd_HWLW_all : TYPE
        DESCRIPTION.

    Returns
    -------
    dict_tidalindicators : TYPE
        DESCRIPTION.

    """
    if hasattr(data_pd_HWLW_all.index[0],'tz'): #timezone present in index
        data_pd_HWLW_all.index = data_pd_HWLW_all.index.tz_localize(None)
    if len(data_pd_HWLW_all['HWLWcode'].unique()) > 2: #aggers are present
        data_pd_HWLW_12 = calc_HWLW12345to21(data_pd_HWLW_all) #convert 12345 to 12 by taking minimum of 345 as 2 (laagste laagwater) #TODO: this drops first/last value if it is a LW, should be fixed
    else:
        data_pd_HWLW_12 = data_pd_HWLW_all.copy()
    
    #split to HW and LW separately, also groupby year
    data_pd_HW = data_pd_HWLW_12.loc[data_pd_HWLW_12['HWLWcode']==1]
    data_pd_LW = data_pd_HWLW_12.loc[data_pd_HWLW_12['HWLWcode']==2]
    
    #yearmean HWLW from HWLW values
    HW_mean_peryear = data_pd_HW.groupby(pd.PeriodIndex(data_pd_HW.index, freq="y"))[['values']].mean()
    LW_mean_peryear = data_pd_LW.groupby(pd.PeriodIndex(data_pd_LW.index, freq="y"))[['values']].mean()
    
    #derive GHHW/GHWS (gemiddeld hoogwater springtij)
    HW_max_permonth = data_pd_HW.groupby(pd.PeriodIndex(data_pd_HW.index, freq="m"))[['values']].max() #proxy for HW at spring tide
    HW_monthmax_peryear = HW_max_permonth.groupby(pd.PeriodIndex(HW_max_permonth.index, freq="y"))[['values']].mean()
    LW_min_permonth = data_pd_LW.groupby(pd.PeriodIndex(data_pd_LW.index, freq="m"))[['values']].min() #proxy for LW at spring tide
    LW_monthmin_peryear = LW_min_permonth.groupby(pd.PeriodIndex(LW_min_permonth.index, freq="y"))[['values']].mean()
   
    dict_HWLWtidalindicators = {'HW_mean':data_pd_HW['values'].mean(), #GHW 
                                'LW_mean':data_pd_LW['values'].mean(), #GLW 
                                'HW_mean_peryear':HW_mean_peryear, #GHW peryear
                                'LW_mean_peryear':LW_mean_peryear, #GLW peryear
                                'HW_monthmax_mean':HW_max_permonth['values'].mean(), #GHHW/GHWS
                                'LW_monthmin_mean':LW_min_permonth['values'].mean(), #GLLW/GLWS
                                'HW_monthmax_mean_peryear':HW_monthmax_peryear, #GHHW/GHWS peryear
                                'LW_monthmin_mean_peryear':LW_monthmin_peryear, #GLLW/GLWS peryear
                                }
    for key in dict_HWLWtidalindicators.keys():
        if not hasattr(dict_HWLWtidalindicators[key],'index'):
            continue
        dict_HWLWtidalindicators[key].index = dict_HWLWtidalindicators[key].index.to_timestamp()
        
    return dict_HWLWtidalindicators

