# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 17:12:42 2022

@author: veenstra
"""


import numpy as np
import statsmodels.api as sm
import pandas as pd
import datetime as dt
from hatyan.timeseries import calc_HWLW12345to12
from hatyan.analysis_prediction import HatyanSettings, prediction #PydanticConfig
#from pydantic import validate_arguments #TODO: enable validator (first add pydantic as dependency, plus how to validate comp df (columns A/phi, then maybe classed should be used instead)


def calc_HWLWtidalindicators(data_pd_HWLW_all, tresh_yearlyHWLWcount=None):
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
        data_pd_HWLW_12 = calc_HWLW12345to12(data_pd_HWLW_all) #convert 12345 to 12 by taking minimum of 345 as 2 (laagste laagwater) #TODO: this drops first/last value if it is a LW, should be fixed
    else:
        data_pd_HWLW_12 = data_pd_HWLW_all.copy()
    
    #split to HW and LW separately, also groupby year
    data_pd_HW = data_pd_HWLW_12.loc[data_pd_HWLW_12['HWLWcode']==1]
    data_pd_LW = data_pd_HWLW_12.loc[data_pd_HWLW_12['HWLWcode']==2]
    
    #count HWLW values per year/month
    HWLW_count_peryear = data_pd_HWLW_12.groupby(pd.PeriodIndex(data_pd_HWLW_12.index, freq="y"))['values'].count()
    HWLW_count_permonth = data_pd_HWLW_12.groupby(pd.PeriodIndex(data_pd_HWLW_12.index, freq="m"))['values'].count()
    
    #yearmean HWLW from HWLW values #maybe also add *_mean_permonth
    HW_mean_peryear = data_pd_HW.groupby(pd.PeriodIndex(data_pd_HW.index, freq="y"))[['values']].mean()
    LW_mean_peryear = data_pd_LW.groupby(pd.PeriodIndex(data_pd_LW.index, freq="y"))[['values']].mean()
    
    #derive GHHW/GHWS (gemiddeld hoogwater springtij) per month
    HW_monthmax_permonth = data_pd_HW.groupby(pd.PeriodIndex(data_pd_HW.index, freq="m"))[['values']].max() #proxy for HW at spring tide
    LW_monthmin_permonth = data_pd_LW.groupby(pd.PeriodIndex(data_pd_LW.index, freq="m"))[['values']].min() #proxy for LW at spring tide
    
    #replace invalids with nan (in case of too less values per month or year)
    if tresh_yearlyHWLWcount is not None:
        tresh_monthlyHWLWcount = tresh_yearlyHWLWcount/13 #not 13 but 12, to also make the threshold valid in short months
        HW_mean_peryear.loc[HWLW_count_peryear<tresh_yearlyHWLWcount] = np.nan
        LW_mean_peryear.loc[HWLW_count_peryear<tresh_yearlyHWLWcount] = np.nan
        HW_monthmax_permonth.loc[HWLW_count_permonth<tresh_monthlyHWLWcount] = np.nan
        LW_monthmin_permonth.loc[HWLW_count_permonth<tresh_monthlyHWLWcount] = np.nan
    
    #derive GHHW/GHWS (gemiddeld hoogwater springtij)
    HW_monthmax_peryear = HW_monthmax_permonth.groupby(pd.PeriodIndex(HW_monthmax_permonth.index, freq="y"))[['values']].mean()
    LW_monthmin_peryear = LW_monthmin_permonth.groupby(pd.PeriodIndex(LW_monthmin_permonth.index, freq="y"))[['values']].mean()
    
    dict_HWLWtidalindicators = {'HW_mean':data_pd_HW['values'].mean(), #GHW
                                'LW_mean':data_pd_LW['values'].mean(), #GLW
                                'HW_mean_peryear':HW_mean_peryear['values'], #GHW peryear
                                'LW_mean_peryear':LW_mean_peryear['values'], #GLW peryear
                                'HW_monthmax_permonth':HW_monthmax_permonth['values'], #GHHW/GHWS permonth
                                'LW_monthmin_permonth':LW_monthmin_permonth['values'], #GLLW/GLWS permonth
                                'HW_monthmax_mean_peryear':HW_monthmax_peryear['values'], #GHHW/GHWS peryear
                                'LW_monthmin_mean_peryear':LW_monthmin_peryear['values'], #GLLW/GLWS peryear
                                }
    
    for key in dict_HWLWtidalindicators.keys():
        if not hasattr(dict_HWLWtidalindicators[key],'index'):
            continue
        dict_HWLWtidalindicators[key].index = dict_HWLWtidalindicators[key].index.to_timestamp()
        
    return dict_HWLWtidalindicators


def calc_wltidalindicators(data_wl_pd, tresh_yearlywlcount=None):
    """
    computes monthly and yearly means from waterlevel timeseries

    Parameters
    ----------
    data_wl_pd : TYPE
        DESCRIPTION.

    Returns
    -------
    dict_wltidalindicators : TYPE
        DESCRIPTION.

    """
    if hasattr(data_wl_pd.index[0],'tz'): #timezone present in index
        data_wl_pd.index = data_wl_pd.index.tz_localize(None)
    
    #count wl values per year/month
    wl_count_peryear = data_wl_pd.groupby(pd.PeriodIndex(data_wl_pd.index, freq="y"))['values'].count()
    wl_count_permonth = data_wl_pd.groupby(pd.PeriodIndex(data_wl_pd.index, freq="m"))['values'].count()
    
    #yearmean wl from wl values
    wl_mean_peryear = data_wl_pd.groupby(pd.PeriodIndex(data_wl_pd.index, freq="y"))[['values']].mean()
    wl_mean_permonth = data_wl_pd.groupby(pd.PeriodIndex(data_wl_pd.index, freq="m"))[['values']].mean()
    
    #replace invalids with nan (in case of too less values per month or year)
    if tresh_yearlywlcount is not None:
        tresh_monthlyywlcount = tresh_yearlywlcount/12
        wl_mean_peryear.loc[wl_count_peryear<tresh_yearlywlcount] = np.nan
        wl_mean_permonth.loc[wl_count_permonth<tresh_monthlyywlcount] = np.nan
        
    dict_wltidalindicators = {'wl_mean_peryear':wl_mean_peryear['values'], #yearly mean wl
                              'wl_mean_permonth':wl_mean_permonth['values'], #monthly mean wl
                              }
    
    for key in dict_wltidalindicators.keys():
        if not hasattr(dict_wltidalindicators[key],'index'):
            continue
        dict_wltidalindicators[key].index = dict_wltidalindicators[key].index.to_timestamp()
        
    return dict_wltidalindicators


#@validate_arguments(config=PydanticConfig)
def calc_HAT_LAT_fromcomponents(comp: pd.DataFrame, hatyan_settings: HatyanSettings = None) -> tuple:
    """
    Derive highest and lowest astronomical tide (HAT/LAT) from a component set.
    The component set is used to make a tidal prediction for an arbitrary period of 19 years with a 1 minute interval. The max/min values of the predictions of all years are the HAT/LAT values.
    The HAT/LAT is very dependent on the A0 of the component set. Therefore, the HAT/LAT values are relevant for the same year as the slotgemiddelde that is used to replace A0 in the component set. For instance, if the slotgemiddelde is valid for 2021.0, HAT and LAT are also relevant for that year.
    The HAT/LAT values are also very dependent on the hatyan_settings used, in general it is important to use the same settings as used to derive the tidal components.
    
    Parameters
    ----------
    comp : pd.DataFrame
        DESCRIPTION.
    hatyan_settings : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    tuple
        DESCRIPTION.

    """
    
    min_vallist_allyears = pd.Series(dtype=float)
    max_vallist_allyears = pd.Series(dtype=float)
    for year in range(2020,2039): # 19 arbitrary consequtive years to capture entire nodal cycle
        times_pred_all = pd.date_range(start=dt.datetime(year,1,1), end=dt.datetime(year+1,1,1), freq='1min')
        ts_prediction = prediction(comp=comp, hatyan_settings=hatyan_settings, times_pred_all=times_pred_all)
        
        min_vallist_allyears.loc[year] = ts_prediction['values'].min()
        max_vallist_allyears.loc[year] = ts_prediction['values'].max()
    #vallist_allyears.plot()
    #print(vallist_allyears)
    #vallist_allyears.to_csv('LAT_HAT_indication_19Y_%s.csv'%(current_station))
    HAT = max_vallist_allyears.max()
    LAT = min_vallist_allyears.min()
    return HAT, LAT


def fit_models(mean_array_todate: pd.Series) -> pd.DataFrame:
    """
    Fit linear model over yearly means in mean_array_todate, including five years in the future.

    Parameters
    ----------
    mean_array_todate : pd.Series
        DESCRIPTION.

    Returns
    -------
    pred_pd : TYPE
        DESCRIPTION.

    """
    
    
    # We'll just use the years. This assumes that annual waterlevels are used that are stored left-padded, the mean waterlevel for 2020 is stored as 2020-1-1. This is not logical, but common practice.
    allyears_DTI = pd.date_range(mean_array_todate.index.min(),mean_array_todate.index.max()+dt.timedelta(days=5*360),freq='AS')
    mean_array_allyears = pd.Series(mean_array_todate,index=allyears_DTI)
    
    df = pd.DataFrame({'year':mean_array_allyears.index.year, 'height':mean_array_allyears.values}) #TODO: make functions accept mean_array instead of df as argument?
    
    # below methods are copied from https://github.com/openearth/sealevel/blob/master/slr/slr/models.py #TODO: install slr package as dependency or keep separate?
    fit, names, X = linear_model(df, with_wind=False, with_nodal=False)
    pred_linear_nonodal = fit.predict(X)
    fit, names, X = linear_model(df, with_wind=False)
    pred_linear_winodal = fit.predict(X)
    
    pred_pd = pd.DataFrame({'pred_linear_nonodal':pred_linear_nonodal,
                            'pred_linear_winodal':pred_linear_winodal},
                            index=allyears_DTI)
    return pred_pd


# copied from https://github.com/openearth/sealevel/blob/master/slr/slr/models.py
def broken_linear_model(df, with_wind=True, quantity='height', start_acceleration=1993):
    """This model fits the sea-level rise has started to rise faster in 1993."""
    y = df[quantity]
    X = np.c_[
        df['year']-1970,
        (df['year'] > start_acceleration),# * (df['year'] - start_acceleration),
        np.cos(2*np.pi*(df['year']-1970)/18.613),
        np.sin(2*np.pi*(df['year']-1970)/18.613)
    ]
    names = ['Constant', 'Trend', f'+trend ({start_acceleration})', 'Nodal U', 'Nodal V']
    if with_wind:
        X = np.c_[
            X,
            df['u2'],
            df['v2']
        ]
        names.extend(['Wind $u^2$', 'Wind $v^2$'])
    X = sm.add_constant(X)
    model_broken_linear = sm.GLSAR(y, X, rho=1, missing='drop')
    fit = model_broken_linear.iterative_fit(cov_type='HC0', missing='drop')
    return fit, names, X


# copied from https://github.com/openearth/sealevel/blob/master/slr/slr/models.py
def linear_model(df, with_wind=True, with_ar=True, with_nodal=True, quantity='height'):
    """Define the linear model with optional wind and autoregression.
    See the latest report for a detailed description.
    """

    y = df[quantity]
    X = np.c_[df['year']-1970,
              ]
    #month = np.mod(df['year'], 1) * 12.0
    names = ['Constant', 'Trend']
    if with_nodal:
        X = np.c_[X,
                  np.cos(2*np.pi*(df['year']-1970)/18.613),
                  np.sin(2*np.pi*(df['year']-1970)/18.613)
                  ]
        names.extend(['Nodal U', 'Nodal V'])
    if with_wind:
        X = np.c_[
            X,
            df['u2'],
            df['v2']
        ]
        names.extend(['Wind $u^2$', 'Wind $v^2$'])
    X = sm.add_constant(X)
    if with_ar:
        model = sm.GLSAR(y, X, missing='drop', rho=1)
    else:
        model = sm.OLS(y, X, missing='drop')
    fit = model.fit(cov_type='HC0')
    return fit, names, X


