# -*- coding: utf-8 -*-
"""
analysis_prediction.py contains hatyan definitions related to tidal analysis and prediction. 

hatyan is a Python program for tidal analysis and prediction, based on the FORTRAN version. 
Copyright (C) 2019-2021 Rijkswaterstaat.  Maintained by Deltares, contact: Jelmer Veenstra (jelmer.veenstra@deltares.nl). 
Source code available at: https://github.com/Deltares/hatyan

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""

import numpy as np
import pandas as pd
import datetime as dt
import warnings
warnings.filterwarnings(action='always', category=DeprecationWarning)

from hatyan.hatyan_core import get_const_list_hatyan, sort_const_list, robust_timedelta_sec, robust_daterange_fromtimesextfreq
from hatyan.hatyan_core import get_freqv0_generic, get_uf_generic
from hatyan.timeseries import check_ts, nyquist_folding, check_rayleigh


class PydanticConfig:
    #https://docs.pydantic.dev/usage/model_config/
    arbitrary_types_allowed = True #necessary to allow for pd.DataFrame etc


class MatrixConditionTooHigh(Exception):
    pass


class HatyanSettings:
    """
    Settings class containing default hatyan settings, to be overwritten by input, initiate with:
    hatyan_settings = hatyan.HatyanSettings(nodalfactors=False)

    source : TYPE, optional
        DESCRIPTION. The default is 'schureman'.
    nodalfactors : bool, optional
        Whether or not to apply nodal factors. The default is True.
    fu_alltimes : bool, optional
        Whether to calculate nodal factors in middle of the analysis/prediction period (default) or on every timestep. The default is True.
    xfac : bool, optional
        Whether or not to apply x-factors. The default is False.
    
    #following are only for analysis
    CS_comps : pandas.DataFrame, optional
        contains the from/derive component lists for components splitting, as well as the amplitude factor and the increase in degrees. The default is None.
    
    #following are only for get_components_from_ts
    analysis_perperiod : False or Y/Q/W, optional
        caution, it tries to analyse each year/quarter/month, but skips if it fails. The default is False.
    return_allperiods : bool, optional
        DESCRIPTION. The default is False.
    return_prediction : bool, optional
        Whether to generate a prediction for the ts time array. The default is False.
    
    """
    #TODO: analysis_peryear,analysis_permonth,return_allyears only for get_components_from_ts, return_prediction only for analysis. Merge analysis and get_components_from_ts? Remove some from HatyanSettings class or maybe split? Add const_list to HatyanSettings?
    
    def __init__(self, source='schureman', nodalfactors=True, fu_alltimes=True, xfac=False, #prediction/analysis 
                 CS_comps=None, analysis_perperiod=False, return_allperiods=False, 
                 analysis_peryear=None, analysis_permonth=None, return_allyears=None,  #TODO: should be phased out, replaced by analysis_perperiod and return_allperiods
                 return_prediction=False,
                 xTxmat_condition_max=12): #analysis only
        if not isinstance(source,str):
            raise Exception('invalid source type, should be str')
        source = source.lower()
        if source not in ['schureman','foreman']:
            raise Exception('invalid source {source}, should be schureman or foreman)')
        
        #for arguments that will be phased out, are now written to newer arguments and overwritten with None #TODO: phase out arguments
        if analysis_peryear is not None:
            if analysis_peryear:
                warnings.warn("WARNING: argument analysis_peryear will be phased out, use analysis_perperiod='Y' instead", category=DeprecationWarning)
                analysis_perperiod = 'Y'
            else:
                warnings.warn("WARNING: argument analysis_peryear will be phased out, use analysis_perperiod=False instead", category=DeprecationWarning)
            analysis_peryear = None
        if analysis_permonth is not None:
            if analysis_permonth:
                warnings.warn("WARNING: argument analysis_permonth will be phased out, use analysis_perperiod='M' instead", category=DeprecationWarning)
                analysis_perperiod = 'M'
            else:
                warnings.warn("WARNING: argument analysis_permonth will be phased out, use analysis_perperiod=False instead", category=DeprecationWarning)
            analysis_permonth = None
        if return_allyears is not None:
            warnings.warn("WARNING: argument return_allyears will be phased out, use return_allperiods instead", category=DeprecationWarning)
            return_allperiods = return_allyears
            return_allyears = None
        
        for var_in in [nodalfactors,fu_alltimes,return_allperiods,return_prediction]:
            if not isinstance(var_in,bool):
                raise Exception(f'invalid {var_in} type, should be bool')
        
        if not (isinstance(xfac,bool) or isinstance(xfac,dict)):
            raise Exception(f'invalid xfac={xfac} type, should be bool or dict')
        
        if not ((analysis_perperiod is False) or (analysis_perperiod in ['Y','Q','M'])):
            raise Exception(f'invalid analysis_perperiod={analysis_perperiod} type, should be False or Y/Q/M')
        
        if CS_comps is not None:
            if not isinstance(CS_comps,(dict,pd.DataFrame)):
                raise Exception('invalid CS_comps type, should be dict')
            CS_comps = pd.DataFrame(CS_comps) #TODO: convert all to dict or pd.DataFrame
            CS_comps_expectedkeys = ['CS_comps_derive', 'CS_comps_from', 'CS_ampfacs', 'CS_degincrs']
            for CS_comps_key in CS_comps_expectedkeys:
                if CS_comps_key not in CS_comps.keys():
                    raise Exception(f'CS_comps does not contain {CS_comps_key}')
            CS_comps_lenvals = [len(CS_comps[key]) for key in CS_comps]
            if len(np.unique(CS_comps_lenvals)) != 1:
                raise Exception(f'CS_comps keys do not have equal lengths:\n{CS_comps}')
        
        self.source = source
        self.nodalfactors = nodalfactors
        self.fu_alltimes = fu_alltimes
        self.xfac = xfac
        self.CS_comps = CS_comps
        self.analysis_perperiod = analysis_perperiod
        self.return_allperiods = return_allperiods
        self.return_prediction = return_prediction
        self.xTxmat_condition_max = xTxmat_condition_max
        
    def __str__(self):
        self_dict = vars(self)
        str_append = ''
        for key,val in self_dict.items():
            if key=='CS_comps' and self.CS_comps is not None:
                str_append += f'{key:20s} = \n{val}\n'
            else:
                str_append += f'{key:20s} = {val}\n'
        return str_append


def vectoravg(A_all, phi_deg_all):
    """
    calculates the vector average of A and phi per constituent, 
    it vector averages over values resulting from multiple periods.
    A regular average is calculated for the amplitude of A0 (middenstand)
    
    Parameters
    ----------
    A_i_all : TYPE
        DESCRIPTION.
    phi_i_deg_all : TYPE
        DESCRIPTION.

    Returns
    -------
    A_i_mean : TYPE
        DESCRIPTION.
    phi_i_deg_mean : TYPE
        DESCRIPTION.

    """
        
    phi_rad_all = np.deg2rad(phi_deg_all)
    v_cos = A_all*np.cos(phi_rad_all)
    v_sin = A_all*np.sin(phi_rad_all)
    mean_v_cos = np.mean(v_cos,axis=1)
    mean_v_sin = np.mean(v_sin,axis=1)
    A_mean = np.sqrt(mean_v_cos**2 + mean_v_sin**2)
    phi_rad_mean = np.arctan2(mean_v_sin,mean_v_cos)
    phi_rad_mean[phi_rad_mean<0] = phi_rad_mean[phi_rad_mean<0]+(2*np.pi)
    
    #if phases of all years are exactly 0, it is the A0 component. Overwrite this A0 with mean amplitude and zero phase if present, otherwise negative values will become positive with 180 phase
    idx_A0 = np.where((phi_deg_all==0).any(axis=1))[0]
    A_mean[idx_A0] = np.mean(A_all[idx_A0,:])
    phi_rad_mean[idx_A0] = 0
    
    phi_deg_mean = np.rad2deg(phi_rad_mean)
    
    return A_mean, phi_deg_mean


def get_components_from_ts(ts, const_list, hatyan_settings=None, **kwargs):#nodalfactors=True, xfac=False, fu_alltimes=True, CS_comps=None, analysis_peryear=False, analysis_permonth=False, source='schureman'):
    """
    Wrapper around the analysis() function, 
    it optionally processes a timeseries per year and vector averages the results afterwards, 
    passes the rest of the arguments on to analysis function
    The timezone of the timeseries, will also be reflected in the phases of the resulting component set, so the resulting component set can be used to make a prediction in the original timezone.
    
    Parameters
    ----------
    ts : pandas.DataFrame
        The DataFrame should contain a 'values' column and a pd.DatetimeIndex as index, it contains the timeseries to be analysed, as obtained from e.g. readts_*.
    const_list : list, pandas.Series or str
        list or pandas.Series: contains the tidal constituent names for which to analyse the provided timeseries ts. 
        str: a predefined name of a component set for hatyan_core.get_const_list_hatyan()
    hatyan_settings : hatyan.HatyanSettings()
        Contains the used settings

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    COMP_mean_pd : pandas.DataFrame
        The DataFrame contains the component data with component names as index, and colums 'A' and 'phi_deg'.
    COMP_all_pd : pandas.DataFrame, optional
        The same as COMP_mean_pd, but with all years added with MultiIndex
    """
    ts_pd = ts #TODO: this is not necessary
    
    if hatyan_settings is None:
        hatyan_settings = HatyanSettings(**kwargs)
    elif len(kwargs)>0:
        raise Exception('both arguments hatyan_settings and other settings (e.g. nodalfactors) are provided, this is not valid')
    
    print('running: get_components_from_ts')
    
    if type(const_list) is str:
        const_list = get_const_list_hatyan(const_list)
    elif type(const_list) is not list:
        const_list = const_list.tolist()
    if hatyan_settings.CS_comps is None:
        n_const = len(const_list)
    else:
        n_const = len(const_list) + len(hatyan_settings.CS_comps)

    if hatyan_settings.analysis_perperiod:
        period = hatyan_settings.analysis_perperiod
        print(f'analysis_perperiod={hatyan_settings.analysis_perperiod}, separate periods are automatically determined from timeseries')
        ts_periods_dt = ts_pd.index.to_period(period).unique() # TODO: to_period is not limited to Y/Q/M, there are more options that are now blocked: https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#offset-aliases
        ts_periods_strlist = [str(x) for x in ts_periods_dt]
        
        n_periods = len(ts_periods_dt)
        A_i_all = np.zeros((n_const,n_periods))*np.nan
        phi_i_deg_all = np.zeros((n_const,n_periods))*np.nan
        for iP, period_dt in enumerate(ts_periods_dt):
            print('analyzing %s of sequence %s'%(period_dt,ts_periods_strlist))
            ts_oneperiod_pd = ts_pd[ts_pd.index.to_period(period)==period_dt]
            try:
                COMP_one = analysis(ts_oneperiod_pd, const_list=const_list, hatyan_settings=hatyan_settings)
                A_i_all[:,iP] = COMP_one.loc[:,'A']
                phi_i_deg_all[:,iP] = COMP_one.loc[:,'phi_deg']
            except MatrixConditionTooHigh: # accept exception if matrix condition is too high, since some years can then be skipped
                print(f'WARNING: analysis of {period_dt} failed because MatrixConditionTooHigh, check if const_list is appropriate for timeseries lenght.')
        if np.isnan(A_i_all).all():
            raise Exception('analysis peryear or permonth failed for all years/months, check warnings above')
        
        COMP_all_pd = pd.DataFrame(data=np.hstack([A_i_all,phi_i_deg_all]), columns=pd.MultiIndex.from_product([['A','phi_deg'],ts_periods_dt]), index=COMP_one.index)
        print('vector averaging analysis results')
        A_i_mean, phi_i_deg_mean = vectoravg(A_all=A_i_all, phi_deg_all=phi_i_deg_all)
        COMP_mean_pd = pd.DataFrame({ 'A': A_i_mean, 'phi_deg': phi_i_deg_mean},index=COMP_one.index)

    else: #dummy values, COMP_years should be equal to COMP_mean
        COMP_mean_pd = analysis(ts_pd, const_list=const_list, hatyan_settings=hatyan_settings)
        COMP_all_pd = None
    
    if hatyan_settings.return_allperiods:
        return COMP_mean_pd, COMP_all_pd
    else:
        return COMP_mean_pd


def analysis(ts, const_list, hatyan_settings=None, **kwargs):#nodalfactors=True, xfac=False, fu_alltimes=True, CS_comps=None, return_prediction=False, source='schureman'):
    """
    harmonic analysis with matrix transformations (least squares fit), optionally with component splitting
    for details about arguments and return variables, see get_components_from_ts() definition
    
    """
    
    if hatyan_settings is None:
        hatyan_settings = HatyanSettings(**kwargs)
    elif len(kwargs)>0:
        raise Exception('both arguments hatyan_settings and other settings (e.g. nodalfactors) are provided, this is not valid')

    print('ANALYSIS initializing')
    print(hatyan_settings)
        
    #drop duplicate times
    bool_ts_duplicated = ts.index.duplicated(keep='first')
    ts_pd = ts.copy() #TODO: this is not necessary
    if not (isinstance(ts.index[0],pd.Timestamp) or isinstance(ts.index[0],dt.datetime)): #works better than isinstance(ts.index,pd.DatetimeIndex) since 1018 time indexes are Index instead of DatetimeIndex
        raise TypeError(f'ts.index is not of expected type ({type(ts.index[0])} instead of pd.Timestamp or dt.datetime)')
    if bool_ts_duplicated.any():
        raise Exception(f'ERROR: {bool_ts_duplicated.sum()} duplicate timesteps in provided timeseries, remove them e.g. with: ts = ts[~ts.index.duplicated(keep="first")]')
    print(f'#timesteps           = {len(ts)}')
    print(f'tstart               = {ts.index[0].strftime("%Y-%m-%d %H:%M:%S")}')
    print(f'tstop                = {ts.index[-1].strftime("%Y-%m-%d %H:%M:%S")}')
    if hasattr(ts.index,'freq'):
        print(f'timestep             = {ts.index.freq}')
    
    #retrieving and sorting const_list
    if type(const_list) is str:
        const_list = get_const_list_hatyan(const_list)
    elif type(const_list) is not list:
        const_list = const_list.tolist()
    const_list = sort_const_list(const_list)
    print(f'components analyzed  = {len(const_list)}')
    
    #check for duplicate components (results in singular matrix)
    if len(const_list) != len(np.unique(const_list)):
        const_list_uniq, const_list_uniq_counts = np.unique(const_list,return_counts=True)
        const_list_counts = pd.DataFrame({'constituent':const_list_uniq,'occurences':const_list_uniq_counts})
        raise Exception('remove duplicate constituents from const_list:\n%s'%(const_list_counts.loc[const_list_counts['occurences']>1]))
    
    #check for length
    if len(ts_pd)<2:
        raise Exception('provided timeseries is less than 2 timesteps long, analysis not possible')

    #remove nans
    ts_pd_nonan = ts_pd[~ts_pd['values'].isna()]
    if len(ts_pd_nonan)==0:
        raise Exception('provided timeseries only contains nan values, analysis not possible')
    times_pred_all_pdDTI = ts_pd_nonan.index.copy() #pd.DatetimeIndex(ts_pd_nonan.index) #TODO: this will not work for OutOfBoundsDatetime
    percentage_nan = 100-len(ts_pd_nonan['values'])/len(ts_pd['values'])*100
    print(f'percentage_nan in values_meas_sel: {percentage_nan:.2f}%')
    
    #get times and time array
    dood_date_mid = pd.Index([ts_pd.index[len(ts_pd.index)//2]]) #middle of analysis period (2july in case of 1jan-1jan), zoals bij hatyan #TODO: this is incorrect in case of e.g. more missings in first half of year than second half
    dood_date_start = ts_pd.index[[0]] #first date (for v0, also freq?)
    if hatyan_settings.fu_alltimes:
        dood_date_fu = times_pred_all_pdDTI
    else:
        dood_date_fu = dood_date_mid
    times_from0_s = robust_timedelta_sec(ts_pd_nonan.index,refdate_dt=dood_date_start[0])
    times_from0_s = times_from0_s[:,np.newaxis]
    
    #get frequency and v0
    t_const_freq_pd, v_0i_rad = get_freqv0_generic(hatyan_settings, const_list, dood_date_mid, dood_date_start)
    omega_i_rads = t_const_freq_pd[['freq']].values.T*(2*np.pi)/3600 #angular frequency, 2pi/T, in rad/s, https://en.wikipedia.org/wiki/Angular_frequency (2*np.pi)/(1/x*3600) = 2*np.pi*x/3600
    u_i_rad, f_i = get_uf_generic(hatyan_settings, const_list, dood_date_fu)
    v_u = v_0i_rad.values + u_i_rad.values
    
    #check rayleigh frequency after nyquist frequency folding process.
    freq_rem = nyquist_folding(ts_pd,t_const_freq_pd)
    check_rayleigh(ts_pd,freq_rem) #TODO: maybe sometimes valuable to not fold with nyquist (eg with strongly varying time interval), in that case: check_rayleigh(ts_pd,t_const_freq_pd)

    #### TIMESERIES ANALYSIS
    N = len(const_list)
    m = len(ts_pd_nonan['values'])
    
    # get xmat and make dot product
    xmat = np.zeros((m,2*N))
    xmat[:,:N] = f_i.values * np.cos(omega_i_rads*times_from0_s+v_u)
    xmat[:,N:] = f_i.values * np.sin(omega_i_rads*times_from0_s+v_u)
    
    xTmat = xmat.T
    print('calculating xTx matrix')
    tic = dt.datetime.now()
    xTxmat = np.dot(xTmat,xmat)
    
    if 'A0' in const_list: #correct center value [N,N] for better matrix condition
        xTxmat_condition = np.linalg.cond(xTxmat)
        print('condition of xTx matrix before center adjustment for A0: %.2f'%(xTxmat_condition))
        xTxmat[N,N] = m
    xTxmat_condition = np.linalg.cond(xTxmat)
    print('condition of xTx matrix: %.2f'%(xTxmat_condition))
    if xTxmat_condition > hatyan_settings.xTxmat_condition_max:#10:#100: #random treshold
        raise MatrixConditionTooHigh(f'ERROR: condition of xTx matrix is too high ({xTxmat_condition:.2f}), check your timeseries length, try different (shorter) component set or componentsplitting.\nAnalysed {check_ts(ts_pd)}')
    xTymat = np.dot(xTmat,ts_pd_nonan['values'].values)
    
    #solve matrix to get beta_roof_mat (and thus a, b)
    beta_roof_mat = np.linalg.solve(xTxmat,xTymat)
    print('matrix system solved, elapsed time: %s'%(dt.datetime.now()-tic))

    phi_i_rad = np.arctan2(beta_roof_mat[N:],beta_roof_mat[:N]) #(a,b) arctan_ab
    A_i = np.sqrt(beta_roof_mat[N:]**2 + beta_roof_mat[:N]**2) #(a,b) sqsqrt_ab

    COMP_pd = pd.DataFrame({'A': A_i, 'phi_deg': np.rad2deg(phi_i_rad%(2*np.pi))}, index=const_list)
    if 'A0' in COMP_pd.index: #correct 180 degrees A0 phase by making amplitude value negative
        if COMP_pd.loc['A0','phi_deg']==180:
            COMP_pd.loc['A0','A'] = -COMP_pd.loc['A0','A']
            COMP_pd.loc['A0','phi_deg'] = 0
    
    if hatyan_settings.CS_comps is not None:
        COMP_pd = split_components(comp=COMP_pd, dood_date_mid=dood_date_mid, hatyan_settings=hatyan_settings)
        
    print('ANALYSIS finished')
    
    if hatyan_settings.return_prediction:
        print('immediately generating a prediction for the same time array as the input ts')
        ts_prediction = prediction(comp=COMP_pd, times_pred_all=ts_pd.index, hatyan_settings=hatyan_settings)
        return COMP_pd, ts_prediction
    else:
        return COMP_pd


def split_components(comp, dood_date_mid, hatyan_settings=None, **kwargs):
    """
    component splitting function
    for details about arguments and return variables, see get_components_from_ts() definition

    """
    
    if hatyan_settings is None:
        hatyan_settings = HatyanSettings(**kwargs)
    elif len(kwargs)>0:
        raise Exception('both arguments hatyan_settings and other settings (e.g. nodalfactors) are provided, this is not valid')
        
    #create sorted and complete component list
    const_list_inclCS_raw = comp.index.tolist() + hatyan_settings.CS_comps['CS_comps_derive'].tolist()
    const_list_inclCS = sort_const_list(const_list=const_list_inclCS_raw)

    #retrieve freq and speed
    t_const_freq_pd, CS_v_0i_rad = get_freqv0_generic(hatyan_settings, const_list=const_list_inclCS, dood_date_mid=dood_date_mid, dood_date_start=dood_date_mid) # with split_components, v0 is calculated on the same timestep as u and f (middle of original series)
    CS_u_i_rad, CS_f_i = get_uf_generic(hatyan_settings, const_list=const_list_inclCS, dood_date_fu=dood_date_mid)
    
    comp_inclCS = pd.DataFrame(comp,index=const_list_inclCS,columns=comp.columns)
    #comp_inclCS_preCS = comp_inclCS.copy()
        
    for comp_main in np.unique(hatyan_settings.CS_comps['CS_comps_from']):
        bool_CS_maincomp = hatyan_settings.CS_comps['CS_comps_from'] == comp_main #boolean of which rows of CS_comps dataframe corresponds to a main constituent, also makes it possible to select two rows
        CS_comps_formain = hatyan_settings.CS_comps.loc[bool_CS_maincomp]
        comp_slave_list = CS_comps_formain['CS_comps_derive'].tolist()
        print(f'splitting component {comp_main} into {comp_slave_list}')
        
        #first update main components based on degincrs/ampfacs of all components that are to be derived
        DBETA = 0
        for iR, CS_comps_row in CS_comps_formain.iterrows():
            comp_slave = CS_comps_row['CS_comps_derive']
            #code from resuda.f, line 440 to 455
            DMU = CS_f_i.loc[0,comp_slave]/CS_f_i.loc[0,comp_main]
            DTHETA = CS_comps_row['CS_ampfacs']
            DGAMMA = np.deg2rad(CS_comps_row['CS_degincrs'])-DBETA-(CS_v_0i_rad.loc[0,comp_slave]+CS_u_i_rad.loc[0,comp_slave])+(CS_v_0i_rad.loc[0,comp_main]+CS_u_i_rad.loc[0,comp_main]) #in FORTRAN code, CS_f_i slave/main is also added, this seems wrong
            DREEEL = 1+DMU*DTHETA*np.cos(DGAMMA)
            DIMAGI = DMU*DTHETA*np.sin(DGAMMA)  
            DALPHA = np.sqrt(DREEEL*DREEEL+DIMAGI*DIMAGI)
            if DALPHA < 1e-50:
                raise Exception('ERROR: DALPHA too small, component splitting failed?')
            DBETA = np.arctan2(DIMAGI,DREEEL)
            
            comp_inclCS.loc[comp_main,'A'] = comp_inclCS.loc[comp_main,'A']/DALPHA
            comp_inclCS.loc[comp_main,'phi_deg'] = (comp_inclCS.loc[comp_main,'phi_deg']-np.rad2deg(DBETA))%360 
        
        #updating slave components after updating main components, this makes a difference when splitting a component into more than two
        for iR, CS_comps_row in CS_comps_formain.iterrows():
            comp_slave = CS_comps_row['CS_comps_derive']
            comp_inclCS.loc[comp_slave,'A'] = comp_inclCS.loc[comp_main,'A'] * CS_comps_row['CS_ampfacs']
            comp_inclCS.loc[comp_slave,'phi_deg'] = (comp_inclCS.loc[comp_main,'phi_deg'] + CS_comps_row['CS_degincrs'])%360
            
    return comp_inclCS


def prediction(comp, times_pred_all=None, times_ext=None, timestep_min=None, hatyan_settings=None, **kwargs):
    """
    generates a tidal prediction from a set of components A and phi values.
    The component set has the same timezone as the timeseries used to create it, therefore the resulting prediction will also be in that original timezone.
    
    Parameters
    ----------
    comp : pandas.DataFrame
        The DataFrame contains the component data with component names as index, and colums 'A' and 'phi_deg'.
    times_pred_all : pandas.DatetimeIndex, optional
        Prediction timeseries. The default is None.
    times_ext : list of datetime.datetime, optional
        Prediction time extents (list of start time and stop time). The default is None.
    timestep_min : int, optional
        Prediction timestep in minutes. The default is None.
    hatyan_settings : hatyan.HatyanSettings()
        Contains the used settings

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    ts_prediction_pd : pandas.DataFrame
        The DataFrame should contain a 'values' column and a pd.DatetimeIndex as index, it contains the prediction times and values.
    
    """
    
    if hatyan_settings is None:
        hatyan_settings = HatyanSettings(**kwargs)
    elif len(kwargs)>0:
        raise Exception('both arguments hatyan_settings and other settings (e.g. nodalfactors) are provided, this is not valid')
    
    print('PREDICTION initializing')
    print(hatyan_settings)
    
    if times_pred_all is None:
        if times_ext is None or timestep_min is None:
            raise Exception('if argument times_pred_all is not provided, the arguments times_ext and timestep_min are obligatory')
        else:
            times_pred_all = robust_daterange_fromtimesextfreq(times_ext,timestep_min)
    else:
        if times_ext is not None or timestep_min is not None:
            raise Exception('if argument times_pred_all is provided, the arguments times_ext and timestep_min are not allowed')
    
    if not len(times_pred_all) > 1:
        raise Exception('ERROR: requested prediction period is not more than one timestep_min')
    
    if isinstance(times_pred_all, pd.core.indexes.datetimes.DatetimeIndex) or isinstance(times_pred_all, pd.core.indexes.base.Index): #TODO: this is probably not necessary, maybe only in case of Index (year 1600 compatibility)
        times_pred_all_pdDTI = times_pred_all
    else:
        times_pred_all_pdDTI = pd.DatetimeIndex(times_pred_all)
    
    print('%-20s = %s'%('components used',len(comp)))
    print('%-20s = %s'%('tstart',times_pred_all_pdDTI[0].strftime('%Y-%m-%d %H:%M:%S')))
    print('%-20s = %s'%('tstop',times_pred_all_pdDTI[-1].strftime('%Y-%m-%d %H:%M:%S')))
    if hasattr(times_pred_all_pdDTI,'freq'):
        print('%-20s = %s'%('timestep',times_pred_all_pdDTI.freq))
    
    dood_date_mid = pd.Index([times_pred_all_pdDTI[len(times_pred_all_pdDTI)//2]]) #middle of analysis period (2july in case of 1jan-1jan), zoals bij hatyan.
    dood_date_start = times_pred_all_pdDTI[:1] #first date (for v0, also freq?)
    
    #sort component list and component dataframe
    if np.isnan(comp.values).any():
        raise Exception('provided component set contains nan values, prediction not possible')
    const_list = sort_const_list(comp.index.tolist())
    COMP = comp.loc[const_list]
    A = np.array(COMP['A'])
    phi_rad = np.array(np.deg2rad(COMP['phi_deg']))

    t_const_freq_pd, v_0i_rad = get_freqv0_generic(hatyan_settings, const_list, dood_date_mid, dood_date_start)
    t_const_speed_all = t_const_freq_pd['freq'].values[:,np.newaxis]*(2*np.pi)

    if hatyan_settings.fu_alltimes:
        dood_date_fu = times_pred_all_pdDTI
    else:
        dood_date_fu = dood_date_mid
    u_i_rad, f_i = get_uf_generic(hatyan_settings, const_list, dood_date_fu)

    print('PREDICTION started')
    omega_i_rads = t_const_speed_all.T/3600 #angular frequency, 2pi/T, in rad/s, https://en.wikipedia.org/wiki/Angular_frequency (2*np.pi)/(1/x*3600) = 2*np.pi*x/3600
    if not isinstance(times_pred_all_pdDTI,pd.DatetimeIndex): #support for years<1677, have to use Index instead of DatetimeIndex (DatetimeIndex is also Index, so isinstance(times_pred_all_pdDTI,pd.Index) does not work
        times_from0allpred_s_orig = (times_pred_all_pdDTI-dood_date_start).total_seconds().values
    else:
        times_from0allpred_s_orig = (times_pred_all_pdDTI-dood_date_start[0]).total_seconds().values
    times_from0allpred_s = np.transpose(times_from0allpred_s_orig[np.newaxis])
    
    f_A = np.multiply(f_i.values,A)
    omeg_t = np.multiply(times_from0allpred_s,omega_i_rads)#_td)
    v_u_phi = np.subtract(np.add(v_0i_rad.values,u_i_rad.values),phi_rad)
    omeg_t_v_u_phi = np.add(omeg_t,v_u_phi)
    ht_res = np.sum(np.multiply(f_A,np.cos(omeg_t_v_u_phi)),axis=1) #not necessary to add A0, since it is already part of the component list
    
    ts_prediction_pd = pd.DataFrame({'values': ht_res},index=times_pred_all_pdDTI)
    print('PREDICTION finished')
    
    return ts_prediction_pd


def prediction_peryear(comp_allyears, timestep_min, hatyan_settings=None, **kwargs):
    raise Exception('ERROR: prediction_peryear() is deprecated, use prediction_perperiod() instead')


def prediction_perperiod(comp_allperiods, timestep_min, hatyan_settings=None, **kwargs):
    """
    Wrapper around prediction(), to use component set of multiple years/months to generate multi-year/month timeseries.

    Parameters
    ----------
    comp_allperiods : TYPE
        DESCRIPTION.
    timestep_min : TYPE
        DESCRIPTION.
    hatyan_settings : hatyan.HatyanSettings()
        Contains the used settings
        
    Returns
    -------
    ts_prediction_peryear : TYPE
        DESCRIPTION.

    """
    
    if hatyan_settings is None:
        hatyan_settings = HatyanSettings(**kwargs)
    elif len(kwargs)>0:
        raise Exception('both arguments hatyan_settings and other settings (e.g. nodalfactors) are provided, this is not valid')

    ts_periods_dt = comp_allperiods.columns.levels[1]
    ts_periods_strlist = [str(x) for x in ts_periods_dt]
    
    ts_prediction_perperiod = pd.DataFrame()
    for period_dt in ts_periods_dt:
        print(f'generating prediction {period_dt} of sequence {ts_periods_strlist}')
        comp_oneyear = comp_allperiods.loc[:,(slice(None),period_dt)]
        comp_oneyear.columns = comp_oneyear.columns.droplevel(1)
        if period_dt.freqstr in ['A-DEC']: #year frequency
            times_ext = [dt.datetime(period_dt.year,1,1),dt.datetime(period_dt.year+1,1,1)-dt.timedelta(minutes=timestep_min)]
        elif period_dt.freqstr in ['M']: #month frequency
            times_ext = [period_dt.to_timestamp().to_pydatetime(),period_dt.to_timestamp().to_pydatetime()+dt.timedelta(days=period_dt.days_in_month)-dt.timedelta(minutes=timestep_min)]
        else:
            raise Exception(f'unknown freqstr: {period_dt.freqstr}')
        ts_prediction_oneperiod = prediction(comp=comp_oneyear,times_ext=times_ext, timestep_min=timestep_min, hatyan_settings=hatyan_settings)
        ts_prediction_perperiod = pd.concat([ts_prediction_perperiod,ts_prediction_oneperiod])
    return ts_prediction_perperiod
