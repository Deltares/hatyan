# -*- coding: utf-8 -*-
"""
analysis_prediction.py contains hatyan definitions related to tidal analysis and prediction. 
"""

import numpy as np
import pandas as pd
import datetime as dt
import logging

from hatyan.hatyan_core import get_const_list_hatyan, sort_const_list, robust_timedelta_sec
from hatyan.hatyan_core import get_freqv0_generic, get_uf_generic
from hatyan.timeseries import Timeseries_Statistics, nyquist_folding, check_rayleigh
from hatyan.metadata import metadata_from_obj, metadata_add_to_obj
from hatyan.deprecated import (deprecated_python_option,
                               DEPRECATED_OPTIONS_PREDICTION_DICT,
                               DEPRECATED_OPTIONS_ANALYSIS_DICT)

__all__ = ["analysis",
           "prediction",
           ]

logger = logging.getLogger(__name__)


class PydanticConfig:
    #https://docs.pydantic.dev/usage/model_config/
    arbitrary_types_allowed = True #necessary to allow for pd.DataFrame etc


class MatrixConditionTooHigh(Exception):
    pass


class HatyanSettings:
    """
    Settings class containing default hatyan settings, to be overwritten by input, initiate with:
    hatyan_settings = hatyan.HatyanSettings(nodalfactors=False)

    """
    
    def __init__(self, 
                 nodalfactors, fu_alltimes, xfac, source, #prediction/analysis 
                 cs_comps=None, analysis_perperiod=None, return_allperiods=None, #analysis only
                 max_matrix_condition=None): #analysis only
        
        if not isinstance(nodalfactors,bool):
            raise TypeError(f'invalid nodalfactors={nodalfactors} type, should be bool')
        self.nodalfactors = nodalfactors
        
        if not isinstance(fu_alltimes,bool):
            raise TypeError(f'invalid fu_alltimes={fu_alltimes} type, should be bool')
        self.fu_alltimes = fu_alltimes
        
        if not (isinstance(xfac,bool) or isinstance(xfac,dict)):
            raise TypeError(f'invalid xfac={xfac} type, should be bool or dict')
        self.xfac = xfac
        
        if not isinstance(source,str):
            raise TypeError(f'invalid source={source} type, should be str')
        source = source.lower()
        if source not in ['schureman','foreman']:
            raise TypeError(f'invalid source {source}, should be "schureman" or "foreman"')
        self.source = source
                
        if return_allperiods is not None:
            if not isinstance(return_allperiods,bool):
                raise TypeError(f'invalid return_allperiods={return_allperiods} type, should be bool')
            self.return_allperiods = return_allperiods
            if return_allperiods and not analysis_perperiod:
                raise TypeError(
                    f"return_allperiods={return_allperiods}, but analysis_perperiod="
                    f"{analysis_perperiod}, this is not supported.")
        
        if analysis_perperiod is not None:
            if not ((analysis_perperiod is False) or (analysis_perperiod in ['Y','Q','M'])):
                raise TypeError(f'invalid analysis_perperiod={analysis_perperiod} type, should be False or Y/Q/M')
            self.analysis_perperiod = analysis_perperiod
        
        if cs_comps is not None:
            if not isinstance(cs_comps,(dict,pd.DataFrame)):
                raise TypeError('invalid cs_comps type, should be dict')
            cs_comps = pd.DataFrame(cs_comps) #TODO: convert all to dict or pd.DataFrame
            cs_comps_expectedkeys = ['CS_comps_derive', 'CS_comps_from', 'CS_ampfacs', 'CS_degincrs']
            for cs_comps_key in cs_comps_expectedkeys:
                if cs_comps_key not in cs_comps.keys():
                    raise KeyError(f'cs_comps does not contain {cs_comps_key}')
            cs_comps_lenvals = [len(cs_comps[key]) for key in cs_comps]
            if len(np.unique(cs_comps_lenvals)) != 1:
                raise ValueError(f'cs_comps keys do not have equal lengths:\n{cs_comps}')
            self.cs_comps = cs_comps
        
        if max_matrix_condition is not None:
            if not isinstance(max_matrix_condition,int) or isinstance(max_matrix_condition,float):
                raise TypeError(f'invalid {max_matrix_condition} type, should be int or float')
            self.max_matrix_condition = max_matrix_condition
        
    def __str__(self):
        self_dict = vars(self)
        str_append = ''
        for key,val in self_dict.items():
            if key=='cs_comps':
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
    idx_A0 = np.nonzero((phi_deg_all==0).any(axis=1))[0]
    A_mean[idx_A0] = np.mean(A_all[idx_A0,:])
    phi_rad_mean[idx_A0] = 0
    
    phi_deg_mean = np.rad2deg(phi_rad_mean)
    
    return A_mean, phi_deg_mean


@deprecated_python_option(**DEPRECATED_OPTIONS_ANALYSIS_DICT)
def analysis(ts, const_list, 
             nodalfactors=True, fu_alltimes=True, xfac=False, 
             source='schureman',
             cs_comps=None,
             analysis_perperiod=False, return_allperiods=False,
             max_matrix_condition=12):
    """
    Analysis of timeseries.
    Optionally processes a timeseries per year and vector averages the results afterwards.
    The timezone of the timeseries, will also be reflected in the phases of the resulting component set, 
    so the resulting component set can be used to make a prediction in the original timezone.
    
    Parameters
    ----------
    ts : pandas.DataFrame
        The DataFrame should contain a 'values' column and a pd.DatetimeIndex as index, it contains the timeseries to be analysed, as obtained from e.g. readts_*.
    const_list : list, pandas.Series or str
        list or pandas.Series: contains the tidal constituent names for which to analyse the provided timeseries ts. 
        str: a predefined name of a component set for hatyan_core.get_const_list_hatyan()
    nodalfactors : bool
        Whether or not to apply nodal factors. The default is True.
    fu_alltimes : bool
        Whether to calculate nodal factors in middle of the analysis/prediction period (default) or on every timestep. The default is True.
    xfac : bool
        Whether or not to apply x-factors. The default is False.
    source : TYPE
        DESCRIPTION. The default is 'schureman'.
    cs_comps : pandas.DataFrame, optional
        contains the from/derive component lists for components splitting, as well as the amplitude factor and the increase in degrees. Only relevant for analysis. The default is None.
    max_matrix_condition: float or int
        the maximum condition of the xTx matrix. The default is 12.
    analysis_perperiod : False or Y/Q/W, optional
        caution, it tries to analyse each year/quarter/month, but skips if it fails. The default is False.
    return_allperiods : bool, optional
        Only relevant if analysis_perperiod is not None. The default is False.
    

    Returns
    -------
    COMP_mean_pd : pandas.DataFrame
        The DataFrame contains the component data with component names as index, and colums 'A' and 'phi_deg'.
    COMP_all_pd : pandas.DataFrame, optional
        The same as COMP_mean_pd, but with all years added with MultiIndex
    """
    # check ts.index type, this works better than isinstance(ts.index,pd.DatetimeIndex) since 1018 time indexes are Index instead of DatetimeIndex
    if not (isinstance(ts.index[0],pd.Timestamp) or isinstance(ts.index[0],dt.datetime)):
        raise TypeError(f'ts.index is not of expected type ({type(ts.index[0])} instead of pd.Timestamp or dt.datetime)')
    
    # remove timezone from timeseries (is added to components dataframe after analysis)
    ts_pd = ts.copy()
    ts_pd.index = ts_pd.index.tz_localize(None)
    
    # validate settings
    hatyan_settings = HatyanSettings(source=source, nodalfactors=nodalfactors, fu_alltimes=fu_alltimes, xfac=xfac,
                                     analysis_perperiod=analysis_perperiod, return_allperiods=return_allperiods,
                                     cs_comps=cs_comps, max_matrix_condition=max_matrix_condition)
    
    logger.info(f'ANALYSIS initializing\n{hatyan_settings}')
        
    #retrieving and sorting const_list
    if type(const_list) is str:
        const_list = get_const_list_hatyan(const_list)
    const_list = sort_const_list(const_list)
    logger.info(f'n components analyzed  = {len(const_list)}')
    
    #check for duplicate components (results in singular matrix)
    if len(const_list) != len(np.unique(const_list)):
        const_list_uniq, const_list_uniq_counts = np.unique(const_list,return_counts=True)
        const_list_counts = pd.DataFrame({'constituent':const_list_uniq,'occurences':const_list_uniq_counts})
        raise ValueError('remove duplicate constituents from const_list:\n%s'%(const_list_counts.loc[const_list_counts['occurences']>1]))
    
    n_const = len(const_list)
    if hasattr(hatyan_settings, "cs_comps"):
        n_const = len(const_list) + len(hatyan_settings.cs_comps)

    if hatyan_settings.analysis_perperiod:
        period = hatyan_settings.analysis_perperiod
        logger.info(f'analysis_perperiod={period}, separate periods are automatically determined from timeseries')
        ts_periods_dt_all = ts_pd.index.to_period(period)
        ts_periods_dt = ts_periods_dt_all.unique() # TODO: to_period is not limited to Y/Q/M, there are more options that are now blocked: https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#offset-aliases
        ts_periods_strlist = [str(x) for x in ts_periods_dt]
        
        n_periods = len(ts_periods_dt)
        A_i_all = np.zeros((n_const,n_periods))*np.nan
        phi_i_deg_all = np.zeros((n_const,n_periods))*np.nan
        for iP, period_dt in enumerate(ts_periods_dt):
            logger.info('analyzing %s of sequence %s'%(period_dt,ts_periods_strlist))
            ts_oneperiod_pd = ts_pd[ts_periods_dt_all==period_dt]
            try:
                COMP_one = analysis_singleperiod(ts_oneperiod_pd, const_list=const_list, hatyan_settings=hatyan_settings)
                A_i_all[:,iP] = COMP_one.loc[:,'A']
                phi_i_deg_all[:,iP] = COMP_one.loc[:,'phi_deg']
            except MatrixConditionTooHigh: # accept exception if matrix condition is too high, since some years can then be skipped
                logger.warning(f'analysis of {period_dt} failed because MatrixConditionTooHigh, check if const_list is appropriate for timeseries lenght.')
        if np.isnan(A_i_all).all():
            raise ValueError('all nans: analysis perperiod failed for all periods, check warnings above')
        
        COMP_all_pd = pd.DataFrame(data=np.hstack([A_i_all,phi_i_deg_all]), columns=pd.MultiIndex.from_product([['A','phi_deg'],ts_periods_dt]), index=COMP_one.index)
        logger.info('vector averaging analysis results')
        A_i_mean, phi_i_deg_mean = vectoravg(A_all=A_i_all, phi_deg_all=phi_i_deg_all)
        COMP_mean_pd = pd.DataFrame({ 'A': A_i_mean, 'phi_deg': phi_i_deg_mean},index=COMP_one.index)
    else: #dummy values, COMP_years should be equal to COMP_mean
        COMP_mean_pd = analysis_singleperiod(ts_pd, const_list=const_list, hatyan_settings=hatyan_settings)
        COMP_all_pd = None
    
    #add metadata
    metadata = metadata_from_obj(ts_pd)
    metadata['nodalfactors'] = hatyan_settings.nodalfactors
    metadata['xfac'] = hatyan_settings.xfac
    metadata['fu_alltimes'] = hatyan_settings.fu_alltimes
    metadata['source'] = hatyan_settings.source
    metadata['tstart'] = ts.index.min().tz_localize(None)
    metadata['tstop'] = ts.index.max().tz_localize(None)
    metadata['tzone'] = ts.index.tz
    COMP_mean_pd = metadata_add_to_obj(COMP_mean_pd, metadata)
    
    if hatyan_settings.return_allperiods:
        COMP_all_pd = metadata_add_to_obj(COMP_all_pd, metadata)
        return COMP_mean_pd, COMP_all_pd
    
    logger.info('ANALYSIS finished')
    
    return COMP_mean_pd


def analysis_singleperiod(ts, const_list, hatyan_settings):
    """
    harmonic analysis with matrix transformations (least squares fit), optionally with component splitting
    for details about arguments and return variables, see analysis() definition
    
    """
    
    #drop duplicate times
    bool_ts_duplicated = ts.index.duplicated(keep='first')
    ts_pd = ts.copy() #TODO: this is not necessary
    if bool_ts_duplicated.any():
        ts_dupl = ts[ts.index.duplicated(keep=False)]
        raise ValueError(f'{bool_ts_duplicated.sum()} duplicate timesteps in provided '
                         'timeseries, can be removed with: '
                         f'ts = ts[~ts.index.duplicated(keep="first")]:\n{ts_dupl}')
    message = (f'#timesteps    = {len(ts)}\n'
               f'tstart        = {ts.index.min().strftime("%Y-%m-%d %H:%M:%S")}\n'
               f'tstop         = {ts.index.max().strftime("%Y-%m-%d %H:%M:%S")}\n'
               f'timestep      = {ts.index.freq}')
    logger.info(message)

    #remove nans
    ts_pd_nonan = ts_pd[~ts_pd['values'].isna()]
    if len(ts_pd_nonan)<2:
        raise ValueError(f'provided timeseries is less than 2 timesteps long (after dropping potential nans), analysis not possible:\n{ts_pd_nonan}')
    times_pred_all_pdDTI = ts_pd_nonan.index.copy() #pd.DatetimeIndex(ts_pd_nonan.index) #TODO: this will not work for OutOfBoundsDatetime
    percentage_nan = 100-len(ts_pd_nonan['values'])/len(ts_pd['values'])*100
    logger.info(f'percentage_nan in values_meas_sel: {percentage_nan:.2f}%')
    
    #get times and time array
    dood_date_mid = ts_pd.index[[len(ts_pd.index)//2]] #middle of analysis period (2july in case of 1jan-1jan), zoals bij hatyan #TODO: this is incorrect in case of e.g. more missings in first half of year than second half
    dood_date_start = ts_pd.index[[0]] #first date (for v0, also freq?)
    if hatyan_settings.fu_alltimes:
        dood_date_fu = times_pred_all_pdDTI
    else:
        dood_date_fu = dood_date_mid
    times_from0_s = robust_timedelta_sec(ts_pd_nonan.index,refdate_dt=dood_date_start[0])
    times_from0_s = times_from0_s[:,np.newaxis]
    
    #get frequency and v0
    t_const_freq_pd, v_0i_rad = get_freqv0_generic(const_list, dood_date_mid, dood_date_start, hatyan_settings.source)
    omega_i_rads = t_const_freq_pd[['freq']].values.T*(2*np.pi)/3600 #angular frequency, 2pi/T, in rad/s, https://en.wikipedia.org/wiki/Angular_frequency (2*np.pi)/(1/x*3600) = 2*np.pi*x/3600
    u_i_rad, f_i = get_uf_generic(const_list, dood_date_fu, hatyan_settings.nodalfactors, hatyan_settings.xfac, hatyan_settings.source)
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
    logger.info('calculating xTx matrix')
    tic = dt.datetime.now()
    xTxmat = np.dot(xTmat,xmat)
    
    if 'A0' in const_list: #correct center value [N,N] for better matrix condition
        xTxmat_condition = np.linalg.cond(xTxmat)
        logger.debug('condition of xTx matrix before center adjustment for A0: %.2f'%(xTxmat_condition))
        xTxmat[N,N] = m
    xTxmat_condition = np.linalg.cond(xTxmat)
    logger.info('condition of xTx matrix: %.2f'%(xTxmat_condition))
    if xTxmat_condition > hatyan_settings.max_matrix_condition:#10:#100: #random treshold
        raise MatrixConditionTooHigh(f'ERROR: condition of xTx matrix is too high ({xTxmat_condition:.2f}), check your timeseries length, try different (shorter) component set or componentsplitting.\nAnalysed {Timeseries_Statistics(ts_pd)}')
    xTymat = np.dot(xTmat,ts_pd_nonan['values'].values)
    
    #solve matrix to get beta_roof_mat (and thus a, b)
    beta_roof_mat = np.linalg.solve(xTxmat,xTymat)
    logger.info('matrix system solved, elapsed time: %s'%(dt.datetime.now()-tic))

    phi_i_rad = np.arctan2(beta_roof_mat[N:],beta_roof_mat[:N]) #(a,b) arctan_ab
    A_i = np.sqrt(beta_roof_mat[N:]**2 + beta_roof_mat[:N]**2) #(a,b) sqsqrt_ab

    COMP_pd = pd.DataFrame({'A': A_i, 'phi_deg': np.rad2deg(phi_i_rad%(2*np.pi))}, index=const_list)
    if 'A0' in COMP_pd.index: #correct 180 degrees A0 phase by making amplitude value negative
        if COMP_pd.loc['A0','phi_deg']==180:
            COMP_pd.loc['A0','A'] = -COMP_pd.loc['A0','A']
            COMP_pd.loc['A0','phi_deg'] = 0
    
    if hasattr(hatyan_settings, "cs_comps"):
        COMP_pd = split_components(comp=COMP_pd, dood_date_mid=dood_date_mid, hatyan_settings=hatyan_settings)
        
    return COMP_pd


def split_components(comp, dood_date_mid, hatyan_settings):
    """
    component splitting function
    for details about arguments and return variables, see analysis() definition

    """
    
    for cs_compname in hatyan_settings.cs_comps['CS_comps_derive']:
        if cs_compname in comp.index:
            raise ValueError(f"component {cs_compname} requested via component splitting, but already present in const_list")
    
    #create sorted and complete component list
    const_list_inclcs_raw = comp.index.tolist() + hatyan_settings.cs_comps['CS_comps_derive'].tolist()
    const_list_inclcs = sort_const_list(const_list=const_list_inclcs_raw)

    #retrieve freq and speed
    _, cs_v_0i_rad = get_freqv0_generic(const_list=const_list_inclcs, dood_date_mid=dood_date_mid, dood_date_start=dood_date_mid, source=hatyan_settings.source) # with split_components, v0 is calculated on the same timestep as u and f (middle of original series)
    cs_u_i_rad, cs_f_i = get_uf_generic(const_list=const_list_inclcs, dood_date_fu=dood_date_mid, nodalfactors=hatyan_settings.nodalfactors, xfac=hatyan_settings.xfac, source=hatyan_settings.source)
    
    comp_inclcs = pd.DataFrame(comp,index=const_list_inclcs,columns=comp.columns)
        
    for comp_main in np.unique(hatyan_settings.cs_comps['CS_comps_from']):
        # boolean of which rows of cs_comps dataframe corresponds to a main constituent, also makes it possible to select two rows
        bool_cs_maincomp = hatyan_settings.cs_comps['CS_comps_from'] == comp_main
        cs_comps_formain = hatyan_settings.cs_comps.loc[bool_cs_maincomp]
        comp_slave_list = cs_comps_formain['CS_comps_derive'].tolist()
        logger.info(f'splitting component {comp_main} into {comp_slave_list}')
        
        #first update main components based on degincrs/ampfacs of all components that are to be derived
        DBETA = 0
        for iR, cs_comps_row in cs_comps_formain.iterrows():
            comp_slave = cs_comps_row['CS_comps_derive']
            #code from resuda.f, line 440 to 455
            DMU = cs_f_i.loc[0,comp_slave]/cs_f_i.loc[0,comp_main]
            DTHETA = cs_comps_row['CS_ampfacs']
            DGAMMA = (np.deg2rad(cs_comps_row['CS_degincrs']) - 
                      DBETA - 
                      (cs_v_0i_rad.loc[0,comp_slave] + cs_u_i_rad.loc[0,comp_slave]) + 
                      (cs_v_0i_rad.loc[0,comp_main] + cs_u_i_rad.loc[0,comp_main])) #in FORTRAN code, CS_f_i slave/main is also added, this seems wrong
            DREEEL = 1+DMU*DTHETA*np.cos(DGAMMA)
            DIMAGI = DMU*DTHETA*np.sin(DGAMMA)  
            DALPHA = np.sqrt(DREEEL*DREEEL+DIMAGI*DIMAGI)
            if DALPHA < 1e-50:
                raise Exception('ERROR: DALPHA too small, component splitting failed?')
            DBETA = np.arctan2(DIMAGI,DREEEL)
            
            comp_inclcs.loc[comp_main,'A'] = comp_inclcs.loc[comp_main,'A']/DALPHA
            comp_inclcs.loc[comp_main,'phi_deg'] = (comp_inclcs.loc[comp_main,'phi_deg']-np.rad2deg(DBETA))%360 
        
        #updating slave components after updating main components, this makes a difference when splitting a component into more than two
        for iR, cs_comps_row in cs_comps_formain.iterrows():
            comp_slave = cs_comps_row['CS_comps_derive']
            comp_inclcs.loc[comp_slave,'A'] = comp_inclcs.loc[comp_main,'A'] * cs_comps_row['CS_ampfacs']
            comp_inclcs.loc[comp_slave,'phi_deg'] = (comp_inclcs.loc[comp_main,'phi_deg'] + cs_comps_row['CS_degincrs'])%360
            
    return comp_inclcs


def prediction_singleperiod(comp:pd.DataFrame, times:pd.DatetimeIndex, hatyan_settings) -> pd.DataFrame:
    
    metadata_comp = metadata_from_obj(comp)
    tzone_comp = metadata_comp.pop('tzone', None)
    
    tzone_pred = times.tz
    
    if tzone_pred is None and tzone_comp is None:
        tzone_convert = False
    elif tzone_pred is not None and tzone_comp is not None:
        tzone_convert = True
    else:
        raise ValueError("provided times and components should both be timezone-aware "
                         "or timezone-naive, not mixed.")
    
    # remove timezone from prediction times: first convert times to tzone of components, then make timezone naive
    if tzone_convert:
        times = times.tz_convert(tzone_comp)
        times = times.tz_localize(None)
        # TODO: workaround to maintain the frequency dropped by tz_localize
        # https://github.com/pandas-dev/pandas/issues/36575
        times.freq = times.inferred_freq
    
    logger.info(f'components used = {len(comp)}\n'
                f'tstart = {times[0].strftime("%Y-%m-%d %H:%M:%S")}\n'
                f'tstop = {times[-1].strftime("%Y-%m-%d %H:%M:%S")}\n'
                f'timestep = {times.freq}')
    
    # middle of analysis period (2july in case of 1jan-1jan), zoals bij hatyan.
    dood_date_mid = times[[len(times)//2]]
    # first date (for v0, also freq?)
    dood_date_start = times[:1]
    
    # sort component list and component dataframe
    if np.isnan(comp.values).any():
        raise Exception('provided component set contains nan values, prediction not possible')
    const_list = sort_const_list(comp.index.tolist())
    COMP = comp.loc[const_list]
    A = np.array(COMP['A'])
    phi_rad = np.array(np.deg2rad(COMP['phi_deg']))

    t_const_freq_pd, v_0i_rad = get_freqv0_generic(const_list, dood_date_mid, dood_date_start, hatyan_settings.source)
    t_const_speed_all = t_const_freq_pd['freq'].values[:,np.newaxis]*(2*np.pi)

    if hatyan_settings.fu_alltimes:
        dood_date_fu = times
    else:
        dood_date_fu = dood_date_mid
    u_i_rad, f_i = get_uf_generic(const_list, dood_date_fu, hatyan_settings.nodalfactors, hatyan_settings.xfac, hatyan_settings.source)

    logger.info('PREDICTION started')
    omega_i_rads = t_const_speed_all.T/3600 #angular frequency, 2pi/T, in rad/s, https://en.wikipedia.org/wiki/Angular_frequency (2*np.pi)/(1/x*3600) = 2*np.pi*x/3600
    if not isinstance(times,pd.DatetimeIndex): #support for years<1677, have to use Index instead of DatetimeIndex (DatetimeIndex is also Index, so isinstance(times_pred_all_pdDTI,pd.Index) does not work
        tdiff = pd.TimedeltaIndex(times-dood_date_start) #pd.TimedeltaIndex is around it to avoid it being an Index in case of outofbounds timesteps (necessary from pandas 2.0.0)
    else:
        tdiff = pd.TimedeltaIndex(times-dood_date_start[0]) #pd.TimedeltaIndex is not necessary here, but for conformity with above
    times_from0allpred_s_orig = tdiff.total_seconds().values
    times_from0allpred_s = np.transpose(times_from0allpred_s_orig[np.newaxis])
    
    f_A = np.multiply(f_i.values,A)
    omeg_t = np.multiply(times_from0allpred_s,omega_i_rads)
    v_u_phi = np.subtract(np.add(v_0i_rad.values,u_i_rad.values),phi_rad)
    omeg_t_v_u_phi = np.add(omeg_t,v_u_phi)
    ht_res = np.sum(np.multiply(f_A,np.cos(omeg_t_v_u_phi)),axis=1) #not necessary to add A0, since it is already part of the component list
    
    ts_prediction_pd = pd.DataFrame({'values': ht_res}, index=times)
    logger.info('PREDICTION finished')
    
    # add timezone to prediction: first interpret times as tzone of components, then convert to timezone of prediction
    if tzone_convert:
        ts_prediction_pd = ts_prediction_pd.tz_localize(tzone_comp)
        # TODO: workaround to maintain the frequency dropped by tz_localize
        # https://github.com/pandas-dev/pandas/issues/36575
        ts_prediction_pd.index.freq = ts_prediction_pd.index.inferred_freq
        ts_prediction_pd = ts_prediction_pd.tz_convert(tzone_pred)
    
    return ts_prediction_pd


@deprecated_python_option(**DEPRECATED_OPTIONS_PREDICTION_DICT)
def prediction(comp, times=None, timestep=None):
    """
    generates a tidal prediction from a set of components A and phi values.
    The component set has the same timezone as the timeseries used to create it.
    If times is timezone-naive the resulting prediction will be in component timezone.
    If times is timezone-aware the resulting prediction will be converted to that timezone.
    If a components dataframe contains multiple column levels (multiple periods),
    The prediction is a concatenation of predictions of all periods (based on the respective A/phi values).

    Parameters
    ----------
    comp : pd.DataFrame
        The DataFrame contains the component data with component names as index, and colums 'A' and 'phi_deg'.
    times : (pd.DatetimeIndex,slice), optional
        pd.DatetimeIndex with prediction timeseries or slice(tstart,stop,timestep) to construct it from. 
        If None, pd.DatetimeIndex is constructed from the tstart/tstop/timestep metadata attrs of the comp object. 
        Only allowed/relevant for component dataframes with single-level columns (single period). The default is None.
    timestep : str
        Only allowed/relevant for component dataframes with multi-level columns (different periods). 
        The string is parsed with pandas.tseries.frequencies.to_offset(). The default is None.

    Returns
    -------
    ts_prediction : TYPE
        The DataFrame contains a 'values' column and a pd.DatetimeIndex as index, it contains the prediction times and values.

    """
    
    # get settings from component attribute and validate their values
    settings_kwargs = {}
    for setting in ['nodalfactors', 'xfac', 'fu_alltimes', 'source']:
        settings_kwargs[setting] = comp.attrs[setting]
    hatyan_settings = HatyanSettings(**settings_kwargs)
    
    logger.info(f'PREDICTION initializing\n{hatyan_settings}')
    
    metadata_comp = metadata_from_obj(comp)
    tzone_comp = metadata_comp.pop('tzone', None)
    
    if hasattr(comp.columns,"levels"):
        logger.info('prediction() per period due to levels in component dataframe columns')
        if timestep is None:
            raise TypeError("prediction() per period, so 'timestep' argument should not be None")
        if times is not None:
            raise TypeError("prediction() per period, so 'times' argument not allowed")
        # convert timestep to tstep of proper type
        tstep = pd.tseries.frequencies.to_offset(timestep)
        
        ts_periods_dt = comp.columns.levels[1]
        ts_periods_strlist = [str(x) for x in ts_periods_dt]
        
        ts_prediction_perperiod_list = []
        for period_dt in ts_periods_dt:
            logger.info(f'generating prediction {period_dt} of sequence {ts_periods_strlist}')
            comp_oneyear = comp.loc[:,(slice(None),period_dt)]
            comp_oneyear.columns = comp_oneyear.columns.droplevel(1)
            if period_dt.freqstr in ['A-DEC','Y-DEC']: #year frequency
                tstart = pd.Timestamp(period_dt.year, 1, 1)
                tstop = pd.Timestamp(period_dt.year+1, 1, 1) - pd.Timedelta(tstep)
            elif period_dt.freqstr in ['M']: #month frequency
                tstart = period_dt.to_timestamp()
                tstop = period_dt.to_timestamp() + pd.Timedelta(days=period_dt.days_in_month) - pd.Timedelta(tstep)
            else:
                raise Exception(f'unknown freqstr: {period_dt.freqstr}')
            # generate date range and do prediction
            times_pred = pd.date_range(start=tstart, end=tstop, freq=tstep, unit="us", tz=tzone_comp)
            ts_prediction_oneperiod = prediction_singleperiod(comp=comp_oneyear, times=times_pred, hatyan_settings=hatyan_settings)
            ts_prediction_perperiod_list.append(ts_prediction_oneperiod)
        ts_prediction = pd.concat(ts_prediction_perperiod_list)
    else:
        logger.info('prediction() atonce')
        if timestep is not None:
            raise TypeError("prediction() atonce, so 'timestep' argument not allowed")
        if times is None:
            raise TypeError("prediction() atonce, so 'times' argument should not be None")
        if isinstance(times,slice):
            tstart = pd.Timestamp(times.start)
            tstop = pd.Timestamp(times.stop)
            tstep = pd.tseries.frequencies.to_offset(times.step)
            times = pd.date_range(start=tstart, end=tstop, freq=tstep, unit="us")
        if not isinstance(times, pd.DatetimeIndex):
            raise TypeError(f'times argument can be of type pd.DatetimeIndex or slice, not {type(times)}')

        # backwards compatibility for timezone-naive times
        tzone_pred = times.tz
        if tzone_pred is None and tzone_comp is not None:
            times = times.tz_localize(tzone_comp)
            # TODO: workaround to maintain the frequency dropped by tz_localize
            # https://github.com/pandas-dev/pandas/issues/36575
            times.freq = times.inferred_freq
            logger.warning("provided times are timezone-naive and provided components are "
                            "timezone-aware. The times are being interpreted as if they would "
                            f"have the same timezone as the components: {tzone_comp}")
        
        ts_prediction = prediction_singleperiod(comp=comp, times=times, hatyan_settings=hatyan_settings)
    
    # add metadata (and update grootheid)
    metadata_comp = metadata_from_obj(comp)
    if 'grootheid' in metadata_comp:
        # update metadata
        if metadata_comp['grootheid'] == 'WATHTE':
            metadata_comp['grootheid'] = 'WATHTBRKD'
    # prevent adding time metadata from component dataframe to prediction dataframe
    if 'tzone' in metadata_comp.keys():
        metadata_comp.pop('tzone')
    if 'tstart' in metadata_comp.keys():
        metadata_comp.pop("tstart")
    if 'tstop' in metadata_comp.keys():
        metadata_comp.pop("tstop")
    ts_prediction = metadata_add_to_obj(ts_prediction, metadata_comp)

    return ts_prediction
