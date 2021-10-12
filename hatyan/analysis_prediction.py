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
    
    import numpy as np
    
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


def get_components_from_ts(ts, const_list, nodalfactors=True, xfac=False, fu_alltimes=True, CS_comps=None, analysis_peryear=False, analysis_permonth=False, return_allyears=False, source='schureman'):
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
    nodalfactors : bool/int, optional
        Whether or not to apply nodal factors. The default is True.
    xfac : bool/int, optional
        Whether or not to apply x-factors. The default is False.
    fu_alltimes : bool/int, optional
        determines whether to calculate nodal factors in middle of analysis period (default) or on every timestep. The default is True.
    analysis_peryear : bool/int, optional
        DESCRIPTION. The default is False.
    analysis_permonth : bool/int, optional
        caution, it tries to analyse each month, but skips it if it fails. analysis_peryear argument has priority. The default is False.
    return_allyears : bool/int, optional
        DESCRIPTION. The default is False.
    CS_comps : pandas.DataFrame, optional
        contains the from/derive component lists for components splitting, as well as the amplitude factor and the increase in degrees. The default is None.

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
    ts_pd = ts
    
    import pandas as pd
    import numpy as np
    
    from hatyan.analysis_prediction import analysis
    from hatyan.analysis_prediction import vectoravg
    from hatyan.hatyan_core import get_const_list_hatyan
       
    print('-'*100)
    print('running: get_components_from_ts')
    
    if type(const_list) is str:
        const_list = get_const_list_hatyan(const_list)
    elif type(const_list) is not list:
        const_list = const_list.tolist()
    n_const = len(const_list)

    if analysis_peryear or analysis_permonth:
        if analysis_peryear:
            print('analysis_peryear=True, separate years are automatically determined from unique calendar years in timeseries')
            ts_years_dt = ts_pd.index.year.unique()
            ts_years = ts_pd.index.year.unique()
        else:
            print('analysis_permonth=True, separate month/year combinations are automatically determined from unique calendar months/years in timeseries')
            ts_years_dt = pd.date_range(start=ts_pd.index.iloc[0], end=ts_pd.index.iloc[-1], freq='M')
            ts_years = ['%d-%02d'%(x.year,x.month) for x in ts_years_dt]

        n_years = len(ts_years)
        if CS_comps is None:
            A_i_all = np.zeros((n_const,n_years))
            phi_i_deg_all = np.zeros((n_const,n_years))
        else:
            A_i_all = np.zeros((n_const+len(CS_comps),n_years))
            phi_i_deg_all = np.zeros((n_const+len(CS_comps),n_years))
        for iY, year_dt in enumerate(ts_years_dt):
            if analysis_peryear:
                print('analyzing %d of sequence %s'%(year_dt,ts_years))
                ts_oneyear_pd = ts_pd[ts_pd.index.year==year_dt]
                COMP_one = analysis(ts_oneyear_pd, const_list=const_list, nodalfactors=nodalfactors, xfac=xfac, fu_alltimes=fu_alltimes, CS_comps=CS_comps, source=source)
                A_i_all[:,iY] = COMP_one.loc[:,'A']
                phi_i_deg_all[:,iY] = COMP_one.loc[:,'phi_deg']
            else:
                print('analyzing %d-%02d of sequence [%s]'%(year_dt.year, year_dt.month, ', '.join(ts_years)))
                ts_oneyear_pd = ts_pd[(ts_pd.index.dt.year==year_dt.year) & (ts_pd.index.dt.month==year_dt.month)]
                try:
                    COMP_one = analysis(ts_oneyear_pd, const_list=const_list, nodalfactors=nodalfactors, xfac=xfac, fu_alltimes=fu_alltimes, CS_comps=CS_comps, source=source)
                    A_i_all[:,iY] = COMP_one.loc[:,'A']
                    phi_i_deg_all[:,iY] = COMP_one.loc[:,'phi_deg']
                except Exception as e:
                    print('WARNING: analysis of %d-%02d failed, error message: "%s'%(year_dt.year,year_dt.month,e))
        
        COMP_all_pd = pd.DataFrame(data=np.hstack([A_i_all,phi_i_deg_all]), columns=pd.MultiIndex.from_product([['A','phi_deg'],ts_years]), index=COMP_one.index)
        print('vector averaging analysis results')
        A_i_mean, phi_i_deg_mean = vectoravg(A_all=A_i_all, phi_deg_all=phi_i_deg_all)
        COMP_mean_pd = pd.DataFrame({ 'A': A_i_mean, 'phi_deg': phi_i_deg_mean},index=COMP_one.index)

    else: #dummy values, COMP_years should be equal to COMP_mean
        COMP_mean_pd = analysis(ts_pd, const_list=const_list, nodalfactors=nodalfactors, xfac=xfac, fu_alltimes=fu_alltimes, CS_comps=CS_comps, source=source)
        COMP_all_pd = None
    
    if return_allyears:
        return COMP_mean_pd, COMP_all_pd
    else:
        return COMP_mean_pd


def analysis(ts, const_list, nodalfactors=True, xfac=False, fu_alltimes=True, CS_comps=None, return_prediction=False, source='schureman'):
    """
    harmonic analysis with matrix transformations (least squares fit), optionally with component splitting
    for details about arguments and return variables, see get_components_from_ts() definition
    
    return_prediction : bool/int, optional
        Whether to generate a prediction for the ts time array. The default is False.
    """
    
    import numpy as np
    import pandas as pd
    import datetime as dt
    
    from hatyan.hatyan_core import get_hatyan_freqs, get_hatyan_v0, get_hatyan_u, get_hatyan_f, get_const_list_hatyan, robust_timedelta_sec
    from hatyan.foreman_core import get_foreman_v0_freq, get_foreman_nodalfactors
    from hatyan.timeseries import check_ts
    
    #drop duplicate times
    ts_pd = ts[~ts.index.duplicated(keep='first')]
    print('-'*100)
    print('ANALYSIS initializing')
    print('%-20s = %s'%('nodalfactors',nodalfactors))
    print('%-20s = %s'%('xfac',xfac))
    print('%-20s = %s'%('fu_alltimes',fu_alltimes))
    if CS_comps is not None:
        print('%-20s = %s'%('CS_comps derive',CS_comps['CS_comps_derive'].tolist()))
        print('%-20s = %s'%('CS_comps from',CS_comps['CS_comps_from'].tolist()))
    else:
        print('%-20s = %s'%('CS_comps', CS_comps))
        
    if len(ts_pd) != len(ts):
        print('WARNING: %i duplicate times of the input timeseries were dropped prior to the analysis'%(len(ts)-len(ts_pd)))
    
    if type(const_list) is str:
        const_list = get_const_list_hatyan(const_list)
    elif type(const_list) is not list:
        const_list = const_list.tolist()
    
    print('%-20s = %s'%('components analyzed',len(const_list)))
    print('%-20s = %s'%('#timesteps',len(ts)))
    print('%-20s = %s'%('tstart',ts.index[0].strftime('%Y-%m-%d %H:%M:%S')))
    print('%-20s = %s'%('tstop',ts.index[-1].strftime('%Y-%m-%d %H:%M:%S')))
    if hasattr(ts.index,'freq'):
        print('%-20s = %s'%('timestep',ts.index.freq))
    
    #check for duplicate components (results in singular matrix)
    if len(const_list) != len(np.unique(const_list)):
        const_list_uniq, const_list_uniq_counts = np.unique(const_list,return_counts=True)
        bool_nonuniq = const_list_uniq_counts>1
        const_list_dupl = pd.DataFrame({'constituent':const_list_uniq[bool_nonuniq],'occurences':const_list_uniq_counts[bool_nonuniq]})
        raise Exception('remove duplicate constituents from const_list:\n%s'%(const_list_dupl))
    
    dood_date_mid = pd.Index([ts_pd.index[len(ts_pd.index)//2]]) #middle of analysis period (2july in case of 1jan-1jan), zoals bij hatyan
    dood_date_start = ts_pd.index[:1] #first date (for v0, also freq?)
    
    ts_pd_nonan = ts_pd[~ts_pd['values'].isna()]
    times_pred_all_pdDTI = pd.DatetimeIndex(ts_pd_nonan.index)
    percentage_nan = 100-len(ts_pd_nonan['values'])/len(ts_pd['values'])*100
    print('percentage_nan in values_meas_sel: %.2f%%'%(percentage_nan))

    #retrieve const_list and frequency in correct order, then retrieve again but now with sorted const_list
    if source.lower()=='schureman':
        t_const_freq_pd = get_hatyan_freqs(const_list)
        const_list = t_const_freq_pd.index.tolist()
        t_const_freq_pd, t_const_speed_all = get_hatyan_freqs(const_list, dood_date=dood_date_mid, return_allraw=True)
    elif source.lower()=='foreman':
        dummy, t_const_freq_pd = get_foreman_v0_freq(const_list=const_list, dood_date=dood_date_start)
        t_const_speed_all = t_const_freq_pd['freq'].values[:,np.newaxis]*(2*np.pi)
        print('WARNING: foreman does not support retrieval of sorted freq list yet, rayleigh check can be wrong')
    else:
        raise Exception('invalid source value (schureman or foreman)')
    t_const_freq = t_const_freq_pd['freq']
            
    #check Rayleigh
    freq_diffs = np.diff(t_const_freq)
    rayleigh_tresh = 0.99
    rayleigh = len(ts_pd['values'])*freq_diffs
    freq_diff_min = rayleigh_tresh/len(ts_pd['values'])
    rayleigh_bool = rayleigh>rayleigh_tresh
    rayleigh_bool_id = np.where(~rayleigh_bool)[0]
    
    if rayleigh_bool.all():
        print('Rayleigh criterion OK (always>%.2f, minimum is %.2f)'%(rayleigh_tresh, np.min(rayleigh)))
        print('Frequencies are far enough apart (always >%.6f, minimum is %.6f)'%(freq_diff_min,np.min(freq_diffs)))
    else:
        print('Rayleigh criterion vandalised (not always>%.2f, minimum is %.2f)'%(rayleigh_tresh, np.min(rayleigh)))
        print('Frequencies with not enough difference (not always >%.6f, minimum is %.6f)'%(freq_diff_min,np.min(freq_diffs)))
        for ray_id in rayleigh_bool_id:
            print(t_const_freq.iloc[[ray_id,ray_id+1]])
    
    #get f and v, only needed after matrix calculations
    if source.lower()=='schureman':
        print('v0 is calculated for start of period: %s'%(dood_date_start[0]))
        v_0i_rad = get_hatyan_v0(const_list, dood_date_start).T #at start of timeseries
    elif source.lower()=='foreman':
        print('v0 is calculated for start of period: %s'%(dood_date_start[0]))
        v_0i_rad, dummy = get_foreman_v0_freq(const_list=const_list, dood_date=dood_date_start)
        v_0i_rad = v_0i_rad.T
    else:
        raise Exception('invalid source value (schureman or foreman)')
    
    if nodalfactors:
        if fu_alltimes:
            print('nodal factors (f and u) are calculated for all timesteps')
            dood_date_fu = times_pred_all_pdDTI
        else:
            print('nodal factors (fu) are calculated for center of period: %s'%(dood_date_mid[0]))
            dood_date_fu = dood_date_mid
        if source.lower()=='schureman':
            f_i = get_hatyan_f(xfac=xfac, const_list=const_list, dood_date=dood_date_fu).T
            u_i_rad = get_hatyan_u(const_list=const_list, dood_date=dood_date_fu).T
        elif source.lower()=='foreman':
            f_i, u_i_rad = get_foreman_nodalfactors(const_list=const_list, dood_date=dood_date_fu)
            f_i, u_i_rad = f_i.T, u_i_rad.T
    else:
        print('no nodal factors (fu) are calculated for (f=1, u=0)')
        f_i = pd.DataFrame(np.ones(len(const_list)),index=const_list).T
        u_i_rad = pd.DataFrame(np.zeros(len(const_list)),index=const_list).T

    v_u = np.add(v_0i_rad.values,u_i_rad.values)
    
    #### TIMESERIES ANALYSIS
    N = len(const_list)
    print('ANALYSIS start (for %i constituents)'%(N))
    
    #times_from0_s = (pd.DatetimeIndex(ts_pd_nonan.index)-dood_date_start[0]).total_seconds().values
    times_from0_s, fancy_pddt = robust_timedelta_sec(ts_pd_nonan.index,refdate_dt=dood_date_start[0])
    times_from0_s = np.transpose(times_from0_s[np.newaxis])
    
    m = len(ts_pd_nonan['values'])
    
    # get xmat and make dot product
    xmat = np.zeros((m,2*N))
    omega_i_rads = t_const_speed_all.T/3600 #angular frequency, 2pi/T, in rad/s, https://en.wikipedia.org/wiki/Angular_frequency (2*np.pi)/(1/x*3600) = 2*np.pi*x/3600
    
    xmat[:,:N] = np.multiply(f_i.values,np.cos(np.multiply(omega_i_rads,times_from0_s)+v_u))
    xmat[:,N:] = np.multiply(f_i.values,np.sin(np.multiply(omega_i_rads,times_from0_s)+v_u))
    xmat_len = xmat.shape[1]
    
    xTmat = xmat.T
    print('calculating xTx matrix')
    tic = dt.datetime.now()
    xTxmat = np.dot(xTmat,xmat)
    print('xTx matrix calculated')
    if 'A0' in const_list: #correct center value for better matrix condition
        xTxmat_condition = np.linalg.cond(xTxmat)
        print('condition of xTx matrix before center adjustment for A0: %.2f'%(xTxmat_condition))
        xTxmat[xmat_len//2,xmat_len//2] = m
    xTxmat_condition = np.linalg.cond(xTxmat)
    print('condition of xTx matrix: %.2f'%(xTxmat_condition))
    if xTxmat_condition > 10:#100: #random treshold
        raise Exception('ERROR: condition of xTx matrix is too high (%.2f), check your timeseries length, try different (shorter) component set or componentsplitting.\nAnalysed %s'%(xTxmat_condition, check_ts(ts_pd)))
    xTymat = np.dot(xTmat,ts_pd_nonan['values'].values)
    
    #solve matrix to get beta_roof_mat (and thus a, b)
    beta_roof_mat = np.linalg.solve(xTxmat,xTymat)
    toc = dt.datetime.now()-tic
    print('matrix system solved, elapsed time: %s'%(toc))

    arctan_ab = np.arctan2(beta_roof_mat[N:],beta_roof_mat[:N]) #(a,b)
    phi_i_rad = arctan_ab
    
    sqsqrt_ab = np.sqrt(np.add(beta_roof_mat[N:]**2,beta_roof_mat[:N]**2)) #(a,b)
    A_i = sqsqrt_ab.flatten()
    phi_i_deg_str = np.rad2deg(phi_i_rad.flatten()%(2*np.pi))
    
    COMP_pd = pd.DataFrame({'A': A_i, 'phi_deg': phi_i_deg_str}, index=const_list)
    if 'A0' in COMP_pd.index: #correct 180 degrees A0 phase by making amplitude value negative
        if COMP_pd.loc['A0','phi_deg']==180:
            COMP_pd.loc['A0','A'] = -COMP_pd.loc['A0','A']
            COMP_pd.loc['A0','phi_deg'] = 0
    
    if CS_comps is not None:
        COMP_pd = split_components(comp=COMP_pd, CS_comps=CS_comps, dood_date_mid=dood_date_mid, xfac=xfac)
        
    print('ANALYSIS finished')
    
    if return_prediction:
        print('immediately generating a prediction for the same time array as the input ts')
        ts_prediction = prediction(comp=COMP_pd, times_pred_all=ts_pd.index, nodalfactors=nodalfactors, xfac=xfac, fu_alltimes=fu_alltimes, source=source)
        return COMP_pd, ts_prediction
    else:
        return COMP_pd


def split_components(comp, CS_comps, dood_date_mid, xfac=False):
    """
    component splitting function
    for details about arguments and return variables, see get_components_from_ts() definition

    """
    
    import numpy as np
    import pandas as pd
    
    from hatyan.hatyan_core import get_hatyan_v0, get_hatyan_u, get_hatyan_f, get_hatyan_freqs

    const_list_inclCS_raw = comp.index.tolist() + CS_comps['CS_comps_derive'].tolist()
    #retrieve const_list and frequency in correct order
    t_const_freq_pd = get_hatyan_freqs(const_list_inclCS_raw)
    const_list_inclCS = t_const_freq_pd.index.tolist()
    #retrieve again but now with sorted const_list
    t_const_freq_pd, t_const_speed_all = get_hatyan_freqs(const_list_inclCS, dood_date=dood_date_mid, return_allraw=True)

    const_list = comp.index.tolist()
    A_i = comp['A'].tolist()
    phi_i_rad_str = np.deg2rad(comp['phi_deg'].tolist())
    
    #if CS_comps is not None: #component splitting
    A_i_inclCS = np.full(shape=(len(const_list_inclCS)),fill_value=np.nan)
    phi_i_rad_str_inclCS = np.full(shape=(len(const_list_inclCS)),fill_value=np.nan)
    CS_v_0i_rad = get_hatyan_v0(const_list=const_list_inclCS, dood_date=dood_date_mid).T.values #with split_components, v0 is calculated on the same timestep as u and f (middle of original series)
    CS_f_i = get_hatyan_f(xfac=xfac, const_list=const_list_inclCS, dood_date=dood_date_mid).T.values
    CS_u_i_rad = get_hatyan_u(const_list=const_list_inclCS, dood_date=dood_date_mid).T.values
    for iC,comp_sel in enumerate(const_list_inclCS):
        if comp_sel in const_list:
            A_i_inclCS[iC] = A_i[const_list.index(comp_sel)]
            phi_i_rad_str_inclCS[iC] = phi_i_rad_str[const_list.index(comp_sel)]
    
    def get_CS_vars(iC_slave, DBETA_in):
        """
        Used to calculate values related to component splitting
        
        """
        #code from resuda.f, line 440 to 455
        DMU = CS_f_i[0,iC_slave]/CS_f_i[0,iC_main]
        DTHETA = ampfac
        DGAMMA = np.deg2rad(degincr)-DBETA_in-(CS_v_0i_rad[0,iC_slave]+CS_u_i_rad[0,iC_slave])+(CS_v_0i_rad[0,iC_main]+CS_u_i_rad[0,iC_main]) #in FORTRAN code, CS_f_i slave/main is also added, this seems wrong
        DREEEL = 1+DMU*DTHETA*np.cos(DGAMMA)
        DIMAGI = DMU*DTHETA*np.sin(DGAMMA)  
        DALPHA = np.sqrt(DREEEL*DREEEL+DIMAGI*DIMAGI)
        if DALPHA < 1e-50:
            raise Exception('ERROR: DALPHA too small, component splitting failed?')
        DBETA = np.arctan2(DIMAGI,DREEEL)
        if np.sign(DIMAGI) == np.sign(DREEEL):
            DBETA = DBETA
        else:
            DBETA = DBETA
        return DTHETA, DALPHA, DBETA
    
    if len(np.unique(CS_comps['CS_comps_derive'])) != len(CS_comps['CS_comps_derive']):
        raise Exception('ERROR: CS_comps_derive contains duplicate components')
        
    for comp_main in np.unique(CS_comps['CS_comps_from']):
        main_ids = np.where(CS_comps['CS_comps_from'] == comp_main)[0]
        comp_slave = CS_comps.loc[main_ids,'CS_comps_derive'].tolist()
        print('splitting component %s into %s'%(comp_main, comp_slave))
        
        iC_main = const_list_inclCS.index(comp_main)
        iC_slave = const_list_inclCS.index(comp_slave[0])
        idslave_CScomp = CS_comps['CS_comps_derive'].tolist().index(comp_slave[0])
        degincr = CS_comps['CS_degincrs'].tolist()[idslave_CScomp]
        ampfac = CS_comps['CS_ampfacs'].tolist()[idslave_CScomp]
        
        DTHETA, DALPHA, DBETA = get_CS_vars(iC_slave,0)
        A_i_inclCS[iC_main] = A_i_inclCS[iC_main]/DALPHA
        phi_i_rad_str_inclCS[iC_main] = (phi_i_rad_str_inclCS[iC_main]-DBETA)%(2*np.pi)
        
        if len(comp_slave) == 1:
            A_i_inclCS[iC_slave] = A_i_inclCS[iC_main]*DTHETA
            phi_i_rad_str_inclCS[iC_slave] = (phi_i_rad_str_inclCS[iC_main]+np.deg2rad(degincr))%(2*np.pi)
        elif len(comp_slave) == 2:
            #T2
            iC_slave2 = const_list_inclCS.index(comp_slave[1])
            idslave_CScomp2 = CS_comps['CS_comps_derive'].tolist().index(comp_slave[1])
            degincr = CS_comps['CS_degincrs'].tolist()[idslave_CScomp2]
            ampfac = CS_comps['CS_ampfacs'].tolist()[idslave_CScomp2]
            
            DTHETA, DALPHA, DBETA = get_CS_vars(iC_slave2, DBETA)
            A_i_inclCS[iC_main] = A_i_inclCS[iC_main]/DALPHA
            phi_i_rad_str_inclCS[iC_main] = (phi_i_rad_str_inclCS[iC_main]-DBETA)%(2*np.pi)
            
            A_i_inclCS[iC_slave2] = A_i_inclCS[iC_main]*DTHETA
            phi_i_rad_str_inclCS[iC_slave2] = (phi_i_rad_str_inclCS[iC_main]+np.deg2rad(degincr))%(2*np.pi)
            
            #revert back to K2
            degincr = CS_comps['CS_degincrs'].tolist()[idslave_CScomp]
            ampfac = CS_comps['CS_ampfacs'].tolist()[idslave_CScomp]
            A_i_inclCS[iC_slave] = A_i_inclCS[iC_main]*ampfac
            phi_i_rad_str_inclCS[iC_slave] = (phi_i_rad_str_inclCS[iC_main]+np.deg2rad(degincr))%(2*np.pi)
        else:
            raise Exception('ERROR: length of comp_slave is invalid (%i)'%(len(comp_slave)))
     
    phi_i_deg_str_inclCS = np.rad2deg(phi_i_rad_str_inclCS)
    
    comp_CS = pd.DataFrame({ 'A': A_i_inclCS, 'phi_deg': phi_i_deg_str_inclCS},index=const_list_inclCS)
    
    return comp_CS


def prediction(comp, times_pred_all=None, times_ext=None, timestep_min=None, nodalfactors=True, xfac=False, fu_alltimes=True, source='schureman'):
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
    nodalfactors : bool/int, optional
        Whether or not to apply nodal factors. The default is True.
    xfac : bool/int, optional
        Whether or not to apply x-factors. The default is False.
    fu_alltimes : bool/int, optional
        determines whether to calculate nodal factors in middle of the prediction period (default) or on every timestep. The default is True.

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    ts_prediction_pd : pandas.DataFrame
        The DataFrame should contain a 'values' column and a pd.DatetimeIndex as index, it contains the prediction times and values.
    
    """
    
    print('-'*100)
    print('PREDICTION initializing')
    print('%-20s = %s'%('nodalfactors',nodalfactors))
    print('%-20s = %s'%('xfac',xfac))
    print('%-20s = %s'%('fu_alltimes',fu_alltimes))
    
    import numpy as np
    import pandas as pd
    from packaging import version
    from hatyan.hatyan_core import get_hatyan_freqs, get_hatyan_v0, get_hatyan_u, get_hatyan_f, robust_daterange_fromtimesextfreq
    from hatyan.foreman_core import get_foreman_v0_freq, get_foreman_nodalfactors
    
    COMP = comp.copy()
    
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
    
    if isinstance(times_pred_all, pd.core.indexes.datetimes.DatetimeIndex) or isinstance(times_pred_all, pd.core.indexes.base.Index):
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
    
    #retrieve const_list and frequency in correct order, then retrieve again but now with sorted const_list
    if source.lower()=='schureman':
        t_const_freq_pd = get_hatyan_freqs(COMP.index.tolist())
        const_list = t_const_freq_pd.index.tolist()
        t_const_freq_pd, t_const_speed_all = get_hatyan_freqs(const_list, dood_date=dood_date_mid, return_allraw=True)
    elif source.lower()=='foreman':
        print('v0 is calculated for start of period: %s'%(dood_date_start[0]))
        dummy, t_const_freq_pd = get_foreman_v0_freq(const_list=COMP.index.tolist(), dood_date=dood_date_start)
        const_list = t_const_freq_pd.index.tolist()
        t_const_speed_all = t_const_freq_pd['freq'].values[:,np.newaxis]*(2*np.pi)
        print('WARNING: foreman does not support retrieval of sorted freq list yet')
    else:
        raise Exception('invalid source value (schureman or foreman)')
    COMP['freq'] = t_const_freq_pd['freq']
    COMP = COMP.sort_values(by='freq')
    
    A = np.array(COMP['A'])
    phi_rad = np.array(np.deg2rad(COMP['phi_deg']))
    
    if source.lower()=='schureman':
        print('v0 is calculated for start of period: %s'%(dood_date_start[0]))
        v_0i_rad = get_hatyan_v0(const_list, dood_date_start).T #at start of timeseries
    elif source.lower()=='foreman':
        print('v0 is calculated for start of period: %s'%(dood_date_start[0]))
        v_0i_rad, dummy = get_foreman_v0_freq(const_list=const_list, dood_date=dood_date_start)
        v_0i_rad = v_0i_rad.T
    else:
        raise Exception('invalid source value (schureman or foreman)')
    
    if nodalfactors:
        if fu_alltimes:
            print('nodal factors (f and u) are calculated for all timesteps')
            dood_date_fu = times_pred_all_pdDTI
        else:
            print('nodal factors (fu) are calculated for center of period: %s'%(dood_date_mid[0]))
            dood_date_fu = dood_date_mid
        if source.lower()=='schureman':
            f_i = get_hatyan_f(xfac=xfac, const_list=const_list, dood_date=dood_date_fu).T
            u_i_rad = get_hatyan_u(const_list=const_list, dood_date=dood_date_fu).T
        elif source.lower()=='foreman':
            f_i, u_i_rad = get_foreman_nodalfactors(const_list=const_list, dood_date=dood_date_fu)
            f_i, u_i_rad = f_i.T, u_i_rad.T
    else:
        print('no nodal factors (fu) are calculated for (f=1, u=0)')
        f_i = pd.DataFrame(np.ones(len(const_list)),index=const_list).T
        u_i_rad = pd.DataFrame(np.zeros(len(const_list)),index=const_list).T
    
    print('PREDICTION started')
    omega_i_rads = t_const_speed_all.T/3600 #angular frequency, 2pi/T, in rad/s, https://en.wikipedia.org/wiki/Angular_frequency (2*np.pi)/(1/x*3600) = 2*np.pi*x/3600
    if ~isinstance(times_pred_all_pdDTI,pd.DatetimeIndex) & (version.parse(pd.__version__) >= version.parse('1.2.0')): #fix for non-backwards compatible change in pandas, pandas version 1.1.2 is used for RWS version.
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
