# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 14:20:24 2022

@author: veenstra
"""
import pandas as pd
import datetime as dt
import numpy as np
import hatyan

shallow_eqs_pd_for = hatyan.get_foreman_shallowrelations(pd_series=True)
shallow_eqs_pd_hat = hatyan.get_hatyan_shallowrelations()
pd.concat([shallow_eqs_pd_hat,shallow_eqs_pd_for],axis=1)

freqs_pd_hat = hatyan.get_hatyan_freqs(const_list='all',dood_date=pd.DatetimeIndex([dt.datetime(2014,1,1)]))

const_list = hatyan.get_const_list_hatyan('all')
for const_notav in ['SA_IHO1','SA_IHO2','MFM','MSQM','SIGMA1','RO1',
                    'M1B','M1B_IHO1','M1B_IHO2','M1D','M1A','M1','M1_IHO1','M1_IHO2','M1_IHO3',
                    'S1_IHO1','S1_IHO2','S1_IHO3','K1_IHO1','K1_IHO2','FI1','THETA1','MA2','MB2','L2A','L2B']:
    const_list.remove(const_notav)
freqs_pd_for = hatyan.get_foreman_v0_freq(const_list=const_list,dood_date=pd.DatetimeIndex([dt.datetime(2014,1,1)]))

freqs_pd = pd.concat([freqs_pd_hat['freq'],freqs_pd_for[1]],axis=1)
freqs_pd['const'] = freqs_pd.index
freqs_pd['bigdiff'] = (np.abs(freqs_pd.iloc[:,0]-freqs_pd.iloc[:,1])>10e-6) | (np.isnan(freqs_pd.iloc[:,0]-freqs_pd.iloc[:,1]))

