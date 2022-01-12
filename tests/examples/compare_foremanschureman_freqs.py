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

dood_date = pd.DatetimeIndex([dt.datetime(2014,1,1)])

freqs_pd_hat = hatyan.get_hatyan_freqs(const_list='all',dood_date=dood_date)

v_0i_rad_harmonic_pd = hatyan.get_foreman_v0freq_fromfromharmonicdood(dood_date=dood_date, mode=None) #list with only harmonic components with more precision than file
foreman_shallowrelations = hatyan.get_foreman_shallowrelations()
const_list_foreman = v_0i_rad_harmonic_pd.index.tolist() + foreman_shallowrelations.index.tolist()

#TODO: add translate table for different naming conventions >> or add duplicates to foreman/schureman table?
freqs_pd_for = hatyan.get_foreman_v0_freq(const_list=const_list_foreman,dood_date=dood_date)

freqs_pd = pd.concat([freqs_pd_hat['freq'],freqs_pd_for[1]],axis=1)
freqs_pd['const'] = freqs_pd.index
freqs_pd['bigdiff'] = (np.abs(freqs_pd.iloc[:,0]-freqs_pd.iloc[:,1])>10e-7)
freqs_pd['nan'] = (np.isnan(freqs_pd.iloc[:,0]-freqs_pd.iloc[:,1]))

print(freqs_pd.sort_values('bigdiff').tail(50))
