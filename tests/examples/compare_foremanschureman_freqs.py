# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 14:20:24 2022

@author: veenstra
"""
import pandas as pd
import datetime as dt
import numpy as np
import hatyan


dood_date = pd.DatetimeIndex([dt.datetime(2014,1,1)])

freqs_pd_schu = hatyan.get_schureman_freqs(const_list='all',dood_date=dood_date)
v0_pd_schu = hatyan.get_schureman_v0(const_list=freqs_pd_schu.index,dood_date=dood_date)

v_0i_rad_harmonic_pd = hatyan.get_foreman_v0freq_fromfromharmonicdood(dood_date=dood_date, mode=None) #list with only harmonic components with more precision than file
foreman_shallowrelations = hatyan.get_foreman_shallowrelations()
const_list_foreman = v_0i_rad_harmonic_pd.index.tolist() + foreman_shallowrelations.index.tolist()

#TODO: add translate table for different naming conventions >> or add duplicates to foreman/schureman table?
v0_pd_for,freqs_pd_for = hatyan.get_foreman_v0_freq(const_list=const_list_foreman,dood_date=dood_date,)

freqs_pd = pd.concat([freqs_pd_schu['freq'],freqs_pd_for],axis=1)
freqs_pd['const'] = freqs_pd.index
freqs_pd['freq_absdiff'] = np.abs(freqs_pd.iloc[:,0]-freqs_pd.iloc[:,1])
freqs_pd['freq_bigdiff'] = freqs_pd['freq_absdiff']>10e-9
#freqs_pd['freq_nan'] = (np.isnan(freqs_pd.iloc[:,0]-freqs_pd.iloc[:,1]))

v0_pd = pd.concat([v0_pd_schu[0],v0_pd_for[0]],axis=1)
v0_pd['v0_absdiff'] = (np.abs(v0_pd.iloc[:,0]-v0_pd.iloc[:,1])+0.5*np.pi)%np.pi-+0.5*np.pi
v0_pd['v0_bigdiff'] = v0_pd['v0_absdiff']>10e-2 #in radians

print(freqs_pd.sort_values('freq_bigdiff').tail(20))
print()
print(v0_pd.sort_values('v0_bigdiff').tail(20))
