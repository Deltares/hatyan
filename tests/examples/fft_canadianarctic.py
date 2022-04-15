# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 11:43:31 2022

@author: veenstra
"""

import hatyan
hatyan.close('all')
import pandas as pd

file_predcan = r'p:\1230882-emodnet_hrsm\fromAmey\Regional_canada_model\bathymetry_checks\CHSdata\Churchill.wl'
ts_pred_can = pd.read_csv(file_predcan,comment='#',delim_whitespace=True, names=['times','values'], parse_dates=['times'])
ts_pred_can['times'] = ts_pred_can['times'].dt.tz_localize(None)
ts_pred_can = ts_pred_can.set_index('times')

source='foreman'
hatyan_settings = hatyan.HatyanSettings(nodalfactors=True,return_prediction=True,source=source,fu_alltimes=False)

#fft analysis
peak_freq, hatyan_freqs_suggestions, hatyan_freqs_matches = hatyan.timeseries_fft(ts_pred_can, prominence=.5**1, plot_fft=True, source=source)
const_list = hatyan_freqs_suggestions.index.unique().tolist()+['A0','SA','SSA','H1','H2']

#harmonic analysis and reprediction
comp,reprediction = hatyan.analysis(ts=ts_pred_can,const_list=const_list,hatyan_settings=hatyan_settings)

#plotting timeseries
fig1, (ax1,ax2) = hatyan.plot_timeseries(ts=reprediction, ts_validation=ts_pred_can)
ax2.set_ylim(-0.4,0.4)
#fig2, (ax3,ax4) = hatyan.plot_components(comp=comp)
