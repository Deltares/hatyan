# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 23:37:25 2021

@author: veenstra

Data van:
https://waterberichtgeving.rws.nl/wbviewer/map.php?set=scheveningen
Station NZB_N, CSV downloaden via tabel:
https://waterberichtgeving.rws.nl/wbviewer/maak_grafiek.php?loc=NZB_N&set=scheveningen&nummer=3&page_type=tabel

"""

import os
import numpy as np
import datetime as dt
import pandas as pd
import matplotlib.pyplot as plt
plt.close('all')
import hatyan

#dir_testdata = 'P:\\1209447-kpp-hydraulicaprogrammatuur\\hatyan\\hatyan_data_acceptancetests'
dir_testdata = 'C:\\DATA\\hatyan_data_acceptancetests'

const_list = 'springneap' #'year' 'halfyear' 'month' etc

file_stroming = os.path.join(dir_testdata,'other','NZB_N - Stroomsnelheid [cm_s].csv')
data_stroming = pd.read_csv(file_stroming,sep=';',skipfooter=2, engine='python', parse_dates=[0], index_col=0)
times_ext = [data_stroming.index[0],data_stroming.index[-1]+dt.timedelta(days=7)]

meas_mag = data_stroming['stroomsnelheid']
meas_dir = data_stroming['stroomrichting']
meas_x = meas_mag*np.sin(np.deg2rad(meas_dir))
meas_y = meas_mag*np.cos(np.deg2rad(meas_dir))

data_x = pd.DataFrame({'values':meas_x})
data_y = pd.DataFrame({'values':meas_y})
data_x.dummymeta = 'dummy'
data_y.dummymeta = 'dummy'

comp_x = hatyan.analysis(ts=data_x, const_list=const_list)
comp_y = hatyan.analysis(ts=data_y, const_list=const_list)
pred_x = hatyan.prediction(comp=comp_x, times_ext=times_ext, timestep_min=10)
pred_y = hatyan.prediction(comp=comp_y, times_ext=times_ext, timestep_min=10)

pred_mag = np.sqrt(pred_x**2+pred_y**2)
pred_dir = np.rad2deg(np.arctan2(pred_x,pred_y))%360

fig, (ax1,ax2) = plt.subplots(2,1,figsize=(12,6),sharex=True)
ax1.plot(meas_mag,label='measured magnitude')
ax1.plot(pred_mag,label='predicted magnitude')
ax2.plot(meas_dir,label='measured direction')
ax2.plot(pred_dir,label='predicted direction')
ax1.legend()
ax2.legend()
ax1.grid()
ax2.grid()
fig.tight_layout()
fig.savefig('analysis_tidalcurrent.png')
