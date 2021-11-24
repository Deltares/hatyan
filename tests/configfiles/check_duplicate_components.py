# -*- coding: utf-8 -*-
"""
Created on Wed Nov 24 18:50:15 2021

@author: veenstra

hatyan is available on Github: https://github.com/Deltares/hatyan/
"""

import hatyan
import numpy as np
import datetime as dt
import pandas as pd

const_list_FES2014 = ['2N2','EPS2','J1','K1','K2','L2','LA2','M2','M3','M4','M6','M8','MF','MKS2','MM','MN4','MS4','MSF','MSQM','MTM','MU2','N2','N4','NU2','O1','P1','Q1','R2','S1','S2','S4','SA','SSA','T2']
const_list_FES2014 = ['2N2','EPS2','J1','K1','K2','L2','LABDA2','M2','M3','M4','M6','M8','MF','MKS2','MM','MN4','MS4','MSF','MSQM','MFM','MU2','N2','N4','NU2','O1','P1','Q1','R2','S1','S2','S4','SA','SSA','T2']
const_list_RWSyear = hatyan.get_const_list_hatyan('year')
const_list_combined = np.unique(const_list_RWSyear+const_list_FES2014)

dood_date = pd.DatetimeIndex([dt.datetime(1900,1,1)]) #dummy value

hat_freqs_FES2014 = hatyan.get_hatyan_freqs(const_list=const_list_FES2014,dood_date=dood_date,sort_onfreq=False)
hat_freqs_RWSyear = hatyan.get_hatyan_freqs(const_list=const_list_RWSyear,dood_date=dood_date,sort_onfreq=False)
hat_freqs_combined = hatyan.get_hatyan_freqs(const_list=const_list_combined,dood_date=dood_date,sort_onfreq=False)
hat_v0 = hatyan.get_hatyan_v0(const_list=const_list_combined,dood_date=dood_date)
hat_u = hatyan.get_hatyan_u(const_list=const_list_combined,dood_date=dood_date)
hat_f = hatyan.get_hatyan_f(const_list=const_list_combined,dood_date=dood_date,xfac=False)


hat_all = pd.DataFrame()#{'freq_FES2014':[hat_freqs_FES2014['freq'].values],'freqs_RWSyear':[hat_freqs_FES2014['freq'].values]},index=hat_freqs.index).sort_values('freqs_RWSyear')
hat_all['freq_combined'] = hat_freqs_combined['freq']
hat_all['freq_dupl'] = hat_all['freq_combined'].duplicated(keep=False)
hat_all['freq_RWSyear'] = hat_freqs_RWSyear['freq']
hat_all['freq_FES2014'] = hat_freqs_FES2014['freq']
hat_all['v0'] = hat_v0
hat_all['u'] = hat_u
hat_all['f'] = hat_f
hat_all = hat_all.sort_values('freq_combined').round(6)



