# -*- coding: utf-8 -*-
"""
Created on Fri Feb  4 15:14:48 2022

@author: veenstra
"""

import os
import pandas as pd
import datetime as dt
import scipy.io as sio
#import numpy as np
import matplotlib.pyplot as plt
plt.close('all')
import hatyan

file_config = os.path.realpath(__file__) #F9 doesnt work, only F5 (F5 also only method to reload external definition scripts)
dir_output, timer_start = hatyan.init_RWS(file_config, interactive_plots=False)
#dir_testdata = 'P:\\1209447-kpp-hydraulicaprogrammatuur\\hatyan\\hatyan_data_acceptancetests'
dir_testdata = 'C:\\DATA\\hatyan_data_acceptancetests'

points_loop = [0]#,100,200]#range(len(datablock_list))#

for iF in points_loop:
    print('iF: %d'%(iF))
    #HATYAN
    
    consts_FES = ['SA','2N2','LABDA2','MF','MFM','P1','SSA','MNS2','M2','MKS2','MU2','Q1','T2','J1','M3','MM','N2','R2','K1','M4','MN4','N4','S1','K2','M6','MS4','NU2','S2','L2','M8','MSF','O1','S4','A0']
    consts_FES = ['A0', 'SA', 'SSA', 'MM', 'MSF', 'MF', 'Q1', 'O1', 'P1', 'S1', 'K1', 'J1', '2N2', 'MU2', 'N2', 'NU2', 'M2', 'MKS2', 'LABDA2', 'L2', 'T2', 'S2', 'R2', 'K2', 'M3', 'N4', 'MN4', 'M4', 'MS4', 'S4', 'M6', 'M8']
    #no N4
    #consts_FES = ['SA','2N2','LABDA2','MF','MFM','P1','SSA','MNS2','M2','MKS2','MU2','Q1','T2','J1','M3','MM','N2','R2','K1','M4','MN4','S1','K2','M6','MS4','NU2','S2','L2','M8','MSF','O1','S4','A0']
    
    file_tim_noSLR = os.path.join(dir_testdata,'other',f'DCSM-FM_OB_all_20181108_{iF+1:04d}.tim') # r'p:\11205235-014-getijde-energie\GTSM_runs\run_noSLR\nesthd2\DCSM-FM_OB_all_20181108_%04d.tim'%(iF+1)
    data_tim_noSLR = pd.read_csv(file_tim_noSLR, comment='*', delim_whitespace=True, names=['times_min','values'])
    data_tim_noSLR.index = dt.datetime(2024,12,15)+pd.to_timedelta(data_tim_noSLR['times_min'],unit='minute')
    comphat_noSLR = hatyan.analysis(ts=data_tim_noSLR, const_list=consts_FES, nodalfactors=True, fu_alltimes=False,source='schureman')
    #print(comp_example[comp_example['const_list']=='M2'])
    
    #file_tim = r'p:\11205235-014-getijde-energie\GTSM_runs\run_2025_10cm\nesthd2\DCSM-FM_OB_all_20181108_%04d.tim'%(iF+1)
    #data_tim = pd.read_csv(file_tim, comment='*', delim_whitespace=True, names=['times_min','values'])
    #data_tim['times'] = dt.datetime(2024,12,15)+pd.to_timedelta(data_tim['times_min'],unit='minute')
    #comphat_2025 = analysis(ts=data_tim, const_list=consts_FES, nodalfactors=True, fu_alltimes=False)
    
    #Timeseries.plot_timeseries(ts=data_tim_noSLR,ts_validation=data_tim)
    
    #comphat_2025diff = pd.DataFrame({'A':comphat_2025['A']-comphat_noSLR['A'],
    #                                 'phi_deg':comphat_2025['phi_deg']-comphat_noSLR['phi_deg']})

    #TTIDE
    file_const_ttide = os.path.join(dir_testdata,'other','Tide_noSLR.mat') #r'p:\11205235-014-getijde-energie\Nesting\Tide_noSLR.mat'
    data_ttide_raw = sio.loadmat(file_const_ttide, mat_dtype=True)['Tidestruct'][0][0]
    #data_ttide.dtype = dtype([('FREQ', 'O'), ('NAME', 'O'), ('TIDECONamp', 'O'), ('TIDECONampnew', 'O'), ('TIDECONphs', 'O'), ('TIDECONphsnew', 'O'), ('TIDECONdiff_amp', 'O'), ('TIDECONdiff_phs', 'O')])
    comp_ttide_noSLR = pd.DataFrame({'A':data_ttide_raw['TIDECONamp'][:,iF],
                                     'phi_deg':data_ttide_raw['TIDECONphs'][:,iF]},
                                    index=[x.strip() for x in data_ttide_raw['NAME']])
    #data_ttide_noSLR.rename(index={'A0':'H0'},inplace=True)
    comp_ttide_noSLR.rename(index={'LDA2':'LABDA2'},inplace=True)
    #file_const_ttide = r'p:\11205235-014-getijde-energie\Nesting\Tide_10cmSLR.mat'
    #data_ttide_raw = sio.loadmat(file_const_ttide, mat_dtype=True)['Tidestruct'][0][0]
    #data_ttide_2025  = pd.DataFrame({'A':data_ttide_raw['TIDECONamp'][:,iF],
    #                                 'phi_deg':data_ttide_raw['TIDECONphs'][:,iF]},
    #                                index=[x.strip() for x in data_ttide_raw['NAME']])
    #data_ttide_2025.rename(index={'A0':'H0'},inplace=True)
    #data_ttide_2025.rename(index={'LDA2':'LABDA2'},inplace=True)
    #data_ttide_2025_diff = pd.DataFrame({'A':data_ttide_2025['A']-data_ttide_noSLR['A'],
    #                                    'phi_deg':data_ttide_2025['phi_deg']-data_ttide_noSLR['phi_deg']})

    #COMPARE
    #Components.plot_components(comp=comphat_noSLR, comp_validation=data_ttide_noSLR)
    #data_ttide_noSLR_noN4 = data_ttide_noSLR#[data_ttide_noSLR.index!='N4']
    hatyan.plot_components(comp=comphat_noSLR, comp_validation=comp_ttide_noSLR)
    """
    data_ttide_2025_noN4 = data_ttide_2025[data_ttide_2025.index!='N4']
    Components.plot_components(comp=comphat_2025, comp_validation=data_ttide_2025_noN4)
    data_ttide_2025_diff_noN4 = data_ttide_2025_diff[data_ttide_2025_diff.index!='N4']
    Components.plot_components(comp=comphat_2025diff, comp_validation=data_ttide_2025_diff_noN4)
    """

hatyan.exit_RWS(timer_start) #provides footer to outputfile when calling this script with python




