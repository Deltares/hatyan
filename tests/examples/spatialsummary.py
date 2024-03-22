# -*- coding: utf-8 -*-
"""
spatialsummary_test.py
provides a spatial summary of several values retrieved from analysis results (A, phi, A0, A-factor, phi-difference)
it plots the values in colored dots on a map

"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')
import hatyan

#dir_testdata = 'P:\\1209447-kpp-hydraulicaprogrammatuur\\hatyan\\hatyan_data_acceptancetests'
dir_testdata = 'C:\\DATA\\hatyan_data_acceptancetests'

stats_all = ['ABDN','AMLAHVN','BAALHK','BATH','BERGSDSWT','BORSSLE','BOURNMH','BRESKS','BROUWHVSGT02','BROUWHVSGT08','CADZD','CROMR','CUXHVN','DELFZL','DENHDR','DENOVBTN','DEVPT','DORDT','DOVR','EEMHVN','EEMSHVN','EURPFM','EURPHVN','FELSWE','FISHGD','GEULHVN','GOIDSOD','GOUDBG','HAGSBNDN','HANSWT','HARLGN','HARMSBG','HARTBG','HARTHVN','HARVT10','HEESBN','HELLVSS','HOEKVHLD','HOLWD','HUIBGT','IJMDBTHVN','IJMDSMPL','IMMHM','KATSBTN','KEIZVR','KINLBVE','KORNWDZBTN','KRAMMSZWT','KRIMPADIJSL','KRIMPADLK','K13APFM','LAUWOG','LEITH','LICHTELGRE','LITHDP','LLANDNO','LOWST','MAASMSMPL','MAASSS','MAESLKRZZDE','MARLGT','MOERDK','NES','NEWHVN','NEWLN','NIEUWSTZL','NORTHSS','OOSTSDE04','OOSTSDE11','OOSTSDE14','OUDSD','OVLVHWT','PARKSS','PETTZD','PORTSMH','RAKND','ROOMPBNN','ROOMPBTN','ROTTDM','ROZBSSNZDE','ROZBSSZZDE','SCHAARVDND','SCHEVNGN','SCHIERMNOG','SCHOONHVN','SHEERNS','SINTANLHVSGR','SPIJKNSE','STAVNSE','STELLDBTN','STORNWY','SUURHBNZDE','TENNSHVN','TERNZN','TERSLNZE','TEXNZE','VLAARDGN','VLAKTVDRN','VLIELHVN','VLISSGN','VURN','WALSODN','WERKDBTN','WESTKPLE','WESTTSLG','WEYMH','WHITBY','WICK','WIERMGDN','YERSKE','ZALTBML','A12','AUKFPFM','AWGPFM','D15','F16','F3PFM','J6','K14PFM','L9PFM','NORTHCMRT','Q1']
stats_all.remove('ZALTBML')
stats_subset = ['ABDN','AUKFPFM','CADZD','CUXHVN','D15','DEVPT','DORDT','DOVR','EURPHVN','F3PFM','K14PFM','VLISSGN']

selected_stations = stats_all

case_list = ['A0','P1','K1','M2','S2','K2','M4','P1_K1','NU2_N2','LABDA2_2MN2','K2_S2','T2_S2']

#ldb file
file_ldb = os.path.join(dir_testdata,'other','wvs_coastline3.ldb') #WGS84 ldb is converted to RD, but does not change anything wrt to matlab converted ldb, which is good
ldb_pd_wgs = pd.read_csv(file_ldb, sep="\s+",skiprows=4,names=['x','y'],na_values=[999.999])
x_out, y_out = hatyan.convert_coordinates(coordx_in=ldb_pd_wgs['x'].values, coordy_in=ldb_pd_wgs['y'].values, epsg_in=4326, epsg_out=28992)
ldb_pd = pd.DataFrame({'RDx':x_out/1000, 'RDy':y_out/1000})


#get coordinates
stats_x = []
stats_y = []
for current_station in selected_stations:
    file_data_comp = os.path.join(dir_testdata,'predictie2019','%s_ana.txt'%(current_station))
    file_data_pred = os.path.join(dir_testdata,'predictie2019','%s_pre.txt'%(current_station)) #to get station coordinates
    
    diablocks_pd_extra = hatyan.get_diablocks(filename=file_data_pred)
    stat_x, stat_y = hatyan.convert_coordinates(coordx_in=diablocks_pd_extra.loc[0,'x'], coordy_in=diablocks_pd_extra.loc[0,'y'], epsg_in=diablocks_pd_extra.loc[0,'epsg'], epsg_out=28992)
    stats_x.append(stat_x/1000) #from m to km
    stats_y.append(stat_y/1000) #from m to km

#get data and plot
for case in case_list:
    print('-'*50)
    print('%-45s = %s'%('case',case))
    print('-'*5)
    if '_' in case:
        const_list = case.split('_')
    else:
        const_list = [case]
    
    A_list = []
    phi_list = []
    for current_station in selected_stations:
        file_data_comp = os.path.join(dir_testdata,'predictie2019','%s_ana.txt'%(current_station))
        
        COMP_fromfile = hatyan.read_components(filename=file_data_comp)
        const_list_sel = [x for x in const_list if x in COMP_fromfile.index] #should not be necessary, but T2 is not always available, this causes wrong results for specific stations
        COMP_sel = COMP_fromfile.loc[const_list_sel,:]
        if len(COMP_sel) == 1:
            A_list.append(COMP_sel.loc[const_list_sel[0],'A'])
            phi_list.append(COMP_sel.loc[const_list_sel[0],'phi_deg']%360)
        else:
            A_list.append(COMP_sel.loc[const_list_sel[0],'A']/COMP_sel.loc[const_list_sel[1],'A'])
            phi_list.append(COMP_sel.loc[const_list_sel[0],'phi_deg']-COMP_sel.loc[const_list_sel[1],'phi_deg'])
        
    
    if len(COMP_sel) == 1:
        if case=='A0':
            v_val_A = [-0.05,0.5]
        else:
            v_val_A = [0,np.max(A_list)]
        v_val_phi = [0,360]
        cmmap_phi = 'hsv'
    else:
        v_val_A = [0,np.max(A_list)]
        v_val_phi = [-20,20]
        cmmap_phi = 'jet'
        
    fig, (ax1, ax2) = plt.subplots(1,2,figsize=(14,7), sharex=True, sharey=True)
    if len(COMP_sel) == 1:
        ax1.set_title('amplitude: %s'%(const_list_sel[0]))
        ax2.set_title('phase: %s'%(const_list_sel[0]))
    else:
        ax1.set_title('amplitude factor: %s/%s'%(const_list_sel[0],const_list_sel[1]))
        ax2.set_title('phase difference: %s-%s'%(const_list_sel[0],const_list_sel[1]))
    ax1.plot(ldb_pd['RDx'],ldb_pd['RDy'],'-k',linewidth=0.4)
    ax2.plot(ldb_pd['RDx'],ldb_pd['RDy'],'-k',linewidth=0.4)
    

    ax1_sc = ax1.scatter(stats_x, stats_y, 20, c=A_list, vmin=v_val_A[0], vmax=v_val_A[1], cmap='jet', edgecolors='face')
    ax2_sc = ax2.scatter(stats_x, stats_y, 20, c=phi_list, vmin=v_val_phi[0], vmax=v_val_phi[1], cmap=cmmap_phi, edgecolors='face')
    ax1_cbar = fig.colorbar(ax1_sc, ax=ax1)
    ax2_cbar = fig.colorbar(ax2_sc, ax=ax2, ticks=np.linspace(v_val_phi[0],v_val_phi[1],9))

    ax1_cbar.ax.set_ylabel('amplitude [m]')
    ax2_cbar.ax.set_ylabel('phase [deg]')
    for ax in (ax1,ax2):
        ax.set_xlabel('RD x [km]')
        ax.set_ylabel('RD y [km]')
        ax.grid()
        ax.set_xlim([-400,360])
        ax.set_ylim([200,1200])
        ax.set_aspect('equal', adjustable='box')
    fig.tight_layout()
    fig.savefig('summary_%s.png'%(case))
