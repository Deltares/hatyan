# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 23:40:17 2021

@author: veenstra
"""

import sys
import datetime as dt
import os
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')
import pandas as pd
import hatyan

file_config = os.path.realpath(__file__) #F9 doesnt work, only F5 (F5 also only method to reload external definition scripts)
dir_output, timer_start = hatyan.init_RWS(file_config, sys.argv, interactive_plots=False) #provides header to outputfile when calling this script with python
#dir_testdata = 'P:\\1209447-kpp-hydraulicaprogrammatuur\\hatyan\\hatyan_data_acceptancetests'
dir_testdata = 'C:\\DATA\\hatyan_data_acceptancetests'

const_list_hatyan195_orig = hatyan.get_const_list_hatyan('all_originalorder')

times_doodsonplot = pd.date_range(start=dt.datetime(2018,12,1),end=dt.datetime(2025,1,1),freq='60min')

const_list_freqv0uf = ['A0','SA','SSA','MSM','MM','MSF','MF','2Q1','Q1','O1','TAU1','M1C','CHI1','PI1','P1','S1','K1','PSI1','J1','OO1','2N2','MU2','N2','NU2','M2','LABDA2','L2','T2','S2','R2','K2','ETA2','M3']
#const_list_freqv0uf = hatyan.get_const_list_hatyan('year')
#const_list_freqv0uf = hatyan.get_const_list_hatyan('springneap')
#const_list_freqv0uf = ['MU2','N2','NU2','M2','2MN2','S2','M4','MS4','M6','2MS6','M8','3MS8'] #xfac list, should also entail constituents for which f==1
#const_list_freqv0uf = ['M2','K1','S2','2MN2'] 
#const_list_freqv0uf = ['MM','MF','Q1','O1','K1','M2','N2','K2','S2','J1','OO1'] #M2=2N2 MU2 NU2 N2 M2. list also for SLS and R
year_v0 = 2006
dood_date_v0freq = pd.date_range(start=dt.datetime(year_v0,1,1,0,0),end=dt.datetime(year_v0,1,1,2,0),freq='10min')
times_freqdiff = pd.date_range(start=dt.datetime(2000,1,1),end=dt.datetime(2040,1,1),freq='1D')

const_list_freqv0uf_hat55 = const_list_freqv0uf[:]
if 'A0' in const_list_freqv0uf_hat55:
    const_list_freqv0uf_hat55.remove('A0')

dood_date_fu = pd.date_range(start=dt.datetime(1980,1,1), end=dt.datetime(2020,1,1), freq='60D')

treshold_freq = 1e-6
treshold_v0uf = 1e-4

dir_output_vuf = os.path.join('vuf_files_peryear')
if not os.path.exists(dir_output_vuf):
    os.makedirs(dir_output_vuf)
dir_output_ufplots = os.path.join('uf_plots_timeseries')
if not os.path.exists(dir_output_ufplots):
    os.makedirs(dir_output_ufplots)
dir_output_v0uplots = os.path.join('v0u_plots_timeseries')
if not os.path.exists(dir_output_v0uplots):
    os.makedirs(dir_output_v0uplots)

print('')
print('#####################################################################')
print('#### GET hatyan55 values ############################################')
print('#####################################################################')


def get_hatyan55_values(file_hatyan55):
    
    #######################
    colnames_freq = ['RANGNR.','NAAM','HOEKSNELHEID[deg/hr]']
    hatyan55_freq_raw = pd.read_csv(file_hatyan55, names=colnames_freq, skiprows=21, nrows=219, delim_whitespace=True, error_bad_lines=False, warn_bad_lines=False)
    bool_droprowsbegin = np.where(hatyan55_freq_raw.iloc[:,1]=='RWS-CIV')[0] #start of this year
    bool_droprowsend = np.where(hatyan55_freq_raw.iloc[:,1]=='GRADEN/UUR')[0] #start of this year
    drop_idx = []
    for x,y in zip(bool_droprowsbegin,bool_droprowsend):
        drop_idx = drop_idx+list(range(x,y+1))
    
    hatyan55_freq = hatyan55_freq_raw.drop(drop_idx).reset_index(drop=True)
    hatyan55_freq = hatyan55_freq.apply(pd.to_numeric,errors='ignore')
    hatyan55_freq = hatyan55_freq.set_index('NAAM')
    hatyan55_freq.index.name=None
    
    #######################
    hatyan55_v0u = pd.DataFrame()
    hatyan55_f = pd.DataFrame()
    colnames_v0uf = ['RANGNR.','VU-FAKTOR','F-FAKTOR']
    hatyan55_v0uf_raw_allyears = pd.read_csv(file_hatyan55, names=colnames_v0uf, skiprows=270, delim_whitespace=True, error_bad_lines=False, warn_bad_lines=False)
    line_startyears = np.where(hatyan55_v0uf_raw_allyears.iloc[:,0]=='JAAR')[0]
    yearnos = hatyan55_v0uf_raw_allyears.iloc[line_startyears,2].astype(int).tolist()
    dood_date_hatyan55 = pd.DatetimeIndex([dt.datetime(x,7,2) for x in yearnos])
    dood_date_hatyan55_v0 = pd.DatetimeIndex([dt.datetime(x,1,1) for x in yearnos])
    for line_no, year in zip(line_startyears, yearnos):
        hatyan55_v0uf_raw_1y = hatyan55_v0uf_raw_allyears.iloc[line_no+3:line_no+222]
        
        bool_droprowsbegin = hatyan55_v0uf_raw_1y[hatyan55_v0uf_raw_1y.iloc[:,1]=='RWS-CIV'].index
        bool_droprowsend = hatyan55_v0uf_raw_1y[hatyan55_v0uf_raw_1y.iloc[:,1]=='GRADEN'].index
        drop_idx = []
        for x,y in zip(bool_droprowsbegin,bool_droprowsend):
            drop_idx = drop_idx+list(range(x,y+1))
        
        hatyan55_v0uf_1y = hatyan55_v0uf_raw_1y.drop(drop_idx).reset_index(drop=True)
        hatyan55_v0uf_1y = hatyan55_v0uf_1y.apply(pd.to_numeric,errors='ignore')
        hatyan55_v0uf_1y.index = hatyan55_freq.index
        
        hatyan55_v0u[year] = hatyan55_v0uf_1y['VU-FAKTOR']
        hatyan55_f[year] = hatyan55_v0uf_1y['F-FAKTOR']
    
    return hatyan55_freq, hatyan55_v0u, hatyan55_f, dood_date_hatyan55, dood_date_hatyan55_v0


file_hatyan55 = os.path.join(dir_testdata,'other','hatyan55_output.txt')
hatyan55_freq, hatyan55_v0u, hatyan55_f, dood_date_hatyan55, dood_date_hatyan55_v0 = get_hatyan55_values(file_hatyan55)

#######################
file_fwithxfac = os.path.join(dir_testdata,'other','vuf2016.txt')
v0uffile_2016 = pd.read_csv(file_fwithxfac, skiprows=3, delim_whitespace=True, names=['Nr.','Hoeksnelh.','V0+u','f(xfac1)','Naam'])
v0uffile_2016 = v0uffile_2016.set_index('Naam')

hatyan55_freq_conv = np.array([0]+(hatyan55_freq.loc[const_list_freqv0uf_hat55,'HOEKSNELHEID[deg/hr]'].values/360).tolist())
v_0i_rad_H_195 = hatyan.get_schureman_v0(const_list_hatyan195_orig, dood_date_hatyan55_v0)
u_i_rad_195 = hatyan.get_schureman_u(const_list=const_list_hatyan195_orig, dood_date=dood_date_hatyan55)
v0u_deg_195 = np.rad2deg(v_0i_rad_H_195+u_i_rad_195)%360
hatyan55py_v0u_diff = ((v0u_deg_195.iloc[1:,:]-hatyan55_v0u.values)+180)%360-180 #first remove A0 and use .values to ignore column names
hatyan55_ucorr = (hatyan55_v0u-np.rad2deg(v_0i_rad_H_195.iloc[1:,:].values)+180)%360-180

print('')
print('#####################################################################')
print('#### READ FOREMAN TABLE TEST ########################################')
print('#####################################################################')

foreman_doodson_harmonic, foreman_nodal_harmonic = hatyan.get_foreman_doodson_nodal_harmonic()
foreman_shallowrelations = hatyan.get_foreman_shallowrelations()

foreman_harmonic_doodson_all_list = foreman_doodson_harmonic.index.tolist()
foreman_harmonic_nodal_all_list = foreman_nodal_harmonic.index.unique().tolist()
foreman_shallowrelations_list = foreman_shallowrelations.index.tolist()
foreman_shallow_dependencies_list = list(set(foreman_shallowrelations.loc[:,[3,5,7,9]].values.reshape(580)))

print( '#### FOREMAN FILE TEST, independent of const_list ############' )
nodal_nodoodson = [x for x in set(foreman_harmonic_nodal_all_list+foreman_shallow_dependencies_list) if x not in foreman_harmonic_doodson_all_list] # doodson is needed. werkt niet want eerste constituent regel wordt altijd naar doodson geschreven..
print( 'nodal_nodoodson, ERROR: provide doodson for:', nodal_nodoodson )
doodson_nonodal = [x for x in foreman_harmonic_doodson_all_list if x not in foreman_harmonic_nodal_all_list] # ok, f=1, u=0
print( 'doodson_nonodal, OK: f=1 and u=0 for:       ', doodson_nonodal )
print( '#### END FOREMAN FILE TEST ###################################' )

pltR_lat_deg = np.arange(-90,90+1,dtype=float) #nan instead of 0, so no division by 0 in R1
pltR_lat_deg[pltR_lat_deg==0] = np.nan
pltR_lat_rad = np.deg2rad(pltR_lat_deg)
R1 = 0.36309*(1.-5.*np.sin(pltR_lat_rad)*np.sin(pltR_lat_rad))/np.sin(pltR_lat_rad)
R2 = 2.59808*np.sin(pltR_lat_rad)

fig, ax1 = plt.subplots(1,1)
ax1.plot(pltR_lat_deg,R1)
ax1.plot(pltR_lat_deg,R2)
ax1.legend(['R1','R2'])
ax1.set_xlabel('latitude in degrees')
ax1.grid()
fig.tight_layout()
fig.savefig('test_foreman_R1R2.png')

print('')
print('#####################################################################')
print('#### CALCULATING AND PLOTTING DOODSON NUMBERS #######################')
print('#####################################################################')

#get and plot doodson values
dood_T_HAT, dood_S_HAT, dood_H_HAT, dood_P_HAT, dood_N_HAT, dood_P1_HAT = hatyan.get_doodson_eqvals(times_doodsonplot)

plt.figure(figsize=(15,5))
plt.title('doodson all')
plt.plot(times_doodsonplot, dood_T_HAT%(2*np.pi),'-',linewidth=1,markersize=0.5, label='T')
plt.plot(times_doodsonplot, dood_S_HAT%(2*np.pi),'-',linewidth=1,markersize=0.5, label='S')#, alpha=0.5)
plt.plot(times_doodsonplot, dood_H_HAT%(2*np.pi),'-',linewidth=1,markersize=0.5, label='H')
plt.plot(times_doodsonplot, dood_P_HAT%(2*np.pi),'-',linewidth=1,markersize=0.5, label='P')
plt.plot(times_doodsonplot, dood_N_HAT%(2*np.pi),'-',linewidth=1,markersize=0.5, label='N')
plt.plot(times_doodsonplot, dood_P1_HAT%(2*np.pi),'-',linewidth=1,markersize=0.5, label='P1')
plt.grid()
plt.legend(loc=1)
#plt.yticks(np.arange(0, 2*np.pi+1, np.pi/2))
plt.savefig(os.path.join('%s_%s_doodson_all.png'%(times_doodsonplot[0].strftime('%Y%m%d'), times_doodsonplot[-1].strftime('%Y%m%d'))))

print('')
print('#####################################################################')
print('#### V0, FREQ TEST ##################################################')
print('#####################################################################')

v_0i_rad_F, t_const_freq_F = hatyan.get_foreman_v0_freq(const_list_freqv0uf, dood_date_v0freq)
v_0i_deg_F = np.rad2deg(v_0i_rad_F)

t_const_freq_H = hatyan.get_schureman_freqs(const_list_freqv0uf, dood_date_v0freq) #dood_date is not so relevant
v_0i_rad_H = hatyan.get_schureman_v0(const_list_freqv0uf, dood_date_v0freq)
v_0i_deg_H = np.rad2deg(v_0i_rad_H)

#export hatyan constituents and frequencies
t_const_freq_pd = hatyan.get_schureman_freqs(const_list_hatyan195_orig)
t_const_freq_pd.to_csv('hatyan_frequencies.csv')

freq_diff = t_const_freq_H['freq']-t_const_freq_F['freq']
v_0i_diff = v_0i_deg_H%(360)-v_0i_deg_F%(360)
freq_diff_toolarge_bool = freq_diff.abs()>treshold_freq
v_0i_diff_toolarge_bool = v_0i_diff[0].abs()>treshold_v0uf

freq_diff_hat55 = t_const_freq_H['freq']-hatyan55_freq_conv

freq_test = pd.DataFrame({'freq(hatyan)':t_const_freq_H['freq'], 'freq(foreman)':t_const_freq_F['freq'], 'freq(H-F)':freq_diff, 'diff_toolarge':freq_diff_toolarge_bool}, index=t_const_freq_H.index)
v0_test = pd.DataFrame({'v0(hatyan)':v_0i_deg_H.loc[:,0], 'v0(foreman)':v_0i_deg_F.loc[:,0], 'v0(H-F)':v_0i_diff.loc[:,0], 'diff_toolarge':freq_diff_toolarge_bool}, index=v_0i_deg_H.index)

print('#### FREQ is not date dependent in this implementation of foreman')
print(freq_test)
print('components with significant (%e) differences for freq:\n%s'%(treshold_freq, freq_diff[freq_diff_toolarge_bool]))
print('#### V0 for %s'%(dood_date_fu[0]))
print(v0_test)
print('components with significant (%e) differences for v0:\n%s'%(treshold_v0uf, v_0i_diff.loc[v_0i_diff_toolarge_bool,0]))
print('#### END V0, FREQ TEST ##################')

print('plotting...')
fig, (ax1,ax2) = plt.subplots(2,1, figsize=(15,8),sharex=True)
ax1.plot(t_const_freq_F,'-o', label='t_const_freq_F')
ax1.plot(t_const_freq_H['freq'],'-o', label='t_const_freq_H')
ax1.plot(freq_diff.values,'-o', label='t_const_freq: H-F')
ax1.plot(hatyan55_freq_conv,'x', markersize=8, label='freq hatyan55')
ax1.plot(freq_diff_hat55.values,'x', markersize=8, label='t_const_freq: H-hatyan55')
ax1.legend(loc=1)
ax1.set_xlim(-0.5,len(const_list_freqv0uf)+0.5)
ax1.set_ylim(-.1,.3)
ax1.grid()
ax1.set_xticks(range(len(const_list_freqv0uf)))
ax1.set_xticklabels(const_list_freqv0uf,rotation=90)#,minor=True) #set minor xticklabels, each constituent
ax1.set_ylabel('freq')

ax2.plot(v_0i_deg_F.loc[:,0]%(360),'-o', label='v_0i_deg_F')
ax2.plot(v_0i_deg_H.loc[:,0]%(360),'-o', label='v_0i_deg_H')
ax2.plot(v_0i_diff.loc[:,0],'-o', label='v_0i_deg: H-F')
if 1: #add hatyan55 v0_corr values
    hatyan55_v0corr = (hatyan55_v0u-np.rad2deg(u_i_rad_195.iloc[1:,:].values))%360
    v0_diff_hat55 = (np.rad2deg(v_0i_rad_H).iloc[1:,0].values-hatyan55_v0corr.loc[const_list_freqv0uf_hat55,year_v0]+180)%360-180
    ax2.plot(hatyan55_v0corr.loc[const_list_freqv0uf_hat55,year_v0],'x', markersize=8, label='v0corr hatyan55')
    ax2.plot(v0_diff_hat55,'x', markersize=8, label='v0: H-hatyan55')

ax2.legend(loc=1)
ax2.set_xlim(-0.5,len(const_list_freqv0uf)+0.5)
ax2.set_ylim(-180,360)
ax2.grid()
ax2.set_yticks(np.arange(-180,360+1,90))
ax2.set_xticks(range(len(const_list_freqv0uf)))
ax2.set_xticklabels(const_list_freqv0uf,rotation=90)#,minor=True) #set minor xticklabels, each constituent
ax2.set_ylabel('v0 [deg]')

fig.tight_layout()
fig.savefig('freq_v0.png')

#difference in frequency over time
print('calculating and plotting frequency difference over time')
times_freqdiff_forplot = times_freqdiff[:1].append(times_freqdiff)
const_list_forplot = ['']+const_list_freqv0uf #to make pcolormesh possible

t_const_freq, t_const_speed_all = hatyan.get_schureman_freqs(const_list_freqv0uf, dood_date=times_freqdiff, sort_onfreq=False, return_allraw=True)
t_const_freq_all = t_const_speed_all/(2*np.pi) #aantal rotaties per uur, freq
np.seterr(divide='ignore') #suppress divide by 0 warning
t_const_perds_all = 1/t_const_freq_all #period [hr]

fig,(ax1) = plt.subplots(1,1,figsize=(11,9),sharex=True)
data_speeddiff = t_const_speed_all-t_const_speed_all[:,:1]
#data_freqdiff = t_const_freq_all-t_const_freq_all[:,:1]
#data_perdsdiff = t_const_perds_all-t_const_perds_all[:,:1]
pc1 = ax1.pcolormesh(times_freqdiff_forplot, const_list_forplot, data_speeddiff)
#pc2 = ax2.pcolormesh(times_freqdiff_forplot, const_list_forplot, data_freqdiff)
#pc3 = ax3.pcolormesh(times_freqdiff_forplot, const_list_forplot, data_perdsdiff)
ax1.set_title('angvelo [rad/hr] difference w.r.t. first timestep')
#ax2.set_title('frequency difference over time')
#ax3.set_title('period difference over time')
fig.colorbar(pc1, ax=ax1)
#fig.colorbar(pc2, ax=ax2)
#fig.colorbar(pc3, ax=ax3)
fig.tight_layout()
fig.savefig('component_speed_difference.png')

print('...done')

#vuf-files per year
print('creating v0uf files per year...')
for year_sel in dood_date_fu.year.unique():
    dood_date_start = pd.DatetimeIndex([dt.datetime(year_sel,1,1)])
    dood_date_mid = pd.DatetimeIndex([dt.datetime(year_sel,7,2)])
    print(year_sel)
    
    v_0i_rad = hatyan.get_schureman_v0(const_list=const_list_freqv0uf, dood_date=dood_date_start) #v0 on beginning of year
    u_i_rad = hatyan.get_schureman_u(const_list=const_list_freqv0uf, dood_date=dood_date_mid)
    f_i_xfac0 = hatyan.get_schureman_f(const_list=const_list_freqv0uf, dood_date=dood_date_mid,xfac=False) #helemaal geen xfac toepassen, is voor eg S2 anders en dit is hatyan55 berekening
    f_i_xfac1 = hatyan.get_schureman_f(const_list=const_list_freqv0uf, dood_date=dood_date_mid,xfac=True)
    
    vuf_v0_deg = np.remainder(np.rad2deg(v_0i_rad),360)
    vuf_u_deg = np.remainder(np.rad2deg(u_i_rad)+180,360)-180
    vuf_v0u_deg = np.remainder(np.rad2deg(v_0i_rad+u_i_rad),360)
    
    const_no = [const_list_hatyan195_orig.index(x) for x in const_list_freqv0uf]
    dir_file = os.path.join(dir_output_vuf,'vuf_%i.txt'%(year_sel))
    
    with open(dir_file,'w') as f:
        #f.write('#### created by Python prototype of HATYAN 2.0 ####\n')
        f.write('  Nr. Hoeksnelh.  V0+u   f(xfac1)  Naam    V0     u      f(xfac0)   --- V0+u en f voor %s ---  \n'%(year_sel))
        f.write('      gr./uur     gr                       gr     gr              \n')
        f.write(' ---- ----------  -----  --------  ------  -----  -----  -------- \n')
        for iC, const in enumerate(const_list_freqv0uf):
            f.write("%5i %10.6f  %5.1f  %8.6f  %-6s  %5.1f  %5.1f  %8.6f \n" % (const_no[iC], t_const_freq_H.loc[const,'angfreq [deg/hr]'], vuf_v0u_deg.loc[const,0], f_i_xfac1.loc[const,0], const, vuf_v0_deg.loc[const,0], vuf_u_deg.loc[const,0], f_i_xfac0.loc[const,0]))

    outpd = pd.DataFrame({'constno':const_no, 'angfreq':t_const_freq_H['angfreq [deg/hr]'].round(6), 'v0+u':vuf_v0u_deg.loc[:,0].round(1), 'f(xfac1)':f_i_xfac1.loc[:,0].round(6), 'v0':vuf_v0_deg.loc[:,0].round(1), 'u':vuf_u_deg.loc[:,0].round(1), 'f(xfac0)':f_i_xfac0.loc[:,0].round(6)}, index=t_const_freq_H.index)
    outpd.to_csv(dir_file.replace('.txt','.csv'))
print('...done')

print('')
print('#####################################################################')
print('#### NODAL FACTORS TEST (F, U) ######################################')
print('#####################################################################')

# NODAL FACTORS (F, U) FOREMAN
f_i_FOR, u_i_rad_FOR = hatyan.get_foreman_nodalfactors(const_list_freqv0uf, dood_date_fu)
u_i_deg_FOR = np.rad2deg(u_i_rad_FOR)

# NODAL FACTORS (F, U) HATYAN
u_i_rad_HAT = hatyan.get_schureman_u(const_list_freqv0uf, dood_date_fu)
u_i_rad_HAT = (u_i_rad_HAT+np.pi)%(2*np.pi)-np.pi
u_i_deg_HAT = np.rad2deg(u_i_rad_HAT)
f_i_HAT_xfac0 = hatyan.get_schureman_f(const_list_freqv0uf, dood_date_fu, xfac=False)
f_i_HAT_xfac1 = hatyan.get_schureman_f(const_list_freqv0uf, dood_date_fu, xfac=True)

f_i_diff = f_i_HAT_xfac0-f_i_FOR
u_i_deg_diff = u_i_deg_HAT-u_i_deg_FOR
f_i_diff_toolarge_bool = f_i_diff[0].abs()>treshold_v0uf
u_i_deg_diff_toolarge_bool = u_i_deg_diff[0].abs()>treshold_v0uf

f_test = pd.DataFrame({'f(hatyan)':f_i_HAT_xfac0.loc[:,0], 'f(foreman)':f_i_FOR.loc[:,0], 'f(H-F)':f_i_diff.loc[:,0], 'diff_toolarge':f_i_diff_toolarge_bool}, index=f_i_HAT_xfac0.index)
u_test = pd.DataFrame({'u(hatyan)':u_i_deg_HAT.loc[:,0], 'u(foreman)':u_i_deg_FOR.loc[:,0], 'u(H-F)':u_i_deg_diff.loc[:,0], 'diff_toolarge':u_i_deg_diff_toolarge_bool}, index=u_i_deg_HAT.index)

print('#### NODAL FACTOR F for %s'%(dood_date_fu[0]))
print(f_test)
print('components with significant (%e) differences for f:\n%s'%(treshold_v0uf, f_i_diff.loc[f_i_diff_toolarge_bool,0]))
print('#### NODAL FACTOR U for %s'%(dood_date_fu[0]))
print(u_test)
print('components with significant (%e) differences for u:\n%s'%(treshold_v0uf, u_i_deg_diff.loc[u_i_deg_diff_toolarge_bool,0]))
print('#### END NODAL FACTOR TEST ##################')

print('plotting...')
fig, (ax1,ax2) = plt.subplots(2,1, figsize=(15,8))
for iC,const in enumerate(const_list_freqv0uf):
    print(const)
    if const in const_list_hatyan195_orig:
        iC_hat = const_list_hatyan195_orig.index(const)
    else:
        iC_hat = 999
    ax1.cla()
    ax1.set_title('f values [-] %03d %s'%(iC_hat, const))
    ax1.plot([dood_date_fu[0],dood_date_fu[-1]],[1,1],'k',linewidth=0.7)
    ax1.plot(dood_date_fu,f_i_FOR.loc[const,:], label='f_FOR')
    ax1.plot(dood_date_fu,f_i_HAT_xfac0.loc[const,:], label='f_HAT_xfac0')
    if const != 'A0':
        ax1.plot(dood_date_hatyan55,hatyan55_f.loc[const,:],'o',label='f hatyan55')
    if const in ['MU2','N2','NU2','M2','2MN2','S2','M4','MS4','M6','2MS6','M8','3MS8']: #xfac list
        ax1.plot(dood_date_fu,f_i_HAT_xfac1.loc[const,:], ':', label='f_HAT_xfac1')
        #if const in const_list_freqv0uf:
        ax1.plot(dt.datetime(2016,7,2),v0uffile_2016.loc[const,'f(xfac1)'],'o',label='f v0uffile 2016')
    ax1.set_xlim(dood_date_fu[0],dood_date_fu[-1])
    ax1.set_ylim(0.8,1.2)
    ax1.legend(loc=1)
    ax1.grid()
    
    ax2.cla()
    ax2.set_title('u values [deg] %03d %s'%(iC_hat, const))
    ax2.plot([dood_date_fu[0],dood_date_fu[-1]],[0,0],'k',linewidth=0.7)
    ax2.plot(dood_date_fu,u_i_deg_FOR.loc[const,:], label='u_FOR')
    ax2.plot(dood_date_fu,u_i_deg_HAT.loc[const,:], label='u_HAT')
    if const != 'A0':
        ax2.plot(dood_date_hatyan55,hatyan55_ucorr.loc[const,:],'o',label='v0+u hatyan55 - v0_HAT')
    ax2.set_xlim(dood_date_fu[0],dood_date_fu[-1])
    ax2.set_ylim(-12,12)
    ax2.legend(loc=1)
    ax2.grid()
    
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output_ufplots,'fu_%03d_%s.png'%(iC_hat,const)))
print('...done')

print('plotting v0u...')
fig, (ax1,ax2) = plt.subplots(2,1, figsize=(15,8))
for iC,const in enumerate(const_list_freqv0uf):
    if const in const_list_hatyan195_orig:
        iC_hat = const_list_hatyan195_orig.index(const)
    else:
        iC_hat = 999
    if const != 'A0':
        print(const)
        ax1.cla()
        ax1.set_title('v0+u values [deg] %03d %s'%(iC_hat, const))
        ax1.plot(dood_date_hatyan55_v0,v0u_deg_195.loc[const,:],'o',label='v0+u hatyan_python')
        ax1.plot(dood_date_hatyan55_v0,hatyan55_v0u.loc[const,:],'^',label='v0+u hatyan55')
        ax1.set_xlim(dood_date_hatyan55_v0[0],dood_date_hatyan55_v0[-1])
        ax1.set_ylim(-5,365)
        ax1.legend(loc=1)
        ax1.grid()

        ax2.cla()
        ax2.set_title('v0+u values diff [deg] %03d %s'%(iC_hat, const))
        ax2.plot([dood_date_hatyan55_v0[0],dood_date_hatyan55_v0[-1]],[0,0],'k',linewidth=0.7)
        ax2.plot(dood_date_hatyan55_v0,hatyan55py_v0u_diff.loc[const,:],'o',label='v0+u hatyan_python-hatyan55')
        ax2.set_xlim(dood_date_hatyan55_v0[0],dood_date_hatyan55_v0[-1])
        ax2.set_ylim(-12,12)
        ax2.legend(loc=1)
        ax2.grid()

        fig.tight_layout()
        fig.savefig(os.path.join(dir_output_v0uplots,'v0u_%03d_%s.png'%(iC_hat,const)))
print('...done')

hatyan.exit_RWS(timer_start)
