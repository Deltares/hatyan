# -*- coding: utf-8 -*-
"""
astrog_test.py

bereken en vergelijk met FORTRAN:
    - culminaties
    - opkomst en ondergang maan en zon
    - maanfasen (schijngestalten)
    - maananomalieën
    - astronomische seizoenen

schrijf resultaat (python) weg als csv (in MET) met pd.to_csv()

"""

import os, pytz
import numpy as np
import datetime as dt
import pandas as pd
import matplotlib.pyplot as plt
plt.close('all')
import hatyan

#dir_testdata = 'P:\\1209447-kpp-hydraulicaprogrammatuur\\hatyan\\hatyan_data_acceptancetests'
dir_testdata = 'C:\\DATA\\hatyan_data_acceptancetests'

# script settings
compare2fortran = True #requires validation data
timeStart = dt.datetime(2000,1,1)
timeEnd   = dt.datetime(2012,1,1)
dT_fortran = True #True is best comparison to fortran, False is more precise 
tz_GMT = 'UTC' # UTC/GMT timezone
try: #only works with pandas 1.2.0 or higher, but 1.1.5 is highest available pandas version in Python 3.6.12 (which is the highest available Python in RHEL)
    tz_MET = 'UTC+01:00' #for timezone conversion to UTC+01:00 (Dutch wintertime)
    pytz.timezone(tz_MET)
except:
    tz_MET = dt.timezone(dt.timedelta(hours=1)) #for timezone conversion to UTC+01:00 (Dutch wintertime)
tz_EurAms = 'Europe/Amsterdam' # for conversion to local timezone, including daylight saving time (DST)

pdtocsv_kwargs = dict(index=False, sep=',', date_format='%Y-%m-%d %H:%M:%S %Z', float_format='%9.5f', na_rep='--:--          ')

#%% calculate astrog arrays
# lunar culmination times, parallax, declination
culminations_python = hatyan.astrog_culminations(tFirst=timeStart, tLast=timeEnd, dT_fortran=dT_fortran, tzone=tz_GMT)
culminations_python.to_csv('moon_culminations.csv',**pdtocsv_kwargs)

# lunar phases
phases_python = hatyan.astrog_phases(tFirst=timeStart, tLast=timeEnd, dT_fortran=dT_fortran, tzone=tz_MET)
phases_python.to_csv('moon_phases.csv',**pdtocsv_kwargs)

# moonrise and -set
moonriseset_python = hatyan.astrog_moonriseset(tFirst=timeStart, tLast=timeEnd, dT_fortran=dT_fortran, tzone=tz_MET)
moonriseset_python.to_csv('moon_riseset.csv',**pdtocsv_kwargs)
moonriseset_python_perday = hatyan.convert2perday(moonriseset_python)
moonriseset_python_perday.to_csv('moon_riseset_perday.csv',**pdtocsv_kwargs)

# sunrise and -set
sunriseset_python = hatyan.astrog_sunriseset(tFirst=timeStart, tLast=timeEnd, dT_fortran=dT_fortran, tzone=tz_MET)
sunriseset_python.to_csv('sun_riseset.csv',**pdtocsv_kwargs)
sunriseset_python_perday = hatyan.convert2perday(sunriseset_python)
sunriseset_python_perday.to_csv('sun_riseset_perday.csv',**pdtocsv_kwargs)

# lunar anomalies
anomalies_python = hatyan.astrog_anomalies(tFirst=timeStart, tLast=timeEnd, dT_fortran=dT_fortran, tzone=tz_MET)
anomalies_python.to_csv('anomalies.csv',**pdtocsv_kwargs)

# astronomical seasons
seasons_python = hatyan.astrog_seasons(tFirst=timeStart, tLast=timeEnd, dT_fortran=dT_fortran, tzone=tz_MET)
seasons_python.to_csv('seasons.csv',**pdtocsv_kwargs)

if compare2fortran:
    #%% load fortran results
    pkl_culm = os.path.join(dir_testdata,'other','astrog20_2000_2011.pkl')
    pkl_phas = os.path.join(dir_testdata,'other','astrog30_phases_2000_2011.pkl')
    txt_phas = os.path.join(dir_testdata,'other','maanfasen.txt')
    pkl_moon = os.path.join(dir_testdata,'other','astrog30_moon_2000_2011.pkl')
    pkl_sun  = os.path.join(dir_testdata,'other','astrog30_sun_2000_2011.pkl')
    pkl_anom = os.path.join(dir_testdata,'other','astrog30_anomalies_2000_2011.pkl')
    pkl_seas = os.path.join(dir_testdata,'other','astrog30_seasons_2000_2011.pkl')
    
    culminations_fortran = pd.read_pickle(pkl_culm)
    culminations_fortran = culminations_fortran[np.logical_and(culminations_fortran['datetime']>=timeStart,culminations_fortran['datetime']<=timeEnd)].reset_index(drop=True)
    
    phases_fortran = pd.read_pickle(pkl_phas)
    phases_fortran = phases_fortran[np.logical_and(phases_fortran['datetime']>=timeStart,phases_fortran['datetime']<=timeEnd)].reset_index(drop=True)
    
    phases_long_fortran = pd.read_csv(txt_phas, sep=';', names=['date','time','type_str'], skiprows=1) # long time series 2021-2035 (Koos Doekes)
    phases_long_fortran['datetime']=pd.to_datetime(phases_long_fortran['date'].astype(str)+phases_long_fortran['time'].astype(str).str.zfill(4))
    phases_long_fortran['type'] = phases_long_fortran['type_str'].replace('EK',1).replace('VM',2).replace('LK',3).replace('NM',4)
    phases_long_python = hatyan.astrog_phases(phases_long_fortran['datetime'].iloc[0]-dt.timedelta(days=5), phases_long_fortran['datetime'].iloc[-1]+dt.timedelta(days=5), dT_fortran=dT_fortran)
    phases_long_python['datetime'] = phases_long_python['datetime'].dt.tz_convert(tz_EurAms) #convert to local timezone
    
    moonriseset_fortran = pd.read_pickle(pkl_moon)
    moonriseset_fortran = moonriseset_fortran[np.logical_and(moonriseset_fortran['datetime']>=timeStart,moonriseset_fortran['datetime']<=timeEnd)]
    
    sunriseset_fortran = pd.read_pickle(pkl_sun)
    sunriseset_fortran = sunriseset_fortran[np.logical_and(sunriseset_fortran['datetime']>=timeStart,sunriseset_fortran['datetime']<=timeEnd)]
    pyth_index         = sunriseset_python['datetime'].dt.date.isin(sunriseset_fortran['datetime'].dt.date.unique())
    sunriseset_python_somedays = sunriseset_python.loc[pyth_index].reset_index(drop=True)
    
    anomalies_fortran = pd.read_pickle(pkl_anom)
    anomalies_fortran = anomalies_fortran[np.logical_and(anomalies_fortran['datetime']>=timeStart,anomalies_fortran['datetime']<=timeEnd)].reset_index(drop=True)
    
    seasons_fortran = pd.read_pickle(pkl_seas)
    seasons_fortran = seasons_fortran[np.logical_and(seasons_fortran['datetime']>=timeStart,seasons_fortran['datetime']<=timeEnd)].reset_index(drop=True)
    
    #%% plot results (differences)
    fig, (ax1,ax2,ax3) = hatyan.plot_astrog_diff(culminations_python, culminations_fortran, typeLab=['lower','upper'], timeBand=[-.18,.18])
    fig.savefig('culmination_differences.png')
    
    fig, (ax1,ax2,ax3) = hatyan.plot_astrog_diff(culminations_python, culminations_fortran, typeCol="parallax", typeUnit='degrees', timeBand=[-.18,.18], typeBand=[-.000005,.000005])
    fig.savefig('culmination_differences_parallax.png')
    
    fig, (ax1,ax2,ax3) = hatyan.plot_astrog_diff(culminations_python, culminations_fortran, typeCol="declination", typeUnit='degrees', timeBand=[-.18,.18], typeBand=[-.0005,.0005])
    fig.savefig('culmination_differences_declination.png')

    fig, (ax1,ax2,ax3) = hatyan.plot_astrog_diff(pd_python=phases_python, pd_fortran=phases_fortran, typeLab=['FQ','FM','LQ','NM'], timeBand=[-30,30])
    fig.savefig('phase_differences.png')

    fig, (ax1,ax2,ax3) = hatyan.plot_astrog_diff(pd_python=phases_long_python, pd_fortran=phases_long_fortran[['datetime','type']], typeLab=['FQ','FM','LQ','NM'], timeBand=[-30,30])
    fig.savefig('phase_differences_longperiod.png')

    fig, (ax1,ax2,ax3) = hatyan.plot_astrog_diff(pd_python=moonriseset_python, pd_fortran=moonriseset_fortran, typeLab=['rise','set'], timeBand=[-30,30])
    fig.savefig('moonRiseSet_differences.png')
    
    fig, (ax1,ax2,ax3) = hatyan.plot_astrog_diff(sunriseset_python_somedays, sunriseset_fortran, typeLab=['rise','set'], timeBand=[-30,30])
    fig.savefig('sunRiseSet_differences.png')
    
    fig, (ax1,ax2,ax3) = hatyan.plot_astrog_diff(pd_python=anomalies_python, pd_fortran=anomalies_fortran, typeLab=['perigeum','apogeum'], timeBand=[-1800,1800])
    fig.savefig('anomaly_differences.png')
    
    fig, (ax1,ax2,ax3) = hatyan.plot_astrog_diff(pd_python=seasons_python, pd_fortran=seasons_fortran, typeLab=['spring','summer','autumn','winter'], timeBand=[-30,30])
    fig.savefig('season_differences.png')

