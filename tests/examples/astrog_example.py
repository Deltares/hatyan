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

import os
import datetime as dt
import pandas as pd
import matplotlib.pyplot as plt
plt.close('all')
import hatyan

#dir_testdata = 'P:\\1209447-kpp-hydraulicaprogrammatuur\\hatyan\\hatyan_data_acceptancetests'
dir_testdata = 'C:\\DATA\\hatyan_data_acceptancetests'

compare2fortran = True #requires validation data

start_date_utc = pd.Timestamp(2000, 1, 1, tz="UTC")
end_date_utc = pd.Timestamp(2012, 1, 1, tz="UTC")
start_date_naive = start_date_utc.tz_convert(None)
end_date_naive = end_date_utc.tz_convert(None)
start_date_met = start_date_utc.tz_convert('UTC+01:00')
end_date_met = end_date_utc.tz_convert('UTC+01:00')

tz_EurAms = 'Europe/Amsterdam' # for conversion to local timezone, including daylight saving time (DST)

dT_fortran = True #True is best comparison to fortran, False is more precise 

pdtocsv_kwargs = dict(index=False, sep=',', date_format='%Y-%m-%d %H:%M:%S %Z', float_format='%9.5f', na_rep='--:--')

# calculate astrog lunar culmination times, parallax, declination
culminations_python = hatyan.astrog_culminations(tFirst=start_date_utc, tLast=end_date_utc, dT_fortran=dT_fortran)
culminations_python.to_csv('moon_culminations.csv',**pdtocsv_kwargs)

# calculate astrog lunar phases
phases_python = hatyan.astrog_phases(tFirst=start_date_met, tLast=end_date_met, dT_fortran=dT_fortran)
phases_python.to_csv('moon_phases.csv',**pdtocsv_kwargs)

# calculate astrog moonrise and -set
moonriseset_python = hatyan.astrog_moonriseset(tFirst=start_date_met, tLast=end_date_met, dT_fortran=dT_fortran)
moonriseset_python.to_csv('moon_riseset.csv',**pdtocsv_kwargs)
moonriseset_python_perday = hatyan.convert2perday(moonriseset_python)
moonriseset_python_perday.to_csv('moon_riseset_perday.csv',**pdtocsv_kwargs)

# calculate astrog sunrise and -set
sunriseset_python = hatyan.astrog_sunriseset(tFirst=start_date_met, tLast=end_date_met, dT_fortran=dT_fortran)
sunriseset_python.to_csv('sun_riseset.csv',**pdtocsv_kwargs)
sunriseset_python_perday = hatyan.convert2perday(sunriseset_python)
sunriseset_python_perday.to_csv('sun_riseset_perday.csv',**pdtocsv_kwargs)

# calculate astrog lunar anomalies
anomalies_python = hatyan.astrog_anomalies(tFirst=start_date_met, tLast=end_date_met, dT_fortran=dT_fortran)
anomalies_python.to_csv('anomalies.csv',**pdtocsv_kwargs)

# calculate astrog astronomical seasons
seasons_python = hatyan.astrog_seasons(tFirst=start_date_met, tLast=end_date_met, dT_fortran=dT_fortran)
seasons_python.to_csv('seasons.csv',**pdtocsv_kwargs)

if compare2fortran:
    # load fortran results
    pkl_culm = os.path.join(dir_testdata,'other','astrog20_2000_2011.pkl')
    pkl_phas = os.path.join(dir_testdata,'other','astrog30_phases_2000_2011.pkl')
    txt_phas = os.path.join(dir_testdata,'other','maanfasen.txt')
    pkl_moon = os.path.join(dir_testdata,'other','astrog30_moon_2000_2011.pkl')
    pkl_sun  = os.path.join(dir_testdata,'other','astrog30_sun_2000_2011.pkl')
    pkl_anom = os.path.join(dir_testdata,'other','astrog30_anomalies_2000_2011.pkl')
    pkl_seas = os.path.join(dir_testdata,'other','astrog30_seasons_2000_2011.pkl')
    
    culminations_fortran = pd.read_pickle(pkl_culm)
    culminations_fortran['datetime'] = pd.to_datetime(culminations_fortran['datetime'])
    culminations_fortran = culminations_fortran.set_index('datetime')
    culminations_fortran = culminations_fortran.loc[start_date_naive:end_date_naive]
    
    phases_fortran = pd.read_pickle(pkl_phas).set_index('datetime')
    phases_fortran = phases_fortran.loc[start_date_naive:end_date_naive]
    
    phases_long_fortran = pd.read_csv(txt_phas, sep=';', names=['date','time','type_str'], skiprows=1) # long time series 2021-2035 (Koos Doekes)
    phases_long_fortran['datetime'] = pd.to_datetime(phases_long_fortran['date'].astype(str)+phases_long_fortran['time'].astype(str).str.zfill(4))
    phases_long_fortran = phases_long_fortran.set_index('datetime')
    phases_long_fortran['type'] = phases_long_fortran['type_str'].str.replace('EK','1').str.replace('VM','2').str.replace('LK','3').str.replace('NM','4').astype(int)
    phases_long_python = hatyan.astrog_phases(phases_long_fortran.index[0]-dt.timedelta(days=5), phases_long_fortran.index[-1]+dt.timedelta(days=5), dT_fortran=dT_fortran)
    phases_long_python = phases_long_python.tz_convert(tz_EurAms) #convert to local timezone
    
    moonriseset_fortran = pd.read_pickle(pkl_moon).set_index('datetime')
    moonriseset_fortran = moonriseset_fortran.loc[start_date_naive:end_date_naive]
    
    sunriseset_fortran = pd.read_pickle(pkl_sun).set_index('datetime')
    sunriseset_fortran = sunriseset_fortran.loc[start_date_naive:end_date_naive]
    selected_dates = sunriseset_fortran.index.date
    sunriseset_python_somedays = sunriseset_python.loc[pd.Index(sunriseset_python.index.date).isin(selected_dates)]
    
    anomalies_fortran = pd.read_pickle(pkl_anom).set_index('datetime')
    anomalies_fortran = anomalies_fortran.loc[start_date_naive:end_date_naive]
    
    seasons_fortran = pd.read_pickle(pkl_seas).set_index('datetime')
    seasons_fortran = seasons_fortran.loc[start_date_naive:end_date_naive]
    
    # plot results (differences)
    fig, (ax1,ax2,ax3) = hatyan.plot_astrog_diff(culminations_python, culminations_fortran, typeLab=['lower','upper'], timeBand=[-.18,.18])
    fig.savefig('culmination_differences.png')
    
    fig, (ax1,ax2,ax3) = hatyan.plot_astrog_diff(culminations_python, culminations_fortran, typeCol="parallax", typeUnit='degrees', timeBand=[-.18,.18], typeBand=[-.000005,.000005])
    fig.savefig('culmination_differences_parallax.png')
    
    fig, (ax1,ax2,ax3) = hatyan.plot_astrog_diff(culminations_python, culminations_fortran, typeCol="declination", typeUnit='degrees', timeBand=[-.18,.18], typeBand=[-.0005,.0005])
    fig.savefig('culmination_differences_declination.png')

    fig, (ax1,ax2,ax3) = hatyan.plot_astrog_diff(pd_python=phases_python, pd_fortran=phases_fortran, typeLab=['FQ','FM','LQ','NM'], timeBand=[-30,30])
    fig.savefig('phase_differences.png')

    fig, (ax1,ax2,ax3) = hatyan.plot_astrog_diff(pd_python=phases_long_python, pd_fortran=phases_long_fortran[['type']], typeLab=['FQ','FM','LQ','NM'], timeBand=[-30,30])
    fig.savefig('phase_differences_longperiod.png')

    fig, (ax1,ax2,ax3) = hatyan.plot_astrog_diff(pd_python=moonriseset_python, pd_fortran=moonriseset_fortran, typeLab=['rise','set'], timeBand=[-30,30])
    fig.savefig('moonRiseSet_differences.png')
    
    fig, (ax1,ax2,ax3) = hatyan.plot_astrog_diff(sunriseset_python_somedays, sunriseset_fortran, typeLab=['rise','set'], timeBand=[-30,30])
    fig.savefig('sunRiseSet_differences.png')
    
    fig, (ax1,ax2,ax3) = hatyan.plot_astrog_diff(pd_python=anomalies_python, pd_fortran=anomalies_fortran, typeLab=['perigeum','apogeum'], timeBand=[-1800,1800])
    fig.savefig('anomaly_differences.png')
    
    fig, (ax1,ax2,ax3) = hatyan.plot_astrog_diff(pd_python=seasons_python, pd_fortran=seasons_fortran, typeLab=['spring','summer','autumn','winter'], timeBand=[-30,30])
    fig.savefig('season_differences.png')
