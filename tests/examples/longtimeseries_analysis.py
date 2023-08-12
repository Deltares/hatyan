# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 14:33:50 2021

@author: veenstra
"""


import os
import datetime as dt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')
from netCDF4 import Dataset, num2date #TODO: move to xarray
import hatyan

dir_data = r'p:\archivedprojects\11205259-004-dcsm-fm\waterlevel_data\all\ncFiles'
selected_stations = ['A2','A12','AUKFPFM','AWGPFM','D15','EURPFM','F3PFM','F16','K13APFM','K14PFM','L9PFM'] #all platforms I can find
selected_stations = ['A2','AUKFPFM','EURPFM','F3PFM','K13APFM'] #more than 17 years of data
selected_stations = ['EURPFM']

do_analysis = True #set to false to check only timeseries length and tstart/tstop
#dir_output = os.path.join(dir_base,case2,'cotidal_FESGTSM_divfM2')
#if not os.path.exists(dir_output):
#    os.mkdir(dir_output)

for current_station in selected_stations:
    print('')
    print(current_station)
    file_nc = os.path.join(dir_data,f'{current_station}.nc')
    
    data_nc = Dataset(file_nc)
    
    data_nc_timevar = data_nc.variables['time']
    
    # get time ids for start and stop time ([0,-1] would also work, but ok)
    time_length = data_nc_timevar.shape[0]
    retrieve_ids_time = [0,-1] #start and stoptime
    retrieve_ids_time = list(np.unique(np.array(range(time_length))[retrieve_ids_time]))

    #get start/stop time and convert to pandas Series of datetime
    data_nc_times = data_nc_timevar[retrieve_ids_time]
    data_nc_datetimes = num2date(data_nc_times, units=data_nc_timevar.units, only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    data_nc_datetimes_pd = pd.Series(data_nc_datetimes,index=retrieve_ids_time).dt.round(freq='S')
    time_length_yr = (data_nc_datetimes_pd.iloc[-1]-data_nc_datetimes_pd.iloc[0]).total_seconds()/3600/24/365.25
    print('ts length is %.2f years'%(time_length_yr))
    print(data_nc_datetimes_pd)

    var_wl = data_nc.variables['waterlevel']

    if do_analysis:
        #constituent list
        #const_list = hatyan.get_const_list_hatyan('year') #this is not a good year when you analyse 35 year of data with 60 or 10 minute interval, but is useful after subsampling
        const_list = hatyan.get_const_list_hatyan('day') # choose from: dict_keys(['all', 'year', 'halfyear', 'month', 'month_deepwater', 'springneap', 'day', 'tidalcycle'])
        #const_list = ['M2','S2']
        
        print('reading data from netcdf')
        meas_times_raw = data_nc_timevar[:]
        meas_times = num2date(meas_times_raw, units=data_nc_timevar.units, only_use_cftime_datetimes=False, only_use_python_datetimes=True)
        meas_vals = var_wl[:]
        print('finished reading netcdf')
        ts_measurements = pd.DataFrame({'values': meas_vals},index=meas_times)
        
        #analysis
        comp_frommeasurements = hatyan.analysis(ts=ts_measurements, const_list=const_list, nodalfactors=True, xfac=False, fu_alltimes=True)
        print(comp_frommeasurements)
        
        fig,(ax1,ax2) = hatyan.plot_components(comp_frommeasurements)
        #fig.savefig('components_%s_4Y.png'%(current_station))
        
        #prediction and validation
        times_pred = slice(dt.datetime(2008,1,1),dt.datetime(2010,1,1), 10)
        ts_prediction = hatyan.prediction(comp=comp_frommeasurements, nodalfactors=True, xfac=False, fu_alltimes=True, times=times_pred)
        fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction, ts_validation=ts_measurements)
        ax1.set_ylim(-2.5,2.5)
        ax2.set_ylim(-0.5,0.5)
        #fig.savefig('prediction_%im_%s_measurements'%(times_step_pred, current_station))

