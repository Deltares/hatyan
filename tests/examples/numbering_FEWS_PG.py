# -*- coding: utf-8 -*-
"""

"""

import os
import datetime as dt
import pandas as pd
import matplotlib.pyplot as plt
plt.close('all')
from netCDF4 import Dataset, num2date
import hatyan

#dir_testdata = 'P:\\1209447-kpp-hydraulicaprogrammatuur\\hatyan\\hatyan_data_acceptancetests'
dir_testdata = 'C:\\DATA\\hatyan_data_acceptancetests'

analyse_ts_bool = False

file_meas = os.path.join(dir_testdata,'other','FEWS_202010221200_testdata_S_2.nc')
data_nc_meas = Dataset(file_meas)
data_meas = data_nc_meas.variables['water_level'][:].flatten()
data_times = num2date(data_nc_meas.variables['time'],units=data_nc_meas.variables['time'].units, only_use_cftime_datetimes=False, only_use_python_datetimes=True)
ts_meas = pd.DataFrame({'values':data_meas},index=data_times)

file_pred = os.path.join(dir_testdata,'other','FEWS_202010221200_testdata_S_4.nc')
data_nc_pred = Dataset(file_pred)
data_pred = data_nc_pred.variables['water_level'][:].flatten()
data_times = num2date(data_nc_pred.variables['time'],units=data_nc_pred.variables['time'].units, only_use_cftime_datetimes=False, only_use_python_datetimes=True)
ts_prediction = pd.DataFrame({'values':data_pred},index=data_times)

file_comp = os.path.join(dir_testdata,'predictie2019','HOEKVHLD_ana.txt')

#station_name = data_pred.var_stations.loc[0,'node_id']
station_name = 'HOEKVHLD'

file_ncout = '%s_getijnummers_new.nc'%(station_name)
file_ncout_nosidx = '%s_getijnummers_nosidx.nc'%(station_name)

if analyse_ts_bool:
    COMP_merged = hatyan.analysis(ts=ts_meas, const_list='year')
    ts_prediction_fromcomp = hatyan.prediction(comp=COMP_merged, times_pred_all=ts_prediction.index, xfac=True)#, fu_alltimes=False)
else:
    COMP_merged = hatyan.read_components(filename=file_comp)
    ts_prediction_fromcomp_2019 = hatyan.prediction(comp=COMP_merged, times_ext=[dt.datetime(2019,1,1),dt.datetime(2019,12,31,23,50)], timestep_min=10, xfac=True, fu_alltimes=False)
    ts_prediction_fromcomp_2020 = hatyan.prediction(comp=COMP_merged, times_ext=[dt.datetime(2020,1,1),dt.datetime(2020,12,31,23,50)], timestep_min=10, xfac=True, fu_alltimes=False)
    ts_prediction_fromcomp = ts_prediction_fromcomp_2019.append(ts_prediction_fromcomp_2020)
    ts_prediction_fromcomp.index = ts_prediction_fromcomp.index-dt.timedelta(hours=1) #convert MET prediction to GMT

print(hatyan.check_ts(ts_prediction))
ts_ext_prediction = hatyan.calc_HWLW(ts=ts_prediction)#, calc_HWLWlocal=True)
ts_ext_prediction_nos = hatyan.calc_HWLWnumbering(ts_ext=ts_ext_prediction, station=station_name)

fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction, ts_validation=ts_prediction_fromcomp, ts_ext=ts_ext_prediction)
ax2.set_ylim(-0.2,0.2)
fig.savefig(file_ncout.replace('.nc','.png'))

fig, (ax1,ax2) = plt.subplots(2,1,figsize=(14,5))
ax1.plot(ts_ext_prediction_nos.index,ts_ext_prediction_nos['HWLWno'],'o')
ax1.set_ylim(ts_ext_prediction_nos['HWLWno'].min()-100,ts_ext_prediction_nos['HWLWno'].max()+100)
ax1.set_ylabel('tide_number_dummies')
ax2.plot(ts_ext_prediction_nos.index,ts_ext_prediction_nos['HWLWno'].diff(),'o')
ax2.set_ylim(-1,2)
fig.savefig(file_ncout.replace('.nc','_nrs.png'))

"""
hatyan.write_tsnetcdf(ts=ts_prediction, station=station_name, vertref='NAP', filename=file_ncout, ts_ext=ts_ext_prediction_nos, tzone_hr=0, nosidx=False)
#from dfm_tools.get_nc_helpers import get_ncvardimlist
#vars_pd, dims_pd = get_ncvardimlist(file_nc=file_ncout)
#data_nc_checkLWval = get_ncmodeldata(file_nc=file_ncout,varname='time_LW',timestep='all',station=0)
data_ncout = Dataset(file_ncout)
data_ncout.variables['waterlevel_astro_LW_numbers']
data_ncout.variables['waterlevel_astro_HW_numbers']
"""
hatyan.write_tsnetcdf(ts=ts_prediction, station=station_name, vertref='NAP', filename=file_ncout_nosidx, ts_ext=ts_ext_prediction_nos, tzone_hr=0, nosidx=True)
hatyan.write_tsnetcdf(ts=ts_prediction, station='HOEKVHLD_copy', vertref='NAP', filename=file_ncout_nosidx, ts_ext=ts_ext_prediction_nos, tzone_hr=0, nosidx=True, mode='a')
data_ncout = Dataset(file_ncout_nosidx)
data_ncout.variables['waterlevel_astro_HW']
data_ncout.variables['HWLWno']
times_all_out = num2date(data_ncout.variables['time'],units=data_ncout.variables['time'].units, only_use_cftime_datetimes=False, only_use_python_datetimes=True)
