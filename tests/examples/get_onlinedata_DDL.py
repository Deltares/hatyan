# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 17:03:52 2021

@author: veenstra
"""

import ddlpy # TODO: we require ddlpy from main/master branch (>0.1.0) >> pip install git+https://github.com/openearth/ddlpy
import hatyan # available via `pip install hatyan` or at https://github.com/Deltares/hatyan
import datetime as dt
import matplotlib.pyplot as plt
plt.close("all")

# input parameters
tstart_dt = dt.datetime(2019,12,24) #period begins with Gecontroleerd and ends with Ongecontroleerd for HOEKVHLD
tstop_dt = dt.datetime(2020,1,5)
#tstart_dt = dt.datetime(2020,11,25,9,47,0) #quite recent period
#tstop_dt = dt.datetime(2021,1,30,9,50,0)
#tstart_dt = dt.datetime(1993,8,25,9,47,0) #VLISSGN got new Waardebepalingsmethode in this year
#tstop_dt = dt.datetime(1994,11,30,9,50,0)
#tstart_dt = dt.datetime(2009,1,1) #common RWS retrieval period
#tstop_dt = dt.datetime(2012,12,31,23,50)

dir_testdata = 'C:\\DATA\\hatyan_data_acceptancetests'


######### online waterlevel data retrieval for one station
if 1: #for RWS
    include_extremes = True

    locations = ddlpy.locations()
    locations["Code"] = locations.index
    # locations = locations_ddlpy.reset_index(drop=False).set_index('Locatie_MessageID')
    
    bool_hoedanigheid = locations['Hoedanigheid.Code'].isin(['NAP'])
    # bool_groepering = locations['Groepering.code'].isin(['NVT']) # TODO: we cannot subset on Groepering (NVT/GETETM2) yet: https://github.com/openearth/ddlpy/issues/21
    # bool_typering = locations['Typering.code'].isin(['GETETTPE']) # TODO: we cannot subset on Typering (NVT/GETETTPE) yet: https://github.com/openearth/ddlpy/issues/21
    
    # get wathte locations (ts and extremes) >> measured waterlevel
    bool_grootheid = locations['Grootheid.Code'].isin(['WATHTE'])
    locs_wathte = locations.loc[bool_grootheid & bool_hoedanigheid]
    
    # get WATHTBRKD locations (ts and extremes) >> computed astronomical waterlevel
    bool_grootheid = locations['Grootheid.Code'].isin(['WATHTBRKD'])
    locs_wathtbrkd = locations.loc[bool_grootheid & bool_hoedanigheid]

    # get types locations (ts and extremes) #TODO: replace with filter for Typering
    bool_grootheid = locations['Grootheid.Code'].isin(['NVT'])
    bool_eenheid = locations['Eenheid.Code'].isin(['DIMSLS'])
    locs_types = locations.loc[bool_grootheid & bool_eenheid]
    
    for current_station in ['HOEKVHLD']:
        locs_wathte_one = locs_wathte.loc[locs_wathte['Code'].isin([current_station])]
        locs_wathtbrkd_one = locs_wathtbrkd.loc[locs_wathtbrkd['Code'].isin([current_station])]
        locs_types_one = locs_types.loc[locs_types['Code'].isin([current_station])]
        
        # TODO: no support for multiple rows, create issue from example code in https://github.com/Deltares/hatyan/issues/187
        if len(locs_wathte_one) > 1:
            raise Exception("duplicate stations for wathte")
        else:
            locs_wathte_one = locs_wathte_one.iloc[0]
        if len(locs_wathtbrkd_one) > 1:
            raise Exception("duplicate stations for wathtbrkd")
        else:
            locs_wathtbrkd_one = locs_wathtbrkd_one.iloc[0]
        if len(locs_types_one) > 1:
            raise Exception("duplicate stations for types")
        else:
            locs_types_one = locs_types_one.iloc[0]
        
        meas_wathte = ddlpy.measurements(locs_wathte_one, start_date=tstart_dt, end_date=tstop_dt)
        meas_wathtbrkd = ddlpy.measurements(locs_wathtbrkd_one, start_date=tstart_dt, end_date=tstop_dt)
        meas_types = ddlpy.measurements(locs_types_one, start_date=tstart_dt, end_date=tstop_dt)
        
        # TODO: rename lowercase code to uppercase Code
        meas_wathte.columns = [x.replace(".code",".Code") for x in meas_wathte.columns]
        meas_wathtbrkd.columns = [x.replace(".code",".Code") for x in meas_wathtbrkd.columns]
        meas_types.columns = [x.replace(".code",".Code") for x in meas_types.columns]
        
        # timeseries
        meas_wathte_ts = meas_wathte.loc[meas_wathte['Groepering.Code'].isin(['NVT'])]
        ts_measwl = hatyan.ddlpy_to_hatyan(meas_wathte_ts)
        meas_wathtbrkd_ts = meas_wathtbrkd.loc[meas_wathtbrkd['Groepering.Code'].isin(['NVT'])]
        ts_astro = hatyan.ddlpy_to_hatyan(meas_wathtbrkd_ts)
        
        if include_extremes:
            #extremes
            meas_wathte_ext = meas_wathte.loc[meas_wathte['Groepering.Code'].isin(['GETETM2'])]
            ts_measwlHWLW = hatyan.ddlpy_to_hatyan(meas_wathte_ext)
            meas_wathtbrkd_ext = meas_wathtbrkd.loc[meas_wathtbrkd['Groepering.Code'].isin(['GETETBRKD2'])]
            ts_astroHWLW = hatyan.ddlpy_to_hatyan(meas_wathtbrkd_ext)
            # extreme types
            meas_wathte_exttype = meas_types.loc[meas_types['Groepering.Code'].isin(['GETETM2'])]
            ts_measwlHWLWtype = hatyan.ddlpy_to_hatyan(meas_wathte_exttype)
            meas_wathtbrkd_exttype = meas_types.loc[meas_types['Groepering.Code'].isin(['GETETBRKD2'])]
            ts_astroHWLWtype = hatyan.ddlpy_to_hatyan(meas_wathtbrkd_exttype)
            
            ts_measwlHWLW = hatyan.convert_HWLWstr2num(ts_measwlHWLW,ts_measwlHWLWtype)
            ts_astrolHWLW = hatyan.convert_HWLWstr2num(ts_astroHWLW,ts_astroHWLWtype)

    stat_name = locs_wathte_one["Naam"]
    stat_code = locs_wathte_one["Code"]
    if include_extremes:
        fig,(ax1,ax2) = hatyan.plot_timeseries(ts=ts_astro,ts_validation=ts_measwl,ts_ext=ts_astroHWLW,ts_ext_validation=ts_measwlHWLW)
        ax1.set_title(f'waterlevels for {stat_code} ({stat_name})')
        ax2.set_ylim(-0.5,0.5)
    else:
        fig,(ax1,ax2) = hatyan.plot_timeseries(ts=ts_astro,ts_validation=ts_measwl)
        ax1.set_title(f'waterlevels for {stat_code} ({stat_name})')
        ax2.set_ylim(-0.5,0.5)


######### simple waterlevel data retrieval for all waterlevel stations or all stations
if 1: #for CMEMS

    locations = ddlpy.locations()
    locations["Code"] = locations.index #TODO: maybe check in index directly?
    
    bool_hoedanigheid = locations['Hoedanigheid.Code'].isin(['NAP'])
    # bool_groepering = locations['Groepering.code'].isin(['NVT']) # TODO: we cannot subset on Groepering (NVT/GETETM2) yet: https://github.com/openearth/ddlpy/issues/21

    # get wathte locations (ts and extremes) >> measured waterlevel
    bool_grootheid = locations['Grootheid.Code'].isin(['WATHTE'])
    locs_wathte = locations.loc[bool_grootheid & bool_hoedanigheid]
    
    for current_station in ['HOEKVHLD']:
        locs_wathte_one = locs_wathte.loc[locs_wathte['Code'].isin([current_station])]
        
        # TODO: no support for multiple rows, create issue from example code in https://github.com/Deltares/hatyan/issues/187
        if len(locs_wathte_one) > 1:
            raise Exception("duplicate stations for wathte")
        else:
            locs_wathte_one = locs_wathte_one.iloc[0]
        
        meas_wathte = ddlpy.measurements(locs_wathte_one, start_date=tstart_dt, end_date=tstop_dt)
        
        # TODO: rename lowercase code to uppercase Code
        meas_wathte.columns = [x.replace(".code",".Code") for x in meas_wathte.columns]
        
        # timeseries
        meas_wathte_ts = meas_wathte.loc[meas_wathte['Groepering.Code'].isin(['NVT'])]
        ts_measwl = hatyan.ddlpy_to_hatyan(meas_wathte_ts)
        
        stat_name = locs_wathte_one["Naam"]
        stat_code = locs_wathte_one["Code"]
        fig, (ax1,ax2) = plt.subplots(2,1, figsize=(8,6), sharex=True)
        ts_measwl["values"].plot(ax=ax1)
        ts_measwl["QC"].plot(ax=ax2)
        ax1.set_title(f'waterlevels for {stat_code} ({stat_name})')
        ax2.set_title(f'QC for {stat_code} ({stat_name})')
        fig.tight_layout()
