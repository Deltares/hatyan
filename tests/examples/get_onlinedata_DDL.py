# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 17:03:52 2021

@author: veenstra
"""

import datetime as dt
import matplotlib.pyplot as plt
plt.close("all")
import ddlpy # TODO: we require ddlpy from main/master branch (>0.1.0) >> pip install git+https://github.com/openearth/ddlpy

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


def ddlpy_to_hatyan(ddlpy_meas):
    cols_mustbe_unique = ['Grootheid.Code','Groepering.Code','Typering.Code','Hoedanigheid.Code']
    for col in cols_mustbe_unique:
        if len(ddlpy_meas[col].drop_duplicates()) != 1:
            raise Exception(f"ddlpy_meas['{col}'] is not unique")
    
    import pandas as pd
    #TODO: below copied from hatyan.getonlinedata.py (more TODO in that script). also in dfm_tools.observations.py
    key_numericvalues = 'Meetwaarde.Waarde_Numeriek'
    isnumeric = True
    if not key_numericvalues in ddlpy_meas.columns: #alfanumeric values for 'Typering.Code':'GETETTPE' #DDL IMPROVEMENT: also include numeric values for getijtype. Also, it is quite complex to get this data in the first place, would be convenient if it would be a column when retrieving 'Groepering.Code':'GETETM2' or 'GETETBRKD2'
        key_numericvalues = 'Meetwaarde.Waarde_Alfanumeriek'
        isnumeric = False
    ts_pd = pd.DataFrame({'values':ddlpy_meas[key_numericvalues].values,
                         'QC':pd.to_numeric(ddlpy_meas['WaarnemingMetadata.KwaliteitswaardecodeLijst'].str[0],downcast='integer').values, # DDL IMPROVEMENT: should be possible with .astype(int), but pd.to_numeric() is necessary for HARVT10 (eg 2019-09-01 to 2019-11-01) since QC contains None values that cannot be ints (in that case array of floats with some nans is returned) >> replace None with int code
                         'Status':ddlpy_meas['WaarnemingMetadata.StatuswaardeLijst'].str[0].values,
                         #'Bemonsteringshoogte':measurements_wathte['WaarnemingMetadata.BemonsteringshoogteLijst'].str[0].astype(int).values, 
                         #'Referentievlak':measurements_wathte['WaarnemingMetadata.ReferentievlakLijst'].str[0].values,
                         #'OpdrachtgevendeInstantie':measurements_wathte['WaarnemingMetadata.OpdrachtgevendeInstantieLijst'].str[0].values,
                         },
                        index=pd.to_datetime(ddlpy_meas['Tijdstip']))
    
    if isnumeric:
        ts_pd['values'] /= 100 #convert from cm to m

    # sort on time values # TODO: do this in ddlpy or in ddl
    ts_pd = ts_pd.sort_index()
    
    return ts_pd


######### oneline waterlevel data retrieval for one station
if 0: #for RWS
    import hatyan # available via `pip install hatyan` or at https://github.com/Deltares/hatyan
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
        ts_measwl = ddlpy_to_hatyan(meas_wathte_ts)
        meas_wathtbrkd_ts = meas_wathtbrkd.loc[meas_wathtbrkd['Groepering.Code'].isin(['NVT'])]
        ts_astro = ddlpy_to_hatyan(meas_wathtbrkd_ts)
        
        if include_extremes:
            #extremes
            meas_wathte_ext = meas_wathte.loc[meas_wathte['Groepering.Code'].isin(['GETETM2'])]
            ts_measwlHWLW = ddlpy_to_hatyan(meas_wathte_ext)
            meas_wathtbrkd_ext = meas_wathtbrkd.loc[meas_wathtbrkd['Groepering.Code'].isin(['GETETBRKD2'])]
            ts_astroHWLW = ddlpy_to_hatyan(meas_wathtbrkd_ext)
            # extreme types
            meas_wathte_exttype = meas_types.loc[meas_types['Groepering.Code'].isin(['GETETM2'])]
            ts_measwlHWLWtype = ddlpy_to_hatyan(meas_wathte_exttype)
            meas_wathtbrkd_exttype = meas_types.loc[meas_types['Groepering.Code'].isin(['GETETBRKD2'])]
            ts_astroHWLWtype = ddlpy_to_hatyan(meas_wathtbrkd_exttype)
            
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
        ts_measwl = ddlpy_to_hatyan(meas_wathte_ts)
        
        stat_name = locs_wathte_one["Naam"]
        stat_code = locs_wathte_one["Code"]
        fig, (ax1,ax2) = plt.subplots(2,1, figsize=(8,6), sharex=True)
        ts_measwl["values"].plot(ax=ax1)
        ts_measwl["QC"].plot(ax=ax2)
        ax1.set_title(f'waterlevels for {stat_code} ({stat_name})')
        ax2.set_title(f'QC for {stat_code} ({stat_name})')
        fig.tight_layout()
