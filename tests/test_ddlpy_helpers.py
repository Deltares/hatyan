# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 14:24:07 2024

@author: veenstra
"""

import pytest
import hatyan
import ddlpy
import datetime as dt
import pandas.api.types as ptypes
import numpy as np


@pytest.mark.acceptance
def test_DDL_QCvalues():
    tstart_dt = dt.datetime(2019,10,1)
    tstop_dt = dt.datetime(2019,10,10)
    tzone = 'UTC+00:00' #'UTC+00:00' for GMT and 'UTC+01:00' for MET
    
    catalog_dict = hatyan.get_DDL_catalog(catalog_extrainfo=['WaardeBepalingsmethoden','MeetApparaten','Typeringen'])
    
    #HARVT10
    cat_locatielijst = catalog_dict['LocatieLijst']#.set_index('Locatie_MessageID',drop=True)
    station_dict = cat_locatielijst[cat_locatielijst['Code']=='HARVT10'].iloc[0]
    ts_meas_pd, metadata, stationdata = hatyan.get_DDL_data(station_dict=station_dict,tstart_dt=tstart_dt,tstop_dt=tstop_dt,tzone=tzone,
                                                            meta_dict={'Grootheid.Code':'WATHTE','Groepering.Code':'NVT','WaardeBewerkingsmethode.Code':'NVT'})
    uniqueQC = ts_meas_pd['QC'].unique() #array([ 0., 99., nan, 25.]), but should be integers without nan array([0, 99, 25])
    assert uniqueQC.dtype=='float64' #this should be int in the future, if None/nan is not in QC list anymore
    assert np.isnan(uniqueQC).sum() == 1 #this one should become 0 in the future and then the second assertion should be valid without indexing
    assert (uniqueQC[~np.isnan(uniqueQC)] == np.array([ 0, 99, 25])).all()

    #STELLDBTN
    station_dict = cat_locatielijst[cat_locatielijst['Code']=='STELLDBTN'].iloc[0]
    ts_meas_pd, metadata, stationdata = hatyan.get_DDL_data(station_dict=station_dict,tstart_dt=tstart_dt,tstop_dt=tstop_dt,tzone=tzone,
                                                            meta_dict={'Grootheid.Code':'WATHTE','Groepering.Code':'NVT','WaardeBewerkingsmethode.Code':'NVT'})
    uniqueQC = ts_meas_pd['QC'].unique() #array([ 0, 25], dtype=int8)
    assert uniqueQC.dtype=='int8'
    assert (uniqueQC == np.array([ 0, 25])).all()


@pytest.mark.unittest
def test_ddlpy_to_hatyan():
    # input parameters
    tstart_dt = dt.datetime(2019,12,24) #period begins with Gecontroleerd and ends with Ongecontroleerd for HOEKVHLD
    tstop_dt = dt.datetime(2020,1,5)

    locations = ddlpy.locations()

    bool_station = locations.index.isin(['HOEKVHLD'])
    bool_hoedanigheid = locations['Hoedanigheid.Code'].isin(['NAP'])
    bool_grootheid = locations['Grootheid.Code'].isin(['WATHTE'])
    locs_wathte = locations.loc[bool_station & bool_grootheid & bool_hoedanigheid]

    meas_wathte = ddlpy.measurements(locs_wathte.iloc[0], start_date=tstart_dt, end_date=tstop_dt)

    # TODO: rename lowercase code to uppercase Code
    meas_wathte.columns = [x.replace(".code",".Code") for x in meas_wathte.columns]

    # timeseries
    meas_wathte_ts = meas_wathte.loc[meas_wathte['Groepering.Code'].isin(['NVT'])]
    ts_measwl = hatyan.ddlpy_to_hatyan(meas_wathte_ts)

    assert ts_measwl.columns.tolist() == ['values', 'QC', 'Status']
    assert ts_measwl.index.name == 'Tijdstip'

    assert ptypes.is_float_dtype(ts_measwl["values"])
    assert ptypes.is_integer_dtype(ts_measwl["QC"])
    assert ptypes.is_object_dtype(ts_measwl["Status"])
    