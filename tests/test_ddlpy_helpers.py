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


@pytest.mark.unittest
def test_ddlpy_to_hatyan():
    # input parameters
    tstart_dt = dt.datetime(2023,12,24)
    tstop_dt = dt.datetime(2024,1,5)

    locations = ddlpy.locations()

    bool_station = locations.index.isin(['HOEKVHLD'])
    bool_hoedanigheid = locations['Hoedanigheid.Code'].isin(['NAP'])
    bool_grootheid = locations['Grootheid.Code'].isin(['WATHTE'])
    bool_groepering = locations['Groepering.Code'].isin(['NVT'])
    locs_wathte = locations.loc[bool_station & bool_grootheid &
                                bool_groepering & bool_hoedanigheid]

    meas_wathte = ddlpy.measurements(locs_wathte.iloc[0], start_date=tstart_dt, end_date=tstop_dt)
    
    # timeseries
    ts_measwl = hatyan.ddlpy_to_hatyan(meas_wathte)

    assert ts_measwl.columns.tolist() == ['values', 'QC', 'Status']
    assert ts_measwl.index.name == 'time'

    assert ptypes.is_float_dtype(ts_measwl["values"])
    assert ptypes.is_integer_dtype(ts_measwl["QC"])
    assert ptypes.is_object_dtype(ts_measwl["Status"])
    