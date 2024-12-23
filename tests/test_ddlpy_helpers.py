# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 14:24:07 2024

@author: veenstra
"""

import pytest
import hatyan
import ddlpy
import datetime as dt
import numpy as np
import pandas.api.types as ptypes


@pytest.fixture(scope="session")
def locations():
    """return all locations"""
    locations = ddlpy.locations()
    return locations


@pytest.mark.unittest
@pytest.mark.parametrize("grootheid", [pytest.param(gr, id=gr) for gr in ['WATHTE','WATHTBRKD']])
def test_ddlpy_to_hatyan(locations, grootheid):
    # input parameters
    tstart_dt = dt.datetime(2023,12,24)
    tstop_dt = dt.datetime(2024,1,5)

    bool_station = locations.index.isin(['HOEKVHLD'])
    bool_hoedanigheid = locations['Hoedanigheid.Code'].isin(['NAP'])
    bool_grootheid = locations['Grootheid.Code'].isin([grootheid])
    bool_groepering = locations['Groepering.Code'].isin(['NVT'])
    locs_wathte = locations.loc[bool_station & bool_grootheid &
                                bool_groepering & bool_hoedanigheid]

    meas_wathte = ddlpy.measurements(locs_wathte.iloc[0], start_date=tstart_dt, end_date=tstop_dt)
    
    # timeseries
    ts_measwl = hatyan.ddlpy_to_hatyan(meas_wathte)

    assert ts_measwl.columns.tolist() == ['values', 'qualitycode', 'status']
    assert ts_measwl.index.name == 'time'

    assert ptypes.is_float_dtype(ts_measwl["values"])
    assert ptypes.is_integer_dtype(ts_measwl["qualitycode"])
    assert ptypes.is_object_dtype(ts_measwl["status"])

    # check if metadata is complete and correct
    meta_fromts = ts_measwl.attrs
    meta_expected = {
        'grootheid': grootheid,
        'groepering': 'NVT',
        'eenheid': 'm',
        'vertref': 'NAP',
        'station': 'HOEKVHLD',
        'origin': 'ddlpy',
        }
    assert meta_fromts == meta_expected


@pytest.mark.unittest
def test_convert_hwlwstr2num(locations):
    tstart_dt = "2020-12-19"
    tstop_dt = "2020-12-31"
    
    bool_grootheid_meas = locations['Grootheid.Code'].isin(['WATHTE'])
    bool_hoedanigheid = locations['Hoedanigheid.Code'].isin(['NAP'])
    bool_groepering_ext = locations['Groepering.Code'].isin(['GETETM2','GETETMSL2'])
    bool_grootheid_exttypes = locations['Grootheid.Code'].isin(['NVT'])
    bool_groepering_ext_meas = locations['Groepering.Code'].isin(['GETETM2','GETETMSL2'])
    bool_station = locations.index.isin(["HOEKVHLD"])
    locs_wathte_ext = locations.loc[bool_grootheid_meas & bool_hoedanigheid & bool_groepering_ext & bool_station]
    locs_exttypes_wathte = locations.loc[bool_grootheid_exttypes & bool_groepering_ext_meas & bool_station]
    
    # no support for multiple rows, so pass one at a time
    meas_wathte_ext = ddlpy.measurements(locs_wathte_ext.iloc[0], start_date=tstart_dt, end_date=tstop_dt)
    meas_wathte_exttypes = ddlpy.measurements(locs_exttypes_wathte.iloc[0], start_date=tstart_dt, end_date=tstop_dt)
    # hatyan timeseries
    ts_measwlHWLW = hatyan.ddlpy_to_hatyan(meas_wathte_ext, meas_wathte_exttypes)
    
    hwlwcode_expected = np.array([2, 1, 3, 4, 5, 1, 2, 1, 3, 4, 5, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1,
           2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 3, 4,
           5, 1, 3, 4, 5, 1, 3, 4, 5, 1, 3, 4, 5, 1, 2])
    
    assert "HWLWcode" in ts_measwlHWLW.columns
    assert np.array_equal(ts_measwlHWLW["HWLWcode"].values, hwlwcode_expected)
    
    # check if metadata is complete and correct
    meta_fromts = ts_measwlHWLW.attrs
    meta_expected = {
        'grootheid': 'WATHTE',
        'groepering': 'GETETM2',
        'eenheid': 'm',
        'vertref': 'NAP',
        'station': 'HOEKVHLD',
        'origin': 'ddlpy',
        }
    assert meta_fromts == meta_expected


@pytest.mark.systemtest
def test_ddlpy_to_components(tmp_path, locations):
    """
    this test shows the metadata of the hatyan timeseries is properly created and passed
    such that the component file can be written without altering or adding metadata.
    This was achieved by fixing issue https://github.com/Deltares/hatyan/issues/357 
    and issue https://github.com/Deltares/hatyan/issues/358.
    """
    start_dt  = "2022-01-01 00:00:00 +01:00"
    end_dt  = "2022-03-31 23:50:00 +01:00"

    bool_grootheid = locations['Grootheid.Code'].isin(['WATHTE']) # measured waterlevels (not astro)
    bool_groepering = locations['Groepering.Code'].isin(['NVT']) # timeseries (not extremes)
    bool_hoedanigheid = locations['Hoedanigheid.Code'].isin(['NAP']) # vertical reference, only NAP
    locs_ddl = locations.loc[bool_grootheid & bool_groepering & bool_hoedanigheid]

    donar_loccode = "VLISSGN"
    locs_ddl_one = locs_ddl.loc[donar_loccode]
    ddl_df = ddlpy.measurements(locs_ddl_one, start_date=start_dt, end_date=end_dt)
    df_meas = hatyan.ddlpy_to_hatyan(ddl_df)
    
    # check if metadata is complete and correct
    meta_fromts = df_meas.attrs
    meta_expected = {
        'grootheid': 'WATHTE',
        'groepering': 'NVT',
        'eenheid': 'm',
        'vertref': 'NAP',
        'station': 'VLISSGN',
        'origin': 'ddlpy'
        }
    assert meta_fromts == meta_expected
    
    comp = hatyan.analysis(ts=df_meas, const_list='month', nodalfactors=True, fu_alltimes=False, xfac=True, analysis_perperiod=False)
    file_comp = tmp_path / f'components_{donar_loccode}.ana'
    hatyan.write_components(comp, filename=file_comp)
