# -*- coding: utf-8 -*-
"""
ddlpy_helpers.py contains functions to convert ddlpy timeseries dataframes to hatyan timeseries dataframes.
The package ddlpy is available at https://github.com/Deltares/ddlpy
"""

import pandas as pd
from hatyan.metadata import metadata_from_ddlpy

__all__ = ["ddlpy_to_hatyan"]


def ddlpy_to_hatyan(ddlpy_meas, ddlpy_meas_exttyp=None):
    """
    Convert ddlpy measurements to hatyan timeseries dataframe.

    Parameters
    ----------
    ddlpy_meas : pd.DataFrame
        ddlpy measurements dataframe. It is assumed that it contains numeric values that 
        represent waterlevel timeseries or waterlevel extremes (measured or astro).
    ddlpy_meas_exttyp : pd.DataFrame, optional
        ddlpy measurements dataframe. If it is supplied it is assumed that it contains 
        alfanumeric values and these represent tidal extreme types (high- and low waters).
        The default is None.

    Returns
    -------
    pd.DataFrame
        hatyan timeseries DataFrame with values/qualitycode/status columns. If ddlpy_meas_typ 
        is supplied, this DataFrame will also include a HWLWcode column.

    """
    
    ts_pd = ddlpy_to_hatyan_plain(ddlpy_meas, isnumeric=True)
    metadata = metadata_from_ddlpy(ddlpy_meas)
    # conver units from cm to meters
    assert metadata['eenheid'] == 'cm'
    ts_pd['values'] /= 100 #convert from cm to m
    metadata['eenheid'] = 'm'
    ts_pd.attrs = metadata
    if ddlpy_meas_exttyp is None:
        return ts_pd
    
    # check if the contents of this dataframe is in deed extreme types
    assert len(ddlpy_meas) == len(ddlpy_meas_exttyp)
    typering_codes = ddlpy_meas_exttyp["Typering.Code"].drop_duplicates().values
    assert len(typering_codes) == 1
    assert typering_codes[0] == "GETETTPE"
    
    ts_pd_typ = ddlpy_to_hatyan_plain(ddlpy_meas_exttyp, isnumeric=False)
    
    ts_pd_combined = convert_exttype_str2num(ts_pd,ts_pd_typ)
    
    return ts_pd_combined
        

def ddlpy_to_hatyan_plain(ddlpy_meas, isnumeric=True):
    
    if ddlpy_meas.empty:
        raise ValueError("supplied dataframe is empty")
    
    cols_mustbe_unique = ['Grootheid.Code','Groepering.Code','Typering.Code','Hoedanigheid.Code']
    for col in cols_mustbe_unique:
        if len(ddlpy_meas[col].drop_duplicates()) != 1:
            raise Exception(f"ddlpy_meas['{col}'] is not unique")
    
    if isnumeric: 
        key_numericvalues = 'Meetwaarde.Waarde_Numeriek'
    else:
        # alfanumeric values for 'Typering.Code':'GETETTPE'
        key_numericvalues = 'Meetwaarde.Waarde_Alfanumeriek'
    
    ts_pd = pd.DataFrame({'values':ddlpy_meas[key_numericvalues],
                          'qualitycode':pd.to_numeric(ddlpy_meas['WaarnemingMetadata.KwaliteitswaardecodeLijst'],downcast='integer'),
                          'status':ddlpy_meas['WaarnemingMetadata.StatuswaardeLijst'].str[0],
                          })
    
    return ts_pd


def convert_exttype_str2num(ts_measwl_ext, ts_measwl_exttype):
    """
    TVL;1;1;hoogwater
    TVL;1;2;laagwater
    TVL;1;3;laagwater 1
    TVL;1;4;topagger
    TVL;1;5;laagwater 2
    """
    
    ts_measwl_ext.loc[ts_measwl_exttype['values']=='hoogwater','HWLWcode'] = 1
    ts_measwl_ext.loc[ts_measwl_exttype['values']=='laagwater','HWLWcode'] = 2
    ts_measwl_ext.loc[ts_measwl_exttype['values']=='laagwater 1','HWLWcode'] = 3
    ts_measwl_ext.loc[ts_measwl_exttype['values']=='topagger','HWLWcode'] = 4
    ts_measwl_ext.loc[ts_measwl_exttype['values']=='laagwater 2','HWLWcode'] = 5
    ts_measwl_ext['HWLWcode'] = ts_measwl_ext['HWLWcode'].astype(int)
    return ts_measwl_ext

