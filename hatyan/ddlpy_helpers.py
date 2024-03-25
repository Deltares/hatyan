# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 17:03:03 2021

@author: veenstra
"""

import pandas as pd

__all__ = ["ddlpy_to_hatyan",
           "convert_HWLWstr2num"]


def ddlpy_to_hatyan(ddlpy_meas):
    
    if ddlpy_meas.empty:
        raise ValueError("supplied dataframe is empty")
    
    cols_mustbe_unique = ['Grootheid.Code','Groepering.Code','Typering.Code','Hoedanigheid.Code']
    for col in cols_mustbe_unique:
        if len(ddlpy_meas[col].drop_duplicates()) != 1:
            raise Exception(f"ddlpy_meas['{col}'] is not unique")
    
    if 'Meetwaarde.Waarde_Numeriek' in ddlpy_meas.columns: 
        key_numericvalues = 'Meetwaarde.Waarde_Numeriek'
        isnumeric = True
    else:
        #alfanumeric values for 'Typering.Code':'GETETTPE' #DDL IMPROVEMENT: also include numeric values for getijtype. Also, it is quite complex to get this data in the first place, would be convenient if it would be a column when retrieving 'Groepering.Code':'GETETM2' or 'GETETBRKD2'
        key_numericvalues = 'Meetwaarde.Waarde_Alfanumeriek'
        isnumeric = False
    
    # DDL IMPROVEMENT: qc conversion should be possible with .astype(int), but pd.to_numeric() is necessary for HARVT10 (eg 2019-09-01 to 2019-11-01) since QC contains None values that cannot be ints (in that case array of floats with some nans is returned) >> now flattened by ddlpy
    ts_pd = pd.DataFrame({'values':ddlpy_meas[key_numericvalues],
                         'QC':pd.to_numeric(ddlpy_meas['WaarnemingMetadata.KwaliteitswaardecodeLijst'],downcast='integer'),
                         'Status':ddlpy_meas['WaarnemingMetadata.StatuswaardeLijst'],
                         })
    
    if isnumeric:
        ts_pd['values'] /= 100 #convert from cm to m
    
    return ts_pd


def convert_HWLWstr2num(ts_measwlHWLW,ts_measwlHWLWtype):
    """
    TVL;1;1;hoogwater
    TVL;1;2;laagwater
    TVL;1;3;laagwater 1
    TVL;1;4;topagger
    TVL;1;5;laagwater 2
    """
    ts_measwlHWLW.loc[ts_measwlHWLWtype['values']=='hoogwater','HWLWcode'] = 1
    ts_measwlHWLW.loc[ts_measwlHWLWtype['values']=='laagwater','HWLWcode'] = 2
    ts_measwlHWLW.loc[ts_measwlHWLWtype['values']=='laagwater 1','HWLWcode'] = 3
    ts_measwlHWLW.loc[ts_measwlHWLWtype['values']=='topagger','HWLWcode'] = 4
    ts_measwlHWLW.loc[ts_measwlHWLWtype['values']=='laagwater 2','HWLWcode'] = 5
    ts_measwlHWLW['HWLWcode'] = ts_measwlHWLW['HWLWcode'].astype(int)
    return ts_measwlHWLW

