# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 20:48:46 2023

@author: veenstra
"""

import numpy as np
import pytz


def metadata_add_to_obj(obj,metadata):
    """
    add metadata to components/timeseries object (is pd.DataFrame)
    """
    for key in metadata.keys():
        obj.__setattr__(key,metadata[key])
    
    return obj


def metadata_from_diablocks(diablocks_pd, block_id):
    diablocks_pd_onerow = diablocks_pd.iloc[block_id]
    
    metadata_keys = ['station', 'grootheid', 'eenheid', 
                     'vertref', 
                     'waarnemingssoort', 
                     'tstart', 'tstop',
                     'timestep_min', 'timestep_unit',
                     'TYP', 'groepering']
    
    #TODO: align with metadata from hatyan.read_components()
    metadata = {key:diablocks_pd_onerow[key] for key in metadata_keys}
    
    # add hardcodedtzone
    #TODO: is this documented in the dia file?
    metadata['tzone'] = pytz.FixedOffset(60)
    
    # add origin
    metadata['origin'] = 'from timeseries dia file'
    
    # replace nan with None (otherwise metadata_compare fails)
    #TODO: avoid nan in metadata (timestep for hoek_har.dia)
    if np.isnan(metadata['timestep_min']):
        metadata['timestep_min'] = None
        metadata['timestep_unit'] = None
    return metadata


def metadata_from_ddlmeta(metadata_ddl, stationdata_ddl):
    
    metadata = {}
    metadata['station'] = stationdata_ddl.iloc[0]['Code']
    metadata['grootheid'] = metadata_ddl.iloc[0]['Grootheid.Code']
    metadata['eenheid'] = metadata_ddl.iloc[0]['Eenheid.Code']
    metadata['vertref'] = metadata_ddl.iloc[0]['Hoedanigheid.Code']
    # metadata['waarnemingssoort'] = metadata_ddl.iloc[0]['']
    # metadata['tstart'] = metadata_ddl.iloc[0]['']
    # metadata['tstop'] = metadata_ddl.iloc[0]['']
    # metadata['timestep_min'] = metadata_ddl.iloc[0]['']
    # metadata['timestep_unit'] = metadata_ddl.iloc[0]['']
    # metadata['TYP'] = metadata_ddl.iloc[0]['Typering.Code']
    metadata['groepering'] = metadata_ddl.iloc[0]['Groepering.Code']
    
    return metadata


def metadata_from_obj(obj):
    obj_vars = vars(obj)
    metadata = {key:obj_vars[key] for key in obj_vars.keys() if not key.startswith('_')}
    if len(metadata) == 0:
        raise ValueError('no metadata found on object')
        # metadata = {}
    return metadata


def metadata_compare(metadata_list):
    
    # remove tstart/tstop since they cannot be compared on equality in case of multifile dia
    metadata_list_notstartstop = []
    for meta in metadata_list:
        meta_new = meta.copy()
        meta_new.pop('tstart')
        meta_new.pop('tstop')
        metadata_list_notstartstop.append(meta_new)
    
    nmeta = len(metadata_list_notstartstop)
    for i in range(1,nmeta):
        meta1 = metadata_list_notstartstop[i-1]
        meta2 = metadata_list_notstartstop[i]
        if meta1!=meta2:
            raise ValueError(f'metadata for two datasets is not equal, cannot be merged:\n{meta1}\n{meta2}')
