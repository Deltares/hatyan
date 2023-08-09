# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 20:48:46 2023

@author: veenstra
"""

import numpy as np


def metadata_add_to_obj(obj,metadata):
    """
    add metadata to components/timeseries object (is pd.DataFrame)
    """
    for key in metadata.keys():
        obj.__setattr__(key,metadata[key])
    
    return obj


def metadata_from_diablocks(diablocks_pd, block_id):
    diablocks_pd_onerow = diablocks_pd.iloc[block_id]
    metadata_keys = ['vertref', 
                     'station', 'TYP', 'groepering', 'grootheid', 'eenheid', 
                     #'tstart', 'tstop', # cannot be compared on equality in case of multifile dia
                     'timestep_min', 'timestep_unit']
    metadata = {key:diablocks_pd_onerow[key] for key in metadata_keys}
    
    # replace nan with None (otherwise metadata_compare fails)
    #TODO: avoid nan in metadata (timestep for hoek_har.dia)
    if np.isnan(metadata['timestep_min']):
        metadata['timestep_min'] = None
        metadata['timestep_unit'] = None
    return metadata


def metadata_from_obj(obj):
    obj_vars = vars(obj)
    metadata = {key:obj_vars[key] for key in obj_vars.keys() if not key.startswith('_')}
    if len(metadata) == 0:
        raise ValueError('no metadata found on object')
    return metadata


def metadata_compare(metadata_list):
    nmeta = len(metadata_list)
    for i in range(1,nmeta):
        meta1 = metadata_list[i-1]
        meta2 = metadata_list[i]
        if meta1!=meta2:
            raise ValueError(f'metadata for two datasets is not equal, cannot be merged:\n{meta1}\n{meta2}')
