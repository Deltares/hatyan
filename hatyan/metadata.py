# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 20:48:46 2023

@author: veenstra
"""

def metadata_add_to_obj(obj,metadata):
    for key in metadata.keys():
        obj.__setattr__(key,metadata[key])
    
    return obj

def metadata_from_diablocks(diablocks_pd, block_id):
    diablocks_pd_onerow = diablocks_pd.iloc[block_id]
    metadata_keys = ['vertref', 
                     'station', 'TYP', 'groepering', 'grootheid', 'eenheid', 
                     'tstart', 'tstop', 'timestep_min', 'timestep_unit']
    metadata = {key:diablocks_pd_onerow[key] for key in metadata_keys}
    return metadata

def metadata_from_obj(obj):
    obj_vars = vars(obj)
    metadata = {key:obj_vars[key] for key in obj_vars.keys() if not key.startswith('_')}
    return metadata
