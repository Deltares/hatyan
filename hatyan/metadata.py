# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 20:48:46 2023

@author: veenstra
"""

import pandas as pd


def metadata_add_to_obj(obj,metadata):
    """
    add metadata to components/timeseries object (is pd.DataFrame)
    via attrs, but this function prevents dropping existing attrs that are not in metadata dict
    """
    for key in metadata.keys():
        obj.attrs[key] = metadata[key]
    
    return obj


def metadata_from_diablocks(diablocks_pd, block_id):
    diablocks_pd_onerow = diablocks_pd.iloc[block_id]
    
    metadata_keys = ['station', 'grootheid', 'eenheid', 
                     'vertref', 'TYP', 'groepering']
    
    #TODO: align with metadata from hatyan.read_components()
    metadata = {key:diablocks_pd_onerow[key] for key in metadata_keys}
    
    # add origin
    metadata['origin'] = 'from timeseries dia file'
    return metadata


def metadata_from_obj(obj):
    metadata = obj.attrs.copy()
    return metadata


def metadata_compare(metadata_list):
    nmeta = len(metadata_list)
    for i in range(1,nmeta):
        meta1 = metadata_list[i-1]
        meta2 = metadata_list[i]
        if meta1!=meta2:
            meta_12_df = pd.concat([pd.Series(meta1), pd.Series(meta2)],axis=1)
            meta_12_df["equal"] = meta_12_df[0]==meta_12_df[1]
            raise ValueError(f'metadata for two datasets is not equal, cannot be merged:\n{meta_12_df}')


def wns_from_metadata(metadata):
    """
    Gets waarnemingssoort from metadata grootheid (quantity), eenheid (unit) and vertref (vertical reference).
    Idea from: https://repos.deltares.nl/repos/lib_tide/trunk/src/hatyan_fortran/HATYAN00/wnstab.f
    Deliberately left out domeincode (type), I(int) and F(float) since we always have floats 
    that are only converted to int upon file writing
    """
    
    meta_sel = {key:metadata[key] for key in ['grootheid','eenheid','vertref']}
    
    if meta_sel == {'grootheid':'WATHTE', 'eenheid':'cm', 'vertref':'NAP'}:
        wns = 1
    elif meta_sel == {'grootheid':'WATHTE', 'eenheid':'cm', 'vertref':'MSL'}:
        wns = 54
    elif meta_sel == {'grootheid':'WATHTBRKD', 'eenheid':'cm', 'vertref':'NAP'}:
        wns = 18
    elif meta_sel == {'grootheid':'WATHTBRKD', 'eenheid':'cm', 'vertref':'MSL'}:
        wns = 55
    else:
        raise ValueError(f'combination of quantity/unit/vertref not found available in wns_from_metadata():\n{meta_sel}')
    
    return wns
