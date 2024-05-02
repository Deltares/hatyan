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
    via attrs, but this function prevents dropping existing attrs that are not in metadata dict
    """
    for key in metadata.keys():
        obj.attrs[key] = metadata[key]
    
    return obj


def metadata_from_diablocks(diablocks_pd, block_id):
    diablocks_pd_onerow = diablocks_pd.iloc[block_id]
    
    metadata_keys = ['station', 'grootheid', 'eenheid', 
                     'vertref', 
                     'tstart', 'tstop',
                     'timestep_min', 'timestep_unit',
                     'TYP', 'groepering']
    
    #TODO: align with metadata from hatyan.read_components()
    metadata = {key:diablocks_pd_onerow[key] for key in metadata_keys}
    
    # add origin
    metadata['origin'] = 'from timeseries dia file'
    
    # replace nan with None (otherwise metadata_compare fails)
    #TODO: avoid nan in metadata (timestep for hoek_har.dia)
    if np.isnan(metadata['timestep_min']): #non-equidistant, nan in py38 and none in py39 (pandas 2.1.2)
        metadata['timestep_min'] = None
        metadata['timestep_unit'] = None
    return metadata


def metadata_from_obj(obj):
    metadata = obj.attrs.copy()
    return metadata


def metadata_compare(metadata_list):
    
    # remove tstart/tstop since they cannot be compared on equality in case of multifile dia
    metadata_list_notstartstop = []
    for meta in metadata_list:
        meta_new = meta.copy()
        if 'tstart' in meta_new:
            meta_new.pop('tstart')
        if 'tstop' in meta_new:
            meta_new.pop('tstop')
        metadata_list_notstartstop.append(meta_new)
    
    nmeta = len(metadata_list_notstartstop)
    for i in range(1,nmeta):
        meta1 = metadata_list_notstartstop[i-1]
        meta2 = metadata_list_notstartstop[i]
        if meta1!=meta2:
            raise ValueError(f'metadata for two datasets is not equal, cannot be merged:\n{meta1}\n{meta2}')


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
