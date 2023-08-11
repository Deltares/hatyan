# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 04:29:26 2023

@author: veenstra
"""


def wns_from_metadata(metadata):
    """
    Gets waarnemingssoort from metadata grootheid (quantity), eenheid (unit) and vertref (vertical reference).
    Originally done in: https://repos.deltares.nl/repos/lib_tide/trunk/src/hatyan_fortran/HATYAN00/wnstab.f
    Deliberately left out domeincode (type), I(int) and F(float) since we always have floats 
    that are only converted to int upon file writing
    """
    
    meta_sel = {key:metadata[key] for key in ['grootheid','eenheid','vertref']}
    
    if meta_sel == {'grootheid':'WATHE', 'eenheid':'cm', 'vertref':'NAP'}:
        wns = 1
    elif meta_sel == {'grootheid':'WATHE', 'eenheid':'cm', 'vertref':'MSL'}:
        wns = 54
    elif meta_sel == {'grootheid':'WATHTBRKD', 'eenheid':'cm', 'vertref':'NAP'}:
        wns = 18
    elif meta_sel == {'grootheid':'WATHTBRKD', 'eenheid':'cm', 'vertref':'MSL'}:
        wns = 55
    else:
        raise ValueError(f'combination of quantity/unit/vertref not found available in wns_from_metadata():\n{wns_from_metadata}')
    
    return wns
    


