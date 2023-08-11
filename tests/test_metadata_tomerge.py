# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 04:44:30 2023

@author: veenstra
"""

import pytest
from hatyan.metadata_tomerge import wns_from_metadata


@pytest.mark.unittest
def test_wns_from_metadata():
    metadata_1 = {'grootheid':'WATHE', 'eenheid':'cm', 'vertref':'NAP'}
    metadata_54 = {'grootheid':'WATHE', 'eenheid':'cm', 'vertref':'MSL'}
    metadata_18 = {'grootheid':'WATHTBRKD', 'eenheid':'cm', 'vertref':'NAP'}
    metadata_55 = {'grootheid':'WATHTBRKD', 'eenheid':'cm', 'vertref':'MSL'}
    
    wns_1 = wns_from_metadata(metadata_1)
    wns_54 = wns_from_metadata(metadata_54)
    wns_18 = wns_from_metadata(metadata_18)
    wns_55 = wns_from_metadata(metadata_55)
    
    assert wns_1 == 1
    assert wns_54 == 54
    assert wns_18 == 18
    assert wns_55 == 55
