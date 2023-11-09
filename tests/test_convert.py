# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 12:13:40 2023

@author: veenstra
"""

import pytest
import numpy as np
import hatyan
from hatyan.convert import convert_tzone2tzinfo


@pytest.mark.unittest
def test_convert_tzone2tzinfo():
    tzone = "UTC+01:00"
    tzinfo = convert_tzone2tzinfo(tzone)
    assert tzinfo.zone == 'Etc/GMT-1'

    
@pytest.mark.unittest
def test_convert_coordinates():
    values_x = np.array([576917.66978449, 571218.73553018, 541425.98321488, 661021.58550526])
    values_y = np.array([5759136.15818497, 5742305.0566163 , 5699181.90968435,
           5894519.40967015])
    
    values_x_RD_exp = np.array([ 67930.00003341,  61679.99979896,  30479.9991833 , 156480.00176811])
    values_y_RD_exp = np.array([444000.00275723, 427360.00284566, 385220.0032781 , 576550.00145627])
    
    values_x_RD_int, values_y_RD_int = hatyan.convert_coordinates(coordx_in=values_x, coordy_in=values_y, epsg_in=25831,epsg_out=28992)
    values_x_RD_str, values_y_RD_str = hatyan.convert_coordinates(coordx_in=values_x, coordy_in=values_y, epsg_in=25831,epsg_out="RD")
    
    assert np.allclose(values_x_RD_int, values_x_RD_str, values_x_RD_exp)
    assert np.allclose(values_y_RD_int, values_y_RD_str, values_y_RD_exp)
