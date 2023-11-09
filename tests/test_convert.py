# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 12:13:40 2023

@author: veenstra
"""

import os
import pytest
import numpy as np
import hatyan
from hatyan.convert import convert_tzone2tzinfo

dir_tests = os.path.dirname(__file__) #F9 doesnt work, only F5 (F5 also only method to reload external definition scripts)
dir_testdata = os.path.join(dir_tests,'data_unitsystemtests')


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


@pytest.mark.unittest
def test_convert_coordinates_diafile():
    file_pred = os.path.join(dir_testdata, "VLISSGN_pre.txt")

    diablocks_pd_extra = hatyan.get_diablocks(filename=file_pred)
    coordx_in = diablocks_pd_extra.loc[0,'x']
    coordy_in = diablocks_pd_extra.loc[0,'y']
    epsg_in = diablocks_pd_extra.loc[0,'epsg']
    WGS84x_int, WGS84y_int = hatyan.convert_coordinates(coordx_in=coordx_in, coordy_in=coordy_in, epsg_in=epsg_in, epsg_out=4326)
    WGS84x_str, WGS84y_str = hatyan.convert_coordinates(coordx_in=coordx_in, coordy_in=coordy_in, epsg_in="RD", epsg_out="W84")

    assert np.isclose(WGS84x_int, WGS84x_str, 3.596056794834692)
    assert np.isclose(WGS84y_int, WGS84y_str, 51.44231093287052)
