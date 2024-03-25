# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 22:50:18 2022

@author: veenstra
"""

from pyproj import Transformer

__all__ = ["convert_coordinates"]


def convert_coordinates(coordx_in, coordy_in, epsg_in, epsg_out=28992):
    
    transformer = Transformer.from_crs(f'epsg:{epsg_in}', f'epsg:{epsg_out}', always_xy=True)
    coordx_out, coordy_out = transformer.transform(coordx_in, coordy_in)
    
    return coordx_out, coordy_out

