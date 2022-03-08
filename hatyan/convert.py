# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 22:50:18 2022

@author: veenstra
"""
from pyproj import Transformer
import pytz

def convert_coordinates(coordx_in, coordy_in, epsg_in, epsg_out=28992):

    epsg_dict = {'RD':28992,'W84':4326,'E50':4230}
    
    if isinstance(epsg_in,str):
        if epsg_in not in epsg_dict.keys():
            raise Exception('when providing epsg_in as a string, the options are: %s'%(list(epsg_dict.keys())))
        else:
            epsgcode_in = epsg_dict[epsg_in]
    else:
        epsgcode_in = epsg_in
    
    if isinstance(epsg_out,str):
        if epsg_out not in epsg_dict.keys():
            raise Exception('when providing epsg_out as a string, the options are: %s'%(list(epsg_dict.keys())))
        else:
            epsgcode_out = epsg_dict[epsg_out]
    else:
        epsgcode_out = epsg_out
        
    transformer = Transformer.from_crs('epsg:%i'%(epsgcode_in), 'epsg:%i'%(epsgcode_out), always_xy=True)
    coordx_out, coordy_out = transformer.transform(coordx_in, coordy_in)
    
    return coordx_out, coordy_out


def convert_tzone2tzinfo(tzone):

    if tzone.startswith('UTC+') or tzone.startswith('UTC-'): #parse to fixed offset like 'Etc/GMT-1'. +/- are counter intuitive but it works: https://pvlib-python.readthedocs.io/en/stable/timetimezones.html#fixedoffsets)
        if len(tzone)!=9 or not tzone.endswith(':00'):
            raise Exception('if tzone starts with UTC+ or UTC-, the string should be 9 characters long and have 0 minutes, like "UTC+01:00"')
        tzone_hr = int(tzone[4:6])
        if tzone[3]=='+':
            tzone = 'Etc/GMT-%d'%(tzone_hr)
        else:
            tzone = 'Etc/GMT+%d'%(tzone_hr)
    tzinfo = pytz.timezone(tzone)
    return tzinfo


