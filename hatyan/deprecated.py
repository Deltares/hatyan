# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 16:29:16 2023

@author: veenstra
"""

__all__ = ["get_components_from_ts",
           "check_ts"]


def get_components_from_ts(**kwargs):
    raise DeprecationWarning("hatyan.get_components_from_ts() was deprecated, use hatyan.analysis() instead")

def check_ts(**kwargs):
    raise DeprecationWarning("hatyan.check_ts() was deprecated, use hatyan.Timeseries_Statistics() instead")
