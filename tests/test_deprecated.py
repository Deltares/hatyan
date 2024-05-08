# -*- coding: utf-8 -*-
"""
Created on Wed May  8 13:37:10 2024

@author: veenstra
"""

import pytest
import hatyan

list_deprecated_funcs = hatyan.deprecated.__all__

def test_deprecated_functions():
    for deprecated_func in list_deprecated_funcs:
        with pytest.raises(DeprecationWarning):
            exec(f"hatyan.{deprecated_func}()")
