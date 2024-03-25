# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 09:25:17 2024

@author: veenstra
"""

import pytest
import subprocess

@pytest.mark.unittest
def test_cli():
    p = subprocess.Popen("python -m hatyan --help")
