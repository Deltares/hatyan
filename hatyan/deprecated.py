# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 16:29:16 2023

@author: veenstra
"""

from typing import Any, Callable, Dict
import functools
import logging

__all__ = ["get_components_from_ts",
           "check_ts"]

logger = logging.getLogger(__name__)


DEPRECATED_OPTIONS_PREDICTION_DICT = {
    'times_pred_all':'times',
    'times_ext':'times',
    'timestep_min':'timestep',
    'nodalfactors':'as attribute of the component dataframe (comp.attrs["nodalfactors"])',
    'xfac':'as attribute of the component dataframe (comp.attrs["xfac"])',
    'fu_alltimes':'as attribute of the component dataframe (comp.attrs["fu_alltimes"])',
    'source':'as attribute of the component dataframe (comp.attrs["source"])',
    }


def deprecated_python_option(**aliases: str) -> Callable:
    def deco(f: Callable):
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            check_old_kwargs(f.__name__, kwargs, aliases)
            return f(*args, **kwargs)

        return wrapper

    return deco


def check_old_kwargs(
    func_name: str, kwargs: Dict[str, Any], aliases: Dict[str, str]
):
    for alias, new in aliases.items():
        if alias in kwargs:
            message = f"Argument '{alias}' has been deprecated for hatyan.{func_name}()"
            if new is not None:
                message += f", use '{new}' instead"
            raise DeprecationWarning(message)


def get_components_from_ts(**kwargs):
    raise DeprecationWarning("hatyan.get_components_from_ts() was deprecated, use hatyan.analysis() instead")


def check_ts(**kwargs):
    raise DeprecationWarning("hatyan.check_ts() was deprecated, use hatyan.Timeseries_Statistics() instead")
