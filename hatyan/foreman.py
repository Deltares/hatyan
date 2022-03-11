# -*- coding: utf-8 -*-
"""
foreman.py contains all foreman definitions now embedded in hatyan. The dataset is derived from "M.G.G. Foreman (2004), Manual for Tidal Heights Analysis and Prediction, Institute of Ocean Sciences (Patricia Bay, Sidney B.C. Canada)"

hatyan is a Python program for tidal analysis and prediction, based on the FORTRAN version. 
Copyright (C) 2019-2021 Rijkswaterstaat.  Maintained by Deltares, contact: Jelmer Veenstra (jelmer.veenstra@deltares.nl). 
Source code available at: https://github.com/Deltares/hatyan

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""

import os
import pandas as pd
import numpy as np
import functools
import datetime as dt

file_path = os.path.realpath(__file__)


#################################################
################# FILECONTENTS ##################
#################################################


@functools.lru_cache() #only caching this already makes foreman slightly faster
def get_foreman_doodson_nodal_harmonic(lat_deg=51.45):
    """
    Omzetten van het tweede deel van de foremantabel in een pandas DataFrame met harmonische (+satellite) componenten.

    Parameters
    ----------
    const_list : TYPE
        DESCRIPTION.
    lat_deg : TYPE, optional
         0 in degrees from equator (pos N, neg S). For R1 and R2, maar variatie heeft niet erg veel invloed op A en phi. The default is 51.45.

    Returns
    -------
    foreman_harmonic_doodson_all : TYPE
        DESCRIPTION.
    foreman_harmonic_nodal_all : TYPE
        DESCRIPTION.
    foreman_harmonic_doodson : TYPE
        DESCRIPTION.
    foreman_harmonic_nodal : TYPE
        DESCRIPTION.

    """
    
    foreman_file = os.path.join(os.path.dirname(file_path),'data','data_foreman_harmonic.txt')
    foreman_harmonic_raw = pd.read_csv(foreman_file, comment='#', names=range(16), index_col=0, skip_blank_lines=True, delim_whitespace=True)
    bool_dupl_index = foreman_harmonic_raw.index.duplicated(keep='first')
    
    lat_rad = np.deg2rad(lat_deg)
    
    R1 = 0.36309*(1.-5.*np.sin(lat_rad)*np.sin(lat_rad))/np.sin(lat_rad) #-1 #for lat=50N
    R2 = 2.59808*np.sin(lat_rad) #2 #for lat=50N
    
    #get all forman harmonic doodson values (first occurences of duplicated component names)
    foreman_doodson_harmonic_lun = foreman_harmonic_raw.loc[~bool_dupl_index].dropna(axis=1)
    foreman_doodson_harmonic_lun[5] = foreman_doodson_harmonic_lun[5].astype(float) #for some reason, column 5/N remains dtype object
    if not len(foreman_doodson_harmonic_lun.columns) == 8:
        raise Exception('unexpected amount of columns for harmonics, inconsistent file')
    foreman_doodson_harmonic_lun.columns = ['T','S','H','P','N','P1','EDN','nsats']
    foreman_doodson_harmonic_lun.index.name = None
    
    #convert from lunar to solar and reverse automatic column sorting
    omega1 = foreman_doodson_harmonic_lun.loc[:,'T']
    corr_array = pd.DataFrame({'S':-omega1,'H':omega1})
    foreman_doodson_harmonic = foreman_doodson_harmonic_lun.add(corr_array,fill_value=0)
    foreman_doodson_harmonic = foreman_doodson_harmonic[foreman_doodson_harmonic_lun.columns]
    
    #get all foreman harmonic nodal values (non-first occurences of duplicated component names)
    foreman_nodal_harmonic_wide = foreman_harmonic_raw.loc[bool_dupl_index]
    #reshape from 15 to 5 columns and drop nan lines
    foreman_nodal_harmonic = pd.DataFrame(foreman_nodal_harmonic_wide.values.reshape(207,5),
                                          index=foreman_nodal_harmonic_wide.index.repeat(3),
                                          columns=['P','N','P1','EDN','factor']).dropna(axis=0)
    
    #multiply with R1/R2 if applicable and convert entire dataframe to floats
    bool_R1 = foreman_nodal_harmonic['factor'].str.contains('R1')
    bool_R2 = foreman_nodal_harmonic['factor'].str.contains('R2')
    foreman_nodal_harmonic.loc[bool_R1,'factor'] = foreman_nodal_harmonic.loc[bool_R1,'factor'].str.replace('R1','').astype(float)*R1
    foreman_nodal_harmonic.loc[bool_R2,'factor'] = foreman_nodal_harmonic.loc[bool_R2,'factor'].str.replace('R2','').astype(float)*R2
    foreman_nodal_harmonic = foreman_nodal_harmonic.astype(float)
    foreman_nodal_harmonic.index.name = None
    
    return foreman_doodson_harmonic, foreman_nodal_harmonic


@functools.lru_cache() #caching useful since reading file and internal consistency check once is sufficient
def get_foreman_shallowrelations():
    """
    Omzetten van het derde deel van de foremantabel in een pandas DataFrame met shallow water relations.

    Returns
    -------
    foreman_shallowrelations : TYPE
        DESCRIPTION.

    """
    
    foreman_file = os.path.join(os.path.dirname(file_path),'data','data_foreman_shallowrelations.txt')
    foreman_shallowrelations = pd.read_csv(foreman_file, comment='#', names=range(10), index_col=0, delim_whitespace=True)
    foreman_shallowrelations.index.name = None
    
    #check for internal numdependencies consistency
    for const in foreman_shallowrelations.index:
        foreman_shallow_const = foreman_shallowrelations.loc[const].dropna().values
        num_dependencies = foreman_shallow_const[0]
        list_shallow_facs = foreman_shallow_const[1::2]
        list_shallow_deps = foreman_shallow_const[2::2]
        if not list_shallow_deps.shape == list_shallow_facs.shape == (num_dependencies,):
            raise Exception(f'ERROR: shallow relations not internally consistent:\n{foreman_shallowrelations.loc[const]}')
        
    #check whether all dependencies are available as harmonic components
    foreman_doodson_harmonic, foreman_nodal_harmonic = get_foreman_doodson_nodal_harmonic()
    list_shallowdependencies = foreman_shallowrelations[[3,5,7,9]].melt()['value'].dropna().unique()
    bool_shallowdependencies_isin_harmonics = pd.Series(list_shallowdependencies).isin(foreman_doodson_harmonic.index)
    if not bool_shallowdependencies_isin_harmonics.all():
        raise Exception(f'ERROR: not all required shallow dependency components are available:\n{list_shallowdependencies[~bool_shallowdependencies_isin_harmonics]}')

    return foreman_shallowrelations, list_shallowdependencies


#################################################
#################### FREQ V0 ####################
#################################################


def get_foreman_v0_freq(const_list, dood_date=pd.DatetimeIndex([dt.datetime(1900,1,1)])):
    """
    Zoekt voor iedere component uit de lijst de v op basis van harmonische doodson getallen en de frequentie rechtstreeks uit de foreman tabel.
    Shallow water componenten worden afgeleid met de relaties beschreven in de foreman tabel.
    """
    
    from hatyan.hatyan_core import get_doodson_eqvals, check_requestedconsts # local import since otherwise cross-dependency
    
    check_requestedconsts(tuple(const_list),source='foreman') #TODO: move check to central location when part of hatyan_settings()?

    #get freq/v0 for harmonic components
    foreman_doodson_harmonic, foreman_nodal_harmonic = get_foreman_doodson_nodal_harmonic()

    doodson_pd = get_doodson_eqvals(dood_date=pd.DatetimeIndex([dt.datetime(1900,1,1)]), mode='freq') #TODO: get freq on multiple dates?
    multiply_variables = doodson_pd.loc[['T','S','H','P','P1'],:]
    t_const_freq_dood = np.dot(foreman_doodson_harmonic.loc[:,['T','S','H','P','P1']],multiply_variables) / (2*np.pi)
    foreman_freqs = pd.DataFrame({'freq':t_const_freq_dood[:,0]},index=foreman_doodson_harmonic.index)
    
    doodson_pd = get_doodson_eqvals(dood_date=dood_date, mode=None)
    multiply_variables = doodson_pd.loc[['T','S','H','P','N','P1'],:]
    v_0i_rad = np.dot(foreman_doodson_harmonic.loc[:,['T','S','H','P','N','P1']],multiply_variables) + 2*np.pi*foreman_doodson_harmonic.loc[:,['EDN']].values
    v_0i_rad_harmonic_pd = pd.DataFrame(v_0i_rad,index=foreman_doodson_harmonic.index)
    
    #derive freq/v0 for shallow water components
    foreman_shallowrelations, list_shallowdependencies = get_foreman_shallowrelations()    
    v_0i_rad = pd.DataFrame(np.zeros((len(const_list),len(dood_date))),index=const_list)
    t_const_freq = pd.DataFrame({'freq':np.zeros((len(const_list)))},index=const_list)
    
    #v and freq for harmonic and shallow constituents
    for iC,const in enumerate(const_list):
        if const in v_0i_rad_harmonic_pd.index:
            v_0i_rad.loc[const] = v_0i_rad_harmonic_pd.loc[const]
            t_const_freq.loc[const,'freq'] = foreman_freqs.loc[const,'freq']
        elif const in foreman_shallowrelations.index: #or is not in foreman_harmonic_doodson_all_list
            v_0i_rad_temp = 0
            t_const_freq_temp = 0
            foreman_shallow_const = foreman_shallowrelations.loc[const].dropna().values
            list_shallow_facs = foreman_shallow_const[1::2]
            list_shallow_deps = foreman_shallow_const[2::2]
            for harm_factor,harm_const in zip(list_shallow_facs,list_shallow_deps):
                v_dependency = v_0i_rad_harmonic_pd.loc[harm_const].values
                freq_dependency = foreman_freqs.loc[harm_const,'freq'] #should be dependent on harmonic doodson numbers (make foreman_freqs_dood_all variable in foreman.py, in freq or harmonic definition)
                v_0i_rad_temp += harm_factor*v_dependency
                t_const_freq_temp += harm_factor*freq_dependency
            v_0i_rad.loc[const,:] = v_0i_rad_temp
            t_const_freq.loc[const,'freq'] = t_const_freq_temp
    
    return v_0i_rad, t_const_freq


#################################################
################# NODALFACTORS ##################
#################################################

def get_foreman_nodalfactors(const_list, dood_date):
    """
    Zoekt voor iedere component uit de lijst de u en f (nodal factors) op basis van satellite doodson getallen uit de foreman tabel.
    Shallow water componenten worden afgeleid met de relaties beschreven in de foreman tabel.
    """
    
    from hatyan.hatyan_core import get_doodson_eqvals, check_requestedconsts # local import since otherwise cross-dependency
    
    check_requestedconsts(tuple(const_list),source='foreman') #TODO: move check to central location when part of hatyan_settings()?
    
    doodson_pd = get_doodson_eqvals(dood_date)
    foreman_doodson_harmonic, foreman_nodal_harmonic = get_foreman_doodson_nodal_harmonic()
    foreman_shallowrelations, list_shallowdependencies = get_foreman_shallowrelations()
    
    #extent const_list with missing but required harmonic components
    bool_shallowrequired = pd.Series(foreman_shallowrelations.index).isin(pd.Series(const_list))
    bool_dependencyalready = pd.Series(list_shallowdependencies).isin(pd.Series(const_list))
    if bool_shallowrequired.any() and not bool_dependencyalready.all(): #if any shallow constituent is requested and if not all dependent constituents are in const_list
        const_list_extra = list_shallowdependencies[~bool_dependencyalready]
        const_list_inclshallow = const_list + const_list_extra.tolist() # const_list + not yet available shallowdependencies
    else:
        const_list_inclshallow = const_list
    
    #allocate dataframes including potential exta constituents (are removed before return)
    #TODO: creating empty dataframe with {} and columns/index takes way more time than below, check if this is done somewhere in the code
    nodal_shape = (len(const_list_inclshallow),len(dood_date))
    f_i_FOR_inclshallow = pd.DataFrame(np.ones(nodal_shape), index=const_list_inclshallow)
    u_i_rad_FOR_inclshallow = pd.DataFrame(np.zeros(nodal_shape), index=const_list_inclshallow)
        
    #f and u for harmonic constituents
    for const in const_list_inclshallow:
        if const in foreman_nodal_harmonic.index.unique(): # if harmonic constituent has nodal factors
            foreman_harmonic_nodal_const = foreman_nodal_harmonic.loc[[const]]
            fore_delta_jk_rad_all = np.dot(foreman_harmonic_nodal_const.loc[:,['P','N','P1']],doodson_pd.loc[['P','N','P1'],:])
            fore_alpha_jk_all = foreman_harmonic_nodal_const.loc[:,['EDN']].values * 2*np.pi #phase correction satellite. 0.5=90 voor M2 en N2, 0=0 voor S2
            fore_r_jk_all = foreman_harmonic_nodal_const.loc[:,['factor']].values #amplitude ratio for satellite. 0.0373 voor M2 en N2, 0.0022 voor S2
            fore_fj_left_all = 1 * fore_r_jk_all * np.cos(fore_delta_jk_rad_all + fore_alpha_jk_all) #should be sum for n sattelites
            fore_fj_right_all = 1 * fore_r_jk_all * np.sin(fore_delta_jk_rad_all + fore_alpha_jk_all) #should be sum for n sattelites
            fore_fj_left = fore_fj_left_all.sum(axis=0)
            fore_fj_right = fore_fj_right_all.sum(axis=0)
            f_i_FOR_inclshallow.loc[const,:] = ( (1+fore_fj_left)**2 + (fore_fj_right)**2)**(1/2.)
            u_i_rad_FOR_inclshallow.loc[const,:] = -np.arctan2(fore_fj_right,1+fore_fj_left) #TODO: added minus to get sign comparable to hatyan
        elif const in foreman_doodson_harmonic.index: #if harmonic constituent has no nodal factors, default values 1 and 0 are applicable
            continue
    
    #f and u for shallow constituents
    for const in const_list_inclshallow:
        if const in foreman_shallowrelations.index: # component has satellites based on shallow water relations
            f_i_FOR_temp = 1.0
            u_i_rad_FOR_temp = 0.0
            foreman_shallow_const = foreman_shallowrelations.loc[const].dropna().values
            list_shallow_facs = foreman_shallow_const[1::2]
            list_shallow_deps = foreman_shallow_const[2::2]
            for harm_factor,harm_const in zip(list_shallow_facs,list_shallow_deps):
                f_i_FOR_temp *= f_i_FOR_inclshallow.loc[harm_const,:].values ** abs(harm_factor)
                u_i_rad_FOR_temp += u_i_rad_FOR_inclshallow.loc[harm_const,:].values * harm_factor
            f_i_FOR_inclshallow.loc[const,:] = f_i_FOR_temp
            u_i_rad_FOR_inclshallow.loc[const,:] = u_i_rad_FOR_temp
    
    #drop extra added constituents
    f_i_FOR = f_i_FOR_inclshallow.loc[const_list,:]
    u_i_rad_FOR = u_i_rad_FOR_inclshallow.loc[const_list,:]
    return f_i_FOR, u_i_rad_FOR


