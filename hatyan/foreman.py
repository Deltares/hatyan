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
    foreman_doodson_harmonic = foreman_harmonic_raw.loc[~bool_dupl_index,:8].astype(float)
    foreman_doodson_harmonic.columns = ['T','S','H','P','N','P1','EDN','nsats']
    
    #get all foreman harmonic nodal values (non-first occurences of duplicated component names)
    foreman_nodal_harmonic_wide = foreman_harmonic_raw.loc[bool_dupl_index]
    #reshape from 15 to 5 columns
    foreman_nodal_harmonic_withna = pd.DataFrame(foreman_nodal_harmonic_wide.values.reshape(207,5),
                                                 index=foreman_nodal_harmonic_wide.index.repeat(3),
                                                 columns=['P','N','P1','EDN','factor'])
    #drop nan lines
    foreman_nodal_harmonic = foreman_nodal_harmonic_withna.loc[~foreman_nodal_harmonic_withna.isnull().all(axis=1)]
    
    #multiply with R1/R2 if applicable and convert entire dataframe to floats
    bool_R1 = foreman_nodal_harmonic['factor'].str.contains('R1')
    bool_R2 = foreman_nodal_harmonic['factor'].str.contains('R2')
    foreman_nodal_harmonic.loc[bool_R1,'factor'] = foreman_nodal_harmonic.loc[bool_R1,'factor'].str.replace('R1','').astype(float)*R1
    foreman_nodal_harmonic.loc[bool_R2,'factor'] = foreman_nodal_harmonic.loc[bool_R2,'factor'].str.replace('R2','').astype(float)*R2
    foreman_nodal_harmonic = foreman_nodal_harmonic.astype(float)

    return foreman_doodson_harmonic, foreman_nodal_harmonic


def get_foreman_shallowrelations(pd_series=False):
    """
    Omzetten van het derde deel van de foremantabel in een pandas DataFrame met shallow water relations.

    Returns
    -------
    foreman_shallowrelations : TYPE
        DESCRIPTION.

    """
    
    foreman_file = os.path.join(os.path.dirname(file_path),'data','data_foreman_shallowrelations.txt')
    foreman_shallowrelations_raw = pd.read_csv(foreman_file, comment='#', names=[0], skip_blank_lines=True)[0]
    
    foreman_shallowrelations = foreman_shallowrelations_raw.str.split(' ', expand=True)
    foreman_shallowrelations = foreman_shallowrelations.set_index(0, drop=True)
    foreman_shallowrelations.index.name = None
    
    if not pd_series:
        return foreman_shallowrelations
    
    foreman_shallowrelations_pd = pd.Series()
    for iC,const in enumerate(foreman_shallowrelations.index):
        foreman_shallow_const = foreman_shallowrelations.loc[const].tolist()
        num_dependencies = int(foreman_shallow_const[0])
        shallowrelation_str = ''
        for iD in range(num_dependencies):
            id_factor = iD*2+1
            id_constname = iD*2+2
            harm_factor = float(foreman_shallow_const[id_factor])
            #if harm_factor%1==0:
            #    harm_factor = int(harm_factor)
            harm_const = foreman_shallow_const[id_constname]
            if harm_factor>=0:
                shallowrelation_str += (f' +{harm_factor}*{harm_const}')
            else:
                shallowrelation_str += (f' {harm_factor}*{harm_const}')
            foreman_shallowrelations_pd[const] = shallowrelation_str
    return foreman_shallowrelations_pd
    


#################################################
#################### FREQ V0 ####################
#################################################

@functools.lru_cache() #only caching this already makes foreman slightly faster
def get_foreman_table(): #TODO: only harmonic and only v0
    
    foreman_doodson_harmonic, foreman_nodal_harmonic = get_foreman_doodson_nodal_harmonic()
    t_const_doodson_lun = foreman_doodson_harmonic.copy()
    omega1 = t_const_doodson_lun.loc[:,'T']
    corr_array = pd.DataFrame({'S':-omega1,'H':omega1})
    t_const_doodson_sol = t_const_doodson_lun.add(corr_array,fill_value=0)
    v0_baseT_for = t_const_doodson_sol[['T','S','H','P','N','P1','EDN']]
    
    return v0_baseT_for

    
def get_foreman_v0_freq(const_list, dood_date=pd.DatetimeIndex([dt.datetime(1900,1,1)])):
    """
    Zoekt voor iedere component uit de lijst de v op basis van harmonische doodson getallen en de frequentie rechtstreeks uit de foreman tabel.
    Shallow water componenten worden afgeleid met de relaties beschreven in de foreman tabel.
    """
    
    #get freq/v0 for harmonic components
    from hatyan.hatyan_core import get_doodson_eqvals # local import since otherwise cross-dependency
    
    t_const_doodson_sol = get_foreman_table()
    
    doodson_pd = get_doodson_eqvals(dood_date=pd.DatetimeIndex([dt.datetime(1900,1,1)]), mode='freq') #TODO: get freq on multiple dates?
    multiply_variables = doodson_pd.loc[['T','S','H','P','P1'],:]
    t_const_freq_dood = np.dot(t_const_doodson_sol.loc[:,['T','S','H','P','P1']],multiply_variables) / (2*np.pi)
    foreman_freqs = pd.DataFrame({'freq':t_const_freq_dood[:,0]},index=t_const_doodson_sol.index)
    
    doodson_pd = get_doodson_eqvals(dood_date=dood_date, mode=None)
    multiply_variables = doodson_pd.loc[['T','S','H','P','N','P1'],:]
    v_0i_rad = np.dot(t_const_doodson_sol.loc[:,['T','S','H','P','N','P1']],multiply_variables) + 2*np.pi*t_const_doodson_sol.loc[:,['EDN']].values
    v_0i_rad_harmonic_pd = pd.DataFrame(v_0i_rad,index=t_const_doodson_sol.index)
    
    #derive freq/v0 for shallow water components
    foreman_shallowrelations = get_foreman_shallowrelations()
    #foreman_shallowrelations_pd = get_foreman_shallowrelations(pd_series=True)
    
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
            foreman_shallow_const = foreman_shallowrelations.loc[const].tolist()
            num_dependencies = int(foreman_shallow_const[0])
            for iD in range(num_dependencies):
                id_factor = iD*2+1
                id_constname = iD*2+2
                harm_factor = float(foreman_shallow_const[id_factor])
                harm_const = foreman_shallow_const[id_constname]
                v_dependency = v_0i_rad_harmonic_pd.loc[harm_const].values
                freq_dependency = foreman_freqs.loc[harm_const,'freq'] #should be dependent on harmonic doodson numbers (make foreman_freqs_dood_all variable in foreman.py, in freq or harmonic definition)
                v_0i_rad_temp += harm_factor*v_dependency
                t_const_freq_temp += harm_factor*freq_dependency
            v_0i_rad.loc[const] = v_0i_rad_temp
            t_const_freq.loc[const,'freq'] = t_const_freq_temp
        else:
            raise Exception('ERROR: constituent %s is not in v_0i_rad_harmonic_pd.index and foreman_shallowrelations.index, this is invalid.'%(const))
    
    return v_0i_rad, t_const_freq


#################################################
################# NODALFACTORS ##################
#################################################
def get_foreman_nodalfactors_fromharmonic_oneconst(foreman_harmonic_nodal_const, dood_date):

    from hatyan.hatyan_core import get_doodson_eqvals # local import since otherwise cross-dependency

    doodson_pd = get_doodson_eqvals(dood_date)
    
    fore_delta_jk_rad_all = np.dot(foreman_harmonic_nodal_const.loc[:,['P','N','P1']],doodson_pd.loc[['P','N','P1'],:])
    fore_alpha_jk_all = foreman_harmonic_nodal_const.loc[:,['EDN']].values * 2*np.pi #phase correction satellite. 0.5=90 voor M2 en N2, 0=0 voor S2
    fore_r_jk_all = foreman_harmonic_nodal_const.loc[:,['factor']].values #amplitude ratio for satellite. 0.0373 voor M2 en N2, 0.0022 voor S2
    fore_fj_left_all = 1 * fore_r_jk_all * np.cos(fore_delta_jk_rad_all + fore_alpha_jk_all) #should be sum for n sattelites
    fore_fj_right_all = 1 * fore_r_jk_all * np.sin(fore_delta_jk_rad_all + fore_alpha_jk_all) #should be sum for n sattelites
    fore_fj_left = fore_fj_left_all.sum(axis=0)
    fore_fj_right = fore_fj_right_all.sum(axis=0)
    
    f_i_FOR = ( (1+fore_fj_left)**2 + (fore_fj_right)**2)**(1/2.)
    u_i_rad_FOR = -np.arctan2(fore_fj_right,1+fore_fj_left) #TODO: added minus to get sign comparable to hatyan

    return f_i_FOR, u_i_rad_FOR


def get_foreman_nodalfactors(const_list, dood_date):
    """
    Zoekt voor iedere component uit de lijst de u en f (nodal factors) op basis van satellite doodson getallen uit de foreman tabel.
    Shallow water componenten worden afgeleid met de relaties beschreven in de foreman tabel.
    """
    
    foreman_shallowrelations = get_foreman_shallowrelations()
    #foreman_shallowrelations_pd = get_foreman_shallowrelations(pd_series=True)
    foreman_doodson_harmonic, foreman_nodal_harmonic = get_foreman_doodson_nodal_harmonic()
       
    f_i_FOR = pd.DataFrame(np.ones((len(const_list),len(dood_date))), index=const_list)
    u_i_rad_FOR = pd.DataFrame(np.zeros((len(const_list),len(dood_date))), index=const_list)
    
    #f and u for harmonic constituents
    for iC,const in enumerate(const_list):
        if const in foreman_doodson_harmonic.index:
            if const not in foreman_nodal_harmonic.index.unique(): # if harmonic constituent has no nodal factors
                continue
            foreman_harmonic_nodal_const = foreman_nodal_harmonic.loc[[const]]
            f_i_FOR.loc[const,:], u_i_rad_FOR.loc[const,:] = get_foreman_nodalfactors_fromharmonic_oneconst(foreman_harmonic_nodal_const, dood_date)
        elif const in foreman_shallowrelations.index: # component has satellites based on shallow water relations
            f_i_FOR_temp = 1.0
            u_i_rad_FOR_temp = 0.0
            #temp_nodal_df = pd.DataFrame()
            foreman_shallow_const = foreman_shallowrelations.loc[const].tolist()
            #foreman_shallow_const_pd = foreman_shallowrelations_pd.loc[const]
            num_dependencies = int(foreman_shallow_const[0])
            for iD in range(num_dependencies):
                id_factor = iD*2+1
                id_constname = iD*2+2
                harm_factor = float(foreman_shallow_const[id_factor])
                harm_const = foreman_shallow_const[id_constname]
                if harm_const not in foreman_nodal_harmonic.index.unique():
                    raise Exception('ERROR: harmonic component %s for shallow water component %s is not available in the harmonic nodal factors (foreman_nodal_harmonic)'%(harm_const,const))
                foreman_harmonic_nodal_const = foreman_nodal_harmonic.loc[[harm_const]]
                #temp_nodal_df_onestep = pd.concat([harm_factor*foreman_nodal_all.loc[harm_const,:3],foreman_nodal_all.loc[harm_const,[4]]],axis=1)
                #temp_nodal_df = temp_nodal_df.append(temp_nodal_df_onestep)
                f_i_dependency, u_i_rad_dependency = get_foreman_nodalfactors_fromharmonic_oneconst(foreman_harmonic_nodal_const, dood_date)#foreman_harmonic_nodal_all[foreman_harmonic_nodal_all_list.index()][iS]
                f_i_FOR_temp *= f_i_dependency**abs(harm_factor)
                u_i_rad_FOR_temp += harm_factor*u_i_rad_dependency
            f_i_FOR.loc[const,:] = f_i_FOR_temp
            u_i_rad_FOR.loc[const,:] = u_i_rad_FOR_temp
            #temp_nodal_df.index = ['M4']*len(temp_nodal_df.index)
            #f_i_FOR2[iC,:], u_i_rad_FOR2[iC,:] = get_foreman_nodalfactors_fromharmonic_oneconst(temp_nodal_df, dood_date)
        else:
            raise Exception('ERROR: constituent %s is not in foreman_doodson_harmonic.index and foreman_shallowrelations.index, this is invalid.'%(const))
    
    #print(f_i_FOR-f_i_FOR2)
    #print(u_i_rad_FOR-u_i_rad_FOR2)
    return f_i_FOR, u_i_rad_FOR


