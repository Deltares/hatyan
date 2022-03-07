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
    foreman_harmonic_raw = pd.read_csv(foreman_file, comment='#', names=[0], skip_blank_lines=True)[0]
    
    lat_rad = np.deg2rad(lat_deg)
    
    R1 = 0.36309*(1.-5.*np.sin(lat_rad)*np.sin(lat_rad))/np.sin(lat_rad) #-1 #for lat=50N
    R2 = 2.59808*np.sin(lat_rad) #2 #for lat=50N

    foreman_harmonic_new = foreman_harmonic_raw.str.split(' ', expand=True)
    foreman_harmonic_new = foreman_harmonic_new.set_index(0, drop=True)
    foreman_harmonic_new.index.name = None
    foreman_harmonic_list_new = list(set(foreman_harmonic_new.index.tolist()))

    #get all foreman harmonic doodson and nodal values
    foreman_harmonic_doodson_all_new = pd.DataFrame()
    foreman_harmonic_nodal_all_new = pd.DataFrame()
    for iC, const in enumerate(foreman_harmonic_list_new):
        foreman_harmonic_sel = foreman_harmonic_new.loc[[const]]
        
        row1 = foreman_harmonic_sel.iloc[0,:]
        row1 = row1[~row1.isnull()] #remove nans
        if not len(row1) == 8:
            raise Exception('ERROR: first row of component %s does not have length 8, so corrupt foreman file, line:\n%s'%(const,row1))
        foreman_harmonic_doodson_all_new[const] = row1
        
        if len(foreman_harmonic_sel.index) > 1:
            row_other_raw = foreman_harmonic_sel.iloc[1:,:]
            row_other = row_other_raw.stack(dropna=True)
            row_otherlen = len(row_other)
            if (row_otherlen%5)!=0:
                raise Exception('ERROR: length of list of nodal values should be a multiple of 5, corrupt foreman file for %s'%(const))
            row_other_reshape = row_other.values.reshape(row_otherlen//5,5)
            row_other_reshape_pd = pd.DataFrame(row_other_reshape,index=[const]*(row_otherlen//5))
            foreman_harmonic_nodal_all_new = foreman_harmonic_nodal_all_new.append(row_other_reshape_pd)
    foreman_doodson_harmonic = foreman_harmonic_doodson_all_new.astype(float).T
    
    bool_R1 = foreman_harmonic_nodal_all_new[4].str.contains('R1')
    bool_R2 = foreman_harmonic_nodal_all_new[4].str.contains('R2')
    foreman_harmonic_nodal_all_new.loc[bool_R1,4] = foreman_harmonic_nodal_all_new.loc[bool_R1,4].str.replace('R1','').astype(float)*R1
    foreman_harmonic_nodal_all_new.loc[bool_R2,4] = foreman_harmonic_nodal_all_new.loc[bool_R2,4].str.replace('R2','').astype(float)*R2
    foreman_nodal_harmonic = foreman_harmonic_nodal_all_new.astype(float)

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
    const_list = foreman_doodson_harmonic.index
    t_const_doodson_lun = np.concatenate([np.zeros((len(const_list),1)),foreman_doodson_harmonic.loc[:,:7].values],axis=1)
    omega1 = t_const_doodson_lun[:,1:2]
    corr_array = np.concatenate([omega1,np.zeros((len(const_list),1)),-omega1,omega1,np.zeros((len(const_list),4))],axis=1)
    t_const_doodson_sol = t_const_doodson_lun+corr_array
    t_const_doodson_sol[:,1] = 0
    
    index_v0 = ['T','S','H','P','N','P1','EDN']
    index_v0_incldummy = ['T','dummy','S','H','P','N','P1','EDN']
    v0_baseT_for_raw = pd.DataFrame(t_const_doodson_sol,index=const_list,columns=index_v0_incldummy)
    v0_baseT_for = v0_baseT_for_raw.loc[:,index_v0]
    
    return v0_baseT_for


def get_foreman_v0freq_fromfromharmonicdood(dood_date=None, mode=None):
    """
    Zoekt de frequentie of v0 voor alle harmonische componenten, in geval van v0 op de gegeven datum (dood_date). Hiervoor zijn de harmonische doodson getallen
    (foreman_harmonic_doodson_all) nodig, afkomstig uit get_foreman_harmonic uit foreman.py 
    """
    
    import numpy as np
    import pandas as pd
    import datetime as dt
    
    from hatyan.schureman import get_doodson_eqvals
    
    if dood_date is None: #in case of frequency
        dood_date = pd.DatetimeIndex([dt.datetime(1900,1,1)]) #dummy value
    
    t_const_doodson_sol = get_foreman_table()
    const_list = t_const_doodson_sol.index

    dood_T_rad, dood_S_rad, dood_H_rad, dood_P_rad, dood_N_rad, dood_P1_rad = get_doodson_eqvals(dood_date=dood_date, mode=mode)

    if mode=='freq':
        dood_rad_array = np.stack([dood_T_rad,dood_S_rad,dood_H_rad,dood_P_rad,np.zeros((dood_N_rad.shape)),dood_P1_rad,np.zeros((dood_T_rad.shape))])
        t_const_freq_dood = np.dot(t_const_doodson_sol,dood_rad_array)/(2*np.pi)
        t_const_freqv0_dood_pd = pd.DataFrame({'freq':t_const_freq_dood[:,0]},index=const_list)
        freqv0_dood_pd = t_const_freqv0_dood_pd
    else:
        dood_rad_array = np.stack([dood_T_rad,dood_S_rad,dood_H_rad,dood_P_rad,dood_N_rad,dood_P1_rad,np.zeros((dood_T_rad.shape))+2*np.pi])
        v_0i_rad = np.dot(t_const_doodson_sol,dood_rad_array)
        v_0i_rad_pd = pd.DataFrame(v_0i_rad,index=const_list)
        freqv0_dood_pd = v_0i_rad_pd
    
    return freqv0_dood_pd


def get_foreman_v0_freq(const_list, dood_date):
    """
    Zoekt voor iedere component uit de lijst de v op basis van harmonische doodson getallen en de frequentie rechtstreeks uit de foreman tabel.
    Shallow water componenten worden afgeleid met de relaties beschreven in de foreman tabel.
    """
    import numpy as np
    import pandas as pd
    
    from hatyan.schureman import get_const_list_hatyan
    
    if type(const_list) is str:
        const_list = get_const_list_hatyan(const_list)
    elif type(const_list) is not list:
        const_list = const_list.tolist()

    foreman_freqs = get_foreman_v0freq_fromfromharmonicdood(dood_date=None, mode='freq') #list with only harmonic components with more precision than file
    v_0i_rad_harmonic_pd = get_foreman_v0freq_fromfromharmonicdood(dood_date=dood_date, mode=None)
    
    foreman_shallowrelations = get_foreman_shallowrelations()
    foreman_harmonic_list = v_0i_rad_harmonic_pd.index.tolist()
    foreman_shallowrelations_list = foreman_shallowrelations.index.tolist()
    
    v_0i_rad = np.zeros((len(const_list),len(dood_date)))
    t_const_freq = np.zeros((len(const_list)))
    
    #v and freq for harmonic and shallow constituents
    for iC,const in enumerate(const_list):
        if const in foreman_harmonic_list:
            v_0i_rad[iC,:] = v_0i_rad_harmonic_pd.loc[const]
            t_const_freq[iC] = foreman_freqs.loc[const,'freq']
        elif const in foreman_shallowrelations_list: #or is not in foreman_harmonic_doodson_all_list
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
                v_0i_rad_temp = v_0i_rad_temp + harm_factor*v_dependency
                t_const_freq_temp = t_const_freq_temp + harm_factor*freq_dependency
            v_0i_rad[iC,:] = v_0i_rad_temp
            t_const_freq[iC] = t_const_freq_temp
        else:
            raise Exception('ERROR: constituent %s is not in foreman_harmonic_doodson_all_list and foreman_shallow_all_list, this is invalid.'%(const))
    v_0i_rad_pd = pd.DataFrame(v_0i_rad,index=const_list)
    t_const_freq_pd = pd.DataFrame({'freq':t_const_freq},index=const_list)
    
    return v_0i_rad_pd, t_const_freq_pd


#################################################
################# NODALFACTORS ##################
#################################################
def get_foreman_nodalfactors_fromharmonic_oneconst(foreman_harmonic_nodal_const, dood_date):
    import numpy as np
        
    from hatyan.hatyan_core import get_doodson_eqvals
    dood_T_rad, dood_S_rad, dood_H_rad, dood_P_rad, dood_N_rad, dood_P1_rad = get_doodson_eqvals(dood_date)
    
    fore_delta_jk_rad_all = np.dot(foreman_harmonic_nodal_const.loc[:,0:2],np.stack([dood_P_rad, dood_N_rad, dood_P1_rad]))
    fore_alpha_jk_all = foreman_harmonic_nodal_const.loc[:,3:3].values * 2*np.pi #phase correction satellite. 0.5=90 voor M2 en N2, 0=0 voor S2
    fore_r_jk_all = foreman_harmonic_nodal_const.loc[:,4:4].values #amplitude ratio for satellite. 0.0373 voor M2 en N2, 0.0022 voor S2
    fore_A_jk = 1
    fore_fj_left_all = fore_A_jk * fore_r_jk_all * np.cos(fore_delta_jk_rad_all + fore_alpha_jk_all) #should be sum for n sattelites
    fore_fj_right_all = fore_A_jk * fore_r_jk_all * np.sin(fore_delta_jk_rad_all + fore_alpha_jk_all) #should be sum for n sattelites
    fore_fj_left = fore_fj_left_all.sum(axis=0)
    fore_fj_right = fore_fj_right_all.sum(axis=0)
    
    f_i_FOR = ( (1+ fore_fj_left)**2 + (fore_fj_right)**2)**(1/2.)
    u_i_rad_FOR = -np.arctan2(fore_fj_right,1+fore_fj_left) #added minus to get sign comparable to hatyan

    return f_i_FOR, u_i_rad_FOR


def get_foreman_nodalfactors(const_list, dood_date):
    """
    Zoekt voor iedere component uit de lijst de u en f (nodal factors) op basis van satellite doodson getallen uit de foreman tabel.
    Shallow water componenten worden afgeleid met de relaties beschreven in de foreman tabel.
    """
    import numpy as np
    import pandas as pd
    from hatyan.hatyan_core import get_const_list_hatyan
    
    if type(const_list) is str:
        const_list = get_const_list_hatyan(const_list)
    elif type(const_list) is not list:
        const_list = const_list.tolist()

    foreman_shallowrelations = get_foreman_shallowrelations()
    foreman_doodson_harmonic, foreman_nodal_harmonic = get_foreman_doodson_nodal_harmonic()
    #foreman_nodal_all = foreman_nodal_harmonic.copy()
    
    foreman_harmonic_doodson_all_list = foreman_doodson_harmonic.index.tolist()
    foreman_harmonic_nodal_all_list = foreman_nodal_harmonic.index.unique().tolist()
    foreman_shallowrelations_list = foreman_shallowrelations.index.tolist()
    
    f_i_FOR = np.ones((len(const_list),len(dood_date)))
    u_i_rad_FOR = np.zeros((len(const_list),len(dood_date)))
    #f_i_FOR2 = np.ones((len(const_list),len(dood_date)))
    #u_i_rad_FOR2 = np.zeros((len(const_list),len(dood_date)))
    
    #f and u for harmonic constituents
    for iC,const in enumerate(const_list):
        if const in foreman_harmonic_doodson_all_list:
            if const in foreman_harmonic_nodal_all_list:
                foreman_harmonic_nodal_const = foreman_nodal_harmonic.loc[[const]]
                f_i_FOR[iC,:], u_i_rad_FOR[iC,:] = get_foreman_nodalfactors_fromharmonic_oneconst(foreman_harmonic_nodal_const, dood_date)
                #f_i_FOR2[iC,:], u_i_rad_FOR2[iC,:] = get_foreman_nodalfactors_fromharmonic_oneconst(foreman_harmonic_nodal_const, dood_date)
        elif const in foreman_shallowrelations_list: # component has satellites based on shallow water relations
            f_i_FOR_temp = 1.0
            u_i_rad_FOR_temp = 0.0
            #temp_nodal_df = pd.DataFrame()
            foreman_shallow_const = foreman_shallowrelations.loc[const].tolist()
            num_dependencies = int(foreman_shallow_const[0])
            for iD in range(num_dependencies):
                id_factor = iD*2+1
                id_constname = iD*2+2
                harm_factor = float(foreman_shallow_const[id_factor])
                harm_const = foreman_shallow_const[id_constname]
                if harm_const not in foreman_harmonic_nodal_all_list:
                    raise Exception('ERROR: harmonic component %s for shallow water component %s is not available in the harmonic nodal factors (foreman_nodal_harmonic)'%(harm_const,const))
                foreman_harmonic_nodal_const = foreman_nodal_harmonic.loc[[harm_const]]
                #temp_nodal_df_onestep = pd.concat([harm_factor*foreman_nodal_all.loc[harm_const,:3],foreman_nodal_all.loc[harm_const,[4]]],axis=1)
                #temp_nodal_df = temp_nodal_df.append(temp_nodal_df_onestep)
                f_i_dependency, u_i_rad_dependency = get_foreman_nodalfactors_fromharmonic_oneconst(foreman_harmonic_nodal_const, dood_date)#foreman_harmonic_nodal_all[foreman_harmonic_nodal_all_list.index()][iS]
                f_i_FOR_temp = f_i_FOR_temp * f_i_dependency**abs(harm_factor)
                u_i_rad_FOR_temp = u_i_rad_FOR_temp + harm_factor*u_i_rad_dependency
            f_i_FOR[iC,:] = f_i_FOR_temp
            u_i_rad_FOR[iC,:] = u_i_rad_FOR_temp
            #temp_nodal_df.index = ['M4']*len(temp_nodal_df.index)
            #f_i_FOR2[iC,:], u_i_rad_FOR2[iC,:] = get_foreman_nodalfactors_fromharmonic_oneconst(temp_nodal_df, dood_date)
        else:
            raise Exception('ERROR: constituent %s is not in foreman_harmonic_doodson_all_list and foreman_shallow_all_list, this is invalid.'%(const))
    
    f_i_FOR_pd = pd.DataFrame(f_i_FOR, index=const_list)
    u_i_rad_FOR_pd = pd.DataFrame(u_i_rad_FOR, index=const_list)
    #f_i_FOR2_pd = pd.DataFrame(f_i_FOR2, index=const_list)
    #u_i_rad_FOR2_pd = pd.DataFrame(u_i_rad_FOR2, index=const_list)
    #print(f_i_FOR_pd-f_i_FOR2_pd)
    #print(u_i_rad_FOR_pd-u_i_rad_FOR2_pd)
    return f_i_FOR_pd, u_i_rad_FOR_pd


