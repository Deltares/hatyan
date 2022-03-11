# - * - coding: utf - 8 - * - 
"""
schureman.py contains definitions with data for the schureman constituents.

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
import functools #to Memoize v0uf table (https://en.wikipedia.org/wiki/Memoization)
import numpy as np
import datetime as dt

file_path = os.path.realpath(__file__)


@functools.lru_cache()
def get_schureman_table():
    """
    Calculate all schureman constituents

    Returns
    -------
    v0uf_all : TYPE
        DESCRIPTION.
    
    """
    
    index_v0 = ['T','S','H','P','N','P1','EDN']
    index_u = ['DKSI','DNU','DQ','DQU','DR','DUK1','DUK2']
    index_f = ['DND73','DND74','DND75','DND76','DND77','DND78','DND79','DFM1','DFK1','DFL2','DFK2','DFM1C']
    #index_fstr =['f_eqs']

    file_schureman_harmonic = os.path.join(os.path.dirname(file_path),'data','data_schureman_harmonic.csv')
    v0uf_baseT = pd.read_csv(file_schureman_harmonic,comment='#',skipinitialspace=True,index_col='component')
    v0uf_base = v0uf_baseT.T
    #v0uf_base.index = index_v0 + index_u + index_f + index_fstr
    
    file_schureman_shallowrelations = os.path.join(os.path.dirname(file_path),'data','data_schureman_shallowrelations.csv')
    shallow_eqs_pd = pd.read_csv(file_schureman_shallowrelations,comment='#',skipinitialspace=True,index_col=0,names=['shallow_eq'])
    shallow_eqs_pd['shallow_eq'] = shallow_eqs_pd['shallow_eq'].str.strip() #remove spaces after
    shallow_eqs_pd['shallow_const'] = shallow_eqs_pd.index
    
    shallow_eqs_pd.index = 'comp_'+shallow_eqs_pd.index.str.replace('(','',regex=False).str.replace(')','',regex=False)#brackets are temporarily removed in order to evaluate functions
    shallow_eqs_pd_str = '\n'.join(f'{key} = {val}' for key, val in shallow_eqs_pd['shallow_eq'].iteritems()) 
    
    #calculate shallow water components and rename back to original component name
    v0uf_base_forv0u = v0uf_base.loc[index_v0+index_u,:].astype(int)
    v0uf_base_forv0u.eval(shallow_eqs_pd_str, inplace=True)
    v0uf_base_forf = v0uf_base.loc[index_f,:].astype(float)
    v0uf_base_forf.eval(shallow_eqs_pd_str.replace('-','+'), inplace=True) #for f only multiplication is applied, never division
    v0uf_all = pd.concat([v0uf_base_forv0u,v0uf_base_forf])
    v0uf_all.rename(columns=shallow_eqs_pd['shallow_const'],inplace=True)
    v0uf_allT = v0uf_all.T
    
    #from hatyan.hatyan_core import get_lunarSLSIHO_fromsolar # local import since otherwise cross-dependency
    #v0uf_allT_lunar, v0uf_allT_lunar_SLS, v0uf_allT_lunar_IHO = get_lunarSLSIHO_fromsolar(v0uf_all)
    
    return v0uf_allT


def get_schureman_freqs(const_list, dood_date=pd.DatetimeIndex([dt.datetime(1900,1,1)]), return_allraw=False):
    """
    Returns the frequencies of the requested list of constituents. Source: beromg.f

    Parameters
    ----------
    const_list : list or pandas.Series
        contains the tidal constituent names.

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    t_const_freq : TYPE
        DESCRIPTION.

    """
    
    from hatyan.hatyan_core import get_doodson_eqvals, check_requestedconsts # local import since otherwise cross-dependency
    
    check_requestedconsts(tuple(const_list),source='schureman') #TODO: move check to central location when part of hatyan_settings()?
    
    doodson_pd = get_doodson_eqvals(dood_date=dood_date, mode='freq') #N is not used here
    multiply_variables = doodson_pd.loc[['T','S','H','P','P1'],:]
    
    v0uf_allT = get_schureman_table()
    v0uf_sel = v0uf_allT.loc[const_list]
    v0uf_sel_freq = v0uf_sel[['T','S','H','P','P1']]
    
    DOMEGA_speed = np.dot(v0uf_sel_freq,multiply_variables)
    if return_allraw: #return array of speeds for each component/timestep
        return DOMEGA_speed
    
    #return dataframe of freq/angvelo/period/angfreq for all components and first timestep
    t_const_speed = DOMEGA_speed[:,0]
    t_const_freq = t_const_speed/(2*np.pi) #aantal rotaties per uur, freq
    np.seterr(divide='ignore') #suppress divide by 0 warning
    t_const_perds = 1/t_const_freq #period [hr]
    t_const_angfreqs = 360/t_const_perds #angfreq [deg/hr]
    freq_pd = pd.DataFrame({'freq':t_const_freq, 'angvelo [rad/hr]':t_const_speed, 'period [hr]':t_const_perds, 'angfreq [deg/hr]':t_const_angfreqs}, index=const_list)
    return freq_pd


def get_schureman_v0(const_list, dood_date):
    """
    Returns the v-values of the requested list of constituents for the requested date(s)
    
    Parameters
    ----------
    const_list : list or pandas.Series
        contains the tidal constituent names.
    dood_date : TYPE
        DESCRIPTION.

    Returns
    -------
    v_0i_rad : TYPE
        DESCRIPTION.

    """
    
    from hatyan.hatyan_core import get_doodson_eqvals, check_requestedconsts # local import since otherwise cross-dependency
    
    check_requestedconsts(tuple(const_list),source='schureman') #TODO: move check to central location when part of hatyan_settings()?

    doodson_pd = get_doodson_eqvals(dood_date=dood_date) #N is not used here
    multiply_variables = doodson_pd.loc[['T','S','H','P','P1'],:]
    
    v0uf_allT = get_schureman_table()
    v0uf_sel = v0uf_allT.loc[const_list]
    DV0 = np.dot(v0uf_sel.loc[:,['T','S','H','P','P1']],multiply_variables) + np.deg2rad(v0uf_sel.loc[:,['EDN']]).values
    DV0_pd = pd.DataFrame(DV0,index=const_list)
    
    return DV0_pd


def get_schureman_constants(dood_date):
    """
    get_schureman_constants

    Parameters
    ----------
    dood_date : TYPE
        DESCRIPTION.

    Returns
    -------
    DOMEGA : TYPE
        DESCRIPTION.
    DIKL : TYPE
        DESCRIPTION.
    DC1681 : TYPE
        DESCRIPTION.
    DC5023 : TYPE
        DESCRIPTION.
    DC0365 : TYPE
        DESCRIPTION.

    """
    
    #TODO: robust_timedelta_sec might also be necesary in other definitions, but is not there yet. Align? (also check if newer pandas version do not have this problem anymore)
    from hatyan.hatyan_core import robust_timedelta_sec # local import since otherwise cross-dependency
    
    #bercon.f: HET BEREKENEN VAN DE 'CONSTANTEN' .0365, .1681 EN .5023, DIE GEBRUIKT WORDEN BIJ DE BEREKENING VAN DE U- EN F-FACTOREN VAN DE GETIJCOMPONENTEN K1 EN K2
    #327932: ratio of mass of sun to combined mass of earth and moon
    #81.53: ratio of mass of earth to mass of moon
    DAGC = 0.01657                            #A/C >> ZIE BLZ. 162 VAN SCHUREMAN, mean lunar parallax in respect to mean radius [rad]
    DAGC1 = 0.00004261                        #A/C1 >> ZIE BLZ. 162 VAN SCHUREMAN, mean solar parallax in respect to mean radius [rad]
    DE = 0.054900489                          #E >> ZIE BLZ. 162 VAN SCHUREMAN, eccentricity of moons orbit
    DMGE = 1./81.53                           #M/E >> ZIE BLZ. 162 VAN SCHUREMAN, mass of moon /mass of earth
    DSGE = 82.53/81.53*327932                 #S/E >> ZIE BLZ. 162 VAN SCHUREMAN, mass of sun / mass of earth
    
    DU = DMGE*DAGC**3                         #U >> ZIE BLZ. 162 VAN SCHUREMAN, basic factor
    DU1 = DSGE*DAGC1**3                       #U1 >> ZIE BLZ. 162 VAN SCHUREMAN, solar coefficient
    DSACCE = DU1/DU                           #S' >> ZIE BLZ. 162 VAN SCHUREMAN, solar factor
    DOMEGA = np.deg2rad(23+27./60.+8.26/3600) #OMEGA >> ZIE BLZ. 162 VAN SCHUREMAN, Obliquity of the Ecliptic, epoch 1 Jan 1900, 23.45229 graden [rad]
    DIKL = np.deg2rad(5+8./60.+43.3546/3600)  #i >> ZIE BLZ. 162 VAN SCHUREMAN, Inclination of moons orbit to plane of ecliptic, epoch 1 Jan 1900 (?), 5.1453762 graden [rad]
    
    DC5023 = 0.5+0.75*DE*DE

    #time dependent
    dood_tstart_sec, fancy_pddt = robust_timedelta_sec(dood_date)    
    dood_Tj = (dood_tstart_sec/3600+12)/(24*36525)
    DE1 = 0.01675104-0.0000418*dood_Tj        #E1 >> ZIE BLZ. 162 VAN SCHUREMAN, eccentricity of earths orbit, epoch 1 Jan 1900
    DCOFSI = (0.5+0.75*DE1**2)*DSACCE         #COEFFICIENT VAN DE SINUSTERMEN IN (217) EN (219) OP BLZ. 45 VAN SCHUREMAN

    DC0365 = DCOFSI*np.sin(DOMEGA)**2
    DC1681 = DCOFSI*np.sin(2*DOMEGA)
    
    return DOMEGA, DIKL, DC1681, DC5023, DC0365


def get_schureman_u(const_list, dood_date):
    """
    Returns the u-values of the requested list of constituents for the requested date(s)

    Parameters
    ----------
    const_list : list or pandas.Series
        contains the tidal constituent names.
    dood_date : TYPE
        DESCRIPTION.

    Returns
    -------
    u_i_rad_HAT : TYPE
        DESCRIPTION.

    """
    
    from hatyan.hatyan_core import get_doodson_eqvals, check_requestedconsts # local import since otherwise cross-dependency
    
    check_requestedconsts(tuple(const_list),source='schureman') #TODO: move check to central location when part of hatyan_settings()?

    doodson_pd = get_doodson_eqvals(dood_date=dood_date)
    N_rad = doodson_pd.loc['N',:].values
    P_rad = doodson_pd.loc['P',:].values
    
    #list of dependencies for U (only P and N are used here)
    DOMEGA, DIKL, DC1681, DC5023, DC0365 = get_schureman_constants(dood_date)
    DHOMI = (DOMEGA-DIKL)*0.5
    DHOPI = (DOMEGA+DIKL)*0.5
    
    DTHN = np.tan(N_rad*0.5)
    DATC = np.arctan2(np.cos(DHOMI)*DTHN,np.cos(DHOPI))
    DATS = np.arctan2(np.sin(DHOMI)*DTHN,np.sin(DHOPI))
    DIH = np.arccos(np.cos(DIKL)*np.cos(DOMEGA)-np.sin(DIKL)*np.sin(DOMEGA)*np.cos(N_rad))
    DKSI = (N_rad-DATC-DATS)
    DNU = DATC-DATS
    DPMKSI = P_rad-DKSI
    DC2PMK = np.cos(DPMKSI+DPMKSI)
    DS2PMK = np.sin(DPMKSI+DPMKSI)
    DCIH = np.cos(DIH)
    DSIH = np.sin(DIH)
    DCHIH = np.cos(DIH*0.5)
    DCTHIH = 1./np.tan(DIH*0.5)
    DS2IH = np.sin(DIH+DIH)
    DQU = np.arctan2(DS2PMK,(3*DCIH/(DCHIH*DCHIH)+DC2PMK))
    DQ = np.arctan2((5*DCIH-1)*np.tan(DPMKSI),(7*DCIH+1))
    DNUACC = np.arctan2(DS2IH*np.sin(DNU),(DS2IH*np.cos(DNU)+DC1681/DC5023))
    DR = np.arctan2(DS2PMK,(DCTHIH*DCTHIH/6.-DC2PMK))
    D2NU2A = np.arctan2(DSIH**2*np.sin(2*DNU),(DSIH**2*np.cos(2*DNU)+DC0365/DC5023))
    DUK1   = -DNUACC 
    DUK2   = -D2NU2A
    multiply_variables = np.stack([DKSI, DNU, DQ, DQU, DR, DUK1, DUK2])
    
    v0uf_allT = get_schureman_table()
    v0uf_sel = v0uf_allT.loc[const_list]
    v0uf_sel_u = v0uf_sel[['DKSI','DNU', 'DQ', 'DQU', 'DR', 'DUK1', 'DUK2']]

    DU = np.dot(v0uf_sel_u.values,multiply_variables)
    DU = np.remainder(DU+np.pi, 2*np.pi) - np.pi

    #create dataframe
    DU_pd = pd.DataFrame(DU, index=const_list)
    
    return DU_pd


def get_schureman_f(const_list, dood_date, xfac):
    """
    Returns the f-values of the requested list of constituents for the requested date(s)

    Parameters
    ----------
    const_list : list or pandas.Series
        contains the tidal constituent names.
    dood_date : TYPE
        DESCRIPTION.

    Returns
    -------
    f_i_HAT : TYPE
        DESCRIPTION.

    """
    
    from hatyan.hatyan_core import get_doodson_eqvals, check_requestedconsts # local import since otherwise cross-dependency
    
    check_requestedconsts(tuple(const_list),source='schureman') #TODO: move check to central location when part of hatyan_settings()?

    doodson_pd = get_doodson_eqvals(dood_date=dood_date)
    N_rad = doodson_pd.loc['N',:].values
    P_rad = doodson_pd.loc['P',:].values
    
    #list of dependencies F (only P and N are used here)
    DOMEGA, DIKL, DC1681, DC5023, DC0365 = get_schureman_constants(dood_date)
    DHOMI = (DOMEGA-DIKL)*0.5
    DHOPI = (DOMEGA+DIKL)*0.5
    DSOMEG = np.sin(DOMEGA)
    DSIKL = np.sin(DIKL)
    DTHN = np.tan(N_rad*0.5)
    DATC = np.arctan2(np.cos(DHOMI)*DTHN,np.cos(DHOPI))
    DATS = np.arctan2(np.sin(DHOMI)*DTHN,np.sin(DHOPI))
    DIH = np.arccos(np.cos(DIKL)*np.cos(DOMEGA)-DSIKL*DSOMEG*np.cos(N_rad))
    DKSI = (N_rad-DATC-DATS)%(2*np.pi)
    DNU = DATC-DATS#
    DCHOM = np.cos(DOMEGA*0.5) 
    DSHOM = np.sin(DOMEGA*0.5) 
    DSOM2 = DSOMEG*DSOMEG
    DCHOM2 = DCHOM*DCHOM 
    DSHOM2 = DSHOM*DSHOM 
    DCHIKL = np.cos(DIKL*0.5)  
    DCHIK4 = DCHIKL*DCHIKL*DCHIKL*DCHIKL
    D132S2 = 1-3./2.*DSIKL*DSIKL 
    DMOF65 = (2./3.-DSOM2)*D132S2  
    DMOF66 = DSOM2*DCHIK4
    DMOF67 = DSOMEG*DCHOM2*DCHIK4
    DMOF68 = np.sin(DOMEGA+DOMEGA)*D132S2
    DMOF69 = DSOMEG*DSHOM2*DCHIK4
    DMOF70 = DCHOM2*DCHOM2*DCHIK4
    DMOF71 = DSOM2*D132S2
    DSIH = np.sin(DIH)#
    DCHIH = np.cos(DIH*0.5)
    DSHIH = np.sin(DIH*0.5)
    DS2IH = np.sin(DIH+DIH)
    DTHIH2 = DSHIH**2/DCHIH**2
    DIGHI2 = np.cos(DIH)/DCHIH**2  
    DPMKSI = P_rad-DKSI
    DC2PMK = np.cos(DPMKSI+DPMKSI)
    DQA = 1./np.sqrt(0.25+1.5*DIGHI2*DC2PMK+2.25*DIGHI2*DIGHI2)
    DRA = 1./np.sqrt(1-12*DTHIH2*DC2PMK+36*DTHIH2*DTHIH2)
    
    #actual values used in the component list below (12 uniques)
    DND73 = (2./3.-DSIH**2)/DMOF65
    DND74 = DSIH**2/DMOF66 
    DND75 = DSIH*DCHIH**2/DMOF67 #DFQ1/DFO1
    DND76 = DS2IH/DMOF68 #DFJ1
    DND77 = DSIH*DSHIH**2/DMOF69 
    DND78 = DCHIH**2*DCHIH**2/DMOF70 #DFN2/DFNU2/DFM2/DFLAB2
    DND79 = DSIH**2/DMOF71 
    DFM1 = DND75/DQA #DFO1/DQA
    DFK1 = np.sqrt(DC5023*DC5023*DS2IH*DS2IH+2*DC5023*DC1681*DS2IH*np.cos(DNU)+DC1681*DC1681)/(DC5023*DMOF68+DC1681) 
    DFL2 = DND78/DRA #DFM2/DRA
    DFK2 = np.sqrt(DC5023*DC5023*DSIH**2*DSIH**2+2*DC5023*DC0365*DSIH**2*np.cos(DNU+DNU)+DC0365*DC0365)/(DC5023*DMOF71+DC0365)
    DFM1C = (1-10*DSHIH**2+15*DSHIH**2*DSHIH**2)*DCHIH**2/((1-10*DSHOM2+15*DSHOM2*DSHOM2)*DCHOM2)
    
    multiply_variables = np.stack([DND73,DND74,DND75,DND76,DND77,DND78,DND79,DFM1,DFK1,DFL2,DFK2,DFM1C])
    sel_cols = ['DND73','DND74','DND75','DND76','DND77','DND78','DND79','DFM1','DFK1','DFL2','DFK2','DFM1C']
    
    v0uf_allT = get_schureman_table()
    v0uf_sel = v0uf_allT.loc[const_list]
    f_i = np.ones(shape=(len(const_list), len(dood_date)))
    for variable, colname in zip(multiply_variables, sel_cols): #this loop is faster than array multiplications, since it only calculates the necessary factors (saves a lot of overhead)
        power = v0uf_sel[colname].values
        idnozero = power!=0
        f_i[idnozero,:] *= variable**power[idnozero][np.newaxis].T
    
    v0uf_M2 = v0uf_allT.loc[['M2']]
    f_i_M2 = np.ones(shape=(len(['M2']), len(dood_date)))
    for variable, colname in zip(multiply_variables, sel_cols): 
        power = v0uf_M2[colname].values
        idnozero = power!=0
        f_i_M2[idnozero,:] *= variable**power[idnozero][np.newaxis].T
    
    #create dataframe
    f_i_pd = pd.DataFrame(f_i, index=const_list)
    f_i_M2_pd = pd.DataFrame(f_i_M2, index = ['M2'])
    
    if xfac: #if variable is not: None, False, 0, more?
        f_i_pd = correct_fwith_xfac(f_i_pd, f_i_M2_pd, xfac=xfac)
    
    return f_i_pd


def correct_fwith_xfac(f_i_pd, f_i_M2_pd, xfac):
    """
    Correct f-values with xfactor, this definition is only ran when xfac=True.

    Parameters
    ----------
    f_i_pd : TYPE
        DESCRIPTION.
    f_i_M2_pd : TYPE
        DESCRIPTION.

    Returns
    -------
    f_i_pd : TYPE
        DESCRIPTION.

    """
    
    if isinstance(xfac,dict):
        print('xfac dictionary provided: %s'%(xfac))
        xfac_values = xfac
    else: #use default xfactor values
        xfac_values = {'MU2': 0.00, # 0 betekent knoopfactor uit, 1 is helemaal aan
                       'N2':  0.00,
                       'NU2': 0.80,
                       'M2':  0.53,
                       '2MN2':0.20,
                       'S2': -0.82,
                       'M4':  0.70,
                       'MS4': 0.00,
                       'M6':  0.75,
                       '2MS6':0.20,
                       'M8':  0.70,
                       '3MS8':0.60}
    
    for xfac_const in xfac_values.keys():
        if not xfac_const in f_i_pd.index:
            continue
        if all(f_i_pd.loc[xfac_const] == 1): # HET IS EEN NIET-KNOOPAFHANKELIJKE COMPONENT (F=1) (like S2)
            f_i_pd.loc[[xfac_const]] = xfac_values[xfac_const]*(f_i_M2_pd.loc[['M2'],:].values-1)+1 # # hvufea.f line 176. uit v0uf_M2 ipv v0uf_sel, want daar is M2 waarde nog niet gecorrigeerd met xfac (kan ook uit f_i_HAT komen maar daar zit M2 niet per definitie in)
        else: # HET IS EEN KNOOPAFHANKELIJKE COMPONENT (F#1)
            f_i_pd.loc[[xfac_const]] = xfac_values[xfac_const]*(f_i_pd.loc[[xfac_const],:]-1)+1 # uit document tom bogaard, en hvufea.f line 181. staat ook in hatyan gebruikershandleiding bij knoopfactoren (pag2-5)

    return f_i_pd



