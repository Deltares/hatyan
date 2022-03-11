# -*- coding: utf-8 -*-
"""
hatyan_core.py contains core components that wrap around schureman/foreman definitions

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

import pandas as pd
import datetime as dt
import numpy as np
import functools

from hatyan.schureman import get_schureman_freqs, get_schureman_v0, get_schureman_u, get_schureman_f, get_schureman_table
from hatyan.foreman import get_foreman_v0_freq, get_foreman_doodson_nodal_harmonic, get_foreman_shallowrelations, get_foreman_nodalfactors


@functools.lru_cache()
def check_requestedconsts(const_list_tuple,source):
    #TODO: move check to central location when part of hatyan_settings()?
    if source=='schureman':
        const_list_allforsource = get_schureman_table().index
    elif source=='foreman':
        foreman_doodson_harmonic, foreman_nodal_harmonic = get_foreman_doodson_nodal_harmonic()
        foreman_shallowrelations, list_shallowdependencies = get_foreman_shallowrelations()
        const_list_allforsource = pd.Series(foreman_doodson_harmonic.index.append(foreman_shallowrelations.index))
    
    bool_constavailable = pd.Series(const_list_tuple).isin(const_list_allforsource)   
    if not bool_constavailable.all():
        raise Exception(f'ERROR: not all requested components available in schureman harmonics or shallowrelations:\n{pd.Series(const_list_tuple).loc[~bool_constavailable]}')


def get_freqv0_generic(hatyan_settings, const_list, dood_date_mid, dood_date_start):

    if type(const_list) is str:
        const_list = get_const_list_hatyan(const_list)

    #retrieve frequency and v0
    print(f'freq is calculated at mid of period: {dood_date_mid[0]}')
    print(f'v0 is calculated for start of period: {dood_date_start[0]}')
    if hatyan_settings.source=='schureman':
        t_const_freq_pd = get_schureman_freqs(const_list, dood_date=dood_date_mid)
        v_0i_rad = get_schureman_v0(const_list, dood_date_start).T #at start of timeseries
    elif hatyan_settings.source=='foreman': #TODO: this is probably quite slow since foreman is not cached
        dummy, t_const_freq_pd = get_foreman_v0_freq(const_list=const_list, dood_date=dood_date_mid) #TODO: does this really matter, maybe just retrieve on dood_date_start? Otherwise maybe split definition?
        v_0i_rad, dummy = get_foreman_v0_freq(const_list=const_list, dood_date=dood_date_start)
        v_0i_rad = v_0i_rad.T
    
    return t_const_freq_pd, v_0i_rad


def get_uf_generic(hatyan_settings, const_list, dood_date_fu):

    #get f and u
    if hatyan_settings.nodalfactors:
        if len(dood_date_fu)==1: #dood_date_mid
            print('nodal factors (f and u) are calculated for center of period: %s'%(dood_date_fu[0]))
        else: #times_pred_all_pdDTI
            print('nodal factors (f and u) are calculated for all timesteps')
        if hatyan_settings.source=='schureman':
            f_i = get_schureman_f(xfac=hatyan_settings.xfac, const_list=const_list, dood_date=dood_date_fu).T
            u_i_rad = get_schureman_u(const_list=const_list, dood_date=dood_date_fu).T
        elif hatyan_settings.source=='foreman': #TODO: this is probably quite slow since foreman is not cached
            f_i, u_i_rad = get_foreman_nodalfactors(const_list=const_list, dood_date=dood_date_fu)
            f_i, u_i_rad = f_i.T, u_i_rad.T
    else:
        print('no nodal factors (f and u) correction applied (f=1, u=0)')
        f_i = pd.DataFrame(np.ones(len(const_list)),index=const_list).T
        u_i_rad = pd.DataFrame(np.zeros(len(const_list)),index=const_list).T
    return u_i_rad, f_i


def get_doodson_eqvals(dood_date, mode=None):
    """
    Berekent de doodson waardes T, S, H, P, N en P1 voor het opgegeven tijdstip.

    Parameters
    ----------
    dood_date : datetime.datetime or pandas.DateTimeIndex
        Date(s) on which the doodson values should be calculated.
    mode : str, optional
        Calculate doodson values with Schureman (hatyan default) or Sea Level Science by Pugh. The default is False.

    Returns
    -------
    dood_T_rad : TYPE
        Bij hatyan 180+, maar in docs niet. In hatyan fortran code is +/-90 in v ook vaak omgedraaid in vergelijking met documentation.
    dood_S_rad : TYPE
        Middelbare lengte van de maan. Mean longitude of moon
    dood_H_rad : TYPE
        Middelbare lengte van de zon. mean longitude of sun
    dood_P_rad : TYPE
        Middelbare lengte van het perieum van de maan. Longitude of lunar perigee
    dood_N_rad : TYPE
        Afstand van de klimmende knoop van de maan tot het lentepunt. Longitude of moons node
    dood_P1_rad : TYPE
        Middelbare lengte van het perieum van de zon. Longitude of solar perigee

    """
    
    DNUJE = 24*36525
    dood_tstart_sec, fancy_pddt = robust_timedelta_sec(dood_date)    
    
    dood_Tj = (dood_tstart_sec/3600+12)/(24*36525) #DTIJJE, #Aantal Juliaanse eeuwen ( = 36525 dagen) die zijn verstreken sinds 31 december 1899 12.00 uur GMT. Number of Julian centuries (36525 days) with respect to Greenwich mean noon, 31 December 1899 (Gregorian calendar)
    if mode=='freq':
        #speed in radians per hour, de afgeleiden van onderstaande functies
        dood_T_rad = np.array(len(dood_Tj)*[np.deg2rad(15)]) #360/24=15 degrees of earth rotation per hour
        dood_S_rad  = (8399.7092745 + 0.0000346*dood_Tj*2)/DNUJE
        dood_H_rad  = ( 628.3319500 + 0.0000052*dood_Tj*2)/DNUJE
        dood_P_rad  = (  71.0180412 - 0.0001801*dood_Tj*2)/DNUJE
        dood_N_rad  = np.ones(shape=dood_date.shape)*np.nan #np.array([np.nan])
        dood_P1_rad = (   0.0300053 + 0.0000079*dood_Tj*2)/DNUJE
    else: #for everything else
        #hatyan documentation. Zie ook p162 of Schureman book
        if fancy_pddt:
            dood_T_rad = np.array(np.deg2rad(180 + dood_date.hour*15.0+dood_date.minute*15.0/60).values)
        else:
            dood_T_rad = np.array([np.deg2rad(180 + x.hour*15.0+x.minute*15.0/60) for x in dood_date])
        dood_S_rad =  (4.7200089 + 8399.7092745*dood_Tj + 0.0000346*dood_Tj**2)
        dood_H_rad =  (4.8816280 + 628.3319500*dood_Tj  + 0.0000052*dood_Tj**2)
        dood_P_rad =  (5.8351526 + 71.0180412*dood_Tj   - 0.0001801*dood_Tj**2)
        dood_N_rad =  (4.5236016 - 33.7571463*dood_Tj   + 0.0000363*dood_Tj**2)
        dood_P1_rad = (4.9082295 + 0.0300053*dood_Tj    + 0.0000079*dood_Tj**2)
    
    doodson_pd = pd.DataFrame(np.stack([dood_T_rad, dood_S_rad, dood_H_rad, dood_P_rad, dood_N_rad, dood_P1_rad]),
                              index=['T','S','H','P','N','P1'])
    return doodson_pd


def robust_daterange_fromtimesextfreq(times_ext,timestep_min):
    """
    Generate daterange. Pandas pd.date_range and pd.DatetimeIndex only support times between 1677-09-21 and 2262-04-11, because of its ns accuracy.
    For dates outside this period, a list is generated and converted to a pd.Index instead.

    Parameters
    ----------
    times_ext : list of datetime.datetime
        DESCRIPTION.
    timestep_min : int
        DESCRIPTION.

    Returns
    -------
    times_pred_all_pdDTI : pd.DatetimeIndex or pd.Index
        DESCRIPTION.

    """
    
    if (np.array(times_ext) > pd.Timestamp.min).all() and (np.array(times_ext) < pd.Timestamp.max).all():
        times_pred_all = pd.date_range(start=times_ext[0], end=times_ext[-1], freq='%imin'%(timestep_min))
        times_pred_all_pdDTI = pd.DatetimeIndex(times_pred_all)
    else:
        print('WARNING: OutOfBoundsDatetime: Out of bounds nanosecond timestamp, so reverting to less fancy (slower) datetime ranges. Fancy ones are possible between %s and %s'%(pd.Timestamp.min, pd.Timestamp.max))
        td_mins = (times_ext[-1]-times_ext[0]).total_seconds()/60
        nsteps = int(td_mins/timestep_min)
        times_pred_all = pd.Series([times_ext[0]+dt.timedelta(minutes=x*timestep_min) for x in range(nsteps+1)])
        times_pred_all_pdDTI = pd.Index(times_pred_all)
    return times_pred_all_pdDTI


def robust_timedelta_sec(dood_date,refdate_dt=None):
    """
    Generate timedelta. Pandas pd.DatetimeIndex subtraction only supports times between 1677-09-21 and 2262-04-11, because of its ns accuracy.
    For dates outside this period, a list subtraction is generated.
    
    Parameters
    ----------
    dood_date : pd.Index or pd.DatetimeIndex
        DESCRIPTION.
    refdate_dt : dt.datetime, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    dood_tstart_sec : np.array
        DESCRIPTION.
    fancy_pddt : bool
        DESCRIPTION.

    """
    #TODO: check if this is still necessary in newer pandas versions
    #TODO: merging with robust_daterange_fromtimesextfreq() possible?
    
    if refdate_dt is None:
        refdate_dt = dt.datetime(1900,1,1)
        
    if (dood_date[[0,-1]] > pd.Timestamp.min).all() and (dood_date[[0,-1]] < pd.Timestamp.max).all():
        fancy_pddt = True
        dood_tstart_sec = ( dood_date-refdate_dt ).total_seconds().values
    else:
        fancy_pddt = False
        #print('WARNING: OutOfBoundsDatetime: Out of bounds nanosecond timestamp, so reverting to less fancy (slower) timedelta arrays. For fancy ones: min=%s and max=%s'%(pd.Timestamp.min, pd.Timestamp.max))
        dood_tstart_sec = np.array([(x-refdate_dt).total_seconds() for x in dood_date])
    return dood_tstart_sec, fancy_pddt


def get_lunarSLSIHO_fromsolar(v0uf_base): #TODO: iets simpeler implementatie in foreman.get_foreman_doodson_nodal_harmonic(), maar dit is ook prima
    
    #conversion to lunar for comparison with SLS and IHO
    v0uf_baseT_solar = v0uf_base.loc[['T','S','H','P','N','P1','EDN']].T
    v0uf_baseT_lunar = v0uf_baseT_solar.copy()
    v0uf_baseT_lunar['S'] = v0uf_baseT_solar['S'] + v0uf_baseT_solar['T'] #ib with relation ω1 =ω0 − ω2 +ω3 (stated in SLS book)
    v0uf_baseT_lunar['H'] = v0uf_baseT_solar['H'] - v0uf_baseT_solar['T'] #ic with relation ω1 =ω0 − ω2 +ω3 (stated in SLS book)
    #lunar IHO (compare to Sea Level Science book from Woodsworth and Pugh)
    v0uf_baseT_lunar_SLS = v0uf_baseT_lunar.copy()
    v0uf_baseT_lunar_SLS['EDN'] = -v0uf_baseT_lunar['EDN']%360 #klopt niet allemaal met tabel 4.1 uit SLS boek, moet dit wel?
    #lunar IHO (compare to c
    v0uf_baseT_lunar_IHO = v0uf_baseT_lunar.copy()
    v0uf_baseT_lunar_IHO['EDN'] = -v0uf_baseT_lunar['EDN']/90 + 5 # (-90 lijkt 6 in IHO lijst, 90 is 4, 180 is 7)
    v0uf_baseT_lunar_IHO.loc[v0uf_baseT_lunar_IHO['EDN']==3,'EDN'] = 7 # convert -180 (3) to +180 (7)
    v0uf_baseT_lunar_IHO[['S','H','P','N','P1']] += 5
    return v0uf_baseT_lunar, v0uf_baseT_lunar_SLS, v0uf_baseT_lunar_IHO


@functools.lru_cache() #TODO: get_foreman_v0_freq is way slower than get_schureman_freqs because it is live, therefore caching this entire function (can also cache hatyan.foreman.get_foreman_v0_freq() but then const_list should be a tuple instead of a list)
def get_full_const_list_withfreqs():

    dood_date = pd.DatetimeIndex([dt.datetime(1900,1,1)]) #dummy value
    
    v0uf_allT = get_schureman_table()
    freqs_pd_schu = get_schureman_freqs(const_list=v0uf_allT.index.tolist(),dood_date=dood_date)
    
    foreman_doodson_harmonic, foreman_nodal_harmonic = get_foreman_doodson_nodal_harmonic()
    foreman_shallowrelations, list_shallowdependencies = get_foreman_shallowrelations()
    const_list_foreman = foreman_doodson_harmonic.index.tolist() + foreman_shallowrelations.index.tolist()
    v0_pd_for,freqs_pd_for = get_foreman_v0_freq(const_list=const_list_foreman,dood_date=dood_date) #TODO: this is slower than schureman, cache it instead of this definition. But then first always retrieve everything (currently foreman only retrieves requested components)
    
    freqs_pd_combined = pd.concat([freqs_pd_schu[['freq']],freqs_pd_for],axis=0)
    bool_duplicated = freqs_pd_combined.index.duplicated(keep='first')
    freqs_pd = freqs_pd_combined.loc[~bool_duplicated] #drop duplicates for OQ2 and M7, different in foreman and schureman, but that close to each other that it does not matter for ordering
    full_const_list = freqs_pd
    
    return full_const_list


def sort_const_list(const_list):
    full_const_list_withfreqs = get_full_const_list_withfreqs()
    const_list_sorted = full_const_list_withfreqs.loc[const_list].sort_values('freq').index.tolist()
    return const_list_sorted


def get_const_list_hatyan(listtype, return_listoptions=False):
    """
    Definition of several hatyan components lists, taken from the tidegui initializetide.m code, often originating from corresponcence with Koos Doekes

    Parameters
    ----------
    listtype : str
        The type of the components list to be retrieved, options:
            - 'all': all available components in hatyan_python
            - 'all_originalorder': all 195 hatyan components in original hatyan-FORTRAN order
            - 'year': default list of 94 hatyan components
            - 'halfyear': list of 88 components to be used when analyzing approximately half a year of data
            - 'month': list of 21 components to be used when analyzing one month of data. If desired, the K1 component can be splitted in P1/K1, N2 in N2/Nu2, S2 in T2/S2/K2 and 2MN2 in Labda2/2MN2.
            - 'month_deepwater': list of 21 components to be used when analyzing one month of data for deep water (from tidegui).
            - 'springneap': list of 14 components to be used when analyzing one spring neap period (approximately 15 days) of data
            - 'day': list of 10 components to be used when analyzing one day
            - 'day_tidegui': list of 5 components to be used when analyzing one day ir two tidal cycles (from tidegui)
            - 'tidalcycle': list of 6 components to be used when analyzing one tidal cycle (approximately 12 hours and 25 minutes)
    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    const_list_hatyan : list of str
        A list of component names.

    """
    
    v0uf_allT = get_schureman_table()
    const_list_all = v0uf_allT.index.tolist()
    
    #TODO: add all_foreman list, first optimize it via caching
    const_lists_dict = {'all_schureman':
                            #alle binnen hatyan beschikbare componenten
                            #A0 en 195 componenten (plus potentially added components)
                            const_list_all,
                        'all_schureman_originalorder': #for writing numbers in component file
                            #alle binnen hatyan beschikbare componenten
                            #A0 en 195 componenten
                            ['A0','SA','SSA','MSM','MM','MSF','SM','MF','SNU2','SN','MFM','2SM','2SMN',
                            '2Q1','NJ1','SIGMA1','NUJ1','Q1','RO1','NUK1','O1','TAU1','MP1','M1B','M1C','M1D','M1A','M1','NO1','CHI1','LP1','PI1','TK1','P1','S1','K1','PSI1','RP1','FI1','KP1','THETA1','LABDAO1','J1',
                            '2PO1','SO1','OO1','KQ1','3MKS2','3MS2','OQ2','MNK2','MNS2','2ML2S2','2MS2K2','NLK2','2N2','MU2','2MS2','SNK2','N2','NU2','2KN2S2','OP2','MSK2','MPS2','M2','MSP2','MKS2','M2(KS)2','2SN(MK)2','LABDA2','2MN2','L2','L2A','L2B','2SK2','T2','S2','R2','K2','MSN2','ETA2','KJ2','MKN2','2KM(SN)2','2SM2','SKM2','2SNU2','3(SM)N2','2SN2','SKN2',
                            'MQ3','NO3','MO3','2MK3','2MP3','M3','SO3','MK3','2MQ3','SP3','SK3','K3','2SO3',
                            '4MS4','2MNS4','3MK4','MNLK4','3MS4','MSNK4','MN4','2MLS4','2MSK4','M4','2MKS4','SN4','3MN4','2SMK4','MS4','MK4','2SNM4','2MSN4','SL4','S4','SK4','2SMN4','3SM4','2SKM4',
                            'MNO5','3MK5','3MP5','M5','MNK5','2MP5','3MO5','MSK5','3KM5',
                            '2(MN)S6','3MNS6','2NM6','4MS6','2MSNK6','2MN6','2MNU6','3MSK6','M6','MSN6','4MN6','MNK6','MKNU6','2(MS)K6','2MS6','2MK6','2SN6','3MSN6','MKL6','2SM6','MSK6','S6',
                            '2MNO7','2NMK7','M7','2MSO7','MSKO7',
                            '2(MN)8','3MN8','3MNKS8','M8','2MSN8','2MNK8','3MS8','3MK8','2SNM8','MSNK8','2(MS)8','2MSK8','3SM8','2SMK8','S8',
                            '2(MN)K9','3MNK9','4MK9','3MSK9',
                            '4MN10','M10','3MSN10','4MS10','2(MS)N10','2MNSK10','3M2S10',
                            '4MSK11',
                            'M12','4MSN12','5MS12','3MNKS12','4M2S12'],
                    
                        'year':
                            #Bij analyse van een jaar wordt gebruik gemaakt
                            #van de 'standaardset' van 94 componenten aanbevolen
                            #Althans voor data van de Noordzee en omgeving.
                            #A0 en 94 componenten
                            #of which shallow water components: ['2MS6', '2SM2', '3MS4', '3MS8', '4MS10', 'M4', 'M6', 'M8', 'MS4']
                            ['A0','SA','SM',
                            'Q1','O1','M1C','P1','S1','K1',
                            '3MKS2','3MS2','OQ2','MNS2','2ML2S2','NLK2','MU2','N2','NU2','MSK2','MPS2','M2','MSP2','MKS2','LABDA2','2MN2','T2','S2','K2','MSN2','2SM2','SKM2',
                            'NO3','2MK3','2MP3','SO3','MK3','SK3',
                            '4MS4','2MNS4','3MS4','MN4','2MLS4','2MSK4','M4','3MN4','MS4','MK4','2MSN4','S4',
                            'MNO5','3MK5','2MP5','3MO5','MSK5','3KM5',
                            '3MNS6','2NM6','4MS6','2MN6','2MNU6','3MSK6','M6','MSN6','MKNU6','2MS6','2MK6','3MSN6','2SM6','MSK6',
                            '2MNO7','M7','2MSO7',
                            '2(MN)8','3MN8','M8','2MSN8','2MNK8','3MS8','3MK8','2(MS)8','2MSK8',
                            '3MNK9','4MK9','3MSK9',
                            '4MN10','M10','3MSN10','4MS10','2(MS)N10','3M2S10',
                            '4MSK11',
                            'M12','4MSN12','5MS12','4M2S12'],
                    
                        'halfyear':
                            #Bij analyse van een halfjaar wordt 
                            #aanbevolen 88 componenten te gebruiken, 
                            #nl. de standaardset minus
                            #S1, MSK2, MPS2, MSP2, MKS2, en T2.
                            #Hierbij nog enkele opmerkingen :
                            #De middelste vier componenten geven 
                            #eigenlijk de seizoensmodulatie van M2
                            #weer. Het is daarom logischer om bij
                            #analyse op een half jaar ook MSK2 en
                            #MKS2 weg te laten. 
                            #De schatting van de jaarlijkse component
                            #Sa is zo natuurlijk onbetrouwbaar, maar
                            #nauwelijks onbetrouwbaarder dan bij
                            #analyse op een jaar.
                            #Toepassing van componentsplitsing bij
                            #analyse op een half jaar geeft geheel
                            #onjuiste resultaten.
                            #Als de beschikbare periode korter is
                            #dan een half jaar, is het aan te bevelen
                            #deze in perioden van 1 maand op te splitsen. 
                            #A0 en 88 componenten
                            ['A0','SA','SM',
                            'Q1','O1','M1C','P1',
                            'K1','3MKS2','3MS2','OQ2','MNS2','2ML2S2','NLK2','MU2','N2','NU2','M2','LABDA2','2MN2','S2','K2','MSN2','2SM2','SKM2',
                            'NO3','2MK3','2MP3','SO3','MK3','SK3',
                            '4MS4','2MNS4','3MS4','MN4','2MLS4','2MSK4','M4','3MN4','MS4','MK4','2MSN4','S4',
                            'MNO5','3MK5','2MP5','3MO5','MSK5','3KM5',
                            '3MNS6','2NM6','4MS6','2MN6','2MNU6','3MSK6','M6','MSN6','MKNU6','2MS6','2MK6','3MSN6','2SM6','MSK6',
                            '2MNO7','M7','2MSO7',
                            '2(MN)8','3MN8','M8','2MSN8','2MNK8','3MS8','3MK8','2(MS)8','2MSK8',
                            '3MNK9','4MK9','3MSK9',
                            '4MN10','M10','3MSN10','4MS10','2(MS)N10','3M2S10',
                            '4MSK11',
                            'M12','4MSN12','5MS12','4M2S12'],
                    
                        'month':
                            #Bij analyse van 1 maand zijn de aanbevolen componenten:
                            #Q1, O1, K1
                            #3MS2, MNS2, Mu2, N2, M2, 2MN2, S2, 2SM2
                            #3MS4, MN4, M4, MS4
                            #2MN6, M6, 2MS6
                            #M8, 3MS8
                            #4MS10
                            #Als er in het programma een functie voor componentsplitsing is ingebouwd, kan
                            #men hiermee K1 splitsen in P1/K1, N2 in N2/Nu2, S2 in T2/S2/K2 en 2MN2 in Labda2/2MN2.
                            #A0 en 21 componenten
                            ['A0','Q1', 'O1', 'K1',
                            '3MS2', 'MNS2', 'MU2', 'N2', 'M2', '2MN2', 'S2', '2SM2',
                            '3MS4', 'MN4', 'M4', 'MS4',
                            '2MN6', 'M6', '2MS6',
                            'M8', '3MS8',
                            '4MS10'],
                            
                        'month_deepwater':
                            # for deep water L2 is better connected than 2MN2 to LABDA2
                            ['A0','Q1', 'O1', 'K1',
                            '3MS2', 'MNS2', 'MU2', 'N2', 'M2', 'L2', 'S2', '2SM2',
                            '3MS4', 'MN4', 'M4', 'MS4',
                            '2MN6', 'M6', '2MS6',
                            'M8', '3MS8',
                            '4MS10'],
                    
                        'springneap':
                            #Bij analyse van 1 spring-doodtijcyclus zijn de aanbevolen componenten:
                            #O1, K1
                            #Mu2, M2, S2, 2SM2
                            #3MS4, M4, MS4
                            #M6, 2MS6
                            #M8, 3MS8
                            #4MS10
                            #De resultaten bij analyse van 15 dagen zijn
                            #erg pover. Men kan zo niet eens N2, en de
                            #samenstellingen daarvan als MN4 etc.,
                            #afsplitsen, zodat de schattingen van M2 en
                            #S2 onbetrouwbaar zijn. 
                            #A0 en 14 componenten
                            ['A0','O1','K1',
                            'MU2','M2','S2','2SM2',
                            '3MS4','M4','MS4',
                            'M6','2MS6',
                            'M8','3MS8',
                            '4MS10'],
                            
                        'day':
                            #harmonic constituents for 1-day analysis
                            #i.e. 2 tidal periods ['A0','M2','M4','M6','M8','O1'],
                            #A0 en 10 componenten
                            ['A0','M1','M2','M3','M4','M5','M6','M7','M8','M10','M12'],
                            
                        'tidalcycle':
                            #harmonic constituents for 1-tide analysis
                            #A0 en 6 componenten
                            ['A0','M2','M4','M6','M8','M10','M12']
                        }
    
    const_list_options = const_lists_dict.keys()
    if listtype in const_list_options:
        const_list_hatyan = const_lists_dict[listtype]
    else:
        raise Exception('ERROR: listtype "%s" is not implemented in get_const_list_hatyan, choose from: %s'%(str(listtype),const_list_options))
        
    if return_listoptions:
        return const_list_hatyan, const_list_options
    else:
        return const_list_hatyan


