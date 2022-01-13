# - * - coding: utf - 8 - * - 
"""
schureman_core.py contains definitions with data for the schureman constituents.

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
file_path = os.path.realpath(__file__)


def get_v0uf_sel(const_list):
    """
    get_v0uf_sel

    Parameters
    ----------
    const_list : TYPE
        DESCRIPTION.
    
    Returns
    -------
    v0uf_sel : TYPE
        DESCRIPTION.

    """
    
    import pandas as pd

    v0uf_allT = get_schureman_table()
    
    const_list_pd = pd.Series(const_list,index=const_list)
    const_list_avaibool = const_list_pd.isin(v0uf_allT.index)
    
    if not const_list_avaibool.all():
        const_list_notavailable = const_list_avaibool.loc[~const_list_avaibool]
        raise Exception('ERROR: not all requested constituents are available:\n%s'%(const_list_notavailable))
    else:
        v0uf_sel = v0uf_allT.loc[const_list]
    
    return v0uf_sel


@functools.lru_cache() #TODO: foreman is way slower, so caching this entire function
def full_const_list():
    import pandas as pd
    import datetime as dt

    dood_date = pd.DatetimeIndex([dt.datetime(1900,1,1)]) #dummy value
    
    freqs_pd_schu = get_schureman_freqs(const_list='all',dood_date=dood_date)
    
    from hatyan.foreman_core import get_foreman_v0freq_fromfromharmonicdood, get_foreman_shallowrelations, get_foreman_v0_freq
    v_0i_rad_harmonic_pd = get_foreman_v0freq_fromfromharmonicdood() #list with only harmonic components with more precision than file
    foreman_shallowrelations = get_foreman_shallowrelations()
    const_list_foreman = v_0i_rad_harmonic_pd.index.tolist() + foreman_shallowrelations.index.tolist()
    v0_pd_for,freqs_pd_for = get_foreman_v0_freq(const_list=const_list_foreman,dood_date=dood_date) #TODO: this is slower than schureman, cache it instead of this definition. But then first always retrieve everything (currently foreman only retrieves requested components)
    
    freqs_pd_combined = pd.concat([freqs_pd_schu[['freq']],freqs_pd_for],axis=0)
    bool_duplicated = freqs_pd_combined.index.duplicated(keep='first')
    freqs_pd = freqs_pd_combined.loc[~bool_duplicated] #drop duplicates for OQ2 and M7, different in foreman and schureman, but that close to each other that it does not matter for ordering
    full_const_list = freqs_pd
    
    return full_const_list


def sort_const_list(const_list):
    from hatyan.schureman_core import full_const_list #import necessary since it is a cached function
    
    full_const_list = full_const_list()
    const_list_sorted = full_const_list.loc[const_list].sort_values('freq').index.tolist()
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
    
    const_lists_dict = {'all':
                            #alle binnen hatyan beschikbare componenten
                            #A0 en 195 componenten (plus potentially added components)
                            const_list_all,
                        'all_originalorder':
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
        raise Exception('ERROR: listtype %s is not implemented in get_const_list_hatyan, choose from: %s'%(str(listtype),const_list_options))
        
    if return_listoptions:
        return const_list_hatyan, const_list_options
    else:
        return const_list_hatyan


def anapred_get_freqv0uf(hatyan_settings, const_list, dood_date_start, dood_date_mid, times_pred_all_pdDTI):
    import numpy as np
    import pandas as pd
    from hatyan.schureman_core import get_schureman_freqs, get_schureman_v0, get_schureman_u, get_schureman_f
    from hatyan.foreman_core import get_foreman_v0_freq, get_foreman_nodalfactors
    
    #retrieve frequency and v0
    print('v0 is calculated for start of period: %s'%(dood_date_start[0]))
    if hatyan_settings.source=='schureman':
        t_const_freq_pd, t_const_speed_all = get_schureman_freqs(const_list, dood_date=dood_date_mid, return_allraw=True)
        v_0i_rad = get_schureman_v0(const_list, dood_date_start).T #at start of timeseries
    elif hatyan_settings.source=='foreman':
        v_0i_rad, t_const_freq_pd = get_foreman_v0_freq(const_list=const_list, dood_date=dood_date_start)
        #t_const_speed_all = t_const_freq_pd['freq'].values[:,np.newaxis]*(2*np.pi)
        v_0i_rad = v_0i_rad.T
    
    #get f and u
    if hatyan_settings.nodalfactors:
        if hatyan_settings.fu_alltimes:
            print('nodal factors (f and u) are calculated for all timesteps')
            dood_date_fu = times_pred_all_pdDTI
        else:
            print('nodal factors (f and u) are calculated for center of period: %s'%(dood_date_mid[0]))
            dood_date_fu = dood_date_mid
        if hatyan_settings.source=='schureman':
            f_i = get_schureman_f(xfac=hatyan_settings.xfac, const_list=const_list, dood_date=dood_date_fu).T
            u_i_rad = get_schureman_u(const_list=const_list, dood_date=dood_date_fu).T
        elif hatyan_settings.source=='foreman':
            f_i, u_i_rad = get_foreman_nodalfactors(const_list=const_list, dood_date=dood_date_fu)
            f_i, u_i_rad = f_i.T, u_i_rad.T
    else:
        print('no nodal factors (f and u) correction applied (f=1, u=0)')
        f_i = pd.DataFrame(np.ones(len(const_list)),index=const_list).T
        u_i_rad = pd.DataFrame(np.zeros(len(const_list)),index=const_list).T
    return t_const_freq_pd, v_0i_rad, u_i_rad, f_i


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
    #import datetime as dt
    import numpy as np
    
    DNUJE = 24*36525
    dood_tstart_sec, fancy_pddt = robust_timedelta_sec(dood_date)    
    
    dood_Tj = (dood_tstart_sec/3600+12)/(24*36525) #DTIJJE, #Aantal Juliaanse eeuwen ( = 36525 dagen) die zijn verstreken sinds 31 december 1899 12.00 uur GMT. Number of Julian centuries (36525 days) with respect to Greenwich mean noon, 31 December 1899 (Gregorian calendar)
    if mode=='freq':
        #speed in radians per hour, de afgeleiden van onderstaande functies
        dood_T_rad = np.array(len(dood_Tj)*[np.deg2rad(15)]) #360/24=15 degrees of earth rotation per hour
        dood_S_rad  = (8399.7092745 + 0.0000346*dood_Tj*2)/DNUJE
        dood_H_rad  = ( 628.3319500 + 0.0000052*dood_Tj*2)/DNUJE
        dood_P_rad  = (  71.0180412 - 0.0001801*dood_Tj*2)/DNUJE
        dood_N_rad  = np.array([np.nan])
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
    
    return dood_T_rad, dood_S_rad, dood_H_rad, dood_P_rad, dood_N_rad, dood_P1_rad


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
    
    #import datetime as dt
    import numpy as np
    
    dood_tstart_sec, fancy_pddt = robust_timedelta_sec(dood_date)    
    dood_Tj = (dood_tstart_sec/3600+12)/(24*36525)
    
    #bercon.f: HET BEREKENEN VAN DE 'CONSTANTEN' .0365, .1681 EN .5023, DIE GEBRUIKT WORDEN BIJ DE BEREKENING VAN DE U- EN F-FACTOREN VAN DE GETIJCOMPONENTEN K1 EN K2
    #327932: ratio of mass of sun to combined mass of earth and moon
    #81.53: ratio of mass of earth to mass of moon
    DAGC = 0.01657                            #A/C >> ZIE BLZ. 162 VAN SCHUREMAN, mean lunar parallax in respect to mean radius [rad]
    DAGC1 = 0.00004261                        #A/C1 >> ZIE BLZ. 162 VAN SCHUREMAN, mean solar parallax in respect to mean radius [rad]
    DE = 0.054900489                          #E >> ZIE BLZ. 162 VAN SCHUREMAN, eccentricity of moons orbit
    DE1 = 0.01675104-0.0000418*dood_Tj        #E1 >> ZIE BLZ. 162 VAN SCHUREMAN, eccentricity of earths orbit, epoch 1 Jan 1900
    DMGE = 1./81.53                           #M/E >> ZIE BLZ. 162 VAN SCHUREMAN, mass of moon /mass of earth
    DSGE = 82.53/81.53*327932                 #S/E >> ZIE BLZ. 162 VAN SCHUREMAN, mass of sun / mass of earth
    
    DU = DMGE*DAGC**3                         #U >> ZIE BLZ. 162 VAN SCHUREMAN, basic factor
    DU1 = DSGE*DAGC1**3                       #U1 >> ZIE BLZ. 162 VAN SCHUREMAN, solar coefficient
    DSACCE = DU1/DU                           #S' >> ZIE BLZ. 162 VAN SCHUREMAN, solar factor
    DOMEGA = np.deg2rad(23+27./60.+8.26/3600) #OMEGA >> ZIE BLZ. 162 VAN SCHUREMAN, Obliquity of the Ecliptic, epoch 1 Jan 1900, 23.45229 graden [rad]
    DIKL = np.deg2rad(5+8./60.+43.3546/3600)  #i >> ZIE BLZ. 162 VAN SCHUREMAN, Inclination of moons orbit to plane of ecliptic, epoch 1 Jan 1900 (?), 5.1453762 graden [rad]
    
    DCOFSI = (0.5+0.75*DE1**2)*DSACCE         #COEFFICIENT VAN DE SINUSTERMEN IN (217) EN (219) OP BLZ. 45 VAN SCHUREMAN

    DC0365 = DCOFSI*np.sin(DOMEGA)**2
    DC1681 = DCOFSI*np.sin(2*DOMEGA)
    DC5023 = 0.5+0.75*DE*DE
    
    return DOMEGA, DIKL, DC1681, DC5023, DC0365


def get_lunarSLSIHO_fromsolar(v0uf_base):
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

    file_schureman_harmonic = os.path.join(os.path.dirname(file_path),'data_schureman_harmonic.csv')
    v0uf_baseT = pd.read_csv(file_schureman_harmonic,comment='#',skipinitialspace=True,index_col='component')
    v0uf_base = v0uf_baseT.T
    #v0uf_base.index = index_v0 + index_u + index_f + index_fstr
    
    file_schureman_shallowrelations = os.path.join(os.path.dirname(file_path),'data_schureman_shallowrelations.csv')
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
    
    #v0uf_allT_lunar, v0uf_allT_lunar_SLS, v0uf_allT_lunar_IHO = get_lunarSLSIHO_fromsolar(v0uf_all)
    
    return v0uf_allT


def get_schureman_freqs(const_list, dood_date=None, sort_onfreq=True, return_allraw=False):
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
    import pandas as pd
    import numpy as np
    import datetime as dt
    
    if type(const_list) is str:
        const_list = get_const_list_hatyan(const_list)
    elif type(const_list) is not list:
        const_list = const_list.tolist()

    const_list_hatyan195, const_listoptions = get_const_list_hatyan('all', return_listoptions=True)
    
    if dood_date is None:
        dood_date = pd.DatetimeIndex([dt.datetime(1900,1,1)]) #dummy value
    
    T_rad_freq, S_rad_freq, H_rad_freq, P_rad_freq, N_rad_freq, P1_rad_freq = get_doodson_eqvals(dood_date=dood_date, mode='freq') #N is not used here
    multiply_variables = np.stack([T_rad_freq,S_rad_freq, H_rad_freq, P_rad_freq, P1_rad_freq])
    
    v0uf_sel = get_v0uf_sel(const_list=const_list)
    v0uf_sel_freq = v0uf_sel[['T','S','H','P','P1']]
    
    DOMEGA_speed = np.dot(v0uf_sel_freq.values,multiply_variables)
    
    t_const_speed = DOMEGA_speed[:,0]
    t_const_freq = t_const_speed/(2*np.pi) #aantal rotaties per uur, freq
    np.seterr(divide='ignore') #suppress divide by 0 warning
    t_const_perds = 1/t_const_freq #period [hr]
    t_const_angfreqs = 360/t_const_perds #angfreq [deg/hr]
    freq_pd = pd.DataFrame({'freq':t_const_freq, 'angvelo [rad/hr]':t_const_speed, 'period [hr]':t_const_perds, 'angfreq [deg/hr]':t_const_angfreqs}, index=const_list)
    if sort_onfreq:
        freq_pd = freq_pd.sort_values(by='freq')
    
    if return_allraw:
        return freq_pd, DOMEGA_speed
    else:
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
    import numpy as np
    import pandas as pd
    
    T_rad, S_rad, H_rad, P_rad, N_rad, P1_rad = get_doodson_eqvals(dood_date=dood_date) #N is not used here
    multiply_variables = np.stack([T_rad,S_rad, H_rad, P_rad, P1_rad])
    
    v0uf_sel = get_v0uf_sel(const_list=const_list)
    v0uf_sel_v = v0uf_sel[['T','S','H','P','P1']]
    
    DV0 = np.dot(v0uf_sel_v.values,multiply_variables) + np.deg2rad(v0uf_sel['EDN']).values[np.newaxis].T
    
    DV0_pd = pd.DataFrame(DV0)
    DV0_pd.index = const_list
    
    return DV0_pd


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
    import numpy as np
    import pandas as pd
    
    if isinstance(const_list, pd.Series) or isinstance(const_list, pd.core.indexes.base.Index):
        const_list = const_list.tolist()
    
    T_rad, S_rad, H_rad, P_rad, N_rad, P1_rad = get_doodson_eqvals(dood_date=dood_date)
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
    
    v0uf_sel = get_v0uf_sel(const_list=const_list)
    v0uf_sel_u = v0uf_sel[['DKSI','DNU', 'DQ', 'DQU', 'DR', 'DUK1', 'DUK2']]

    DU = np.dot(v0uf_sel_u.values,multiply_variables)
    DU += np.pi
    DU = np.remainder(DU,2*np.pi)
    DU -= np.pi

    #create dataframe
    DU_pd = pd.DataFrame(DU)
    DU_pd.index = const_list
    
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
    import numpy as np
    import pandas as pd
    
    if isinstance(const_list, pd.Series) or isinstance(const_list, pd.core.indexes.base.Index):
        const_list = const_list.tolist()

    T_rad, S_rad, H_rad, P_rad, N_rad, P1_rad = get_doodson_eqvals(dood_date=dood_date)

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
    
    v0uf_sel = get_v0uf_sel(const_list=const_list)
    f_i = np.ones(shape=(len(const_list), len(dood_date)))
    for variable, colname in zip(multiply_variables, sel_cols): #this loop is faster than array multiplications, since it only calculates the necessary factors (saves a lot of overhead)
        power = v0uf_sel[colname].values
        idnozero = power!=0
        f_i[idnozero,:] *= variable**power[idnozero][np.newaxis].T
    
    v0uf_M2 = get_v0uf_sel(const_list=['M2'])
    f_i_M2 = np.ones(shape=(len(['M2']), len(dood_date)))
    for variable, colname in zip(multiply_variables, sel_cols): 
        power = v0uf_M2[colname].values
        idnozero = power!=0
        f_i_M2[idnozero,:] *= variable**power[idnozero][np.newaxis].T
    
    #create dataframe
    f_i_pd = pd.DataFrame(f_i)
    f_i_pd.index = const_list
    f_i_M2_pd = pd.DataFrame(f_i_M2)
    f_i_M2_pd.index = ['M2']
    
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
    
    const_list = f_i_pd.index.tolist()
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
        #print(xfac_const)
        if xfac_const in const_list:
            if all(f_i_pd.loc[xfac_const] == 1): # HET IS EEN NIET-KNOOPAFHANKELIJKE COMPONENT (F=1) (like S2)
                f_i_pd.loc[[xfac_const]] = xfac_values[xfac_const]*(f_i_M2_pd.loc[['M2'],:].values-1)+1 # # hvufea.f line 176. uit v0uf_M2 ipv v0uf_sel, want daar is M2 waarde nog niet gecorrigeerd met xfac (kan ook uit f_i_HAT komen maar daar zit M2 niet per definitie in)
            else: # HET IS EEN KNOOPAFHANKELIJKE COMPONENT (F#1)
                f_i_pd.loc[[xfac_const]] = xfac_values[xfac_const]*(f_i_pd.loc[[xfac_const],:]-1)+1 # uit document tom bogaard, en hvufea.f line 181. staat ook in hatyan gebruikershandleiding bij knoopfactoren (pag2-5)

    return f_i_pd


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
    import numpy as np
    import pandas as pd
    import datetime as dt
    
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
    import pandas as pd
    import datetime as dt
    import numpy as np
    
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


