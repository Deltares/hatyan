# - * - coding: utf - 8 - * - 
"""
hatyan_core.py contains definitions with data for the hatyan constituents.

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
import functools #to Memoize v0uf table (https://en.wikipedia.org/wiki/Memoization)


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

    v0uf_allT = calcwrite_baseforv0uf()
    
    const_list_pd = pd.Series(const_list,index=const_list)
    const_list_avaibool = const_list_pd.isin(v0uf_allT.index)
    
    if not const_list_avaibool.all():
        const_list_notavailable = const_list_avaibool.loc[~const_list_avaibool]
        raise Exception('ERROR: not all requested constituents are available:\n%s'%(const_list_notavailable))
    else:
        v0uf_sel = v0uf_allT.loc[const_list]
    
    return v0uf_sel


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
    v0uf_allT = calcwrite_baseforv0uf()
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


def get_hatyan_constants(dood_date):
    """
    get_hatyan_constants

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


@functools.lru_cache() #TODO: phase out pickle file
def calcwrite_baseforv0uf():
    """
    Calculate and write table for all constituents. Alternatively (faster), this data can be read from data_components_hatyan.pkl. Whether to calculate or read the table is controlled with the v0uf_calculatewrite boolean in the top of this script 

    Returns
    -------
    v0uf_all : TYPE
        DESCRIPTION.
    
    """
        
    #TODO: v0uf_base.T naar csv file en van daaruit inlezen. feqsstr als comment en ook comments meegeven
    #TODO: hatyan naar schureman hernoemen
    # in solar days T=Cs=w0 (s is van solar, er is ook een Cl=w1=ia for lunar), S=w2=ib, H=w3=ic, P=w4=id, N=w5=ie, P1=w6=if, EDN=ExtendedDoodsonNumber
    # N zit is hier altijd 0 (regression of moons node), want dat wordt apart als knoopfactor gedaan
    #                         comp       T, S, H, P, N,P1,EDN,     u eqs                    f eqs 73-79           f M1 K1 L2 K2 M1C feqsstr          #comparison of v0 columns with lunar conversion to SLS and IHO
    v0uf_base = pd.DataFrame({'A0':     [0, 0, 0, 0, 0, 0,  0,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0,[       ]],  #  0
                              'SA':     [0, 0, 1, 0, 0, 0,  0,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0,[       ]],  #  1 #SA 2nd option from IHO
                              'SA_IHO1':[0, 0, 1, 0, 0,-1,  0,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0,[       ]],  #    #newly added: SA 1st option from IHO (used in SLS and t_tide). Nulpunt v0 ligt dicht bij jaarwisseling (2015,1,3,21,50,0)
                              'SA_IHO2':[0, 0, 1, 0, 0, 0,  0,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0,[       ]],  #    #newly added: SA 2nd option from IHO (used in hatyan). Nulpunt v0 ligt dicht bij equinox/lentepunt (2015,3,22,19,50,0)
                              'SSA':    [0, 0, 2, 0, 0, 0,  0,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0,[       ]],  #  2
                              'MSM':    [0, 1,-2, 1, 0, 0,  0,     0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0,['DND73']],  #  3 #not in SLS. IHO also called MNum
                              'MM':     [0, 1, 0,-1, 0, 0,  0,     0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0,['DND73']],  #  4 #not in SLS
                              'MSF':    [0, 2,-2, 0, 0, 0,  0,     0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0,['DND73']],  #  5 #not in SLS
                              'MF':     [0, 2, 0, 0, 0, 0,  0,    -2, 0, 0, 0, 0, 0, 0,     0, 1, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0,['DND74']],  #  7 
                              'MFM':    [0, 3, 0,-1, 0, 0,  0,    -2, 0, 0, 0, 0, 0, 0,     0, 1, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0,['DND74']],  # 10 #corresponds to MTM in FES2014?
                              'MSQM':   [0, 4,-2, 0, 0, 0,  0,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0,[       ]],  #    #newly added: MSQM from IHO, unsure about nodal factor
                              '2Q1':    [1,-4, 1, 2, 0, 0, 90,     2,-1, 0, 0, 0, 0, 0,     0, 0, 1, 0, 0, 0, 0,     0, 0, 0, 0, 0,['DND75']],  # 13
                              'SIGMA1': [1,-4, 3, 0, 0, 0, 90,     2,-1, 0, 0, 0, 0, 0,     0, 0, 1, 0, 0, 0, 0,     0, 0, 0, 0, 0,['DND75']],  # 15
                              'Q1':     [1,-3, 1, 1, 0, 0, 90,     2,-1, 0, 0, 0, 0, 0,     0, 0, 1, 0, 0, 0, 0,     0, 0, 0, 0, 0,['DND75']],  # 17
                              'RO1':    [1,-3, 3,-1, 0, 0, 90,     2,-1, 0, 0, 0, 0, 0,     0, 0, 1, 0, 0, 0, 0,     0, 0, 0, 0, 0,['DND75']],  # 18 #also called RHO1
                              'O1':     [1,-2, 1, 0, 0, 0, 90,     2,-1, 0, 0, 0, 0, 0,     0, 0, 1, 0, 0, 0, 0,     0, 0, 0, 0, 0,['DND75']],  # 20
                              'TAU1':   [1,-2, 3, 0, 0, 0,-90,     0,-1, 0, 0, 0, 0, 0,     0, 0, 0, 1, 0, 0, 0,     0, 0, 0, 0, 0,['DND76']],  # 21 #not in SLS
                              'M1B':    [1,-1, 1,-1, 0, 0,-90,     2,-1, 0, 0, 0, 0, 0,     0, 0, 1, 0, 0, 0, 0,     0, 0, 0, 0, 0,['DND75']],  # 23 #not in SLS. M1B 2nd option from IHO
                             'M1B_IHO1':[1,-1, 1,-1, 0, 0, 90,     2,-1, 0, 0, 0, 0, 0,     0, 0, 1, 0, 0, 0, 0,     0, 0, 0, 0, 0,['DND75']],  #    #newly added: M1B 1st option from IHO
                             'M1B_IHO2':[1,-1, 1,-1, 0, 0,-90,     2,-1, 0, 0, 0, 0, 0,     0, 0, 1, 0, 0, 0, 0,     0, 0, 0, 0, 0,['DND75']],  #    #newly added: M1B 2nd option from IHO (used in hatyan)
                              'M1C':    [1,-1, 1, 0, 0, 0,  0,     1,-1, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 1,['DFM1C']],  # 24 #not in SLS
                              'M1D':    [1,-1, 1, 0, 0, 0,-90,     1,-1, 1, 0, 0, 0, 0,     0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,['DFM1' ]],  # 25 #not in SLS
                              'M1A':    [1,-1, 1, 1, 0, 0,-90,     0,-1, 0, 0, 0, 0, 0,     0, 0, 0, 1, 0, 0, 0,     0, 0, 0, 0, 0,['DND76']],  # 26 #not in SLS
                              'M1':     [1,-1, 1, 1, 0, 0,-90,     0,-1, 0,-1, 0, 0, 0,     0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,['DFM1' ]],  # 27 #M1 3rd option from IHO
                              'M1_IHO1':[1,-1, 1, 0, 0, 0,-90,     0,-1, 0,-1, 0, 0, 0,     0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,['DFM1' ]],  #    #newly added: M1 1st option from IHO
                              'M1_IHO2':[1,-1, 1, 0, 0, 0,180,     0,-1, 0,-1, 0, 0, 0,     0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,['DFM1' ]],  #    #newly added: M1 2nd option from IHO (used in SLS
                              'M1_IHO3':[1,-1, 1, 1, 0, 0,-90,     0,-1, 0,-1, 0, 0, 0,     0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,['DFM1' ]],  #    #newly added: M1 3rd option from IHO (used in hatyan)
                              'CHI1':   [1,-1, 3,-1, 0, 0,-90,     0,-1, 0, 0, 0, 0, 0,     0, 0, 0, 1, 0, 0, 0,     0, 0, 0, 0, 0,['DND76']],  # 29 #lunar conversion EDN different than SLS
                              'PI1':    [1, 0,-2, 0, 0, 1, 90,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0,[       ]],  # 31
                              'P1':     [1, 0,-1, 0, 0, 0, 90,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0,[       ]],  # 33
                              'S1':     [1, 0, 0, 0, 0, 0,  0,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0,[       ]],  # 34 #S1 1st option from IHO
                              'S1_IHO1':[1, 0, 0, 0, 0, 0,  0,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0,[       ]],  #    #newly added: S1 1st option from IHO (used in hatyan)
                              'S1_IHO2':[1, 0, 0, 0, 0, 0,180,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0,[       ]],  #    #newly added: S1 2nd option from IHO 
                              'S1_IHO3':[1, 0, 0, 0, 0, 0,-90,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0,[       ]],  #    #newly added: S1 3rd option from IHO (used in SLS)
                              'K1':     [1, 0, 1, 0, 0, 0,-90,     0, 0, 0, 0, 0, 1, 0,     0, 0, 0, 0, 0, 0, 0,     0, 1, 0, 0, 0,['DFK1' ]],  # 35 #S=w2=ib is 0 instead of 1 in SLS (mistake in SLS?). K1 2nd option from IHO. 
                              'K1_IHO1':[1, 0, 1, 0, 0, 0,  0,     0, 0, 0, 0, 0, 1, 0,     0, 0, 0, 0, 0, 0, 0,     0, 1, 0, 0, 0,['DFK1' ]],  #    #newly added: K1 1st option from IHO
                              'K1_IHO2':[1, 0, 1, 0, 0, 0,-90,     0, 0, 0, 0, 0, 1, 0,     0, 0, 0, 0, 0, 0, 0,     0, 1, 0, 0, 0,['DFK1' ]],  #    #newly added: K1 2nd option from IHO (used in hatyan)
                              'PSI1':   [1, 0, 2, 0, 0,-1,-90,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0,[       ]],  # 36 #EDN 270 instead of 90 in SLS (mistake in SLS?)
                              'FI1':    [1, 0, 3, 0, 0, 0,-90,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0,[       ]],  # 38 #EDN 270 instead of 90 in SLS (mistake in SLS?). IHO also called PHI1
                              'THETA1': [1, 1,-1, 1, 0, 0,-90,     0,-1, 0, 0, 0, 0, 0,     0, 0, 0, 1, 0, 0, 0,     0, 0, 0, 0, 0,['DND76']],  # 40 #EDN 0 instead of 90 in SLS (mistake in SLS?).
                              'J1':     [1, 1, 1,-1, 0, 0,-90,     0,-1, 0, 0, 0, 0, 0,     0, 0, 0, 1, 0, 0, 0,     0, 0, 0, 0, 0,['DND76']],  # 42 #EDN 270 instead of 90 in SLS (mistake in SLS?).
                              'OO1':    [1, 2, 1, 0, 0, 0,-90,    -2,-1, 0, 0, 0, 0, 0,     0, 0, 0, 0, 1, 0, 0,     0, 0, 0, 0, 0,['DND77']],  # 45
                              'EPS2':   [2,-5, 4, 1, 0, 0,  0,     2,-2, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 1, 0,     0, 0, 0, 0, 0,['DND78']],  #    #newly added: EPS2 (same v0/freq as MNS2 and same u/f as M2)
                              '2N2':    [2,-4, 2, 2, 0, 0,  0,     2,-2, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 1, 0,     0, 0, 0, 0, 0,['DND78']],  # 55
                              'MU2':    [2,-4, 4, 0, 0, 0,  0,     2,-2, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 1, 0,     0, 0, 0, 0, 0,['DND78']],  # 56
                              'N2':     [2,-3, 2, 1, 0, 0,  0,     2,-2, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 1, 0,     0, 0, 0, 0, 0,['DND78']],  # 59
                              'NU2':    [2,-3, 4,-1, 0, 0,  0,     2,-2, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 1, 0,     0, 0, 0, 0, 0,['DND78']],  # 60
                              'MA2':    [2,-2, 1, 0, 0, 0,  0,     2,-2, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 1, 0,     0, 0, 0, 0, 0,['DND78']],  #    #newly added from IHO table: MA2, α2, H1. Corresponds to MPS2 (except for EDN)
                              'M2':     [2,-2, 2, 0, 0, 0,  0,     2,-2, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 1, 0,     0, 0, 0, 0, 0,['DND78']],  # 65
                              'MB2':    [2,-2, 3, 0, 0, 0,  0,     2,-2, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 1, 0,     0, 0, 0, 0, 0,['DND78']],  #    #newly added from IHO table: MB2, β2, H2, Ma2 and MA2*. Corresponds to MSP2 (except for EDN)
                              'LABDA2': [2,-1, 0, 1, 0, 0,180,     2,-2, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 1, 0,     0, 0, 0, 0, 0,['DND78']],  # 70 IHO also called LAMBDA2, FES2014 calls it LA2
                              'L2':     [2,-1, 2,-1, 0, 0,180,     2,-2, 0, 0,-1, 0, 0,     0, 0, 0, 0, 0, 0, 0,     0, 0, 1, 0, 0,['DFL2' ]],  # 72
                              'L2A':    [2,-1, 2,-1, 0, 0,180,     2,-2, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 1, 0,     0, 0, 0, 0, 0,['DND78']],  # 73 #not in SLS
                              'L2B':    [2,-1, 2, 1, 0, 0,  0,     0,-2, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 0, 1,     0, 0, 0, 0, 0,['DND79']],  # 74 #not in SLS
                              'T2':     [2, 0,-1, 0, 0, 1,  0,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0,[       ]],  # 76
                              'S2':     [2, 0, 0, 0, 0, 0,  0,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0,[       ]],  # 77
                              'R2':     [2, 0, 1, 0, 0,-1,180,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0,[       ]],  # 78
                              'K2':     [2, 0, 2, 0, 0, 0,  0,     0, 0, 0, 0, 0, 0, 1,     0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 1, 0,['DFK2' ]],  # 79
                              'ETA2':   [2, 1, 2,-1, 0, 0,  0,     0,-2, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0, 0, 1,     0, 0, 0, 0, 0,['DND79']],  # 81 #not in SLS
                              'M3':     [3,-3, 3, 0, 0, 0,  0,     3,-3, 0, 0, 0, 0, 0,     0, 0, 0, 0, 0,1.5,0,     0, 0, 0, 0, 0,['78f1p5']]})# 96 #EDN 180 instead of 0 in SLS and IHO (mistake in hatyan?)
    index_v0 = ['T','S','H','P','N','P1','EDN']
    index_u = ['DKSI','DNU','DQ','DQU','DR','DUK1','DUK2']
    index_f = ['DND73','DND74','DND75','DND76','DND77','DND78','DND79','DFM1','DFK1','DFL2','DFK2','DFM1C']
    index_fstr =['f_eqs']
    v0uf_base.index = index_v0 + index_u + index_f + index_fstr
    v0uf_base = v0uf_base.loc[index_v0 + index_u + index_f].astype(float)
    
    #TODO: definitie eruit halen
    #conversion to lunar for comparison with SLS and IHO
    def get_lunarSLSIHO_fromsolar(v0uf_base):
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
    
    v0uf_baseT_lunar, v0uf_baseT_lunar_SLS, v0uf_baseT_lunar_IHO = get_lunarSLSIHO_fromsolar(v0uf_base)
    
    #v0uf_test = pd.DataFrame([[0,4,-2,0,0,0,0]],columns=['T','S','H','P','N','P1','EDN'],index=['MSQM']).T #solar hatyan convention
    #v0uf_testT_lunar, v0uf_testT_lunar_SLS, v0uf_testT_lunar_IHO = get_lunarSLSIHO_fromsolar(v0uf_test)
    
    #TODO: tabel naar textfile en samenvoegen met foreman
    #shallow water components
    shallow_eqs = {'SM': 'S2 - M2',                        #[ 6]      
                   'SNU2': 'S2 - NU2',                     #[ 8]      
                   'SN': 'S2 - N2',                        #[ 9]      
                   '2SM': '2*S2 - 2*M2',                   #[ 11]     
                   '2SMN': '2*S2 - M2 - N2',               #[ 12]     
                   'NJ1': 'N2 - J1',                       #[ 14]     
                   'NUJ1': 'NU2 - J1',                     #[ 16]     
                   'NUK1': 'NU2 - K1',                     #[ 19]     
                   'MP1': 'M2 - P1',                       #[ 22]     
                   'NO1': 'N2 - O1',                       #[ 28]     
                   'LP1': 'L2 - P1',                       #[ 30]     
                   'TK1': 'T2 - K1',                       #[ 32]     
                   'RP1': 'R2 - P1',                       #[ 37]     
                   'KP1': 'K2 - P1',                       #[ 39]     
                   'LABDAO1': 'LABDA2 - O1',               #[ 41]     
                   '2PO1': 'P1 + P1 - O1',                 #[ 43]     
                   'SO1': 'S2 - O1',                       #[ 44]     
                   'KQ1': 'K2 - Q1',                       #[ 46]     
                   '3MKS2': '3*M2 - K2 - S2',              #[ 47]     
                   '3MS2': '3*M2 - 2*S2',                  #[ 48]     
                   'OQ2': 'O1 + Q1',                       #[ 49]     
                   'MNK2': 'M2 + N2 - K2',                 #[ 50]     
                   'MNS2': 'M2 + N2 - S2',                 #[ 51]     
                   '2ML2S2': '2*M2 + L2 - 2*S2',           #[ 52]     
                   '2MS2K2': '2*M2 + S2 - 2*K2',           #[ 53]     
                   'NLK2': 'N2 + L2 - K2',                 #[ 54]     
                   '2MS2': '2*M2 - S2',                    #[ 57]     
                   'SNK2': 'S2 + N2 - K2',                 #[ 58]     
                   '2KN2S2': '2*K2 + N2 - 2*S2',           #[ 61]     
                   'OP2': 'O1 + P1',                       #[ 62]     
                   'MSK2': 'M2 + S2 - K2',                 #[ 63]     
                   'MPS2': 'M2 + P1 - S1',                 #[ 64]     
                   'MSP2': 'M2 - P1 + S1',                 #[ 66]     
                   'MKS2': 'M2 + K2 - S2',                 #[ 67]     
                   'M2(KS)2': 'M2 + 2*K2 - 2*S2',            #[ 68]     
                   '2SN(MK)2': '2*S2 + N2 - (M2 + K2)',      #[ 69]     
                   '2MN2': '2*M2 - N2',                    #[ 71]     
                   '2SK2': '2*S2 - K2',                    #[ 75]     
                   'MSN2': 'M2 + S2 - N2',                 #[ 80]     
                   'KJ2': 'K1 + J1',                       #[ 82]     
                   'MKN2': 'M2 + K2 - N2',                 #[ 83]     
                   '2KM(SN)2': '2*K2 + M2 - (S2 + N2)',      #[ 84]     
                   '2SM2': '2*S2 - M2',                    #[ 85]     
                   'SKM2': 'S2 + K2 - M2',                 #[ 86]     
                   '2SNU2': '2*S2 - NU2',                  #[ 87]     
                   '3(SM)N2': '3*S2 + N2 - 3*M2',            #[ 88]     
                   '2SN2': '2*S2 - N2',                    #[ 89]     
                   'SKN2': 'S2 + K2 - N2',                 #[ 90]     
                   'MQ3': 'M2 + Q1',                       #[ 91]     
                   'NO3': 'N2 + O1',                       #[ 92]     
                   'MO3': 'M2 + O1',                       #[ 93]     
                   '2MK3': '2*M2 - K1',                    #[ 94]     
                   '2MP3': '2*M2 - P1',                    #[ 95]     
                   'SO3': 'S2 + O1',                       #[ 97]     
                   'MK3': 'M2 + K1',                       #[ 98]     
                   '2MQ3': '2*M2 - Q1',                    #[ 99]     
                   'SP3': 'S2 + P1',                       #[100]     
                   'SK3': 'S2 + K1',                       #[101]     
                   'K3': 'K2 + K1',                        #[102]     
                   '2SO3': '2*S2 - O1',                    #[103]     
                   '4MS4': '4*M2 - 2*S2',                  #[104]     
                   '2MNS4': '2*M2 + N2 - S2',              #[105]     
                   '3MK4': '3*M2 - K2',                    #[106]     
                   'MNLK4': 'M2 + N2 + L2 - K2',           #[107]     
                   '3MS4': '3*M2 - S2',                    #[108]     
                   'MSNK4': 'M2 + S2 + N2 - K2',           #[109]     
                   'MN4': 'M2 + N2',                       #[110]     
                   '2MLS4': '2*M2 + L2 - S2',              #[111]     
                   '2MSK4': '2*M2 + S2 - K2',              #[112]     
                   'M4': '2*M2',                           #[113]     
                   '2MKS4': '2*M2 + K2 - S2',              #[114]     
                   'SN4': 'S2 + N2',                       #[115]     
                   '3MN4': '3*M2 - N2',                    #[116]     
                   '2SMK4': '2*S2 + M2 - K2',              #[117]     
                   'MS4': 'M2 + S2',                       #[118]     
                   'MK4': 'M2 + K2',                       #[119]     
                   '2SNM4': '2*S2 + N2 - M2',              #[120]     
                   '2MSN4': '2*M2 + S2 - N2',              #[121]     
                   'SL4': 'S2 + L2',                       #[122]     
                   'S4': '2*S2',                           #[123]     
                   'SK4': 'S2 + K2',                       #[124]     
                   '2SMN4': '2*S2 + M2 - N2',              #[125]     
                   '3SM4': '3*S2 - M2',                    #[126]     
                   '2SKM4': '2*S2 + K2 - M2',              #[127]     
                   'MNO5': 'M2 + N2 + O1',                 #[128]     
                   '3MK5': '3*M2 - K1',                    #[129]     
                   '3MP5': '3*M2 - P1',                    #[130]     
                   'M5': '2*M2 + M1',                      #[131]     
                   'MNK5': 'M2 + N2 + K1',                 #[132]     
                   '2MP5': '2*M2 + P1',                    #[133]     
                   '3MO5': '3*M2 - O1',                    #[134]     
                   'MSK5': 'M2 + S2 + K1',                 #[135]     
                   '3KM5': '3*K1 + M2',                    #[136]     
                   '2(MN)S6': '2*M2 + 2*N2 - S2',            #[137]     
                   '3MNS6': '3*M2 + N2 - S2',              #[138]     
                   '2NM6': '2*N2 + M2',                    #[139]     
                   '4MS6': '4*M2 - S2',                    #[140]     
                   '2MSNK6': '2*M2 + S2 + N2 - K2',        #[141]     
                   '2MN6': '2*M2 + N2',                    #[142]     
                   '2MNU6': '2*M2 + NU2',                  #[143]     
                   '3MSK6': '3*M2 + S2 - K2',              #[144]     
                   'M6': '3*M2',                           #[145]     
                   'MSN6': 'M2 + S2 + N2',                 #[146]     
                   '4MN6': '4*M2 - N2',                    #[147]     
                   'MNK6': 'M2 + N2 + K2',                 #[148]     
                   'MKNU6': 'M2 + K2 + NU2',               #[149]     
                   '2(MS)K6': '2*M2 + 2*S2 - K2',            #[150]     
                   '2MS6': '2*M2 + S2',                    #[151]     
                   '2MK6': '2*M2 + K2',                    #[152]     
                   '2SN6': '2*S2 + N2',                    #[153]     
                   '3MSN6': '3*M2 + S2 - N2',              #[154]     
                   'MKL6': 'M2 + K2 + L2',                 #[155]     
                   '2SM6': '2*S2 + M2',                    #[156]     
                   'MSK6': 'M2 + S2 + K2',                 #[157]     
                   'S6': '3*S2',                           #[158]     
                   '2MNO7': '2*M2 + N2 + O1',              #[159]     
                   '2NMK7': '2*N2 + M2 + K1',              #[160]     
                   'M7': '3*M2 + M1',                      #[161]     
                   '2MSO7': '2*M2 + S2 + O1',              #[162]     
                   'MSKO7': 'M2 + S2 + K2 + O1',           #[163]     
                   '2(MN)8': '2*M2 + 2*N2',                  #[164]     
                   '3MN8': '3*M2 + N2',                    #[165]     
                   '3MNKS8': '3*M2 + N2 + K2 - S2',        #[166]     
                   'M8': '4*M2',                           #[167]     
                   '2MSN8': '2*M2 + S2 + N2',              #[168]     
                   '2MNK8': '2*M2 + N2 + K2',              #[169]     
                   '3MS8': '3*M2 + S2',                    #[170]     
                   '3MK8': '3*M2 + K2',                    #[171]     
                   '2SNM8': 'M2 + 2*S2 + N2',              #[172]     
                   'MSNK8': 'M2 + S2 + N2 + K2',           #[173]     
                   '2(MS)8': '2*M2 + 2*S2',                  #[174]     
                   '2MSK8': '2*M2 + S2 + K2',              #[175]     
                   '3SM8': '3*S2 + M2',                    #[176]     
                   '2SMK8': '2*S2 + M2 + K2',              #[177]     
                   'S8': '3*S2 + S2',                      #[178]     
                   '2(MN)K9': '2*M2 + 2*N2 + K1',            #[179]     
                   '3MNK9': '3*M2 + N2 + K1',              #[180]     
                   '4MK9': '4*M2 + K1',                    #[181]     
                   '3MSK9': '3*M2 + S2 + K1',              #[182]     
                   '4MN10': '4*M2 + N2',                   #[183]     
                   'M10': '5*M2',                          #[184]     
                   '3MSN10': '3*M2 + S2 + N2',             #[185]     
                   '4MS10': '4*M2 + S2',                   #[186]     
                   '2(MS)N10': '2*M2 + 2*S2 + N2',           #[187]     
                   '2MNSK10': '2*M2 + N2 + S2 + K2',       #[188]     
                   '3M2S10': '3*M2 + 2*S2',                #[189]     
                   '4MSK11': '4*M2 + S2 + K1',             #[190]     
                   'M12': '5*M2 + M2',                     #[191]     
                   '4MSN12': '4*M2 + S2 + N2',             #[192]     
                   '5MS12': '5*M2 + S2',                   #[193]     
                   '3MNKS12': '3*M2 + N2 + K2 + S2',       #[194]     
                   '4M2S12': '4*M2 + 2*S2',                #[195]   
                   'N4': '2*N2',                           #extra added   
                   }
    
    shallow_eqs_pd = pd.Series(shallow_eqs)
    
    v0uf_base_forv0u = v0uf_base.loc[index_v0+index_u,:].astype(int)
    v0uf_base_forf = v0uf_base.loc[index_f,:].astype(float)
    shallow_eqs_sel = shallow_eqs_pd#[const_list_shallow_sel]
    shallow_eqs_sel_str = '\n'.join('comp_%s = %s'%(key.replace('(','').replace(')',''), val) for key, val in shallow_eqs_sel.iteritems()) #brackets are temporarily removed in order to evaluate functions
    v0uf_base_forv0u.eval(shallow_eqs_sel_str, inplace=True)
    v0uf_base_forf.eval(shallow_eqs_sel_str.replace('-','+'), inplace=True) #for f only multiplication is applied, never division
    
    #TODO: kan deels weg?
    if 1:
        v0uf_all = pd.concat([v0uf_base_forv0u,v0uf_base_forf])
    else: #option does not work at f definition yet
        v0uf_base_forfstr = v0uf_base.loc[index_fstr,:]
        shallow_eqs_sel_str_nomult = shallow_eqs_sel_str[:]
        for mult in [2,3,4,5,6]:
            for comp in ['M2','N2','K2','K1','S2']:
                shallow_eqs_sel_str_nomult = shallow_eqs_sel_str_nomult.replace('%i*%s'%(mult,comp),'+'.join(mult*[comp]))
        v0uf_base_forfstr.eval(shallow_eqs_sel_str_nomult.replace('-','+'), inplace=True) #for f only multiplication is applied, never division
        v0uf_all = pd.concat([v0uf_base_forv0u,v0uf_base_forfstr])
    
    #remove component_ prefix and restore brackets in component names (was necessary since equation cannot start with number)
    new_cols = [x.replace('comp_','') for x in v0uf_all.columns]
    for comp in ['M2(KS)2','2SN(MK)2','2KM(SN)2','3(SM)N2','2(MN)S6','2(MS)K6','2(MN)8','2(MS)8','2(MN)K9','2(MS)N10']:
        new_cols = [x.replace(comp.replace('(','').replace(')',''),comp) for x in new_cols]
    v0uf_all.columns.values[:] = new_cols
    
    #TODO: kan deels weg?
    if 1:
        v0uf_allT = v0uf_all.T
    else: #option does not work at f definition yet
        v0uf_allT_obj = v0uf_all.T
        v0uf_allT = v0uf_allT_obj.loc[:,index_v0+index_u].astype(int)
        v0uf_allT[index_fstr] = v0uf_allT_obj.loc[:,index_fstr]
    
    v0uf_allT_lunar, v0uf_allT_lunar_SLS, v0uf_allT_lunar_IHO = get_lunarSLSIHO_fromsolar(v0uf_all)
    
    return v0uf_allT


def get_hatyan_freqs(const_list, dood_date=None, sort_onfreq=True, return_allraw=False):
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


def get_hatyan_v0(const_list, dood_date):
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


def get_hatyan_u(const_list, dood_date):
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
    DOMEGA, DIKL, DC1681, DC5023, DC0365 = get_hatyan_constants(dood_date)
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


def get_hatyan_f(const_list, dood_date, xfac):
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
    DOMEGA, DIKL, DC1681, DC5023, DC0365 = get_hatyan_constants(dood_date)
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


