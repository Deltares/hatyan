# -*- coding: utf-8 -*-
"""
astrog.py contains all astro-related definitions, previously embedded in a separate program but now part of hatyan.

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
import functools
import pandas as pd
import numpy as np
import datetime as dt
import requests
import warnings
import matplotlib.pyplot as plt

from hatyan.schureman import get_schureman_freqs
    

file_path = os.path.realpath(__file__)


def astrog_culminations(tFirst,tLast,dT_fortran=False,tzone='UTC'): #TODO: add simple lon correction at end of definition, currenty calculating culmination at lon=0
    """
    Makes use of the definitions dT, astrab and astrac.
    Calculates lunar culminations, parallax and declination. By default the lunar culmination is calculated at coordinates lon=0 (Greenwich), since EHMOON is used to calculate it. Possible to add lon-correction at end of definition.

    Parameters
    ----------
    tFirst : pd.Timestamp, datetime.datetime or string ("yyyymmdd")
        Start of timeframe for output.
    tLast : pd.Timestamp, datetime.datetime or string ("yyyymmdd")
        End of timeframe for output.
    dT_fortran : boolean, optional
        Reproduce fortran difference between universal time and terrestrial time (dT). Can be True (for latest fortran reproduction) or False (international definition). The default is False.
    tzone : string/dt.timezone, optional
        Timezone to convert the output dataset to. The default is 'UTC'.

    Raises
    ------
    Exception
        Checks input times tFirst and tLast.

    Returns
    -------
    astrog_df : pandas DataFrame
        datetime:    lunar culmination at Greenwich in UTC (datetime)
        type:        type of culmination (1=lower, 2=upper)
        parallax:    lunar parallax (degrees)
        declination: lunar declination (degrees)

    """

    # check input times (datetime or string)
    [tFirst,tLast] = convert_str2datetime(datetime_in_list=[tFirst,tLast])
    
    # constants
    EHMINC       = 346.8 # increment of ephemeris hour angle of moon (deg/day)
    M2_period_hr = get_schureman_freqs(['M2']).loc['M2','period [hr]'] # interval between lunar culminations (hours)
    
    # first and last datetime in calculation (add enough margin, and an extra day for timezone differences)
    date_first = tFirst-dt.timedelta(hours=M2_period_hr+1*24)
    date_last = tLast+dt.timedelta(hours=M2_period_hr+1*24)

    # estimate culminations (time and type)
    astrabOutput = astrab(date_first,dT_fortran=dT_fortran)
    EHMOON = astrabOutput['EHMOON']
    ICUL = np.floor(EHMOON[0]%360/180).astype(int)+1 # ICUL=1: next culmination is lower culmination, ICUL=2: next culmination is upper culmination
    CULEST = pd.date_range(start=date_first+dt.timedelta(days=(180.*ICUL-EHMOON[0])/EHMINC), end=date_last, freq='%iN'%(M2_period_hr*3600*1e9)) #defined freq as M2_period in nanoseconds
    CULTYP = np.zeros(len(CULEST),dtype=int)
    CULTYP[::2] = ICUL
    CULTYP[1::2] = (ICUL%2)+1
    
    # calculate exact time of culminations
    CULTIM = astrac(CULEST, dT_fortran=dT_fortran, mode=CULTYP) #TODO: no correction with dT necessary?
    astrabOutput = astrab(CULTIM, dT_fortran=dT_fortran)
    PAR = astrabOutput['PARLAX']/3600 # conversion from arcseconds to degrees
    DEC = astrabOutput['DECMOO']
    
    # make dataframe
    astrog_df = pd.DataFrame({'datetime':CULTIM,'type':CULTYP,'parallax':PAR,'declination':DEC}) #CULTIM.round('S') decreases fortran reproduction
    astrog_df['type_str'] = astrog_df['type'].astype(str).replace('1','lowerculmination').replace('2','upperculmination')

    #set timezone, check datetime order and filter datetimerange
    astrog_df = check_crop_dataframe(astrog_df, tFirst, tLast, tzone)
    
    return astrog_df


def astrog_phases(tFirst,tLast,dT_fortran=False,tzone='UTC'):
    """
    Makes use of the definitions dT, astrab and astrac.
    Calculates lunar phases. The lunar phases are independent of coordinates.

    Parameters
    ----------
    tFirst : datetime.datetime or string ("yyyymmdd")
        Start of timeframe for output.
    tLast : datetime.datetime or string ("yyyymmdd")
        End of timeframe for output.
    dT_fortran : boolean, optional
        Reproduce fortran difference between universal time and terrestrial time (dT). Can be True (for latest fortran reproduction) or False (international definition). The default is False.
    tzone : string/dt.timezone, optional
        Timezone to convert the output dataset to. The default is 'UTC'.

    Raises
    ------
    Exception
        Checks input times tFirst and tLast.

    Returns
    -------
    astrog_df : pandas DataFrame
        datetime:  lunar phase in UTC (datetime)
        type:      type of phase (1=FQ, 2=FM, 3=LQ, 4=NM)

    """

    # check input times (datetime or string)
    [tFirst,tLast] = convert_str2datetime(datetime_in_list=[tFirst,tLast])
        
    # constants
    ELOINC = 12.2           # increment of ELONG ecliptic elongation of moon-sun (deg/day)
    FASINT = 29.530587981/4 # quarter of a lunar synodic month (days)

    # first and last datetime in calculation (add enough margin (done later), and an extra day for timezone differences)
    date_first = tFirst-dt.timedelta(days=FASINT+1)
    date_last = tLast+dt.timedelta(days=FASINT+1)

    # estimate first lunar phase (time and type), correct first date (FAEST_first to 45 deg from there)
    astrabOutput = astrab(date_first,dT_fortran=dT_fortran)
    ELONG = astrabOutput['ELONG']
    FAEST_first = date_first - pd.TimedeltaIndex((ELONG-45)/ELOINC, unit='D')

    # use the first date to create a new daterange from the correct starting time. The frequency is 29 days, 12 hours and 44 minutes, following from dood_S-dood_H
    date = pd.date_range(start=FAEST_first[0],end=date_last,freq='%iN'%(FASINT*24*3600*1e9))

    # estimate all lunar phases (time and type)
    astrabOutput = astrab(date,dT_fortran=dT_fortran)
    ELONG = astrabOutput['ELONG']
    FATYP = (np.floor(ELONG/90).astype(int)+3)%4+1 #make sure the next phase is searched for (modulus to use 'FATYP-1')
    FAEST = date - pd.TimedeltaIndex((90*FATYP-ELONG%360)/ELOINC, unit='D')

    # calculate exact time of phase, loop until date_last
    TIMDIF = pd.TimedeltaIndex(-dT(FAEST,dT_fortran=dT_fortran),unit='S')
    FATIM = astrac(FAEST,dT_fortran=dT_fortran,mode=FATYP+2) + TIMDIF

    # make dataframe
    astrog_df = pd.DataFrame({'datetime':FATIM.round('S'),'type':FATYP})
    astrog_df['type_str'] = astrog_df['type'].astype(str).replace('1','FQ').replace('2','FM').replace('3','LQ').replace('4','NM')

    #set timezone, check datetime order and filter datetimerange
    astrog_df = check_crop_dataframe(astrog_df, tFirst, tLast, tzone)

    return astrog_df


def astrog_sunriseset(tFirst,tLast,dT_fortran=False,tzone='UTC',lon=5.3876,lat=52.1562):
    """
    Makes use of the definitions dT, astrab and astrac.
    Calculates sunrise and -set at requested location.

    Parameters
    ----------
    tFirst : datetime.datetime or string ("yyyymmdd")
        Start of timeframe for output.
    tLast : datetime.datetime or string ("yyyymmdd")
        End of timeframe for output.
    dT_fortran : boolean, optional
        Reproduce fortran difference between universal time and terrestrial time (dT). Can be True (for latest fortran reproduction) or False (international definition). The default is False.
    tzone : string/dt.timezone, optional
        Timezone to convert the output dataset to. The default is 'UTC'.
    lon : float, optional
        Longitude, defined positive eastward. The default is -5.3876 (Amersfoort).
    lat : float, optional
        Latitude, defined positive northward, cannot exceed 59 (too close to poles). The default is 52.1562 (Amersfoort).

    Raises
    ------
    Exception
        Checks input times tFirst and tLast.

    Returns
    -------
    astrog_df : pandas DataFrame
        datetime: time of rise or set in UTC (datetime)
        type:     type  (1=sunrise, 2=sunset)

    """

    # check input times (datetime or string)
    [tFirst,tLast] = convert_str2datetime(datetime_in_list=[tFirst,tLast])

    # first and last datetime in calculation (add enough margin, and an extra day for timezone differences)
    date_first = tFirst - dt.timedelta(days=1)
    date_last = tLast + dt.timedelta(days=1)
    
    # --- sunrise and -set ---
    # estimate times: starting at tFirst, 0h local solar time
    DAYEST = pd.date_range(start=date_first,end=date_last,freq='%iN'%(24*3600*1e9))
    OPEST  = DAYEST + dt.timedelta(days=-lon/360.+.25) # correct for longitude and 'floor' date to 00:00 +6h
    ONEST  = DAYEST + dt.timedelta(days=-lon/360.+.75) # correct for longitude and 'floor' date to 00:00 +18h

    # calculate exact times
    TIMDIF = pd.TimedeltaIndex(-dT(OPEST,dT_fortran=dT_fortran),unit='S')
    OPTIM  = astrac(OPEST,dT_fortran=dT_fortran,mode=np.array(9),lon=lon,lat=lat) + TIMDIF
    ONTIM  = astrac(ONEST,dT_fortran=dT_fortran,mode=np.array(10),lon=lon,lat=lat) + TIMDIF

    # make dataframe
    astrog_df = pd.DataFrame({'datetime':np.concatenate((OPTIM.round('S'),ONTIM.round('S'))),'type':np.repeat([1,2],len(OPTIM))})
    astrog_df = astrog_df.sort_values('datetime').reset_index(drop=True)
    astrog_df['type_str'] = astrog_df['type'].astype(str).replace('1','sunrise').replace('2','sunset')

    #set timezone, check datetime order and filter datetimerange
    astrog_df = check_crop_dataframe(astrog_df, tFirst, tLast, tzone)
    
    return astrog_df


def astrog_moonriseset(tFirst,tLast,dT_fortran=False,tzone='UTC',lon=5.3876,lat=52.1562):
    """
    Makes use of the definitions dT, astrab and astrac.
    Calculates moonrise and -set at requested location.

    Parameters
    ----------
    tFirst : datetime.datetime or string ("yyyymmdd")
        Start of timeframe for output.
    tLast : datetime.datetime or string ("yyyymmdd")
        End of timeframe for output.
    dT_fortran : boolean, optional
        Reproduce fortran difference between universal time and terrestrial time (dT). Can be True (for latest fortran reproduction) or False (international definition). The default is False.
    tzone : string/dt.timezone, optional
        Timezone to convert the output dataset to. The default is 'UTC'.
    lon : float, optional
        Longitude, defined positive eastward. The default is -5.3876 (Amersfoort).
    lat : float, optional
        Latitude, defined positive northward, cannot exceed 59 (too close to poles). The default is 52.1562 (Amersfoort).

    Raises
    ------
    Exception
        Checks input times tFirst and tLast.

    Returns
    -------
    astrog_df : pandas DataFrame
        datetime: time of rise or set in UTC (datetime)
        type:     type  (1=moonrise, 2=moonset)

    """

    # check input times (datetime or string)
    [tFirst,tLast] = convert_str2datetime(datetime_in_list=[tFirst,tLast])

    # constants
    from hatyan.schureman import get_schureman_freqs
    M2_period_hr = get_schureman_freqs(['M2']).loc['M2','period [hr]'] # CULINT
    EHMINC = 346.8 # increment of EHMOON ephemeris hour angle of moon (deg/day) TODO: 360/(M2_period_hr*2/24) = 347.8092506037627

    # first and last datetime in calculation (add enough margin, and an extra day for timezone differences)
    date_first = tFirst-dt.timedelta(hours=M2_period_hr+1*24)
    date_last = tLast+dt.timedelta(hours=M2_period_hr+1*24)

    # --- moonrise and -set ---
    # estimate times
    astrabOutput = astrab(date_first,dT_fortran=dT_fortran,lon=lon,lat=lat)
    ALTMOO = astrabOutput['ALTMOO']
    EHMOON = astrabOutput['EHMOON']

    if ALTMOO < -(0.5667+(0.08+0.2725*astrabOutput['PARLAX'])/3600): # first phenomenon is moonrise #the PARLAX output is in arcseconds, so /3600 is conversion to degrees
        OPEST = pd.date_range(start=date_first+dt.timedelta(days=(270-EHMOON[0])/EHMINC), end=date_last, freq='%iN'%(M2_period_hr*2*3600*1e9))
        ONEST = OPEST + dt.timedelta(hours=M2_period_hr)
    else: # first phenomenon is moonset
        ONEST = pd.date_range(start=date_first+dt.timedelta(days=(90-EHMOON[0])/EHMINC), end=date_last, freq='%iN'%(M2_period_hr*2*3600*1e9))
        OPEST = ONEST + dt.timedelta(hours=M2_period_hr)

    # calculate exact times
    TIMDIF = pd.TimedeltaIndex(-dT(OPEST,dT_fortran=dT_fortran),unit='S')
    OPTIM  = astrac(OPEST,dT_fortran=dT_fortran,mode=np.array(7),lon=lon,lat=lat) + TIMDIF
    ONTIM  = astrac(ONEST,dT_fortran=dT_fortran,mode=np.array(8),lon=lon,lat=lat) + TIMDIF

    # make dataframe
    astrog_df = pd.DataFrame({'datetime':np.concatenate((OPTIM.round('S'),ONTIM.round('S'))),'type':np.repeat([1,2],len(OPTIM))})
    astrog_df = pd.DataFrame(astrog_df).sort_values('datetime').reset_index(drop=True)
    astrog_df['type_str'] = astrog_df['type'].astype(str).replace('1','moonrise').replace('2','moonset')

    #set timezone, check datetime order and filter datetimerange
    astrog_df = check_crop_dataframe(astrog_df, tFirst, tLast, tzone)

    return astrog_df


def astrog_anomalies(tFirst,tLast,dT_fortran=False,tzone='UTC'):
    """
    Makes use of the definitions dT, astrab and astrac.
    Calculates lunar anomalies. The lunar anomalies are independent of coordinates.

    Parameters
    ----------
    tFirst : datetime.datetime or string ("yyyymmdd")
        Start of timeframe for output.
    tLast : datetime.datetime or string ("yyyymmdd")
        End of timeframe for output.
    dT_fortran : boolean, optional
        Reproduce fortran difference between universal time and terrestrial time (dT). Can be True (for latest fortran reproduction) or False (international definition). The default is False.
    tzone : string/dt.timezone, optional
        Timezone to convert the output dataset to. The default is 'UTC'.

    Raises
    ------
    Exception
        Checks input times tFirst and tLast.

    Returns
    -------
    astrog_df : pandas DataFrame
        datetime:   lunar anomaly in UTC (datetime)
        type:       type of anomaly (1=perigeum, 2=apogeum)

    """

    # check input times (datetime or string)
    [tFirst,tLast] = convert_str2datetime(datetime_in_list=[tFirst,tLast])

    # constants
    ANMINC = 13.06       # increment of ANM anomaly of moon (deg/day)
    ANOINT = 27.554551/2 # half of a lunar anomalistic month (days)

    # first and last datetime in calculation (add enough margin, and an extra day for timezone differences)
    date_first = tFirst-dt.timedelta(days=ANOINT+1)
    date_last = tLast+dt.timedelta(days=ANOINT+1)

    # estimate first lunar anomaly (time and type)
    astrabOutput = astrab(date_first,dT_fortran=dT_fortran)
    DPAXDT = astrabOutput['DPAXDT']
    ANM = astrabOutput['ANM']
    if DPAXDT>0: # ANOTYP=1: perigeum first
        ref_deg = 360    
        IANO = 1
    elif DPAXDT<=0: # ANOTYP=2: apogeum first
        ref_deg = 180    
        IANO = 2
    ANOEST = pd.date_range(start=date_first+dt.timedelta(days=(ref_deg-ANM[0])/ANMINC),end=date_last,freq='%iN'%(ANOINT*24*3600*1e9))
    ANOTYP = np.zeros(len(ANOEST),dtype=int)
    ANOTYP[::2] = IANO
    ANOTYP[1::2] = (IANO%2)+1

    # calculate exact times
    TIMDIF = pd.TimedeltaIndex(-dT(ANOEST,dT_fortran=dT_fortran),unit='S')
    ANOTIM = astrac(ANOEST,dT_fortran=dT_fortran,mode=ANOTYP+14) + TIMDIF

    # make dataframe
    astrog_df = pd.DataFrame({'datetime':ANOTIM.round('S'),'type':ANOTYP})
    astrog_df['type_str'] = astrog_df['type'].astype(str).replace('1','perigeum').replace('2','apogeum')

    #set timezone, check datetime order and filter datetimerange
    astrog_df = check_crop_dataframe(astrog_df, tFirst, tLast, tzone)
    
    return astrog_df


def astrog_seasons(tFirst,tLast,dT_fortran=False,tzone='UTC'):
    """
    Makes use of the definitions dT, astrab and astrac.
    Calculates astronomical seasons. The seasons are independent of coordinates.

    Parameters
    ----------
    tFirst : datetime.datetime or string ("yyyymmdd")
        Start of timeframe for output.
    tLast : datetime.datetime or string ("yyyymmdd")
        End of timeframe for output.
    dT_fortran : boolean, optional
        Reproduce fortran difference between universal time and terrestrial time (dT). Can be True (for latest fortran reproduction) or False (international definition). The default is False.
    tzone : string/dt.timezone, optional
        Timezone to convert the output dataset to. The default is 'UTC'.

    Raises
    ------
    Exception
        Checks input times tFirst and tLast.

    Returns
    -------
    astrog_df : pandas DataFrame
        datetime:   start of astronomical season in UTC (datetime)
        type:       type of astronomical season (1=spring, 2=summer, 3=autumn, 4=winter)

    """

    # check input times (datetime or string)
    [tFirst,tLast] = convert_str2datetime(datetime_in_list=[tFirst,tLast])

    # estimate start of seasons (time and type)
    SEIEST = pd.date_range(start=dt.datetime(tFirst.year,int(np.ceil(tFirst.month/3)*3),1),end=tLast+dt.timedelta(days=1),freq='%iMS'%(3))+dt.timedelta(days=20)
    SEITYP = np.floor(SEIEST.month/3).astype(int)

    # calculate exact times, loop until tLast
    TIMDIF = pd.TimedeltaIndex(-dT(SEIEST,dT_fortran=dT_fortran),unit='S')
    SEITIM = astrac(SEIEST,dT_fortran=dT_fortran,mode=SEITYP+10) + TIMDIF

    # make dataframe
    astrog_df = pd.DataFrame({'datetime':SEITIM.round('S'),'type':SEITYP})
    astrog_df['type_str'] = astrog_df['type'].astype(str).replace('1','spring').replace('2','summer').replace('3','autumn').replace('4','winter')
    
    #set timezone, check datetime order and filter datetimerange
    astrog_df = check_crop_dataframe(astrog_df, tFirst, tLast, tzone)
    
    return astrog_df


def astrab(date,dT_fortran=False,lon=5.3876,lat=52.1562):
    """
    Python version of astrab.f in FORTRAN 77
    Calculates 18 astronomical parameters at requested time.

    Parameters
    ----------
    date : datetime.datetime or pandas.DatetimeIndex
        Requested time for calculation.
    dT_fortran : boolean, optional
        Reproduce fortran difference between universal time and terrestrial time (dT). Can be True (for latest fortran reproduction) or False (international definition). The default is False.
    lon : float, optional
        Longitude for altitudes, defined positive eastward. The default is 5.3876 (Amersfoort).
    lat : float, optional
        Latitude for altitudes, defined positive northward. The default is 52.1562 (Amersfoort).

    Raises
    ------
    Exception
        Checks if input is valid.

    Returns
    -------
    astrabOutput : dictionary
        18 astronomical variables:
            EHMOON, lunar ephemeris hour angle (degrees between +90 and +450)
            DECMOO, lunar declination (degrees)
            PARLAX, lunar horizontal parallax (arcseconds)
            DPAXDT, time derivative of parallax (arcseconds/day)
            ALTMOO, lunar altitude (degrees, negative below horizon)
            ELONG,  ecliptic elongation moon-sun (degrees, between +45 and +405)
            ALTSUN, solar altitude (degrees, negative under horizon)
            LONSUN, solar longitude (degrees, between 45 and 405)
            EQELON, equatorial elongaton moon-sun (degrees, between 0 and 360)
            DECSUN, solar declination (degrees)
            DISSUN, relative distance earth-sun (astronomical units)
            EHARI,  ephemeris hour angle vernal equinox (degrees, between 0 and 360)
            RASUN,  solar right ascension (degrees, between 0 and 360)
            EHSUN,  solar ephemeris hour angle (degrees, between 0 and 360)
            LONMOO, lunar longitude (degrees, between 0 and 360)
            LATMOO, lunar latitude (degrees)
            RAMOON, lunar right ascension (degrees, between 0 and 360)
            ANM,    mean lunar anomaly (degrees, between 0 and 360)
    
    """
    
    dT_TT_days = dT(date,dT_fortran=dT_fortran)/3600/24 #Difference between terrestrial and universal time in days.
    
    # check input
    if isinstance(date, dt.datetime):
        date = pd.DatetimeIndex([date])
    elif not isinstance(date, pd.DatetimeIndex):
        raise Exception('Input variable date should be datetime or pd.DateTimeIndex')
    
    if np.abs(lon)>180:
        raise Exception('Input variable longitude larger than 180deg')
    if np.abs(lat)>90:
        raise Exception('Input variable latitude larger than 90deg')
    
    if (date<dt.datetime(1900,1,1)).any() or (date>dt.datetime(2091,1,1)).any():
        import matplotlib.pyplot as plt
        fig,ax = plt.subplots()
        ax.plot(date)
        raise Exception('Requested time out of range (1900-2091)')

    # constants - general
    EPOCH  = dt.datetime(1899, 12, 31, 12, 0, 0) # 1900.0 # -12h shift because julian date 0 is at noon?
    # The average orbital elements of the celestial bodies are calculated for the epoch 1900.0.
    # The values are corrected for the year 1990. Intitial values are from the vernal equinox. #TODO: should correction be recalculated?
    
    #TODO: some constants correspond with schureman.get_schureman_constants(), merge constants and put in dictionary?
    # constants - sun
    LABOS  = 4.8816237     # longitude sun (rad)
    NSUN   = 0.01720279153 # increment longitude sun (rad/day)
    PEROS  = 4.9082230     # longitude perigeum sun (rad)
    BETSUN = 8.2188819E-7  # increment longitude perigeum sun (rad/day)

    # constants - moon
    LABOM  = 4.7199860     # longitude moon (rad)
    NMOON  = 0.2299715020  # increment longitude moon (rad/day)
    PEROM  = 5.8352992     # longitude perigeum moon (rad)
    BETMOO = 0.0019443591  # increment longitude perigeum moon (rad/day)
    NODOM  = 4.523572      # longitude lunar orbital node lunar (rad)
    GAMMOO = -9.2421851E-4  # increment longitude lunar orbital node lunar (rad/day)
    INMOON = 0.089804108   # inclination lunar orbit (rad) >> DIKL
    PARMEA = 3422.608      # mean horizontal lunar parallax (arcseconds) >> DAGC=np.deg2rad(PARMEA/3600)

    # constants - planets
    VENTZE = 1.10079       # elongation Venus-Earth (rad)
    VENTIN = 0.010760328   # increment elongation Venus-Earth (rad/day)
    TJUPZE = 3.86848       # elongation Earth-Jupiter (rad)
    TJUPIN = 0.015751909   # increment elongation Earth-Jupiter (rad/day)
    TMARZE = 2.89636       # elongation Earth-Mars (rad)
    TMARIN = 0.008056024   # increment elongation Earth-Mars (rad/day)
    TSATZE = 3.37079       # elongation Earth-Saturnus (rad)
    TSATIN = 0.016618143   # increment elongation Earth-Saturnus (rad/day)

    # constants - ecliptic
    OBZERO = 0.40931977    # inclination of ecliptic (rad) >> DOMEGA
    OBINC = -6.21937E-9    # increment inclination of ecliptic (rad/day)

    # constants - vernal equinox
    ARZERO = 4.881523      # ephemeris hour angle of vernal equinox (rad)
    NARIES = 6.30038809878 # increment ephemeris hour angle of vernal equinox (rad/day)

    #TODO: tabellen naar file (of aparte definitie)
    # constants - lunar orbital disturbances
    # selected from Brown's Tables of the Motion of the Moon (1909)
    # storingen in longitude
    distP = np.array([[      0,        0,     1,       1,         1,         1,       1,       0,        0,        0,      0,        0,      2,       2,        2,       2,      1,        1,        1,      1,      1,       1,      1,      0,      0,      0,        0,       0,      1,      1,     1,      0,     3,      3,       3,      3,      2,      2,      2,     2,     2,      2,      1,      1,     1,     1,     1,       1,      1,      1,     1,      0,      0,     2,     2,     1,      1,     4,     4,      2,      2],
                      [      0,        0,     0,       0,         0,         0,       0,       1,        1,        1,      1,        0,      0,       0,        0,       0,      1,        1,        1,      1,     -1,      -1,     -1,      2,      2,      0,        0,       0,      0,      0,     0,      1,     0,      0,       0,      0,      1,      1,      1,    -1,    -1,     -1,      2,      2,    -2,    -2,     0,       0,      0,      0,     0,      1,      1,     0,     0,     1,     -1,     0,     0,      0,      0],
                      [      0,        0,     0,       0,         0,         0,       0,       0,        0,        0,      0,        0,      0,       0,        0,       0,      0,        0,        0,      0,      0,       0,      0,      0,      0,      2,        2,       2,      0,      0,     0,      0,     0,      0,       0,      0,      0,      0,      0,     0,     0,      0,      0,      0,     0,     0,     2,       2,     -2,     -2,    -2,      2,     -2,     0,     0,     0,      0,     0,     0,      2,     -2],
                      [      4,        2,     4,       2,         0,        -2,      -4,       2,        0,       -2,     -4,        1,      2,       0,       -2,      -4,      2,        0,       -2,     -4,      2,       0,     -2,      0,     -2,      2,        0,      -2,      1,     -1,    -3,      1,     2,      0,      -2,     -4,      0,     -2,     -4,     2,     0,     -2,      0,     -2,     0,    -2,     2,       0,      2,      0,    -2,     -2,      2,    -1,    -3,     1,     -1,     0,    -2,      0,      0],
                      [ 13.902, 2369.902, 1.979, 191.953, 22639.500, -4586.426, -38.428, -24.420, -666.608, -164.773, -1.877, -125.154, 14.387, 769.016, -211.656, -30.773, -2.921, -109.420, -205.499, -4.391, 14.577, 147.361, 28.475, -7.486, -8.096, -5.741, -411.608, -55.173, -8.466, 18.609, 3.215, 18.023, 1.060, 36.124, -13.193, -1.187, -7.649, -8.627, -2.740, 1.181, 9.703, -2.494, -1.167, -7.412, 2.580, 2.533, -0.992, -45.099, -6.382, 39.532, 9.366, -2.152, -1.440, 1.750, 1.225, 1.267, -1.089, 1.938, -0.952, -3.996, -1.298]])
    # storingen in latitude
    distC = np.array([[     0,     2,     3,     0,     0,     1,     1,    -1,    -1],
                      [     0,     0,     0,     1,     2,     1,     1,     1,     1],
                      [     0,     0,     0,     0,     0,     0,     0,     0,     0],
                      [     1,    -2,     0,     0,    -2,     2,    -2,     0,    -2],
                      [-0.725, 5.679,-1.300,-1.302,-0.740, 0.787, 2.056, 0.679,-1.540]])
    
    distS = np.array([[       0,       0,       0,       1,       1,       1,       1,    1,      1,     2,       2,       2,       2,     3,      3,      0,     0,       0,       0,       0,       0,       1,       1,       1,       1,      -1,      -1,     -1,      2,    2,     2,     1,      0,     1,     -1],
                      [       0,       0,       0,       0,       0,       0,       0,    0,      0,     0,       0,       0,       0,     0,      0,      1,     1,       1,       1,       1,       2,       1,       1,       1,       1,       1,       1,      1,      1,    1,    -1,     2,      0,     0,      0],
                      [       0,       0,       0,       0,       0,       0,       0,    0,      0,     0,       0,       0,       0,     0,      0,      0,     0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,      0,      0,    0,     0,     0,      2,     2,      2],
                      [       1,       2,       4,       4,       2,       0,      -2,   -3,     -4,     2,       0,      -2,      -4,     0,     -2,      2,     1,       0,      -2,      -4,      -2,       2,       0,      -2,      -4,       2,       0,     -2,      0,   -2,     0,    -2,     -2,    -2,      0],
                      [ -112.79, 2373.36,   14.06,    6.98,  192.72,22609.07,-4578.13, 5.44, -38.64, 14.78,  767.96, -152.53,  -34.07, 50.64, -16.40, -25.10, 17.93, -126.98, -165.06,   -6.46,  -16.35,  -11.75, -115.18, -182.36,   -9.66,  -23.59, -138.76, -31.70, -10.56,-7.59, 11.67, -6.12, -52.14, -9.52, -85.13]])
    
    distN = np.array([[       0,      0,      1,      1,     -1,      -1,      -2,     -2,       0,      0],
                      [       0,      0,      0,      0,      0,       0,       0,      0,       1,     -1],
                      [       1,      1,      1,      1,      1,       1,       1,      1,       1,      1],
                      [      -2,     -4,     -2,     -4,      0,      -2,       0,     -2,      -2,     -2],
                      [-526.069, -3.352, 44.297, -6.000, 20.599, -30.598, -24.649, -2.000, -22.571, 10.985]])
    # storingen in parallax
    distR = np.array([[      0,       0,      1,        1,       1,      1,       0,       0,      0,       0,      2,       2,       2,      2,       1,      1,      1,      1,       1,       0,       1,      0,      3,       3,       2,      2,       1],
                      [      0,       0,      0,        0,       0,      0,       1,       1,      1,       0,      0,       0,       0,      0,       1,      1,     -1,     -1,      -1,       0,       0,      1,      0,       0,       1,     -1,       0],
                      [      0,       0,      0,        0,       0,      0,       0,       0,      0,       0,      0,       0,       0,      0,       0,      0,      0,      0,       0,       2,       0,      0,      0,       0,       0,      0,      -2],
                      [      4,       2,      2,        0,      -2,     -4,       2,       0,     -2,       1,      2,       0,      -2,     -4,       0,     -2,      2,      0,      -2,      -2,       1,      1,      0,      -2,       0,      0,       0],
                      [ 0.2607, 28.2333, 3.0861, 186.5398, 34.3117, 0.6008, -0.2993, -0.3988, 1.9135, -0.9781, 0.2833, 10.1657, -0.3039, 0.3722, -0.9469, 1.4404, 0.2297, 1.1502, -0.2252, -0.1052, -0.1093, 0.1494, 0.6215, -0.1187, -0.1038, 0.1268, -0.7136]])

    # process input data
    TIME  = np.array((date-EPOCH).total_seconds()/86400) # days since 1900
    RLONG = -np.deg2rad(lon)
    RLATI = np.deg2rad(lat)

    # mean orbital elements
    # of sun, including additive correction
    LABSUN = (LABOS - np.deg2rad(6.22/3600) + NSUN*TIME) % (2*np.pi)
    PERSUN = PEROS+BETSUN*TIME
    # of moon
    LABMOO = (LABOM+NMOON *TIME) % (2*np.pi)
    PERMOO = (PEROM+BETMOO*TIME) % (2*np.pi)
    NODMOO = NODOM+GAMMOO*TIME
    # additive corrections for moon
    LABMOO = LABMOO+np.deg2rad(( 0.79 +14.27 *np.sin(6.0486+6.34913E-5*TIME)+7.261*np.sin(NODMOO))/3600)
    PERMOO = PERMOO+np.deg2rad((-1.966- 2.076*np.sin(NODMOO))/3600)
    NODMOO = NODMOO+np.deg2rad(( 0.59 +95.96 *np.sin(NODMOO)+15.58*np.sin(NODMOO+4.764400))/3600)
    CINMOO = np.deg2rad((-8.636*np.cos(NODMOO)-1.396*np.cos(NODMOO+4.764400))/3600)

    # mean heliocentric elongation Venus-Earth, Earth-Jupiter, Earth-Mars and Earth-Saturnus
    VENTER = VENTZE+TIME*VENTIN
    TERJUP = TJUPZE+TIME*TJUPIN
    TERMAR = TMARZE+TIME*TMARIN
    TERSAT = TSATZE+TIME*TSATIN

    # Nutation in longitude
    CNULON = -17.248*np.sin(NODMOO)-1.273*np.sin(2.*LABSUN)
    # Nutation in inclination of ecliptic
    CNUTOB =   9.21 *np.cos(NODMOO)+0.55 *np.cos(2.*LABSUN)

    # inclination of ecliptic, sine and cosine of inclination
    OBLIQ = OBZERO+TIME*OBINC+np.deg2rad(CNUTOB/3600)
    sin_OBLIQ = np.sin(OBLIQ)
    cos_OBLIQ = np.cos(OBLIQ)

    # parameters of sun and moon for disturbance equations
    ANM = (LABMOO-PERMOO) # output value 18: mean lunar anomaly (rad)
    ANM_arr = ANM[:,np.newaxis] # output value 18: mean lunar anomaly (rad)
    ANS = (LABSUN-PERSUN)
    ANS_arr = ANS[:,np.newaxis]
    FNO = (LABMOO-NODMOO)
    FNO_arr = FNO[:,np.newaxis]
    ELO = (LABMOO-LABSUN)
    ELO_arr = ELO[:,np.newaxis]
    DAM = NMOON-BETMOO
    DAS = NSUN-BETSUN
    DFN = NMOON-GAMMOO
    DEL = NMOON-NSUN

    # disturbance equations
    # disturbances in longitude sun
    # equation of center
    CENTR=(6910.10 -17.33*TIME/36525)*np.sin(ANS)+72.01*np.sin(2.*ANS) + 1.05*np.sin(3.*ANS)
    # planetary disturbances
    PLANET = (4.838 * np.cos(  VENTER        + 1.5708) +
              5.526 * np.cos(2*VENTER        + 1.5723) +
              0.666 * np.cos(3*VENTER        + 4.7195) +
              2.497 * np.cos(2*VENTER  - ANS + 4.4986) +
              1.559 * np.cos(3*VENTER  - ANS + 1.3607) +
              1.024 * np.cos(3*VENTER -2*ANS + 0.8875) +
              7.208 * np.cos(  TERJUP        + 1.5898) +
              2.731 * np.cos(2*TERJUP        + 4.7168) +
              2.600 * np.cos(  TERJUP  - ANS + 3.0503) +
              1.610 * np.cos(2*TERJUP  - ANS + 5.1068) +
              0.556 * np.cos(3*TERJUP  - ANS + 3.0946) +
              2.043 * np.cos(2*TERMAR        + 1.5660) +
              1.770 * np.cos(2*TERMAR  - ANS + 5.3454) +
              0.585 * np.cos(4*TERMAR -2*ANS + 3.2432) +
              0.500 * np.cos(4*TERMAR  - ANS + 5.5317) +
              0.425 * np.cos(3*TERMAR  - ANS + 5.5449) +
              0.419 * np.cos(  TERSAT        + 1.5767) +
              0.320 * np.cos(  TERSAT  - ANS + 4.5242) )

    # geometric disturbance by the moon
    GEOM = 6.454*np.sin(ELO)-0.424*np.sin(ELO-ANM)
    
    # aberration (correction for optical path)
    ABER = -(20.496+.344*np.cos(ANS))
    
    # longitude sun with all disturbances
    LONSUN = LABSUN + np.deg2rad((CENTR+PLANET+GEOM+CNULON+ABER)/3600) # output value 8: solar longitude (rad)
    
    # relative distance from sun
    DISSUN = 1.000140 - 0.016712*np.cos(ANS) - 0.000140*np.cos(2.*ANS) # output value 11: relative distance earth-sun (astronomical units)
    
    # longitude moon with all disturbances
    CLONM = (np.sin(distP[0,:]*ANM_arr + distP[1,:]*ANS_arr + distP[2,:]*FNO_arr + distP[3,:]*ELO_arr) * distP[4,:]).sum(axis=1)
    LONMOO = LABMOO+np.deg2rad((CNULON+CLONM)/3600) # output value 15: lunar longitude (rad)
    
    # latitude moon with all disturbances
    CLM = (np.cos(distC[0,:]*ANM_arr + distC[1,:]*ANS_arr + distC[2,:]*FNO_arr + distC[3,:]*ELO_arr) * distC[4,:]).sum(axis=1)
    SLM = (np.sin(distS[0,:]*ANM_arr + distS[1,:]*ANS_arr + distS[2,:]*FNO_arr + distS[3,:]*ELO_arr) * distS[4,:]).sum(axis=1)
    SF=FNO + np.deg2rad(SLM/3600)
    NLM = (np.sin(distN[0,:]*ANM_arr + distN[1,:]*ANS_arr + distN[2,:]*FNO_arr + distN[3,:]*ELO_arr) * distN[4,:]).sum(axis=1)
    LATMOO = ((18519.7+CLM)*np.sin(SF) - 6.241*np.sin(3.*SF) + NLM)*np.deg2rad((1.+CINMOO/INMOON)/3600) # output value 16: lunar latitude (rad)

    # lunar parallax. output value 3: lunar horizontal parallax (arcseconds)
    PARLAX = PARMEA + (np.cos(distR[0,:]*ANM_arr + distR[1,:]*ANS_arr + distR[2,:]*FNO_arr + distR[3,:]*ELO_arr) * distR[4,:]).sum(axis=1)
    
    #derrivative of lunar parallax. output value 4: time derivative of parallax (arcseconds/day)
    DPAXDT_sub1 = np.sin(distR[0,:]*ANM_arr + distR[1,:]*ANS_arr + distR[2,:]*FNO_arr + distR[3,:]*ELO_arr)
    DPAXDT_sub2 = distR[0,:]*DAM + distR[1,:]*DAS + distR[2,:]*DFN + distR[3,:]*DEL
    DPAXDT = -(DPAXDT_sub1 * DPAXDT_sub2 * distR[4,:]).sum(axis=1)
    
    # ecliptic elongation moon-sun
    ELONG = LONMOO-LONSUN # output value 6: ecliptic elongation moon-sun (rad)

    # transformation to equatorial coordinates
    sin_LONSUN = np.sin(LONSUN)
    cos_LONSUN = np.cos(LONSUN)
    RASUN = np.arctan2(sin_LONSUN*cos_OBLIQ, cos_LONSUN)           # output value 13: solar right ascension (rad)
    TEMP3 = sin_LONSUN*sin_OBLIQ
    DECSUN = np.arcsin(TEMP3)                                      # output value 10: solar declination (rad)
    sin_LONMOO = np.sin(LONMOO)
    cos_LONMOO = np.cos(LONMOO)
    sin_LATMOO = np.sin(LATMOO)
    cos_LATMOO = np.cos(LATMOO)
    RAMOON = np.arctan2(cos_LATMOO*sin_LONMOO*cos_OBLIQ - sin_LATMOO*sin_OBLIQ, cos_LATMOO*cos_LONMOO) # output value 17: lunar right ascension (rad)
    TEMP8 = cos_LATMOO*sin_LONMOO*sin_OBLIQ + sin_LATMOO*cos_OBLIQ
    DECMOO = np.arcsin(TEMP8)                                      # output value 2: lunar declination (rad)
    EQELON = RAMOON-RASUN                                          # output value 9: equatorial elongaton moon-sun (rad)

    # uurhoeken
    EHARI = (ARZERO+NARIES*TIME+cos_OBLIQ*np.deg2rad(CNULON/3600)) % (2*np.pi) # output value 12: ephemeris hour angle of vernal equinox (rad)
    LHARI = EHARI-dT_TT_days*NARIES-RLONG                         # local hour angle of vernal equinox (rad)
    EHSUN = EHARI-RASUN                                           # output value 14: solar ephemeris hour angle (rad)
    EHMOON = EHARI-RAMOON                                         # output value  1: lunar ephemeris hour angle (rad)
    LHSUN = LHARI-RASUN                                           # local solar hour angle (rad)
    LHMOON = LHARI-RAMOON                                         # local lunar hour angle (rad)

    # transformation to local coordinates (astrab.f regel 558)
    # LOCALE HOOGTE VAN DE ZON (ZONDER PARALLAX-CORRECTIE)
    ARGUM = TEMP3*np.sin(RLATI)+np.cos(DECSUN)*np.cos(RLATI)*np.cos(LHSUN)
    ALTSUN = np.arcsin(ARGUM)#TODO: np.nan_to_num(np.arcsin(ARGUM),nan=np.copysign(np.pi/2,ARGUM))  # output value 7: solar altitude (rad). TODO: np.nan_to_num seems not necesary, filling in condition from if-statement
    # LOCALE HOOGTE VAN DE MAAN (MET PARALLAX-CORRECTIE)
    ARGUM = TEMP8*np.sin(RLATI)+np.cos(DECMOO)*np.cos(RLATI)*np.cos(LHMOON)
    ALTMOO = np.arcsin(ARGUM)
    ALTMOO = ALTMOO-np.cos(ALTMOO)*np.deg2rad(PARLAX/3600)# TODO: np.nan_to_num(ALTMOO-np.cos(ALTMOO)*np.deg2rad(PARLAX/3600),nan=np.copysign(np.pi/2,ARGUM)) # output value 5: lunar altitude (rad). TODO: np.nan_to_num seems not necesary, filling in condition from if-statement
    
    # summarize in dataframe and convert output to degrees
    astrabOutput = {'EHMOON': (np.rad2deg(EHMOON)-90) % 360 + 90,
                    'DECMOO': np.rad2deg(DECMOO),
                    'PARLAX': PARLAX,
                    'DPAXDT': DPAXDT,
                    'ALTMOO': np.rad2deg(ALTMOO),
                    'ELONG' : (np.rad2deg(ELONG)-45) % 360 + 45,
                    'ALTSUN': np.rad2deg(ALTSUN),
                    'LONSUN': (np.rad2deg(LONSUN)-45) % 360 + 45,
                    'EQELON': np.rad2deg(EQELON) % 360,
                    'DECSUN': np.rad2deg(DECSUN),
                    'DISSUN': DISSUN,
                    'EHARI' : EHARI,
                    'RASUN' : np.rad2deg(RASUN) % 360,
                    'EHSUN' : np.rad2deg(EHSUN) % 360,
                    'LONMOO': np.rad2deg(LONMOO) % 360,
                    'LATMOO': np.rad2deg(LATMOO),
                    'RAMOON': np.rad2deg(RAMOON) % 360,
                    'ANM'   : np.rad2deg(ANM) % 360}

    return astrabOutput


def astrac(timeEst,mode,dT_fortran=False,lon=5.3876,lat=52.1562):
    """
    Python version of astrac.f in FORTRAN 77.
    Calculates exact time of requested astronomical phenomenon.

    Parameters
    ----------
    timeEst : datetime.datetime or pandas.DatetimeIndex
        Estimated time for iteration.
    mode : numpy.array of integer(s)
        Requested phenomenon:
            1:  lunar lower culmination (EHMOON=180 deg.)
            2:  lunar upper culmination (EHMOON=360 deg.)
            3:  lunar first quarter (ELONG=90 deg.)
            4:  full moon (ELONG=180 deg.)
            5:  lunar last quarter (ELONG=270 deg.)
            6:  new moon (ELONG=360 deg.)
            7:  moonrise (ALTMOO=-34 BOOGMIN-SEMIDIAM., ascending)
            8:  moonset (ALTMOO=-34 BOOGMIN-SEMIDIAM., descending)
            9:  sunrise (ALTSUN=-50 arcseconds, ascending)
            10: sunset (ALTSUN=-50 arcseconds, descending)
            11: vernal equinox (LONSUN=360 deg.)
            12: summer solstice (LONSUN=90 deg.)
            13: autumnal equinox (LONSUN=180 deg.)
            14: winter solstice (LONSUN=270 deg.)
            15: perigeum (DPAXDT=0, descending)
            16: apogeum (DPAXDT=0, ascending)
    dT_fortran : boolean, optional
        Reproduce fortran difference between universal time and terrestrial time (dT). Can be True (for latest fortran reproduction) or False (international definition). The default is False.
    lon : float, optional
        Longitude for rise and set, defined positive eastward. The default is 5.3876 (Amersfoort).
    lat : float, optional
        Latitude for rise and set, defined positive northward, cannot exceed 59 for modes 7 to 10 (too close to poles). The default is 52.1562 (Amersfoort).

    Raises
    ------
    Exception
        Checks if latitude is not too close to poles.

    Returns
    -------
    TIMOUT : datetime
        Exact time after iteration.

    """

    if isinstance(timeEst, pd.DatetimeIndex):
        pass
    elif isinstance(timeEst, dt.datetime):
        timeEst = pd.DatetimeIndex([timeEst])
    else:
        raise Exception('Input variable date should be datetime or pd.DateTimeIndex')

    # constants - iteration targets
    itertargets_pd = pd.DataFrame({'IPAR':  ['EHMOON']*2 +['ELONG']*4 +          ['ALTMOO']*2 +  ['ALTSUN']*2 +   ['LONSUN']*4 +          ['DPAXDT']*2 ,
                                   'ANGLE':[180,   360,   90,   180,  270,  360, -0.5667,-0.5667,-0.8333,-0.8333, 360,  90,   180,  270,  0,     0],
                                   #'ANGLE': [180,   360,   90,   180,  270,  360, -34/60, -34/60, -50/60, -50/60,  360,  90,   180,  270,  0,     0], #TODO: probably usefull to add more accuracy, but astrac testbank has to be redefined
                                   'CRITER':[1e-3,  1e-3,  1e-4, 1e-4, 1e-4, 1e-4, 2e-4,   2e-4,   2e-4,   2e-4,   1e-5, 1e-5, 1e-5, 1e-5, 2e-3,  2e-3],
                                   'RAT':   [346.8, 346.8, 12.2, 12.2, 12.2, 12.2, 346.8, -346.8,  360,   -360,    1,    1,    1,    1,   -10.08, 10.08]})
    itertargets_pd.index = range(1,len(itertargets_pd)+1)
    
    ANG = itertargets_pd.loc[mode,'ANGLE'] # required value after iteration
    CRIT = itertargets_pd.loc[mode,'CRITER'] # allowed difference between ANG and iteration result
    RATE = itertargets_pd.loc[mode,'RAT'] # estimated change per day for iteration. [deg/day]?
    IPAR_all = itertargets_pd.loc[mode,'IPAR'] # define astrab output parameter corresponding to requested mode
    if len(np.unique(IPAR_all))!=1:
        raise Exception('incorrectly mixed modes requested, results in more than one IPAR')
    IPAR = np.unique(IPAR_all)[0]
    if IPAR in ['ALTMOO','ALTSUN']: # correct RATE in case of rise and set for latitude
        if np.abs(lat)>59:
            raise Exception('Latitude to close to poles (>59deg), cannot take polar days and nights into account')
        RATE=RATE*np.cos(np.deg2rad(lat))

    # calculate value at start of iteration
    ITER=1
    TNEW = timeEst
    astrabOutput = astrab(TNEW,dT_fortran=dT_fortran,lon=lon,lat=lat)
    PNEW = astrabOutput[IPAR]
    # iterate until criterium is reached or max 20 times (including catch for sunrise and moonrise)
    while (abs(ANG-PNEW) > CRIT).any():# and ITER <=20:
        TOLD = TNEW.copy()
        POLD = PNEW.copy()
        if IPAR=='ALTMOO': # correction for semidiameter moon #TODO, might not be necessary to put in loop since PARLAX is not that variable over iterations
            ANG = itertargets_pd.loc[mode,'ANGLE']-(0.08+0.2725*astrabOutput['PARLAX'])/3600. # ANGLE in degrees and PARLAX in arcseconds (/3600 gives degrees)
        addtime = pd.TimedeltaIndex(np.nan_to_num((ANG-POLD)/RATE,0),unit='D') #nan_to_num to make sure no NaT output in next iteration
        #print(f'{ITER} timediff in hours:\n%s'%(np.array(addtime.total_seconds()/3600).round(2)))
        if IPAR in ['ALTMOO','ALTSUN'] and (np.abs(np.array(addtime.total_seconds()/3600)) > 24).any(): #catch for ALTMOO and ALTSUN to let the iteration process stop before it escalates
            raise Exception(f'Iteration step resulted in time changes larger than 24 hours (max {np.abs(addtime.total_seconds()).max()/3600:.2f} hours), try using a lower latitude')
        TNEW = TOLD + addtime
        ITER = ITER+1
        astrabOutput = astrab(TNEW,dT_fortran=dT_fortran,lon=lon,lat=lat)
        PNEW = astrabOutput[IPAR]
        RATE = np.array((PNEW-POLD)/((TNEW-TOLD).total_seconds()/86400))
        if ITER>20: #iteration catch from fortran, mibht not be necessary anymore
            raise Exception('Stopped after %s iterations, datetime=\n%s' %(ITER,TNEW))
    TIMOUT = TNEW#.round('S') # rounding everything to seconds reduces the accuracy of the reporduction of FORTRAN culmination times

    return TIMOUT


@functools.lru_cache() #caching this prevents retrieval from internet every second
def get_leapsecondslist_fromurlorfile():
    """
    

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    dict_leap_seconds : TYPE
        DESCRIPTION.
    expirydate : TYPE
        DESCRIPTION.

    """
    
    refdate = dt.datetime(1900,1,1)

    url_leap_seconds_list = 'https://raw.githubusercontent.com/eggert/tz/main/leap-seconds.list' #previously, https://www.ietf.org/timezones/data/leap-seconds.list was used but this was outdated on 24-01-2022. #TODO: get most up to date source from somewhere.
    file_leap_seconds_list = os.path.join(os.path.dirname(file_path),'data','leap-seconds.list')
    
    #retrieve leap-seconds.list via url and write to file in hatyan sourcecode folder. If it fails, an old version of the file is used.
    try:
        resp = requests.get(url_leap_seconds_list)
        if resp.status_code==404:
            resp.raise_for_status()
        with open(file_leap_seconds_list, 'wb') as f:
            f.write(resp.content)
    except requests.HTTPError as e: #catch raised 404 error
        print(f'WARNING: leap-seconds.list not retrieved, using local copy. Error message: "{e}"')
    except (FileNotFoundError, PermissionError) as e: #catch problems with writing to file
        print(f'WARNING: leap-seconds.list not written, using local copy. Error message: "{e}"')
        
    #get expiry date from file
    with open(file_leap_seconds_list) as f:
        resp_pd_all = pd.Series(f.readlines())
    expirydate_linestart = '#@'
    expirydate_line = resp_pd_all.loc[resp_pd_all.str.startswith(expirydate_linestart)]
    if len(expirydate_line) != 1:
        raise Exception('less or more than 1 expirydate line found with "{expirydate_linestart}": {expirydate_line}')
    expirydate = refdate + dt.timedelta(seconds=int(expirydate_line.iloc[0].split()[1]))
    
    #get leapsecond list from file
    resp_pd = pd.read_csv(file_leap_seconds_list,comment='#',names=['seconds_since_19000101','leap_seconds'],delim_whitespace=True)
    resp_pd['datetime'] = refdate + pd.to_timedelta(resp_pd['seconds_since_19000101'],unit='S')
    resp_pd['datetime_str'] = resp_pd['datetime'].dt.strftime('%Y-%m-%d')
    leap_seconds_pd = resp_pd.set_index('datetime_str')
    
    return leap_seconds_pd, expirydate


def dT(dateIn,dT_fortran=False):
    """
    Calculates difference between terrestrial time and universal time. Uses a leap-second file that automatically updates via hatyan.astrog.get_leapsecondslist_fromurlorfile(). Dates after the expiry date of the file are also corrected with the last available value. Dates befor 1972 are corrected with the first available value.
    Background is available on https://astro.ukho.gov.uk/nao/miscellanea/DeltaT/
    TT = TAI + 32.184 seconds. TAI - UTC = 37 s (latest value from leap-seconds.list). So UTC = TAI-37s = TT-32.184s-37s
    Uses the international definition unless dT_fortran=True. 
    
    Parameters
    ----------
    dateIn : datetime.datetime or pandas.DatetimeIndex
        Date for correction. Definition makes use of provided year.
    dT_fortran : boolean, optional
        When True, use latest fortran dT and increment value instead of default international definition. The default is False.
    Raises
    ------
    Warning
        Checks if hard-coded values can still be used.

    Returns
    -------
    dT_TT : float
        Difference dT between terrestrial time (TT) and universal time (UT1) in seconds.

    """
    
    if isinstance(dateIn, dt.datetime):
        dateIn = pd.DatetimeIndex([dateIn])
    elif not isinstance(dateIn, pd.DatetimeIndex):
        raise Exception('Input variable date should be datetime or pd.DateTimeIndex')
    
    if dT_fortran: # reproduce fortran dT_TT values with latest dT and increment values from fortran code
        warnings.warn('WARNING: If dT_fortran=True, the last values for dT and its increment are used from the fortran code, this is not accurate for years that are far away from 2012. Use dT_fortran=False since the default approach uses an automatically updated list of leap-seconds instead')
        # historical hard-coded values (taken from FORTRAN comments) from Astronomical Almanac - Reduction of time scales (only last ones are used to reproduce fortran code)
        dT_TTyear     = [ 1980,  1993,  2002,  2012 ] # year of used dT_TTval value
        dT_TTval      = [50.97, 59.35, 64.90, 67.184] # difference between TT and UT1 (32.184s + leap seconds)
        dT_TTinc      = [0.998,  0.70,  0.42,  0.676] # yearly increment of dT curve: (dT_last-dT_5yBefore)/5
        
        # approximation of dT in requested year (dT = TT-UT1)
        dT_TT = (dT_TTval[-1]+dT_TTinc[-1]*(np.array(dateIn.year)-dT_TTyear[-1]))#/86400 

    else: #use most exact approximation
        # get list with leap seconds
        if (dateIn<dt.datetime(1972,1,1)).any():
            warnings.warn('WARNING: The current definition of the relationship between UTC and TAI dates from 1 January 1972. This first dT value is also applied before that date even though this might not be accurate.')
        leap_seconds_pd, expirydate = get_leapsecondslist_fromurlorfile()

        NTP_date = leap_seconds_pd['datetime'].values#tolist()
        leap_sec = leap_seconds_pd['leap_seconds'].values#.tolist()
        #import matplotlib.pyplot as plt
        #fig,ax=plt.subplots()
        #ax.plot(leap_seconds_pd['datetime'], leap_seconds_pd['leap_seconds'])
        ind = np.zeros(shape=dateIn.shape,dtype=int)
        for iD, NTP_date_one in enumerate(NTP_date):
            ind[dateIn>NTP_date_one] = iD
        dT_TT = leap_sec[(ind)] + 32.184 #voor conversie van Terrestrial Time naar UTC (TAI = TT-32.184 = UTC+dT  >> UTC = TT-32.184-dT)

    return dT_TT


def check_crop_dataframe(astrog_df, tFirst, tLast, tzone):
    
    #set timezone, check datetime order and filter datetimerange
    astrog_df['datetime'] = pd.to_datetime(astrog_df['datetime']).dt.tz_localize('UTC',ambiguous=False,nonexistent='shift_forward') # set timezone (UTC)
    astrog_df['datetime'] = astrog_df['datetime'].dt.tz_convert(tzone) #convert timezone to tzone
    if (np.diff(astrog_df.sort_values('datetime').index)!=1).any():
        raise Exception('something went wrong which resulted in off ordering of the dataframe')
    astrog_df_dtnaive = astrog_df['datetime'].dt.tz_localize(None)
    astrog_df = astrog_df[np.logical_and(astrog_df_dtnaive>=tFirst,astrog_df_dtnaive<=tLast)].reset_index(drop=True)
    return astrog_df


def convert_str2datetime(datetime_in_list):
    """
    Tries to convert datetime_in_list (list of str or datetime.datetime) to list of datetime.datetime

    Parameters
    ----------
    datetime_in_list : list of str/dt.datetime/pd.Timestamp
        DESCRIPTION.

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    datetime_out_list : list of pd.Timestamp
        DESCRIPTION.

    """

    datetime_out_list = datetime_in_list
    for iDT, datetime_in in enumerate(datetime_in_list):
        if isinstance(datetime_in,pd._libs.tslibs.timestamps.Timestamp):
            datetime_out = datetime_in
        elif isinstance(datetime_in, dt.datetime):
            datetime_out = pd.Timestamp(datetime_in)
            if hasattr(datetime_out,'tz'):
                if datetime_out.tz != None:
                    raise Exception('tFirst and tLast should be timezone naive dt.datetime or "yyyymmdd" str')
        else:
            try:
                datetime_out = dt.datetime.strptime(datetime_in,'%Y%m%d')
            except:
                raise Exception('date_input should be timezone naive dt.datetime or "yyyymmdd" str')
                
    return datetime_out_list


def convert2perday(dataframeIn, timeformat='%H:%M %Z'):
    """
    converts normal astrog pd.DataFrame to one with the same information restructured per day

    Parameters
    ----------
    dataframeIn : pd.DataFrame
        with columns 'datetime' and 'type_str'.
    timeformat : str, optional
        format of the timestrings in dataframeOut. The default is '%H:%M %Z'.

    Returns
    -------
    dataframeOut : pd.DataFrame 
        The 'datetime' column contains dates, with columns containing all unique 'type_str' values.

    """
    
    dataframeOut = dataframeIn.copy()
    dataframeOut.index = dataframeOut['datetime'].dt.date
    for type_sel in dataframeOut['type_str'].unique():
        dataframeOut[type_sel] = dataframeOut['datetime'][dataframeOut['type_str']==type_sel].dt.strftime(timeformat)
    dataframeOut.drop(['type','type_str'],axis='columns',inplace=True)
    dataframeOut = dataframeOut[~dataframeOut.index.duplicated(keep='first')]
    dataframeOut['datetime'] = dataframeOut.index #overwrite datetime with dates
    
    return dataframeOut


def plot_astrog_diff(pd_python, pd_fortran, typeUnit='-', typeLab=None, typeBand=None, timeBand=None):
    """
    Plots results of FORTRAN and python verison of astrog for visual inspection.
    Top plot shows values or type, middle plot shows time difference, bottom plot shows value/type difference.

    Parameters
    ----------
    pd_python : pandas DataFrame
        DataFrame from astrog (python) with times (UTC).
    pd_fortran : pandas DataFrame
        DataFrame from astrog (FORTRAN) with times (UTC).
    typeUnit : string, optional
        Unit of provided values/types. The default is '-'.
    typeLab : TYPE, optional
        Labels of provided types. The default is ['rise','set'].
    typeBand : list of floats, optional
        Expected bandwith of accuracy of values/types. The default is [-.5,.5].
    timeBand : list of floats, optional
        Expected bandwith of accuracy of times (seconds). The default is [0,60].
    timeLim : list of floats, optional
        Time limits of x-axis. The default is None (takes limits from pd_python).

    Returns
    -------
    fig : figure handle
        Output figure.
    axs : axis handles
        Axes in figure.

    """
    
    if hasattr(pd_python['datetime'].dtype,'tz'):
        pd_python = pd_python.copy() #do not overwrite original dataframe, so make a copy
        pd_python['datetime'] = pd_python['datetime'].dt.tz_localize(None) #Passing None will remove the time zone information preserving local time. https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.Series.dt.tz_localize.html#pandas.Series.dt.tz_localize
    
    typeName = pd_python.columns[1]
    fig, (ax1,ax2,ax3) = plt.subplots(3,1,figsize=(15,9),sharex=True)
    ax1.set_title('%s'%(typeName))
    ax1.plot(pd_python['datetime'], pd_python[typeName], label='python')
    ax1.plot(pd_fortran['datetime'],pd_fortran[typeName],label='FORTRAN',linestyle='dashed')
    if typeName == 'type':
        if typeLab is not None:
            ax1.set_ylim(1,len(typeLab))
            ax1.set_yticks(np.arange(1,len(typeLab)+1,step=1))
            ax1.set_yticklabels(typeLab)
    ax1.set_ylabel('%s [%s]'%(typeName, typeUnit))
    ax1.legend(loc=1)

    ax2.plot(pd_python['datetime'],(pd_python['datetime'] - pd.to_datetime(pd_fortran['datetime'])).dt.total_seconds())
    if timeBand is not None:
        ax2.plot(pd_python['datetime'].iloc[[0,-1]],[timeBand[0],timeBand[0]],color='k',linestyle='dashed')
        ax2.plot(pd_python['datetime'].iloc[[0,-1]],[timeBand[1],timeBand[1]],color='k',linestyle='dashed')
    ax2.set_ylabel('python - FORTRAN [seconds]')
    ax2.set_title('time difference [seconds]')

    ax3.plot(pd_python['datetime'],pd_python[typeName] - pd_fortran[typeName])
    if typeBand is not None:
        ax3.plot(pd_python['datetime'].iloc[[0,-1]],[typeBand[0],typeBand[0]],color='k',linestyle='dashed')
        ax3.plot(pd_python['datetime'].iloc[[0,-1]],[typeBand[1],typeBand[1]],color='k',linestyle='dashed')
    ax3.set_ylabel('python - FORTRAN [%s]'%(typeUnit))
    ax3.set_title('%s difference [%s]'%(typeName, typeUnit))

    ax1.set_xlim(pd_python['datetime'].iloc[[0,-1]])
    fig.tight_layout()

    axs = (ax1,ax2,ax3)
    return fig, axs


