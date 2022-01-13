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


def astrog_culminations(tFirst,tLast,mode_dT='exact',tzone='UTC'):
    """
    Makes use of the definitions dT, astrab and astrac.
    Calculates lunar culminations, parallax and declination. By default the lunar culmination is calculated at coordinates 52,0 (Netherlands,Greenwich).

    Parameters
    ----------
    tFirst : pd.Timestamp, datetime.datetime or string ("yyyymmdd")
        Start of timeframe for output.
    tLast : pd.Timestamp, datetime.datetime or string ("yyyymmdd")
        End of timeframe for output.
    mode_dT : string, optional
        Method to calculate difference between universal time and terrestrial time (dT). Can be 'last' (for fortran reproduction), 'historical' or 'exact' (most accurate). The default is 'exact'.
    tzone : string/dt.timezone, optional
        Timezone to convert the output dataset to. The default is 'UTC'.

    Raises
    ------
    Exception
        Checks input times tFirst and tLast.

    Returns
    -------
    dataCulminations : pandas DataFrame
        datetime:    lunar culmination at Greenwich in UTC (datetime)
        type:        type of culmination (1=lower, 2=upper)
        parallax:    lunar parallax (degrees)
        declination: lunar declination (degrees)

    """
    import pandas as pd
    import numpy as np
    import datetime as dt
    
    from hatyan.hatyan_core import get_hatyan_freqs
    
    # check input times (datetime or string)
    [tFirst,tLast] = convert_str2datetime(datetime_in_list=[tFirst,tLast])
    
    # constants
    EHMINC       = 346.8                                            # increment of ephemeris hour angle of moon (deg/day)
    M2_period_hr = get_hatyan_freqs(['M2']).loc['M2','period [hr]'] # interval between lunar culminations (days)
    
    # first and last datetime in calculation (add enough margin, and an extra day for timezone differences)
    date_first = tFirst-dt.timedelta(hours=M2_period_hr+1*24)
    date_last = tLast+dt.timedelta(hours=M2_period_hr+1*24)

    # estimate culminations (time and type)
    astrabOutput = astrab(date_first,dT(date_first,mode_dT=mode_dT))
    EHMOON = astrabOutput['EHMOON']
    EHMOON[EHMOON>=360] -= 360. #subtract 360 if larger than 360
    # ICUL=1: next culmination is lower culmination
    # ICUL=2: next culmination is upper culmination
    ICUL = (EHMOON[0]/180.).astype(int)+1
    CULEST=pd.date_range(start=date_first+dt.timedelta(days=(180.*ICUL-EHMOON[0])/EHMINC),end=date_last,freq='%iN'%(M2_period_hr*3600*1e9)) #defined freq as M2_period in nanoseconds
    CULTYP = np.empty((len(CULEST),)).astype(int)
    CULTYP[::2]  = ICUL
    CULTYP[1::2] = (ICUL%2)+1
    
    # calculate exact time of culminations
    CULTIM = astrac(CULEST,dT(CULEST,mode_dT=mode_dT),CULTYP)
    astrabOutput = astrab(CULTIM,dT(CULTIM,mode_dT=mode_dT))
    PAR = astrabOutput['PARLAX']/3600.
    DEC = astrabOutput['DECMOO']
    
    # make dataframe and crop for requested timeframe
    dataCulminations = pd.DataFrame({'datetime':CULTIM,'type':CULTYP,'parallax':PAR,'declination':DEC}) #CULTIM.round('S') decreases fortran reproduction
    dataCulminations['type_str'] = dataCulminations['type'].astype(str).replace('1','lowerculmination').replace('2','upperculmination')
    dataCulminations['datetime'] = pd.to_datetime(dataCulminations['datetime']).dt.tz_localize('UTC',ambiguous=False,nonexistent='shift_forward') # set timezone to UTC
    dataCulminations['datetime'] = dataCulminations['datetime'].dt.tz_convert(tzone) #convert timezone to tzone
    
    #filter datetimerange
    dataCulminations_dtnaive = dataCulminations['datetime'].dt.tz_localize(None)
    dataCulminations = dataCulminations[np.logical_and(dataCulminations_dtnaive>=tFirst,dataCulminations_dtnaive<=tLast)].reset_index(drop=True)
    
    return dataCulminations


def astrog_phases(tFirst,tLast,mode_dT='exact',tzone='UTC'):
    """
    Makes use of the definitions dT, astrab and astrac.
    Calculates lunar phases. The lunar phases are independent of coordinates.

    Parameters
    ----------
    tFirst : datetime.datetime or string ("yyyymmdd")
        Start of timeframe for output.
    tLast : datetime.datetime or string ("yyyymmdd")
        End of timeframe for output.
    mode_dT : string, optional
        Method to calculate difference between universal time and terrestrial time (dT). Can be 'last' (for fortran reproduction), 'historical' or 'exact' (most accurate). The default is 'exact'.
    tzone : string/dt.timezone, optional
        Timezone to convert the output dataset to. The default is 'UTC'.

    Raises
    ------
    Exception
        Checks input times tFirst and tLast.

    Returns
    -------
    dataPhases : pandas DataFrame
        datetime:  lunar phase in UTC (datetime)
        type:      type of phase (1=FQ, 2=FM, 3=LQ, 4=NM)

    """
    import pandas as pd
    import numpy as np
    import datetime as dt

    # check input times (datetime or string)
    [tFirst,tLast] = convert_str2datetime(datetime_in_list=[tFirst,tLast])
        
    # constants
    ELOINC = 12.2           # increment of ecliptic elongation of moon-sun (deg/day)
    FASINT = 29.530587981/4 # quarter of a lunar synodic month (days)

    # first and last datetime in calculation (add enough margin (done later), and an extra day for timezone differences)
    date_first = tFirst-dt.timedelta(days=FASINT+1)
    date_last = tLast+dt.timedelta(days=FASINT+1)

    # estimate first lunar phase (time and type), correct first date (FAEST_first to 45 deg from there)
    astrabOutput = astrab(date_first,dT(date_first,mode_dT=mode_dT))
    ELONG = astrabOutput['ELONG']
    ELONG[ELONG>=360] -= 360.
    FAEST_first = date_first - pd.TimedeltaIndex((ELONG-45)%360/ELOINC, unit='D')

    # use the first date to create a new daterange from the correct starting time. The frequency is 29 days, 12 hours and 44 minutes, following from dood_S-dood_H
    date = pd.date_range(start=FAEST_first[0],end=date_last,freq='%iN'%(FASINT*24*3600*1e9))

    # estimate all lunar phases (time and type)
    astrabOutput = astrab(date,dT(date,mode_dT=mode_dT))
    ELONG = astrabOutput['ELONG']
    ELONG[ELONG>=360] -= 360.
    FATYP=(np.array(ELONG/90.).astype(int)+3)%4+1 #make sure the next phase is searched for (modulus to use 'FATYP-1')
    FAEST=date-pd.TimedeltaIndex((90.*FATYP-ELONG)/ELOINC, unit='D')

    # calculate exact time of phase, loop until date_last
    TIMDIF = pd.TimedeltaIndex(-dT(FAEST,mode_dT=mode_dT)+1./2880.,unit='D')
    FATIM = astrac(FAEST,dT(FAEST,mode_dT=mode_dT),FATYP+2)+TIMDIF

    # make dataframe and crop for requested timeframe
    dataPhases = pd.DataFrame({'datetime':FATIM.round('S'),'type':FATYP})
    dataPhases['type_str'] = dataPhases['type'].astype(str).replace('1','FQ').replace('2','FM').replace('3','LQ').replace('4','NM')
    dataPhases['datetime'] = pd.to_datetime(dataPhases['datetime']).dt.tz_localize('UTC',ambiguous=False,nonexistent='shift_forward') # set timezone (UTC)
    dataPhases['datetime'] = dataPhases['datetime'].dt.tz_convert(tzone) #convert timezone to tzone
    if (np.diff(dataPhases.sort_values('datetime').index)!=1).any():
        raise Exception('something went wrong with moonphases which resulted in off ordering of the dataframe, check FAEST_first degree correction')
    
    #filter datetimerange
    dataPhases_dtnaive = dataPhases['datetime'].dt.tz_localize(None)
    dataPhases = dataPhases[np.logical_and(dataPhases_dtnaive>=tFirst,dataPhases_dtnaive<=tLast)].reset_index(drop=True)

    return dataPhases


def astrog_sunriseset(tFirst,tLast,mode_dT='exact',tzone='UTC',lon=5.3876,lat=52.1562):
    """
    Makes use of the definitions dT, astrab and astrac.
    Calculates sunrise and -set at requested location.

    Parameters
    ----------
    tFirst : datetime.datetime or string ("yyyymmdd")
        Start of timeframe for output.
    tLast : datetime.datetime or string ("yyyymmdd")
        End of timeframe for output.
    mode_dT : string, optional
        Method to calculate difference between universal time and terrestrial time (dT). Can be 'last' (for fortran reproduction), 'historical' or 'exact' (most accurate). The default is 'exact'.
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
    dataSun : pandas DataFrame
        datetime: time of rise or set in UTC (datetime)
        type:     type  (1=sunrise, 2=sunset)

    """
    import pandas as pd
    import numpy as np
    import datetime as dt

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
    TIMDIF = pd.TimedeltaIndex(-dT(OPEST,mode_dT=mode_dT)+1./2880.,unit='D')
    OPTIM  = astrac(OPEST,dT(OPEST,mode_dT=mode_dT),np.array( 9),lon=lon,lat=lat)+TIMDIF
    ONTIM  = astrac(ONEST,dT(ONEST,mode_dT=mode_dT),np.array(10),lon=lon,lat=lat)+TIMDIF

    # make dataframe and crop for requested timeframe
    dataSun = pd.DataFrame({'datetime':np.concatenate((OPTIM.round('S'),ONTIM.round('S'))),'type':np.concatenate((np.full(len(OPTIM),1),np.full(len(OPTIM),2)))})
    dataSun = dataSun.sort_values('datetime').reset_index(drop=True)
    dataSun['type_str'] = dataSun['type'].astype(str).replace('1','sunrise').replace('2','sunset')
    dataSun['datetime'] = pd.to_datetime(dataSun['datetime']).dt.tz_localize('UTC',ambiguous=False,nonexistent='shift_forward') # set timezone (UTC)
    dataSun['datetime'] = dataSun['datetime'].dt.tz_convert(tzone) #convert timezone to tzone

    #filter datetimerange
    dataSun_dtnaive = dataSun['datetime'].dt.tz_localize(None)
    dataSun = dataSun[np.logical_and(dataSun_dtnaive>=tFirst,dataSun_dtnaive<=tLast)].reset_index(drop=True)
    
    return dataSun


def astrog_moonriseset(tFirst,tLast,mode_dT='exact',tzone='UTC',lon=5.3876,lat=52.1562):
    """
    Makes use of the definitions dT, astrab and astrac.
    Calculates moonrise and -set at requested location.

    Parameters
    ----------
    tFirst : datetime.datetime or string ("yyyymmdd")
        Start of timeframe for output.
    tLast : datetime.datetime or string ("yyyymmdd")
        End of timeframe for output.
    mode_dT : string, optional
        Method to calculate difference between universal time and terrestrial time (dT). Can be 'last' (for fortran reproduction), 'historical' or 'exact' (most accurate). The default is 'exact'.
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
    dataMoon : pandas DataFrame
        datetime: time of rise or set in UTC (datetime)
        type:     type  (1=moonrise, 2=moonset)

    """
    import pandas as pd
    import numpy as np
    import datetime as dt

    # check input times (datetime or string)
    [tFirst,tLast] = convert_str2datetime(datetime_in_list=[tFirst,tLast])

    # constants
    from hatyan.hatyan_core import get_hatyan_freqs
    M2_period_hr = get_hatyan_freqs(['M2']).loc['M2','period [hr]'] # CULINT
    EHMINC = 346.8 # increment of ephemeris hour angle of moon (deg/day)

    # first and last datetime in calculation (add enough margin, and an extra day for timezone differences)
    date_first = tFirst-dt.timedelta(hours=M2_period_hr+1*24)
    date_last = tLast+dt.timedelta(hours=M2_period_hr+1*24)

    # --- moonrise and -set ---
    # estimate times
    astrabOutput=astrab(date_first,dT(date_first,mode_dT=mode_dT),lon=lon,lat=lat)
    ALTMOO=astrabOutput['ALTMOO']
    EHMOON=astrabOutput['EHMOON']
    # first phenomenon is moonrise
    if ALTMOO < -(0.5667+(0.08+0.2725*astrabOutput['PARLAX'])/3600.):
        OPEST=pd.date_range(start=date_first+dt.timedelta(days=(270.-EHMOON[0])/EHMINC),end=date_last,freq='%iN'%(M2_period_hr*2*3600*1e9))
        ONEST=OPEST+dt.timedelta(hours=M2_period_hr)
    # first phenomenon is moonset
    else:
        ONEST=pd.date_range(start=date_first+dt.timedelta(days=(90.-EHMOON[0])/EHMINC),end=date_last,freq='%iN'%(M2_period_hr*2*3600*1e9))
        OPEST=ONEST+dt.timedelta(hours=M2_period_hr)

    # calculate exact times
    TIMDIF = pd.TimedeltaIndex(-dT(OPEST,mode_dT=mode_dT)+1./2880.,unit='D')
    OPTIM  = astrac(OPEST,dT(OPEST,mode_dT=mode_dT),np.array(7),lon=lon,lat=lat)+TIMDIF
    ONTIM  = astrac(ONEST,dT(ONEST,mode_dT=mode_dT),np.array(8),lon=lon,lat=lat)+TIMDIF

    # make dataframe and crop for requested timeframe
    dataMoon = {'datetime':np.concatenate((OPTIM.round('S'),ONTIM.round('S'))),'type':np.concatenate((np.full(len(OPTIM),1),np.full(len(OPTIM),2)))}
    dataMoon = pd.DataFrame(dataMoon).sort_values('datetime').reset_index(drop=True)
    dataMoon['type_str'] = dataMoon['type'].astype(str).replace('1','moonrise').replace('2','moonset')
    dataMoon['datetime'] = pd.to_datetime(dataMoon['datetime']).dt.tz_localize('UTC',ambiguous=False,nonexistent='shift_forward') # set timezone (UTC)
    dataMoon['datetime'] = dataMoon['datetime'].dt.tz_convert(tzone) #convert timezone to tzone

    #filter datetimerange
    dataMoon_dtnaive = dataMoon['datetime'].dt.tz_localize(None)
    dataMoon = dataMoon[np.logical_and(dataMoon_dtnaive>=tFirst,dataMoon_dtnaive<=tLast)].reset_index(drop=True)

    return dataMoon


def astrog_anomalies(tFirst,tLast,mode_dT='exact',tzone='UTC'):
    """
    Makes use of the definitions dT, astrab and astrac.
    Calculates lunar anomalies. The lunar anomalies are independent of coordinates.

    Parameters
    ----------
    tFirst : datetime.datetime or string ("yyyymmdd")
        Start of timeframe for output.
    tLast : datetime.datetime or string ("yyyymmdd")
        End of timeframe for output.
    mode_dT : string, optional
        Method to calculate difference between universal time and terrestrial time (dT). Can be 'last' (for fortran reproduction), 'historical' or 'exact' (most accurate). The default is 'exact'.
    tzone : string/dt.timezone, optional
        Timezone to convert the output dataset to. The default is 'UTC'.

    Raises
    ------
    Exception
        Checks input times tFirst and tLast.

    Returns
    -------
    dataAnomaly : pandas DataFrame
        datetime:   lunar anomaly in UTC (datetime)
        type:       type of anomaly (1=perigeum, 2=apogeum)

    """
    import pandas as pd
    import numpy as np
    import datetime as dt

    # check input times (datetime or string)
    [tFirst,tLast] = convert_str2datetime(datetime_in_list=[tFirst,tLast])

    # constants
    ANMINC   = 13.06       # increment of anomaly of moon (deg/day)
    ANOINT   = 27.554551/2 # half of a lunar anomalistic month (days)

    # first and last datetime in calculation (add enough margin, and an extra day for timezone differences)
    date_first = tFirst-dt.timedelta(days=ANOINT+1)
    date_last = tLast+dt.timedelta(days=ANOINT+1)

    # estimate first lunar anomaly (time and type)
    astrabOutput=astrab(date_first,dT(date_first,mode_dT=mode_dT))
    DPAXDT=astrabOutput['DPAXDT']
    ANM   =astrabOutput['ANM']
    if DPAXDT>0.:
        # ANOTYP=1: perigeum first
        if ANM<90.:
            ANM=ANM+360.
        ANOEST = pd.date_range(start=date_first+dt.timedelta(days=(360.-ANM[0])/ANMINC),end=date_last,freq='%iN'%(ANOINT*24*3600*1e9))
        ANOTYP = np.empty((len(ANOEST),)).astype(int)
        ANOTYP[::2]  = 1
        ANOTYP[1::2] = 2

    elif DPAXDT<=0.:
        # ANOTYP=2: apogeum first
        if ANM>270.:
            ANM=ANM-360.
        ANOEST = pd.date_range(start=date_first+dt.timedelta(days=(180.-ANM[0])/ANMINC),end=date_last,freq='%iN'%(ANOINT*24*3600*1e9))
        ANOTYP = np.empty((len(ANOEST),)).astype(int)
        ANOTYP[::2]  = 2
        ANOTYP[1::2] = 1

    # calculate exact times
    TIMDIF = pd.TimedeltaIndex(-dT(ANOEST,mode_dT=mode_dT)+1./48., unit='D')
    ANOTIM = astrac(ANOEST,dT(ANOEST,mode_dT=mode_dT),ANOTYP+14)+TIMDIF

    # make dataframe and crop for requested timeframe
    dataAnomaly = pd.DataFrame({'datetime':ANOTIM.round('S'),'type':ANOTYP})
    dataAnomaly['type_str'] = dataAnomaly['type'].astype(str).replace('1','perigeum').replace('2','apogeum')
    dataAnomaly['datetime'] = pd.to_datetime(dataAnomaly['datetime']).dt.tz_localize('UTC',ambiguous=False,nonexistent='shift_forward') # set timezone (UTC)
    dataAnomaly['datetime'] = dataAnomaly['datetime'].dt.tz_convert(tzone) #convert timezone to tzone

    #filter datetimerange
    dataAnomaly_dtnaive = dataAnomaly['datetime'].dt.tz_localize(None)
    dataAnomaly = dataAnomaly[np.logical_and(dataAnomaly_dtnaive>=tFirst,dataAnomaly_dtnaive<=tLast)].reset_index(drop=True)
    
    return dataAnomaly


def astrog_seasons(tFirst,tLast,mode_dT='exact',tzone='UTC'):
    """
    Makes use of the definitions dT, astrab and astrac.
    Calculates astronomical seasons. The seasons are independent of coordinates.

    Parameters
    ----------
    tFirst : datetime.datetime or string ("yyyymmdd")
        Start of timeframe for output.
    tLast : datetime.datetime or string ("yyyymmdd")
        End of timeframe for output.
    mode_dT : string, optional
        Method to calculate difference between universal time and terrestrial time (dT). Can be 'last' (for fortran reproduction), 'historical' or 'exact' (most accurate). The default is 'exact'.
    tzone : string/dt.timezone, optional
        Timezone to convert the output dataset to. The default is 'UTC'.

    Raises
    ------
    Exception
        Checks input times tFirst and tLast.

    Returns
    -------
    dataSeasons : pandas DataFrame
        datetime:   start of astronomical season in UTC (datetime)
        type:       type of astronomical season (1=spring, 2=summer, 3=autumn, 4=winter)

    """
    import pandas as pd
    import numpy as np
    import datetime as dt

    # check input times (datetime or string)
    [tFirst,tLast] = convert_str2datetime(datetime_in_list=[tFirst,tLast])

    # estimate start of seasons (time and type)
    SEIEST = pd.date_range(start=dt.datetime(tFirst.year,int(np.ceil(tFirst.month/3)*3),1),end=tLast+dt.timedelta(days=1),freq='%iMS'%(3))+dt.timedelta(days=20)
    SEITYP = (SEIEST.month/3).astype(int)

    # calculate exact times, loop until tLast
    TIMDIF = pd.TimedeltaIndex(-dT(SEIEST,mode_dT=mode_dT)+1./2880., unit='D') # conversion to UTC
    SEITIM = astrac(SEIEST,dT(SEIEST,mode_dT=mode_dT),SEITYP+10)+TIMDIF

    # make dataframe and crop for requested timeframe
    dataSeasons = pd.DataFrame({'datetime':SEITIM.round('S'),'type':SEITYP})
    dataSeasons['type_str'] = dataSeasons['type'].astype(str).replace('1','spring').replace('2','summer').replace('3','autumn').replace('4','winter')
    dataSeasons['datetime'] = pd.to_datetime(dataSeasons['datetime']).dt.tz_localize('UTC',ambiguous=False,nonexistent='shift_forward') # set timezone (UTC)
    dataSeasons['datetime'] = dataSeasons['datetime'].dt.tz_convert(tzone) #convert timezone to tzone

    #filter datetimerange
    dataSeasons_dtnaive = dataSeasons['datetime'].dt.tz_localize(None)
    dataSeasons = dataSeasons[np.logical_and(dataSeasons_dtnaive>=tFirst,dataSeasons_dtnaive<=tLast)].reset_index(drop=True)
    
    return dataSeasons


def astrab(date,dT_TT,lon=5.3876,lat=52.1562):
    """
    Python version of astrab.f in FORTRAN 77
    Calculates 18 astronomical parameters at requested time.

    Parameters
    ----------
    date : datetime.datetime or pandas.DatetimeIndex
        Requested time for calculation.
    dT_TT : float
        Difference between terrestrial and universal time in days.
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
    import pandas as pd
    import numpy as np
    import datetime as dt
    
    # check input
    if isinstance(date, pd.DatetimeIndex):
        pass
    elif isinstance(date, dt.datetime):
        date = pd.DatetimeIndex([date])
    else:
        raise Exception('Input variable date should be datetime or pd.DateTimeIndex')
    
    if np.abs(lon)>180:
        raise Exception('Input variable longitude larger than 180deg')
    if np.abs(lat)>90:
        raise Exception('Input variable latitude larger than 90deg')

    # constants - general
    EPOCH  = dt.datetime(1899, 12, 31, 12, 0, 0) # 1900.0 # -12h shift because julian date 0 is at noon?
    # The average orbital elements of the celestial bodies are calculated for the epoch 1900.0.
    # The values are corrected for the year 1990. Intitial values are from the vernal equinox.

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
    GAMMOO =-9.2421851E-4  # increment longitude lunar orbital node lunar (rad/day)
    INMOON = 0.089804108   # inclination lunar orbit (rad)
    PARMEA = 3422.608      # mean horizontal lunar parallax (arcseconds)

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
    OBZERO = 0.40931977    # inclination of ecliptic (rad)
    OBINC  =-6.21937E-9    # increment inclination of ecliptic (rad/day)

    # constants - vernal equinox
    ARZERO = 4.881523      # ephemeris hour angle of vernal equinox (rad)
    NARIES = 6.30038809878 # increment ephemeris hour angle of vernal equinox (rad/day)

    # constants - lunar orbital disturbances
    # selected from Brown's Tables of the Motion of the Moon (1909)
    # in longitude
    #TODO: tabellen naar file
    distP = {'col1':[        0,            0,            1,            1,        1,        1,            1,            0,            0,            0,            0,            0,            2,            2,            2,            2,            1,            1,            1,            1,            1,            1,            1,            0,            0,            0,            0,            0,            1,            1,            1,            0,            3,            3,            3,            3,            2,            2,            2,            2,            2,            2,            1,            1,            1,            1,            1,            1,            1,            1,            1,            0,            0,            2,            2,            1,            1,            4,            4,            2,            2],
             'col2':[        0,            0,            0,            0,        0,        0,            0,            1,            1,            1,            1,            0,            0,            0,            0,            0,            1,            1,            1,            1,           -1,           -1,           -1,            2,            2,            0,            0,            0,            0,            0,            0,            1,            0,            0,            0,            0,            1,            1,            1,           -1,           -1,           -1,            2,            2,           -2,           -2,            0,            0,            0,            0,            0,            1,            1,            0,            0,            1,           -1,            0,            0,            0,            0],
             'col3':[        0,            0,            0,            0,        0,        0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            2,            2,            2,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            2,            2,           -2,           -2,           -2,            2,           -2,            0,            0,            0,            0,            0,            0,            2,           -2],
             'col4':[        4,            2,            4,            2,        0,       -2,           -4,            2,            0,           -2,           -4,            1,            2,            0,           -2,           -4,            2,            0,           -2,           -4,            2,            0,           -2,            0,           -2,            2,            0,           -2,            1,           -1,           -3,            1,            2,            0,           -2,           -4,            0,           -2,           -4,            2,            0,           -2,            0,           -2,            0,           -2,            2,            0,            2,            0,           -2,           -2,            2,           -1,           -3,            1,           -1,            0,           -2,            0,            0],
             'col5':[   13.902,     2369.902,        1.979,      191.953,22639.500,-4586.426,      -38.428,      -24.420,     -666.608,     -164.773,       -1.877,     -125.154,       14.387,      769.016,     -211.656,      -30.773,       -2.921,     -109.420,     -205.499,       -4.391,       14.577,      147.361,       28.475,       -7.486,       -8.096,       -5.741,     -411.608,      -55.173,       -8.466,       18.609,        3.215,       18.023,        1.060,       36.124,      -13.193,       -1.187,       -7.649,       -8.627,       -2.740,        1.181,        9.703,       -2.494,       -1.167,       -7.412,        2.580,        2.533,        -.992,      -45.099,       -6.382,       39.532,        9.366,       -2.152,       -1.440,        1.750,        1.225,        1.267,       -1.089,        1.938,        -.952,       -3.996,       -1.298],}
    distP = pd.DataFrame(distP)
    # in latitude
    distC = {'col1':[     0,     2,     3,     0,     0,     1,     1,    -1,    -1],
             'col2':[     0,     0,     0,     1,     2,     1,     1,     1,     1],
             'col3':[     0,     0,     0,     0,     0,     0,     0,     0,     0],
             'col4':[     1,    -2,     0,     0,    -2,     2,    -2,     0,    -2],
             'col5':[ -.725, 5.679,-1.300,-1.302, -.740,  .787, 2.056,  .679,-1.540]}
    distC = pd.DataFrame(distC)
    distS = {'col1':[       0,       0,       0,       1,       1,       1,       1,       1,       1,       2,       2,       2,       2,       3,       3,       0,       0,       0,       0,       0,       0,       1,       1,       1,       1,      -1,      -1,      -1,       2,       2,       2,       1,       0,       1,      -1],
             'col2':[       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       1,       1,       1,       1,       1,       2,       1,       1,       1,       1,       1,       1,       1,       1,       1,      -1,       2,       0,       0,       0],
             'col3':[       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       2,       2,       2],
             'col4':[       1,       2,       4,       4,       2,       0,      -2,      -3,      -4,       2,       0,      -2,      -4,       0,      -2,       2,       1,       0,      -2,      -4,      -2,       2,       0,      -2,      -4,       2,       0,      -2,       0,      -2,       0,      -2,      -2,      -2,       0],
             'col5':[ -112.79, 2373.36,   14.06,    6.98,  192.72,22609.07,-4578.13,    5.44,  -38.64,   14.78,  767.96, -152.53,  -34.07,   50.64,  -16.40,  -25.10,   17.93, -126.98, -165.06,   -6.46,  -16.35,  -11.75, -115.18, -182.36,   -9.66,  -23.59, -138.76,  -31.70,  -10.56,   -7.59,   11.67,   -6.12,  -52.14,   -9.52,  -85.13]}
    distS = pd.DataFrame(distS)
    distN = {'col1':[       0,       0,       1,       1,      -1,      -1,      -2,      -2,       0,       0],
             'col2':[       0,       0,       0,       0,       0,       0,       0,       0,       1,      -1],
             'col3':[       1,       1,       1,       1,       1,       1,       1,       1,       1,       1],
             'col4':[      -2,      -4,      -2,      -4,       0,      -2,       0,      -2,      -2,      -2],
             'col5':[-526.069,  -3.352,  44.297,  -6.000,  20.599, -30.598, -24.649,  -2.000, -22.571,  10.985]}
    distN = pd.DataFrame(distN)
    # in parallax
    distR = {'col1':[        0,        0,        1,        1,        1,        1,        0,        0,        0,        0,        2,        2,        2,        2,        1,        1,        1,        1,        1,        0,        1,        0,        3,        3,        2,        2,        1],
             'col2':[        0,        0,        0,        0,        0,        0,        1,        1,        1,        0,        0,        0,        0,        0,        1,        1,       -1,       -1,       -1,        0,        0,        1,        0,        0,        1,       -1,        0],
             'col3':[        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        2,        0,        0,        0,        0,        0,        0,       -2],
             'col4':[        4,        2,        2,        0,       -2,       -4,        2,        0,       -2,        1,        2,        0,       -2,       -4,        0,       -2,        2,        0,       -2,       -2,        1,        1,        0,       -2,        0,        0,        0],
             'col5':[    .2607,  28.2333,   3.0861, 186.5398,  34.3117,    .6008,   -.2993,   -.3988,   1.9135,   -.9781,    .2833,  10.1657,   -.3039,    .3722,   -.9469,   1.4404,    .2297,   1.1502,   -.2252,   -.1052,   -.1093,    .1494,    .6215,   -.1187,   -.1038,    .1268,   -.7136]}
    distR = pd.DataFrame(distR)

    # process input data
    TIME  = np.array((date-EPOCH).total_seconds()/86400) # days since 1900
    if (TIME<0).any() or (TIME>70000).any():
        raise Exception('Requested time out of range (1900-2091)')
    RLONG=-np.deg2rad(lon)
    RLATI= np.deg2rad(lat)

    # mean orbital elements
    # of sun, including additive correction
    LABSUN=(LABOS -np.deg2rad(6.22/3600) +NSUN*TIME) % (2*np.pi)
    PERSUN=PEROS+BETSUN*TIME
    # of moon
    LABMOO=(LABOM+NMOON *TIME) % (2*np.pi)
    PERMOO=(PEROM+BETMOO*TIME) % (2*np.pi)
    NODMOO=NODOM+GAMMOO*TIME
    # additive corrections for moon
    LABMOO=LABMOO+np.deg2rad(( 0.79 +14.27 *np.sin(6.0486+6.34913E-5*TIME)+7.261*np.sin(NODMOO))/3600)
    PERMOO=PERMOO+np.deg2rad((-1.966- 2.076*np.sin(NODMOO))/3600)
    NODMOO=NODMOO+np.deg2rad(( 0.59 +95.96 *np.sin(NODMOO)+15.58*np.sin(NODMOO+4.764400))/3600)
    CINMOO=np.deg2rad((-8.636*np.cos(NODMOO)-1.396*np.cos(NODMOO+4.764400))/3600)

    # mean heliocentric elongation Venus-Earth, Earth-Jupiter, Earth-Mars and Earth-Saturnus
    VENTER=VENTZE+TIME*VENTIN
    TERJUP=TJUPZE+TIME*TJUPIN
    TERMAR=TMARZE+TIME*TMARIN
    TERSAT=TSATZE+TIME*TSATIN

    # Nutation in longitude
    CNULON=-17.248*np.sin(NODMOO)-1.273*np.sin(2.*LABSUN)
    # Nutation in inclination of ecliptic
    CNUTOB=  9.21 *np.cos(NODMOO)+0.55 *np.cos(2.*LABSUN)

    # inclination of ecliptic, sine and cosine of inclination
    OBLIQ=OBZERO+TIME*OBINC+np.deg2rad(CNUTOB/3600)
    SINOB=np.sin(OBLIQ)
    COSOB=np.cos(OBLIQ)

    # parameters of sun and moon for disturbance equations
    ANM=LABMOO-PERMOO # output value 18: mean lunar anomaly (rad)
    ANS=LABSUN-PERSUN
    FNO=LABMOO-NODMOO
    ELO=LABMOO-LABSUN
    DAM=NMOON-BETMOO
    DAS=NSUN-BETSUN
    DFN=NMOON-GAMMOO
    DEL=NMOON-NSUN

    # disturbance equations
    # disturbances in longitude sun
    # equation of center
    CENTR=(6910.10 -17.33*TIME/36525)*np.sin(ANS)+72.01*np.sin(2.*ANS) + 1.05*np.sin(3.*ANS)
    # planetary disturbances
    PLANET= (4.838*np.cos(   VENTER       +1.5708)+
             5.526*np.cos(2.*VENTER       +1.5723)+
             0.666*np.cos(3.*VENTER       +4.7195)+
             2.497*np.cos(2.*VENTER   -ANS+4.4986)+
             1.559*np.cos(3.*VENTER   -ANS+1.3607)+
             1.024*np.cos(3.*VENTER-2.*ANS+0.8875)+
             7.208*np.cos(   TERJUP       +1.5898)+
             2.731*np.cos(2.*TERJUP       +4.7168)+
             2.600*np.cos(   TERJUP   -ANS+3.0503)+
             1.610*np.cos(2.*TERJUP   -ANS+5.1068)+
             0.556*np.cos(3.*TERJUP   -ANS+3.0946)+
             2.043*np.cos(2.*TERMAR       +1.5660)+
             1.770*np.cos(2.*TERMAR   -ANS+5.3454)+
             0.585*np.cos(4.*TERMAR-2.*ANS+3.2432)+
             0.500*np.cos(4.*TERMAR   -ANS+5.5317)+
             0.425*np.cos(3.*TERMAR   -ANS+5.5449)+
             0.419*np.cos(   TERSAT       +1.5767)+
             0.320*np.cos(   TERSAT   -ANS+4.5242))
    # geometric disturbance by the moon
    GEOM=6.454*np.sin(ELO)-0.424*np.sin(ELO-ANM)
    # aberration (correction for optical path)
    ABER=-(20.496+.344*np.cos(ANS))
    # longitude sun with all disturbances
    LONSUN=LABSUN+np.deg2rad((CENTR+PLANET+GEOM+CNULON+ABER)/3600) # output value 8: solar longitude (rad)
    # relative distance from sun
    DISSUN=1.000140 - 0.016712*np.cos(ANS) - 0.000140*np.cos(2.*ANS) # output value 11: relative distance earth-sun (astronomical units)
    # longitude moon with all disturbances
    CLONM=0
    # TODO: for loopjes vervangen door array multiplicatie
    for i in range(0,len(distP)):
        CLONM=CLONM+np.sin(distP['col1'][i]*ANM+distP['col2'][i]*ANS+distP['col3'][i]*FNO+distP['col4'][i]*ELO) * distP['col5'][i]
    LONMOO=LABMOO+np.deg2rad((CNULON+CLONM)/3600) # output value 15: lunar longitude (rad)
    # latitude moon with all disturbances
    CLM=0
    for i in range(0,len(distC)):
        CLM=CLM+np.cos(distC['col1'][i]*ANM+distC['col2'][i]*ANS+distC['col3'][i]*FNO+distC['col4'][i]*ELO) * distC['col5'][i]
    SLM=0
    for i in range(0,len(distS)):
        SLM=SLM+np.sin(distS['col1'][i]*ANM+distS['col2'][i]*ANS+distS['col3'][i]*FNO+distS['col4'][i]*ELO) * distS['col5'][i]
    SF=FNO+np.deg2rad(SLM/3600)
    NLM=0
    for i in range(0,len(distN)):
        NLM=NLM+np.sin(distN['col1'][i]*ANM+distN['col2'][i]*ANS+distN['col3'][i]*FNO+distN['col4'][i]*ELO) * distN['col5'][i]
    LATMOO=((18519.7+CLM)*np.sin(SF) - 6.241*np.sin(3.*SF) + NLM)*np.deg2rad((1.+CINMOO/INMOON)/3600) # output value 16: lunar latitude (rad)

    # lunar parallax
    PARLAX=PARMEA # output value 3: lunar horizontal parallax (arcseconds)
    for i in range(0,len(distR)):
        PARLAX=PARLAX+np.cos(distR['col1'][i]*ANM+distR['col2'][i]*ANS+distR['col3'][i]*FNO+distR['col4'][i]*ELO) * distR['col5'][i]
    # derrivative of lunar parallax
    DPAXDT=0. # output value 4: time derivative of parallax (arcseconds/day)
    for i in range(0,len(distR)):
        DPAXDT=(DPAXDT-np.sin(distR['col1'][i]*ANM+distR['col2'][i]*ANS+distR['col3'][i]*FNO+distR['col4'][i]*ELO) *
                (distR['col1'][i]*DAM+distR['col2'][i]*DAS+distR['col3'][i]*DFN+distR['col4'][i]*DEL) * distR['col5'][i])

    # ecliptic elongation moon-sun
    ELONG=LONMOO-LONSUN # output value 6: ecliptic elongation moon-sun (rad)

    # transformation to equatorial coordinates
    TEMP1=np.sin(LONSUN)
    TEMP2=np.cos(LONSUN)
    RASUN=np.arctan2(TEMP1*COSOB,TEMP2)                          # output value 13: solar right ascension (rad)
    TEMP3=TEMP1*SINOB
    DECSUN=np.arcsin(TEMP3)                                      # output value 10: solar declination (rad)
    TEMP4=np.sin(LONMOO)
    TEMP5=np.cos(LONMOO)
    TEMP6=np.sin(LATMOO)
    TEMP7=np.cos(LATMOO)
    RAMOON=np.arctan2(TEMP7*TEMP4*COSOB-TEMP6*SINOB,TEMP7*TEMP5) # output value 17: lunar right ascension (rad)
    TEMP8=TEMP7*TEMP4*SINOB+TEMP6*COSOB
    DECMOO=np.arcsin(TEMP8)                                      # output value 2: lunar declination (rad)
    EQELON=RAMOON-RASUN                                          # output value 9: equatorial elongaton moon-sun (rad)

    # uurhoeken
    EHARI  =(ARZERO+NARIES*TIME+COSOB*np.deg2rad(CNULON/3600)) % (2*np.pi)    # output value 12: ephemeris hour angle of vernal equinox (rad)
    LHARI  =EHARI-dT_TT*NARIES-RLONG                             # local hour angle of vernal equinox (rad)
    EHSUN  =EHARI-RASUN                                          # output value 14: solar ephemeris hour angle (rad)
    EHMOON =EHARI-RAMOON                                         # output value  1: lunar ephemeris hour angle (rad)
    LHSUN  =LHARI-RASUN                                          # local solar hour angle (rad)
    LHMOON =LHARI-RAMOON                                         # local lunar hour angle (rad)

    # transformation to local coordinates
    ARGUM=TEMP3*np.sin(RLATI)+np.cos(DECSUN)*np.cos(RLATI)*np.cos(LHSUN)
    ALTSUN=np.nan_to_num(np.arcsin(ARGUM),nan=np.copysign(np.pi/2,ARGUM))  # output value 7: solar altitude (rad). Makes use of np.nan_to_num, filling in condition from if-statement
    ARGUM=TEMP8*np.sin(RLATI)+np.cos(DECMOO)*np.cos(RLATI)*np.cos(LHMOON)
    ALTMOO=np.arcsin(ARGUM)
    ALTMOO=np.nan_to_num(ALTMOO-np.cos(ALTMOO)*np.deg2rad(PARLAX/3600),nan=np.copysign(np.pi/2,ARGUM)) # output value 5: lunar altitude (rad). Makes use of np.nan_to_num, filling in condition from if-statement

    # summarize in dataframe and convert output to degrees
    #TODO: kan met minder output en zonder dict?
    astrabOutput = {'EHMOON':((np.rad2deg(EHMOON)+ 720.-90.) % 360.)+90.,
                    'DECMOO':  np.rad2deg(DECMOO),
                    'PARLAX':  PARLAX,
                    'DPAXDT':  DPAXDT,
                    'ALTMOO':  np.rad2deg(ALTMOO),
                    'ELONG' :((np.rad2deg(ELONG )+1080.-45.) % 360.)+45.,
                    'ALTSUN':  np.rad2deg(ALTSUN),
                    'LONSUN':((np.rad2deg(LONSUN)+ 720.-45.) % 360.)+45.,
                    'EQELON': (np.rad2deg(EQELON)+ 720.)     % 360.,
                    'DECSUN':  np.rad2deg(DECSUN),
                    'DISSUN':  DISSUN,
                    'EHARI' :  EHARI,
                    'RASUN' : (np.rad2deg(RASUN )+ 360.)     % 360.,
                    'EHSUN' : (np.rad2deg(EHSUN )+ 720.)     % 360.,
                    'LONMOO': (np.rad2deg(LONMOO)+ 720.)     % 360.,
                    'LATMOO':  np.rad2deg(LATMOO),
                    'RAMOON': (np.rad2deg(RAMOON)+ 360.)     % 360.,
                    'ANM'   : (np.rad2deg(ANM)   + 720.)     % 360.}

    return astrabOutput


def astrac(timeEst,dT_TT,mode,lon=5.3876,lat=52.1562):
    """
    Python version of astrac.f in FORTRAN 77.
    Calculates exact time of requested astronomical phenomenon.

    Parameters
    ----------
    timeEst : datetime.datetime or pandas.DatetimeIndex
        Estimated time for iteration.
    dT_TT : float
        Difference between terrestrial and universal time in days.
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
    import pandas as pd
    import numpy as np
    import datetime as dt

    if isinstance(timeEst, pd.DatetimeIndex):
        pass
    elif isinstance(timeEst, dt.datetime):
        timeEst = pd.DatetimeIndex([timeEst])
    else:
        raise Exception('Input variable date should be datetime or pd.DateTimeIndex')

    # constants - iteration targets
    ANGLE  = np.array([180.   , 360.   , 90.    , 180.    ,   270.  , 360.    ,   -.5667,    -.5667,     -.8333,     -.8333, 360.     , 90.     , 180.     , 270.     ,   0.   ,  0.   ])
    CRITER = np.array([  0.001,   0.001,  0.0001,   0.0001,   0.0001,   0.0001,   0.0002,    0.0002,     0.0002,     0.0002,   0.00001,  0.00001,   0.00001,   0.00001,   0.002,  0.002])
    RAT    = np.array([346.8  , 346.8  , 12.2   ,  12.2   ,    12.2 ,  12.2   , 346.8   , -346.8   ,   360.    ,  -360.    ,   1.     ,  1.     ,   1.     ,   1.     , -10.08 , 10.08 ])
    ANG = ANGLE[mode-1] # required value after iteration
    CRIT = CRITER[mode-1] # allowed difference between ANG and iteration result
    RATE = RAT[mode-1] # estimated change per day for iteration
    if (mode>=7).any() and (mode<=10).any(): # correct RATE in case of rise and set for latitude
        if np.abs(lat)>59:
            raise Exception('Latitude to close to poles (>59deg), cannot take polar days and nights into account')
        RATE=RATE*np.cos(np.deg2rad(lat))

    # define astrab output parameter corresponding to requested mode
    #TODO: mode omschrijven naar leesbare naam en code verwerken?
    if ((mode== 1) | (mode== 2)).all():
        IPAR = 'EHMOON'
    elif ((mode>= 3) & (mode<= 6)).all():
        IPAR = 'ELONG'
    elif ((mode>= 7) & (mode<= 8)).all():
        IPAR = 'ALTMOO'
    elif ((mode>= 9) & (mode<=10)).all():
        IPAR = 'ALTSUN'
    elif ((mode>=11) & (mode<=14)).all():
        IPAR = 'LONSUN'
    elif ((mode>=15) & (mode<=16)).all():
        IPAR = 'DPAXDT'
    else:
        raise Exception('Requested mode (%s) not recognized' % mode)

    # calculate value at start of iteration
    TNEW=timeEst
    astrabOutput = astrab(TNEW,dT_TT,lon=lon,lat=lat)
    PNEW=astrabOutput[IPAR]

    # iterate until criterium is reached or max 20 times
    ITER=1
    while (abs(ANG-PNEW) > CRIT).any():# and ITER <=20:
        if ITER==11:
            print(ITER)
            print('')
        TOLD=TNEW
        POLD=PNEW
        if (mode==7).any() or (mode==8).any(): # correction for semidiameter moon
            ANG=ANGLE[mode-1]-(0.08+0.2725*astrabOutput['PARLAX'])/3600.
        TNEW=TOLD+pd.TimedeltaIndex(np.nan_to_num((ANG-POLD)/RATE),unit='D') #nan_to_num to make sure no NaT output in next iteration
        astrabOutput = astrab(TNEW,dT_TT,lon=lon,lat=lat)
        ITER=ITER+1
        PNEW=astrabOutput[IPAR]
        RATE=np.array((PNEW-POLD)/((TNEW-TOLD).total_seconds()/86400))
        if ITER>20:
            raise Exception('Stopped after %s iterations, datetime=%s' %(ITER-1,TNEW))
    TIMOUT=TNEW#.round('S') # rounding everything to seconds reduces the accuracy of the reporduction of FORTRAN culmination times

    return TIMOUT


def dT(dateIn,mode_dT='exact'):
    """
    Calculates difference between terrestrial time and universal time.
    Current hard-coded values valid until 2023, update arrays afterwards.

    Parameters
    ----------
    dateIn : datetime.datetime or pandas.DatetimeIndex
        Date for correction. Definition makes use of provided year.
    mode : string, optional
        'last': using the last hard-coded value (as last FORTRAN version)
        'historical': using all (previous) hard-coded values (historical FORTRAN versions)
        'exact' (default): determine dT based on number of leap seconds (follows international definition)

    Raises
    ------
    Warning
        Checks if hard-coded values can still be used.

    Returns
    -------
    dT_TT : float
        Difference dT between terrestrial time (TT) and universal time (UT1) in seconds

    """
    import pandas as pd
    import numpy as np
    import datetime as dt

    if isinstance(dateIn, pd.DatetimeIndex):
        pass
    elif isinstance(dateIn, dt.datetime):
        dateIn = pd.DatetimeIndex([dateIn])
    else:
        raise Exception('Input variable date should be datetime or pd.DateTimeIndex')

    if mode_dT=='last' or mode_dT=='historical': # use approximation of dT based on hard-coded values
        # historical hard-coded values (taken from FORTRAN comments) from Astronomical Almanac - Reduction of time scales
        dT_TTyear     = [ 1980,  1993,  2002,  2012 ] # year of used dT_TTval value
        dT_TTval      = [50.97, 59.35, 64.90, 67.184] # difference between TT and UT1 (32.184s + leap seconds)
        dT_TTinc      = [0.998,  0.70,  0.42,  0.676] # yearly increment of dT curve: (dT_last-dT_5yBefore)/5
        if (dateIn.year>dT_TTyear[-1]+11).any(): # check if the last hard-coded value can still be used
            print('WARNING: update hard-coded arrays in definition dT to continue using astrog for modes "last" and "historical"')
        #SCL: changed way to set ind to work with pandas DatetimeIndex
        ind = np.full(len(dateIn),-1) # use the last hard-coded value (same result as last available FORTRAN version)
        if mode_dT=='historical': # use the historical hard-coded values to reproduce results from older FORTRAN versions
            ind[dateIn<dt.datetime(2013,10, 1)]=2
            ind[dateIn<dt.datetime(2001,11,22)]=1
            ind[dateIn<dt.datetime(1994, 8,23)]=0
        dT_TTyear = list(map(dT_TTyear.__getitem__, ind))
        dT_TTval  = list(map(dT_TTval.__getitem__,  ind))
        dT_TTinc  = list(map(dT_TTinc.__getitem__,  ind))

        # approximation of dT in requested year (dT = TT-UT1)
        dT_TT = (dT_TTval+dT_TTinc*(np.array(dateIn.year)-dT_TTyear))/86400 # approximation of dT_TT in [year]

    elif mode_dT=='exact':
        # hard-coded list with leap seconds
        NTP = {'1972-01-01':10,
               '1972-07-01':11,
               '1973-01-01':12,
               '1974-01-01':13,
               '1975-01-01':14,
               '1976-01-01':15,
               '1977-01-01':16,
               '1978-01-01':17,
               '1979-01-01':18,
               '1980-01-01':19,
               '1981-07-01':20,
               '1982-07-01':21,
               '1983-07-01':22,
               '1985-07-01':23,
               '1988-01-01':24,
               '1990-01-01':25,
               '1991-01-01':26,
               '1992-07-01':27,
               '1993-07-01':28,
               '1994-07-01':29,
               '1996-01-01':30,
               '1997-07-01':31,
               '1999-01-01':32,
               '2006-01-01':33,
               '2009-01-01':34,
               '2012-07-01':35,
               '2015-07-01':36,
               '2017-01-01':37,
               }
        NTP_valid = dt.datetime(2021, 12, 28) #display warning after this date

        NTP_date = [dt.datetime.strptime(x,'%Y-%m-%d') for x in [*NTP]]
        leap_sec = list(NTP.values())
        ind = np.full(len(dateIn),-1)
        for iD in range(len(NTP_date)-2,-1,-1):
            ind[dateIn<NTP_date[iD]]=iD
        dT_TT  = (np.array(32.184) + list(map(leap_sec.__getitem__,  ind)) )/86400

        if (dateIn>NTP_valid).any():
            print('WARNING: hard-coded leap-second dataset is officially valid for dates up to %s for mode "exact". After this date, Astrog extrapolates based on last available leap-second values. To update the dataset, check: https://www.ietf.org/timezones/data/leap-seconds.list'%(NTP_valid.date()))
            leap_sec_5yBefore = leap_sec[np.where(np.array(NTP_date)<NTP_valid-dt.timedelta(days=370*5))[0][-1]]
            dT_TTyear = NTP_valid.year
            dT_TTval  = leap_sec[-1]
            dT_TTinc  = (dT_TTval-leap_sec_5yBefore)/5
            # approximation of dT in requested year (dT = TT-UT1)
            dT_TT_extrapolated = (np.array(32.184) + dT_TTval+dT_TTinc*(np.array(dateIn.year)-dT_TTyear))/86400
            dT_TT[np.array(dateIn.year)-dT_TTyear>0] = dT_TT_extrapolated[np.array(dateIn.year)-dT_TTyear>0]

    else:
        raise Exception('mode=%s not recognized'%(mode_dT))

    return dT_TT


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
    import datetime as dt
    import pandas as pd

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
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    
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


