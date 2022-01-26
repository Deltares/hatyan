# -*- coding: utf-8 -*-
"""
astrog_RWS.py
voor productie data conform astrog30
Jan Rolf Hendriks
20210107

"""

import os, sys, pytz
import numpy as np
import datetime as dt
import pandas as pd
import matplotlib.pyplot as plt
plt.close('all')
import hatyan

#file_config = os.path.realpath(__file__)
#dir_output, timer_start = hatyan.init_RWS(file_config, sys.argv, interactive_plots=False)

# script settings
timeStart = dt.datetime(2022,1,20)
timeEnd   = dt.datetime(2022,1,27)
dT_fortran = False #True is best comparison to fortran, False is more precise 
tz_GMT = 'UTC' # UTC/GMT timezone
try: #only works with pandas 1.2.0 or higher, but 1.1.5 is highest available pandas version in Python 3.6.12 (which is the highest available Python in RHEL)
    tz_MET = 'UTC+01:00' #for timezone conversion to UTC+01:00 (Dutch wintertime)
    pytz.timezone(tz_MET)
except:
    tz_MET = dt.timezone(dt.timedelta(hours=1)) #for timezone conversion to UTC+01:00 (Dutch wintertime)
tz_EurAms = 'Europe/Amsterdam' # for conversion to local timezone, including daylight saving time (DST)

pdtocsv_kwargs = dict(index=False, sep=';', date_format='%Y%m%d %H%M%S', float_format='%9.5f', na_rep='--:--')

# sunrise and -set
# ook uitgegeven door KNMI voor lon=5, lat=52 (ronde coordinaten, rotonde bij IJsselstein) voor tabel met zon op/onder en begin seizoenen:
# https://cdn.knmi.nl/ckeditor_assets/attachments/170/Tijden_van_zonopkomst_en_-ondergang_2022.pdf

print('sun')
for lon in []:#np.arange(-180,180+1,45): #30 degrees is 30/360*24=2 hours
    print()
    print(lon)
    tz_LOCAL = dt.timezone(dt.timedelta(hours=lon/360*24))
    sunriseset_python = hatyan.astrog_sunriseset(tFirst=timeStart, tLast=timeEnd, dT_fortran=dT_fortran, tzone=tz_GMT, lon=lon, lat=52.0)
    print(sunriseset_python.iloc[-2:])
    
lat,lon =  54.7165, 135.3084 #Vladivostok #this one crashes for longer time periods (rates of increase go off track). Sort of solved by switching signs of ALTMOO RATE in astrac, but moonrise/set are then switched and it crashes for lat,lon=50,45
lat,lon =  50, 45.3084 #fake
lat,lon = -33.8688, 151.2093 #sydney
#lat,lon =  52.1561,   5.3878 #amersfoort
#lat,lon =  51.47869,  -0.01080 #greenwich
sunriseset_python = hatyan.astrog_moonriseset(tFirst=timeStart, tLast=timeEnd, dT_fortran=dT_fortran, tzone='UTC+10:00', lon=lon,lat=lat)
#sunriseset_python = hatyan.astrog_sunriseset(tFirst=timeStart, tLast=timeEnd, dT_fortran=dT_fortran, tzone='UTC', lon=lon,lat=lat)
datesmoonrise = pd.DatetimeIndex(sunriseset_python.loc[:,'datetime'].dt.tz_localize(None))
print(sunriseset_python)
astrabOutput = hatyan.astrab(datesmoonrise,dT_fortran=dT_fortran,lon=lon,lat=lat)
#print(astrabOutput['ALTSUN'])
#print(astrabOutput['ALTMOO'])
#print(astrabOutput['ALTMOO'])

print('moon')
for lon in []:#np.arange(-180,180+1,45): #30 degrees is 30/360*24=2 hours suntime
    print()
    print(lon)
    tz_LOCAL = dt.timezone(dt.timedelta(hours=lon/360*24))
    try:
        lat=50
        tzone_local = dt.timezone(dt.timedelta(hours=lon/360*24))
        sunriseset_python = hatyan.astrog_moonriseset(tFirst=timeStart, tLast=timeEnd, dT_fortran=dT_fortran, tzone=tz_GMT, lon=lon, lat=lat) 
        #print('SUCCES')
        print(sunriseset_python.iloc[:3])
        datesmoonrise = pd.DatetimeIndex(sunriseset_python.loc[:,'datetime'].dt.tz_localize(None))
        astrabOutput = hatyan.astrab(datesmoonrise,dT_fortran=dT_fortran,lon=lon,lat=lat)
        print(astrabOutput['ALTMOO'])
    except Exception as e:
        print(f'FAILED: {e}')
#%%
#hatyan.exit_RWS(timer_start) #provides footer to outputfile when calling this script with python
