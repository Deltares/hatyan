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
import pandas as pdvla
import matplotlib.pyplot as plt
plt.close('all')
import hatyan

#file_config = os.path.realpath(__file__)
#dir_output, timer_start = hatyan.init_RWS(file_config, sys.argv, interactive_plots=False)

# script settings
timeStart = dt.datetime(2022,1,15)
timeEnd   = dt.datetime(2022,1,17)
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

#sunriseset_python = hatyan.astrog_moonriseset(tFirst=dt.datetime(2022,1,15), tLast=dt.datetime(2022,1,20), dT_fortran=dT_fortran, tzone='Australia/Sydney', lon=151.2093, lat=-33.8688)
sunriseset_python = hatyan.astrog_moonriseset(tFirst=timeStart, tLast=timeEnd, dT_fortran=dT_fortran, tzone='Vladivostok', lon=135.3084, lat=54.7165)
print(sunriseset_python)


print('moon')
for lon in []:#np.arange(-180,180+1,45): #30 degrees is 30/360*24=2 hours suntime
    print()
    print(lon)
    tz_LOCAL = dt.timezone(dt.timedelta(hours=lon/360*24))
    try:
        sunriseset_python = hatyan.astrog_moonriseset(tFirst=timeStart, tLast=timeEnd, dT_fortran=dT_fortran, tzone=tz_GMT, lon=lon, lat=52.0)
        #print('SUCCES')
        print(sunriseset_python.iloc[:3])
    except Exception as e:
        print(f'FAILED: {e}')
#%%
#hatyan.exit_RWS(timer_start) #provides footer to outputfile when calling this script with python
