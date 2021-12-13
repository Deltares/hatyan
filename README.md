[![pytest-devenv](https://github.com/Deltares/hatyan/actions/workflows/pytest-devenv.yml/badge.svg?branch=main)](https://github.com/Deltares/hatyan/actions/workflows/pytest-devenv.yml)
[![pytest-py39](https://github.com/Deltares/hatyan/actions/workflows/pytest-py39.yml/badge.svg?branch=main)](https://github.com/Deltares/hatyan/actions/workflows/pytest-py39.yml)
[![sigrid-publish](https://github.com/Deltares/hatyan/actions/workflows/sigrid-publish.yml/badge.svg?branch=main)](https://github.com/Deltares/hatyan/actions/workflows/sigrid-publish.yml)
[![rpm-build-core](https://github.com/Deltares/hatyan/actions/workflows/rpm-build-core.yml/badge.svg?event=release)](https://github.com/Deltares/hatyan/actions/workflows/rpm-build-core.yml)

# hatyan

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


Installation
--------

Install hatyan OPTION 1: Install from github in an existing or new environment:

- optional: download Anaconda 64 bit Python 3.7 (or higher) from https://www.anaconda.com/distribution/#download-section (miniconda should also be sufficient, but this is not yet tested). Install it with the recommended settings, but check 'add Anaconda3 to my PATH environment variable' if you want to use conda from the windows command prompt instead of anaconda prompt
- open command window (or anaconda prompt)
- optional: ``conda create --name hatyan_env -c conda-forge python=3.7 git spyder -y`` (or higher python version)
- optional: ``conda activate hatyan_env``
- ``python -m pip install git+https://github.com/Deltares/hatyan`` (this command installs hatyan and all required packages, add a tag like ``@v2.3.0`` if you require a specific version)
- to update hatyan: ``python -m pip install --upgrade git+https://github.com/Deltares/hatyan``
- optional: ``conda deactivate``

Install hatyan OPTION 2: get and install RPM on CentOS/RHEL

- get the latest rpm file (see developer information for building procedure)
- install hatyan on CentOS: ``rpm -i hatyan_python-2.2.30-1.x86_64.rpm``
- upgrade hatyan on CentOS: ``rpm -U hatyan_python-2.2.31-1.x86_64.rpm``
- installing the RPM results in a hatyan command in linux, this activates a Python virtual environment and sets necessary Qt environment variables. It creates a folder with a python environment hatyan_env, doc en tests (/opt/hatyan_python/hatyan_env/) and a file that provides the hatyan command (/usr/bin/hatyan)
- check version: ``hatyan --version``
- test installation: ``hatyan /opt/hatyan_python/tests/examples/predictie_2019_19Ycomp4Ydia_VLISSGN_interactive.py`` (or use the ``hatyan --test`` shortcut)
- this should result in several interactive figures popping up, described in chapter 5 (Quick start guide) of the hatyan user manual (gebruikershandleiding).
- if you see the message "RuntimeError: Invalid DISPLAY variable", restart the MobaXterm connection and try again.
- the followning warning can be ignored: "QXcbConnection: XCB error: 145 (Unknown), sequence: 171, resource id: 0, major code: 139 (Unknown), minor code: 20". To avoid it, disable the extension RANDR in Mobaxterm settings (Settings > Configuration > X11)


Getting started
--------

[Documentation is available on Github](https://htmlpreview.github.io/?https://github.com/Deltares/hatyan/blob/main/doc/hatyan/index.html) (replace 'main' in the url with any tagname to view older versions) and there is background information in [the doc folder](https://github.com/Deltares/hatyan/tree/main/doc). Copy the code below to your own script to get started. For more examples, check [the examples folder](https://github.com/Deltares/hatyan/tree/main/tests/examples).

```python
import datetime as dt
import pandas as pd
from netCDF4 import Dataset, num2date
import hatyan
hatyan.close('all')

#defining a list of the components to be analysed (can also be 'half_year' and others, 'year' contains 94 components and the mean H0)
const_list = hatyan.get_const_list_hatyan('year')

#reading and editing time series, results in a pandas DataFrame a 'values' column (water level in meters) and a pd.DatetimeIndex as index (timestamps as datetime.datetime)
file_data_meas = 'http://uhslc.soest.hawaii.edu:80/opendap/rqds/global/hourly/h825a.nc' #Cuxhaven dataset from UHSLC database #os.path.join(r'n:\\Deltabox\\Bulletin\\veenstra\\VLISSGN_waterlevel_20101201_20140101.noos')
times_ext = [dt.datetime(2017,1,1),dt.datetime(2018,12,31)]
timestep_min = 10
ts_data = Dataset(file_data_meas)
ts_data_values = ts_data['sea_level'][0,-18000:]/1000-5 #correct from mm to meters and for 5m offset
ts_data_times = num2date(ts_data['time'][-18000:],units=ts_data['time'].units, only_use_cftime_datetimes=False)
ts_meas = pd.DataFrame({'values':ts_data_values},index=ts_data_times)
#ts_meas = hatyan.resample_timeseries(ts_meas, timestep_min=60) #resampling only works well when timesteps are rounded to seconds
ts_meas = hatyan.crop_timeseries(ts=ts_meas, times_ext=times_ext)

#tidal analysis and plotting of results. commented: saving figure  
comp_frommeas, comp_allyears = hatyan.get_components_from_ts(ts=ts_meas, const_list=const_list, nodalfactors=True, return_allyears=True, fu_alltimes=True, analysis_peryear=True)
fig,(ax1,ax2) = hatyan.plot_components(comp=comp_frommeas, comp_allyears=comp_allyears)
#fig.savefig('components.png')

#tidal prediction and plotting of results. commented: saving figure and writing to netCDF
ts_prediction = hatyan.prediction(comp=comp_frommeas, nodalfactors=True, fu_alltimes=True, times_ext=times_ext, timestep_min=timestep_min)
fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction, ts_validation=ts_meas)
ax1.legend(['prediction','measurement','difference','mean of prediction'])
ax2.set_ylim(-0.5,0.5)
#fig.savefig('prediction.png')

#calculation of HWLW and plotting of results. commented: saving figure
ts_ext_meas = hatyan.calc_HWLW(ts=ts_meas)
ts_ext_prediction = hatyan.calc_HWLW(ts=ts_prediction)
fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction, ts_validation=ts_meas, ts_ext=ts_ext_prediction, ts_ext_validation=ts_ext_meas)
ax1.set_xlim([dt.datetime(2018,6,1),dt.datetime(2018,7,1)])
ax2.set_ylim(-1,1)
#fig.savefig('prediction_HWLW.png')

fig, ax = hatyan.plot_HWLW_validatestats(ts_ext=ts_ext_prediction, ts_ext_validation=ts_ext_meas)
#fig.savefig('prediction_HWLW_validatestats.png')
#hatyan.write_tsnetcdf(ts=ts_prediction, ts_ext=ts_ext_prediction, station='Cuxhaven', vertref='MSL', filename='prediction.nc')
```

Information for developers
--------

Create python environment hatyan_env and install hatyan in developer mode:

- download Anaconda 64 bit Python 3.7 (or higher) from https://www.anaconda.com/distribution/#download-section (miniconda should also be sufficient, but this is not yet tested)
- install it with the recommended settings, but check 'add Anaconda3 to my PATH environment variable' if you want to use conda from the windows command prompt instead of anaconda prompt
- Download git from https://git-scm.com/download/win, install with default settings
- create a branch called work_yourname on https://github.com/Deltares/hatyan
- open command window in a folder where you want to clone the hatyan github repo, e.g. C:\\DATA
- open git bash window where you want to checkout (e.g. C:\\DATA\\)
- ``git remote update origin --prune`` (update local branch list)
- ``git clone -b work_yourname https://github.com/Deltares/hatyan hatyan_github`` (repos gets cloned in C:\\DATA\\hatyan_github, this is a checkout of the work_yourname branch)
- update your branch if main has been updated: add+commit+push everything in branch first, ``git checkout main``, ``git pull``, ``git checkout development``, ``git merge main -m ''``, ``git push``
- open command line and navigate to hatyan local folder, e.g. ``C:\\DATA\\hatyan_github``
- ``conda env create -f environment.yml`` (This yml file installs Python 3.6.12 since that is the latest available Python on RHEL6)
- ``conda info --envs`` (should show hatyan_env virtual environment in the list)
- ``conda activate hatyan_env``
- ``python -m pip install -e . -r requirements_dev.txt`` (pip developer mode, also install all packages in requirements_dev.txt containing CentOS tested libraries, linked via setup.py)
- ``conda deactivate``
- to remove hatyan_env when necessary: ``conda remove -n hatyan_env --all``

Increase the hatyan version number:

- open command line and navigate to hatyan local folder, e.g. ``C:\\DATA\\hatyan_github``
- ``conda activate hatyan_env``
- ``bumpversion major`` or ``bumpversion minor`` or ``bumpversion patch``
- the hatyan version number of all relevant files will be updated, as stated in setup.cfg

Running the testbank:

- open command line and navigate to hatyan local folder, e.g. ``C:\\DATA\\hatyan_github``
- ``conda activate hatyan_env``
- ``pytest`` (runs all tests)
- ``pytest -m unittest``
- ``pytest -m systemtest``
- ``pytest -m acceptance`` (runs the acceptance tests, which are the scripts in [the examples folder](https://github.com/Deltares/hatyan/tree/main/tests/examples))
- ``pytest -m "not acceptance"`` (excludes all acceptance tests)
- the following arguments are automatically provided via pytest.ini: ``-v --tb=short``, add ``--cov=hatyan`` for a coverage summary

Generate documentation:

- open command line and navigate to hatyan local folder, e.g. ``C:\\DATA\\hatyan_github``
- ``conda activate hatyan_env``
- ``python scripts/generate_documentation.py``

Generate RPM (RHEL/CentOS installer):

- use the script in scripts/hatyan_rpmbuild.sh (for instance on the CentOS7 Deltares buildserver)
- preparation: activate environment, run testbank and check acceptancetest output, update history.rst, git add+commit, bumpversion minor, (run testbank and) backup acceptancetest output, generate documentation, git add+commit+push, merge branch with main, create tag+release on github (e.g. v2.3.0)
- this script uses the rpmbuild command and the specfile to generate an RPM on a CentOS/RHEL machine with the correct dependencies installed
- rpmbuild uses the specfile scripts/hatyan_python-latest.spec as input (set the versiontag variable to the newly created github tag)
- the dependencies for the RPM are documented in the specfile
- the required Python libraries are documented in requirements_dev.txt: these are fixed versions, which is at least relevant for sip, since it needs to be compatible with pyqt5==5.7.1 for Qt5 plots
- additionally, the library pyqt5==5.7.1 (hardcoded in specfile) is for interative QT5 plots. There is a newer version but it requires glibc >2.14, while 2.12 is the highest version available on CentOS/RedHat 6)
- to test hatyan on CentOS without installing an RPM: use the script scripts/hatyan_rpmbuild_nobinaries.sh, this creates a comparable setup in the home directory and a ~/hatyan_fromhome.sh file comparable to hatyan command. If you get an error about X11-forwarding, first try the xterm command.
