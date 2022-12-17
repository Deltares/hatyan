[![pytest-devenv](https://github.com/Deltares/hatyan/actions/workflows/pytest-devenv.yml/badge.svg?branch=main)](https://github.com/Deltares/hatyan/actions/workflows/pytest-devenv.yml)
[![pytest-py39](https://github.com/Deltares/hatyan/actions/workflows/pytest-py39.yml/badge.svg?branch=main)](https://github.com/Deltares/hatyan/actions/workflows/pytest-py39.yml)
[![sigrid-publish](https://github.com/Deltares/hatyan/actions/workflows/sigrid-publish.yml/badge.svg?branch=main)](https://github.com/Deltares/hatyan/actions/workflows/sigrid-publish.yml)
[![pypi-upload](https://github.com/Deltares/hatyan/actions/workflows/pypi-upload.yml/badge.svg?event=release)](https://github.com/Deltares/hatyan/actions/workflows/pypi-upload.yml)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=Deltares_hatyan&metric=alert_status)](https://sonarcloud.io/dashboard?id=Deltares_hatyan)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Deltares/hatyan/HEAD)

# hatyan

hatyan is a Python program for tidal analysis and prediction, based on the FORTRAN version and developed for Rijkswaterstaat.


Information and examples
--------
- [docs folder](https://github.com/Deltares/hatyan/tree/main/docs) with background information
- [online documentation](https://htmlpreview.github.io/?https://github.com/Deltares/hatyan/blob/main/docs/hatyan/index.html) generated from docstrings (replace 'main' in the url with any tagname to view older versions)
- [jupyter notebooks](https://github.com/Deltares/hatyan/blob/main/notebooks) with example code
- [use binder](https://mybinder.org/v2/gh/Deltares/hatyan/HEAD) to run these notebooks interactively (loading takes a while)
- [github folder](https://github.com/Deltares/hatyan/tree/main/tests/examples) with more example scripts


Installation
--------

Install hatyan OPTION 1: Install from pip/github in an existing or new environment:

- optional: download and install Anaconda 64 bit Python 3.7 (or higher) from https://www.anaconda.com/distribution/#download-section
- open anaconda prompt
- optional: ``conda create --name hatyan_env -c conda-forge python=3.7 git spyder -y`` (or higher python version)
- optional: ``conda activate hatyan_env``
- ``python -m pip install hatyan`` (this installs hatyan and all required packages from PyPI, add a version like ``==2.3.0`` if you require a specific version. Optionally add ``--upgrade``)
- alternatively: ``python -m pip install git+https://github.com/Deltares/hatyan`` (this installs hatyan and all required packages from github, add a tag like ``@v2.3.0``, ``@main`` or ``@development`` if you require a specific release/branch. Optionally add ``--upgrade``)

Install hatyan OPTION 2: create python venv in your Linux home directory and install from zipfile (no internet required for installation)

- might be necessary to activate python3: ``module load anaconda3`` (or any other python3 installation)
- ``python -m venv ~/venv_hatyan``
- ``source ~/venv_hatyan/bin/activate``
- optional (requires internet): ``python -m pip install --upgrade pip setuptools``
- ``pip install ~/PyQt5-5.15.7-cp37-abi3-manylinux1_x86_64.whl`` (download from https://pypi.org/project/PyQt5/#files)
- ``pip install ~/hatyan-2.5.64.zip`` (download from https://github.com/Deltares/hatyan/releases/tag/v2.5.64 or any other release)
- ``deactivate``
- to test, copy example script to e.g. home directory (https://github.com/Deltares/hatyan/blob/main/tests/examples/validate_astro_DDL.py)
- ``source ~/venv_hatyan/bin/activate``
- ``python -m venv ~/venv_hatyan``
- ``python ~/validate_astro_DDL.py``
- ``deactivate``

Install hatyan OPTION 3: get and install RPM on CentOS/RHEL (this will be deprecated soon)

- get the latest rpm file (see developer information for building procedure)
- install hatyan on CentOS: ``rpm -i hatyan_python-2.2.30-1.x86_64.rpm``
- upgrade hatyan on CentOS: ``rpm -U hatyan_python-2.2.31-1.x86_64.rpm``
- installing the RPM results in a hatyan command in linux, this activates a Python virtual environment and sets necessary Qt environment variables. It creates a folder with a python environment hatyan_env, doc en tests (/opt/hatyan_python/hatyan_env/) and a file that provides the hatyan command (/usr/bin/hatyan)
- check version: ``hatyan --version``
- test installation: ``hatyan /opt/hatyan_python/tests/examples/predictie_2019_19Ycomp4Ydia_VLISSGN_interactive.py`` (or use the ``hatyan --test`` shortcut)
- this should result in several interactive figures popping up, described in chapter 5 (Quick start guide) of the hatyan user manual (gebruikershandleiding).
- if you see the message "RuntimeError: Invalid DISPLAY variable", restart the MobaXterm connection and try again.
- the followning warning can be ignored: "QXcbConnection: XCB error: 145 (Unknown), sequence: 171, resource id: 0, major code: 139 (Unknown), minor code: 20". To avoid it, disable the extension RANDR in Mobaxterm settings (Settings > Configuration > X11)


