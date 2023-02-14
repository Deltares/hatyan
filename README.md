[![generate-documentation](https://github.com/Deltares/hatyan/actions/workflows/generate-documentation.yml/badge.svg)](https://github.com/Deltares/hatyan/actions/workflows/generate-documentation.yml)
[![pytest-hmcenv](https://github.com/Deltares/hatyan/actions/workflows/pytest-hmcenv.yml/badge.svg?branch=main)](https://github.com/Deltares/hatyan/actions/workflows/pytest-devenv.yml)
[![pytest-py39](https://github.com/Deltares/hatyan/actions/workflows/pytest-py39.yml/badge.svg?branch=main)](https://github.com/Deltares/hatyan/actions/workflows/pytest-py39.yml)
[![sigrid-publish](https://github.com/Deltares/hatyan/actions/workflows/sigrid-publish.yml/badge.svg?branch=main)](https://github.com/Deltares/hatyan/actions/workflows/sigrid-publish.yml)
[![pypi-upload](https://github.com/Deltares/hatyan/actions/workflows/pypi-upload.yml/badge.svg?event=release)](https://github.com/Deltares/hatyan/actions/workflows/pypi-upload.yml)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=Deltares_hatyan&metric=alert_status)](https://sonarcloud.io/dashboard?id=Deltares_hatyan)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Deltares/hatyan/HEAD)


# hatyan

Abca Python package for harmonic tidal analysis and prediction, based on the FORTRAN version and developed for Rijkswaterstaat. Hatyan contains the methods to derive water level extremes (high and low waters) and several other water level indicators (Kenmerkende Waarden). Furthermore, hatyan provides easier access to Rijkswaterstaat data via their data distribution layer (DataDistributieLaag, DDL).


Information and examples
--------
- [docs folder](https://github.com/Deltares/hatyan/tree/main/docs) with background information
- [online documentation](https://htmlpreview.github.io/?https://github.com/Deltares/hatyan/blob/main/docs/hatyan/index.html) generated from docstrings (replace 'main' in the url with any tagname to view older versions)
- [jupyter notebooks](https://github.com/Deltares/hatyan/blob/main/notebooks) with example code
- [use binder](https://mybinder.org/v2/gh/Deltares/hatyan/HEAD) to run these notebooks interactively (loading takes a while)
- [github folder](https://github.com/Deltares/hatyan/tree/main/tests/examples) with more example scripts


Installation
--------

- optional: download and install Anaconda 64 bit Python 3.8 (or higher) from https://www.anaconda.com/distribution/#download-section
- open anaconda prompt
- optional: ``conda create --name hatyan_env -c conda-forge python=3.8 git spyder -y`` (or higher python version)
- optional: ``conda activate hatyan_env``
- ``python -m pip install hatyan`` (this installs hatyan and all required packages from PyPI, add a version like ``==2.3.0`` if you require a specific version. Optionally add ``--upgrade``)
- alternatively: ``python -m pip install git+https://github.com/Deltares/hatyan`` (this installs hatyan and all required packages from github, add a tag like ``@v2.3.0`` if you require a specific release/branch.)
