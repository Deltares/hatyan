Information for developers
--------

## Checkout dfm_tools git repository

- this is just a suggestion, feel free to work with VScode or any other git-compatible workflow
- download git from [git-scm.com](https://git-scm.com/download/win), install with default settings
- open git bash window where you want to clone the dfm_tools github repository (e.g. ``C:\DATA\``)
- git clone https://github.com/deltares/hatyan (creates a folder hatyan with the checked out repository)
- ``cd hatyan``
- optional: ``git config --global user.email [emailaddress]``
- optional: ``git config --global user.name [username]``

## Setup local developer environment

- download and install Anaconda 64 bit Python 3.9 (or higher) from [anaconda.com](https://www.anaconda.com/distribution/#download-section) (miniconda should also be sufficient, but this is not yet tested). Install it with the recommended settings.
- open anaconda prompt and navigate to hatyan checkout folder, e.g. ``C:\DATA\hatyan``
- ``conda env create -f environment_hmc.yml`` (installs environment with HMC fixed python/package versions)
- ``conda activate hatyan_env``
- ``python -m pip install -e .[test]`` (pip developer mode, any updates to the local folder are immediately available in your python. It also installs all requirements via pip, ``[test]`` installs also the developer requirements)
- ``conda deactivate``
- to remove hatyan_env when necessary: ``conda remove -n hatyan_env --all``



Increase the hatyan version number:

- open anaconda prompt and navigate to hatyan local folder, e.g. ``C:\\DATA\\hatyan_github``
- ``conda activate hatyan_env``
- ``bumpversion major`` or ``bumpversion minor`` or ``bumpversion patch``
- the hatyan version number of all relevant files will be updated, as stated in setup.cfg

Running the testbank:

- open anaconda prompt and navigate to hatyan local folder, e.g. ``C:\\DATA\\hatyan_github``
- ``conda activate hatyan_env``
- ``pytest`` (runs all tests)
- ``pytest -m unittest``
- ``pytest -m systemtest``
- ``pytest -m acceptance`` (runs the acceptance tests, which are the scripts in [the examples folder](https://github.com/Deltares/hatyan/tree/main/tests/examples))
- ``pytest -m "not acceptance"`` (excludes all acceptance tests)
- the following arguments are automatically provided via pytest.ini: ``-v --tb=short``, add ``--cov=hatyan`` for a coverage summary

Generate documentation (automatically runs via Github Actions upon push to main):

- open anaconda prompt and navigate to hatyan local folder, e.g. ``C:\\DATA\\hatyan_github``
- ``conda activate hatyan_env``
- ``python scripts/generate_documentation.py``

Publish to PyPI (automatically runs via Github Actions upon release creation):

- open anaconda promptand navigate to hatyan local folder, e.g. ``C:\\DATA\\hatyan_github``
- ``conda activate hatyan_env``
- ``python setup.py sdist bdist_wheel``
- to check before uploading: ``twine check dist/*``
- ``twine upload dist/*``
- a new version should now be available on https://pypi.org/manage/project/hatyan/releases/
