Information for developers
--------

Create python environment hatyan_env and install hatyan in developer mode:

- download and install Anaconda 64 bit Python 3.8 (or higher) from https://www.anaconda.com/distribution/#download-section (miniconda should also be sufficient, but this is not yet tested). Install it with the recommended settings, but check 'add Anaconda3 to my PATH environment variable' if you want to use conda from the windows command prompt instead of anaconda prompt
- download git from https://git-scm.com/download/win, install with default settings
- create a branch called work_yourname on https://github.com/Deltares/hatyan
- open git bash window where you want to clone the hatyan github repository (e.g. C:\\DATA\\)
- optional: ``git config --global user.email [emailaddress]``
- optional: ``git config --global user.name [username]``
- optional: ``git remote update origin --prune`` (update local branch list)
- ``git clone -b work_yourname https://github.com/Deltares/hatyan hatyan_github`` (repo gets cloned in C:\\DATA\\hatyan_github, this is a checkout of the work_yourname branch)
- update your branch if main has been updated: add+commit+push everything in branch first, ``git checkout main``, ``git pull``, ``git checkout development``, ``git merge main -m ''``, ``git push``
- open anaconda prompt and navigate to hatyan local folder, e.g. ``C:\\DATA\\hatyan_github``
- ``conda env create -f environment.yml`` (This yml file installs Python 3.6.12 since that is the latest available Python on RHEL6)
- ``conda info --envs`` (should show hatyan_env virtual environment in the list)
- ``conda activate hatyan_env``
- ``python -m pip install -e .`` (pip developer mode, also install all packages in requirements.txt containing CentOS tested libraries, linked via setup.py) >> maybe add ``test`` to install also test requirements [like this](https://stackoverflow.com/questions/15422527/best-practices-how-do-you-list-required-dependencies-in-your-setup-py)
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
