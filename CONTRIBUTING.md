Information for developers
--------

Create python environment hatyan_env and install hatyan in developer mode:

- download Anaconda 64 bit Python 3.7 (or higher) from https://www.anaconda.com/distribution/#download-section (miniconda should also be sufficient, but this is not yet tested). Install it with the recommended settings, but check 'add Anaconda3 to my PATH environment variable' if you want to use conda from the windows command prompt instead of anaconda prompt
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

Generate documentation (automatically runs via Github Actions upon push to main):

- open command line and navigate to hatyan local folder, e.g. ``C:\\DATA\\hatyan_github``
- ``conda activate hatyan_env``
- ``python scripts/generate_documentation.py``

Generate RPM (RHEL/CentOS installer, automatically runs via Github Actions upon release creation):

- use the script in scripts/hatyan_rpmbuild.sh (for instance on the CentOS7 Deltares buildserver)
- preparation: activate environment, run testbank and check acceptancetest output, update history.rst, git add+commit, bumpversion minor, (run testbank and) backup acceptancetest output, generate documentation, git add+commit+push, merge branch with main, create tag+release on github (e.g. v2.3.0)
- this script uses the rpmbuild command and the specfile to generate an RPM on a CentOS/RHEL machine with the correct dependencies installed
- rpmbuild uses the specfile scripts/hatyan_python-latest.spec as input (set the versiontag variable to the newly created github tag)
- the dependencies for the RPM are documented in the specfile
- the required Python libraries are documented in requirements_dev.txt: these are fixed versions, which is at least relevant for sip, since it needs to be compatible with pyqt5==5.7.1 for Qt5 plots
- additionally, the library pyqt5==5.7.1 (hardcoded in specfile) is for interative QT5 plots. There is a newer version but it requires glibc >2.14, while 2.12 is the highest version available on CentOS/RedHat 6)
- to test hatyan on CentOS without installing an RPM: use the script scripts/hatyan_rpmbuild_nobinaries.sh, this creates a comparable setup in the home directory and a ~/hatyan_fromhome.sh file comparable to hatyan command. If you get an error about X11-forwarding, first try the xterm command.

Publish to PyPI (automatically runs via Github Actions upon release creation):

- open command line and navigate to hatyan local folder, e.g. ``C:\\DATA\\hatyan_github``
- ``conda activate hatyan_env``
- ``python setup.py sdist bdist_wheel``
- ``twine upload dist/*``
- a new version should now be available on https://pypi.org/manage/project/hatyan/releases/
