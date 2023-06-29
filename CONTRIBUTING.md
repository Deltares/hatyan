Information for developers
--------

## Checkout hatyan git repository

- this is just a suggestion, feel free to work with VScode or any other git-compatible workflow
- download git from [git-scm.com](https://git-scm.com/download/win), install with default settings
- open git bash window where you want to clone the hatyan github repository (e.g. ``C:\DATA\``)
- git clone https://github.com/deltares/hatyan (creates a folder hatyan with the checked out repository)
- ``cd hatyan``
- optional: ``git config --global user.email [emailaddress]``
- optional: ``git config --global user.name [username]``

## Setup local developer environment

- download and install Anaconda 64 bit Python 3.9 (or higher) from [anaconda.com](https://www.anaconda.com/distribution/#download-section) (miniconda should also be sufficient, but this is not yet tested). Install it with the recommended settings.
- open anaconda prompt and navigate to hatyan checkout folder, e.g. ``C:\DATA\hatyan``
- ``conda create --name hatyan_hmcenv python=3.8.13 git spyder -y`` (``git`` and ``spyder``, python version should be the one at HMC)
- ``conda activate hatyan_hmcenv``
- ``pip install -r environment_hmc.txt`` (installs fixed python/package versions like on HMC)
- ``python -m pip install -e .[test]`` (pip developer mode, any updates to the local folder are immediately available in your python. ``[test]`` installs also the developer requirements)
- ``conda deactivate``
- to remove hatyan_hmcenv when necessary: ``conda remove -n hatyan_hmcenv --all``

## Contributing

- open an existing issue or create a new issue at https://github.com/Deltares/hatyan/issues
- create a branch via ``Development`` on the right. This branch is now linked to the issue and the issue will be closed once the branch is merged with main again
- open git bash window in local hatyan folder (e.g. ``C:\DATA\hatyan``)
- ``git fetch origin`` followed by ``git checkout [branchname]``
- make your local changes to the hatyan code (including docstrings and unittests for functions), after each subtask do ``git commit -am 'description of what you did'`` (``-am`` adds all changed files to the commit)
- check if all edits were committed with ``git status``, if there are new files created also do ``git add [path-to-file]`` and commit again
- ``git push`` to push your committed changes your branch on github
- open a pull request at the branch on github, there you can see what you just pushed and the automated checks will show up (testbank and code quality analysis).
- optionally make additional local changes (+commit+push) untill you are done with the issue and the automated checks have passed
- optionally increase the hatyan version with: ``bumpversion patch``
- request a review on the pull request
- after review, squash+merge the branch into main

## Increase the hatyan version number

- open anaconda prompt and navigate to hatyan local folder, e.g. ``C:\\DATA\\hatyan``
- ``conda activate hatyan_hmcenv``
- ``bumpversion major`` or ``bumpversion minor`` or ``bumpversion patch``
- the hatyan version number of all relevant files will be updated, as stated in setup.cfg

## Running the testbank

- open anaconda prompt and navigate to hatyan local folder, e.g. ``C:\\DATA\\hatyan``
- ``conda activate hatyan_hmcenv``
- ``pytest`` (runs all tests)
- ``pytest -m unittest``
- ``pytest -m systemtest``
- ``pytest -m acceptance`` (runs the acceptance tests, which are the scripts in [the examples folder](https://github.com/Deltares/hatyan/tree/main/tests/examples))
- ``pytest -m "not acceptance"`` (excludes all acceptance tests)
- the following arguments are automatically provided via pytest.ini: ``-v --tb=short``, add ``--cov=hatyan`` for a coverage summary

## Generate documentation (automatically runs via Github Actions upon push to main)

- open anaconda prompt and navigate to hatyan local folder, e.g. ``C:\\DATA\\hatyan``
- ``conda activate hatyan_hmcenv``
- ``python scripts/generate_documentation.py``

## Publish to PyPI (automatically runs via Github Actions upon release creation)

- open anaconda promptand navigate to hatyan local folder, e.g. ``C:\\DATA\\hatyan``
- ``conda activate hatyan_hmcenv``
- ``python setup.py sdist bdist_wheel``
- to check before uploading: ``twine check dist/*``
- ``twine upload dist/*``
- a new version should now be available on https://pypi.org/manage/project/hatyan/releases/
