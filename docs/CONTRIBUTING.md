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
- ``python -m pip install -e .[dev,docs,examples]`` (pip developer mode, any updates to the local folder are immediately available in your python. It also installs all requirements via pip, square brackets are to install optional dependency groups)
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

## Running the testbank

- open anaconda prompt and navigate to hatyan local folder (e.g. ``C:\\DATA\\hatyan``)
- ``conda activate hatyan_hmcenv``
- ``pytest`` (runs all tests)
- ``pytest -m "not acceptance"``
- ``pytest -m acceptance`` (runs the acceptance tests, which are the scripts in [the examples folder](https://github.com/Deltares/hatyan/tree/main/tests/examples))

## Generate documentation (automatically runs via Github Actions upon push to main)

- open anaconda prompt and navigate to hatyan local folder, e.g. ``C:\\DATA\\hatyan``
- ``conda activate hatyan_hmcenv``
```
cp README.md docs
mkdocs build
```

## Increase the hatyan version number

- commit all changes via git
- open anaconda prompt and navigate to hatyan local folder, e.g. ``C:\\DATA\\hatyan``
- ``conda activate hatyan_hmcenv``
- ``bumpversion major`` or ``bumpversion minor`` or ``bumpversion patch``
- the hatyan version number of all relevant files will be updated, as stated in setup.cfg

## Create release

- bump the versionnumber with ``bumpversion minor``
- update ``docs/whats-new.md`` and add a date to the current release heading
- run local testbank
- local check with: ``python -m build`` and ``twine check dist/*`` ([does not work on WCF](https://github.com/pypa/setuptools/issues/4133))
- make sure the remote ``main`` branch is up to date (important issues solved, all pullrequests and branches closed, commit+push all local changes)
- copy the hatyan version from https://github.com/Deltares/hatyan/blob/main/setup.cfg (e.g. ``0.11.0``)
- go to https://github.com/Deltares/hatyan/releases/new
- click ``choose a tag`` and type v+versionnumber (e.g. ``v0.11.0``), click ``create new tag: v0.11.0 on publish``
- set the release title to the tagname (e.g. ``v0.11.0``)
- click `Generate release notes` and replace the `What's Changed` info by a tagged link to ``docs/whats-new.md``
- if all is set, click ``Publish release``
- a release is created and the github action publishes it on PyPI (https://pypi.org/project/hatyan/)
