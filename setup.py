#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='hatyan',
    version='2.5.42',
    author="Jelmer Veenstra",
    author_email='Jelmer.Veenstra@Deltares.nl',
    url='https://repos.deltares.nl/repos/lib_tide/trunk/src/hatyan_python',
    description="hatyan is a tidal analysis and prediction tool of Rijkswaterstaat",
    long_description=readme,
    long_description_content_type='text/markdown',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Environment :: X11 Applications :: Qt',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Information Analysis',
    ],
    platforms=['platform independent'],
    keywords='hatyan',
    license='LGPL',
    python_requires='>=3.6',
    install_requires=requirements,
    packages=find_packages(include=['hatyan']),
    include_package_data=True,
    zip_safe=False,
)
