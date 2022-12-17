#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

setup(
    name='hatyan',
    version='2.5.93',
    author="Jelmer Veenstra",
    author_email='Jelmer.Veenstra@Deltares.nl',
    url='https://github.com/Deltares/hatyan',
    description="hatyan is a tidal analysis and prediction tool of Rijkswaterstaat",
    long_description=readme,
    long_description_content_type='text/markdown',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Environment :: X11 Applications :: Qt',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Information Analysis',
    ],
    platforms=['platform independent'],
    keywords='hatyan',
    license='LGPL',
    python_requires='>=3.6',
    install_requires=['scipy>=1.5.4', 'numpy>=1.18.4', 'matplotlib>=3.2.1', 'pandas>=1.1.2', 'netCDF4>=1.5.3', 'pyproj>=2.2.0', 'sip>=4.19.8', 'packaging>=21.0', 'requests', 'statsmodels'],
    packages=find_packages(include=['hatyan']),
    test_suite='tests',
    tests_require=['pytest>=5.0.1','bump2version>=0.5.11','pytest-cov','pdoc3','pandoc'], #pyqt5,twine
    include_package_data=True,
    zip_safe=False,
)
