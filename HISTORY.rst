=======
History
=======

* Thu Dec 14 2021 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.4.4
- bugfix in diawriting for last entries of timeseries
- bugfix in diafile writing that solves rounding issue. Made dia writing more efficient by writing values with pd.to_csv instead of in loop
- added DDL to hatyan (might be merged with ddlpy package on github later)

* Fri Oct 22 2021 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.3.3
- improvement in get_components_from_ts() to support analysis_peryear=True and CS_comps=!None at the same time. Now possible to remove stats_anaperyear from configfiles
- updated example in readme to work with online UHSLC file
- added prediction_peryear() definition

* Fri Oct 8 2021 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.3.0
- simplified rpmbuild scripts and specfile, make more directly dependent on github

* Thu Oct 7 2021 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.2.89
- generalized imports, now only 'import hatyan' necessary, all modules can be loaded from there

* Mon Oct 4 2021 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.2.88
- bugfix for foreman frequencies, were actually retrieved from schureman

* Fri Oct 1 2021 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.2.87
- fix for non-backwards compatible change in pandas, to make it work in pandas versions both older and newer than 1.2.0

* Wed Sep 29 2021 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.2.85
- first version of write_tsnoos() definition added
- fixed extreme calculation for Lith Dorp, the algorithm is now more robust and less sensitive
- fixed agger calculation for Lith Dorp

* Thu Sep 9 2021 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.2.81
- added foreman to analysis(), get_components_from_ts() and prediction(), and added testcase and extended configfile for schureman/foreman comparison
- Converted code for multiindex multiyear component dataframe from 4 to 1 line.
- added saveguard for duplicate indices when using Timeseries.resample_timeseries(). improved error message
- Simplified Timeseries.resample_timeseries() and extended it with tstart/tstop arguments.

* Thu Jul 22 2021 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.2.78
- Made hatyan analysis/prediction more robust by also supporting times outside pandas.date_range extremes (years 1677 to 2262)
- updated user error message in Timeseries.crop_timeseries()

* Sat Jul 17 2021 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.2.77
- astrog fix of periodicity bug, by removing -1 from range-loops (results in extra precision). Also updated expected values in systemtests.
- cleaned up small things in astrog code
- improved astrog_test.py by providing pdtocsv_kwargs as dict and made output more readable (timezones and commas)
- updated metadata in setup.py based on https://docs.python.org/3/distutils/setupscript.html
- added scripts/generate_documentation.py to copy contents of README.md to hatyan/__init__.py automatically, followed by generation of html+pdf documentation. Updated (reduced) README.md text about documentation generation.
- updated html/pdf documentation

* Thu Jul 15 2021 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.2.75
- major improvements of astrog definitions and interactions, cleaned up code, added internal timezones (default='UTC'), fixed moonphases bug, made leapsecond default choice, made datasets/csvwriting/plotting more generic
- astrog_test.py is now more readable, so can be easily used by users
- updated history.rst, updated astrog docstrings, updated documentation

* Fri Jul 2 2021 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.2.70
- updated documentation, history.rst, description of timeseries DataFrame in docstrings
- B&O: cleaned get_diablocks and made retrieval a bit simpler, added netcdf read to testbank (just as an example), added documentation for Foreman and a reference in the foreman_core.py docstring, printed more parameter...
- added dir_testdata P-drive in comment in all configfiles
- replaced continue_process with continue in master configfiles
- updated numbering_FEWS_PG.py (replaced dfm_tools with netCDF4)
- cleaned testbank, splitted unitsystem and acceptance testdata in order to reduce repos size and updated configfiles accordingly. Removed acceptance testdata from repos
- Removed some duplicate codelines from foreman_core.py and astrog.py
- added fft analysis and translation to to hatyan

* Fri Jun 11 2021 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.2.64
- Added MSQM to hatyan_core.py
- code/configfile overhaul, ts dataframes now have timestamps as index

* Tue May 4 2021 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.2.62
- improved rayleigh treshold feedback
- converted astrog to Astrog in testbank
- reduced amount of data in 'testdata_unitsystemtests' folder
- Removed foreman.py, old version of foreman_core.py
- Some astrog improvements (more general)
- removed as_ex_nld.dia, since it was not used
- updated docstring for calc_HWLWnumbering() and updated documentation accordingly
- removed --use-feature=in-tree-build from specfile again (was just as test)
- updated readme and documentation

* Fri Apr 30 2021 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.2.59
- updated html and pdf documentation
- updated documentation
- updated README.md and __init__.py docstring
- specfile rollback to requirements_dev.txt
- --use-feature=in-tree-build  in specfile to test with new pip version
- updated pkl inclusion in MANIFEST.in
- Renamed foreman and hatyan data to data_components_*.*, renamed phasediff dataset to data_M2phasediff_perstation.txt, updated all relevant files. test: used requirements.txt in specfile (instead of requirements_dev.tx...
- replaced _middenstanden_predictie2019.txt by _slotgemiddelden_predictie2019.txt in configfiles
- updated readme.md to contain more installation/building info for RPM (removed from user manual)
- added pyqt5 installation in specfile (was via requirements_CentOS.txt previously, now using requirements_dev.txt)
- included pkl file via MANIFEST.in instead of setup.py (2nd test)
- added pkl file to setup.py (include), as a test
- updated specfile: requirements_CentOS.txt to requirements_dev.txt
- added writing of figure to configfiles/predictie_2019_frommergedcomp_WSdwarsstroming_test.py

* Mon Apr 12 2021 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.2.56
- Added EPS2 component to hatyan_core and data_components.pkl (same freq/v0 as MNS2 and same u/f as M2)
- updated testbank and configfiles accordingly
- restructured data folder
- merged test_analysis_settings and test_analysis
- moved station_M2phasediff.txt to hatyan code folder, splitted configfiles/acceptancetests from main testbank script
- redirected test_hatyan_main.py to new testdata_predictie2019 folder and removed *_all.py testcases.
- removed *_all.py configfiles and changed testdata_predictie2019 location
- moved testdata_predictie2019 to separate folder
- removed rest of testdata_analysis and analysis_components_test.py
- added ``-v --tb=short`` arguments to pytest.ini
- interactive plots back to False
- bugfix in configfile export_freq_v0uf_data_test.py, v0 and v0+u difference plots are now 0 straight lines

* Sat Apr 3 2021 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.2.55
- removed analysis_components_all.py
- added analysis testcase
- removed faulty datasets
- improved hatyan55 and v0uf2016 file validation data to export_freq_v0uf_data_test.py figures, toned down analysis_components_test.py (removed faulty datasets)
- corrected foreman shallow water relation for 2MSN4
- added hatyan55 and v0uf2016 file validation data to export_freq_v0uf_data_test.py figures
- bugfix in plot_components (now also possible to plot nonexistent components)
- bugfix in plot_components, diff is now between -180 and 180 instead of 0 and 360
- cleaned up commented code in hatyan/components.py

* Wed Mar 31 2021 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.2.53
- Improved components_plot (sorting and difference now better implemented). Added timeshift for component set. Added test for available constituents to hatyan_core.py (with new get_v0uf_sel() definition)
- small updates in configfile
- renamed foreman.py to foreman_core.py, cleaned up a bit and replaced a for-loop with matrices.
- merged export_freq_v0uf_data_test.py and foreman_test.py and improved output

* Fri Mar 12 2021 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.2.52
- Removed commented parts of code throughout entire code
- cleaned up foreman file, small updates on configfiles
- foreman: Z0 replaced by A0 and made more generic (now works for v0freq as well as uf). Made frequency calculation the default (over reading from foreman file) and removed some bugs there (now higher accuracy and more in line with v0 calcula...

* Thu Mar 11 2021 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.2.51
- updated and improved foreman.py and foreman_test.py to pandas and arrays, more efficient and better usable.
- Put foreman in main folder again, repaired everything, made it faster and (started with) removing unnecessary parts of code, updated header conform other hatyan scripts
- foreman weer werkend krijgen, alle scripts gecheckt en sneller gemaakt (dood_date als array ipv loop over losse datums)
- hatyan BO: updated export_freq_v0uf_data_test.py configfile to way more efficient and more functionalities (but less unnecessary output)
- removed waterkeringen_normtrajecten_20160613.ldb from datafolder since it is not used
- made plot_HWLW_validatestats() more robust (fail with try/except instead of crash)

* Wed Mar 10 2021 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.2.49
- General B&O: removed some duplicate code from several definitions
- wrapper_RWS.py: replaced FILE_CONFIG.txt' by %s'%(os.path.basename(file_config))
- merged export_doodnum_test.py in export_freq_v0uf_data_test.py

* Tue Mar 9 2021 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.2.48
- Removed old+unused+slow HWLW statistics calculation in timeseries.py(). Improved unique timestep calculation in check_ts() definition.

* Tue Mar 9 2021 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.2.47
- It is now possible to read diafiles that contain multiple diablocks for one station (and append and sort them automatically), this was an issue at the kenmerkende waarden project but has now been taken care of. Due to this change, the entire dia related code also had to be better structured, hopefully also improving the SIG score. Testbank is extended (read multi diafile, multi diablock) and testconfigfiles are updated to work with the new code. Updated all readts_dia_HWLW() to readts_dia()
- improved dia reading in timeseries.py, made more efficient and less duplicated (was selected due to SIG violations)
- cleaned up Rayleigh prints in analysis_prediction.py, added comments to hatyan_core.py
- bugfix for strptime in astrog (added .datetime and proper string value)
- removed unnecessary line of code
- switched order of N and P1 doodson numbers again to correspond with 'normal' order
- added comments to hatyan_core.py

* Thu Feb 25 2021 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.2.43
- Added MA2, MB2 and alternatives for SA and S1 (for research purposes). systemtests still all pass. updated data_components.pkl file
- added fstr column, but not functional in f definition, testbank does work again.
- added comments in hatyan_core.py with differences compared to IHO/SLS
- Tested addition of N column (works) and added option to recalculate v0uf_all table in hatyan_core.py. also added lunar table for comparison.
- added extra stations (including -360 to testscript, not updated M2phasediff file yet)
- switched order of P1 and N. Removed N column from pkl file, since it has no function
- added N4 in data_components.pkl

* Wed Jan 20 2021 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.2.30
- cleaned up several parts of code and testbank, based on SIG score/points
- added number_HWLW() function and improved calc_HWLW, including testcases
- added HWLW 345 code calculation and numbering
- added first version of astrog (moonphases and such)
- updated documentation (pdf to 2.2.28, html to 2.2.30)
- made RPM creation more generic (on Deltares buildserver), created some necessary scripts in hatyan_python/scripts folder instead of copy-paste codelines from manual/readme
- removed VM files from repos, renamed RWS folder to build

* Thu Oct 28 2020 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.2.22
- added README.md and the same text to __init__.py (includes example code, installation guide and developer information that were previously in user manual)
- restructured specfile to make update via zip possible and make installation more according to standards

* Thu Sep 18 2020 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.2.20
- reprogrammed extremes calculation, with boolean for local extremes output
- bugfix in netCDF extremes writing
- added catch for singular matrix

* Thu Jul 23 2020 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.2.16
- renamed RPM from hatyan to hatyan_python, command stays hatyan
- H0 as component instead of separate

* Wed Jul 22 2020 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.2.12
- hatyan and venv are now moved to /opt/hatyan_python, since /opt/hatyan is occupied by fortran hatyan. name of program is still hatyan.
- removed readts_mat and corresponding data
- moved get_outputfoldername() to wrapper_RWS.py
- added more documentation to docstrings of hatyan functions
- improved components dataframe, for easier sorting and differences
- better error for singular matrix
- now phi_deg instead of phi_rad

* Fri May 22 2020 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.2.0
- changed name from hatyan2 to hatyan, increased version to hatyan-2.2.0 (previous official release was hatyan2-1.0)
- restructured and slimmed down testbank
- slimmed down datafiles in RPM
- completed overhaul to new direct call instead of old configfiles

* Tue May 19 2020 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.1.8.1
- cleaned svn structure, moved settings files to data folders
- moved vectoravg() outside of get_components_from_ts(), in order to remove get_components_from_ts() in the future (and ts_ids, ts_years)
- extended bumpversion to also update version numbers in RWS/hatyan-rpm.spec and RWS/hatyan_commands.sh file
- converted hatyan to new interaction (configfiles converted to python scripts that call hatyan)
- replaced Timeseries and Components classes with pandas DataFrame

* Fri Mar 15 2020 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.1.4
- made requirements.txt more flexible, but hardcoded matplotlib, pyqt5 and sip files in spec file to avoid "ImportError: Failed to import any qt binding" and "ImportError: Cannot load backend 'Qt5Agg' which requires the 'qt5' interactive framework, as 'headless' is currently running" and "TypeError: float() argument must be a string or a number, not 'Timestamp'"

* Fri Feb 5 2020 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.1.3
- RPM's merged into one (code and venv), venv is now moved to /opt/hatyan_python/env/
- dependencies for hatyan code are now installed via setup.py>>requirements.txt (pip install -e hatyan)

* Fri Feb 4 2020 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.1.2
- includes post and preun added, for pip install of python program (no internet should be required)
- upgrade for pyproj to 1.9.6, since windows venv did not support 1.9.5.1
- added pytest==5.0.1 pytest-cov teamcity-messages for testbank

* Fri Jan 3 2020 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.1.1
- new name for python environment (hatyan_venv instead of hatyan_py_env), it also which fixes more software versions and contains netcdf
- upgrade pip in the building process
- pip install sip==4.19.8 toegevoegd, met (automatisch) nieuwere versie hiervan of van dependencies was koppeling naar qt niet meer goed.
- netCDF4==1.5.3 toegevoegd
- made rh-python36 version dependency minimal instead of fixed

* Thu Jun 20 2019 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.1.0
- bugfix in component plotting
- finalised component splitting, now correct and more robust
- fix in dia-file output format, now compatible with DONAR
- added test for minimial dia-inputfile contents, including coordinate check
- updated component output file
- replaced LDA2 to LABDA2 (removed exception and replacement), to avoid errors
- improved spatial summary programming, incl ldb coordinate conversion. WGS84 and RD supported
- added vertical reference checks, file_station checks icm stations_strict setting. added testcases with wrong data from koos to show the functionality.
- implemented block read for dia files, more structured and stable. Also makes it possible to select a specific block from a file and prevents reading in wrong data.
- renamed hatyan_py to hatyan
- final release for RWS for June 2019

* Fri May 1 2019 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.0.10
- added component splitting
- added all necessary functionality
- added numerous configfiles for almost all 121 donar stations
- moved all individual script tests to configfiles
- added spatial summary plotting functionality with coordinate conversion
- restructured code and made more stable
- pre-final release for RWS for 1 June 2019

* Fri Aug 17 2018 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.0.7
- better and more output written to output (screen and file), to facilitate debugging
- added expect package to requirements, facilitates line-buffered tee-output instead of blocks
- replaced component numbers by names in figures
- catch hiaat-values in dia files (999999999/99) and replace by nan

* Fri Jun 22 2018 Jelmer Veenstra <jelmer.veenstra@deltares.nl> 2.0.6
- final first RWS RPM, delivered and installed in June 2018