# What's new

## UNRELEASED

### Feat
- `datetime.timezone` support in `hatyan.write_components()` to support components from ddlpy timeseries in [#359](https://github.com/Deltares/hatyan/pull/359)
- adding metdata to timeseries from ddlpy in `hatyan.ddlpy_to_hatyan()` in [#360](https://github.com/Deltares/hatyan/pull/360)

### Fix
- aligned (conversion of) units in [#364](https://github.com/Deltares/hatyan/pull/364)


## 2.9.0 (2024-09-11)

### Feat
- flexible timezone support in astrog in [#333](https://github.com/Deltares/hatyan/pull/333)
- support for missing `comp` timezone attribute in `hatyan.prediction()` in [#335](https://github.com/Deltares/hatyan/pull/335)
- read status from diafile in `hatyan.read_dia()` and aligned status/qualitycode columns in `hatyan.ddlpy_to_hatyan()` in [#337](https://github.com/Deltares/hatyan/pull/337) and [#353](https://github.com/Deltares/hatyan/pull/353)
- astrog datetimes as index instead of column in [#347](https://github.com/Deltares/hatyan/pull/347)

### Fix
- repaired support for equidistant multiblock diafiles with varying timesteps in [#314](https://github.com/Deltares/hatyan/pull/314)
- fixed passing of timezones in case of mixed timezones in comp and times in `hatyan.prediction()` in [#336](https://github.com/Deltares/hatyan/pull/336) and [#349](https://github.com/Deltares/hatyan/pull/349)
- cleaner metadata header by `hatyan.write_dia()` and conserving index frequency in [#354](https://github.com/Deltares/hatyan/pull/354)


## 2.8.0 (2024-05-08)
This release contains many changes of which several also impact the user. These changes benefit the user friendliness and robustness of hatyan.

### Feat
- reading/writing analysis settings (xfac, nodalfactors, fu_alltimes) from/to component file in [#175](https://github.com/Deltares/hatyan/pull/175)
- made hatyan callable via `hatyan /path/to/script.py` (deprecates `init_RWS()` and `exit_RWS()`) in [#113](https://github.com/Deltares/hatyan/pull/113) and [#226](https://github.com/Deltares/hatyan/issues/226)
- retained `freq` attribute of `ts.index` in case of multifile equidistant dia in [#118](https://github.com/Deltares/hatyan/pull/118)
- added support for file patterns in `hatyan.readts_dia()` in [#120](https://github.com/Deltares/hatyan/pull/120)
- uniform usage of `analysis` (deprecates `hatyan.get_components_from_ts()`) in [#125](https://github.com/Deltares/hatyan/pull/125)
- added metadata to timeseries and components and passing it between those objects in the hatyan process in [#131](https://github.com/Deltares/hatyan/pull/131)
- improved metadata in component files in [#131](https://github.com/Deltares/hatyan/pull/131)
- replaced `times_ext` `timestep_min` and `times_pred_all` for `hatyan.prediction()` with `times` argument in [#143](https://github.com/Deltares/hatyan/pull/143)
- integrated ddlpy to in ddl example script in [#202](https://github.com/Deltares/hatyan/pull/202) and [#213](https://github.com/Deltares/hatyan/pull/213)
- simplified ddlpy to hatyan timeseries dataframe conversoin in `hatyan.ddlpy_to_hatyan()` in https://github.com/Deltares/hatyan/pull/284
- prevent writing components with non-standard settings to file and prevent prediction with different settings than provided components in [#218](https://github.com/Deltares/hatyan/pull/218)
- use `pandas.DataFrame.attrs` instead of separate metadata attrs in [#219](https://github.com/Deltares/hatyan/pull/219)
- moved from prints to logging in [#236](https://github.com/Deltares/hatyan/pull/236)
- made `hatyan.get_diablocks()` private and added public `hatyan.get_diaxycoords()` instead [#244](https://github.com/Deltares/hatyan/pull/244)
- made `hatyan.Timeseries_Statistics()` public [#244](https://github.com/Deltares/hatyan/pull/244)
- support for timezones in hatyan timeseries in [#258](https://github.com/Deltares/hatyan/pull/258)
- simplified `hatyan.write_tsnetcdf()` by reading metadata from timeseries dataframes in [#260](https://github.com/Deltares/hatyan/pull/260)
- made all read/write timeserie function names consistent by dropping "ts" in [#265](https://github.com/Deltares/hatyan/pull/265)
- use prediction settings from `comp` argument instead of settings arguments in `hatyan.prediction()` [#263](https://github.com/Deltares/hatyan/issues/263)
- re-introduce clear arguments for `hatyan.analysis()` in [#270](https://github.com/Deltares/hatyan/pull/270)
- dropped support for `ìnt` timesteps in slices and in `hatyan.prediction()` [#288](https://github.com/Deltares/hatyan/pull/288)
- renamed `CS_comps` and `xTx_condition_max` arguments for `hatyan.prediction()` to lowercase [#299](https://github.com/Deltares/hatyan/pull/292)
- fix for outofbounds nanosecond timesteps in `hatyan.analysis()` in [#296](https://github.com/Deltares/hatyan/pull/296)

### Deprecated
- hatyan ddl functions `get_DDL_catalog()`, `get_DDL_queryserver()`, `get_DDL_data()`, `get_DDL_stationmetasubset()` in [#206](https://github.com/Deltares/hatyan/pull/206)
- moved kenmerkendewaarden functions and examples to https://github.com/Deltares-research/kenmerkendewaarden in [#211](https://github.com/Deltares/hatyan/pull/211)
- deprecated `hatyan.write_tsdia_HWLW()` in favor of `hatyan.write_tsdia()` in [#223](https://github.com/Deltares/hatyan/pull/223)
- deprecated `hatyan.convert_coords()` and `hatyan.check_ts()` in [#244](https://github.com/Deltares/hatyan/pull/244)
- deprecated `return_listoptions` argument for `hatyan.get_const_list_hatyan()` in [#262](https://github.com/Deltares/hatyan/pull/262)
- deprecated public `hatyan.HatyanSettings()` in [#270](https://github.com/Deltares/hatyan/pull/270)
- deprecated `hatyan.components_timeshift()` in [#281](https://github.com/Deltares/hatyan/pull/281)


## 2.7.0 (2023-08-03)

### Feat
- added py38/py310/py311 testbanks in [#89](https://github.com/Deltares/hatyan/pull/89)

### Fix
- introduced `min_width` of two hours to avoid incorrectly marked primary extremes in [#86](https://github.com/Deltares/hatyan/pull/86)
- fixed `pandas>=2.0.0` AttributeError for outofbounds datetimes in [commit aab36ba](https://github.com/Deltares/hatyan/commit/aab36ba6a5adeb4cec255f39c505f397f6a60be5)


## 2.6.0 (2023-02-15)

### Feat
- added Kenmerkende Waarden functions
- moved from analysis_peryear and analysis_permonth to analysis_perperiod
- modernized python package and updated dependencies in several commits
- updated `xTxmat_condition_max=12` (was `10` before) to allow several correct cases
- moved from analysis_peryear and analysis_permonth to analysis_perperiod
- added xtrack and fes2014b constituent list
- added shallow water components 2MK2, M(SK)2, M(KS)2. Added harmonic components S3, M!M, MSTM from tugo (provided by Henrique Guarneri)
- introduced extra component OQ2_tugo
- corrected nodal factor for MSQM based on tugo/HG
- added get_status keyword to `readts_dia()`
- used proper (non)equidistant distinguising method for diafiles
- added nyquist folding method to properly distinguish overlapping frequencies after folding
- introducing hatyan_settings
- improved readability of `split_compoonents()`
- introduced caching of v0uf data and discontinued pickle file, introduced harmonic/shallow csv file
- updated astrog leap-seconds.list for astrog and disable exact dT extrapolation

### Fix
- improved HWLWnumbering error in case of missing station from phasediff textfile
- improved rayleigh calculation method
- astrog bugfix for longitudes away from 0, removed for loops from astrab (faster code)


## 2.5.0 (2021-12-14)

### Feat
- added DDL to hatyan
- added `prediction_peryear()` definition
- support `analysis_peryear=True` and `CS_comps=!None` at the same time

### Fix
- bugfix in diawriting for last entries of timeseries
- bugfix in diafile writing that solves rounding issue. Made dia writing more efficient by writing values with pd.to_csv instead of in loop


## 2.3.0 (2021-10-08)

See removed [HISTORY.rst](https://github.com/Deltares/hatyan/blob/442f1b7b0975e40f29afe6fbf7d252c6271b3741/HISTORY.rst) for older features/fixes
