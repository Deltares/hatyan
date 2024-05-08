# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 17:03:52 2021

@author: veenstra
"""

import ddlpy # available via pip install rws-ddlpy or at https://github.com/Deltares/ddlpy
import matplotlib.pyplot as plt
plt.close("all")

# input parameters
# start_date = "2022-12-19 00:00:00 +01:00"
# end_date = "2022-12-31 00:00:00 +01:00"
start_date = "2020-03-22 00:00:00 +01:00" # DONAR validation data available for hoekvhld
end_date = "2020-03-29 00:00:00 +01:00"
#start_date = "2020-11-25 09:47:00 +01:00" #quite recent period
#end_date = "2021-01-30 09:50:00 +01:00"
#start_date = "1993-08-25 09:47:00 +01:00" #VLISSGN got new Waardebepalingsmethode in this year
#end_date = "1994-11-30 09:50:00 +01:00"
#start_date = "2009-01-01 00:00:00 +01:00" #common RWS retrieval period
#end_date = "2012-12-31 23:50:00 +01:00"


locations = ddlpy.locations()


######### online waterlevel data retrieval for one station
if 1: #for RWS
    import hatyan # available via `pip install hatyan` or at https://github.com/Deltares/hatyan
    include_extremes = True
    
    bool_grootheid_meas = locations['Grootheid.Code'].isin(['WATHTE'])
    bool_grootheid_astro = locations['Grootheid.Code'].isin(['WATHTBRKD'])
    bool_hoedanigheid = locations['Hoedanigheid.Code'].isin(['NAP'])
    bool_groepering_ts = locations['Groepering.Code'].isin(['NVT'])
    bool_groepering_ext = locations['Groepering.Code'].isin(['GETETM2','GETETBRKD2','GETETMSL2','GETETBRKDMSL2'])
    
    # get WATHTE locations (ts) >> measured waterlevel
    locs_wathte = locations.loc[bool_grootheid_meas & bool_hoedanigheid & bool_groepering_ts]
    # get WATHTBRKD locations (ts) >> computed astronomical waterlevel
    locs_wathtbrkd = locations.loc[bool_grootheid_astro & bool_hoedanigheid & bool_groepering_ts]
    
    if include_extremes:
        locs_wathte_ext = locations.loc[bool_grootheid_meas & bool_hoedanigheid & bool_groepering_ext]
        locs_wathtbrkd_ext = locations.loc[bool_grootheid_astro & bool_hoedanigheid & bool_groepering_ext]
        
        # get types locations (ts and extremes)
        # we cannot subset locations on Typering in ddlpy (NVT/GETETTPE), so we use a groepering+grootheid combination
        bool_grootheid_exttypes = locations['Grootheid.Code'].isin(['NVT'])
        bool_groepering_ext_meas = locations['Groepering.Code'].isin(['GETETM2','GETETMSL2'])
        bool_groepering_ext_astro = locations['Groepering.Code'].isin(['GETETBRKD2','GETETBRKDMSL2'])
        locs_exttypes_wathte = locations.loc[bool_grootheid_exttypes & bool_groepering_ext_meas]
        locs_exttypes_wathtbrkd = locations.loc[bool_grootheid_exttypes & bool_groepering_ext_astro]
    
    for current_station in ['AWGPFM']:
        locs_wathte_one = locs_wathte.loc[locs_wathte.index.isin([current_station])]
        locs_wathtbrkd_one = locs_wathtbrkd.loc[locs_wathtbrkd.index.isin([current_station])]
        
        # no support for multiple rows, so pass one at a time
        if len(locs_wathte_one) > 1:
            raise Exception("multiple stations in locs_wathte dataframe")
        if len(locs_wathtbrkd_one) > 1:
            raise Exception("multiple stations in locs_wathtbrkd dataframe")
        meas_wathte = ddlpy.measurements(locs_wathte_one.iloc[0], start_date=start_date, end_date=end_date)
        meas_wathtbrkd = ddlpy.measurements(locs_wathtbrkd_one.iloc[0], start_date=start_date, end_date=end_date)
        # hatyan timeseries
        ts_measwl = hatyan.ddlpy_to_hatyan(meas_wathte)
        ts_astro = hatyan.ddlpy_to_hatyan(meas_wathtbrkd)
        
        if include_extremes:
            locs_wathte_ext_one = locs_wathte_ext.loc[locs_wathte_ext.index.isin([current_station])]
            locs_wathtbrkd_ext_one = locs_wathtbrkd_ext.loc[locs_wathtbrkd_ext.index.isin([current_station])]
            locs_wathte_exttypes_one = locs_exttypes_wathte.loc[locs_exttypes_wathte.index.isin([current_station])]
            locs_wathtbrkd_exttypes_one = locs_exttypes_wathtbrkd.loc[locs_exttypes_wathtbrkd.index.isin([current_station])]
            # no support for multiple rows, so pass one at a time
            if len(locs_wathte_ext_one) > 1:
                raise Exception("multiple stations in locs_wathte_ext dataframe")
            if len(locs_wathtbrkd_ext_one) > 1:
                raise Exception("multiple stations in locs_wathtbrkd_ext dataframe")
            if len(locs_wathte_exttypes_one) > 1:
                raise Exception("multiple stations in locs_wathte_exttypes dataframe")
            if len(locs_wathtbrkd_exttypes_one) > 1:
                raise Exception("multiple stations in locs_wathtbrkd_exttypes dataframe")
            meas_wathte_ext = ddlpy.measurements(locs_wathte_ext_one.iloc[0], start_date=start_date, end_date=end_date)
            meas_wathtbrkd_ext = ddlpy.measurements(locs_wathtbrkd_ext_one.iloc[0], start_date=start_date, end_date=end_date)
            meas_wathte_exttypes = ddlpy.measurements(locs_wathte_exttypes_one.iloc[0], start_date=start_date, end_date=end_date)
            meas_wathtbrkd_exttypes = ddlpy.measurements(locs_wathtbrkd_exttypes_one.iloc[0], start_date=start_date, end_date=end_date)
            # hatyan timeseries
            ts_measwlHWLW = hatyan.ddlpy_to_hatyan(meas_wathte_ext, meas_wathte_exttypes)
            ts_astroHWLW = hatyan.ddlpy_to_hatyan(meas_wathtbrkd_ext, meas_wathtbrkd_exttypes)
    
    stat_name = locs_wathte_one["Naam"].iloc[0]
    stat_code = locs_wathte_one.index[0]
    if include_extremes:
        fig,(ax1,ax2) = hatyan.plot_timeseries(ts=ts_astro,ts_validation=ts_measwl,ts_ext=ts_astroHWLW,ts_ext_validation=ts_measwlHWLW)
        ax1.set_title(f'waterlevels for {stat_code} ({stat_name})')
        ax2.set_ylim(-0.5,0.5)
    else:
        fig,(ax1,ax2) = hatyan.plot_timeseries(ts=ts_astro,ts_validation=ts_measwl)
        ax1.set_title(f'waterlevels for {stat_code} ({stat_name})')
        ax2.set_ylim(-0.5,0.5)


######### simple waterlevel data retrieval for all waterlevel stations or all stations
if 1: #for CMEMS
    
    bool_grootheid = locations['Grootheid.Code'].isin(['WATHTE']) # measured waterlevels (not astro)
    bool_groepering = locations['Groepering.Code'].isin(['NVT']) # timeseries (not extremes)
    bool_hoedanigheid = locations['Hoedanigheid.Code'].isin(['NAP']) # vertical reference
    locs_wathte = locations.loc[bool_grootheid & bool_groepering & bool_hoedanigheid]
    
    for current_station in ['HOEKVHLD']:
        locs_wathte_one = locs_wathte.loc[locs_wathte.index.isin([current_station])]
        
        # no support for multiple rows, so pass one at a time
        if len(locs_wathte_one) > 1:
            raise Exception("duplicate stations for wathte")
        
        meas_wathte = ddlpy.measurements(locs_wathte_one.iloc[0], start_date=start_date, end_date=end_date)
        
        stat_name = locs_wathte_one.iloc[0]["Naam"]
        stat_code = current_station
        fig, (ax1,ax2) = plt.subplots(2,1, figsize=(8,6), sharex=True)
        ax1.plot(meas_wathte["Meetwaarde.Waarde_Numeriek"])
        # QC codes are documented in https://rijkswaterstaatdata.nl/publish/pages/223004/uitleg-attributen-bij-bulkwaarnemingenservice-dd-06-11-2023-.pdf
        ax2.plot(meas_wathte["WaarnemingMetadata.KwaliteitswaardecodeLijst"].astype(int))
        ax1.set_title(f'waterlevels for {stat_code} ({stat_name})')
        ax2.set_title(f'QC for {stat_code} ({stat_name})')
        fig.tight_layout()
