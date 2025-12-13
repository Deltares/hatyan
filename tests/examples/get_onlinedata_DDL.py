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
start_date = "2020-03-22 00:00:00 +01:00" # DONAR validation data available for hoekvanholland
end_date = "2020-03-29 00:00:00 +01:00"
#start_date = "2020-11-25 09:47:00 +01:00" #quite recent period
#end_date = "2021-01-30 09:50:00 +01:00"
#start_date = "1993-08-25 09:47:00 +01:00" #vlissingen got new Waardebepalingsmethode in this year
#end_date = "1994-11-30 09:50:00 +01:00"
#start_date = "2009-01-01 00:00:00 +01:00" #common RWS retrieval period
#end_date = "2012-12-31 23:50:00 +01:00"


locations = ddlpy.locations()

# retrieval for one station including comparison to DONAR validation dataset
if 1: # for RWS
    import os
    import hatyan
    import urllib
    
    # use time range available in validation data, passing winter/summertime date
    current_station = "hoekvanholland"
    
    # first download validation data
    url_hatyan = "https://raw.githubusercontent.com/Deltares/hatyan/main/tests/data_unitsystemtests/"
    file_vali = 'diawia_HOEKVHLD_astro_tijdreeks.dia'
    url = url_hatyan+file_vali
    with urllib.request.urlopen(url) as response:
        data = response.read()
    with open(file_vali, "w") as f:
        f.write(data.decode('utf-8'))
    ts_vali = hatyan.read_dia(file_vali)
    ts_vali = hatyan.crop_timeseries(ts_vali, times=slice(start_date, end_date))
    
    bool_procestype = locations['ProcesType'].isin(['astronomisch']) # astro (not measured)
    bool_grootheid = locations['Grootheid.Code'].isin(['WATHTE']) # waterlevels
    bool_groepering = locations['Groepering.Code'].isin(['']) # timeseries (not extremes)
    bool_hoedanigheid = locations['Hoedanigheid.Code'].isin(['NAP']) # vertical reference
    locs_wathte = locations.loc[bool_procestype & bool_grootheid & bool_groepering & bool_hoedanigheid]
    
    locs_wathte_one = locs_wathte.loc[locs_wathte.index.isin([current_station])]
    
    # no support for multiple rows, so pass one at a time
    if len(locs_wathte_one) > 1:
        raise Exception("duplicate stations for WATHTE")

    # get the measurements
    meas_wathte = ddlpy.measurements(locs_wathte_one.iloc[0], start_date=start_date, end_date=end_date)
    meas_wathte_hat = hatyan.ddlpy_to_hatyan(meas_wathte)
    
    stat_name = locs_wathte_one.iloc[0]["Naam"]
    stat_code = current_station
    fig, ax = plt.subplots(figsize=(8,5))
    ax.plot(ts_vali["values"], label="donar timeseries")
    ax.plot(meas_wathte_hat["values"], "--", label="ddlpy timeseries")
    ax.set_title(f'waterlevels for {stat_code} ({stat_name})')
    ax.legend(loc=1)
    fig.tight_layout()
    os.remove(file_vali)


######### online waterlevel data retrieval for one station
if 1: #for RWS
    import hatyan # available via `pip install hatyan` or at https://github.com/Deltares/hatyan
    include_extremes = True
    
    bool_procestype_meas = locations['ProcesType'].isin(['meting'])
    bool_procestype_astro = locations['ProcesType'].isin(['astronomisch'])
    bool_grootheid = locations['Grootheid.Code'].isin(['WATHTE'])
    bool_hoedanigheid = locations['Hoedanigheid.Code'].isin(['NAP'])
    bool_groepering_ts = locations['Groepering.Code'].isin([''])
    bool_groepering_ext = locations['Groepering.Code'].isin(['GETETM2','GETETBRKD2','GETETMSL2','GETETBRKDMSL2'])
    
    # get meas locations (ts) >> measured waterlevel
    locs_meas = locations.loc[bool_procestype_meas & bool_grootheid & bool_hoedanigheid & bool_groepering_ts]
    # get astro locations (ts) >> computed astronomical waterlevel
    locs_astro = locations.loc[bool_procestype_astro & bool_grootheid & bool_hoedanigheid & bool_groepering_ts]
    
    if include_extremes:
        locs_meas_ext = locations.loc[bool_procestype_meas & bool_grootheid & bool_hoedanigheid & bool_groepering_ext]
        locs_astro_ext = locations.loc[bool_procestype_astro & bool_grootheid & bool_hoedanigheid & bool_groepering_ext]
        
        # get types locations (ts and extremes)
        # we cannot subset locations on Typering in ddlpy (NVT/GETETTPE), so we use a groepering+grootheid combination
        bool_grootheid_exttypes = locations['Grootheid.Code'].isin(['NVT'])
        bool_groepering_ext_meas = locations['Groepering.Code'].isin(['GETETM2','GETETMSL2'])
        bool_groepering_ext_astro = locations['Groepering.Code'].isin(['GETETBRKD2','GETETBRKDMSL2'])
        locs_exttypes_meas = locations.loc[bool_grootheid_exttypes & bool_groepering_ext_meas]
        locs_exttypes_astro = locations.loc[bool_grootheid_exttypes & bool_groepering_ext_astro]
    
    for current_station in ['hoekvanholland']:
        locs_meas_one = locs_meas.loc[locs_meas.index.isin([current_station])]
        locs_astro_one = locs_astro.loc[locs_astro.index.isin([current_station])]
        
        # no support for multiple rows, so pass one at a time
        if len(locs_meas_one) > 1:
            raise Exception("multiple stations in locs_meas dataframe")
        if len(locs_astro_one) > 1:
            raise Exception("multiple stations in locs_astro dataframe")
        meas_meas = ddlpy.measurements(locs_meas_one.iloc[0], start_date=start_date, end_date=end_date)
        meas_astro = ddlpy.measurements(locs_astro_one.iloc[0], start_date=start_date, end_date=end_date)
        # hatyan timeseries
        ts_measwl = hatyan.ddlpy_to_hatyan(meas_meas)
        ts_astro = hatyan.ddlpy_to_hatyan(meas_astro)
        
        if include_extremes:
            locs_meas_ext_one = locs_meas_ext.loc[locs_meas_ext.index.isin([current_station])]
            locs_astro_ext_one = locs_astro_ext.loc[locs_astro_ext.index.isin([current_station])]
            locs_meas_exttypes_one = locs_exttypes_meas.loc[locs_exttypes_meas.index.isin([current_station])]
            locs_astro_exttypes_one = locs_exttypes_astro.loc[locs_exttypes_astro.index.isin([current_station])]
            # no support for multiple rows, so pass one at a time
            if len(locs_meas_ext_one) > 1:
                raise Exception("multiple stations in locs_meas_ext dataframe")
            if len(locs_astro_ext_one) > 1:
                raise Exception("multiple stations in locs_astro_ext dataframe")
            if len(locs_meas_exttypes_one) > 1:
                raise Exception("multiple stations in locs_meas_exttypes dataframe")
            if len(locs_astro_exttypes_one) > 1:
                raise Exception("multiple stations in locs_astro_exttypes dataframe")
            meas_meas_ext = ddlpy.measurements(locs_meas_ext_one.iloc[0], start_date=start_date, end_date=end_date)
            meas_astro_ext = ddlpy.measurements(locs_astro_ext_one.iloc[0], start_date=start_date, end_date=end_date)
            meas_meas_exttypes = ddlpy.measurements(locs_meas_exttypes_one.iloc[0], start_date=start_date, end_date=end_date)
            meas_astro_exttypes = ddlpy.measurements(locs_astro_exttypes_one.iloc[0], start_date=start_date, end_date=end_date)
            # hatyan timeseries
            ts_measwlHWLW = hatyan.ddlpy_to_hatyan(meas_meas_ext, meas_meas_exttypes)
            ts_astroHWLW = hatyan.ddlpy_to_hatyan(meas_astro_ext, meas_astro_exttypes)
    
    stat_name = locs_meas_one["Naam"].iloc[0]
    stat_code = locs_meas_one.index[0]
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
    
    bool_procestype = locations['ProcesType'].isin(['meting']) # measured waterlevels (not astro)
    bool_grootheid = locations['Grootheid.Code'].isin(['WATHTE']) # waterlevels
    bool_groepering = locations['Groepering.Code'].isin(['']) # timeseries (not extremes)
    bool_hoedanigheid = locations['Hoedanigheid.Code'].isin(['NAP']) # vertical reference
    locs_wathte = locations.loc[bool_procestype & bool_grootheid & bool_groepering & bool_hoedanigheid]
    
    for current_station in ['hoekvanholland']:
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
        ax2.plot(meas_wathte["WaarnemingMetadata.Kwaliteitswaardecode"].astype(int))
        ax1.set_title(f'waterlevels for {stat_code} ({stat_name})')
        ax2.set_title(f'QC for {stat_code} ({stat_name})')
        fig.tight_layout()
