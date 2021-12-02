# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 17:03:03 2021

@author: veenstra
"""



def get_DDL_catalog(catalog_extrainfo=[]):
    """
    check get_DDL_data() for details

    Parameters
    ----------
    catalog_filter : TYPE, optional
        DESCRIPTION. The default is []. Possibilities are in https://rijkswaterstaat.github.io/wm-ws-dl/?json#ophalencatalogus, for instance 'MeetApparaten' and 'Parameters'.

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    result_cat_dict : TYPE
        DESCRIPTION.
    """

    import requests
    import pandas as pd
    
    #the webservices 
    url_catalog = 'https://waterwebservices.rijkswaterstaat.nl/METADATASERVICES_DBO/OphalenCatalogus'
    
    #The request for ophalencatalogus
    catalog_filter = ['Compartimenten','Eenheden','Grootheden','Hoedanigheden','Groeperingen']+catalog_extrainfo
    request_cat = {"CatalogusFilter": {x:True for x in catalog_filter}}
    
    # pull catalog from the API and store in json format
    resp = requests.post(url_catalog, json=request_cat)
    if not resp.ok:
        raise Exception('%s for %s: %s'%(resp.reason, resp.url, str(resp.text)))
    result_cat = resp.json()
    #print('json catalog keys: %s'%(result_cat.keys())) #['AquoMetadataLijst', 'LocatieLijst', 'AquoMetadataLocatieLijst', 'BemonsteringshoogteLijst', 'OpdrachtgevendeInstantieLijst', 'KwaliteitswaardecodeLijst', 'ReferentievlakLijst', 'StatuswaardeLijst', 'Succesvol'] or ['Succesvol', 'Foutmelding']
    if not result_cat['Succesvol']:
        raise Exception('catalog query not succesful, Foutmelding: %s'%(result_cat['Foutmelding']))
    
    result_cat_dict = {}
    for catalog_key in result_cat.keys():
        if catalog_key!='Succesvol':
            result_cat_dict[catalog_key] = pd.json_normalize(result_cat[catalog_key])
    return result_cat_dict


def get_DDL_queryserver(query_station,query_metadata,query_tstart,query_tstop,check_available=False):
    """
    check get_DDL_data() for details
    """
    import requests
        
    tzinfo_numraw = query_tstart.strftime('%z') #'+0100'
    tzinfo_numstr = tzinfo_numraw[:3]+':'+tzinfo_numraw[-2:] #'+01:00'
    query_tstart_str = query_tstart.strftime('%Y-%m-%dT%H:%M:%S.000'+tzinfo_numstr) #"2021-01-14T09:47:00.000+01:00"
    query_tstop_str  = query_tstop.strftime('%Y-%m-%dT%H:%M:%S.000'+tzinfo_numstr) #"2021-11-27T10:00:00.000+01:00"
    
    if check_available: # Check if data is available
        url_ddl = 'https://waterwebservices.rijkswaterstaat.nl/ONLINEWAARNEMINGENSERVICES_DBO/CheckWaarnemingenAanwezig'
        request_ddl = {"AquoMetadataLijst" :[query_metadata],
                         "LocatieLijst":[query_station],
                         "Periode":{"Begindatumtijd":query_tstart_str, # TODO IMPROVEMENT: if user accidentally switches start/stop dates in query, OphalenWaarnemingen returns 'Begindatum is groter dan einddatum. (check_available=False)', CheckWaarnemingenAanwezig crashes instead of returning a proper error
                                    "Einddatumtijd":query_tstop_str}
                        }
        # TODO IMPROVEMENT: would be valuable to quickly get available start/stop time and number of available measurements (timesteps). Below is a (slow) example, seems to take as much time as retrieving measurements
        # TODO IMPROVEMENT: welke groeperingsperiodes zijn beschikbaar en kan 'geen' ook? (JsonProcessingException: Can not construct instance of nl.ordina.request.OphalenAantalWaarnemingenRequest$Groepering from String value 'geen': value not one of declared Enum instance names.)
        """
        url_ddl = 'https://waterwebservices.rijkswaterstaat.nl/ONLINEWAARNEMINGENSERVICES_DBO/OphalenAantalWaarnemingen'
        request_ddl = {"AquoMetadataLijst" :[query_metadata],
                         "Groeperingsperiode" : 'Jaar', 
                         "LocatieLijst":[query_station],
                         "Periode":{"Begindatumtijd":query_tstart_str,
                                    "Einddatumtijd":query_tstop_str}
                        }
        print(result['AantalWaarnemingenPerPeriodeLijst'][0]['AantalMetingenPerPeriodeLijst'])
        """
        # TODO IMPROVEMENT: OphalenLaatsteWaarnemingen seems to be valuable to get the end time for a station, however resulted waarnemingenlijst for one station has multiple entries (probably multiple WaardeBepalingsmethode?) but also MetingenLijst sometimes also has multiple entries, how to interpret this?
        """
        url_ddl = 'https://waterwebservices.rijkswaterstaat.nl/ONLINEWAARNEMINGENSERVICES_DBO/OphalenLaatsteWaarnemingen'
        request_ddl = {"AquoPlusWaarnemingMetadataLijst":[{"AquoMetadata":query_metadata}],"LocatieLijst":[query_station]}
        result['WaarnemingenLijst'][2]['MetingenLijst'] 
        """
    else:
        #retrieve data
        url_ddl = 'https://waterwebservices.rijkswaterstaat.nl/ONLINEWAARNEMINGENSERVICES_DBO/OphalenWaarnemingen'
        request_ddl = {"AquoPlusWaarnemingMetadata":{"AquoMetadata":query_metadata},
                      "Locatie":query_station,#{"X":518882.333320247,"Y":5760829.11729589,"Code":"EURPFM"}, # TODO IMPROVEMENT: it seems not not possible to retreive by station Naam/Code/Locatie_MessageID only. X+Y+Code is minimum, so supplying entire dict. Why is this so strict? It seems odd that one needs to supply a six decimal RD coordinate (so micrometer accuracy) while the Code is in itself already unique. Also supplying Locatie_MessageID only is not possible
                      "Periode":{"Begindatumtijd":query_tstart_str, # TODO IMPROVEMENT: longer timeseries (eg 4 years) take a long time or return error (Foutmelding: Het max aantal waarnemingen (157824) is overschreven, beperk uw request.). Can this not be extended? >> now retrieving per year, is that always possible with this limit?
                                 "Einddatumtijd":query_tstop_str} # TODO IMPROVEMENT: recent data for eg HOEKVLD is not available but it is as station HOEK, can these stations not be one, but with a different kaliteitscode/statuswaarde/GrootheidCode/GroeperingCode etc? I was a bit surprised that 'ongecontroleerd' (and HOEK in general) also has kwaliteitscode=0
                      }
    #print(request_ddl)
    resp = requests.post(url_ddl, json=request_ddl)
    if not resp.ok:
        raise Exception('%s for %s: %s'%(resp.reason, resp.url, str(resp.text)))
    result = resp.json()
    if not result['Succesvol']:
        raise Exception('query not succesful, Foutmelding: %s from %s'%(result['Foutmelding'],url_ddl))
    return result


def get_DDL_data(query_station,query_tstart,query_tstop,query_tzone='UTC+01:00',meta_dict={'Grootheid.Code':'WATHTE','Groepering.Code':'NVT'},allow_multipleresultsfor=None):
    """
    ddl tutorial: https://rijkswaterstaat.github.io/wm-ws-dl/?python#tutorial-locations
    normalizing json output: https://towardsdatascience.com/how-to-convert-json-into-a-pandas-dataframe-100b2ae1e0d8

    query_tzone: MET/CET results in Europe/Amsterdam (so including DST), use fixed offset instead
    """
    import pandas as pd
    import numpy as np
    import datetime as dt
    import pytz

    #parse meta_dict to query_metadata dict
    query_metadata = {}
    for metakeypoint in meta_dict:
        metakeymain = metakeypoint.split('.')[0]
        metakeysub = metakeypoint.split('.')[1]
        query_metadata[metakeymain] = {metakeysub:meta_dict[metakeypoint]}

    if query_tzone.startswith('UTC+') or query_tzone.startswith('UTC-'): #parse to fixed offset like 'Etc/GMT-1'. +/- are counter intuitive but it works: https://pvlib-python.readthedocs.io/en/stable/timetimezones.html#fixedoffsets)
        if len(query_tzone)!=9 or not query_tzone.endswith(':00'):
            raise Exception('if query_tzone starts with UTC+ or UTC-, the string should be 9 characters long and have 0 minutes, like "UTC+01:00"')
        query_tzone_hr = int(query_tzone[4:6])
        if query_tzone[3]=='+':
            query_tzone = 'Etc/GMT-%d'%(query_tzone_hr)
        else:
            query_tzone = 'Etc/GMT+%d'%(query_tzone_hr)
    tzinfo = pytz.timezone(query_tzone)
    query_tstart = query_tstart.replace(tzinfo=tzinfo)
    query_tstop = query_tstop.replace(tzinfo=tzinfo)
    year_list = range(query_tstart.year, query_tstop.year+1)
    
    print('processing station %s/%s: %s to %s (%d years)'%(query_station['Code'],query_station['Naam'],query_tstart,query_tstop,len(year_list)))
    result_available = get_DDL_queryserver(query_station,query_metadata,query_tstart,query_tstop,check_available=True)
    if result_available['WaarnemingenAanwezig']!='true': # TODO IMPROVEMENT: WaarnemingenAanwezig is now a 'true' string instead of a True boolean (result_available['Succesvol'] is also a boolean)
        print('WARNING: no values present for this query, returning None')
    else:
        result_wl0_metingenlijst_alldates = pd.DataFrame()
        result_wl0_aquometadata_unique = pd.DataFrame()
        for year in year_list:
            query_tstart_oneyear = np.maximum(query_tstart,dt.datetime(year,1,1,tzinfo=tzinfo))
            query_tstop_oneyear = np.minimum(query_tstop,dt.datetime(year+1,1,1,tzinfo=tzinfo))
    
            result_available = get_DDL_queryserver(query_station,query_metadata,query_tstart_oneyear,query_tstop_oneyear,check_available=True)
            if result_available['WaarnemingenAanwezig']!='true': # TODO IMPROVEMENT: WaarnemingenAanwezig is now a 'true' string instead of a True boolean (result_available['Succesvol'] is also a boolean)
                print('year %d: no values'%(year))
                continue
            print('year %d: retrieving data'%(year))
            
            result_wl = get_DDL_queryserver(query_station,query_metadata,query_tstart_oneyear,query_tstop_oneyear,check_available=False)
            
            #print('json result keys: %s'%(result_wl.keys())) #['WaarnemingenLijst', 'Succesvol'] or ['Succesvol', 'Foutmelding']
            if not result_wl['Succesvol']:
                raise Exception('measurement query not succesful, Foutmelding: %s'%(result_wl['Foutmelding']))
            
            range_nresults = range(len(result_wl['WaarnemingenLijst']))
            for result_idx in range_nresults:
                result_wl0 = result_wl['WaarnemingenLijst'][result_idx]
                #print('json result waarnemingenlijst keys: %s'%(result_wl0.keys())) #['Locatie', 'MetingenLijst', 'AquoMetadata']
                result_wl0_metingenlijst = pd.json_normalize(result_wl0['MetingenLijst']) # the actual waterlevel data for this station
                if not result_wl0_metingenlijst['Tijdstip'].is_monotonic_increasing:
                    #print('WARNING: retrieved timeseries is not monotonic increasing, so it was sorted')
                    result_wl0_metingenlijst = result_wl0_metingenlijst.sort_values('Tijdstip').reset_index(drop=True) # TODO IMPROVEMENT: data in response is not always sorted on time
                last_timestamp_tzaware = pd.to_datetime(result_wl0_metingenlijst['Tijdstip'].iloc[-1]).tz_convert(query_tstart.tzinfo)
                if not last_timestamp_tzaware.isoformat().startswith(str(year)): #need to remove the last data entry if it is 1 january in next year (correct for timezone first and this does not happen for eg extremes since they probably do not have a value on that exact datetime)
                    result_wl0_metingenlijst = result_wl0_metingenlijst.iloc[:-1]
                result_wl0_metingenlijst_alldates = result_wl0_metingenlijst_alldates.append(result_wl0_metingenlijst)
                #result_wl0_locatie = pd.json_normalize(result_wl0['Locatie'])
                result_wl0_aquometadata = pd.json_normalize(result_wl0['AquoMetadata'])
                result_wl0_aquometadata_unique = result_wl0_aquometadata_unique.append(result_wl0_aquometadata).drop_duplicates()
            
            if allow_multipleresultsfor is None:
                result_wl0_aquometadata_uniqueallowed = result_wl0_aquometadata_unique
            else: #allow multiple results for eg ['MeetApparaat', 'WaardeBepalingsmethode']
                try:
                    if not isinstance(allow_multipleresultsfor,list):
                        allow_multipleresultsfor = [allow_multipleresultsfor]
                    list_dropcolumns = ['AquoMetadata_MessageID']+['%s.%s'%(colname,postfix) for colname in allow_multipleresultsfor for postfix in ['Code','Omschrijving']]
                    result_wl0_aquometadata_uniqueallowed = result_wl0_aquometadata_unique.drop(list_dropcolumns,axis=1).drop_duplicates()
                except KeyError as e:
                    metakeys_forCode = [x.replace('.Code','') for x in result_wl0_aquometadata_unique.keys() if x.endswith('.Code')]
                    raise Exception('%s, available are: %s'%(e,metakeys_forCode))
            if len(result_wl0_aquometadata_uniqueallowed)>1:
                bool_nonuniquecols = (result_wl0_aquometadata_unique.iloc[0]!=result_wl0_aquometadata_unique).any(axis=0)
                metakeys_forCode_nonunique = [x.replace('.Code','') for x in result_wl0_aquometadata_unique.loc[:,bool_nonuniquecols].columns if x.endswith('.Code')]
                for iR, result_one in enumerate(result_wl['WaarnemingenLijst']):
                    print(f'Result {iR+1}:')
                    metakey_list = ['Compartiment','Eenheid','Grootheid','Hoedanigheid','Groepering']+metakeys_forCode_nonunique
                    for metakey in metakey_list:
                        print('%28s: %s,'%("'%s'"%metakey,result_one['AquoMetadata'][metakey]))
                raise Exception('query returned more than one result (differences in %s, details above), use more specific query_metadata argument or more extensive allow_multipleresultsfor argument (the latter might result in duplicate timesteps)'%(metakeys_forCode_nonunique))
                
                
        result_wl0_metingenlijst_alldates = result_wl0_metingenlijst_alldates.reset_index(drop=True)
        # TODO IMPROVEMENT: WaarnemingMetadata: all values are nested lists of length 1, can be flattened (they are actually not list/lijst, but statuswaarde instead of statuswaardelijst and kwaliteitswaardecode instead of kwaliteitswaardecodelijst).
        # TODO IMPROVEMENT: WaarnemingMetadata: Bemonsteringshoogte/Referentievlak/OpdrachtgevendeInstantie is probably constant for each query result, so could be added to aquometadata frame instead of a value per timestep (probably makes query faster and can be longer)
        # TODO IMPROVEMENT: WaarnemingMetadata: there seems to be no explanation in the catalog or metadata of the KwaliteitswaardecodeLijst values
        # create improved pandas DataFrame
        ts_meas_pd = pd.DataFrame({'values':result_wl0_metingenlijst_alldates['Meetwaarde.Waarde_Numeriek'].values,
                                   'QC':result_wl0_metingenlijst_alldates['WaarnemingMetadata.KwaliteitswaardecodeLijst'].str[0].astype(int).values, 
                                   'Status':result_wl0_metingenlijst_alldates['WaarnemingMetadata.StatuswaardeLijst'].str[0].values,
                                   #'Bemonsteringshoogte':result_wl0_metingenlijst_alldates['WaarnemingMetadata.BemonsteringshoogteLijst'].str[0].astype(int).values, 
                                   #'Referentievlak':result_wl0_metingenlijst_alldates['WaarnemingMetadata.ReferentievlakLijst'].str[0].values,
                                   #'OpdrachtgevendeInstantie':result_wl0_metingenlijst_alldates['WaarnemingMetadata.OpdrachtgevendeInstantieLijst'].str[0].values,
                                   },
                                  index=pd.to_datetime(result_wl0_metingenlijst_alldates['Tijdstip']))
        #convert timezone from MET to requested timezone
        ts_meas_pd.index = ts_meas_pd.index.tz_convert(query_tstart.tzinfo)

        bool_timeduplicated = ts_meas_pd.index.duplicated()
        if bool_timeduplicated.any():
            print('WARNING: query returned %d duplicate times, use less extensive allow_multipleresultsfor'%bool_timeduplicated.sum())
 
        return ts_meas_pd, result_wl0_aquometadata_unique


def get_DDL_stationmetasubset(catalog_dict, station=None,stationcolumn='Naam',meta_dict=None, error_empty=True):
    import numpy as np
    
    cat_aquometadatalijst = catalog_dict['AquoMetadataLijst'].set_index('AquoMetadata_MessageID')
    cat_locatielijst = catalog_dict['LocatieLijst'].set_index('Locatie_MessageID')
    cat_AquoMetadataLocatieLijst_locidx = catalog_dict['AquoMetadataLocatieLijst'].set_index('Locatie_MessageID')
    cat_AquoMetadataLocatieLijst_metaidx = catalog_dict['AquoMetadataLocatieLijst'].set_index('AquoMetaData_MessageID')
    # TODO IMPROVEMENT: AquoMetadataLocatieLijst is missing AquoMetaData_MessageIDs: [38, 43, 75, 98, 99, 176] (and probably somewhat more), so these will not be present in cat_locatielijst_sel
    
    if meta_dict is None and station is None: #this makes the rest of the function a bit simpler to write
        raise Exception('using this function has no added value when not supplying meta and station')

    if meta_dict is not None:
        #loop over meta_dict keys, append to array and 
        bool_aquometadatalijst_containingmeta_list = []
        for metakey in meta_dict.keys():
            bool_aquometadatalijst_containingmeta_list.append(cat_aquometadatalijst[metakey].str.contains(meta_dict[metakey],case=False))
        bool_aquometadatalijst_containingmeta = np.array(bool_aquometadatalijst_containingmeta_list).all(axis=0) #logical and
        
        cat_aquometadatalijst_containingmeta = cat_aquometadatalijst.loc[bool_aquometadatalijst_containingmeta]
        bool_metaid_intranslatetable = cat_aquometadatalijst.index.isin(cat_AquoMetadataLocatieLijst_metaidx.index)
        if not bool_metaid_intranslatetable.all():
            cat_AquoMetadataLocatieLijst_missingAquoIDs = cat_aquometadatalijst.loc[~bool_metaid_intranslatetable].index
            print('WARNING: AquoMetadataLocatieLijst is missing AquoMetaData_MessageIDs:\n%s\nThese IDs will not be present in cat_locatielijst_sel. This issue can be reported to Servicedesk data via: https://www.rijkswaterstaat.nl/formulieren/contactformulier-servicedesk-data'%(cat_aquometadatalijst.loc[cat_AquoMetadataLocatieLijst_missingAquoIDs,['Grootheid.Code','Grootheid.Omschrijving']]))
            cat_aquometadatalijst_containingmeta = cat_aquometadatalijst.loc[bool_aquometadatalijst_containingmeta & bool_metaid_intranslatetable]
        bool_locatielijst_containingmeta = cat_locatielijst[stationcolumn].str.contains('THISSTRINGISNOTINCOLUMN') #first generate all false bool
        bool_locatielijst_containingmeta_idxtrue = cat_AquoMetadataLocatieLijst_metaidx.loc[cat_aquometadatalijst_containingmeta.index,'Locatie_MessageID']
        bool_locatielijst_containingmeta.loc[bool_locatielijst_containingmeta_idxtrue] = True
    if station is not None:
        bool_locatielijst_containingstation = cat_locatielijst[stationcolumn].str.contains(station,case=False)
        cat_locatielijst_containingstation = cat_locatielijst.loc[bool_locatielijst_containingstation]
        bool_locid_intranslatetable = cat_locatielijst.index.isin(cat_AquoMetadataLocatieLijst_locidx.index)
        if not bool_locid_intranslatetable.all():
            raise Exception('first time that this occurs, check code')
            cat_AquoMetadataLocatieLijst_missingLocIDs = cat_locatielijst.loc[~bool_locid_intranslatetable].index
            print('WARNING: AquoMetadataLocatieLijst is missing Locatie_MessageIDs:\n%s\nThese IDs will not be present in cat_aquometadatalijst_sel. This issue can be reported to Servicedesk data via: https://www.rijkswaterstaat.nl/formulieren/contactformulier-servicedesk-data'%(cat_locatielijst.loc[cat_AquoMetadataLocatieLijst_missingLocIDs]))
            cat_locatielijst_containingstation = cat_locatielijst.loc[bool_locatielijst_containingstation & bool_locid_intranslatetable]
        bool_aquometadatalijst_containingstation = cat_aquometadatalijst['Grootheid.Omschrijving'].str.contains('THISSTRINGISNOTINCOLUMN') #first generate all false bool
        bool_aquometadatalijst_containingstation_idxtrue = cat_AquoMetadataLocatieLijst_locidx.loc[cat_locatielijst_containingstation.index,'AquoMetaData_MessageID']
        bool_aquometadatalijst_containingstation.loc[bool_aquometadatalijst_containingstation_idxtrue] = True
    
    if meta_dict is not None and station is not None: 
        cat_aquometadatalijst_sel = cat_aquometadatalijst.loc[bool_aquometadatalijst_containingmeta & bool_aquometadatalijst_containingstation]
        cat_locatielijst_sel = cat_locatielijst.loc[bool_locatielijst_containingmeta & bool_locatielijst_containingstation]
    elif meta_dict is None and station is not None:
        cat_aquometadatalijst_sel = cat_aquometadatalijst.loc[bool_aquometadatalijst_containingstation]
        cat_locatielijst_sel = cat_locatielijst.loc[bool_locatielijst_containingstation]
    elif meta_dict is not None and station is None:
        cat_aquometadatalijst_sel = cat_aquometadatalijst.loc[bool_aquometadatalijst_containingmeta]
        cat_locatielijst_sel = cat_locatielijst.loc[bool_locatielijst_containingmeta]

    if error_empty and len(cat_aquometadatalijst_sel)==0 and len(cat_locatielijst_sel)==0:
        raise Exception('Catalog query yielded no results')
        
    return cat_aquometadatalijst_sel, cat_locatielijst_sel




