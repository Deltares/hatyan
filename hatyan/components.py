# -*- coding: utf-8 -*-
"""
components.py contains all the definitions related to hatyan components.

hatyan is a Python program for tidal analysis and prediction, based on the FORTRAN version. 
Copyright (C) 2019-2021 Rijkswaterstaat.  Maintained by Deltares, contact: Jelmer Veenstra (jelmer.veenstra@deltares.nl). 
Source code available at: https://github.com/Deltares/hatyan

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
import datetime as dt
import pytz
import warnings

from hatyan.schureman import get_schureman_freqs, get_schureman_v0 #TODO: this is not generic foreman/schureman
from hatyan.hatyan_core import sort_const_list, get_const_list_hatyan
from hatyan.metadata import (metadata_add_to_obj, metadata_from_obj, 
                             metadata_compare, wns_from_metadata)


def plot_components(comp, comp_allperiods=None, comp_validation=None, sort_freqs=True):
    """
    Create a plot with the provided analysis results

    Parameters
    ----------
    comp : TYPE
        DESCRIPTION.
    comp_allyears : TYPE, optional
        DESCRIPTION. The default is None.
    comp_validation : TYPE, optional
        DESCRIPTION. The default is None.
    sort_freqs : BOOL, optional
        Whether to sort the component list on frequency or not, without sorting it is possible to plot components that are not available in hatyan. The default is True.
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The generated figure handle, with which the figure can be adapted and saved.
    axs : (tuple of) matplotlib.axes._subplots.AxesSubplot
        The generated axis handle, whith which the figure can be adapted.

    """
    
    COMP = comp.copy()
    if comp_allperiods is not None:
        comp_allperiods = comp_allperiods.copy()
        comp_legend_labels = comp_allperiods['A'].columns
    
    if comp_validation is not None:
        COMPval = comp_validation.rename(columns={"A": "A_validation", "phi_deg": "phi_deg_validation"})
        COMP = pd.concat([COMP,COMPval],axis=1)
        COMP['A_diff'] = COMP['A']-COMP['A_validation']
        COMP['phi_deg_diff'] = (COMP['phi_deg']-COMP['phi_deg_validation']+180)%360-180
    
    const_list = COMP.index.tolist() #is combined const_list
    if sort_freqs:
        const_list = sort_const_list(const_list=const_list)
        COMP = COMP.loc[const_list] #sorted
    
    size_figure = (15,9)
    size_line_ts = 0.7
    size_line_comp = 1.2
    size_marker_comp = 4
    
    fig, (ax1, ax2) = plt.subplots(2,1,figsize=size_figure, sharex=True)
    ax1.set_title('Amplitudes and Phases per component')
    ax1.plot([0,len(const_list)],[0,0],'-k',linewidth=size_line_ts)
    if comp_allperiods is not None: #and COMP_validation is None:
        ax1.plot(comp_allperiods['A'].values,'o-',color='gray',linewidth=size_line_comp,markersize=size_marker_comp, label=comp_legend_labels)
    ax1.plot(COMP['A'],'-o',linewidth=size_line_comp,markersize=size_marker_comp,label='comp')
    if comp_validation is not None:
        ax1.plot(COMP['A_validation'],'o-',linewidth=size_line_comp,markersize=size_marker_comp,label='comp_validation')
        ax1.plot(COMP['A_diff'],'go-',linewidth=size_line_comp,markersize=size_marker_comp,label='difference')
    ax1.grid(axis='y')#, which='major')
    ax1.set_xlim(-0.5,len(const_list)+0.5)
    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax1.set_ylim(-.05,1)
    ax1.set_ylabel('amplitudes [m]')
    ax1.legend(loc='lower right')
    ax2.plot([0,len(const_list)],[0,0],'-k',linewidth=size_line_ts)
    if comp_allperiods is not None:
        ax2.plot(comp_allperiods['phi_deg'],'o-',color='gray',linewidth=size_line_comp,markersize=size_marker_comp, label=comp_legend_labels)
    ax2.plot(COMP['phi_deg'],'-o',linewidth=size_line_comp,markersize=size_marker_comp,label='comp')
    if comp_validation is not None:
        ax2.plot(COMP['phi_deg_validation'],'o-',linewidth=size_line_comp,markersize=size_marker_comp,label='comp_validation')
        ax2.plot(COMP['phi_deg_diff'],'go-',linewidth=size_line_comp,markersize=size_marker_comp,label='difference')
    ax2.grid(axis='y')#, which='major')
    ax2.set_xlim(-0.5,len(const_list)-0.5)
    #ax2.set_xticklabels([]) #remove major xtick labels
    ax2.set_xticks(np.arange(len(const_list)))#,minor=True) #set minor xticks to 1 (create tick for each constituent, but no grid)
    ax2.set_xticklabels(const_list,rotation=90)#,minor=True) #set minor xticklabels, each constituent
    ax2.set_xlabel('constituents')
    ax2.set_ylim(-5,360+5)
    ax2.set_yticks(np.arange(0,360+.001,90))
    ax2.set_ylabel('phases [degrees]')
    ax2.legend(loc='lower right')
    fig.tight_layout()
    
    axs = (ax1,ax2)
    return fig, axs


def write_components(comp, filename):
    """
    Writes the provided analysis results to a file

    Parameters
    ----------
    comp : TYPE
        DESCRIPTION.
    filename : TYPE
        DESCRIPTION.
    metadata : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    None.

    """
    
    #get metadata before copying DataFrame
    metadata = metadata_from_obj(comp)
    waarnemingssoort = wns_from_metadata(metadata)
    
    station = metadata.pop('station')
    grootheid = metadata.pop('grootheid')
    vertref = metadata.pop('vertref')
    unit = metadata.pop('eenheid')
    
    tstart = metadata.pop('tstart')
    tstop = metadata.pop('tstop')
    tzone_min = metadata.pop('tzone')._minutes
    tstart_str = tstart.strftime("%Y%m%d  %H%M%S")
    tstop_str = tstop.strftime("%Y%m%d  %H%M%S")
    
    nodalfactors = metadata.pop('nodalfactors')
    fu_alltimes = metadata.pop('fu_alltimes')
    if not nodalfactors:
        raise ValueError("nodalfactors=False cannot be distinguished in components file header, so not supported by write_components()")
    if fu_alltimes:
        raise ValueError("fu_alltimes=True cannot be distinguished in components file header, so not supported by write_components()")
        
    if 'A0' in comp.index.tolist():
        midd = comp.loc['A0','A']*100
        comp = comp.drop('A0',axis=0)
    else:
        midd = 0
        comp = comp.copy()

    #sort components by frequency
    t_const_freq_pd = get_schureman_freqs(comp.index.tolist())
    comp['freq'] = t_const_freq_pd['freq']
    comp = comp.sort_values(by='freq')
    
    const_list_hatyan195_orig = get_const_list_hatyan('all_schureman_originalorder')
    comp['const_no'] = [const_list_hatyan195_orig.index(x) for x in comp.index]  
    comp['const_speed'] = t_const_freq_pd['freq'].values*360  
    ncomp = len(comp.index)
    
    with open(filename,'w') as f:
        
        """
        header of component file is described in https://repos.deltares.nl/repos/lib_tide/trunk/src/hatyan_fortran/HATYAN40/anadea.f
        
        STAT:
        Locatiecode station
        Parametercode
        Hoedanigheidcode
        Eenheidcode
        Waarnemingssoort
        
        PERD:
        BEGINDATUM
        BEGINTIJD
        EINDDATUM
        EINDTIJD
        TIJDVERSCHIL T.O.V. GMT
        
        CODE: COMPONENTENCODE (?)
        
        NCOM: AANTAL COMPONENTEN
        """
        
        from hatyan import __version__ as hatyan_version
        f.write(f'* written with hatyan-{hatyan_version}\n')            
        
        if metadata is None: #TODO: maybe remove this or add metadata retrieval in try-except loop
            f.write('* no metadata available\n') #TODO: HATYAN40\anadea.f regel 297 schrijft format van header voor, '60' na start/stop datetime is tijdzone, '1' na eenheid is Waarnemingssoort. Beide belangrijke regels/gegevens
        else:
            for key in metadata.keys():
                f.write(f'* {key} :  {metadata[key]}\n')
        
        f.write(f'STAT  {station}    {grootheid}    {vertref}    {unit}    {waarnemingssoort}\n')
        f.write(f'PERD  {tstart_str}  {tstop_str}     {tzone_min}\n')
        f.write( 'CODE      3\n')
        f.write(f'MIDD {midd:9.2f}\n')
        f.write(f'NCOM {ncomp:5d}\n')
        for compname in comp.index.tolist():
            comp_one = comp.loc[compname]
            f.write("COMP %4i %12.6f %9.3f %7.2f  %-12s\n" % (comp_one['const_no'],
                                                              comp_one['const_speed'], 
                                                              comp_one['A']*100, 
                                                              comp_one['phi_deg']%360, 
                                                              compname))


def merge_componentgroups(comp_main, comp_sec, comp_sec_list=['SA','SM']):
    """
    Merges the provided component groups into one

    Parameters
    ----------
    comp_main : TYPE
        DESCRIPTION.
    comp_sec : TYPE
        DESCRIPTION.
    comp_sec_list : TYPE, optional
        DESCRIPTION. The default is ['SA','SM'].

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    COMP_merged : TYPE
        DESCRIPTION.

    """
    
    comp_main_meta = metadata_from_obj(comp_main).copy()
    comp_sec_meta = metadata_from_obj(comp_sec).copy()
    
    meta_settings_list = ['origin','nodalfactors','xfac','fu_alltimes','groepering','timestep_min','timestep_unit','TYP']
    comp_main_meta_others = {}
    for key in meta_settings_list:
        if key in comp_main_meta:
            comp_main_meta_others[f'{key}'] = comp_main_meta.pop(key)
        if key in comp_sec_meta:
            comp_main_meta_others[f'{key}_sec'] = comp_sec_meta.pop(key)
    metadata_compare([comp_main_meta,comp_sec_meta])
    
    # add metadata for analysis settings
    comp_merged_meta = comp_main_meta.copy()
    for key in meta_settings_list:
        if key in comp_main_meta_others:
            comp_merged_meta[key] = comp_main_meta_others[key]
    
    #add metadata for secondary components (list + origin)
    comp_sec_str = ', '.join(comp_sec_list)
    components_sec_str = f"{comp_sec_str} imported"
    if 'origin_sec' in comp_main_meta_others:
        origin_sec = comp_main_meta_others['origin_sec']
        components_sec_str += ' '+origin_sec
    comp_merged_meta['components_sec'] = components_sec_str
    
    comp_merged = comp_main.copy()
    
    comp_sec_list_sel = [comp for iC,comp in enumerate(comp_merged.index) if comp in comp_sec_list]
    if comp_sec_list_sel != []:
        comp_merged = comp_merged.drop(comp_sec_list_sel)
    comp_merged = pd.concat([comp_sec.loc[comp_sec_list],comp_merged])

    t_const_freq = get_schureman_freqs(comp_merged.index.tolist())
    comp_merged['freq'] = t_const_freq['freq']
    comp_merged = comp_merged.sort_values(by='freq')
    
    # add metadata
    comp_merged = metadata_add_to_obj(comp_merged, comp_merged_meta)
    
    return comp_merged


def read_components(filename):
    """
    Reads analysis results from a file.

    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.
    
    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    
    print('reading file: %s'%(filename))

    xfac = True
    line_compstart = None
    metadata_available = False
    with open(filename) as f:
        for i, line in enumerate(f):
            if line.startswith('*'):
                if 'theoretische' in line:
                    print(f"xfac=False derived from headerline: '{line.strip()}'")
                    xfac = False
                elif 'empirisch' in line:
                    print(f"xfac=True derived from headerline: '{line.strip()}'")
                    xfac = True
            elif line.startswith('STAT'):
                station = line.split()[1]
                print('retrieving data from components file for station %s'%(station))
                grootheid = line.split()[2]
                vertref = line.split()[3]
                print('the vertical reference level in the imported file is: %s'%(vertref))
                eenheid = line.split()[4]
                # waarnemingssoort = int(line.split()[5])
                metadata_available = True
            elif line.startswith('PERD'):
                dateline = line.split()
                tstart = pd.Timestamp(dateline[1]+' '+dateline[2])
                tstop = pd.Timestamp(dateline[3]+' '+dateline[4])
                tzone = pytz.FixedOffset(int(dateline[5]))
                metadata_available = True
            elif line.startswith('MIDD'):
                A0_cm = float(line.split()[1])
            elif line.startswith('COMP'):
                line_compstart = i
                break #break because last line before actual data
    
    
    #retrieve raw data
    if line_compstart is None:
        raise Exception('invalid file, no line that starts with COMP')
    Aphi_datapd_raw_noA0 = pd.read_csv(filename, delimiter=r"\s+", header=line_compstart-1, names=['COMP', 'hat_id', 'freq', 'A', 'phi', 'name'])
    
    if "A0" in Aphi_datapd_raw_noA0["name"].tolist():
        bool_a0 = Aphi_datapd_raw_noA0["name"] == "A0"
        Aphi_datapd_raw_noA0 = Aphi_datapd_raw_noA0.loc[~bool_a0]
        print(UserWarning("A0 component found in components file, the file is probably generated with hatyan 2.7.0 or older. Dropping the component since the data is available in MIDD."))
    
    Aphi_datapd_A0line = pd.DataFrame({'A': [A0_cm], 'phi': [0], 'name': ['A0']})
    Aphi_datapd_raw = pd.concat([Aphi_datapd_A0line,Aphi_datapd_raw_noA0],ignore_index=True)
    comp_pd = pd.DataFrame({'A': Aphi_datapd_raw['A'].values/100, 'phi_deg': Aphi_datapd_raw['phi'].values}, index=Aphi_datapd_raw['name'].values)
    
    # add metadata
    if not metadata_available:
        warnings.warn('No metadata available in component file (STAT/PERD lines), this might cause issues in the rest of the process. You are probably using a component file generated with hatyan 2.7.0 or older.')
        return comp_pd
    
    metadata = {'station':station,
                'grootheid':grootheid, 'eenheid':eenheid,
                'vertref':vertref,
                'tstart':tstart, 'tstop':tstop, 'tzone':tzone, 
                'nodalfactors':True, 'xfac':xfac,
                'origin':'from component file'}
    comp_pd = metadata_add_to_obj(comp_pd, metadata)
    return comp_pd


def components_timeshift(comp,hours):
        
    comp_out = comp.copy()
    
    refdate = dt.datetime(2000,1,1)
    corrdate = refdate+dt.timedelta(hours=hours)
    
    v0_twodates = get_schureman_v0(const_list=comp.index, dood_date=pd.DatetimeIndex([refdate,corrdate]))
    hourcorr_v0_deg = np.rad2deg(v0_twodates[1]-v0_twodates[0])#%360
    
    comp_out['phi_deg'] = (comp_out['phi_deg']+hourcorr_v0_deg)%360
    
    # add metadata
    metadata = metadata_from_obj(comp)
    metadata['tzone'] = pytz.FixedOffset(metadata['tzone']._minutes + hours*60)
    comp_out = metadata_add_to_obj(comp_out,metadata)
    
    return comp_out


