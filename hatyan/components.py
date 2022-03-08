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

from hatyan.schureman import get_schureman_freqs, get_schureman_v0 #TODO: this is not generic foreman/schureman
from hatyan.hatyan_core import sort_const_list, get_const_list_hatyan


def plot_components(comp, comp_allyears=None, comp_validation=None, sort_freqs=True):
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
    if comp_allyears is not None:
        comp_allyears = comp_allyears.copy()
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
    if comp_allyears is not None: #and COMP_validation is None:
        ax1.plot(comp_allyears['A'].values,'o-',color='gray',linewidth=size_line_comp,markersize=size_marker_comp, label='separate years')
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
    if comp_allyears is not None:
        ax2.plot(comp_allyears['phi_deg'],'o-',color='gray',linewidth=size_line_comp,markersize=size_marker_comp, label='separate years')
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
    


    
    
def write_components(comp, filename, metadata=None):
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
    
    COMP = comp.copy()
    
    t_const_freq_pd = get_schureman_freqs(COMP.index.tolist())
    COMP['freq'] = t_const_freq_pd['freq']
    COMP = COMP.sort_values(by='freq')

    const_list_hatyan195_orig = get_const_list_hatyan('all_schureman_originalorder')
    const_no = [const_list_hatyan195_orig.index(x) for x in COMP.index]
    const_speed = t_const_freq_pd['freq'].values*360
    
    with open(filename,'w') as f:
        if metadata is None:
            f.write('* no metadata available\n')
        else:
            for key in metadata.keys():
                f.write('* %-20s: %s\n'%(key, metadata[key]))

        if 'A0' in COMP.index.tolist():
            f.write('MIDD   %.2f cm\n'%(COMP.loc['A0','A']*100))
        else:
            f.write('MIDD   %.2f cm\n'%(np.nan))
        f.write('NCOM   %i\n'%(len(COMP.index)))
        for iC, compname in enumerate(COMP.index.tolist()):
            f.write("COMP %4i %12.6f %9.3f %7.2f  %-12s\n" % (const_no[iC], const_speed[iC], COMP.loc[compname,'A']*100, COMP.loc[compname,'phi_deg']%360, compname))


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
    
    COMP_merged = comp_main.copy()
    
    comp_sec_list_sel = [comp for iC,comp in enumerate(COMP_merged.index) if comp in comp_sec_list]
    if comp_sec_list_sel != []:
        COMP_merged = COMP_merged.drop(comp_sec_list_sel)
    COMP_merged = comp_sec.loc[comp_sec_list].append(COMP_merged)

    t_const_freq = get_schureman_freqs(COMP_merged.index.tolist())
    COMP_merged['freq'] = t_const_freq['freq']
    COMP_merged = COMP_merged.sort_values(by='freq')
    
    return COMP_merged


def read_components(filename, get_metadata=False):
    """
    Reads analysis results from a file.

    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.
    get_metadata : TYPE, optional
        DESCRIPTION. The default is False.

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

    file = open(filename)
    line_compstart = None
    for i, line in enumerate(file):
        if line.startswith('STAT'):
            station_fromfile = line.split()[1]
            print('retrieving data from components file for station %s'%(station_fromfile))
            file_vertref = line.split()[3]
            print('the vertical reference level in the imported file is: %s'%(file_vertref))
        elif line.startswith('PERD'):
            dateline = line.split()
            times_compfile_ext = [dt.datetime.strptime(dateline[1]+dateline[2],'%Y%m%d%H%M'),dt.datetime.strptime(dateline[3]+dateline[4],'%Y%m%d%H%M')]
            times_compfile_step = int(dateline[5])
        elif line.startswith('MIDD'):
            A0_cm = float(line.split()[1])
        elif line.startswith('COMP'):
            line_compstart = i
            break #break because last line before actual data
    #retrieve raw data
    if line_compstart is None:
        raise Exception('invalid file, no line that starts with COMP')
    Aphi_datapd_raw_noA0 = pd.read_csv(filename, delimiter=r"\s+", header=line_compstart-1, names=['COMP', 'hat_id', 'freq', 'A', 'phi', 'name'])
    
    Aphi_datapd_A0line = pd.DataFrame({'A': [A0_cm], 'phi': [0], 'name': ['A0']})
    Aphi_datapd_raw = Aphi_datapd_A0line.append(Aphi_datapd_raw_noA0,ignore_index=True)
    COMP_pd = pd.DataFrame({'A': Aphi_datapd_raw['A'].values/100, 'phi_deg': Aphi_datapd_raw['phi'].values}, index=Aphi_datapd_raw['name'].values)
    if get_metadata:
        meta = {'station': station_fromfile, 'times_ext': times_compfile_ext, 'times_stepmin': times_compfile_step, 'origin':'import', 'vertref':file_vertref}#, 'usedxfac':None}
        return COMP_pd, meta
    else:
        return COMP_pd


def components_timeshift(comp,hours):
        
    comp_out = comp.copy()
    
    refdate = dt.datetime(2000,1,1)
    corrdate = refdate+dt.timedelta(hours=hours)
    
    v0_twodates = get_schureman_v0(const_list=comp.index, dood_date=pd.DatetimeIndex([refdate,corrdate]))
    hourcorr_v0_deg = np.rad2deg(v0_twodates[1]-v0_twodates[0])#%360
    
    comp_out['phi_deg'] = (comp_out['phi_deg']+hourcorr_v0_deg)%360
    
    return comp_out


