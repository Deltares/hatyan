# -*- coding: utf-8 -*-
"""
wrapper_RWS.py contains wrapper functions around the hatyan process for RWS related calculations.

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

import os
import sys
import shutil
import datetime as dt
import matplotlib
import matplotlib.pyplot as plt


def init_RWS(interactive_plots=True):
    """
    Initializes the hatyan process for RWS related calculations. Besides the return variables,
    it prints a header for the print output (shows up in the hatyan diagnostics file)

    Parameters
    ----------
    interactive_plots : bool/int, optional
        sets the correct matplotlib backend so plots are (not) displayed on both RedHat and windows. The default is True.

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    dir_output : path/str
        the output directory for the hatyan process, the current directory is set to this folder.
    timer_start : datetime.datetime
        provides a start time with which exit_RWS calculates the total time of the process.

    """
    
    argvlist = sys.argv
    
    file_config = os.path.realpath(argvlist[0])
    
    if len(argvlist) == 1:
        dir_output = get_outputfoldername(file_config)
    elif len(argvlist) == 2: #for running testbank with command `python configfile.py dir_output`
        dir_output = argvlist[1]
    else:
        raise Exception('ERROR: something wrong with input arguments')
    os.chdir(dir_output)
    
    with open('__NOT_FINISHED__','w') as f:
        f.write('CAUTION, this hatyan process has not yet properly finished, check FILE_DIAGNOSICS.txt for possible errors')

    try:
        import hatyan
        version_no = hatyan.__version__
    except ModuleNotFoundError:
        version_no = None
    
    #set the storage location of interactive plots
    #TODO: this is not necessary since cwd is already dir_output
    import matplotlib
    matplotlib.rcParams["savefig.directory"] = dir_output
    
    #set the matplotlib backend depending on the interactive_plots argument
    import matplotlib.pyplot as plt
    if interactive_plots:
        try:
            plt.switch_backend('Qt5agg')
        except:
            raise Exception('Failed to switch to Qt5agg backend, check if you have X-forwarding enabled (and mesa-libGL and possibly other libraries installed) or use argument interactive_plots=False')
    else:
        plt.switch_backend('Agg')
    
    ##################################################################
    ##### fixed settings #############################################
    ##################################################################
    
    timer_start = dt.datetime.now()
    sys.argv.append(timer_start)
    print('#'*50)
    print('-'*50)
    print('hatyan-%s: RWS tidal analysis and prediction'%(version_no))
    print('-'*50)
    print('INITIALISATION')
    print('started at %s'%(timer_start.strftime('%Y-%m-%d %H:%M:%S')))
    print('%-45s = %s'%('file_config',file_config))
    print('%-45s = %s'%('dir_output',dir_output))
    print('-'*50)
        
    #copy configfile to dir_output
    print('copying configfile to dir_output\\%s'%(os.path.basename(file_config)))
    shutil.copy(file_config,dir_output)
    print('END OF INITIALISATION')
    


def exit_RWS():
    """
    Provides a footer to the print output (shows up in the hatyan diagnostics file)

    Parameters
    ----------
    timer_start : TYPE
        The start time of the hatyan process, which is used to calculate the total time of the process.

    Returns
    -------
    None.

    """

    if matplotlib.get_backend() == 'Qt5Agg':
        print('MESSAGE: interactive plots opened, close them to continue')
        plt.show()
    
    timer_start = sys.argv[-1]
    if not isinstance(timer_start,dt.datetime):
        raise Exception('exit_RWS() can only be called if init_RWS() is called in the same process')
    print('-'*50)
    timer_stop = dt.datetime.now()
    print('calculation finished at %s'%(timer_stop.strftime('%Y-%m-%d %H:%M:%S')))
    timer_elapsed = (timer_stop - timer_start).total_seconds()/60
    print('elapsed time: %.2f minutes (%.2f seconds)'%(timer_elapsed,timer_elapsed*60))
    print('-'*50)
    
    os.remove('__NOT_FINISHED__')


def get_outputfoldername(file_config):
    """
    Creates an output folder based on the start time of the filename of the configfile and the current time. 

    Parameters
    ----------
    file_config : str or path
        path to the configuration file.

    Returns
    -------
    dir_output : str or path
        path to the output directory.

    """
    
    mode = os.path.basename(file_config).split('.')[0]
    time_now = dt.datetime.now()
    dir_output = os.path.join(os.getcwd(),'output__%s__%s'%(time_now.strftime('%Y%m%d_%H%M%S'),mode))
    
    os.makedirs(dir_output, exist_ok=False)
    
    return dir_output


def close(fig=None):
    """
    wrapper around matplotlib.pyplot.close() to enable closing of plots generated by hatyan without explicitly importing matplotlib

    Parameters
    ----------
    fig : None or int or str or .Figure
        The figure to close. There are a number of ways to specify this:
        None: the current figure
        .Figure: the given .Figure instance
        int: a figure number
        str: a figure name
        'all': all figures.

    Returns
    -------
    None.

    """
    
    plt.close(fig)
    