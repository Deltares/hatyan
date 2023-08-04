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
import argparse
import shutil
import datetime as dt
import matplotlib
import matplotlib.pyplot as plt


def init_RWS():
    """
    Initializes the hatyan process by creating a `dir_output` if it does not exist,
    setting current and matplotlib savefig directory to `dir_output`,
    setting the matplotlib backend to Qt5agg if --interactive-plots sysargv is used,
    starting a script timer for exit_RWS(),
    copying the input script to `dir_output`,
    and printing the initialisation header

    Parameters
    ----------
    interactive_plots : bool, optional
        sets the correct matplotlib backend so plots are (not) displayed on both RedHat and windows. The default is None.

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    #parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-D', '--dir-output') # option that takes a value, default is None
    parser.add_argument('-I', '--interactive-plots', action='store_true') # on/off flag, default is False
    args = parser.parse_args()
    
    file_config = os.path.realpath(sys.argv[0])
    
    if not os.path.isfile(file_config): # escape for running with F9 or unsaved script
        print('init_RWS() silently failed, file_config not found')
        return
    
    if args.dir_output is None:
        foldername = os.path.basename(file_config).split('.')[0]
        time_now = dt.datetime.now()
        dir_output = 'output__%s__%s'%(time_now.strftime('%Y%m%d_%H%M%S'),foldername)
        os.makedirs(dir_output, exist_ok=False)
    else: # for running testbank with command `python configfile.py dir_output`
        dir_output = args.dir_output
    os.chdir(dir_output)
    
    #set the storage location of interactive plots
    import matplotlib
    matplotlib.rcParams["savefig.directory"] = dir_output

    #set the matplotlib backend depending on the interactive_plots argument
    import matplotlib.pyplot as plt
    if args.interactive_plots:
        try:
            plt.switch_backend('Qt5agg')
        except:
            raise Exception('Failed to switch to Qt5agg backend, check if you have X-forwarding enabled (and mesa-libGL and possibly other libraries installed) or use argument interactive_plots=False')
    
    # get hatyan version and start time
    import hatyan
    hatyan_version = hatyan.__version__
    timer_start = dt.datetime.now()
    timer_start_str = timer_start.strftime('%Y-%m-%d %H:%M:%S')
    # appending to sysargv list, environment variable would be more suitable that remains after process finishes.
    sys.argv.append(timer_start)
    
    # print initialization
    print('############### HATYAN INITALIZING ###############')
    print('-'*50)
    print(f'hatyan-{hatyan_version}: RWS tidal analysis and prediction')
    print(f'started at:  {timer_start_str}')
    print(f'file_config: {file_config}')
    print(f'dir_output:  {dir_output}')
    #copy configfile to dir_output
    shutil.copy(file_config,dir_output)
    print('-'*50)


def exit_RWS():
    """
    Provides a finishing footer to the print output

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    timer_start = sys.argv[-1]
    if not isinstance(timer_start,dt.datetime):
        raise Exception('exit_RWS() can only be called if init_RWS() is called in the same process')

    timer_stop = dt.datetime.now()
    timer_stop_str = timer_stop.strftime('%Y-%m-%d %H:%M:%S')
    timer_elapsed = (timer_stop - timer_start).total_seconds()/60

    # print de-initialization
    print('-'*50)
    if matplotlib.get_backend() == 'Qt5Agg':
        print('WARNING: close open plots to continue')
        plt.show()
    print(f'finished at: {timer_stop_str}')
    print(f'elapsed time: {timer_elapsed:.2f} minutes {timer_elapsed*60:.2f} seconds')
    print('-'*50)
    print('################# HATYAN FINISHED ################')
    
    with open('FINISHED.txt','w') as f:
        f.write('The presence of this file means the run finished succesfully. In a later stage some logging will be printed in this file.')


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
    