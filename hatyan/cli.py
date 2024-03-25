# -*- coding: utf-8 -*-

"""
Console script for hatyan.
    - ``hatyan --help``
"""
import sys
# import logging
import click
import hatyan
import os
import shutil
import datetime as dt
import matplotlib
import matplotlib.pyplot as plt

@click.command()
@click.argument(
    'filename',
    #help="Python script to run, will be copied to `dir_output`"
)
@click.option(
    "-u", "--unique-outputdir",
    is_flag=True,
    help="add timestamp to `dir_output` so output is never overwritten"
)
@click.option(
    "-i", "--interactive-plots",
    is_flag=True,
    help="show interactive plots at the end of the hatyan process"
)
@click.option(
    "-r", "--redirect-stdout",
    is_flag=True,
    help="redirecting stdout to dir_output/STDOUT.txt, "
    "warnings/errors are still printed to console"
)
# @click.option('-v', '--verbose', is_flag=True, help="this is not used")
@click.version_option(hatyan.__version__)
def cli(filename, unique_outputdir, interactive_plots, redirect_stdout, 
        # verbose
        ):
    """
    Initializes hatyan by creating a `dir_output`, setting current and matplotlib
    savefig directory to `dir_output` and printing the initialisation header.
    Then runs the provided configfile (FILENAME).
    Wraps up hatyan by printing a de-initialisation footer with the script runtime.
    """
    # level = logging.INFO
    # if verbose:
    #     level = logging.DEBUG
    # logging.basicConfig(level=level)
    
    # get file_config and start time
    file_config = os.path.abspath(filename)
    timer_start = dt.datetime.now()
    
    # check for file existence before creating folder
    if not os.path.isfile(file_config):
        raise FileNotFoundError(f"file_config not found: {file_config}")
    
    # create dir_output and chdir
    foldername = os.path.splitext(os.path.basename(file_config))[0]
    if unique_outputdir:
        dir_output = f"{foldername}__{timer_start.strftime('%Y%m%d_%H%M%S')}"
    else:
        dir_output = foldername
    os.makedirs(dir_output, exist_ok=True)
    os.chdir(dir_output)
    
    #copy configfile to cwd (=dir_output)
    shutil.copyfile(file_config,os.path.basename(file_config))
    
    #set the storage location of interactive plots
    matplotlib.rcParams["savefig.directory"] = dir_output
    
    #redirecting stdout, stderr is still printed to console
    if redirect_stdout:
        # set sys.stdout to redirect prints in this command
        sys.stdout = open('STDOUT.txt', 'w')
        # set stdout to redirect prints in actual execution of the config_file
        stdout = sys.stdout
    else:
        stdout = None
    
    # initialization print
    print("############### HATYAN INITALIZING ###############")
    print("--------------------------------------------------")
    print(f"hatyan-{hatyan.__version__}: RWS tidal analysis and prediction")
    print(f"started at:  {timer_start.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"file_config: {file_config}")
    print(f"dir_output:  {dir_output}")
    print("--------------------------------------------------")
    
    # write above prints to file even if buffer it not yet full
    sys.stdout.flush()
    
    # run the configfile
    # we use subprocess instead of exec since this supports oneline list generations
    import subprocess
    file_config_abs = os.path.abspath(file_config)
    p = subprocess.run(f"{sys.executable} {file_config_abs}", stdout=stdout)
    if p.returncode:
        raise RuntimeError("hatyan run failed, check error messages above")
    
    # get stop time
    timer_stop = dt.datetime.now()
    timer_elapsed = (timer_stop - timer_start).total_seconds()/60
    
    # de-initialization print
    print("--------------------------------------------------")
    print(f"finished at: {timer_stop.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"elapsed time: {timer_elapsed:.2f} minutes {timer_elapsed*60:.2f} seconds")
    print("--------------------------------------------------")
    print("################# HATYAN FINISHED ################")
    
    # show all created plots in case of interactive-plots
    if interactive_plots:
        print(f"WARNING: close open plots to continue (mpl.backend='{matplotlib.get_backend()}')")
        plt.show()
    os.chdir("..")

