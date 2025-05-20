# -*- coding: utf-8 -*-

"""
Console script for hatyan.
    - ``hatyan --help``
"""

import logging
import click
import hatyan
import os
import shutil
import datetime as dt
import matplotlib
import matplotlib.pyplot as plt
import importlib

logger = logging.getLogger(__name__)


@click.command()
@click.argument(
    'filename',
    type=click.Path(exists=True),
)
@click.option(
    "-o", "--overwrite",
    is_flag=True,
    help="overwrite `dir_output` if it exists"
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
@click.option(
    '-l',
    '--loglevel',
    type=str,
    help="set logging level to ERROR/WARNING/INFO/DEBUG, default is INFO",
)
@click.version_option(hatyan.__version__)
def cli(filename, overwrite, interactive_plots, redirect_stdout, loglevel):
    """
    Initializes hatyan by creating a `dir_output`, setting current and matplotlib
    savefig directory to `dir_output` and printing the initialisation header.
    Then runs the provided configfile (FILENAME).
    Wraps up hatyan by printing a de-initialisation footer with the script runtime.
    """
    
    # get file_config and start time
    file_config = os.path.abspath(filename)
    timer_start = dt.datetime.now()
    
    # create dir_output and chdir
    dir_output = os.path.splitext(os.path.basename(file_config))[0]
    if os.path.exists(dir_output) and not overwrite:
        raise FileExistsError(f"Directory '{dir_output}' already exists. Use '--overwrite' if you want to overwrite existing results.")
    os.makedirs(dir_output, exist_ok=True)
    os.chdir(dir_output)
    
    #copy configfile to cwd (=dir_output)
    shutil.copyfile(file_config, os.path.basename(file_config))
    
    #set the storage location of interactive plots
    matplotlib.rcParams["savefig.directory"] = dir_output
    
    #redirecting stdout, stderr is still printed to console
    if redirect_stdout:
        stream = open('STDOUT.txt', "w")
    else:
        stream = None
    
    # set logging level and stdout
    if loglevel is None:
        loglevel = "INFO"
    logging.basicConfig(level=loglevel, stream=stream)
    
    # initialization print
    logger.info("############### HATYAN INITALIZING ###############")
    logger.info("--------------------------------------------------")
    logger.info(f"hatyan-{hatyan.__version__}: RWS tidal analysis and prediction")
    logger.info(f"started at:  {timer_start.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"file_config: {file_config}")
    logger.info(f"dir_output:  {dir_output}")
    logger.info("--------------------------------------------------")
    
    # run the configfile
    # exec from within cli somehow does not support oneline list generation with predefined variables
    # with open(file_config) as f:
    #     exec(f.read())
    # therefore we use importlib instead: https://stackoverflow.com/questions/67631/how-can-i-import-a-module-dynamically-given-the-full-path
    spec = importlib.util.spec_from_file_location("arbitrary.name", file_config)
    foo = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(foo)
    
    # get stop time
    timer_stop = dt.datetime.now()
    timer_elapsed = (timer_stop - timer_start).total_seconds()/60
    
    # de-initialization print
    logger.info("--------------------------------------------------")
    logger.info(f"finished at: {timer_stop.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"elapsed time: {timer_elapsed:.2f} minutes {timer_elapsed*60:.2f} seconds")
    logger.info("--------------------------------------------------")
    logger.info("################# HATYAN FINISHED ################")
    
    # show all created plots in case of interactive-plots
    if interactive_plots:
        logger.warning(f"close open plots to continue (mpl.backend='{matplotlib.get_backend()}')")
        plt.show()
    os.chdir("..")

