# -*- coding: utf-8 -*-

"""
Console script for hatyan.
    - ``hatyan --help``
"""
import sys
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
@click.option('-v', '--verbose', is_flag=True, help="set logging level to debug instead of info")
@click.version_option(hatyan.__version__)
def cli(filename, unique_outputdir, interactive_plots, redirect_stdout, 
        verbose
        ):
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
        sys.stdout = open('STDOUT.txt', "w")
    
    level = logging.INFO
    if verbose:
        level = logging.DEBUG
    logging.basicConfig(level=level, stream=sys.stdout)
    
    
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
    # therefore we use subprocess instead, this also requires flushing the print buffer first
    # with open(file_config) as f:
    #     exec(f.read())
    spec = importlib.util.spec_from_file_location("script", file_config)
    foo = importlib.util.module_from_spec(spec)
    sys.modules["module.name"] = foo
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
        # logger.warning(f"close open plots to continue (mpl.backend='{matplotlib.get_backend()}')")
        plt.show()
    os.chdir("..")

