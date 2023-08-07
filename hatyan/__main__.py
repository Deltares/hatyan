# -*- coding: utf-8 -*-
#print('aa')


import os
import sys
import argparse
import shutil
import datetime as dt
import matplotlib
import matplotlib.pyplot as plt
import hatyan


# get start time
timer_start = dt.datetime.now()
timer_start_str = timer_start.strftime('%Y-%m-%d %H:%M:%S')

hatyan_version = hatyan.__version__
hatyan_description = (
    "Initializes hatyan by creating a `dir_output`, "
    "setting current and matplotlib savefig directory to `dir_output` "
    "and printing the initialisation header. "
    "Then runs the provided configfile. "
    "Wraps up hatyan by printing a finishing footer with the script runtime."
)

#parse arguments: https://docs.python.org/3/library/argparse.html
parser = argparse.ArgumentParser(prog='hatyan', description=hatyan_description)
parser.add_argument('filename',
                    help='Python script to run, will be copied to `dir_output`') # positional argument
parser.add_argument('-u', '--unique-outputdir', action='store_true',
                    help='add timestamp to `dir_output` so output is never overwritten')  # on/off flag, default is False
parser.add_argument('-i', '--interactive-plots', action='store_true',
                    help='manually setting the matplotlib backend to Qt5agg') # on/off flag, default is False
parser.add_argument('-r', '--redirect-stdout', action='store_true',
                    help='redirecting stdout to dir_output/STDOUT.txt, warnings/errors are still printed to console') # on/off flag, default is False
parser.add_argument('--version', action='version', version=f'%(prog)s {hatyan_version}')
args = parser.parse_args()


file_config = args.filename # os.path.realpath(sys.argv[0])

# check for file existence before creating folder
if not os.path.isfile(file_config):
    raise FileNotFoundError(f'file_config not found: {file_config}')

# create dir_output and chdir
foldername = os.path.basename(file_config).split('.')[0]
if args.unique_outputdir:
    dir_output = f"{foldername}__{timer_start.strftime('%Y%m%d_%H%M%S')}"
else:
    dir_output = foldername
os.makedirs(dir_output, exist_ok=True)
os.chdir(dir_output)

#set the storage location of interactive plots
matplotlib.rcParams["savefig.directory"] = dir_output

#set the matplotlib backend depending on the interactive_plots argument
if args.interactive_plots:
    try:
        plt.switch_backend('Qt5agg')
    except:
        raise Exception('Failed to switch to Qt5agg backend, check if you have X-forwarding enabled (and mesa-libGL and possibly other libraries installed) or use argument interactive_plots=False')


#copy configfile to cwd (=dir_output)
shutil.copyfile(file_config,os.path.basename(file_config))

#redirecting stdout, stderr is still printed to console
if args.redirect_stdout:
    sys.stdout = open('STDOUT.txt', 'w')

# print initialization
print("############### HATYAN INITALIZING ###############")
print("--------------------------------------------------")
print(f"hatyan-{hatyan_version}: RWS tidal analysis and prediction")
print(f"started at:  {timer_start_str}")
print(f"file_config: {file_config}")
print(f"dir_output:  {dir_output}")
print("--------------------------------------------------")


# run the configfile
with open(file_config) as f:
    exec(f.read())


# get stop time
timer_stop = dt.datetime.now()
timer_stop_str = timer_stop.strftime('%Y-%m-%d %H:%M:%S')
timer_elapsed = (timer_stop - timer_start).total_seconds()/60

# print de-initialization
print("--------------------------------------------------")
print(f"finished at: {timer_stop_str}")
print(f"elapsed time: {timer_elapsed:.2f} minutes {timer_elapsed*60:.2f} seconds")
print("--------------------------------------------------")
print("################# HATYAN FINISHED ################")

if matplotlib.get_backend() == 'Qt5Agg':
    print('WARNING: close open plots to continue')
    plt.show()

