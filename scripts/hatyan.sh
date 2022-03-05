#!/bin/bash
set -e

# Help
display_help() {
    # Display Help
    echo "hatyan is a Python program for tidal analysis and prediction"
    echo "Usage: $0 [option...]" >&2
    echo "Usage: $0 [file]" >&2
    echo
    echo "   -h, --help              Display this help message"
    echo "   -t, --test              Run hatyan with example script"
    echo "   -T, --testmodules       Test if all linux libraries are installed"
    echo "   -v, --version           Display hatyan version"
    echo
    exit 1
}

source /opt/hatyan_python/hatyan_env/bin/activate

while :
do
    case "$1" in
      -h | --help)
          display_help  # Call your function
          exit 0
          ;;
      -t | --test)
          hatyan /opt/hatyan_python/tests/examples/predictie_2019_19Ycomp4Ydia_VLISSGN_interactive.py
          exit 0
          ;;
      -T | --testmodules)
        	echo "testing if all dependencies for hatyan are installed";
        	for LIBRARYNAME in glibc coreutils expect stix-fonts fontconfig freetype libstdc++ jasper libXcursor libXrender xorg-x11-xauth mesa-libGL mesa-libEGL libXi; do 
        		echo "--------------------------"
        		echo "searching for '$LIBRARYNAME' on this system";
        		#export LIBRARYVERSION=$(ldconfig -p | grep $LIBRARYNAME);
        		export LIBRARYVERSION=$(rpm -qa |grep ^$LIBRARYNAME-);
        		if [ -z "$LIBRARYVERSION" ]; then
        			echo "ERROR: '$LIBRARYNAME' is not available on this system";
        		else
        			echo "SUCCESS: '$LIBRARYNAME' is available on this system:";
        			echo "$LIBRARYVERSION";
        		fi;
        	done
          exit 0
          ;;
      -v | --version)
          echo hatyan-$(python -c "import hatyan; print(hatyan.__version__)")
          exit 0
          ;;
      -*)
          echo "Error: Unknown option: $1" >&2
          ## or call function display_help
          exit 1 
          ;;
      *)  # No more options
          break
          ;;
    esac
done

#main program
echo $(python --version)
pythonversion=$(python -c "import sys;print(f'{sys.version_info.major}.{sys.version_info.minor}')")

if [ ! -f $1 ]; then echo "ERROR: configfile not found."
else
	if [[ $pythonversion == 3.6 ]]; then
		export QT_QPA_PLATFORM_PLUGIN_PATH=/opt/hatyan_python/hatyan_env/lib/python3.6/site-packages/PyQt5/Qt/plugins/platforms
	else
		export QT_QPA_PLATFORM_PLUGIN_PATH=/opt/hatyan_python/hatyan_env/lib/python${pythonversion}/site-packages/PyQt5/Qt
	fi
	export QT_XKB_CONFIG_ROOT=/usr/share/X11/xkb
	echo hatyan-$(python -c "import hatyan; print(hatyan.__version__)") started
	dir_output=$(python -c "import sys; from hatyan.wrapper_RWS import init_RWS; dir_output, timer_start = init_RWS(file_config=sys.argv[1],silent=True,interactive_plots=False); print(dir_output)" $1)
	unbuffer python $1 $dir_output 2>&1 |tee $dir_output/FILE_DIAGNOSTICS.txt
fi
