#!/bin/bash
set -e

source /opt/hatyan_python/hatyan_env/bin/activate
if [ -z $1 ]; then echo "ERROR: no (path to) configfile --version or --test provided as argument."
elif [ "$1" == "--version" ]; then echo hatyan-$(python -c "import hatyan; print(hatyan.__version__)") 
elif [ "$1" == "--test" ]; then
	echo "testing if all dependencies for hatyan are installed";
	for LIBRARYNAME in rh-python36-python rh-python36-python-libs glibc coreutils expect stix-fonts fontconfig freetype libstdc++ jasper libXcursor libXrender xorg-x11-xauth mesa-libGL mesa-libEGL libXi; do 
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
elif [ ! -f $1 ]; then echo "ERROR: configfile not found."
else
	export QT_QPA_PLATFORM_PLUGIN_PATH=/opt/hatyan_python/hatyan_env/lib/python3.6/site-packages/PyQt5/Qt/plugins/platforms
	export QT_XKB_CONFIG_ROOT=/usr/share/X11/xkb
	echo hatyan-$(python -c "import hatyan; print(hatyan.__version__)") started
	dir_output=$(python -c "import sys; from hatyan.wrapper_RWS import init_RWS; dir_output, timer_start = init_RWS(file_config=sys.argv[1],silent=True,interactive_plots=False); print(dir_output)" $1)
	unbuffer python $1 $dir_output 2>&1 |tee $dir_output/FILE_DIAGNOSTICS.txt
fi
