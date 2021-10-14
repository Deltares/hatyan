#rpmbuild requires (sudo yum -y install): centos-release-scl-rh, rh-python36-python,  rh-python36-python-virtualenv, rpm-build. Start rpmbuild like this
#rpmbuild -v -bb ~/hatyan_github/scripts/hatyan_python-latest.spec --define "VERSIONTAG main"
#TODO: update rpmbuild requires list (rh is old, replace with python3*?)

Name:        hatyan_python
Version:     2.3.1
Release:     1
BuildArch:   x86_64
URL:         https://github.com/Deltares/hatyan
AutoReq:     no
Summary:     Python version of the hatyan RWS program, packed with relocatable Python env including necessary Python libraries
License:     LGPL
Provides:    hatyan_python
Requires:    python3 python3-libs python3-pip python3-setuptools glibc >= 2.12 coreutils expect stix-fonts fontconfig freetype libstdc++ jasper libXcursor libXrender xorg-x11-xauth mesa-libGL mesa-libEGL libXi
#TODO: update rpminstall requires list (python3 is shipped now?)
#TODO: check on teamcity?

%description
%{summary}

#prevent build_id links (/usr/lib/) on github platform
%define _build_id_links none
#define versiontag (defaults to main if not passed as rpmbuild define flag, there should be a github tag created with that name, e.g. v2.2.86)
%{!?VERSIONTAG: %define VERSIONTAG main}

#install the code into directories on the build machine
%install
#clear build folder, clone specific hatyan versiontag
rm -rf %{_topdir}/BUILD/*
git clone -b %{VERSIONTAG} https://github.com/Deltares/hatyan.git %{_topdir}/BUILD/hatyan_github #the BUILD folder is where it automatically clones to
#create sh script for running hatyan on linux in one command
mkdir -p $RPM_BUILD_ROOT/usr/bin
EXECFILE=$RPM_BUILD_ROOT/usr/bin/hatyan
cp %{_topdir}/BUILD/hatyan_github/scripts/hatyan.sh $EXECFILE
chmod +x $EXECFILE
#create folder for hatyan_env and potentially other folders/files
mkdir -p $RPM_BUILD_ROOT/opt/hatyan_python
cp -r %{_topdir}/BUILD/hatyan_github/doc $RPM_BUILD_ROOT/opt/hatyan_python
cp -r %{_topdir}/BUILD/hatyan_github/tests $RPM_BUILD_ROOT/opt/hatyan_python
#cp -r %{_topdir}/BUILD/hatyan_github/hatyan $RPM_BUILD_ROOT/opt/hatyan_python
# create python3 venv to install virtualenv in # this is necessary on h6-c7, not on Github Actions since default there is python3
python3 -m venv hatyan_setup_venv
. hatyan_setup_venv/bin/activate #Was (but does not work on github): source hatyan_setup_venv/bin/activate
python --version #TODO: this version is used for venv and virtualenv. Python 3.8 requires glibc>2.24 or so, but might not be available on destination machine, how to fix pythonversion? Github has setup-python available for specific version 
python -m pip install --upgrade pip setuptools
python -m pip install virtualenv #==15.1.0 #in this version the relocatable flag still works
#create empty virtualenv (this one should be relocatable, not possible with venv)
#/opt/rh/rh-python36/root/usr/bin/virtualenv $RPM_BUILD_ROOT/opt/hatyan_python/hatyan_env
python -m virtualenv $RPM_BUILD_ROOT/opt/hatyan_python/hatyan_env
deactivate #deactivate venv hatyan_setup_venv
# upgrade pip and setuptools to make sure all dependencies are handled well
$RPM_BUILD_ROOT/opt/hatyan_python/hatyan_env/bin/python -m pip install --upgrade pip setuptools
#install hatyan package from source, also install old library versions to make it work on CentOS (prevent errors related to Qt and others)
$RPM_BUILD_ROOT/opt/hatyan_python/hatyan_env/bin/python -m pip install %{_topdir}/BUILD/hatyan_github -r %{_topdir}/BUILD/hatyan_github/requirements_dev.txt
#install pyqt5==5.7.1 to avoid "Failed to import any qt binding" error. The fixed version is necessary since CentOS/RHEL6 have glibc 2.12 and higher pyqt5 versions require glibc>=2.14
$RPM_BUILD_ROOT/opt/hatyan_python/hatyan_env/bin/python -m pip install pyqt5==5.7.1
#make existing environment relocatable and remove BUILDROOT prefix in activate file:
#/opt/rh/rh-python36/root/usr/bin/virtualenv --relocatable $RPM_BUILD_ROOT/opt/hatyan_python/hatyan_env  #relocatable flag was dropped in newer virtualenv version which comes with non-rh python. It also broke installation entirely (no module hatyan)
sed -i "s#/.*/rpmbuild/BUILDROOT/.*x86_64##g" $RPM_BUILD_ROOT/opt/hatyan_python/hatyan_env/bin/* #replace BUILDROOT prefix in all bin folders, makes sure activate/python/pip can find each other
exit 0 #to prevent compiling

# gathers list of files and packs them to the RPM
%files
/opt/hatyan_python
/usr/bin/hatyan
