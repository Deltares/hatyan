Name:        hatyan_python
Version:     main #2.2.90
Release:     1
#BuildArch:   noarch
#Buildroot:   ~/rpmbuild/%{name}-%{version}-root
URL:         https://github.com/Deltares/hatyan
AutoReq:     no
Summary:     Python version of the hatyan RWS program, packed with relocatable Python env including necessary Python libraries
License:     LGPL
Provides:    hatyan_python
Requires:    rh-python36-python >= 3.6.3 rh-python36-python-libs >= 3.6.3 rh-python36-python-virtualenv >= 3.6.3 glibc >= 2.12 coreutils expect stix-fonts fontconfig freetype libstdc++ jasper libXcursor libXrender xorg-x11-xauth mesa-libGL mesa-libEGL libXi

%description
%{summary}

#rpmbuild requires (sudo yum -y install): centos-release-scl-rh, rh-python36-python,  rh-python36-python-virtualenv, rpm-build
#start rpmbuild like this (default and more strict)
#rpmbuild -v -bb ~/hatyan_github/scripts/hatyan_python-latest.spec
#rpmbuild -v -bb ~/hatyan_github/scripts/hatyan_python-latest.spec --define "_topdir /u/veenstra/rpmbuild" --define "HATYANROOTFOLDER ~/hatyan_github"

#install the code into directories on the build machine
%install
echo TESTTEST
echo %{version}
echo TESTTEST
#make local copy of hatyan sources, to install from later. first all files in root (but not folders), then the hatyan and scripts folder
cp %{HATYANROOTFOLDER}/* %{_topdir}/SOURCES | true
cp -r %{HATYANROOTFOLDER}/hatyan %{_topdir}/SOURCES
#create sh script for running hatyan on linux in one command
mkdir -p $RPM_BUILD_ROOT/usr/bin
EXECFILE=$RPM_BUILD_ROOT/usr/bin/hatyan
cp %{HATYANROOTFOLDER}/scripts/hatyan.sh $EXECFILE
chmod +x $EXECFILE
#create folder for hatyan_env and potentially other folders/files
mkdir -p $RPM_BUILD_ROOT/opt/hatyan_python
cp -r %{HATYANROOTFOLDER}/doc $RPM_BUILD_ROOT/opt/hatyan_python
cp -r %{HATYANROOTFOLDER}/tests $RPM_BUILD_ROOT/opt/hatyan_python
#cp -r %{HATYANROOTFOLDER}/hatyan $RPM_BUILD_ROOT/opt/hatyan_python
# create empty virtual environment
/opt/rh/rh-python36/root/usr/bin/virtualenv $RPM_BUILD_ROOT/opt/hatyan_python/hatyan_env
# upgrade pip and setuptools to make sure all dependencies are handled well
$RPM_BUILD_ROOT/opt/hatyan_python/hatyan_env/bin/python -m pip install --upgrade pip
$RPM_BUILD_ROOT/opt/hatyan_python/hatyan_env/bin/python -m pip install --upgrade setuptools
#install hatyan package from source, also install old library versions to make it work on CentOS (prevent errors related to Qt and others)
$RPM_BUILD_ROOT/opt/hatyan_python/hatyan_env/bin/python -m pip install %{_topdir}/SOURCES -r %{_topdir}/SOURCES/requirements_dev.txt
#install pyqt5==5.7.1 to avoid "Failed to import any qt binding" error. The fixed version is necessary since CentOS/RHEL6 have glibc 2.12 and higher pyqt5 versions require glibc>=2.14
$RPM_BUILD_ROOT/opt/hatyan_python/hatyan_env/bin/python -m pip install pyqt5==5.7.1
#make existing environment relocatable and remove BUILDROOT prefix in activate file:
/opt/rh/rh-python36/root/usr/bin/virtualenv --relocatable $RPM_BUILD_ROOT/opt/hatyan_python/hatyan_env
sed -i "s#/.*/rpmbuild/BUILDROOT/.*x86_64##g" $RPM_BUILD_ROOT/opt/hatyan_python/hatyan_env/bin/activate
exit 0 #to prevent compiling

# gathers list of files and packs them to the RPM
%files
/opt/hatyan_python
/usr/bin/hatyan
